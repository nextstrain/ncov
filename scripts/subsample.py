from augur.utils import AugurException
from augur.filter import run as augur_filter, register_arguments as register_filter_arguments
from augur.index import index_sequences
from augur.io import write_sequences, open_file, read_sequences, read_metadata
from get_distance_to_focal_set import get_distance_to_focal_set # eventually from augur.priorities (or similar)
from priorities import create_priorities # eventually from augur.priorities (or similar)
import yaml
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from os import path
import pandas as pd
from tempfile import NamedTemporaryFile
import jsonschema
# from pkg_resources import resource_string

DESCRIPTION = "Subsample sequences based on user-defined YAML configuration"

def register_arguments(parser):
    parser.add_argument('--scheme', required=True, metavar="YAML", help="subsampling scheme")
    parser.add_argument('--output-dir', required=True, metavar="PATH", help="directory to save intermediate results")
    parser.add_argument('--metadata', required=True, metavar="TSV", help="metadata")
    parser.add_argument('--alignment', required=True, metavar="FASTA", help="alignment to subsample")
    parser.add_argument('--alignment-index', required=False, metavar="INDEX", help="sequence index of alignment")
    parser.add_argument('--reference', required=True, metavar="FASTA", help="reference (which was used for alignment)")
    parser.add_argument('--include-strains-file', required=False, nargs="+", default=None, metavar="TXT", help="strains to force include")
    parser.add_argument('--exclude-strains-file', required=False, nargs="+", default=None, metavar="TXT", help="strains to force exclude")
    parser.add_argument('--output-fasta', required=True, metavar="FASTA", help="output subsampled sequences")
    parser.add_argument('--output-metadata', required=True, metavar="TSV", help="output subsampled metadata")
    parser.add_argument('--output-log', required=False, metavar="TSV", help="log file explaining why strains were excluded / included")
    parser.add_argument('--use-existing-outputs', required=False, action="store_true", help="use intermediate files, if they exist")

def run(args):

    config = parse_scheme(args.scheme)

    generate_sequence_index(args)

    samples = [Sample(name, data, args) for name, data in config.items()]

    graph = make_graph(samples)

    traverse_graph(
        graph,
        lambda s: s.filter()
    )

    combine_samples(args, samples)

def parse_scheme(filename):
    with open(filename) as fh:
        try:
            data = yaml.safe_load(fh)
        except yaml.YAMLError as exc:
            print(exc)
            raise AugurException(f"Error parsing subsampling scheme {filename}")
    validate_scheme(data)
    return data


def validate_scheme(scheme):
    try:
        # When we move this to `augur subsample`, load the schema via:
        # schema = yaml.safe_load(resource_string(__package__, path.join("data", "schema-subsampling.yaml")))
        with open(path.join(path.dirname(path.realpath(__file__)), "subsample_schema.yaml")) as fh:
            schema = yaml.safe_load(fh)
    except yaml.YAMLError as err:
        raise AugurException("Subsampling schema definition is not valid YAML. Error: {}".format(err))
    # check loaded schema is itself valid -- see http://python-jsonschema.readthedocs.io/en/latest/errors/
    try:
        jsonschema.Draft6Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise AugurException("Subsampling schema definition is not valid. Error: {}".format(path, err))

    try:
        jsonschema.Draft6Validator(schema).validate(scheme)
    except jsonschema.exceptions.ValidationError as err:
        print(err)
        raise AugurException("Subsampling scheme failed validation")

class Sample():
    """
    A class to hold information about a sample. A subsampling scheme will consist of multiple
    samples. Each sample may depend on the priorities based off another sample.
    """
    def __init__(self, name, config, cmd_args):
        self.name = name
        self.tmp_dir = cmd_args.output_dir
        self.alignment = cmd_args.alignment
        self.alignment_index = cmd_args.alignment_index
        self.reference = cmd_args.reference
        self.metadata = cmd_args.metadata
        self.initialise_filter_args(config, cmd_args)
        self.priorities = config.get("priorities", None) 
        self.use_existing_outputs = args.use_existing_outputs

    def initialise_filter_args(self, config, subsample_args):
        """
        Currently this method is needed as we need to call `augur filter`'s `run()` with an
        argparse instance. An improvement here would be to expose appropriate filtering
        functions and call them as needed, with the output being returned rather than
        written to disk.
        """
        # create the appropriate command-line arguments for the augur filter run we want
        arg_list = [
            "--metadata", self.metadata,
            "--sequences", self.alignment,
            "--sequence-index", self.alignment_index,
            "--output", path.join(self.tmp_dir, f"sample.{self.name}.fasta"),         # filtered sequences in FASTA forma
            "--output-metadata", path.join(self.tmp_dir, f"sample.{self.name}.tsv"),  # metadata for strains that passed filters
            "--output-strains", path.join(self.tmp_dir, f"sample.{self.name}.txt"),   # list of strains that passed filters (no header)
            "--output-log", path.join(self.tmp_dir, f"sample.{self.name}.log.tsv")
        ]
        # convert the YAML config into the command-line arguments for augur filter
        for name, value in config.items():
            if isinstance(value, dict):
                pass # we explicitly ignore dictionary config entries
            elif isinstance(value, list):
                arg_list.append(f"--{name}")
                arg_list.extend([str(v) for v in value])
            elif isinstance(value, bool):
                if value:
                    arg_list.append(f"--{name}")
            else:
                arg_list.append(f"--{name}")
                arg_list.append(str(value))
        # mock an ArgumentParser so that we can use augur filters interface, avoiding the need to duplicate logic
        parser = ArgumentParser(prog="Mock_Augur_Filter")
        register_filter_arguments(parser)
        self.filter_args, unused_args = parser.parse_known_args(arg_list)
        if unused_args:
            print(f"Warning - the following config parameters are not part of augur filter and may be ignored:")
            print(' '.join(unused_args))

    def calculate_required_priorities(self):
        """
        If computation of this sample requires priority information of another sample
        (the "focus"), then this function will compute those priorities by calling
        a method on the focal sample object.
        """
        if not self.priorities:
            return
        focal_sample = self.priorities.get('sample', None)
        if not focal_sample:
            raise AugurException(f"Cannot calculate priorities needed for {self.name} as the {self.get_priority_focus_name()} sample wasn't linked")
        print(f"Calculating priorities of {focal_sample.name}, as required by {self.name}")
        priorities_file = focal_sample.calculate_priorities()
        print(f"\tSetting {self.name} filter priority file to {priorities_file}")
        self.filter_args.priority = priorities_file

    def calculate_priorities(self):
        """
        Calculate the priorities TSV file for samples in the alignment vs this sample

        Returns the filename of the priorities file (TSV)
        """

        proximity_output_file = path.join(self.tmp_dir, f"proximity_{self.name}.tsv")
        if self.use_existing_outputs and check_outputs_exist(proximity_output_file):
            print(f"Using existing proximity scores for {self.name}")
        else:
            print(f"Calculating proximity of {self.name}")
            get_distance_to_focal_set(
                self.alignment,
                self.reference,
                self.filter_args.output,
                proximity_output_file,
                ignore_seqs=["Wuhan/Hu-1/2019"]  # TODO - use the config to define this?
            )

        priorities_path = path.join(self.tmp_dir, f"priorities_{self.name}.tsv")
        if self.use_existing_outputs and check_outputs_exist(priorities_path):
            print(f"Using existing priorities for {self.name}")
        else:
            print(f"Calculating priorities of {self.name}")
            create_priorities(
                self.alignment_index,
                proximity_output_file,
                priorities_path
            )
        return priorities_path

    def get_priority_focus_name(self):
        if not self.priorities:
            return None
        return self.priorities['focus']

    def set_priority_sample(self, sample):
        if not self.priorities:
            raise AugurException(f"No priorities set for {self.name}")
        self.priorities['sample'] = sample

    def filter(self):
        print("\n---------------------------------\nCONSTRUCTING SAMPLE FOR", self.name, "\n---------------------------------")
        self.calculate_required_priorities()
        if self.use_existing_outputs and check_outputs_exist(self.filter_args.output_metadata, self.filter_args.output_strains, self.filter_args.output_log):
            print(f"Using existing filtering results for {self.name}")
        else:
            print("Calling augur filter")
            print("Filter arguments:")
            for k,v in self.filter_args.__dict__.items():
                if v is not None:
                    print(f"\t{k: <30}{v}")
            augur_filter(self.filter_args)

        # In the future, instead of `augur_filter` saving data to disk, it would return
        # data to the calling process. In lieu of that, we read the data just written.
        try:
            self.sampled_strains = set(pd.read_csv(self.filter_args.output_strains, header=None)[0])
        except pd.errors.EmptyDataError:
            self.sampled_strains = set()
        self.filter_log = pd.read_csv(
            self.filter_args.output_log,
            header=0,
            sep="\t",
            index_col="strain"
        )


def make_graph(samples):
    """"
    Given a config file, construct a graph of samples to perform in an iterative fashion, such that
    priorities 
    This is a DAG, however an extremely simple one which we can construct outselves rather than relying on
    extra libraries.
    Constraints:
    * Each sample can only use priorities of one other sample
    * Acyclic
    Structure:
    tuple: (sample name, list of descendent samples) where a "descendant" sample requires the linked sample to be
    created prior to it's creation. Each entry in the list has this tuple structure.
    """

    included = set() # set of samples added to graph
    graph = (None, [])

    # add all the samples which don't require priorities to the graph
    for sample in samples:
        if not sample.get_priority_focus_name():
            graph[1].append((sample, []))
            included.add(sample.name)

    def add_descendants(level):
        parent_sample = level[0]
        descendants = level[1]
        for sample in samples:
            if sample.name in included:
                continue
            if sample.get_priority_focus_name() == parent_sample.name:
                sample.set_priority_sample(parent_sample)
                descendants.append((sample, []))
                included.add(sample.name)
        for inner_level in descendants:
            add_descendants(inner_level)

    for level in graph[1]:
        add_descendants(level)

    # from pprint import pprint
    # print("\ngraph"); pprint(graph);print("\n")

    if len(samples)!=len(included):
        AugurException("Incomplete graph construction")

    return graph

def traverse_graph(level, callback):
    this_sample, descendents = level
    if this_sample:
        callback(this_sample)
    for child in descendents:
        traverse_graph(child, callback)

def generate_sequence_index(args):
    if args.alignment_index:
        print("Skipping sequence index creation as an index was provided")
        return    
    print("Creating ephemeral sequence index file")
    with NamedTemporaryFile(delete=False) as sequence_index_file:
        sequence_index_path = sequence_index_file.name
        index_sequences(args.alignment, sequence_index_path)
        args.alignment_index = sequence_index_path


def combine_samples(args, samples):
    """Collect the union of strains which are included in each sample and write them to disk.
    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from argparse
    samples : list[Sample]
        list of samples
    """
    print("\n\n")
    ### Form a union of each sample set, which is the subsampled strains list
    sampled_strains = set()
    for sample in samples:
        print(f"Sample \"{sample.name}\" included {len(sample.sampled_strains)} strains")
        sampled_strains.update(sample.sampled_strains)
    print(f"In total, {len(sampled_strains)} strains are included in the resulting subsampled dataset")

    ## Iterate through the input sequences, streaming a subsampled version to disk.
    sequences = read_sequences(args.alignment)
    sequences_written_to_disk = 0
    with open_file(args.output_fasta, "wt") as output_handle:
        for sequence in sequences:
            if sequence.id in sampled_strains:
                sequences_written_to_disk += 1
                write_sequences(sequence, output_handle, 'fasta')
    print(f"{sequences_written_to_disk} sequences written to {args.output_fasta}")

    ## Iterate through the metadata in chunks, writing out those entries which are in the subsample
    metadata_reader = read_metadata(
        args.metadata,
        id_columns=["strain", "name"], # TODO - this should be an argument
        chunk_size=10000 # TODO - argument
    )
    metadata_header = True
    metadata_mode = "w"
    metadata_written_to_disk = 0
    for metadata in metadata_reader:
        df = metadata.loc[metadata.index.intersection(sampled_strains)]
        df.to_csv(
            args.output_metadata,
            sep="\t",
            header=metadata_header,
            mode=metadata_mode,
        )
        metadata_written_to_disk += df.shape[0]
        metadata_header = False
        metadata_mode = "a"
    print(f"{metadata_written_to_disk} metadata entries written to {args.output_metadata}")

    ## Combine the log files (from augur filter) for each sample into a larger log file
    ## Format TBD
    ## TODO

def check_outputs_exist(*paths):
    for p in paths:
        if not (path.exists(p) and path.isfile(p)):
            return False
    return True

if __name__ == "__main__":
    # the format of this block is designed specifically for future transfer of this script
    # into augur in the form of `augur subsample`
    parser = ArgumentParser(
        usage=DESCRIPTION,
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    register_arguments(parser)
    args = parser.parse_args()
    run(args)