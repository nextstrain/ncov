import argparse, os, glob, sys
from Bio import AlignIO, SeqIO, SeqFeature, Seq
import numpy
import pandas as pd
from augur.utils import write_json

# enable importation of bindingcalculator from subdir of this script
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from bindingcalculator import BindingCalculator

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate immune escape from input spike sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignment', type=str, required=False, help="spike AA alignment")
    parser.add_argument('--output', type=str, required=True, help="JSON of escape estimates for nodes")
    parser.add_argument('--output-csv', type=str, required=True, help="CSV of escape data")
    args = parser.parse_args()

    # Trevor added this hardcoded reference: it seems kind of hacky and should
    # maybe be done better.
    # from SouthAfrica/NHLS-UCT-LA-Z957/2022
    reference = "MFVFLVLLPLVSSQCVNLITRTQ---SYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLDVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLGRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLKRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNHNAQALNTLVKQLSSKFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*"
    length = len(reference)
    assert "X" not in reference, "reference has ambiguous amino acids"

    alignment = AlignIO.read(args.alignment, 'fasta')

    # get RBD mutations and mutated sites for each strain
    rbd_start = 331  # start of RBD
    rbd_end = 531  # end of RBD
    records = []
    for record in alignment:
        mutated_rbd_sites = []
        rbd_mutations = []
        for site in range(rbd_start, rbd_end):
            ref_aa = reference[site - 1]
            aa = str(record.seq[site - 1])
            if (ref_aa != aa) and (ref_aa != "-" and aa != "-") and (aa != "X"):
                mutated_rbd_sites.append(site)
                rbd_mutations.append(f"{ref_aa}{site}{aa}")
        records.append((record.name, mutated_rbd_sites, rbd_mutations))
    # convert to data frame
    diffs_df = pd.DataFrame(
        records,
        columns=["strain", "mutated RBD sites", "RBD mutations"],
    )

    # compute escape, which is negative log2 binding retained
    bindcalc = BindingCalculator(
        eliciting_virus="SARS-CoV-2",
        known_to_neutralize="Omicron BA.2",
    )
    log_base = 2
    escape_df = diffs_df.assign(
        binding_retained=lambda x: x["mutated RBD sites"].map(bindcalc.binding_retained),
        escape_score=lambda x: -numpy.log(x["binding_retained"]) / numpy.log(log_base),
    )

    # write to CSV, which is useful for debugging
    escape_df.to_csv(args.output_csv, index=False, float_format="%.3g")

    # output JSON with escape scores for each strain
    escape_dict = (
        escape_df
        .set_index("strain")
        [["escape_score"]]
        .to_dict(orient="index")
    )
    write_json({"nodes": escape_dict}, args.output)
