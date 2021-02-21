"""
Based on a rotation project by Benjamin Otter @benjaminotter in neherlab

https://github.com/neherlab/sequence_distances
"""
import argparse, io, os
from treetime import TreeAnc
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO, Phylo
import subprocess
import json
import csv

A=ord('A')
C=ord('C')
G=ord('G')
T=ord('T')
N=ord('N')
gap=ord('-')

def sequence_to_int_array(s, fill_value=N, fill_gaps=True):
    seq = np.frombuffer(str(s).upper().encode('utf-8'), dtype=np.int8).copy()
    if fill_gaps:
        seq[(seq!=A) & (seq!=C) & (seq!=G) & (seq!=T)] = fill_value
    else:
        seq[(seq!=A) & (seq!=C) & (seq!=G) & (seq!=T) & (seq!=gap)] = fill_value
    return seq


def get_distance(n, s_diffs, Ns):
    """
    Parameters
    ----------
    n : Clade
        Clade of the focal alignment tree.
    s_diffs : set of tuples (int, str)
        differences of input sequence to reference sequence.
    Ns : set of int
        positions in input sequence which contain ambiguous bases or gap characters.

    Returns
    -------
    int
        number of symmetric differences

    """
    d = [x for x in s_diffs.symmetric_difference(n.differences) if x[0] not in Ns]
    return len(d)

def get_differences(s, ref):
    """
    Parameters
    ----------
    s : array of str
        input sequence.
    ref : array of str
        reference sequence.

    Returns
    -------
    set of differences, set of Ns
        set of differences (int, str) between the input sequence and the reference sequence ,
        contains the position and the mutated base of the input sequence.
        set of Ns (int) contains the positions in the input sequence with ambiguous or gap
        characters.

    """
    good_idx = s!=N
    differences = np.where((ref!=s)&good_idx)[0]
    return set(zip(differences, s[differences])), set(np.where(~good_idx)[0])

def remove_close_children(T):
    """

    Parameters
    ----------
    T : Bio.Phylo.Newick.Tree
        phylogenetic tree to prune similar Clades.

    Returns
    -------
    T : Bio.Phylo.Newick.Tree
        phylogenetic tree where all similar Clades are pruned.

    """
    n_kids = 0
    for n in T.find_clades(order='postorder'):
        identical_kids = []
        for c in n.clades:
            #append identical kids to list
            if (c.is_terminal() and len((n.differences).symmetric_difference(c.differences))==0):
                identical_kids.append(c)
        #sort by ascending N values
        identical_kids.sort(key=lambda x: len(x.Ns))
        #remove all but sequence with lowest amount of Ns
        for k in identical_kids[1:]:
            n_kids+=1
            T.prune(k)
    print('Removed %i identical children from tree.'%(n_kids))
    return T

def find_closest(T, s_diffs, Ns, margin):
    """
    Iterative algorithm to find closest Clade in phylogenetic tree of focal-alignment
    given input sequence data and a margin value which defines the search radius.
    Low margin --> greedy search

    Parameters
    ----------
    T : Bio.Phylo.Newick.Tree
        phylogenetic tree of focal-alignment.
    s_diffs : set of tuples (int, str)
        differences of input sequence to reference sequence.
    Ns : set of int
        positions in input sequence which contain ambiguous bases or gap characters.
    margin : int
        search radius.

    Returns
    -------
    best : Bio.Phylo.Newick.Clade
        Clade with lowest distance to input sequence data.
    best_d : int
        distance of best Clade compared to input sequence.

    """
    best = T.root
    best_d = get_distance(best, s_diffs, Ns)
    children = best.clades

    while(len(children)):
        #get difference values for all nodes to compare
        dists = np.array([get_distance(n, s_diffs, Ns) for n in children])
        min_d = dists.min()

        #if there is a node that has better distance value, set as new best
        if min_d < best_d:
            best = children[np.argmin(dists)]
            best_d = min_d

        #if the distance is 0 return best node
        if best_d == 0:
            return best, best_d

        #select nodes that are within margin around current best node
        min_nodes = [children[idx] for idx in np.where(dists<(best_d+margin))[0]]

        #there are nodes within margin --> continue search for closest node
        children = [c for l in [n.clades for n in min_nodes] for c in l]

    #return best child
    return best, best_d

def find_closest_brute(T, s_diffs, Ns, margin):
    """
    brute_force algorithm to find closest Clade in phylogenetic tree of focal-alignment
    given input sequence data

    Parameters
    ----------
    T : Bio.Phylo.Newick.Tree
        phylogenetic tree of focal-alignment.
    s_diffs : set of tuples (int, str)
        differences of input sequence to reference sequence.
    Ns : set of int
        positions in input sequence which contain ambiguous bases or gap characters.
    margin : int
        search radius (not used in this method)

    Returns
    -------
    best : Bio.Phylo.Newick.Clade
        Clade with lowest distance to input sequence data.
    best_d : int
        distance of best Clade compared to input sequence

    """
    best = T.root
    d = get_distance(best, s_diffs, Ns)

    for n in T.find_clades():
        new_d = get_distance(n,s_diffs,Ns)
        if new_d<d:
            d = new_d
            best=n
    return best, d

def export_tree(T, metadata):
    """
    exports a tree and the metadata (s_diffs, Ns) for each Clade to
    a tree.nwk and metadata.json file

    Parameters
    ----------
    T : Bio.Phylo.Newick.Tree
        phylogenetic tree to export.
    metadata : dict
        dictionary that contains the s_diffs and Ns values
        for each Clade in the tree.

    Returns
    -------
    None.

    """
    if not os.path.isdir('exports'):
        os.mkdir('exports')

    Phylo.write(T, 'exports/tree.nwk', 'newick')
    with open('exports/metadata.json', 'w') as out_file:
        json.dump(metadata, out_file)

def build_tree(focal_alignment):
    """
    Parameters
    ----------
    focal_alignment : path-like
        path to a fasta file containing the sequences to build the
        focal-alignment tree.

    Returns
    -------
    Bio.Phylo.Newick.Tree
        pyhlogenetic tree of focal-alignment.

    """
    tree_cmd = ["FastTree", "-fastest", "-nt",  "-noml", "-nome",  "-nosupport", focal_alignment]

    T = Phylo.read(io.StringIO(subprocess.check_output(tree_cmd).decode()), 'newick')
    T.root_with_outgroup('Wuhan/Hu-1/2019')
    #T.root_at_midpoint()

    tt = TreeAnc(tree=T, aln=focal_alignment)
    tt.infer_ancestral_sequences(reconstruct_tip_states=False)
    tt.prune_short_branches()
    tt.optimize_tree()
    return tt.tree

def align(T, alignment, ref, method, margin, output):
    """
    Finds closest Clade in T for each sequence in alignment using the specified
    search function.
    Writes the results to a .tsv file with columns ['strain','distance','closest strain']
    (closest strain) is only saved if it is a terminal node.

    Parameters
    ----------
    T : Bio.Phylo.Newick.Tree
        phylogenetic tree of focal-alignment.
    alignment : path-like
        path to fasta file containing the sequences to search for closest Clade.
    ref : path like
        path to fasta file conatining the reference sequence.
    method : function
        function to search for closest Clade in tree.
    margin : int
        search radius for iterative search function.
    output : path-like
        path for output.

    Returns
    -------
    None.

    """
    i = 0
    print('Finding closest node:')
    with open(output, 'w', newline='') as out_file:

        tsv_writer = csv.writer(out_file, delimiter = '\t')
        header = ['strain','distance','closest strain']
        tsv_writer.writerow(header)

        with open(alignment) as fh:
            for l,s in SimpleFastaParser(fh):
                i+=1
                if i% 100==0:
                    print('%i sequences done'%(i))

                s_diffs, Ns = get_differences(sequence_to_int_array(s), ref)
                n,d = method(T,s_diffs,Ns,margin)

                row = [l,d,n.name]

                tsv_writer.writerow(row)
    print('Table of closest nodes has been exported to %s.'%(output))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on genetic proximity to focal sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--focal-alignment", type = str, required=True, help="focal smaple of sequences")
    parser.add_argument("--reference", type = str, required=True, help="reference sequence")
    parser.add_argument("--alignment", type=str, required=False, help="FASTA file of alignment")
    parser.add_argument("--tree", type=str, required = False, help="nwk file of focal alignment tree")
    parser.add_argument("--metadata", type=str, required = False, help="json file of focal alignment tree metadata")
    parser.add_argument("--export", action="store_true", default=False , required=False, help="export tree and metadata to tree.nwk and metadata.json")
    parser.add_argument("--method", type=str, default='recursive', required=False, choices=['brute', 'recursive'], help="choose between brute force search or recursive search")
    parser.add_argument("--margin", type=int, default = 6 ,required = False, help="margin for iterative search")
    parser.add_argument("--output", type=str, required=False, help="table with closest sequences")
    args = parser.parse_args()


    ref = sequence_to_int_array(SeqIO.read(args.reference, 'fasta').seq)

    #load tree or build tree
    if args.tree is not None:
        tree = Phylo.read(args.tree, "newick")
    else:
        tree = build_tree(args.focal_alignment)

    #load metadata or get metadata from Tree
    if args.metadata is not None:
        with open(args.metadata, 'r') as fh:
            metadata = json.load(fh)
        for n in T.find_clades():
            d, Ns = metadata[n.name]
            n.differences=set([(int(p),b) for p,b in d.items()])
            n.Ns = set(Ns)
    else:
        metadata = {}
        for n in tree.find_clades():
            n.differences, n.Ns = get_differences(sequence_to_int_array("".join(n.sequence)), ref)
            metadata[n.name]=[dict(n.differences),list(n.Ns)]

    #export metadata and tree
    if args.export==True:
        export_tree(tree, metadata)
        print("Tree and metadata saved in 'exports/' folder")

    #find closest node with specified method
    p = remove_close_children(tree)
    if args.method == 'recursive':
        method = find_closest
    elif args.method == 'brute':
        method = find_closest_brute

    align(p, args.alignment, ref, method, args.margin, args.output)
