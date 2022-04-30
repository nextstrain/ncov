import json
import argparse
from collections import defaultdict
from copy import deepcopy
from operator import is_
from augur.dates import numeric_date, SUPPORTED_DATE_HELP_TEXT # requires augur 15.0 or above


def collect_terminals(basal_node):
    """
    Returns a list of nodes that have no children
    """
    terminals = []
    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        if len(node.get("children", []))==0:
            terminals.append(node)
        else:
            for child in node["children"]:
                nodes_to_process.append(child)
    return terminals

def count_descendant_tips(basal_node, all_tips):
    count = 0
    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        if node['name'] in all_tips:
            count+=1
        for child in node.get("children", []):
            nodes_to_process.append(child)
    return count

def prune_orphan_branches(basal_node, terminals):
    """
    Makes a single pass through the tree to remove
    all orphan branches (those which have no children but
    which were internal nodes in the input tree).
    Modifies the tree-structure in place.
    Returns True if any modifications were made.
    This function should be called multiple times (until
    the return value is False) in order to remove all orphans!
    """
    made_changes = False
    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        n_children_before = len(node.get('children', []))
        children = []
        for child in node.get('children', []):
            hanging = len(child.get('children', []))==0
            orphan_branch = hanging and child['name'] not in terminals
            if not orphan_branch:
                children.append(child)
        n_children_after = len(children)
        node['children'] = children
        if n_children_after!=n_children_before:
            made_changes = True
        if len(node['children'])==0:
            del node['children']
        else:
            for child in node["children"]:
                nodes_to_process.append(child)
    return made_changes

def is_node_after_cutoff(node, cutoff_numeric):
    return node['node_attrs']["num_date"]["value"] > cutoff_numeric

def prune_prior_tips(basal_node, tips, cutoff):
    """
    Remove tips before the <cutoff> from the tree, as well as
    removing any orphan branches which result from this.
    Modifies the tree-structure in-place.
    """
    print(f"\nPruning tips prior to {cutoff}")
    cutoff_numeric = numeric_date(cutoff)
    prune_count = [0] # a list as I couldn't get `global` to work as I expected. TODO.
    
    def include(n):
        if len(n.get('children', []))>0: # child is an internal node
            return True
        if is_node_after_cutoff(n, cutoff_numeric):
            return True
        prune_count[0]+=1 ## pruned tips
        return False

    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        if len(node.get("children", []))>0:
            node["children"] = [child for child in node['children'] if include(child)]
            if len(node['children'])==0:
                del node['children']
            else:
                for child in node["children"]:
                    nodes_to_process.append(child)
    print(f"\tRemoved {prune_count[0]}/{len(tips)} tips")

    while prune_orphan_branches(basal_node, tips):
        pass

def get_trait(node, trait):
    """Returns the node-attr (trait) value"""
    return node.get('node_attrs', {}).get(trait, {}).get('value', False)

def collect_descendants(basal_node):
    """
    Returns a list of nodes that have no children
    """
    descendants = []
    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        descendants.append(node)
        for child in node.get("children", []):
            nodes_to_process.append(child)
    return descendants

def prune_from_tree(basal_node, node_name_to_prune):
    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        if node_name_to_prune in [n['name'] for n in node.get('children', [])]:
            node['children'] = [n for n in node['children'] if n['name']!=node_name_to_prune]
            if len(node['children'])==0:
                del node['children']
            break
        for child in node.get("children", []):
            nodes_to_process.append(child)

def get_clades(basal_node, tips, clade_label, cutoff):
    """
    Find the set of nodes which best explain the tips observed in the tree
    based on the provided clade labels.
    1. We first find the nodes labelled under they <clade_label> key
    2. We then consider that some clades are nested within others. If a clade 'c' is
    within another clade ('C'), and 'c' is after the cutoff, then we do not consider 'c' a new clade.
    If 'c' is prior to the cutoff, then it is a separate clade and it is removed from 'C'.
    3. Clades without any tips (i.e. the are all before the cutoff) are removed
  
    Returns a dictionary with keys: clade names, values: root nodes (of each clade)
    """
    print("\nFinding labelled nodes (branches) in pruned tree:")
    cutoff_numeric = numeric_date(cutoff)
    clades = {}
    nodes_to_process = [basal_node]
    while nodes_to_process:
        node = nodes_to_process.pop()
        label_value = node.get("branch_attrs", {}).get("labels", {}).get(clade_label)
        if (label_value):
            clades[label_value]=deepcopy(node)
            print(f"\t{label_value} is {node['name']}")
        for child in node.get("children", []):
            nodes_to_process.append(child)

    descendants = {lab: set([n['name'] for n in collect_descendants(node)]) for lab, node in clades.items()}


    print("\nDropping labelled nodes if they are a descendent of another label we've found")
    print("(and the parent label is after the cutoff):")
    part_of_other_subtrees = set()
    for label_name in clades.keys():
        ## consider other labels ('query') and remove them if they are a descendant of 
        ## the `label_name` AND the query name branches after the cutoff point
        for query_clade, query_node in clades.items():
            if query_clade==label_name:
                continue
            if query_node['name'] in descendants[label_name] and is_node_after_cutoff(query_node, cutoff_numeric):
                print(f"\tRemoving {query_clade} as it's within the subtree of {label_name}")
                part_of_other_subtrees.add(query_clade)
    for c in part_of_other_subtrees:
        del clades[c]
    ## the remaining clades are to be displayed as separate subtrees, so ensure they are pruned from other subtrees
    ## (otherwise, because we've made copies of the tree, some clades would be displayed multiple times)
    for c1 in clades.keys():
        for c2 in [cc for cc in clades.keys() if cc!=c1]:
            prune_from_tree(clades[c1], clades[c2]['name'])

    for label_name, labelled_node in list(clades.items()):
        while prune_orphan_branches(labelled_node, tips):
            pass

    for label_name, labelled_node in list(clades.items()):
        if count_descendant_tips(labelled_node, tips)==0:
            print(f"\tRemoving {label_name} as it has zero tips (after other clades pruned out)")
            del clades[label_name]

    return clades

def label_basal_nodes(clades, label_key):
    """
    Create / update a label with the clade name across the provided nodes
    """
    print(f"\nSetting branch labels under key '{label_key}'")
    for clade_name, node in clades.items():
        if 'branch_attrs' not in node: node['branch_attrs']={}
        if 'labels' not in node['branch_attrs']: node['branch_attrs']['labels']={}
        if label_key in node['branch_attrs']['labels']:
            existing_name = node['branch_attrs']['labels'][label_key]
            if existing_name!=clade_name:
                raise Exception(f"{node['name']} had conflicting branch labels: {existing_name} vs. {clade_name}")
            print(f"\t{clade_name} is still on the same node (as the input tree)")
        else:
            node['branch_attrs']['labels'][label_key] = clade_name
            print(f"\t{clade_name} is now on node {node['name']}")

def combine_mutations(m1, m2):
    """
    Combine mutations from two subsequent branches (nodes).
    Reversions are removed, and multiple mutations are collapsed.
    """
    def pos(mut):
        return int(mut[1:-1])
    m = defaultdict(list)
    for key in set(m1.keys()) | set(m2.keys()):
        muts = m1.get(key, []) + m2.get(key, [])
        # sort via position
        muts.sort(key=pos)
        # prune out reversions
        muts = [mut for mut in muts if mut[0]!=mut[-1]]
        # store in return dict, and combine mutations at
        # same position (e.g. A123B + B123C is A123C)
        if len(muts)==0:
            continue
        m[key].append(muts[0])
        for idx in range(1, len(muts)):
            if pos(m[key][-1])==pos(muts[idx]):
                if m[key][-1][-1]!=muts[idx][0]: ## e.g. A123B + D123C
                    raise Exception(f"Mutations are inconsistent!: {m[key][-1]} vs {muts[idx]}")
                m[key][-1] = f"{m[key][-1][0]}{pos(muts[idx])}{muts[idx][-1]}"
            else:
                m[key].append(muts[idx])
    return m

def transfer_labels(parentLabels, childLabels, muts):
    labels = {k:v for k,v in parentLabels.items() if k!='aa'}
    for k,v in childLabels.items():
        if k!='aa':
            labels[k]=v # child label overwrites parent in case of conflict
    # aa label code copied from `augur export v2` src
    aa = {gene:data for gene, data in muts.items() if len(data) and gene!="nuc"}
    if aa:
        labels['aa'] = '; '.join("{!s}: {!s}".format(key,', '.join(val)) for (key,val) in aa.items())
    return labels

def move_node_down(basal_node, tips):
    """
    Walk down the tree (from the `basal_node`) while the nodes have n=1 children
    (and aren't themselves terminal, or the single child isn't terminal).
    This is useful when a clade definition may have a number of branches prior to the
    first actual bifurcation point.

    Returns a tuple with [0]: the node (either the basal node, or a descendant of it)
                         [1]: the move count
    """
    node = basal_node
    move_count=0
    while len(node.get("children", []))==1 and len(node['children'][0].get("children", []))>0:
        node = node['children'][0]
        move_count+=1
    return (node, move_count)


def collapse_branches(clade_name, basal_node):
    """
    Collapse parent-child branches where the parent has only 1 child.
    Mutations are combined.
    Other branch attrs & node attrs are combined, and in the case of
    conflicts the child value is taken.
    Returns the number of internal branches collapsed
    """
    terminals = set([n['name'] for n in collect_terminals(basal_node)])
    nodes_to_process = [basal_node]
    collapse_count = 0
    while nodes_to_process:
        node = nodes_to_process.pop()
        if len(node.get("children", []))==1 and len(node['children'][0].get("children", []))>0:
            child = node['children'][0]
            ## Move the child's node-attributes to this node, overwriting as necessary
            collapse_count+=1

            node['name'] = child['name'] if child['name'] in terminals else f"{clade_name}_{collapse_count:05}"
            for child_attr, value in child['node_attrs'].items():
                node['node_attrs'][child_attr] = value # overwrites existing trait on node
            ## Transfer branch attrs (mutations, labels)
            node['branch_attrs']['mutations'] = combine_mutations(node['branch_attrs'].get('mutations', {}), child['branch_attrs'].get('mutations', {}))
            labels = transfer_labels(node['branch_attrs'].get('labels', {}), child['branch_attrs'].get('labels', {}), node['branch_attrs']['mutations'])
            if labels:
                node['branch_attrs']['labels'] = labels
            ## Make the grandchildren the new children, thus "removing" the child from the tree
            if child['name'] in terminals:
                del node['children'] # just in case
                continue
            else:
                node['children']=child['children']
            ## add back to the stack...
            nodes_to_process.append(node)
        else:
            for child in node.get("children", []):
                nodes_to_process.append(child)
    return collapse_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Restrict the tips in the tree to those after the provided <cutoff>, and create subtrees\
            based on the clades formed by the provided <trait>",
        epilog="Script designed for Nextstrain SARS-CoV-2 builds, but should be useable for other datasets."
    )
    parser.add_argument("--input", required=True, metavar="JSON", help="Auspice JSON")
    parser.add_argument("--output", required=True, metavar="JSON", help="output (auspice) JSON")
    parser.add_argument("--cutoff", required=True, metavar="DATE", help=f"Date cutoff. Formats: {SUPPORTED_DATE_HELP_TEXT}")
    parser.add_argument("--clade-label", required=True, metavar="LABEL", help="Label key to partition tips into clades")
    args = parser.parse_args()

    with open(args.input, 'r') as fh:
        jsn = json.load(fh)

    tips = set([n['name'] for n in collect_terminals(jsn['tree'])])
    prune_prior_tips(jsn['tree'], tips, args.cutoff)
    clades = get_clades(jsn['tree'], tips, args.clade_label, args.cutoff)

    # for each of the identified clades, collapse (descendant) branches if they have n=1 children
    print("\nMoving labels down the tree if they are on branches with n=1 children & collapsing descendant internal nodes with n=1 children:")
    for clade_name, clade_node in clades.items():
        (new_clade_node, move_count) = move_node_down(clade_node, tips)
        collapse_count = collapse_branches(clade_name, new_clade_node)
        clades[clade_name] = new_clade_node
        print(f"\t{clade_name} moved label down {move_count} branches and collapsed {collapse_count} internal branches")

    label_basal_nodes(clades, args.clade_label)

    ## Change the single-input tree to a list of subtrees
    jsn['tree'] = list(clades.values())

    with open(args.output, 'w') as fh:
        print(f"\nWriting (auspice) JSON to {args.output}")
        jsn = json.dump(jsn, fh, indent=2)