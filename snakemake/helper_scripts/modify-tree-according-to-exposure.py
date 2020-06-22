import argparse
import json
from functools import partial
import copy

def traverse(node, fun):
    """ Traverse the tree calling `fun(node)` on each node """
    fun(node)
    for child in node.get("children", []):
        traverse(child, fun)

def modify_branches(node):
    """ Traverse the tree calling `fun` on each child """
    for idx, child in enumerate(node.get("children", [])):
        node["children"][idx] = modify(child)
        modify_branches(child)


def modify(node):
    """
    Modify the branch to NODE, if NODE is a terminal node and the EXPOSURE_TRAIT
    of NODE differs from the SAMPLING_TRAIT.
    This creates an intermediate branch, with a single child, such that we have
    PARENT-INTERMEDIATE-NODE.
    """
    if "children" in node:
        return node
    trait = node["node_attrs"].get(SAMPLING_TRAIT, {}).get("value", "")
    recoded = node["node_attrs"].get(EXPOSURE_TRAIT, {}).get("value", "")
    if not trait or not recoded or trait == recoded:
        return node


    ## Make a new node which is the INTERMEDIATE (i.e. a single child, which is `node` with slight modifications)
    n = {
        "name": node["name"]+"_travel_history",
        "node_attrs": {
            "div": node["node_attrs"]["div"],
            "num_date": node["node_attrs"]["num_date"],
            EXPOSURE_TRAIT: node["node_attrs"][EXPOSURE_TRAIT]
        },
        "children": [node]
    }

    # Store a backup of the original exposure, to put back in later
    node["node_attrs"]["exposure_backup"] = node["node_attrs"][EXPOSURE_TRAIT]
    # change actual NODE EXPOSURE_TRAIT to have SAMPLING_TRAIT
    node["node_attrs"][EXPOSURE_TRAIT] = node["node_attrs"][SAMPLING_TRAIT]
    return n


def collect_exposures(values, node):
    v = node["node_attrs"].get(SAMPLING_TRAIT, {}).get("value", None)
    if v:
        values.add(v)

def switch_attribute(node):
    """
    For a given node, replace SAMPLING_TRAIT with the info in EXPOSURE_TRAIT (removing the latter)
    """
    exposure = node["node_attrs"].get(EXPOSURE_TRAIT, None)
    sampling = node["node_attrs"].get(SAMPLING_TRAIT, None)
    if (sampling and not exposure):
        raise Exception("Where there's SAMPLING_TRAIT we should always have EXPOSURE_TRAIT")
    if exposure:
        node["node_attrs"][SAMPLING_TRAIT] = copy.copy(exposure)
        # if has a 'real' exposure history, then put this back in place now! Delete all others.
        if "exposure_backup" in node["node_attrs"] and node["node_attrs"]["exposure_backup"]["value"] != node["node_attrs"][EXPOSURE_TRAIT]["value"]:
            node["node_attrs"][EXPOSURE_TRAIT]["value"] = node["node_attrs"]["exposure_backup"]["value"]
            del node["node_attrs"]["exposure_backup"]
        else:
            del node["node_attrs"][EXPOSURE_TRAIT]

def reset_colors(colorings, values_wanted, colors):
    """
    Reset the color scale for SAMPLING_TRAIT so that it's in the correct order
    """
    trait = [x for x in colorings if x["key"]==SAMPLING_TRAIT][0]
    trait["scale"] = []
    values_wanted_lower = {x.lower() for x in values_wanted}
    values_added = set()
    for [name, value] in colors:
        if name.lower() in values_wanted_lower:
            trait["scale"].append([name, value])
            values_added.add(name.lower())
        else:
            print("Heads up: a colour has been set for SAMPLING_TRAIT -> {} but it is not found in the data!".format(name))
    for missing_name in values_wanted_lower - values_added:
        print("WARNING: Colour for {} -> {} is missing & auspice will choose a shade of grey for this.".format(SAMPLING_TRAIT, missing_name))


def update_latlongs(geo_resolutions, values_wanted, latlongs):
    """Add in values to the geo_resolutions section of the JSON from the latlongs"""
    trait = [x for x in geo_resolutions if x["key"]==SAMPLING_TRAIT][0]
    for x in values_wanted:
        if x not in trait["demes"].keys():
            try:
                trait["demes"][x] = latlongs[x]
            except KeyError:
                print("WARNING: The lat/long value for {} -> {} is missing from the lat/longs TSV. Please add it!".format(SAMPLING_TRAIT, x))
                # sys.exit(2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Modify tree based on exposure",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--colors', type=str, metavar="TSV", required=True, help="input colors TSV")
    parser.add_argument('--lat-longs', type=str, metavar="TSV", required=True, help="input lat/longs TSV")
    parser.add_argument('--sampling', type=str, required=True, help="name of sampling location attribute")
    parser.add_argument('--exposure', type=str, required=True, help="name of exposure location attribute")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    SAMPLING_TRAIT = args.sampling
    EXPOSURE_TRAIT = args.exposure

    with open(args.input) as fh:
        input_json = json.load(fh)

    # Read colorings & lat-longs
    colors = []
    with open(args.colors) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) == 3 and fields[0] == SAMPLING_TRAIT:
                colors.append([fields[1], fields[2]])
    latlongs = {}
    with open(args.lat_longs) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) == 4 and fields[0] == SAMPLING_TRAIT:
                latlongs[fields[1]] = {"latitude": float(fields[2]), "longitude": float(fields[3])}

    # Modify branches
    modify_branches(input_json["tree"])

    # And now move EXPOSURE_TRAIT to `sampling` (i.e. what was `sampling` is gone)
    traverse(input_json["tree"], switch_attribute)

    # Collect all EXPOSURE_TRAIT values (we may have created ones which didn't previously exist)
    sampling_values = set()
    traverse(input_json["tree"], partial(collect_exposures, sampling_values))

    # Ensure the colors & lat-longsare up-to-date
    reset_colors(input_json["meta"]["colorings"], sampling_values, colors)
    update_latlongs(input_json["meta"]["geo_resolutions"], sampling_values, latlongs)

    # Remove EXPOSURE_TRAIT from filters - we only included it so that `augur_export` would include the trait
    # on the nodes.
    input_json["meta"]["filters"] = [f for f in input_json["meta"]["filters"] if f!=EXPOSURE_TRAIT]

    with open(args.output, 'w') as fh:
        json.dump(input_json, fh, indent=2)
