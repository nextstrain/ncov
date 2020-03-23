import json
import sys
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
    Modify the branch to NODE, if NODE is a terminal node and the `division_exposure`
    of NODE differs from the `division`.
    This creates an intermediate branch, with a single child, such that we have
    PARENT-INTERMEDIATE-NODE.
    """
    if "children" in node:
        return node
    division = node["node_attrs"].get("division", {}).get("value", "")
    recoded = node["node_attrs"].get("division_exposure", {}).get("value", "")
    if not division or not recoded or division == recoded:
        return node


    ## Make a new node which is the INTERMEDIATE (i.e. a single child, which is `node` with slight modifications)
    n = {
        "name": node["name"]+"_travel_history",
        "node_attrs": {
            "div": node["node_attrs"]["div"],
            "num_date": node["node_attrs"]["num_date"],
            "division_exposure": node["node_attrs"]["division_exposure"]
        },
        "children": [node]
    }

    # Story a backup of the original exposure, to put back in later
    node["node_attrs"]["division_exposure_backup"] = node["node_attrs"]["division_exposure"]
    # change actual NODE division_exposure to have collection division
    node["node_attrs"]["division_exposure"] = node["node_attrs"]["division"]
    return n


def collect_division_exposures(values, node):
    v = node["node_attrs"].get("division", {}).get("value", None)
    if v:
        values.add(v)

def switch_division(node):
    """
    For a given node, replace "division" with the info in "division_exposure" (removing the latter)
    """
    division_exposure = node["node_attrs"].get("division_exposure", None)
    division = node["node_attrs"].get("division", None)
    if (division and not division_exposure):
        raise Exception("Where there's division we should always have division_exposure")
    if division_exposure:
      node["node_attrs"]["division"] = copy.copy(division_exposure)
      # if has a 'real' exposure history, then put this back in place now! Delete all others.
      if "division_exposure_backup" in node["node_attrs"] and node["node_attrs"]["division_exposure_backup"]["value"] != node["node_attrs"]["division_exposure"]["value"]:
        node["node_attrs"]["division_exposure"]["value"] = node["node_attrs"]["division_exposure_backup"]["value"]
        del node["node_attrs"]["division_exposure_backup"]
      else:
        del node["node_attrs"]["division_exposure"]

def reset_colors(colorings, values_wanted, colors):
    """
    Reset the color scale for `division` so that it's in the correct order
    """
    trait = [x for x in colorings if x["key"]=="division"][0]
    trait["scale"] = []
    values_wanted_lower = {x.lower() for x in values_wanted}
    values_added = set()
    for [name, value] in colors:
        if name.lower() in values_wanted_lower:
            trait["scale"].append([name, value])
            values_added.add(name.lower())
        else:
            print("Heads up: a colour has been set for division -> {} but it is not found in the data!".format(name))
    for missing_name in values_wanted_lower - values_added:
        print("WARNING: Colour for division -> {} is missing & auspice will choose a shade of grey for this.".format(missing_name))


def update_latlongs(geo_resolutions, values_wanted, latlongs):
    """Add in values to the geo_resolutions section of the JSON from the latlongs"""
    trait = [x for x in geo_resolutions if x["key"]=="division"][0]
    for x in values_wanted:
        if x not in trait["demes"].keys():
            try:
                trait["demes"][x] = latlongs[x]
            except KeyError:
                print("WARNING: The lat/long value for division -> {} is missing from the lat/longs TSV. Please add it!".format(x))
                # sys.exit(2)

if __name__ == '__main__':

    with open(sys.argv[1]) as fh:
        input_json = json.load(fh)

    # Read colorings & lat-longs
    colors = []
    with open(sys.argv[2]) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields)==3 and fields[0]=="division":
                colors.append([fields[1], fields[2]])
    latlongs = {}
    with open(sys.argv[3]) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields)==4 and fields[0]=="division":
                latlongs[fields[1]] = {"latitude": float(fields[2]), "longitude": float(fields[3])}

    # Modify branches
    modify_branches(input_json["tree"])

    # And now move `division_exposure` to `division` (i.e. what was `division` is gone)
    traverse(input_json["tree"], switch_division)

    # Collect all `division_exposure` values (we may have created ones which didn't previously exist)
    division_values = set()
    traverse(input_json["tree"], partial(collect_division_exposures, division_values))

    # Ensure the colors & lat-longsare up-to-date
    reset_colors(input_json["meta"]["colorings"], division_values, colors)
    update_latlongs(input_json["meta"]["geo_resolutions"], division_values, latlongs)

    # Remove `division_exposure` from filters - we only included it so that `augur_export` would include the trait
    # on the nodes. 
    input_json["meta"]["filters"] = [f for f in input_json["meta"]["filters"] if f!="division_exposure"]

    with open(sys.argv[4], 'w') as fh:
        json.dump(input_json, fh, indent=2)


