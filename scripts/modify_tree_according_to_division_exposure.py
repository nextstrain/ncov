import json
import sys
from functools import partial

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

    # change actual NODE division_exposure to have collection division 
    node["node_attrs"]["division_exposure"] = node["node_attrs"]["division"]
    return n


def collect_division_exposures(values, node):
    v = node["node_attrs"].get("division_exposure", {}).get("value", None)
    if v:
        values.add(v)


def update_colors(colorings, values_wanted, colors):
    trait = [x for x in colorings if x["key"]=="division_exposure"][0]
    values_present = {x[0] for x in trait["scale"]}
    for x in values_wanted:
        if x not in values_present:
            print(x, colors[x])
            trait["scale"].append([x, colors[x]])

def update_latlongs(geo_resolutions, values_wanted, latlongs):
    trait = [x for x in geo_resolutions if x["key"]=="division_exposure"][0]
    for x in values_wanted:
        if x not in trait["demes"].keys():
            print(x, latlongs[x])
            trait["demes"][x] = latlongs[x]

if __name__ == '__main__':

    with open(sys.argv[1]) as fh:
        input_json = json.load(fh)

    # Read colorings & lat-longs
    colors = {}
    with open(sys.argv[2]) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields)==3 and fields[0]=="division_exposure":
                colors[fields[1]] = fields[2]
    latlongs = {}
    with open(sys.argv[3]) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields)==4 and fields[0]=="division":
                latlongs[fields[1]] = {"latitude": float(fields[2]), "longitude": float(fields[3])}

    # Modify branches
    modify_branches(input_json["tree"])

    # Collect all `division_exposure` values (we may have created ones which didn't previously exist)
    division_exposure_values = set()
    traverse(input_json["tree"], partial(collect_division_exposures, division_exposure_values))

    # Ensure the colors & lat-longsare up-to-date
    update_colors(input_json["meta"]["colorings"], division_exposure_values, colors)
    update_latlongs(input_json["meta"]["geo_resolutions"], division_exposure_values, latlongs)

    with open(sys.argv[4], 'w') as fh:
        json.dump(input_json, fh, indent=2)


