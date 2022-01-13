import argparse
import json
import re
from numpy import linspace
from math import floor

def adjust_coloring_for_epiweeks(dataset):
    """
    If an auspice JSON specifies a colouring with the key "epiweek" (case sensitive) then we create a categorical
    colorscale which evenly spaces the canonical nextstrain rainbow across the observed time window.

    NOTE: epiweek must be in CDC format ("YYYYMM") but this may be relaxed to include ISO format in the future.
    """
    EPIKEY="epiweek"
    try:
        (cidx, coloring) = [(i, c) for i, c in enumerate(dataset['meta'].get("colorings", [])) if c['key']==EPIKEY][0]
    except IndexError: # coloring doesn't define an epiweek
        return 

    # remove any duplicate coloring entries in the JSON to ensure the entry we edit is the one used by Auspice
    # (NOTE: this is augur bug https://github.com/nextstrain/augur/issues/719)
    dataset['meta']['colorings'] = [c for i,c in enumerate(dataset['meta']['colorings']) if not (c['key']==EPIKEY and i!=cidx)]

    # delay import to support older setups not using epiweeks package
    from epiweeks import Year, Week
    
    observed_values = set()
    def recurse(node):
        value = node.get("node_attrs", {}).get(EPIKEY, {}).get("value", False)
        if value:
            # we validate using both the epiweeks package and a regex (epiweeks will perform coercion of non-valid data into valid data)
            if not re.match(r'^(\d{4})(\d{2})$', value):
                raise(ValueError(f"Epiweek value {value} was not in format YYYYMM."))
            week = Week.fromstring(value, system="cdc") # raises ValueError if not valid
            observed_values.add(week)
        for child in node.get("children", []):
            recurse(child)
    try:
        recurse(dataset["tree"])
    except ValueError as e:
        print(str(e))
        print("Skipping color scale creation for epiweek.")
        return
    observed_values = sorted(list(observed_values))

    ## generate epiweeks across the entire observed range for color generation
    epiweeks = [ observed_values[0] ]
    while epiweeks[-1] < observed_values[-1]:
        epiweeks.append(epiweeks[-1]+1)
    ## generate rainbow colour scale across epiweeks.
    ## Since a "default" augur install does not include matplotlib, rather than interpolating between values in the scale
    ## we reuse them. This only applies when n(epiweeks)>30, where distinguising between colors is problematic anyway.
    rainbow = ["#511EA8", "#482BB6", "#4039C3", "#3F4ACA", "#3E5CD0", "#416CCE", "#447CCD", "#4989C4", "#4E96BC", "#559FB0", "#5DA8A4", "#66AE96", "#6FB388", "#7AB77C", "#85BA6F", "#91BC64", "#9DBE5A", "#AABD53", "#B6BD4B", "#C2BA46", "#CDB642", "#D6B03F", "#DDA83C", "#E29D39", "#E69036", "#E67F33", "#E56D30", "#E2592C", "#DF4428", "#DC2F24"]
    color_indicies = [floor(x) for x in linspace(0, len(rainbow), endpoint=False, num=len(epiweeks))]
    coloring['scale'] = [
        [epiweek.cdcformat(), rainbow[color_indicies[i]]]
        for i,epiweek in enumerate(epiweeks)
        if epiweek in observed_values
    ]
    ## auspice will order the legend according to the provided color scale, so there is no need to set
    ## `coloring['legend']` unless we want to restrict this for some reason.
    coloring['type'] = 'categorical' # force the scale type to be categorical

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove extraneous colorings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        input_json = json.load(f)

    keys_to_remove = ["genbank_accession", "gisaid_epi_isl"]

    fixed_colorings = []
    for coloring in input_json["meta"]["colorings"]:
        if coloring['key'] not in keys_to_remove:
            fixed_colorings.append(coloring)

    input_json["meta"]["colorings"] = fixed_colorings

    adjust_coloring_for_epiweeks(input_json)

    with open(args.output, 'w') as f:
        json.dump(input_json, f, indent=2)
