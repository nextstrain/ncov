from os import listdir
from pathlib import Path
import argparse
from collections import defaultdict

def bold(s):
    return('\033[1m' + s + '\033[0m')

###################################
# Import functions
###################################

# Utils functions to read and add to very simple files containing no tabs, only listing of one entry in every line
def read_simple_file(name):
    with open(name) as myfile:
        data_file = myfile.readlines()
    return [l.strip() for l in data_file]

def read_dict(name):
    with open(name) as myfile:
        data_file = myfile.readlines()
    return {l.split("\t")[0]:l.split("\t")[1].strip() for l in data_file}

def add_to_simple_file(file_name, line):
    with open(file_name) as myfile:
        data_file = myfile.readlines()
    if not data_file[-1].endswith("\n"):
        line = "\n" + line
    with open(file_name, "a") as myfile:
        myfile.write(line + "\n")

# Return the current metadata as dictionary sorted by gisaid_epi_isl, only containing the sequences in sorted_info
def read_metadata(file_name, additional_info):
    with open(file_name) as myfile:
        data_file = myfile.readlines()

    metadata = {}

    header = data_file[0].strip().split("\t")
    for line in data_file[1:]:
        l = line.strip().split("\t")
        id = l[2]
        if id in additional_info:
            metadata[id] = {}
            for i in range(len(l)):
                type = header[i]
                if type == "gisaid_epi_isl":
                    continue
                metadata[id][type] = l[i]

    for id in additional_info:
        if id not in metadata:
            print(bold("WARNING: " + id + " missing from metadata! Please download most recent metadata before running this script."))
            return None

    return metadata


# Read ordering file and return as hierarchical dictionary:
def read_ordering_file(file_name):
    with open(file_name) as myfile:
        data_file = myfile.readlines()

    data = {"Asia": {}, "Oceania": {}, "Africa": {}, "Europe": {}, "South America": {}, "North America": {}}

    region = ""
    country = ""
    division = ""

    for line in data_file:
        if line == "\n":
            continue
        if line.startswith("###"):
            if len(line.split("### ")) > 1: #country
                place = line.strip().split("### ")[1]
                country = place
                if country not in data[region]:
                    data[region][country] = {}

        else:
            if line.startswith("#"):
                if len(line.split("# ")) > 1: #region or division
                    place = line.strip().split("# ")[1]
                    if place in data:
                        region = place
                    else:
                        division = place
                        if division not in data[region][country]:
                            data[region][country][division] = []

            else:
                l = line.strip().split("\t")
                type = l[0] #location, division etc
                place = l[1]
                if type == "division":
                    division = place
                    if division not in data[region][country]:
                        data[region][country][division] = []
                if type == "location":
                    location = place
                    if location not in data[region][country][division]:
                        data[region][country][division].append(location)

    return data


# cut a string apart at a ":" into key (e.g. "gisaid_epi_isl") and content (e.g. "EPI_ISL_537730")
def cut(s):
    key = s.split(":")[0]
    content = ":".join(s.split(":")[1:])[1:]
    return (key, content)

# Read all files within a given folder that start with "additional-info-changes" and return a dictionary first by
# GISAID epi isl, then by type of information (e.g. "additional_host_info")
def read_data(path):

    additional_info = {}

    for file in sorted(listdir(path)):
        if file == '.DS_Store':
            continue

        id = ""
        if file.startswith("additional-info-changes"):
            with open(path + file) as f:
                data = f.readlines()

            added = False #only consider newly added additional info
            removed = False

            for i in range(len(data)):
                k = data[i].strip("\n")

                if k.endswith("info added"):
                    added = True
                    removed = False
                if k.endswith("info changed"): #skip changed info
                    added = False
                    removed = False
                if k.endswith("info removed"): #consider removed info
                    added = False
                    removed = True

                if ":" in k:
                    if added or removed:
                        (key, content) = cut(k)
                        key = key.strip()

                        if key == "gisaid_epi_isl":
                                id = content
                                if added:
                                    if id in additional_info:
                                        print("WARNING: additional info added two times for same strain! (" + id + ")")
                                    additional_info[id] = {}
                        else:
                            if added:
                                additional_info[id][key] = content
                            if removed:
                                if id in additional_info:
                                    if key in additional_info[id]:
                                        additional_info[id].pop(key)
                                        if additional_info[id] == {}:
                                            additional_info.pop(id)

    return additional_info

# Rearrange the previously read info into a dictionary that stores for every type of unique additional info
# (e.g. "Imported from France") all sequences that came with this info.
def rearrange_additional_info(additional_info):
    sorted_info = {}

    for id in additional_info:
        info_found = 0
        for key in additional_info[id]:

            #Special case:
            if additional_info[id]["additional_location_info"] == "Migrants ship" and "additional_host_info" in additional_info[id] and additional_info[id]["additional_host_info"] != "":
                content = additional_info[id]["additional_location_info"] + " " + additional_info[id]["additional_host_info"]
                print("Merge two additional info (host and location) to " + bold(content))
                if content not in sorted_info:
                    sorted_info[content] = []
                sorted_info[content].append((id, additional_info[id]["strain"]))

            elif key == "additional_host_info" or key == "additional_location_info":
                content = additional_info[id][key]
                if content == "" or content == " ":
                    continue
                info_found += 1
                access = (content, key)
                if access not in sorted_info:
                    sorted_info[access] = []
                sorted_info[access].append((id, additional_info[id]["strain"]))
        if info_found > 1:
            if additional_info[id]["additional_location_info"] != additional_info[id]["additional_host_info"]:
                print("Warning: " + id + " has more than one relevant info (\"" + additional_info[id]["additional_host_info"] + "\" and \"" + additional_info[id]["additional_location_info"] + "\"). Possible conflict!")

    return sorted_info



###################################
# Comparison and search functions
###################################

# Given the variants dictionary stored in variants.txt, apply on given place name
def apply_variant(place, variants):
    if place in variants:
        return variants[place]
    return place


# Turn string into lower case and remove dashes for easier comparison
def break_down(s):
    return s.lower().replace("_", " ").replace("-", " ")


# Given just one name of either a region, country, division or location, search for this within the color_ordering.tsv
# file and return the found hierarchy (region, country, division, location)
def find_place_in_ordering(place, ordering, variants):

    if place == "Unknown":
        return None

    place = apply_variant(place, variants)

    for region in ordering:
        for country in ordering[region]:
            for division in ordering[region][country]:
                for location in ordering[region][country][division]:
                    if break_down(place) == break_down(location):
                        return (region, country, division, location)

    for region in ordering:
        for country in ordering[region]:
            for division in ordering[region][country]:
                if break_down(place) == break_down(division):
                    return (region, country, division, "")

    for region in ordering:
        for country in ordering[region]:
            if break_down(place) == break_down(country):
                # If given name is country, return as division also, leave location empty. (Useful for travel exposure)
                return (region, country, country, "")

    for region in ordering:
        if break_down(place) == break_down(region):
            return (region, region, region, "")


    return None


# Util function to find for a given info the longest matching pattern it can start/end with and return the extracted
# place name(s), or return None
def find_longest_pattern(info, list_of_patterns):
    pattern_found = False
    best_pattern = ("", "", "")
    for pattern in list_of_patterns:
        if "XXX" not in pattern:
            print("Invalid pattern found (missing XXX): " + pattern)
            continue
        pre = pattern.split("XXX")[0]
        post = pattern.split("XXX")[1]

        if info.startswith(pre) and info.endswith(post):
            pattern_found = True
            if len(pattern) > len(best_pattern[0]):
                best_pattern = (pattern, pre, post)
    if pattern_found:
        second_half = info
        if best_pattern[1] != "":
            second_half = info.split(best_pattern[1])[1]
        center = second_half
        if best_pattern[2] != "":
            center = second_half.split(best_pattern[2])[0]

        return center.strip()
    return None



###################################
# Interpret additional info
###################################

# Given an old and a new set of region/country/division/location, produce the necessary annotations for a given seq.
def create_annontation(id, strain, new, old, travel, info, annotations_append, prev = True):
    (region, country, division, location) = new
    (region_original, country_original, division_original, location_original) = old

    info_str = " (" + info + ")"

    t = ""
    if travel:
        t = "_exposure"

    if region_original != region:
        p = ""
        if prev:
            p = " # previously " + region_original + info_str
        annotations_append.append(strain + "\t" + id + "\t" + "region" + t + "\t" + region + p)
        info_str = ""
        print("Adjust region" + t + " " + region_original + " to " + region + " for sequence " + strain)

    if country_original != country:
        p = ""
        if prev:
            p = " # previously " + country_original + info_str
        annotations_append.append(strain + "\t" + id + "\t" + "country" + t + "\t" + country + p)
        info_str = ""
        print("Adjust country" + t + " " + country_original + " to " + country + " for sequence " + strain)

    if division_original != division:
        p = ""
        if prev:
            p = " # previously " + division_original + info_str
        annotations_append.append(strain + "\t" + id + "\t" + "division" + t + "\t" + division + p)
        print("Adjust division" + t + " " + division_original + " to " + division + " for sequence " + strain)
        info_str = ""

    if location != None and location_original != location:
        p = ""
        if prev:
            p = " # previously " + location_original + info_str
        annotations_append.append(strain + "\t" + id + "\t" + "location" + t + "\t" + location + p)
        print("Adjust location" + t + " " + location_original + " to " + location + " for sequence " + strain)
        info_str = ""
    if info_str != "":
        print("No adjustments necessary for sequence " + strain)
    return annotations_append


def interpret_location(place, ordering, variants, pattern_found):
    # Pattern: several hierarchical place names separated by slashes or commas.
    # Remark: Order region/country/division/location is assumed
    for delim in ["/", ","]:
        if delim in place:
            place_list = [elem.strip() for elem in place.split(delim)]

            # Pattern: all 4 hierarchies given
            if len(place_list) == 4:
                print("Interpreted pattern: " + bold("region/country/division/location"))
                region_guess = place_list[0]
                country_guess = place_list[1]
                division_guess = place_list[2]
                location_guess = place_list[3]
                ordering_result = find_place_in_ordering(location_guess, ordering, variants)
                if ordering_result == None:
                    answer = input("Could not correctly interpret " + place + " (unknown location " + location_guess + "). Leave empty to approve of guess " + bold(region_guess + ", " + country_guess + ", " + division_guess + ", " + location_guess) + " or press any key to return to manual processing: ")
                    if answer != "":
                        return None
                    else:
                        return (region_guess, country_guess, division_guess, location_guess)
                print("Location found: " + bold("/".join(ordering_result)))
                if ordering_result[0] == region_guess and ordering_result[1] == country_guess:
                    return ordering_result
                else:
                    print("Found region or country conflicting with given additional info.")
                    return None

            # Pattern: only 3 hierarchies given - determine whether region or location is missing
            elif len(place_list) == 3:
                ordering_result = find_place_in_ordering(place_list[0], ordering, variants)
                if ordering_result == None:
                    print("Interpretation failed. " + bold(place_list[0]) + " is not known")
                    return None
                (region, country, division, location) = ordering_result
                if region == country:
                    print("Interpreted pattern: " + bold("region/country/division"))
                    region_guess = place_list[0]
                    country_guess = place_list[1]
                    division_guess = place_list[2]
                    location_guess = ""

                    ordering_result = find_place_in_ordering(division_guess, ordering, variants)
                    if ordering_result == None:
                        answer = input("Could not correctly interpret " + place + " (unknown division " + division_guess + "). Leave empty to approve of guess " + bold(region_guess + ", " + country_guess + ", " + division_guess) + " or press any key to return to manual processing: ")
                        if answer == "":
                            return (region_guess, country_guess, division_guess, location_guess)
                        else:
                            return None
                    print("Division found: " + bold("/".join(ordering_result)))
                    if ordering_result[0] == region_guess and ordering_result[1] == country_guess:
                        return ordering_result
                    else:
                        print("Found region or country conflicting with given additional info.")
                        return None

                else:
                    print("Interpreted pattern: " + bold("country/division/location"))
                    region_guess = region
                    country_guess = place_list[0]
                    division_guess = place_list[1]
                    location_guess = place_list[2]

                    ordering_result = find_place_in_ordering(location_guess, ordering, variants)
                    if ordering_result == None:
                        answer = input(
                            "Could not correctly interpret " + place + " (unknown location " + location_guess + "). Leave empty to approve of guess " + bold(region_guess + ", " + country_guess + ", " + division_guess + ", " + location_guess) + " or press any key to return to manual processing: ")
                        if answer != "":
                            return None
                        else:
                            return (region_guess, country_guess, division_guess, location_guess)
                    print("Location found: " + bold("/".join(ordering_result)))
                    if ordering_result[0] == region_guess and ordering_result[1] == country_guess:
                        return ordering_result
                    else:
                        print("Found region or country conflicting with given additional info.")
                        return None

            elif len(place_list) == 2:
                ordering_result0 = find_place_in_ordering(place_list[0], ordering, variants)
                ordering_result1 = find_place_in_ordering(place_list[1], ordering, variants)

                if ordering_result0 == None and ordering_result1 == None:
                    answer = input("Unknown places " + bold(place_list[0]) + " and " + bold(place_list[1]) + ". Enter 0 or 1 to select one location to continue, leave empty to return to manual processing: ")
                    if answer == "0":
                        return ("", "", "", place_list[0])
                    elif answer == "1":
                        return ("", "", "", place_list[1])
                    else:
                        return None

                else:
                    if ordering_result0 == None and ordering_result1 != None:
                        return ordering_result1
                    if ordering_result0 != None and ordering_result1 == None:
                        return ordering_result0
                    if ordering_result0 != None and ordering_result1 != None:
                        answer = input("Two locations found: " + bold(place_list[0]) + " and " + bold(place_list[1]) + ". Enter 0 or 1 to select one location to continue, leave empty to return to manual processing: ")
                        if answer == "0":
                            return ordering_result0
                        elif answer == "1":
                            return ordering_result1
                        else:
                            return None

            else:
                print("Interpretation failed. Unknown amount of \"/\" or \",\". Returning to manual processing...\n")
                return None

    # Pattern: no slashes contained in place - assume single place name
    ordering_result = find_place_in_ordering(place, ordering, variants)
    if ordering_result == None:
        if pattern_found:
            answer = input("Unknown location " + bold(place) + ". Press enter to use division provided by strains, press any key to input division manually: ")
            if answer == "":
                return ("", "", "", place)
            else:
                return None
        else:
            print("Unknown location " + bold(place))
            return ("", "", "", place)
    else:
        return ordering_result


def adjust_caps(s):
    if s.upper() == s:
        s = "/".join([" ".join(["'".join([f[0].upper() + f[1:].lower() if len(f) > 1 else f.lower() for f in d.split("'")]) for d in a.split(" ")]) for a in s.split("/")])
    return s

# Check an additional info for patterns fitting additional location info and if possible interpret result
def check_additional_location(info, strain_list, location_pattern, ordering, metadata, annotations_append, variants, auto):
    place = find_longest_pattern(info, location_pattern)

    if place != None:
        return annotations_append, True
    else:
        return annotations_append, False

    if info.endswith(" (interpreted as patient residence)"):
        place = info.split(" (interpreted as patient residence)")[0]

    if place == None:
        # No pattern found - try whether found as location anyway
        pattern_found = False
        place = info
        print("Testing for single location name...")
        ordering_result = find_place_in_ordering(place, ordering, variants)
        if ordering_result == None:
            if auto:
                return annotations_append, True
            print("Location not found. Returning to manual processing...\n")
            return annotations_append, False
    else:
        pattern_found = True
        print("Known " + bold("patient residence") + " pattern found. Extracted location(s): " + bold(place))

    place = "".join(place.split("_h"))
    place = " ".join(place.split("_"))
    place = adjust_caps(place)

    place_interpretation = interpret_location(place, ordering, variants, pattern_found)

    if place_interpretation is None and auto:
        return annotations_append, True
    while place_interpretation == None:
        s = "Could not identify " + bold(place) + ". You have the following options:"
        s += "\n- Enter different name (e.g. in case of typo) for repeated search"
        s += "\n- Enter desired region/country/division/location in this exact format"
        s += "\n- Leave answer empty to return to manual processing"
        s += "\nType your input here: "
        answer = input(s)
        if answer == "":
            print("Return to manual processing...\n")
            return annotations_append, False

        place = answer
        place_interpretation = interpret_location(place, ordering, variants, pattern_found)


    (region, country, division, location) = place_interpretation
    if region == "" and country == "" and division == "":
        if auto:
            return annotations_append, True
        if pattern_found:
            print("Location not found in ordering and no additional geography given. Assume new location " + bold(location))
            regions = []
            countries = []
            divisions = []
            locations = []
            for (id, strain) in strain_list:
                if metadata[id]["region"] not in regions:
                    regions.append(metadata[id]["region"])
                if metadata[id]["country"] not in countries:
                    countries.append(metadata[id]["country"])
                if metadata[id]["division"] not in divisions:
                    divisions.append(metadata[id]["division"])
                if metadata[id]["location"] not in locations:
                    locations.append(metadata[id]["location"])

            if len(regions) == 1 and len(countries) == 1 and len(divisions) == 1:
                answer = input("Suggested pattern (please double check!): " + bold(metadata[id]["region"] + ", " + metadata[id]["country"] + ", " + metadata[id]["division"] + ", " + place) + ". Leave empty to approve, otherwise press any key to return to manual processing: ")
                if answer != "":
                    print("Return to manual processing...\n")
                    return annotations_append, False
                else:
                    region = metadata[id]["region"]
                    country = metadata[id]["country"]
                    division = metadata[id]["division"]
                    location = place
                    for (id, strain) in strain_list:
                        new = (region, country, division, location)
                        old = (metadata[id]["region"], metadata[id]["country"], metadata[id]["division"], metadata[id]["location"])
                        create_annontation(id, strain, new, old, False, info, annotations_append)
                    return annotations_append, True

            else:
                print("Could not correctly interpret " + info + " (unknown description " + place + "). No suggestions possible (several countries/divisions with this info). Returning to manual processing...\n")
                return annotations_append, False

    else:
        answer = "" if auto else input("Interpreted " + bold(info) + " as " + bold(region + "/" + country + "/" + division + "/" + location) + ". Press enter to approve, press any other key to return to manual processing: ")
        if answer != "":
            print("Return to manual processing...\n")
            return annotations_append, False
        for (id, strain) in strain_list:
            new = (region, country, division, location)
            old = (metadata[id]["region"], metadata[id]["country"], metadata[id]["division"], metadata[id]["location"])
            create_annontation(id, strain, new, old, False, info, annotations_append)
        return annotations_append, True


    return annotations_append, False


# Check the given info for matching patterns known to carry travel history information, and apply the necessary
# annotations
def check_travel_history(info, strain_list, travel_pattern, ordering, metadata, annotations_append, variants, auto):

    place = find_longest_pattern(info, travel_pattern)

    if info.endswith(" (interpreted as travel exposure)"):
        place = info.split(" (interpreted as travel exposure)")[0]

    if place != None:

        print("Known " + bold("travel exposure") + " pattern found. Extracted location(s): " + bold(place))


        if "," in place:
            places = [p.strip() for p in place.split(",")]
            print("Several location names found. Process each individually: " + bold("[" + ", ".join(places) + "]"))
        elif " and " in place:
            places = [p.strip() for p in place.split(" and ")]
            print("Several location names found. Process each individually: " + bold("[" + ", ".join(places) + "]"))
        else:
            places = [place]

        results = {"region_exposure": [], "country_exposure": [], "division_exposure": []}
        # In case several travels are listed, check all of them and find overlaps
        for place in places:
            if "/" in place:
                place = place.split("/")[-1].strip()
            ordering_result = find_place_in_ordering(place, ordering, variants)
            if ordering_result is None and auto:
                return annotations_append, True
            while ordering_result == None:
                s = "Could not identify " + bold(place) + ". You have the following options:"
                s += "\n- Enter different name (e.g. in case of typo) for repeated search"
                s += "\n- Enter desired region/country/division in this exact format"
                s += "\n- Leave answer empty to return to manual processing"
                s += "\nType your input here: "
                answer = input(s)
                if answer == "":
                    print("Return to manual processing...\n")
                    return annotations_append, False
                if answer.count("/") == 2:
                    ordering_result = (answer.split("/")[0], answer.split("/")[1], answer.split("/")[2], "")
                else:
                    place = answer
                    ordering_result = find_place_in_ordering(place, ordering, variants)

            (region_exposure, country_exposure, division_exposure, location_exposure) = ordering_result
            if region_exposure not in results["region_exposure"]:
                results["region_exposure"].append(region_exposure)
            if country_exposure not in results["country_exposure"]:
                results["country_exposure"].append(country_exposure)
            if division_exposure not in results["division_exposure"]:
                results["division_exposure"].append(division_exposure)

        if len(results["region_exposure"]) > 1:
            answer = "" if auto else input("Contains exposures from different regions. No annotation possible. Press " + bold("ENTER") + " to skip this info, otherwise press any key to return to manual processing: ")
            if answer == "":
                return annotations_append, True
            else:
                print("Return to manual processing...\n")
                return annotations_append, False
        final_region = results["region_exposure"][0]

        if len(results["country_exposure"]) > 1:
            final_country = final_region
        else:
            final_country = results["country_exposure"][0]

        if len(results["division_exposure"]) > 1:
            final_division = final_country
        else:
            final_division = results["division_exposure"][0]

        answer = "" if auto else input("Interpreted as " + bold(final_region + ", " + final_country + ", " + final_division) + ". Press " + bold("ENTER") + " to approve, otherwise press any key to return to manual processing: ")
        if answer != "":
            print("Return to manual processing...\n")
            return annotations_append, False

        for (id, strain) in strain_list:
            new = (final_region, final_country, final_division, None)
            old = (metadata[id]["region"], metadata[id]["country"], metadata[id]["division"], metadata[id]["location"])
            create_annontation(id, strain, new, old, True, info, annotations_append)

        return annotations_append, True

    return annotations_append, False





# Iterate over every info given and search for known patterns, or provide interactive interface for manual processing
def check_additional_info(additional_info, path_to_config_files, auto):
    metadata = read_metadata(path_to_nextstrain + "ncov/data/downloaded_gisaid.tsv", additional_info)
    if metadata == None:
        return []
    ordering = read_ordering_file(path_to_nextstrain + "ncov/defaults/color_ordering.tsv")
    variants = read_dict(path_to_config_files + "variants.txt")

    sorted_info = rearrange_additional_info(additional_info)

    # Collected patterns
    info_ignore = read_simple_file(path_to_config_files + "info_ignore.txt")
    location_pattern = read_simple_file(path_to_config_files + "location_pattern.txt")
    travel_pattern = read_simple_file(path_to_config_files + "travel_pattern.txt")
    purpose_of_sequencing = read_dict(path_to_config_files + "purpose_of_sequencing.txt")

    print("\n##########################################\n")

    annotations_append = []

    for key in sorted_info:
        (info, info_type) = key
        strain_list = sorted_info[key]
        info_found = False

        if info.startswith("Other: "):
            info = info[7:]
            print("Treating \"Other: " + info + "\" as \"" + info + "\"")

        if info.startswith("other: "):
            info = info[7:]
            print("Treating \"other: " + info + "\" as \"" + info + "\"")

        print("Processing " + bold(info) + ":")

        while True:

            auto_comments = ["travel surveillance"]

            # Special case:
            if (info.startswith("Resident of ") or info.startswith("resident of ")) and " tested in " in info:
                if info.startswith("Resident of "):
                    division = ((info.split(" tested in ")[0].strip(",")).strip(";")).split("Resident of ")[1].strip()
                if info.startswith("resident of "):
                    division = ((info.split(" tested in ")[0].strip(",")).strip(";")).split("resident of ")[1].strip()
                division_exposure = info.split(" tested in ")[1].strip()
                ordering_result = find_place_in_ordering(division, ordering, variants)
                if ordering_result != None:
                    (region, country, division, location) = ordering_result
                    ordering_result = find_place_in_ordering(division_exposure, ordering, variants)
                    if ordering_result != None:
                        (region_exposure, country_exposure, division_exposure, location) = ordering_result

                        answer = "" if auto else input("Interpret " + bold(info) + " as special case with division " + bold(division) + " and exposure " + bold(division_exposure) + ". Leave empty to approve, or press any key to continue with manual processing: ")
                        if answer == "":
                            for (id, strain) in strain_list:
                                new_residence = (region, country, division, None)
                                new_exposure = (region_exposure, country_exposure, division_exposure, None)
                                old = (metadata[id]["region"], metadata[id]["country"], metadata[id]["division"],metadata[id]["location"])
                                create_annontation(id, strain, new_residence, old, False, info, annotations_append) # Change residence
                                create_annontation(id, strain, new_exposure, new_residence, True, info, annotations_append, prev = False) # Change exposure right back
                            info_found = True
            if info_found:
                break

            if info in purpose_of_sequencing:
                answer = "" if auto else input("Identified as \"purpose_of_sequencing\". Press " + bold("ENTER") + " to approve, otherwise press any key: ")
                if answer == "":
                    for (id, strain) in strain_list:
                        annotations_append.append(strain + "\t" + id + "\t" + "sampling_strategy" + "\t" + purpose_of_sequencing[info] + " # " + info)
                        print("purpose_of_sequencing info added as annotation for strain " + id)
                    break

            if info in info_ignore:
                answer = "" if auto else input("Identified as \"Ignore\". Press " + bold("ENTER") + " to approve, otherwise press any key: ")
                if answer == "":
                    break

            annotations_append, info_found = check_travel_history(info, strain_list, travel_pattern, ordering, metadata,
                                                                  annotations_append, variants, auto)
            if info_found:
                break

            annotations_append, info_found = check_additional_location(info, strain_list, location_pattern, ordering, metadata, annotations_append, variants, auto)
            if info_found or auto:
                print("Location pattern found. Skipping entry...")
                break

            s = bold(info) + " did not contain known pattern or could not be interpreted. You have the following options:"
            s += "\n" + bold("l") + " - force interpretation as " + bold("patient residence")
            s += "\n" + bold("t") + " - force interpretation as " + bold("travel exposure")
            s += "\n" + bold("i") + " - add info to " + bold("ignore")
            s += "\n" + bold("a") + " - add to annotations as a " + bold("comment")
            s += "\n" + bold("nl") + " - add new " + bold("patient residence") + " pattern"
            s += "\n" + bold("nt") + " - add new " + bold("travel exposure") + " pattern"
            s += "\n" + bold("ns") + " - add new " + bold("sequencing purpose") + " pattern"
            s += "\n" + bold("s") + " - " + bold("skip") + " this additional info"
            s += "\n" + bold("r") + " - " + bold("repeat") + " as before"
            s += "\nType input here: "


            answer = input(s)
            if answer == "a":
                for (id, strain) in strain_list:
                    annotations_append.append("# " + strain + "\t" + id + "\t" + info_type + ": " + info)
                    print("Add comment for " + id)
                break
            if answer == "i":
                add_to_simple_file(path_to_config_files + "info_ignore.txt", info)
                break
            elif answer == "l":
                print("Process " + bold(info) + " now as " + bold(info + " (interpreted as patient residence)"))
                info = info + " (interpreted as patient residence)"
            elif answer == "t":
                print("Process " + bold(info) + " now as " + bold(info + " (interpreted as travel exposure)"))
                info = info + " (interpreted as travel exposure)"
            elif answer == "nl":
                pattern = input("Type pattern here (don't forget XXX as placeholder): ")
                add_to_simple_file(path_to_config_files + "location_pattern.txt", pattern)
                location_pattern.append(pattern)
            elif answer == "nt":
                pattern = input("Type pattern here (don't forget XXX as placeholder): ")
                add_to_simple_file(path_to_config_files + "travel_pattern.txt", pattern)
                travel_pattern.append(pattern)
            elif answer == "ns":
                pattern1 = input("Type additional info pattern here: ")
                pattern2 = input("Type desired metadata entry here: ")
                add_to_simple_file(path_to_config_files + "purpose_of_sequencing.txt", pattern1 + "\t" + pattern2)
                purpose_of_sequencing[pattern1] = pattern2
            elif answer == "s":
                break
            print("\n")


        print("\n-------\n")

    with open(path_to_outputs + "omicron_additional_info.txt", "w") as out:
        for key in sorted_info:
            (info, info_type) = key
            strain_list = sorted_info[key]
            for (id, strain) in strain_list:
                if metadata[id]["Nextstrain_clade"] == "21K (Omicron)":
                    out.write(id + ": " + info + "\n")

    return annotations_append


path_to_input = "scripts/curate_metadata/inputs_new_sequences/"
path_to_config_files = "scripts/curate_metadata/config_files_additional_info/"
path_to_nextstrain = "../"
path_to_outputs = "scripts/curate_metadata/outputs_new_sequences/"

Path(path_to_outputs).mkdir(parents=True, exist_ok=True)
Path(path_to_input).mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Parse a file containing additional info from GISAID to create annotations from recognized patterns.",
    )
    parser.add_argument("--auto", action="store_true", help="Automatically parse recognized patterns as opposed to asking for user verification.")
    args = parser.parse_args()


    additional_info = read_data(path_to_input)
    if additional_info != {}:

        annotations_append = check_additional_info(additional_info, path_to_config_files, args.auto)

        # Print out necessary annotations
        annotations = read_simple_file(path_to_nextstrain + "ncov-ingest/source-data/gisaid_annotations.tsv")
        annot_lines_to_write = []
        for line in annotations_append:
            if line in annotations:
                continue
            #print(line)
            if "=" not in line:
                annot_lines_to_write.append(line)
            if len(line.split("\t")) == 4:
                number_of_occurences = annotations.count(line.split("\t")[1])
                irrelevant_occurences = sum([(line.split("\t")[1] + "\t" + s) in annotations for s in ["title", "authors", "paper_url"]])
                if number_of_occurences > irrelevant_occurences:
                    print("Warning: " + line.split("\t")[1] + " already exists in annotations!")

        with open(path_to_outputs + "additional_info_annotations.tsv", 'w') as out:
            out.write("\n".join(sorted(annot_lines_to_write)))

    else:
        print("No additional_info files provided.")
