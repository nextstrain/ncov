import os
from difflib import SequenceMatcher
from pathlib import Path
from geopy.geocoders import Nominatim
from collections import Counter

# Used for printing bold in console output
def bold(s):
    return('\033[1m' + s + '\033[0m')

replace_special_char = {
    "é":"e",
    "è":"e",
    "ü":"ue",
    "ä":"ae",
    "ö":"oe",
    "í":"i",
    "ó":"o",
    "ç":"c",
    "á":"a",
    "'":" ",
    "â":"a",
    "š":"s",
    "ť":"t",
    "ñ":"n",
    "ř":"r",
    "ž":"z",
    "ů":"u",
    "ý":"y",
    "ě":"e",
    "ň":"n",
    "ã":"a",
    "ê":"e",
    "č":"c",
    "ô":"o",
    "ı":"i",
    "ú":"u",
    "ś":"s",
    "ą":"q",
    "à":"a",
    "å":"a",
    "ł":"l",
    "-":" ",
    "î":"i",
    "ŕ":"r",
    "ľ":"l",
    "ď":"d",
    "ć": "c",
    "ș": "s",
}

# Compare strings without considering special characters and caps
def clean_string(s):
    s = s.lower()
    for c in replace_special_char:
        s = s.replace(c, replace_special_char[c])
    return s

# Read in files containing simple lists, dictionaries or geoLocationRules
def read_local_file(file_name):
    path_file_name = path_to_config_files + file_name

    with open(path_file_name) as myfile:
        file_content = myfile.readlines()

    list_format = [accepted_additions_file]

    if file_name in list_format:  # simple list
        return [line.strip() for line in file_content[1:]]

    abbreviations_format = [abbreviations_file]

    if file_name in abbreviations_format:  # dictionary, keys seperated from content with tabs
        content = {}
        for line in file_content:
            if line == "\n":
                continue
            if line.startswith("#"):
                key = line.strip().split("# ")[1]
                content[key] = {}
                continue
            l = line.strip().split("\t")
            if l[0] in content[key]:
                print("Attention, duplicate found while reading " + file_name + ": " + l[0] + " -> " + l[1] + ", " + content[key][l[0]])
            content[key][l[0]] = l[1]
        return content

    geoLocationRule_format = [geoLocationRules_file, manualAnnotationRules_file, internationalExceptions_file]

    if file_name in geoLocationRule_format: # Read as simple dictionary
        content = {}
        for line in file_content:
            if line == "\n":
                continue
            l = line.strip().split("\t")
            k = l[0]
            c = l[1]
            if k in content:
                print("Attention, duplicate found while reading " + file_name + ": " + k + " -> " + c + ", " + content[k])
            content[k] = c
        return content

# Read ordering file into simple dictionary based on type (e.g. location, division, pango_lineage)
def read_ordering(path):
    with open(path + ordering_file) as myfile:
        file_content = myfile.readlines()

    ordering = {}
    for line in file_content:
        if line == "\n" or line.startswith("#"):
            continue
        l = line.strip().split("\t")
        type = l[0]
        name = l[1]
        if type not in ordering:
            ordering[type] = []
        ordering[type].append(name)

    return ordering

# Read lat_longs into dictionary based on type (e.g. location, division)
def read_latlongs(path):
    with open(path + latlongs_file) as myfile:
        file_content = myfile.readlines()

    latlongs = {"location": {}, "division": {}, "country": {}, "region": {}}

    for line in file_content:
        if line == "\n" or line.startswith("#"):
            continue
        l = line.strip().split("\t")
        type = l[0]  # location, division etc.
        name = l[1]

        if name not in latlongs[type]:
            latlongs[type][name] = (float(l[2]), float(l[3])) # Store as float to enable sorting by lat_longs
        else:
            print("Duplicate in lat_longs? (" + l[0] + " " + l[1] + ")\n")

    return latlongs

# Read metadata into a multi-level dictionary accessible via data[region][country][division] = list_of_locations
def read_metadata(metadata_filename, data, geo_location_occurences, genbank=False):

    with open(path_to_metadata + metadata_filename) as f:
        header = f.readline().split("\t")
        country_i = header.index("country")
        region_i = header.index("region")
        division_i = header.index("division")
        location_i = header.index("location")

        line = f.readline()
        while line:
            l = line.split("\t")
            country = l[country_i]
            region = l[region_i]
            division = l[division_i]
            location = l[location_i]

            # automatically increment genbank locations to the threshold since we
            # don't want to skip any for now.
            increment = 1 if not genbank else 20

            geo_location_occurences["region"].update({region: increment})
            geo_location_occurences["country"].update({country: increment})
            geo_location_occurences["division"].update({division: increment})
            geo_location_occurences["location"].update({location: increment})

            if region not in data:
                data[region] = {}
            if country not in data[region]:
                data[region][country] = {}
            if division not in data[region][country]:
                data[region][country][division] = []
            if location not in data[region][country][division]:
                data[region][country][division].append(location)

            line = f.readline()

    return data, geo_location_occurences

# Read metadata again, but this time consider region_exposure, country_exposure and division_exposure (e.g. travel info)
# In order to be properly displayed, exposure geographies need also be included in color_ordering and lat_longs
# Only include exposure geography if already specified in accepted exposure file, otherwise print warning
def read_exposure(data, metadata_filename, accepted_additions_file):
    # divisions and countries that are accepted additions to the metadata
    accepted_exposure = read_local_file(accepted_additions_file)

    # Check given accepted exposures and print warning if already included in the data
    # (e.g. if the country was unknown when the exposure was registered, but had sequences added in the meantime)
    for region in data:
        for country in data[region]:
            if country + " (" + region + ")" in accepted_exposure:
                print("Specified exposure " + bold(country + " (" + region + ")") + " is no longer needed and can be removed from " + accepted_additions_file + ".")
            for division in data[region]:
                if division + " (" + country + ", " + region + ")" in accepted_exposure:
                    print("Specified exposure " + bold(division + " (" + country + ", " + region + ")") + " is no longer needed and can be removed from " + accepted_additions_file + ".")

    with open(path_to_metadata + metadata_filename) as f:
        header = f.readline().split("\t")
        region_i = header.index("region_exposure")
        country_i = header.index("country_exposure")
        division_i = header.index("division_exposure")
        epi_i = header.index("gisaid_epi_isl")

        line = f.readline()
        while line:
            l = line.split("\t")
            region = l[region_i]
            country = l[country_i]
            division = l[division_i]
            epi = l[epi_i]

            if region not in data:
                print("Strain " + epi + " has unknown region_exposure " + bold(region) + ". Please correct!")
            else:
                s1 = country + " (" + region + ")"
                country_present = True
                if country not in data[region]:
                    country_present = False
                    if s1 in accepted_exposure or country == region:
                        data[region][country] = {}
                        country_present = True
                    else:
                        print("Strain " + epi + " has unknown country_exposure " + bold(country) + ". Please correct or consider adding " + bold(s1) + " to " + accepted_additions_file + "!")
                if country_present:
                    s2 = division + " (" + country + ", " + region + ")"
                    if division not in data[region][country]:
                        if s2 in accepted_exposure or division == country:
                            data[region][country][division] = [""]
                        else:
                            print("Strain " + epi + " has unknown division_exposure " + bold(division) + ". Please correct or consider adding " + bold(s2) + " to " + accepted_additions_file + "!")
            line = f.readline()
    return data

# Given a set of corrections, adjust the data dictionary accordingly
def correct_data(data, corrections):
    for (region, country, division, location, region_correct, country_correct, division_correct, location_correct) in corrections:
        if country_correct not in data[region_correct]:
            data[region_correct][country_correct] = {}
        if division_correct not in data[region_correct][country_correct]:
            data[region_correct][country_correct][division_correct] = []
        if location_correct not in data[region_correct][country_correct][division_correct]:
            data[region_correct][country_correct][division_correct].append(location_correct)

        # If no entries are contained downstream, assume geography level is now obsolete and delete from the data dictionary
        # (e.g. if a misspelled region is corrected, delete it from data after copying over all downstream countries, divisions and locations)
        data[region][country][division].remove(location)
        if data[region][country][division] == []:
            del data[region][country][division]
        if data[region][country] == {}:
            del data[region][country]
        if data[region] == {}:
            del data[region]
    return data

# Given 3 sets of (region, country, division, location), decide whether the 'given' set corresponds to the
# 'before' set (considering the wildcard *) and formulate the correction needed to change it to the to 'correct' set
def formulate_correction(given, before, correct):
    (region, country, division, location) = given
    (region_before, country_before, division_before, location_before) = before
    (region_correct, country_correct, division_correct, location_correct) = correct
    if region_correct == "*":
        region2 = region
    else:
        region2 = region_correct

    if country_correct == "*":
        country2 = country
    else:
        country2 = country_correct

    if division_correct == "*":
        division2 = division
    else:
        division2 = division_correct

    if location_correct == "*":
        location2 = location
    else:
        location2 = location_correct

    if (region == region_before or region_before == "*") \
            and (country == country_before or country_before == "*") \
            and (division == division_before or division_before == "*") \
            and (location == location_before or location_before == "*"):
        return (region, country, division, location, region2, country2, division2, location2)
    return None

# Given a set of rules, traverse the data dictionary and find and correct cases where rules apply
# For manualAnnotationRules: Enable different delimiter than '/' in case this is included in a location name
def apply_rules(data, ruleSet, delimiter = ["/"], print_rules = True):
    rules = read_local_file(ruleSet)

    applied_rules = {}
    for g in rules:
        for d in delimiter:
            rules_apply = []
            if d not in g or d not in rules[g]:
                continue
            (region_before, country_before, division_before, location_before) = g.split(d)
            (region_correct, country_correct, division_correct, location_correct) = rules[g].split(d)

            # Due to reoccuring bug: Since empty divisions are automatically filled with the country name later in the
            # ncov-ingest pipeline, give a warning when detecting a rule that might be affected
            if country_before == division_before:
                recommended_rule = d.join([region_before, country_before, "", ""])
                if recommended_rule not in rules:
                    print(bold("Attention: Consider automatic division filler applied after geoLocationRules (Hint: add [" + recommended_rule + "\t" + rules[g] + "])"))

            for region in data:
                for country in data[region]:
                    for division in data[region][country]:
                        for location in data[region][country][division]:
                            correction = formulate_correction((region, country, division, location), (region_before, country_before, division_before, location_before), (region_correct, country_correct, division_correct, location_correct))
                            if correction is not None:
                                rules_apply.append(correction)
                                if print_rules:
                                    print("/".join(correction[:4]) + "\t" + "/".join(correction[4:]))
                                applied_rules[correction[:4]] = correction[4:]

            data = correct_data(data, rules_apply)

    # Also return all rules that were applied to the data to make detection of conflicting annotations easier later
    return data, applied_rules

# Check all locations if the appear as division elsewhere (ignore cases where division == location)
def check_division_inconsistency(data):
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                if division != "":
                    for location in data[region][country][division]:
                        if location != "":
                            if location in data[region][country] and location != division:
                                print(bold(location) + " found as both division and location within division " + bold(division) + ".")
                                if list(data[region][country][location]) != [""]: # Locations found below both divisions - needs manual review which division is proper
                                    print("Conflict found: Both divisions contain locations:")
                                    l = data[region][country][location]
                                    if len(l) > 10:
                                        s = ", ".join(l[:10]) + "... (plus " + str(len(l) - 10) + " more)"
                                    else:
                                        s = ", ".join(l)
                                    print("division " + bold(location) + ": location(s) " + s)

                                    l = data[region][country][division]
                                    if len(l) > 10:
                                        s = ", ".join(l[:10]) + "... (plus " + str(len(l) - 10) + " more)"
                                    else:
                                        s = ", ".join(l)
                                    print("division " + bold(division) + ": location(s) " + s)
                                    print("(Template for correction" + "[" + "/".join([region, country, location, "?"]) + "\t" + "/".join([region, country, division, "?"]) + "])\n")

                                else: # No location found below the affected location/division - change to proper level
                                    print("/".join([region, country, location, ""]) + "\t" + "/".join([region, country, division, location]) + "\n")

# Search for duplicates on location and division level
def check_duplicates(data, abbreviations_file):
    abbreviations = read_local_file(abbreviations_file)

    # Collect all locations and their (region, country, division) origin
    # Cut away all present duplicate specifiers (e.g. 'Guadalupe ES' -> 'Guadalupe') - data is treated as if no
    # duplicate adjustment has happened before. This way, changes in duplicates can be detected and properly treated
    # (e.g. if another duplicate appears in the same country but different division, adjust the duplicate specifier
    # from country abbreviation to division) (e.g. 'Guadalupe ES' -> 'Guadalupe (Extremadura)')
    location_origin = {}
    for region in data:
        if region not in region_order:
            continue
        for country in data[region]:
            if country == region:
                continue
            if country not in abbreviations["country"]:
                print("Abbreviation missing for " + country + ". Please add to " + abbreviations_file)
                continue
            for division in data[region][country]:
                if country == "USA" and division not in abbreviations["division"]:
                    print("Abbreviation missing for US state " + division + ". Please add to " + abbreviations_file)
                    continue
                for location in data[region][country][division]:
                    if location == "":
                        continue
                    if location.endswith(" (" + division + ")"): # Cut away already existing duplicate specifiers
                        location = location.split(" (" + division + ")")[0]
                    elif location.endswith(" " + abbreviations["country"][country]):
                        location = location.split(" " + abbreviations["country"][country])[0]
                    elif country == "USA" and location.endswith(" " + abbreviations["division"][division]):
                        location = location.split(" " + abbreviations["division"][division])[0]
                    elif len(location.split(" ")[-1]) > 1 and location.upper().split(" ")[-1] == location.split(" ")[-1] or "(" in location:
                        # If parentheses or caps found at the end of the location, consider potential invalid duplicate specifier
                        if not location.split(" ")[-1].isnumeric():
                            #print(f"{'/'.join([region, country, division, location])}\t{'/'.join([region, country, division, location.split(' ')[0]])}")
                            print("Potential duplicate inconsistent with current rules: " + location)

                    if location not in location_origin:
                        location_origin[location] = []
                    location_origin[location].append((region, country, division))
    print()

    # Filter for duplicates (locations that have more than one (region, country, division) origin set)
    locations_duplicates = {}
    for location in location_origin:
        if len(location_origin[location]) > 1: # more than one (region, country, division) origin
            reduced = []
            countries = []
            for combination in location_origin[location]:
                if combination not in reduced:
                    countries.append(combination[1])
                    reduced.append(combination)
            # If, after reducing, only one set of (region, country, division) is left, that means there was a location
            # unnecessarily specified as duplicate where the location doesn't exist in any other country/division
            # In that case print out warning that this location needs not be considered as a duplicate anymore and the
            # corresponding geoLocationRules and annotations should be corrected
            if len(reduced) == 1 and countries != ["USA"]:
                print("Unnecessary duplicate: " + bold(location) + "\n")
            else:
                # Unless country is USA, then leave these "unneccessary" duplicate specifications, as the automatic
                # county assignment already considers duplicates
                locations_duplicates[location] = reduced

    # Apply duplicates rules
    for location in locations_duplicates:
        printed_message = []
        divisions = {}
        for (region, country, division) in locations_duplicates[location]:
            if country not in divisions:
                divisions[country] = []
            divisions[country].append(division)

        for (region, country, division) in locations_duplicates[location]:
            if country == "USA": # For locations in the USA: always use state abbreviation
                location_new = location + " " + abbreviations["division"][division]
                if location in data[region][country][division]:
                    printed_message.append("/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, location_new]))
                if location + " (" + division + ")" in data[region][country][division]:
                    #printed_message.append("Please update duplicate " + bold(location + " (" + division + ")") + " to " + bold(location_new) + " for consistency.")
                    printed_message.append("/".join([region, country, division, location + " (" + division + ")"]) + "\t" + "/".join([region, country, division, location_new]))

            elif len(divisions[country]) == 1: # Among-country duplicate - use country abbreviation
                location_new = location + " " + abbreviations["country"][country]
                if location in data[region][country][division]:
                    printed_message.append("/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, location_new]))
                if location + " (" + division + ")" in data[region][country][division]:
                    #printed_message.append("Please update duplicate " + bold(location + " (" + division + ")") + " to " + bold(location_new) + " for consistency.")
                    printed_message.append("/".join([region, country, division, location + " (" + division + ")"]) + "\t" + "/".join([region, country, division, location_new]))

            else: # Within-country duplicate - use division as unique identifier
                location_new = location + " (" + division + ")"
                if location in data[region][country][division]:
                    printed_message.append("/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, location_new]))
                if location + " " + abbreviations["country"][country] in data[region][country][division]:
                    #printed_message.append("Please update duplicate " + bold(location + " " + abbreviations["country"][country]) + " to " + bold(location_new) + " for consistency.")
                    printed_message.append("/".join([region, country, division, location + " " + abbreviations["country"][country]]) + "\t" + "/".join([region, country, division, location_new]))
        if printed_message != []:
            #print("Duplicate found: " + bold(location))
            for l in printed_message:
                print(l)
            #print()

    ### DIVISION ###
    print("\n----------\n")
    print("Checking for division duplicates...\n")

    # Collect all divisions and their (region, country) origin
    division_origin = {}
    for region in data:
        if region not in region_order:
            continue
        for country in data[region]:
            if country == region or country not in abbreviations["country"]:
                continue
            for division in data[region][country]:
                if division.endswith(" " + abbreviations["country"][country]):
                    division = division.split(" " + abbreviations["country"][country])[0]
                if division not in division_origin:
                    division_origin[division] = []
                # Special cases where a duplicate specification seems out of place for one of the countries (e.g. US states)
                if (division == "Montana" or division == "Maryland") and country == "USA":
                    print("(Ignoring duplicate division " + division + " in favor of the USA.)")
                elif country == "Luxembourg" and division == "Luxembourg":
                    print("(Ignoring duplicate division " + division + " in favor of the country Luxembourg.)")
                else:
                    division_origin[division].append((region, country))

    # Filter for duplicates
    division_duplicates = {}
    for division in division_origin:
        if len(division_origin[division]) > 1:
            reduced = []
            countries = []
            for combination in division_origin[division]:
                if combination not in reduced:
                    countries.append(combination[1])
                    reduced.append(combination)
            if len(reduced) == 1:
                print("Unnecessary duplicate: " + bold(division) + "\n")
            else:
                division_duplicates[division] = reduced

    print()
    # Apply duplicates rules
    for division in division_duplicates:
        printed_message = []

        for (region, country) in division_duplicates[division]:
            division_new = division + " " + abbreviations["country"][country]
            if division in data[region][country]:
                printed_message.append("/".join([region, country, division, "*"]) + "\t" + "/".join([region, country, division_new, "*"]))

        if printed_message != []:
            print("Duplicate found: " + bold(division))
            for l in printed_message:
                print(l)
            print()

    ### COUNTRY ###
    print("\n----------\n")
    print("Checking for country duplicates...\n")

    country_origin = {}
    for region in data:
        if region not in region_order:
            continue
        for country in data[region]:
            if country not in country_origin:
                country_origin[country] = []
            country_origin[country].append(region)

    for country in country_origin:
        if len(country_origin[country]) > 1:
            print("Duplicate country found: " + bold(country) + " within " + bold(", ".join(country_origin[country])))
            if len(country_origin[country]) == 2:
                print("/".join([country_origin[country][0], country, "*", "*"]) + " <-> " + "/".join([country_origin[country][1], country, "*", "*"]))

    return locations_duplicates

# Search for all locations, divisions, countries and regions that are not present in the lat_longs.tsv file
def missing_coordinates(data, path, geo_location_occurences):
    missing_latlongs = {"region": [], "country": {}, "division": {}, "location": {}}

    latlongs = read_latlongs(path)

    for region in data:
        if region not in region_order:
            missing_latlongs["region"].append(region)

        for country in data[region]:
            if country not in latlongs["country"]:
                if region not in missing_latlongs["country"]:
                    missing_latlongs["country"][region] = []
                missing_latlongs["country"][region].append(country)

            division_threshold_function = lambda division: division not in latlongs["division"] and (geo_location_occurences["division"][division] >= 20)
            for division in filter(division_threshold_function, data[region][country]):
                if region not in missing_latlongs["division"]:
                    missing_latlongs["division"][region] = {}
                if country not in missing_latlongs["division"][region]:
                    missing_latlongs["division"][region][country] = []
                missing_latlongs["division"][region][country].append(division)

    return missing_latlongs

# Print out missing locations, divisions etc. in a sorted manner
def print_missing_places(missing_latlongs):
    ### DIVISION ###
    print("\n----------\n")
    if missing_latlongs['division']:
        print("Missing divisions:")
        for region in missing_latlongs["division"]:
            print("# " + region + " #")
            for country in missing_latlongs["division"][region]:
                print(country)
                for division in missing_latlongs["division"][region][country]:
                    print("\tdivision\t" + bold(division))
            print()
    else:
        print("No missing divisions")

    ### COUNTRY ###
    print("\n----------\n")
    if missing_latlongs['country']:
        print("\nMissing countries:")
        for region in missing_latlongs["country"]:
            print("# " + region + " #")
            for country in missing_latlongs["country"][region]:
                print("\tcountry\t" + bold(country))
    else:
        print("No missing countries")

    ### REGION ###
    if missing_latlongs['region']:
        print("\n----------\n")
        print("\nMissing regions:")
        for region in missing_latlongs["region"]:
            print("\tregion\t" + bold(region))


# For all missing place names, search the data for similarly spelled names
# Always search only within the same ordering level as the missing name, as well as below the same higher-level place
# (e.g. for a missing location, compare only to locations within the same division)
def search_similar_names(data, missing_latlongs, locations_duplicates):

    abbreviations = read_local_file(abbreviations_file)

    ### DIVISION ###
    print("\n----------\n")
    identical = []
    similar = {}

    for region in missing_latlongs["division"]:
        for country in missing_latlongs["division"][region]:
            for division in missing_latlongs["division"][region][country]:
                similarity_score = 0
                identical_hit = False
                best_match = None

                for division2 in data[region][country]:
                    if division2 == division:
                        continue
                    if clean_string(division) == clean_string(division2): # Identical except for alternative chars
                        identical.append("/".join([region, country, bold(division), "*"]) + "\t" + "/".join([region, country, bold(division2), "*"]))
                        identical_hit = True
                        break

                    diff = SequenceMatcher(None, division, division2).ratio() # Similarity score if not perfect hit
                    if diff > 0.6:
                        if diff > similarity_score:
                            similarity_score = diff
                            best_match = division2

                if not identical_hit and best_match is not None:
                    while similarity_score in similar:
                        similarity_score += 0.000000000000001
                    similar[similarity_score] = "/".join([region, country, bold(division), "*"]) + "\t" + "/".join([region, country, bold(best_match), "*"])

    if identical:
        print("Identical divisions:")
        for l in identical:
            print(l)

    if similar:
        print("\nSimilar divisions (sorted by descending similarity):")
        for l in sorted(similar, reverse=True):
            print(similar[l])

    ### COUNTRY ###
    print("\n----------\n")
    identical = []
    similar = {}

    for region in missing_latlongs["country"]:
        for country in missing_latlongs["country"][region]:
            similarity_score = 0
            identical_hit = False
            best_match = None

            for country2 in data[region]:
                if country2 == country:
                    continue
                if clean_string(country) == clean_string(country2):  # Identical except for alternative chars
                    identical.append("/".join([region, bold(country), "*", "*"]) + "\t" + "/".join([region, bold(country2), "*", "*"]))
                    identical_hit = True
                    break

                diff = SequenceMatcher(None, country, country2).ratio()  # Similarity score if not perfect hit
                if diff > 0.6:
                    if diff > similarity_score:
                        similarity_score = diff
                        best_match = country2

            if not identical_hit and best_match is not None:
                while similarity_score in similar:
                    similarity_score += 0.000000000000001
                similar[similarity_score] = "/".join([region, bold(country), "*", "*"]) + "\t" + "/".join([region, bold(best_match), "*", "*"])

    if identical:
        print("Identical countries:")
        for l in identical:
            print(l)

    if similar:
        print("\nSimilar countries (sorted by descending similarity):")
        for l in sorted(similar, reverse=True):
            print(similar[l])


    ### REGION ###
    print("\n----------\n")
    identical = []
    similar = {}

    for region in missing_latlongs["region"]:
        similarity_score = 0
        identical_hit = False
        best_match = None

        for region2 in data:
            if region2 == region:
                continue
            if clean_string(region) == clean_string(region2):  # Identical except for alternative chars
                identical.append("/".join([bold(region), "*", "*", "*"]) + "\t" + "/".join([bold(region2), "*", "*", "*"]))
                identical_hit = True
                break

            diff = SequenceMatcher(None, region, region2).ratio()  # Similarity score if not perfect hit
            if diff > 0.6:
                if diff > similarity_score:
                    similarity_score = diff
                    best_match = region2

        if not identical_hit and best_match is not None:
            while similarity_score in similar:
                similarity_score += 0.000000000000001
            similar[similarity_score] = "/".join([bold(region), "*", "*", "*"]) + "\t" + "/".join([bold(best_match), "*", "*", "*"])

    if identical:
        print("Identical regions:")
        for l in identical:
            print(l)

    if similar:
        print("\nSimilar regions (sorted by descending similarity):")
        for l in sorted(similar, reverse=True):
            print(similar[l])

# Using geoLocator, search for missing coordinates in a supervised manner
# Automatically add all new coordinates to lat_longs.tsv and store in output folder for manual replacing of the old file
def search_missing_latlongs(missing_latlongs):
    geolocator = Nominatim(user_agent="hello@nextstrain.org")
    new_lat_longs = []

    for country in missing_latlongs["location"]:
        print("# " + country + " #")
        for division in missing_latlongs["location"][country]:
            print("\ndivision: " + division)
            for location in missing_latlongs["location"][country][division]:
                full_location = location + ", " + division + ", " + country
                new_lat_longs.append(find_place("location", location, full_location, geolocator))
        print()

    for region in missing_latlongs["division"]:
        for country in missing_latlongs["division"][region]:
            print("# " + country + " #")
            for division in missing_latlongs["division"][region][country]:
                full_division = division + ", " + country
                new_lat_longs.append(find_place("division", division, full_division, geolocator, region))
            print()

    for region in missing_latlongs["country"]:
        for country in missing_latlongs["country"][region]:
            new_lat_longs.append(find_place("country", country, country, geolocator))

    auto_add_lat_longs(new_lat_longs)

# Suggest geoLocator hits to the user and ask for confirmation or alternative spellings
def find_place(geo_level, place, full_place, geolocator, region = "*"):
    typed_place = full_place
    redo = True
    tries = 0
    while redo == True:
        if tries < 5:
            try:
                new_place = geolocator.geocode(typed_place, language='en')
            except:
                tries += 1
                continue
        else:
            new_place = None
            tries = 0

        if str(new_place) == 'None':
            print("\nCurrent place for missing {}:\t".format(geo_level) + full_place)
            print("The place as currently written could not be found.")
            answer = 'n'
        else:
            new_place_string = new_place.address
            full_place_string = full_place
            for level in full_place.split(", "):
                if clean_string(level) in clean_string(new_place_string):
                    new_place_string = bold(level).join(new_place_string.split(level))
                    full_place_string = bold(level).join(full_place_string.split(level))
            for level in new_place.address.split(", "):
                if clean_string(level) in clean_string(full_place_string):
                    full_place_string = bold(level).join(full_place_string.split(level))
                    new_place_string = bold(level).join(new_place_string.split(level))

            print("\nCurrent place for missing {}:\t".format(geo_level) + full_place_string)

            print("Geopy suggestion: "+ new_place_string)

            if geo_level != "division":
                answer = input('Is this the right place [y/n]? ')
            else:
                answer = input('Is this the right place (a - alter division level) [y/n/a]? ')

        if answer.lower() == 'y':
            coordinates = (geo_level + "\t" + place + "\t" + str(new_place.latitude) + "\t" + str(new_place.longitude))
            redo = False

        elif geo_level == "division" and answer.lower() == "a":
            division2 = input("Type correct division to produce corrective rule: ")
            (division, country) = full_place.split(", ")
            print(bold("/".join([region, country, division, ""]) + "\t" + "/".join([region, country, division2, division])))
            redo = False
            coordinates = ("location" + "\t" + place + "\t")

        else:
            # Let the user correct/have more detail for what's typed
            print("For: "+full_place)
            typed_place =  input("Type a more specific place name or 'NA' to leave blank: ")
            if typed_place.lower() == 'na':
                coordinates = (geo_level + "\t" + place + "\t")
                redo = False

    #print(coordinates)
    return coordinates

# Add new coordinates to lat_longs.tsv and sort the file before storing in the output folder
def auto_add_lat_longs(new_lat_longs):
    with open("defaults/lat_longs.tsv") as f:
        lat_longs_old = f.readlines()

    lat_longs = lat_longs_old + [l + "\n" for l in new_lat_longs if len(l.split("\t")) == 4]

    dataset = {"location": [], "division": [], "country": [], "region": []}
    for line in lat_longs:
        if line == "\n":
            continue
        dataset[line.split("\t")[0]].append(line)

    lat_longs_sorted = []

    regions_list = []
    for type in dataset:
        no_special_char = {clean_string(dataset[type][i].split("\t")[1]): i for i in range(len(dataset[type]))}
        for line in sorted(no_special_char):
            i = no_special_char[line]
            line_orig = dataset[type][i]
            if line_orig.startswith("country") and line_orig.split("\t")[1] in region_order:
                regions_list.append(line_orig)
                continue
            lat_longs_sorted.append(line_orig)
        if type == "country":
            lat_longs_sorted.append("\n")
            lat_longs_sorted += regions_list
        lat_longs_sorted.append("\n")

    if lat_longs_sorted != lat_longs_old:
        with open(path_to_output_files + latlongs_file, "w") as f:
            for line in lat_longs_sorted:
                f.write(line)
        print(bold("\nNew lat_longs written out to " + path_to_output_files + latlongs_file + ". Remember to replace the old file in " + path_to_default_files + "."))
    else:
        print("No changes to " + latlongs_file + ".")

# Given either the new lat_longs file (if new coordinates were added in this iteration of the script) or the old file,
# Construct the color_ordering.tsv file based on the data dictionary. Only places with existing coordinates are added.
# Countries and divisions are sorted according to the coordinates, locations are sorted by the alphabet. Regions have
# a fixed ordering.
def build_ordering(data, new_latlongs):
    ordering = read_ordering(path_to_default_files)

    if new_latlongs:
        latlongs_path = path_to_output_files
    else:
        latlongs_path = path_to_default_files
    latlongs = read_latlongs(latlongs_path)

    # Drop all empty locations
    data_clean = {}
    for region in data:
        data_clean[region] = {}
        for country in data[region]:
            data_clean[region][country] = {}
            for division in data[region][country]:
                data_clean[region][country][division] = []
                for location in data[region][country][division]:
                    if location != "":
                        data_clean[region][country][division].append(location)

    with open(path_to_output_files + ordering_file, "w") as out:
        for hierarchy in ordering:
            if hierarchy not in ["region", "country", "division", "location"]:
                for l in ordering[hierarchy]:
                    out.write(hierarchy + "\t" + l + "\n")
            else:
                for region in region_order:

                    if hierarchy == "region":
                        out.write("region\t" + region + "\n")
                    else:

                        out.write("\n# " + region + "\n")
                        for country in sort_by_coordinates(data_clean[region], latlongs["country"]):

                            if hierarchy == "country":
                                out.write("country\t" + country + "\n")
                            else:

                                if hierarchy == "location":
                                    if sum([len(data_clean[region][country][d]) for d in data_clean[region][country]]) > 0:  # only write country as a comment if there is data following it
                                        out.write("\n### " + country)

                                if hierarchy == "division":
                                    if len(data_clean[region][country]) > 0:
                                        out.write("\n### " + country + "\n")

                                for division in sort_by_coordinates(data_clean[region][country], latlongs["division"]):

                                    if hierarchy == "division":
                                        out.write("division\t" + division + "\n")
                                        continue

                                    if len(data_clean[region][country][division]) > 0:  # only write division as a comment if there is data following it
                                        out.write("\n# " + division + "\n")

                                    for location in sorted(data_clean[region][country][division]):
                                        out.write("location\t" + location + "\n")

                    if hierarchy == "location" or hierarchy == "division":
                        out.write("\n################\n")

            out.write("\n################\n\n\n")

    new_ordering = read_ordering(path_to_output_files)

    if not new_ordering == ordering:
        print(bold("Attention: " + ordering_file + " was altered! Remember to replace the old file in " + path_to_default_files + "."))
    else:
        print("No changes to " + ordering_file + ".")

# Sort a list of divisions or countries by latitude or longitude (whichever has the larger range)
def sort_by_coordinates(data, coordinates):
    max_lat = -90
    min_lat = 90
    max_long = -150
    min_long = 150
    for loc in data:
        if loc in coordinates:
            (lat, long) = coordinates[loc]
            max_lat = max(max_lat, lat)
            min_lat = min(min_lat, lat)
            max_long = max(max_long, long)
            min_long = min(min_long, long)

    index = 1
    if (max_lat - min_lat) > (max_long - min_long):
        index = 0

    loc_per_coord = {}
    for loc in data:
        if loc in coordinates:
            coord = coordinates[loc][index]
            if coordinates[loc][index] in loc_per_coord:
                loc_per_coord[coord].append(loc)
            else:
                loc_per_coord[coord] = [loc]
    sorted_locs = []
    for coord in sorted(loc_per_coord):
        sorted_locs.extend(loc_per_coord[coord])
    return sorted_locs

# Collect all stored annotations in a categorized manner, sorted by topics (e.g. geography-related) and annotation type
# (e.g. location, division, country...)
def read_annotations(annotationsFile, gisaid):
    types = {"geography": ["location", "division", "country", "region", "division_exposure", "country_exposure", "region_exposure"], "special": ["sampling_strategy", "date", "host", "strain"], "paper": ["title", "paper_url"]}
    types_inverted = {t:section for section, type in types.items() for t in type}
    annotations = {"comments": [], "geography": {}, "special": {}, "paper": {}}

    with open(path_to_annotations + annotationsFile) as f:
        line = f.readline()
        while line:
            if line.startswith("#"):
                annotations["comments"].append(line.strip())
            else:
                l = line.strip().split("\t")
                if line.endswith("\t\n"):
                    l.append("")
                if gisaid:
                    if len(l) != 4:
                        print("Invalid annotation length (annotation deleted): " + line.strip())
                        line = f.readline()
                        continue
                    else:
                        id = l[0] + "\t" + l[1]
                        type = l[2]
                        content = l[3]
                else:
                    if len(l) != 3:
                        print("Invalid annotation: " + line.strip())
                        line = f.readline()
                        continue
                    else:
                        id = l[0]
                        type = l[1]
                        content = l[2]
                if type not in types_inverted:
                    print("Invalid annotation type (annotation deleted): " + line.strip())
                else:
                    section = types_inverted[type]
                    if id not in annotations[section]:
                        annotations[section][id] = {}
                    if type in annotations[section][id]:
                        print("Duplicate annotation (first annotation deleted): " + line.strip() + " vs. " + type + "\t" + annotations[section][id][type])
                    annotations[section][id][type] = content
            line = f.readline()
    return annotations

# Given applied geoLocationRules and manualAnnotationRules, search the metadata for all affected strains
def create_annotations(metadata_filename, applied_rules_geoLocation, applied_rules_manual, gisaid):
    geoLocationAnnotations = {}
    manualAnnotations = {}

    with open(path_to_metadata + metadata_filename) as f:
        header = f.readline().split("\t")
        country_i = header.index("country")
        region_i = header.index("region")
        division_i = header.index("division")
        location_i = header.index("location")
        strain_i = header.index("strain")
        gisaid_epi_isl_i = header.index("gisaid_epi_isl")
        genbank_accession_i = header.index("genbank_accession")
        host_i = header.index("host")

        line = f.readline()
        while line:
            l = line.split("\t")
            country = l[country_i]
            region = l[region_i]
            division = l[division_i]
            location = l[location_i]
            strain = l[strain_i]

            if gisaid:
                id = l[gisaid_epi_isl_i]
            else:
                id = l[genbank_accession_i]

            if (region, country, division, location) in applied_rules_geoLocation:
                geoLocationAnnotations[id] = (region, country, division, location), applied_rules_geoLocation[(region, country, division, location)], strain

            if (region, country, division, location) in applied_rules_manual:
                manualAnnotations[id] = (region, country, division, location), applied_rules_manual[(region, country, division, location)], strain

            line = f.readline()

    return geoLocationAnnotations, manualAnnotations

# Compare old annotations with new alterations to the data and flag & adjust conflicts
# (Since annotations overwrite geoLocationRules in the ncov-ingest pipeline, there is a need to find all annotations
# that conflict with new rules and adjust the annotation accordingly)
# For geoLocationRules, only test whether there are conflicting annotations that need adjustment.
# For manualAnnotationRules, also produce new annotations.
def find_conflicting_annotations(annotations, geoLocationAnnotations, manualAnnotations, gisaid):
    for id in annotations["geography"]:
        EPI_ISL = id.split("\t")[-1]
        for ruleSet in [geoLocationAnnotations, manualAnnotations]:
            if EPI_ISL in ruleSet:
                (region2, country2, division2, location2) = ruleSet[EPI_ISL][1]
                annotations_correct = {"region": region2, "country": country2, "division": division2, "location": location2}
                for type in annotations_correct:
                    if type in annotations["geography"][id]:
                        name0 = annotations["geography"][id][type]
                        comment = ""
                        if "#" in name0:
                            (name0, comment) = name0.split(" #")
                        if name0 != annotations_correct[type]:
                            print(f"Conflicting annotation: {id}\t{bold(type + ' ' + name0)} will be replaced with {bold(annotations_correct[type])}")
                            annotations["geography"][id][type] = annotations_correct[type]
                            if comment != "":
                                annotations["geography"][id][type] += " #" + comment

    for EPI_ISL in manualAnnotations:
        (region, country, division, location) = manualAnnotations[EPI_ISL][0]
        (region2, country2, division2, location2) = manualAnnotations[EPI_ISL][1]
        strain = manualAnnotations[EPI_ISL][2]
        if gisaid:
            id = strain + "\t" + EPI_ISL
        else:
            id = EPI_ISL
        annotations_correct = {"region": (region, region2), "country": (country, country2), "division": (division, division2), "location": (location, location2)}
        for type in annotations_correct:
            if annotations_correct[type][0] != annotations_correct[type][1]:
                if id not in annotations["geography"]:
                    annotations["geography"][id] = {}
                if type not in annotations["geography"][id]:
                    annotations["geography"][id][type] = annotations_correct[type][1] + " # previously " +  annotations_correct[type][0]

    return annotations

clade_dates = {
    "19A": "2019-12-01",
    "19B": "2019-12-01",
    "20A": "2020-01-20",
    "20A.EU2": "2020-02-15",
    "20B": "2020-02-14",
    "20C": "2020-02-25",
    "20D": "2020-03-12",
    "20E (EU1)": "2020-05-27",
    "20F": "2020-05-24",
    "20G": "2020-06-11",
    "20H (Beta, V2)": "2020-08-10",
    "20I (Alpha, V1)": "2020-09-20",
    "20J (Gamma, V3)": "2020-10-29",
    "21A (Delta)": "2020-10-30",
    "21B (Kappa)": "2020-10-30",
    "21C (Epsilon)": "2020-08-03",
    "21D (Eta)": "2020-11-21",
    "21E (Theta)": "2021-01-10",
    "21F (Iota)": "2020-11-20",
    "21G (Lambda)": "2021-01-05",
    "21H": "2021-01-05",
}
# Check for special cases where annotations need to be introduced, e.g. special characters in strain names, or adjustment to "Mink"
# Also check for sampling dates that are too early for the assigned clade and auto-add to exclude
def special_metadata_checks(metadata_filename, annotations, gisaid):
    special_annotations = {}

    unknown_clades = []
    with open(path_to_metadata + metadata_filename) as f:
        header = f.readline().strip().split("\t")
        strain_i = header.index("strain")
        gisaid_epi_isl_i = header.index("gisaid_epi_isl")
        genbank_accession_i = header.index("genbank_accession")
        host_i = header.index("host")
        clade_i = header.index("Nextstrain_clade")
        date_i = header.index("date")
        clock_deviation_i = header.index("clock_deviation")

        line = f.readline()
        while line:
            l = line.strip().split("\t")
            host = l[host_i]
            strain = l[strain_i]

            if gisaid:
                id = strain + "\t" + l[gisaid_epi_isl_i]
            else:
                id = l[genbank_accession_i]

            # Check for special cases where annotations need to be introduced, e.g. special characters in strain names, or adjustment to "Mink"
            if host == "Neovison vison" or host == "Mustela lutreola":
                print("Adjust host " + host + " to Mink")
                if id not in special_annotations:
                    special_annotations[id] = {}
                special_annotations[id]["host"] = "Mink # previously " + host

            problematic_char = ["'", "`"]
            for c in problematic_char:
                if c in strain:
                    strain2 = strain.replace(c, "-")
                    print("Adjust strain " + strain + " to " + strain2)
                    if id not in special_annotations:
                        special_annotations[id] = {}
                    special_annotations[id]["strain"] = strain2 + " # previously " + strain

            line = f.readline()

    for id in special_annotations:
        if id not in annotations["special"]:
            annotations["special"][id] = {}
        for type in special_annotations[id]:
            if type in annotations["special"][id]:
                if annotations["special"][id][type] != special_annotations[id][type]:
                    print("Conflicting annotation: " + id + "\t" + bold(type + " " + annotations["special"][id][type]) + " will be replaced with " + bold(special_annotations[id][type]))
            annotations["special"][id][type] = special_annotations[id][type]

    return annotations

# Write the adjusted annotation set to the output folder in a sorted manner
def write_annotations(annotations, annotationsFile):
    with open(path_to_output_files + annotationsFile, "w") as out:
        for section in annotations:
            if section == "comments":
                for line in sorted(annotations[section]):
                    out.write(line + "\n")
            else:
                for id in sorted(annotations[section]):
                    for type in sorted(annotations[section][id]):
                        out.write(id + "\t" + type + "\t" + annotations[section][id][type] + "\n")


path_to_metadata = "data/"
path_to_default_files = "defaults/" # Contains color_ordering.tsv and lat_longs.tsv
path_to_config_files = "scripts/curate_metadata/config_curate_metadata/"
path_to_output_files = "scripts/curate_metadata/output_curate_metadata/"
path_to_annotations = "../ncov-ingest/source-data/" # Contains gisaid and genbank annotation files
Path(path_to_output_files).mkdir(parents=True, exist_ok=True)
Path(path_to_config_files).mkdir(parents=True, exist_ok=True)

gisaid_metadata_file = "downloaded_gisaid.tsv"
genbank_metadata_file = "metadata_genbank.tsv"
abbreviations_file = "abbreviations.txt"
accepted_additions_file = "acceptedExposureAdditions.txt"
geoLocationRules_file = "geoLocationRules.txt"
manualAnnotationRules_file = "manualAnnotationRules.txt"
internationalExceptions_file = "internationalExceptions.txt"
gisaidAnnotationsFile = "gisaid_annotations.tsv"
genbankAnnotationsFile = "genbank_annotations.tsv"
ordering_file = "color_ordering.tsv"
latlongs_file = "lat_longs.tsv"
exclude_file = "exclude.txt"

if not os.path.exists(path_to_config_files + geoLocationRules_file):
    with open(path_to_config_files + geoLocationRules_file, 'w'): pass

if not os.path.exists(path_to_config_files + manualAnnotationRules_file):
    with open(path_to_config_files + manualAnnotationRules_file, 'w'): pass

manualAnnotationsDelimiter = ["/", ","]
region_order = ["Asia", "Oceania", "Africa", "Europe", "South America", "North America"]

if __name__ == '__main__':

    print("\n===============================================\n")
    # Collect metadata from gisaid and genbank in a joint, multi-level dictionary accessible via data[region][country][division] = list_of_locations
    data = {}
    # count occurences
    geo_location_occurences = {"region": Counter(), "country": Counter(), "division": Counter(), "location": Counter()}
    print("Reading GISAID metadata...")
    data, geo_location_occurences = read_metadata(gisaid_metadata_file, data, geo_location_occurences)

    print("Reading GenBank metadata...")
    data, geo_location_occurences = read_metadata(genbank_metadata_file, data, geo_location_occurences, genbank=True)

    print("\n===============================================\n")
    # Add countries, regions and divisions that only appear in the exposure info from gisaid metadata
    # (We don't have any exposure info for genbank)
    print("Checking exposure consistency...")
    data = read_exposure(data, gisaid_metadata_file, accepted_additions_file)

    print("\n===============================================\n")
    print("Applying new geoLocationRules...")
    # Apply new general rules that will be copied over into ncov-ingest/source-data/gisaid_geoLocationRules.tsv
    data, applied_rules_geoLocation = apply_rules(data, geoLocationRules_file)

    print("\nApplying manualAnnotationRules...")
    # Apply new specific rules that will be translated into annotations
    # (needed for cases where geoLocationRules won't work, e.g. for locations containing the delimiter char '/')
    data, applied_rules_manual = apply_rules(data, manualAnnotationRules_file, delimiter = manualAnnotationsDelimiter)

    print("\n(Applying adjustments to international cruiseships...)")
    # Some cruiseships need adjustment to one country only so they don't appear multiple times in color_ordering.tsv
    data, applied_rules_exceptions = apply_rules(data, internationalExceptions_file, print_rules = False)

    print("\n===============================================\n")
    # Check for names that appear as division and location
    print("Checking for division inconsistencies...")
    check_division_inconsistency(data)

    print("\n===============================================\n")
    # Check for locations and divisions that appear multiple times
    print("Checking for location duplicates...")
    locations_duplicates = check_duplicates(data, abbreviations_file)

    print("\n===============================================\n")
    # List all locations, divisions, coutnries and regions that have no lat_longs in defaults/lat_longs.tsv
    print("Checking for missing lat_longs...")
    missing_latlongs = missing_coordinates(data, path_to_default_files, geo_location_occurences)
    print_missing_places(missing_latlongs)

    print("\n===============================================\n")
    # For all missing places, search for names that are similar or identical (when not considering special characters)
    # For locations, also consider special cases for USA counties
    print("Checking for similar names...")
    search_similar_names(data, missing_latlongs, locations_duplicates)

    print("\n===============================================\n")
    # Using geopy, search coordinates for all missing places and auto-sort them into lat_longs.tsv
    answer = input("Would you like to search for missing lat_longs [y/n]? ")
    if answer == "y":
        print("\nSearching for lat_longs...")
        search_missing_latlongs(missing_latlongs)
    else:
        print("\nAuto-sorting " + latlongs_file + "...")
        auto_add_lat_longs([])

    print("\n===============================================\n")
    # Reconstruct the color_ordering.tsv file based on the data dictionary and old and new lat_longs
    print("Constructing new color_ordering file...")
    build_ordering(data, answer == "y")


    print("\n===============================================\n")
    # Collect all known annotations and sort by type & strain
    print("\nReading annotations...")
    annotations_gisaid = read_annotations(gisaidAnnotationsFile, gisaid=True)
    annotations_open = read_annotations(genbankAnnotationsFile, gisaid=False)

    answer2 = input("Would you like to check annotations for conflicts with geoLocationRules and produce manualAnnotations [y/n]? ")
    if answer2 == "y":

        # Collect all strains for which new rules apply
        print("\n----------\n")
        print("Applying new geoLocationRules and manualAnnotationRules to the metadata...")
        geoLocationAnnotations_gisaid, manualAnnotations_gisaid = create_annotations(gisaid_metadata_file, applied_rules_geoLocation, applied_rules_manual, gisaid = True)
        geoLocationAnnotations_open, manualAnnotations_open = create_annotations(genbank_metadata_file, applied_rules_geoLocation, applied_rules_manual, gisaid=False)

        # Compare affected strains with annotations and search for conflicts (-> update annotations to fit new rules)
        # Also insert new annotations created by manualAnnotationRules
        print("\n----------\n")
        print("Searching for conflicting annotations and adding manualAnnotationRules...")
        annotations_gisaid = find_conflicting_annotations(annotations_gisaid, geoLocationAnnotations_gisaid, manualAnnotations_gisaid, gisaid = True)
        annotations_open = find_conflicting_annotations(annotations_open, geoLocationAnnotations_open, manualAnnotations_open, gisaid = False)

    print("\n----------\n")
    # Perform special checks on the metadata, e.g. check for Mink host consistency, check if date is consistent with clade...
    answer3 = input("Would you like to perform additional metadata checks (e.g. date, host, strain name) [y/n]? ")
    if answer3 == "y":
        print("Traversing metadata...")
        annotations_gisaid = special_metadata_checks(gisaid_metadata_file, annotations_gisaid, gisaid = True)
        annotations_open = special_metadata_checks(genbank_metadata_file, annotations_open, gisaid = False)


    # Sort and write updated annotation files to output folder
    print("\n----------\n")
    print("Writing updated annotation files to " + path_to_output_files + "...")
    write_annotations(annotations_gisaid, gisaidAnnotationsFile)
    write_annotations(annotations_open, genbankAnnotationsFile)

    with open(path_to_annotations + gisaidAnnotationsFile, "r") as f:
        annot_gisaid_old = f.read()
    with open(path_to_output_files + gisaidAnnotationsFile, "r") as f:
        annot_gisaid_new = f.read()
    if annot_gisaid_old != annot_gisaid_new:
        print(bold("Attention: " + gisaidAnnotationsFile + " was altered! Remember to replace the old file in " + path_to_annotations + "."))
    else:
        print("No changes to " + gisaidAnnotationsFile + ".")

    with open(path_to_annotations + genbankAnnotationsFile, "r") as f:
        annot_open_old = f.read()
    with open(path_to_output_files + genbankAnnotationsFile, "r") as f:
        annot_open_new = f.read()
    if annot_open_old != annot_open_new:
        print(bold("Attention: " + genbankAnnotationsFile + " was altered! Remember to replace the old file in " + path_to_annotations + "."))
    else:
        print("No changes to " + genbankAnnotationsFile + ".")

