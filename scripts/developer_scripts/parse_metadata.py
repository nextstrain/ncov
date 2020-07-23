from os import listdir
from difflib import SequenceMatcher

# Things to make things recogised as Cruise ships & ignored/special treatment
cruise_abbrev = ["Grand Princess", "Cruise", "cruise", "Diamond Princess"]

#path to files used in the script
path_to_script_files = "scripts/developer_scripts/"

def bold(s):
    return('\033[1m' + s + '\033[0m')

################################################################################
# Utils for reading files
################################################################################

# Read files which store duplicates, variants etc.
def read_local_file(file_name): #TODO: how will final file structure look like? Also, combine everything into one file for compactness?

    path_file_name = path_to_script_files + file_name

    with open(path_file_name) as myfile:
        file_content = myfile.readlines()

    first_files = [path_to_script_files+fi for fi in ["sequences_exclude.txt", "duplicates.txt", "accepted_exposure_additions.txt"]]

    if path_file_name in first_files: #simple list
        return [line.strip() for line in file_content[1:]]

    second_files = [path_to_script_files+fi for fi in ["variants.txt", "wrong_regions.txt", "abbreviations.txt", "false_divisions.txt", ] ]

    if path_file_name in second_files: #dictionary, keys seaprated from content with tabs
        content = {}
        for line in file_content[1:]:
            l = line.strip().split("\t")
            if l[0] in content:
                print("Attention, duplicate found while reading " + file_name + ": " + l[0] + " -> " + l[1] + ", " + content[l[0]])
            content[l[0]] = l[1]
        return content

# Read ordering and lat_longs file and return as dictionary:
def read_geography_file(file_name):
    lat_longs = ("lat_longs" in file_name)
    with open(file_name) as myfile:
        data_file = myfile.readlines()

    if lat_longs:
        # dictionary containing all locations, divisions ety. as dict, linking name to coordinates
        data = {"location": {}, "division": {}, "country": {}, "region": {}, "recency": {}}
    else:
        # dictionary containing all locations, divisions etc. as lists
        data = {"location": [], "division": [], "country": [], "region": [], "recency": []}

    for line in data_file:
        if line == "\n":
            continue
        l = line.strip().split("\t")
        if l[0][:1] == "#": #if a comment - ignore!
            continue
        type = l[0] #location, division etc
        name = l[1]
        if name not in data[type]:
            if lat_longs:
                data[type][name] = (float(l[2]), float(l[3]))
            else:
                data[type].append(name)
        else:
            s = "ordering"
            if lat_longs:
                s = "lat_longs"
            print("Duplicate in " + s + "? (" + l[0] + " " + l[1] + ")\n") #if already in the dictionary, print warning

    return data



################################################################################
# Step 1: Collection of data from metadata file in hierarchical manner
################################################################################

##### Step 1.1: Collection of all standard, non-exposure related data

# Read the metadata and return in dictionary format data[region][country][division][location] = list of all strains + GISAID id with this combination
def read_metadata(metadata):
    data = {}

    sequences_exclude = read_local_file("sequences_exclude.txt") #TODO: temporary solution to exclude two inconsistent Thailand samples

    for line in metadata[1:]:
        l = line.split("\t")
        region = l[5]
        country = l[6]
        division = l[7]
        location = l[8]
        id = l[2]
        strain = l[0]

        if location == "Unknown" or location == "UNKNOWN" or location == "unknown": #TODO: separate
            additions_to_annotation.append(strain + "\t" + id + "\tlocation\t# previously " + location)
            location = ""

        if region == "United Kingdom": #TODO: separate this, make it more applicable for other countries
            additions_to_annotation.append(strain + "\t" + id + "\tregion\tEurope")
            additions_to_annotation.append(strain + "\t" + id + "\tcountry\tUnited Kingdom")
            additions_to_annotation.append(strain + "\t" + id + "\tdivision\t" + country)
            if division != country:
                print("UK issue: division " + division + " != country " + country)
            region = "Europe"
            country = "United Kingdom"

        if strain in sequences_exclude:
            continue

        if region not in data:
            data[region] = {}
        if country not in data[region]:
            data[region][country] = {}
        if division not in data[region][country]:
            data[region][country][division] = {}
        if location not in data[region][country][division]:
            data[region][country][division][location] = []
        data[region][country][division][location].append(strain + "\t" + id)  # store strain and id of each seq with this combination of region/country/division/location
    additions_to_annotation.append("\n=============================\n")
    return data


##### Step 1.2: Collection of regions, countries and divisions of exposure

# TODO: many inconsistencies to resolve
def read_exposure(data, metadata):
    # divisions and countries that are accepted additions to the metadata
    accepted_additions = read_local_file("accepted_exposure_additions.txt")

    print("\n=============================\n")
    print("Travel history includes:")

    bad_div = {}
    bad_ctry = {}

    for line in metadata[1:]:
        l = line.split("\t")
        region2 = l[9]
        country2 = l[10]
        division2 = l[11]
        id = l[2]
        strain = l[0]

        s = division2 + " (" + country2 + ", " + region2 + ")"
        s2 = country2 + " (" + region2 + ")"

        if s2 in bad_ctry:
            bad_ctry[s2].append(line.strip())
        else:
            if country2 not in data[region2]:
                if s2 not in accepted_additions and country2 != region2:
                    bad_ctry[s2] = [line.strip()]
                else:
                    data[region2][country2] = {}
                    #print("Added country " + bold(s2) + " to the dataset") #optional confirmation of added countries
        if s in bad_div:
            bad_div[s].append(line.strip())
        else:
            if country2 in data[region2]:
                if division2 not in data[region2][country2]:
                    if s not in accepted_additions and division2 != country2:
                        bad_div[s] = [line.strip()]
                    else:
                        data[region2][country2][division2] = []
                        #print("Added division " + bold(s) + " to the dataset") #optional confirmation of added divisions
    print("\n\nUnchecked travel histories: (consider adding to accepted_exposure_additions.txt)\n")
    for division in bad_div:
        print("Strains with unknown division " + bold(division))
        for l in bad_div[division]:
            print(l)
        print()

    print()
    for country in bad_ctry:
        print("Strains with unknown country " + bold(country))
        for l in bad_ctry[country]:
            print(l)
        print()
    print("\n=============================\n")

    return data


################################################################################
# Utils for manipulating metadata
################################################################################

# Correct the metadata dictionary in a given manner
# e.g. switch all locations and strains from a misspelled division to the correct division
# e.g. turn a certain false division into a location below the correct division, and move all connected strains
def correct_data(data, type, corrections): #TODO: add region correction (e.g. for Turkey, Georgia)

    if type == "country":
        for country in corrections:
            (region, country_correct) = corrections[country]
            if country_correct not in data[region]:
                data[region][country_correct] = {}
            for division in data[region][country]:
                if division not in data[region][country_correct]:
                    data[region][country_correct][division] = {}
                for location in data[region][country][division]:
                    if location not in data[region][country_correct][division]:
                        data[region][country_correct][division][location] = []
                    for strain in data[region][country][division][location]:
                        additions_to_annotation.append(strain + "\tcountry\t" + country_correct + " # previously " + country)
                        data[region][country_correct][division][location].append(strain)
            del data[region][country]

    if type == "division":
        for division in corrections:
            (region, country, division_correct) = corrections[division]
            if division_correct not in data[region][country]:
                data[region][country][division_correct] = {}
            for location in data[region][country][division]:
                if location not in data[region][country][division_correct]:
                    data[region][country][division_correct][location] = []
                for strain in data[region][country][division][location]:
                    additions_to_annotation.append(strain + "\tdivision\t" + division_correct + " # previously " + division)
                    data[region][country][division_correct][location].append(strain)
            del data[region][country][division]

    if type == "location":
        for location in corrections:
            (region, country, division, location_correct) = corrections[location]
            if location_correct not in data[region][country][division]:
                data[region][country][division][location_correct] = []
            for strain in data[region][country][division][location]:
                additions_to_annotation.append(strain + "\tlocation\t" + location_correct + " # previously " + location)
                data[region][country][division][location_correct].append(strain)
            del data[region][country][division][location]

    if type == "div_to_loc":
        for location in corrections:
            (region, country, division) = corrections[location]
            if division not in data[region][country]:
                data[region][country][division] = {}
            for sub_location in data[region][country][location]:
                if sub_location != "":
                    print("Attention, additional location assigned to false division: " + sub_location)
                if location not in data[region][country][division]:
                    data[region][country][division][location] = []
                for strain in data[region][country][location][sub_location]:
                    additions_to_annotation.append(strain + "\tdivision\t" + division + " # previously false division " + location)
                    additions_to_annotation.append(strain + "\tlocation\t" + location)
                    data[region][country][division][location].append(strain)
            del data[region][country][location]

    return data

# Search the ordering file for a similar name as the one given, and return it if the score is above a fixed threshold
def check_similar(ordering, name):
    diff_max = 0
    name_max = ""
    for name0 in ordering:
        diff = SequenceMatcher(None, name, name0).ratio()
        if diff > diff_max:
            diff_max = diff
            name_max = name0

    if diff_max > 0.7:
        return name_max
    return ""


################################################################################
# Step 2: Clean up data
################################################################################


##### Step 2.0:
def adjust_to_database(data): #TODO: temporary solution, needs reworking
    for region in data:
        for country in data[region]:
            if country + ".txt" in listdir(path_to_script_files + "country_ordering/"): #TODO: correct path?

                variants = {}
                with open(path_to_script_files + "country_ordering/Belgium_variants.txt") as myfile: #TODO: this could be prettier...
                    belgium_variants = myfile.readlines()
                for line in belgium_variants:
                    if line == "\n":
                        continue
                    l = line.strip().split("\t")
                    variants[l[0]] = l[1]

                with open(path_to_script_files + "country_ordering/" + country + ".txt") as myfile:
                    country_ordering = myfile.readlines()

                arrondissement_to_location = {}
                location_to_arrondissement = {}
                duplicates = {}

                for line in country_ordering:
                    if line == "\n" or "------" in line:
                        continue
                    if line.startswith("### "):
                        continue
                    if line.startswith("# "):
                        arrondissement = line.strip()[2:]
                        arrondissement_to_location[arrondissement] = []
                        continue

                    location = line.strip()
                    if location not in arrondissement_to_location[arrondissement]:
                        arrondissement_to_location[arrondissement].append(location)
                    if location in location_to_arrondissement:
                        if location_to_arrondissement[location] != arrondissement:
                            duplicates[location] = (arrondissement, location_to_arrondissement[location])
                    location_to_arrondissement[location] = arrondissement

                division_to_correct = {}
                div_to_loc = {}
                for division in data[region][country]:

                    if division == country:
                        continue

                    if division in duplicates:
                        print("Attention duplicate: " + bold(division) + " found in " + bold(duplicates[division][0]) + " and " + bold(duplicates[division][1]))
                        print("Suggestion: select one and adjust database by deleting duplicate (no better solution due to missing additional info")

                    if division in arrondissement_to_location: # given division is actually an arrondissement => no changes necessary
                        continue

                    if division in variants and variants[division] in arrondissement_to_location: # given division is an arrondissement, but missspelled => simple adjustment
                        division_to_correct[division] = (region, country, variants[division])
                        continue

                    if division in location_to_arrondissement:
                        div_to_loc[division] = (region, country, location_to_arrondissement[division])
                        continue

                    if division in variants and variants[division] in location_to_arrondissement:
                        division_to_correct[division] = (region, country, variants[division]) #first correct to properly spelled division
                        div_to_loc[variants[division]] = (region, country, location_to_arrondissement[variants[division]]) #then to location
                        continue
                    print("Missing division in " + country + " database: " + bold(division))

                data = correct_data(data, "division", division_to_correct)
                data = correct_data(data, "div_to_loc", div_to_loc)

    print("\n=============================\n")
    return data


##### Step 2.1: Apply all known variants stored in an external file variants.txt
def apply_variants(data): #TODO: currently, the file variants.txt doesn't distinguish between location or division - what if we want to correct only one type, not the other?
    variants = read_local_file("variants.txt")

    countries_to_switch = {}
    for region in data:
        for country in data[region]:
            if country in variants:
                country_correct = variants[country]
                print("Apply variant (country): " + bold(country) + " -> " + bold(country_correct))
                countries_to_switch[country] = (region, country_correct)

    data = correct_data(data, "country", countries_to_switch)

    divisions_to_switch = {}
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                if division in variants:
                    division_correct = variants[division]
                    print("Apply variant (division): " + bold(division) + " -> " + bold(division_correct))
                    divisions_to_switch[division] = (region, country, division_correct)

    data = correct_data(data, "division", divisions_to_switch)

    locations_to_switch = {}
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                for location in data[region][country][division]:
                    if location in variants:
                        location_correct = variants[location]
                        print("Apply variant (location): " + bold(location) + " -> " + bold(location_correct))
                        locations_to_switch[location] = (region, country, division, location_correct)

    data = correct_data(data, "location", locations_to_switch)

    additions_to_annotation.append("\n=============================\n")
    print("\n=============================\n")
    return data

def apply_typical_errors(data): #TODO: rename, maybe join with UK as region? also use correct_data()
    wrong_regions = read_local_file("wrong_regions.txt")

    countries_to_switch = {}
    for country in wrong_regions:
        correct_region = wrong_regions[country]
        for region in data:
            if region == correct_region:
                continue
            if country in data[region]:
                print("Found incorrect region " + bold(region) + " for country " + bold(country) + " (correct region: " + bold(correct_region) + ")" )

                for division in data[region][country]:
                    if division not in data[correct_region][country]:
                        data[correct_region][country][division] = {}
                    for location in data[region][country][division]:
                        if location not in data[correct_region][country][division]:
                            data[correct_region][country][division][location] = []
                        for strain in data[region][country][division][location]:
                            additions_to_annotation.append(strain + "\tregion\t" + correct_region + " # previously " + region)
                            data[correct_region][country][division][location].append(strain)
                del data[region][country]
    return data

##### Step 2.2 Check for "false" division that appear as location elsewhere (known cases stored in false_divisions.txt as well as checking for new cases)
def check_false_divisions(data):

    # Known false divisions
    div_as_loc_known = {}
    known_false_divisions = read_local_file("false_divisions.txt")
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                if division in known_false_divisions:
                    div_as_loc_known[division] = (region, country, known_false_divisions[division])
                    print("False division corrected: " + bold(division) + " (true division: " + bold(known_false_divisions[division]) + ")")
    data = correct_data(data, "div_to_loc", div_as_loc_known)


    # Check for unknown cases:
    div_as_loc = {}
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                for location in data[region][country][division]:
                    if location in data[region][country] and location != division:
                        div_as_loc[location] = (region, country, division)
                        print("Unknown location found as division: " + bold(location) + " (true division: " + bold(division) + ")")
                        print("(Suggestion: add " + location + " -> " + division + " to false_divisions.txt)")

    additions_to_annotation.append("\n=============================\n")
    print("\n=============================\n")

##### Step 2.3: Check for duplicate divisions/locations in different countries/divisions (known cases stored in duplicates.txt as well as checking for new cases)
def check_duplicate(data):

    #Check known duplicates
    # TODO: Only locations covered properly (divisions: only alert)
    duplicates = read_local_file("duplicates.txt")
    abbreviations = read_local_file("abbreviations.txt")

    duplicate_locations = {}
    for region in data:
        for country in list(data[region].keys()):
            for division in list(data[region][country]):
                for location in list(data[region][country][division]):
                    if location in duplicates:
                        print("Known duplicate detected: " + bold(location))
                        data = correct_data(data, "location", {location: (region, country, division, location + " " + abbreviations[division])})

    #Check for new cases
    division_to_country = {}
    location_to_division = {}
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                if division == "":
                    continue
                if division not in division_to_country:
                    division_to_country[division] = []
                division_to_country[division].append((country, region))
                for location in data[region][country][division]:
                    if location == "":
                        continue
                    if location not in location_to_division:
                        location_to_division[location] = []
                    location_to_division[location].append((division, country, region))

    print()

    cruise_ship_duplicates = 0
    #TODO: a bit chaotic, go over it again
    for division in division_to_country:
        if not any(x in division for x in cruise_abbrev): #ignore cruise ship ones
            if len(division_to_country[division]) > 1:
                if division_to_country[division][0][1] == division_to_country[division][1][1]:
                    s = ", ".join([country for (country, region) in division_to_country[division]])
                else:
                    s = ", ".join([country + " (" + region + ")" for (country, region) in division_to_country[division]])

                print("New duplicate division detected: " + bold(division + " (" + s + ")"))
        else:
            cruise_ship_duplicates = cruise_ship_duplicates + 1

    if cruise_ship_duplicates: print("("+str(cruise_ship_duplicates)+" cruise ship entries ignored for duplicate divisions)")
    

    cruise_ship_duplicates = 0
    cruis_ship_abbrev = 0
    for location in location_to_division:
        if len(location_to_division[location]) > 1:
            if location_to_division[location][0][1] == location_to_division[location][1][1]:
                if location_to_division[location][0][2] == location_to_division[location][1][2]:
                    s = ", ".join([division for (division, country, region) in location_to_division[location]])
                else:
                    s = ", ".join([division + " (" + country + ", " + region + ")" for (division, country, region) in location_to_division[location]])
            else:
                s = ", ".join([division + " (" + country + ")" for (division, country, region) in location_to_division[location]])


            if "Cruise" in location:
                cruise_ship_duplicates = cruise_ship_duplicates + 1
            else:
                print("New duplicate location detected: " + bold(location + " (in both " + s + ")"))
                print("(Suggestion: Add " + location + " to duplicates.txt)")

            
            for (division, country, region) in location_to_division[location]:
                if division not in abbreviations:
                    if not any(x in division for x in cruise_abbrev):
                        print("Attention: Missing abbreviation for " + bold(division) + " (Suggestion: add to abbreviations.txt)")
                    else:
                        cruis_ship_abbrev = cruis_ship_abbrev + 1
                        

    if cruise_ship_duplicates: print("("+str(cruise_ship_duplicates)+" cruise ship entries ignored for duplicate locations)")
    if cruis_ship_abbrev: print("("+ str(cruis_ship_abbrev) + " cruise ship entries ignored for missing state abbreviations)")
    print("\n=============================\n")

##### Step 2.4: Check for missing names in ordering and lat_longs as well as return a clean, reduced version of the metadata
def check_for_missing(data):
    data_clean = {}

    missing = {"country": [], "division": {}, "location": {}}
    clean_missing = {"country": [], "division": {}, "location": {}} # Same as above, but without formatting or notes

    for region in data:
        data_clean[region] = {}

        for country in data[region]:
            if country not in ordering["country"]:
                missing["country"].append(bold(country))
                clean_missing["country"].append(country)
            else:
                data_clean[region][country] = {}

            for division in data[region][country]:
                if division == "":
                    continue

                if division not in ordering["division"] or division not in lat_longs["division"]:
                    s = bold(division)
                    name0 = check_similar(ordering["division"], division)
                    if division not in ordering["division"] and division in lat_longs["division"]:
                        s = s + " (only missing in ordering => auto-added to color_ordering.tsv)"
                        if country not in data_clean[region]:
                            print("Conflict: division " + division + " should be added to color_ordering.tsv, but country " + country + " is missing from dataset")
                        else:
                            data_clean[region][country][division] = []
                    else: #only check for additional hints like "similar name" or "present as location" if not auto-added to color_ordering
                        if division in ordering["division"] and division not in lat_longs["division"]:
                            s = s + " (only missing in lat_longs)"
                        else:
                            if name0 != "":
                                s += " (similar name: " + name0 + " - consider adding to variants.txt)"
                            if division in ordering["location"] or division in lat_longs["location"]:
                                s = s + " (present as location)"
                    if country not in missing["division"]:
                        missing["division"][country] = []
                        clean_missing["division"][country] = []
                    missing["division"][country].append(s)
                    clean_missing["division"][country].append(division)

                else:
                    if country not in data_clean[region]:
                        print("Conflict: division " + division + " should be added to color_ordering.tsv, but country " + country + " is missing from dataset")
                    else:
                        data_clean[region][country][division] = []

                for location in data[region][country][division]:
                    if location == "":
                        continue

                    if location not in ordering["location"] or location not in lat_longs["location"]:
                        s = bold(location)
                        name0 = check_similar(ordering["location"], location)
                        if location not in ordering["location"] and location in lat_longs["location"]:
                            s = s + " (only missing in ordering => auto-added to color_ordering.tsv)"
                            if country not in data_clean[region]:
                                print("Conflict: location " + location + " should be added to color_ordering.tsv, but country " + country + " is missing from dataset")
                            else:
                                if division not in data_clean[region][country]:
                                    if not any(x in location for x in cruise_abbrev) and not any(x in division for x in cruise_abbrev):
                                        print("Conflict: location " + location + " should be added to color_ordering.tsv, but division " + division + " is missing from dataset")
                                else:
                                    data_clean[region][country][division].append(location)
                        else: #only check for additional hints like "similar name" or "present as division" if not auto-added to color_ordering
                            if name0 != "":
                                s += " (similar name: " + name0 + " - consider adding to variants.txt)"
                            if location in ordering["location"] and location not in lat_longs["location"]:
                                s = s + " (only missing in lat_longs)"
                            if location in ordering["division"] or location in lat_longs["division"]:
                                s = s + " (present as division)"

                        if country not in missing["location"]:
                            missing["location"][country] = {}
                            clean_missing["location"][country] = {}
                        if division not in missing["location"][country]:
                            if any(x in division for x in cruise_abbrev):
                                print("Cruise-associated division ignored ("+division+")")
                            else:
                                missing["location"][country][division] = []
                                clean_missing["location"][country][division] = []
                        if not any(x in location for x in cruise_abbrev) and not any(x in division for x in cruise_abbrev):
                            missing["location"][country][division].append(s)
                            clean_missing["location"][country][division].append(location)
                        else:
                            print("Cruise-associated location ignored ("+location+")")

                    else:
                        if country not in data_clean[region]:
                            print("Conflict: location " + location + " should be added to color_ordering.tsv, but country " + country + " is missing from dataset")
                        else:
                            if division not in data_clean[region][country]:
                                if not any(x in location for x in cruise_abbrev) and not any(x in division for x in cruise_abbrev):
                                    print("Conflict: location " + location + " should be added to color_ordering.tsv, but division " + division + " is missing from dataset")
                            else:
                                data_clean[region][country][division].append(location)

    if missing['location']:
        print("\n\nMissing locations:")
        for country in missing["location"]:
            print("# " + country + " #")
            for division in missing["location"][country]:
                print(division)
                for location in missing["location"][country][division]:
                    print("\tlocation\t" + location)
            print()
    else:
        print("No missing locations")

    if missing['division']:
        print("\nMissing divisions:")
        for country in missing["division"]:
            print("# " + country + " #")
            for division in missing["division"][country]:
                print("division\t" + division)
            print()
    else:
        print("No missing divisions")

    if missing['country']:
        print("\nMissing countries:")
        for country in missing["country"]:
            print(country)
    else:
        print("No missing countries")


    ##### Ask user if they want to look for lat-longs now, or end script for time being.

    find_lat_longs = input("\n\nWould you like to look for lat-longs for these places now? y or n \n(it's suggested to make any necessary file additions before this step): ")

    if find_lat_longs.lower() == 'y':

        from geopy.geocoders import Nominatim
        geolocator = Nominatim(user_agent="hello@nextstrain.org")
        new_lat_longs = []

        print("Getting lat-long for missing places:\n")

        for country in clean_missing["location"]:
            print("# " + country + " #")
            for division in clean_missing["location"][country]:
                print("\ndivision: "+division)
                for location in clean_missing["location"][country][division]:
                    if any(x in location for x in cruise_abbrev):
                        print(" One cruise ship location ignored ("+ location +").")
                        continue

                    full_location = location +", "+ division+", "+country

                    new_lat_longs.append(find_place("location", location, full_location, geolocator))
            print()


        for country in clean_missing["division"]:
            print("# " + country + " #")
            for division in clean_missing["division"][country]:
                print("division\t" + division)
                full_division = division+", "+country

                new_lat_longs.append(find_place("division", division, full_division, geolocator))
            print()

        for country in clean_missing["country"]:
            print(country)
            
            new_lat_longs.append(find_place("country", country, country, geolocator))

        print("\nNew locations to be written out: ")
        print(*new_lat_longs, sep='\n')

        with open(path_to_script_files+"new_lat-longs.tsv", 'w') as out:
            out.write("\n".join(new_lat_longs))
        print("New lat-longs written out to "+path_to_script_files+"new_lat-longs.tsv")


    print("\n=============================\n")
    return data_clean


# Get the geo-locator to find a possible location - returns result
# call with ex: 'Dallas, Texas, USA', geolocator
def ask_geocoder(full_unknown_place, geolocator):
    new_place = geolocator.geocode(full_unknown_place, language='en')
    return new_place

# Allows user to try typing different locations to get lat-long, or tell to leave blank
# Call with ex: 'location', 'Dallas', 'Dallas, Texas, USA', geolocator
def find_place(geo_level, place, full_place, geolocator):
    typed_place = full_place
    redo = True                    
    while redo == True:
        print("\nCurrent place for missing {}:\t".format(geo_level) + full_place)
        new_place = ask_geocoder(typed_place, geolocator)

        if str(new_place) == 'None':
            print("The place as currently written could not be found.")
            answer = 'n'
        else:
            print("Geopy suggestion: "+ new_place.address)
            answer = input('Is this the right place? Type y or n: ')

        if answer.lower() == 'y':
            answer = (geo_level + "\t" + place + "\t" + str(new_place.latitude) + "\t" + str(new_place.longitude))
            redo = False

        else:
            # Let the user correct/have more detail for what's typed
            print("For: "+full_place)
            typed_place =  input("Type a more specific place name or 'NA' to leave blank: ")
            if typed_place.lower() == 'na':
                print("Writing out a line with blank lat-long to be filled by hand")
                answer = (geo_level + "\t" + place + "\t")
                redo = False

    return answer

################################################################################
# Step 3: Storage of locations, divisions etc hierarchical manner
################################################################################

# Util needed to sort a given list of locations or divisions by their coordinated stored in lat_longs
# TODO: enable for countries as well
def sort_by_coordinates(data, coordinates):
    loc_per_coord = {}
    for loc in data:
        if loc in coordinates:
            if coordinates[loc][1] in loc_per_coord:
                loc_per_coord[coordinates[loc][1]].append(loc)
            else:
                loc_per_coord[coordinates[loc][1]] = [loc]
        else:
            print("Missing coordinates: " + bold(loc))
    sorted_locs = []
    for coord in sorted(loc_per_coord):
        sorted_locs.extend(loc_per_coord[coord])
    return sorted_locs

# Write a given hierarchy (location, division, country, region, recency) into the new ordering file.
# Sort locations and divisions by coordinates to retain proximity coloring
def write_ordering(data, hierarchy):
    mode = "a"
    if hierarchy == "location":
        mode = "w"

    with open(path_to_script_files+"color_ordering.tsv", mode) as out:
        if hierarchy == "recency":
            out.write("recency\tOlder\nrecency\tOne month ago\nrecency\tOne week ago\nrecency\t3-7 days ago\nrecency\t1-2 days ago\nrecency\tNew\n")
            return

        # Give fixed order of regions to retain the usual coloring order
        region_order = ["Asia", "Oceania", "Europe", "Africa", "South America", "North America"] #TODO: is this the correct order?

        for region in region_order:

            if hierarchy == "region":
                out.write("region\t" + region + "\n")
                continue

            out.write("\n# " + region + "\n")
            for country in sorted(data[region]): #TODO: would be nice to sort this by coordinate too, but would need to add most lat_longs first!

                if hierarchy == "country":
                    out.write("country\t" + country + "\n")
                    continue

                if hierarchy == "location":
                    if sum([len(data[region][country][d]) for d in data[region][country]]) > 0:  # only write country as a comment if there is data following it
                        out.write("\n### " + country)

                if hierarchy == "division":
                    if len(data[region][country]) > 0:
                        out.write("\n### " + country + "\n")

                for division in sort_by_coordinates(data[region][country], lat_longs["division"]):

                    if hierarchy == "division":
                        out.write("division\t" + division + "\n")
                        continue

                    if len(data[region][country][division]) > 0:  # only write division as a comment if there is data following it
                        out.write("\n# " + division + "\n")

                    for location in sort_by_coordinates(data[region][country][division], lat_longs["location"]):
                        out.write("location\t" + location + "\n")

            if hierarchy == "location" or hierarchy == "division":
                out.write("\n################\n")

        out.write("\n################\n\n\n")




################################################################################
# Step 0: Read data
################################################################################

# Read current metadata
#path_to_ncov = "../../" # TODO: adjust file structure properly
with open("data/metadata.tsv") as myfile:
    metadata = myfile.readlines()

# Read orderings and lat_longs
ordering = read_geography_file("defaults/color_ordering.tsv") #TODO: combine with read_local_files()?
lat_longs = read_geography_file("defaults/lat_longs.tsv")

# List that will contain all proposed annotations collected throughout the script
additions_to_annotation = []
with open("../ncov-ingest/source-data/gisaid_annotations.tsv") as myfile:
    annotations = myfile.read()


################################################################################
# Step 1: Collection of data from metadata file in hierarchical manner
################################################################################

##### Step 1.1: Collection of all standard, non-exposure related data

# Hierarchical ordering of all regions, countries, divisions and locations
# Each location (also empty ones) hold a list of all strains & GISAID IDs with this region+country+division+location
data = read_metadata(metadata)


##### Step 1.2: Collection of regions, countries and divisions of exposure
# In case some geographic units are only found in the exposure information of the metadata, iterate again over the metadata and add to the dataset
# Since travel history related entries are prone to errors, check for each entry whether it collides with already existing data.

# TODO: Currently commented out due to numerous inconsistencies
#data = read_exposure(data, metadata)


################################################################################
# Step 2: Clean up data
################################################################################

##### Step 2.0: Adjust the divisions and locations by comparing them to a known database - only accessible for Belgium at the moment
data = adjust_to_database(data)

##### Step 2.1: Apply all known variants stored in an external file variants.txt
data = apply_typical_errors(data) #TODO: do this earlier (before reading metadata), join with UK as region?
data = apply_variants(data)

##### Step 2.2 Check for "false" division that appear as location elsewhere (known cases stored in false_divisions.txt as well as checking for new cases)
check_false_divisions(data)

##### Step 2.3: Check for duplicate divisions/locations in different countries/divisions (known cases stored in duplicates.txt as well as checking for new cases)
check_duplicate(data)

##### Step 2.4: Check for missing names in ordering and lat_longs as well as return a clean, reduced version of the metadata
data = check_for_missing(data) # =====> From here on, strains are dropped, only region/country/division/location remain

################################################################################
# Step 3: Storage of locations, divisions etc hierarchical manner
################################################################################

write_ordering(data, "location")
write_ordering(data, "division")
write_ordering(data, "country")
write_ordering(data, "recency")
write_ordering(data, "region")


##### Bonus step: Print out all collected annotations - if considered correct, they can be copied by the user to annotations.tsv
# Only print line if not yet present
# Print warning if this GISAID ID is already in the file
annot_lines_to_write = []
for line in additions_to_annotation:
    if line in annotations:
        continue
    print(line)
    if "=" not in line:
        annot_lines_to_write.append(line)
    if len(line.split("\t")) == 4:
        if line.split("\t")[1] in annotations:
            print("Warning: " + line.split("\t")[1] + " already exists in annotations!")

with open(path_to_script_files+"new_annotations.tsv", 'w') as out:
    out.write("\n".join(annot_lines_to_write))
print("New annotation additions written out to "+path_to_script_files+"new_annotations.tsv")
