from os import listdir
from difflib import SequenceMatcher
from pathlib import Path


# Things to make things recogised as Cruise ships & ignored/special treatment
cruise_abbrev = ["Grand Princess", "Cruise", "cruise", "Diamond Princess"]

#path to files used in the script
path_to_config_files = "scripts/developer_scripts/config_files_parse_metadata/"
path_to_output_files = "scripts/developer_scripts/output_files_parse_metadata/"
Path(path_to_output_files).mkdir(parents=True, exist_ok=True)

def bold(s):
    return('\033[1m' + s + '\033[0m')

################################################################################
# Utils for reading files
################################################################################

# Read files which store duplicates, variants etc.
def read_local_file(file_name): #TODO: how will final file structure look like? Also, combine everything into one file for compactness?

    path_file_name = path_to_config_files + file_name

    with open(path_file_name) as myfile:
        file_content = myfile.readlines()

    first_files = [path_to_config_files+fi for fi in ["duplicates.txt", "accepted_exposure_additions.txt"]]

    if path_file_name in first_files: #simple list
        return [line.strip() for line in file_content[1:]]

    second_files = [path_to_config_files+fi for fi in ["wrong_regions.txt", "abbreviations.txt", "false_divisions.txt"] ]

    if path_file_name in second_files: #dictionary, keys seperated from content with tabs
        content = {}
        for line in file_content[1:]:
            l = line.strip().split("\t")
            if l[0] in content:
                print("Attention, duplicate found while reading " + file_name + ": " + l[0] + " -> " + l[1] + ", " + content[l[0]])
            content[l[0]] = l[1]
        return content

    third_files = [path_to_config_files+fi for fi in ["variants.txt", "international_exceptions.txt"] ]

    if path_file_name in third_files: #need two level-dict
        if path_file_name == path_to_config_files+"variants.txt":
            content = {'location': {}, 'division': {}, 'country': {}, 'region': {}}
        if path_file_name == path_to_config_files+"international_exceptions.txt":
            content = {'location': {}, 'division': {}}
        for line in file_content[1:]:
            if line == "\n":
                continue
            l = line.strip().split("\t")
            if line.endswith("\t\n"):
                l = [l[0], l[1], ""] #allow empty assignment of hierarchy (e.g. set location to blank)
            entry = l[2]
            if len(l) == 4:
                entry = (l[2],l[3])
            if l[0] not in content:
                content[l[0]] = {}
            if l[1] not in content[l[0]]: # allow duplicates (e.g. multiple "San Rafael" in different divisions)
                content[l[0]][l[1]] = []
            else: #check whether already existing variant has hierarchical ordering or not
                conflict = False
                for c in content[l[0]][l[1]]:
                    if type(c) is not tuple:
                        print("Warning: Variant " + str(entry) + " can not be applied due to the presence of another instance of this name in variants.txt without hierarchical ordering.")
                        conflict = True
                if  conflict:
                    continue
            content[l[0]][l[1]].append(entry)

        return content

    fourth_files = [path_to_config_files + fi for fi in ["manual_adjustments.txt"]]

    if path_file_name in fourth_files: # / and tab as separator
        content = {}
        for line in file_content[1:]:
            if line == "\n":
                continue
            l = line.strip().split("\t")[0].split("/") + line.strip().split("\t")[1].split("/")
            if len(l) < 8:
                for i in range(8-len(l)):
                    l.append("")
            k = "/".join(l[:4])
            c = "/".join(l[4:])
            if k in content:
                print("Attention, duplicate found while reading " + file_name + ": " + k + " -> " + c + ", " + content[k])
            content[k] = c
        return content



# Read ordering and lat_longs file and return as dictionary:
def read_geography_file(file_name, hierarchical = False):
    lat_longs = ("lat_longs" in file_name)
    with open(file_name) as myfile:
        data_file = myfile.readlines()

    if not hierarchical:
        if lat_longs:
            # dictionary containing all locations, divisions ety. as dict, linking name to coordinates
            data = {"location": {}, "division": {}, "country": {}, "region": {}}
        else:
            # dictionary containing all locations, divisions etc. as lists
            data = {"location": [], "division": [], "country": [], "region": []}
            color_ordering_other = {}

        for line in data_file:
            if line == "\n":
                continue
            l = line.strip().split("\t")
            if l[0][:1] == "#": #if a comment - ignore!
                continue
            type = l[0] #location, division etc
            name = l[1]

            if lat_longs:
                if name not in data[type]:
                    data[type][name] = (float(l[2]), float(l[3]))
                else:
                    print("Duplicate in lat_longs? (" + l[0] + " " + l[1] + ")\n")  # if already in the dictionary, print warning
            else:
                if type in data:
                    if name not in data[type]:
                        data[type].append(name)
                    else:
                        print("Duplicate in color_ordering? (" + l[0] + " " + l[1] + ")\n")  # if already in the dictionary, print warning
                else:
                    if type not in color_ordering_other:
                        color_ordering_other[type] = []
                    color_ordering_other[type].append(name)
        if lat_longs:
            return data
        else:
            return data, color_ordering_other
    else: #hierarchical structure of ordering for checking similar names only in the same country
        data = {"Asia": {}, "Oceania": {}, "Africa": {}, "Europe": {}, "South America": {}, "North America": {}}

        region = ""
        country = ""
        division = ""

        for line in data_file:
            if line == "\n":
                continue
            if line.startswith("###"):
                if len(line.split("### ")) > 1:  # country
                    country = line.strip().split("### ")[1]
                    if country not in data[region]:
                        data[region][country] = {}

            else:
                if line.startswith("#"):
                    if len(line.split("# ")) > 1:  # region or division
                        place = line.strip().split("# ")[1]
                        if place in data:
                            region = place
                        else:
                            division = place
                            if division not in data[region][country]:
                                data[region][country][division] = []

                else:
                    l = line.strip().split("\t")
                    type = l[0]  # location, division etc
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
    "ú": "u",
    "ś":"s",
    "ą":"q",
    "à":"a",
    "å":"a",
    "ł":"l",
    "-":" ",
    "î": "i",
    "ŕ": "r",
    "ľ": "l",
    "ď": "d"
}


def clean_string(s):
    s = s.lower()
    for c in replace_special_char:
        s = s.replace(c, replace_special_char[c])
    return s


def pre_sort_lat_longs(lat_longs):
    dataset = {"location": [], "division": [], "country": [], "region": []}
    regions = ["Africa", "Asia", "Europe", "North America", "Oceania", "South America"]
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
            if line_orig.startswith("country") and line_orig.split("\t")[1] in regions:
                regions_list.append(line_orig)
                continue
            lat_longs_sorted.append(line_orig)
        if type == "country":
            lat_longs_sorted.append("\n")
            lat_longs_sorted += regions_list
        lat_longs_sorted.append("\n")

    return lat_longs_sorted


#Function to support supervised addition of new entries into lat_longs. The user must review every new entry and approve it to be written into the lat_longs file. Ground truth lat_longs is not overwritten, but a copy is made in the developer_scripts folder.
def auto_add_lat_longs(new_lat_longs):

    with open("defaults/lat_longs.tsv") as f:
        lat_longs = f.readlines()
    lat_longs = pre_sort_lat_longs(lat_longs)
    for entry in new_lat_longs:
        if len(entry.split("\t")) < 4:
            continue
        correct_hierarchy = False
        for i in range(len(lat_longs)):
            if lat_longs[i] == "\n" and not correct_hierarchy:
                continue
            if lat_longs[i] != "\n" and entry[:4] != lat_longs[i][:4]: #first characters correspond to country, division, location etc.
                continue
            correct_hierarchy = True
            if lat_longs[i] != "\n" and clean_string(entry) > clean_string(lat_longs[i]):
                continue
            print("\n")
            for k in range(3):
                print(lat_longs[i-3+k].strip())
            print(bold(entry))
            for k in range(3):
                print(lat_longs[i+k].strip())
            answer = input("Approve of this new entry (y)?")
            if answer == "y":
                lat_longs = lat_longs[:i] + [entry + "\n" ] + lat_longs[i:]
            break

    local_file = path_to_output_files + "lat_longs.tsv"
    with open(local_file, "w") as f:
        for line in lat_longs:
            f.write(line)



################################################################################
# Step 1: Collection of data from metadata file in hierarchical manner
################################################################################

##### Step 1.1: Collection of all standard, non-exposure related data

# Read the metadata and return in dictionary format data[region][country][division][location] = list of all strains + GISAID id with this combination
def read_metadata(metadata, data, source):

    header = metadata[0].split("\t")
    region_i = header.index("region")
    country_i = header.index("country")
    division_i = header.index("division")
    location_i = header.index("location")

    strain_i = header.index("strain")

    if source == "gisaid":
        epi_i = header.index("gisaid_epi_isl")
    if source == "genbank":
        epi_i = header.index("genbank_accession")

    for line in metadata[1:]:
        l = line.split("\t")

        region = l[region_i]
        country = l[country_i]
        division = l[division_i]
        location = l[location_i]
        strain = l[strain_i]
        id = l[epi_i]

        host = l[14]
        if host == "Neovison vison" or host ==  "Mustela lutreola":
            print("Adjust host " + host + " to Mink")
            if source == "gisaid":
                additions_to_annotation.append(strain + "\t" + id + "\thost\tMink # previously " + host)
            if source == "genbank":
                additions_to_annotation.append(id + "\thost\tMink # previously " + host)
        problematic_char = ["'", "`"]

        for c in problematic_char:
            if c in strain:
                strain2 = strain.replace(c, "-")
                print("Adjust strain " + strain + " to " + strain2)
                if source == "gisaid":
                    additions_to_annotation.append(strain + "\t" + id + "\tstrain\t" + strain2 + " # previously " + strain)
                if source == "genbank":
                    additions_to_annotation.append(id + "\tstrain\t" + strain2 + " # previously " + strain)
            

        if region not in data:
            data[region] = {}
        if country not in data[region]:
            data[region][country] = {}
        if division not in data[region][country]:
            data[region][country][division] = {}
        if location not in data[region][country][division]:
            data[region][country][division][location] = []
        if source == "gisaid":
            data[region][country][division][location].append(strain + "\t" + id)  # store strain and id of each seq with this combination of region/country/division/location
        if source == "genbank":
            data[region][country][division][location].append(id)
    return data


##### Step 1.2: Collection of regions, countries and divisions of exposure

# TODO: many inconsistencies to resolve
def read_exposure(data, metadata):
    # divisions and countries that are accepted additions to the metadata
    accepted_additions = read_local_file("accepted_exposure_additions.txt")

    print("\n=============================\n")
    #print("Travel history includes:")

    bad_div = {}
    bad_ctry = {}

    header = metadata[0].split("\t")
    region_i = header.index("region_exposure")
    country_i = header.index("country_exposure")
    division_i = header.index("division_exposure")

    for line in metadata[1:]:
        l = line.split("\t")
        region2 = l[region_i]
        country2 = l[country_i]
        division2 = l[division_i]

        if region2 == "United Kingdom": #TODO: separate this, make it more applicable for other countries
            region2 = "Europe"
            division2 = country2
            country2 = "United Kingdom"

        if region2 not in data:
            continue

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
                        data[region2][country2][division2] = {}
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
def correct_data(data, type, corrections, add_annotations = True): #TODO: add region correction (e.g. for Turkey, Georgia)

    if type == "region":
        for (region, region_correct) in corrections:
            if region_correct not in data:
                data[region_correct] = {}
            for country in data[region]:
                if country not in data[region_correct]:
                    data[region_correct][country] = {}
                for division in data[region][country]:
                    if division not in data[region_correct][country]:
                        data[region_correct][country][division] = {}
                    for location in data[region][country][division]:
                        if location not in data[region_correct][country][division]:
                            data[region_correct][country][division][location] = []
                        for strain in data[region][country][division][location]:
                            if region != region_correct:
                                if add_annotations:
                                    additions_to_annotation.append(strain + "\tregion\t" + region_correct + " # previously " + region)
                            data[region_correct][country][division][location].append(strain)
            del data[region]

    if type == "country":
        for (region, country, region_correct, country_correct) in corrections:
            if country_correct not in data[region_correct]:
                data[region_correct][country_correct] = {}
            for division in data[region][country]:
                if division not in data[region_correct][country_correct]:
                    data[region_correct][country_correct][division] = {}
                for location in data[region][country][division]:
                    if location not in data[region_correct][country_correct][division]:
                        data[region_correct][country_correct][division][location] = []
                    for strain in data[region][country][division][location]:
                        if country != country_correct:
                            if add_annotations:
                                additions_to_annotation.append(strain + "\tcountry\t" + country_correct + " # previously " + country)
                        if region != region_correct:
                            if add_annotations:
                                additions_to_annotation.append(strain + "\tregion\t" + region_correct + " # previously " + region)
                        data[region_correct][country_correct][division][location].append(strain)
            del data[region][country]

    if type == "division":
        for (region, country, division, region_correct, country_correct, division_correct) in corrections:
            if country_correct not in data[region_correct]:
                data[region_correct][country_correct] = {}
            if division_correct not in data[region_correct][country_correct]:
                data[region_correct][country_correct][division_correct] = {}
            for location in data[region][country][division]:
                if location not in data[region_correct][country_correct][division_correct]:
                    data[region_correct][country_correct][division_correct][location] = []
                for strain in data[region][country][division][location]:
                    if division != division_correct:
                        if add_annotations:
                            additions_to_annotation.append(strain + "\tdivision\t" + division_correct + " # previously " + division)
                    if country != country_correct:
                        if add_annotations:
                            additions_to_annotation.append(strain + "\tcountry\t" + country_correct + " # previously " + country)
                    if region != region_correct:
                        if add_annotations:
                            additions_to_annotation.append(strain + "\tregion\t" + region_correct + " # previously " + region)
                    data[region_correct][country_correct][division_correct][location].append(strain)
            del data[region][country][division]

    if type == "location":
        for (region, country, division, location, region_correct, country_correct, division_correct, location_correct) in corrections:
            if country_correct not in data[region_correct]:
                data[region_correct][country_correct] = {}
            if division_correct not in data[region_correct][country_correct]:
                data[region_correct][country_correct][division_correct] = {}
            if location_correct not in data[region_correct][country_correct][division_correct]:
                data[region_correct][country_correct][division_correct][location_correct] = []
            for strain in data[region][country][division][location]:
                if location != location_correct:
                    if add_annotations:
                        additions_to_annotation.append(strain + "\tlocation\t" + location_correct + " # previously " + location)
                if division != division_correct:
                    if add_annotations:
                        additions_to_annotation.append(strain + "\tdivision\t" + division_correct + " # previously " + division)
                if country != country_correct:
                    if add_annotations:
                        additions_to_annotation.append(strain + "\tcountry\t" + country_correct + " # previously " + country)
                if region != region_correct:
                    if add_annotations:
                        additions_to_annotation.append(strain + "\tregion\t" + region_correct + " # previously " + region)
                data[region_correct][country_correct][division_correct][location_correct].append(strain)
            del data[region][country][division][location]
            if data[region][country][division] == {}:
                del data[region][country][division]
            if data[region][country] == {}:
                del data[region][country]

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
                    if add_annotations:
                        additions_to_annotation.append(strain + "\tdivision\t" + division + " # previously false division " + location)
                        additions_to_annotation.append(strain + "\tlocation\t" + location)
                    data[region][country][division][location].append(strain)
            del data[region][country][location]

    return data

# Search the ordering file for a similar name as the one given, and return it if the score is above a fixed threshold
def check_similar(ordering, name, type):
    diff_max = 0
    name_max = ""
    for name0 in ordering:
        diff = SequenceMatcher(None, name, name0).ratio()
        if name0 in name or name in name0:
            diff = 0.8
        if diff > diff_max:
            diff_max = diff
            name_max = name0

    if diff_max > 0.6:
        return name_max
    return ""


################################################################################
# Step 2: Clean up data
################################################################################


##### Step 2.0:
def adjust_to_database(data): #TODO: temporary solution, needs reworking
    for region in data:
        for country in data[region]:
            if country + ".txt" in listdir(path_to_config_files + "country_ordering/"):

                variants = {}
                with open(path_to_config_files + "country_ordering/" + country + "_variants.txt") as myfile:
                    country_variants = myfile.readlines()
                for line in country_variants:
                    if line == "\n":
                        continue
                    l = line.strip().split("\t")
                    variants[l[0]] = l[1]

                with open(path_to_config_files + "country_ordering/" + country + ".txt") as myfile:
                    country_ordering = myfile.readlines()

                arrondissement_to_location = {}
                location_to_arrondissement = {}
                provinces = []
                duplicates = {}

                for line in country_ordering:
                    if line == "\n" or "------" in line:
                        continue
                    if line.startswith("### "):
                        province = clean_string(line.strip()[4:])
                        provinces.append(province)
                        continue
                    if line.startswith("# "):
                        arrondissement = line.strip()[2:]
                        arrondissement_to_location[clean_string(arrondissement)] = []
                        continue

                    location = line.strip()
                    if location not in arrondissement_to_location[clean_string(arrondissement)]:
                        arrondissement_to_location[clean_string(arrondissement)].append(location)
                    if clean_string(location) in location_to_arrondissement:
                        if location_to_arrondissement[clean_string(location)] != arrondissement:
                            duplicates[clean_string(location)] = (arrondissement, location_to_arrondissement[clean_string(location)])
                    location_to_arrondissement[clean_string(location)] = arrondissement

                for division in data[region][country]:
                    
                    for location in data[region][country][division]:
                        division_c = clean_string(division)

                        if division == country:
                            continue
                        
                        if division_c in provinces and location == "":
                            continue

                        # division appears two times in country ordering - advise to pick one
                        if division_c in duplicates:
                            print("Attention duplicate: " + bold(division) + " found in " + bold(duplicates[division_c][0]) + " and " + bold(duplicates[division_c][1]))
                            print("Suggestion: check additional info for zip code")
                            print("Suggestion: check additional info for zip code")
                            continue

                        ### location given
                        if location != "":

                            # Ignore location, only check for correct division
                            if division_c in arrondissement_to_location:
                                continue

                            '''
                            # consistent with dataset
                            if clean_string(location) in location_to_arrondissement and division == location_to_arrondissement[clean_string(location)]:
                                continue

                            # other way around (in case of duplicates overwriting each other in location_to_arrondissement)
                            if division_c in arrondissement_to_location and location in arrondissement_to_location[division_c]:
                                continue

                            # location given, but with wrong division - adjust to correct division
                            if clean_string(location) in location_to_arrondissement and division != location_to_arrondissement[clean_string(location)]:
                                print("Wrong division " + bold(division) + " given for location " + bold(location))
                                print("Suggestion: add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, location_to_arrondissement[clean_string(location)], location]) + "] to manual_adjustments.txt")
                                continue

                            # location given, but with wrong spelling. Division is correct - adjust to correct location
                            if location in variants and clean_string(variants[location]) in location_to_arrondissement and division == location_to_arrondissement[clean_string(variants[location])]:
                                print("Location " + bold(location) + " should be adjusted to " + bold(variants[location]))
                                print("Suggestion: add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, variants[location]]) + "] to manual_adjustments.txt")
                                continue

                            # location given, but with wrong spelling. Division false - adjust both location and division
                            if location in variants and clean_string(variants[location]) in location_to_arrondissement and division != location_to_arrondissement[clean_string(variants[location])]:
                                print("Location " + bold(location) + " should be adjusted to " + bold(variants[location]) + ". Wrong division " + bold(division) + " given for location " + bold(variants[location]))
                                print("Suggestion: add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, location_to_arrondissement[clean_string(variants[location])], variants[location]]) + "] to manual_adjustments.txt")
                                continue
                            '''


                        ### location empty
                        else:
                            # given division is proper - no changes necessary
                            if division_c in arrondissement_to_location:
                                continue

                            # given division is proper, but misspelled - adjust spelling
                            if division in variants and (clean_string(variants[division]) in provinces or clean_string(variants[division]) in arrondissement_to_location):
                                print("Division " + bold(division) + " should be adjusted to " + bold(variants[division]))
                                print("Suggestion: add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, variants[division], location]) + "] to manual_adjustments.txt")
                                continue

                            # given division is actually a location
                            if division_c in location_to_arrondissement:
                                print("Given division " + bold(division) + " is actually a location within division " + bold(location_to_arrondissement[division_c]))
                                print("Suggestion: add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, location_to_arrondissement[division_c], division]) + "] to manual_adjustments.txt")
                                continue

                            # given division is misspelled and location
                            if division in variants and clean_string(variants[division]) in location_to_arrondissement:
                                print("Given division " + bold(division) + " is a misspelled location " + bold(variants[division]) + " within division " + bold(location_to_arrondissement[clean_string(variants[division])]))
                                print("Suggestion: add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, location_to_arrondissement[clean_string(variants[division])], variants[division]]) + "] to manual_adjustments.txt")
                                continue


                        print("Missing combination in " + country + " database: " + bold(division + ", " + location))


    print("\n=============================\n")
    return data


##### Step 2.05: Apply manual adjustments set in manual_adjustments.txt
def manual_adjustments(data):
    manual_adjustments = read_local_file("manual_adjustments.txt")

    for g in manual_adjustments:
        seqs_to_correct = []
        (region_before, country_before, division_before, location_before) = g.split("/")
        (region_correct, country_correct, division_correct, location_correct) = manual_adjustments[g].split("/")
        for region in data:
            for country in data[region]:
                for division in data[region][country]:
                    for location in data[region][country][division]:
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

                        if (region == region_before or region_before == "*") and (country == country_before or country_before == "*") and (division == division_before or division_before == "*") and (location == location_before or location_before == "*"):
                            seqs_to_correct.append((region, country, division, location, region2, country2, division2, location2))
                            print("Manual adjustment: " + bold("/".join([region, country, division, location])) + " -> " + bold("/".join([region2, country2, division2, location2])))

        data = correct_data(data, "location", seqs_to_correct)
    print("\n=============================\n")
    return data



##### Step 2.1: Apply all known variants stored in an external file variants.txt
def apply_variants(data):
    variants = read_local_file("variants.txt")

    countries_to_switch = []
    for region in data:
        for country in data[region]:
            if country in variants['country']:
                match_found = False
                # if the first entry has no specified hierarchy, all other entries of this place name are ignored
                if type(variants['country'][country][0]) is not tuple:
                    match_found = True
                    country_correct = variants['country'][country][0]
                else:
                    for country_option in variants['country'][country]:
                        if country_option[1] == "(" + region + ")":
                            match_found = True
                            country_correct = country_option[0]
                            break
                if match_found:
                    print("Apply variant (country): " + bold(country) + " -> " + bold(country_correct))
                    countries_to_switch.append((region, country, region, country_correct))

    data = correct_data(data, "country", countries_to_switch)

    divisions_to_switch = []
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                if division in variants['division']:
                    match_found = False
                    if type(variants['division'][division][0]) is not tuple:
                        match_found = True
                        division_correct = variants['division'][division][0]
                    else:
                        for division_option in variants['division'][division]:
                            if division_option[1] == "(" + region + ", " + country + ")":
                                match_found = True
                                division_correct = division_option[0]
                                break
                    if match_found:
                        print("Apply variant (division): " + bold(division) + " -> " + bold(division_correct))
                        divisions_to_switch.append((region, country, division, region, country, division_correct))

    data = correct_data(data, "division", divisions_to_switch)

    locations_to_switch = []
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                for location in data[region][country][division]:
                    if location in variants['location']:
                        match_found = False
                        if type(variants['location'][location][0]) is not tuple:
                            match_found = True
                            location_correct = variants['location'][location][0]
                        else:
                            for location_option in variants['location'][location]:
                                if location_option[1] == "(" + region + ", " + country + ", " + division + ")":
                                    match_found = True
                                    location_correct = location_option[0]
                                    break
                        if match_found:
                            print("Apply variant (location): " + bold(location) + " -> " + bold(location_correct))
                            locations_to_switch.append((region, country, division, location, region, country, division, location_correct))

    data = correct_data(data, "location", locations_to_switch)

    print("\n=============================\n")
    return data

def apply_typical_errors(data): #TODO: rename, maybe join with UK as region? also use correct_data()
    wrong_regions = read_local_file("wrong_regions.txt")

    countries_to_switch = []
    for country in wrong_regions:
        region_correct = wrong_regions[country]
        for region in data:
            if region == region_correct:
                continue
            if country in data[region]:
                print("Found incorrect region " + bold(region) + " for country " + bold(country) + " (correct region: " + bold(region_correct) + ")" )
                countries_to_switch.append((region, country, region_correct, country))

    data = correct_data(data, "country", countries_to_switch)


    print("\nAdjustments made to avoid international duplicates (e.g. cruise ships) for generation of color_ordering.tsv:\n")
    divisions_to_switch = []
    locations_to_switch = []
    international_exceptions = read_local_file("international_exceptions.txt")
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                if division in international_exceptions["division"]:
                    (region_correct, country_correct) = tuple(international_exceptions["division"][division][0].split(", "))
                    if region == region_correct and country == country_correct:
                        continue
                    print("division " + division + ": " + region + ", " + country + " => " + region_correct + ", " + country_correct)
                    divisions_to_switch.append((region, country, division, region_correct, country_correct, division))
                for location in data[region][country][division]:
                    if location in international_exceptions["location"]:
                        (region_correct, country_correct, division_correct) = tuple(international_exceptions["location"][location][0].split(", "))
                        if region_correct == region and country_correct == country and division_correct == division:
                            continue
                        print("location " + location + ": " + region + ", " + country + ", " + division + " => " + region_correct + ", " + country_correct + ", " + division_correct)
                        locations_to_switch.append((region, country, division, location, region_correct, country_correct, division_correct, location))
    data = correct_data(data, "division", divisions_to_switch, add_annotations = False) #Changes only needed for generation of color_ordering to avoid international duplicates, should stay in original metadata
    data = correct_data(data, "location", locations_to_switch, add_annotations = False)
    print()

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
                if division != "":
                    for location in data[region][country][division]:
                        if location != "":
                            if location in data[region][country] and location != division:
                                div_as_loc[location] = (region, country, division)
                                print("Unknown location found as division: " + bold(location) + " (true division: " + bold(division) + ")")
                                print("(Suggestion: add " + "[" + "/".join([region, country, location, ""]) + "\t" + "/".join([region, country, division, location]) + "]" + " to manual_adjustments.txt)")
                                if list(data[region][country][location]) != [""]:
                                    print("Attention: location(s) " + ", ".join(data[region][country][location]) + " would be lost.")

    print("\n=============================\n")

##### Step 2.3: Check for duplicate divisions/locations in different countries/divisions (known cases stored in duplicates.txt as well as checking for new cases)
def check_duplicate(data):

    #Check known duplicates
    # TODO: Only locations covered properly (divisions: only alert)
    duplicates = read_local_file("duplicates.txt")
    abbreviations = read_local_file("abbreviations.txt")

    duplicate_locations = []
    for region in data:
        for country in data[region]:
            for division in data[region][country]:
                for location in data[region][country][division]:
                    if location in duplicates:
                        print("Known duplicate detected: " + bold(location))
                        if abbreviations.get(division) is not None:
                            print("Please add [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, location + " " + abbreviations[division]]) + "] to manual_adjustments.txt")
                            location_correct = location + " " + abbreviations[division]
                            duplicate_locations.append((region, country, division, location, region, country, division, location_correct))
                        else:
                            print("No abbreviation for " + division + ", please add one to abbreviations.txt and rerun.")
    data = correct_data(data, "location", duplicate_locations)

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
        if len(division_to_country[division]) > 1:
            if not any(x in division for x in cruise_abbrev): #ignore cruise ship ones
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
                print("Suggestion: Add " + location + " to duplicates.txt")

            
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
    locations_skipped = {}
    n_skipped_locations = 0

    for region in data:
        data_clean[region] = {}

        for country in data[region]:

            if country not in ordering["country"] or country not in lat_longs["country"]:
                s = bold(country)
                if country not in ordering["country"] and country in lat_longs["country"]:
                    s = s + " (only missing in ordering => auto-added to color_ordering.tsv)"
                    data_clean[region][country] = {}
                else:
                    if country in ordering["country"] and country not in lat_longs["country"]:
                        s = s + " (only missing in lat_longs)"
                    else:
                        if country in ordering["division"] or country in lat_longs["division"]:
                            s = s + " (present as division)"

                missing["country"].append(s)
                if "(only missing in ordering" not in s:
                    clean_missing["country"].append(country)

            else:
                data_clean[region][country] = {}


            for division in data[region][country]:
                if division == "":
                    continue

                if division not in ordering["division"] or division not in lat_longs["division"]:
                    s = bold(division)
                    name0 = ""
                    if country in hierarchical_ordering.get(region, ""):
                        name0 = check_similar(hierarchical_ordering[region][country], division, "division")
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
                                s += " (similar name in same country: " + bold(name0) + " - consider adding " + "[" + "/".join([region, country, division, "*"]) + "\t" + "/".join([region, country, name0, "*"]) + "]" + " to manual_adjustments.txt)"
                            if division in ordering["location"] or division in lat_longs["location"]:
                                s = s + " (present as location)"
                    if country not in missing["division"]:
                        missing["division"][country] = []
                        clean_missing["division"][country] = []
                    missing["division"][country].append(s)
                    if "(only missing in ordering" not in s:
                        clean_missing["division"][country].append(division)

                else:
                    if country not in data_clean[region]:
                        print("Conflict: division " + division + " should be added to color_ordering.tsv, but country " + country + " is missing from dataset")
                    else:
                        data_clean[region][country][division] = []

                for location in data[region][country][division]:
                    if location == "":
                        continue

                    if region == "North America":
                        if location not in ordering["location"] or location not in lat_longs["location"]:
                            s = bold(location)
                            name0 = check_similar(hierarchical_ordering[region][country][division], location, "location") if hierarchical_ordering[region].get(country) is not None and hierarchical_ordering[region][country].get(division) else ""
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
                                    s += " (similar name in same division: " + bold(name0) + " - consider adding " + "[" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, name0]) + "]" + " to manual_adjustments.txt)"
                                if location in ordering["location"] and location not in lat_longs["location"]:
                                    s = s + " (only missing in lat_longs)"
                                if location in ordering["division"] or location in lat_longs["division"]:
                                    s = s + " (present as division)"
                                if country == "USA" and "County" not in location:
                                    ct = "County"
                                    if division == "Louisiana":
                                        ct = "Parish"
                                    s = s + " (correction to " + ct + " might be necessary using [" + "/".join([region, country, division, location]) + "\t" + "/".join([region, country, division, location + " " + ct]) + "]"

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
                                if "(only missing in ordering" not in s:
                                    clean_missing["location"][country][division].append(location)
                            else:
                                print("Cruise-associated location ignored ("+location+")")
                        else:
                            if country not in data_clean[region]:
                                print(
                                            "Conflict: location " + location + " should be added to color_ordering.tsv, but country " + country + " is missing from dataset")
                            else:
                                if division not in data_clean[region][country]:
                                    if not any(x in location for x in cruise_abbrev) and not any(x in division for x in cruise_abbrev):
                                        print("Conflict: location " + location + " should be added to color_ordering.tsv, but division " + division + " is missing from dataset")
                                else:
                                    data_clean[region][country][division].append(location)
                    else:
                        if country in data_clean[region]:
                            if division in data_clean[region][country]:
                                data_clean[region][country][division].append(location)
                        if location not in lat_longs["location"] and location not in ordering["location"]: #only completely new locations are considered for the counting
                            if region not in locations_skipped:
                                locations_skipped[region] = {}
                            if country not in locations_skipped [region]:
                                locations_skipped[region][country] = {}
                            if division not in locations_skipped[region][country]:
                                locations_skipped[region][country][division] = []
                            locations_skipped[region][country][division].append(location)
                            n_skipped_locations += 1

    print("Number of non-North American locations skipped: " + str(n_skipped_locations))


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
            print("country\t" + country)
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

        with open(path_to_output_files+"new_lat-longs.tsv", 'w') as out:
            out.write("\n".join(new_lat_longs))
        print("New lat-longs written out to "+path_to_output_files+"new_lat-longs.tsv")

        answer = input("Would you like to use auto-sort for these lat_longs? y or n")
        if answer == "y":
            auto_add_lat_longs(new_lat_longs)


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
    tries = 0
    while redo == True:
        if tries < 5:
            try:
                new_place = ask_geocoder(typed_place, geolocator)
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
                if level.lower() in new_place_string.lower():
                    new_place_string = bold(level).join(new_place_string.split(level))
                    full_place_string = bold(level).join(full_place_string.split(level))

            print("\nCurrent place for missing {}:\t".format(geo_level) + full_place_string)

            print("Geopy suggestion: "+ new_place_string)
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

    print(answer)
    return answer

################################################################################
# Step 3: Storage of locations, divisions etc hierarchical manner
################################################################################

# Util needed to sort a given list of locations or divisions by their coordinated stored in lat_longs
# TODO: enable for countries as well
def sort_by_coordinates(data, coordinates):
    max_lat = -90
    min_lat = 90
    max_long = -150
    min_long = 150
    for hierarchy in data:
        (lat, long) = coordinates[hierarchy]
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

    with open(path_to_output_files+"color_ordering.tsv", mode) as out:
        if hierarchy not in ["region", "country", "division", "location"]:
            for l in data[hierarchy]:
                out.write(hierarchy + "\t" + l + "\n")
            out.write("\n################\n\n\n")
            return

        # Give fixed order of regions to retain the usual coloring order
        region_order = ["Asia",
                        "Oceania",
                        "Africa",
                        "Europe",
                        "South America",
                        "North America"]

        for region in region_order:

            if hierarchy == "region":
                out.write("region\t" + region + "\n")
                continue

            out.write("\n# " + region + "\n")
            for country in sort_by_coordinates(data[region], lat_longs["country"]): #TODO: would be nice to sort this by coordinate too, but would need to add most lat_longs first!

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

                    for location in sorted(data[region][country][division]):
                        out.write("location\t" + location + "\n")

            if hierarchy == "location" or hierarchy == "division":
                out.write("\n################\n")

        out.write("\n################\n\n\n")


def auto_add_annotations(additions_to_annotation):
    enable_duplicate_check = True

    with open("../ncov-ingest/source-data/gisaid_annotations.tsv") as myfile:
        annotations = myfile.readlines()
    types = {"geography": ["location", "division", "country", "region", "division_exposure", "country_exposure", "region_exposure"], "special": ["sampling_strategy", "date", "host", "strain"], "paper": ["title", "paper_url"], "genbank": ["genbank_accession"]}
    sections = {"comments": [], "geography": [], "special": [], "paper": [], "genbank": []}

    print("The following annotations have unknown type:")
    for list in [annotations, additions_to_annotation]:
        for line in list:
            if not line.endswith("\n"):
                line = line + "\n"
            if line.startswith("#"):
                sections["comments"].append(line)
                continue
            t1 = line.split("\t")[2]
            type_found = False
            for t in types:
                if t1 in types[t]:
                    if line not in sections[t]:
                        sections[t].append(line)
                    type_found = True
                    break
            if not type_found:
                print(line)

    with open(path_to_output_files + "gisaid_annotations.tsv", "w") as out:
        for t in sections:
            for l in sorted(sections[t]):
                out.write(l)
    print("New annotation auto-added to " + path_to_output_files + "gisaid_annotations.tsv")

    if enable_duplicate_check:
        print("The following annotations have duplicate annotations:")
        duplicate_check = {}
        with open(path_to_output_files + "gisaid_annotations.tsv") as list:
            for line in list:
                if not line.startswith("#"):
                    t = line.split("\t")[2]
                    if t not in types["paper"] and t not in types["genbank"]:
                        epi = line.split("\t")[1]
                        if t not in duplicate_check:
                            duplicate_check[t] = []
                        if epi not in duplicate_check[t]:
                            duplicate_check[t].append(epi)
                        else:
                            print("Attention: Duplicate annotation for " + epi + ", " + t)






if __name__ == '__main__':

    ################################################################################
    # Step 0: Read data
    ################################################################################

    # Read current metadata
    #path_to_ncov = "../../" # TODO: adjust file structure properly
    with open("data/downloaded_gisaid.tsv") as myfile:
        metadata_gisaid = myfile.readlines()

    with open("data/metadata_genbank.tsv") as myfile:
        metadata_genbank = myfile.readlines()


    # Read orderings and lat_longs
    ordering, ordering_other = read_geography_file("defaults/color_ordering.tsv") #TODO: combine with read_local_files()?
    hierarchical_ordering = read_geography_file("defaults/color_ordering.tsv", True)
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
    data = {}
    data = read_metadata(metadata_gisaid, data, "gisaid")
    data = read_metadata(metadata_genbank, data, "genbank")



    ##### Step 1.2: Collection of regions, countries and divisions of exposure
    # In case some geographic units are only found in the exposure information of the metadata, iterate again over the metadata and add to the dataset
    # Since travel history related entries are prone to errors, check for each entry whether it collides with already existing data.

    # TODO: Currently commented out due to numerous inconsistencies
    data = read_exposure(data, metadata_gisaid)


    ################################################################################
    # Step 2: Clean up data
    ################################################################################

    ##### Step 2.0: Before checking for any of the other adjustments, apply manually designed adjustments stored in manual_adjustments.txt
    # This includes manually setting the region, country, division and location before and after the adjustment
    data = manual_adjustments(data)

    ##### Step 2.05: Adjust the divisions and locations by comparing them to a known database - only accessible for Belgium at the moment
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
    write_ordering(data, "region")
    for type in ordering_other:
        write_ordering(ordering_other, type)

    ##### Bonus step: Print out all collected annotations - if considered correct, they can be copied by the user to annotations.tsv
    with open(path_to_output_files+"new_annotations.tsv", 'w') as out:
        out.write("\n".join(sorted(additions_to_annotation)))
    print("New annotation additions written out to " + path_to_output_files + "new_annotations.tsv")

    auto_add_annotations(additions_to_annotation)

    # Only print line if not yet present
    # Print warning if this GISAID ID is already in the file
    lines_exclude = ["title", "authors", "paper_url", "genbank_accession", "sampling_strategy"]
    annot_lines_to_write = []
    for line in additions_to_annotation:
        if line in annotations:
            continue
        #print(line)
        if len(line.split("\t")) == 4:
            epi = line.split("\t")[1]
            if epi in annotations:
                number_of_occurences = annotations.count(line.split("\t")[1])
                irrelevant_occurences = sum([(line.split("\t")[1] + "\t" + s) in annotations for s in lines_exclude])
                if number_of_occurences > irrelevant_occurences:
                    for l in annotations.split("\n"):
                        if epi in l:
                            if not l.startswith("#"):
                                print("Warning: " + epi + " already exists in annotations! (" + bold(line.split("\t")[2]) + " " + line.split("\t")[3] + " vs " + bold(l.split("\t")[2]) + " " + l.split("\t")[3] + ")")

