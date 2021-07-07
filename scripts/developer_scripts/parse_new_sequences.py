import os
import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import random
import numpy as np
import math
import datetime
register_matplotlib_converters()



def bold(s):
    return('\033[1m' + s + '\033[0m')

def strike(s):
    return('\u0336'.join(s) + '\u0336')

# Cut a string of the format "key: content" into a tuple (key, content)
def cut(s):
    key = s.split(":")[0]
    content = ":".join(s.split(":")[1:])[1:]
    return (key, content)

# Read all files within given directory that start with "metadata-changes" in a sorted manner. Contents are stored in a
# dictionary first by GISAID epi isl, then by info-type (e.g. date, division, host)
def read_data(path):

    data = {} # added and changed sequences
    list_of_strains = []

    for file in sorted(os.listdir(path)):
        if file == '.DS_Store':
            continue

        if file.startswith("metadata-changes"):
            with open(path + file) as f:
                metadata_changes = f.readlines()

            added = False
            changed = False
            for i in range(len(metadata_changes)):
                k = metadata_changes[i].strip()

                if k.endswith("sequences added") or k.endswith("sequence added"):
                    added = True
                    changed = False
                if k.endswith("sequences changed") or k.endswith("sequence changed"):
                    added = False
                    changed = True
                if k.endswith("sequences removed") or k.endswith("sequence removed"):
                    added = False
                    changed = False

                if k.startswith("gisaid_epi_isl"):
                    (key, id) = cut(k)
                    if added:
                        list_of_strains.append(id)
                        if id in data:
                            print("Attention, same sequence added two times! (" + id + ")")
                        data[id] = {}
                        j = -2
                        while (i+j) < len(metadata_changes) and metadata_changes[i+j] != "\n":
                            k = metadata_changes[i+j].strip()
                            (key, content) = cut(k)
                            if key != "gisaid_epi_isl":
                                data[id][key] = content
                            j += 1

                    elif changed:
                        if id in data:
                            j = 1
                            k = metadata_changes[i+j]
                            while k != "\n":
                                (key, content) = cut(k.strip())
                                data[id][key] = content.split("=>")[1].strip().strip("\"") #overwrite with latest changes
                                j += 1
                                if (i+j) >= len(metadata_changes):
                                    break
                                k = metadata_changes[i+j]

        if file.startswith("metadata-additions"):
            with open(path + file) as f:
                metadata_additions = f.readlines()
            header = metadata_additions[0].strip().split("\t")
            for line in metadata_additions[1:]:
                l = line.strip().split("\t")
                id = l[2]
                if id in data:
                    print("Attention, same sequence added two times! (" + id + ")")
                data[id] = {}
                for i in range(len(header)):
                    if header != "gisaid_epi_isl":
                        data[id][header[i]] = l[i]

    return (data, list_of_strains)

def check_for_recency(counts, list_of_strains, lab_collection, path_to_metadata, table_file_name):
    print("\n----------------------------------------------\n")

    countries = {}
    subm_labs = {}
    origlab_authors = {}
    cutoff_date = datetime.datetime.combine(datetime.date.today(), datetime.datetime.min.time()) - datetime.timedelta(days=30) #last XX days


    with open(path_to_metadata + "downloaded_gisaid.tsv") as f:
        header = f.readline().split("\t")
        country_i = header.index("country")
        subm_date_i = header.index("date_submitted")
        strain_i = header.index("gisaid_epi_isl")
        subm_lab_i = header.index("submitting_lab")
        orig_lab_i = header.index("originating_lab")
        author_i = header.index("authors")


        print("Collecting all labs from the last month from metadata... (this may take a while)")

        line = f.readline()
        while line:
            l = line.split("\t")
            country = l[country_i]
            lab = l[subm_lab_i]
            orig_lab = l[orig_lab_i]
            author = l[author_i]
            if country == "United Kingdom":
                line = f.readline()
                continue
            if country in counts:
                if country in subm_labs and lab in subm_labs[country] and subm_labs[country][lab] > 20 and orig_lab in origlab_authors[country] and origlab_authors[country][orig_lab] > 20 and author in origlab_authors[country] and origlab_authors[country][author] > 20:
                    line = f.readline()
                    continue
                date = datetime.datetime.fromisoformat(l[subm_date_i])
                if date > cutoff_date:
                    strain = l[strain_i]
                    if strain not in list_of_strains:
                        if country not in countries:
                            print(country)
                            countries[country] = 0
                        countries[country] += 1
                        if country not in subm_labs:
                            subm_labs[country] = {}
                            origlab_authors[country] = {}
                        if lab not in subm_labs[country]:
                            subm_labs[country][lab] = 0
                        subm_labs[country][lab] += 1
                        if author not in origlab_authors[country]:
                            origlab_authors[country][author] = 0
                        origlab_authors[country][author] += 1
                        if orig_lab not in origlab_authors[country]:
                            origlab_authors[country][orig_lab] = 0
                        origlab_authors[country][orig_lab] += 1
            line = f.readline()

    print("\nSearching for twitter handles... ")

    rare_countries = [
        c
        for c in counts
        if c != "United Kingdom" and (c not in countries or countries[c] <= 20)
    ]

    lab_dictionary = read_excel_lab_file(table_file_name)
    lab_collection_present = {}
    for country, value in subm_labs.items():
        if country not in lab_collection_present:
            lab_collection_present[country] = {}

        for lab in value:
            n = subm_labs[country][lab]
            if country in lab_dictionary and lab.lower() in lab_dictionary[country]:
                k = lab_dictionary[country][lab.lower()]
                for l in k.split(", "):
                    if l not in lab_collection_present[country]:
                        lab_collection_present[country][l] = 0
                    lab_collection_present[country][l] += n
            else:
                print("Lab " + bold(lab) + " (" + country + ") not found and will be excluded from future steps. Please fill into excel table to include it.")
        for lab in origlab_authors[country]:
            n = origlab_authors[country][lab]
            if country in lab_dictionary and lab.lower() in lab_dictionary[country]:
                k = lab_dictionary[country][lab.lower()]
                for l in k.split(", "):
                    if l not in lab_collection_present[country]:
                        if l not in lab_collection_present[country]:
                            lab_collection_present[country][l] = 0
                        lab_collection_present[country][l] += n

    lab_collection_present["United Kingdom"] = {"@CovidGenomicsUK": 1000}

    rare_labs = {}
    for region in lab_collection:
        for country in lab_collection[region]:
            for lab in lab_collection[region][country]:
                if country not in lab_collection_present or lab not in lab_collection_present[country] or country in rare_countries:
                    if region not in rare_labs:
                        rare_labs[region] = {}
                    if country not in rare_labs[region]:
                        rare_labs[region][country] = []
                    rare_labs[region][country].append(lab)

    print("\nCountries that have submitted < 20 sequences last month (all of these will be included in the tweet):")
    print(rare_countries)

    print("\nSubmitters that have not submitted last month (all of these will be included in the tweet):")
    print(rare_labs)

    return rare_countries, rare_labs


# Double check sample dates for invalid format or unrealistic / impossible date
def check_dates(data, today):
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

    invalid_sample_date = {}
    suspicious_sample_date = {}

    for id in list(data.keys()):
        date = data[id]["date"]
        strain = data[id]["strain"]
        country = data[id]["country"]

        if len(date) != len(today):
            invalid_sample_date[strain] = (date, country)
            data.pop(id)
            continue

        day = int(date[8:])
        month = int(date[5:7])
        year = int(date[:4])

        day_today = int(today[8:])
        month_today = int(today[5:7])
        year_today = int(today[:4])


        #check for past dates
        if year < 2019 or ((year) == 2019 and month < 12):
            invalid_sample_date[strain] = (date, country)
            data.pop(id)
            continue

        # Check for future dates
        if (year > year_today) or (year == year_today and month > month_today) or (year == year_today and month == month_today and day > day_today):
            invalid_sample_date[strain] = (date, country)
            data.pop(id)
            continue

        clade = data[id]["Nextstrain_clade"]
        dev = data[id]["clock_deviation"]
        if clade == "":
            print("Clade missing for sequence " + id)
        elif clade not in clade_dates:
            print("Unknown clade " + clade + " for sequence " + id)
        else:
            clade_day = clade_dates[clade]
            day_clade = int(clade_day[8:])
            month_clade = int(clade_day[5:7])
            year_clade = int(clade_day[:4])

            if (year < year_clade) or (year == year_clade and month < month_clade) or (year == year_clade and month == month_clade and day < day_clade):
                suspicious_sample_date[strain] = date + " (" + clade + ", clock deviation = " + dev + ")"
                data.pop(id)
                continue


    invalid_dates_by_country = {}
    for strain in invalid_sample_date:
        (date, country) = invalid_sample_date[strain]
        if country not in invalid_dates_by_country:
            invalid_dates_by_country[country] = {}
        if date not in invalid_dates_by_country[country]:
            invalid_dates_by_country[country][date] = 0
        invalid_dates_by_country[country][date] += 1



    print("\n----------------------------------------------\n")
    print("Invalid sample dates (automatically excluded from total counts):")
    for country, value_ in invalid_dates_by_country.items():
        print(country)
        for date in invalid_dates_by_country[country]:
            print(date + " (" + str(value_[date]) + ")")
        print("")

    print("\nSample date before clade (automatically excluded from total counts):")
    for strain, value in suspicious_sample_date.items():
        print(strain + ": " + value)

    return data


flagged_properties = {"originating_lab": ["Synlab Haut de France"]}

# Check for certain unique properties and potentially exclude (e.g. all sequences from a certain submission lab)
def check_flagged_properties(data):

    flagged_strains = {
        p: {name: [] for name in flagged_properties[p]}
        for p in flagged_properties
    }

    seqs_found = False
    for id in list(data.keys()):
        strain = data[id]["strain"]
        exclude = False
        for p in flagged_properties:
            prop = data[id][p]
            for name in flagged_properties[p]:
                if prop == name:
                    flagged_strains[p][name].append(strain)
                    seqs_found = True
                    exclude = True
        if exclude:
            data.pop(id)

    if seqs_found:
        print(bold("\nFlagged properties found! Please check outputs_new_sequences/sequences_exclude.txt for strain names to exclude") + " (automatically excluded from total counts).\n")
    else:
        print("\nNo flagged properties found.\n")

    with open(path_to_outputs + "sequences_exclude.txt", "w") as out:
        out.write("\n\nStrains to add to exclude (based on flagged properties):\n")
        for p, value in flagged_strains.items():
            for name in flagged_properties[p]:
                out.write(p + " = \"" + name + "\":\n")
                for strain in value[name]:
                    out.write(strain + "\n")
                out.write("\n")

    return data


# Plot distribution of dates to detect suspicious behaviour, e.g. large bunch of sequences from the same date
def plot_dates(data, path):
    dates_by_country = {}

    for id in data:
        country = data[id]["country"]
        date = datetime.datetime.strptime(data[id]["date"], '%Y-%m-%d')
        if country not in dates_by_country:
            dates_by_country[country] = {}
        if date not in dates_by_country[country]:
            dates_by_country[country][date] = 0
        dates_by_country[country][date] += 1

    #remove old plots
    for f in os.listdir(path):
        if f.startswith("dates_"):
            os.remove(path + f)

    for country in dates_by_country:
        dates = list(dates_by_country[country].keys())
        values = list(dates_by_country[country].values())
        plt.figure()
        plt.bar(dates, values)
        plt.title(country)
        plt.xticks(rotation=45, ha="right", size = 7)
        plt.savefig(path + "dates_" + country.replace(".", ""))
        plt.close()


# Count for every country the number of new sequences
def print_counts(data):
    counts = {}
    for id in data:
        country = data[id]["country"]
        division = data[id]["division"]
        if country not in counts:
            counts[country] = {}
        if division not in counts[country]:
            counts[country][division] = 0
        counts[country][division] += 1

    sum_total = 0
    for country, value_ in counts.items():
        sum_country = sum(value_[division] for division in counts[country])
        sum_total += sum_country

    print("\n----------------------------------------------\n")
    print("Total counts: " + str(sum_total))

    with open(path_to_outputs + "tweet_resources.txt", "w") as out:
        for country, value in counts.items():
            s = country + ": "
            sum_country = 0
            for division in counts[country]:
                sum_country += value[division]
            s = s + str(sum_country)
            if len(counts[country]) == 1:
                s = s + " (" + division + ")"
            else:
                s = (
                    s
                    + " ("
                    + ", ".join(
                        str(counts[country][division]) + " " + division
                        for division in counts[country]
                    )
                    + ")"
                )

            print(s)
            out.write(s + "\n")
        out.write("\n\n\n")
    return counts

def read_excel_lab_file(table_file_name):
    excel_table = pd.read_excel(table_file_name, index_col=0, skiprows=1)
    excel_table = excel_table.fillna("empty?")
    lab_dictionary = {}
    for country, row in excel_table.iterrows():
        description = row["Who"]
        handle = row["Who to tag"]
        if country not in lab_dictionary:
            lab_dictionary[country] = {}
        if description in lab_dictionary[country]:
            print("Warning: lab description is found two times in excel table in same country (" + str(
                country) + ", " + str(description) + ")")
        lab_dictionary[country][description.lower()] = handle
    return lab_dictionary

# Collect all submitting and originating labs as well as authors and try to infer as many twitter handles as possible
# from a given excel file
def collect_labs(data, table_file_name):
    print("\n----------------------------------------------")
    submitting_labs = {}
    originating_labs = {}
    authors = {}
    for id in data:
        region = data[id]["region"]
        country = data[id]["country"]
        submitting_lab = data[id]["submitting_lab"]
        originating_lab = data[id]["originating_lab"]
        author = data[id]["authors"]

        _extracted_from_collect_labs_13(
            region, submitting_labs, country, submitting_lab
        )

        if region not in originating_labs:
            originating_labs[region] = {}
        if country not in originating_labs[region]:
            originating_labs[region][country] = []
        if originating_lab not in originating_labs[region][country] and originating_lab != submitting_lab:
            originating_labs[region][country].append(originating_lab)

        _extracted_from_collect_labs_13(region, authors, country, author)
    lab_dictionary = read_excel_lab_file(table_file_name)
    lab_UK = lab_dictionary["United Kingdom"]["COVID-19 Genomics UK Consortium".lower()]
    lab_collection = {}

    print("\nSubmitting labs:\n(Note: small differences in spelling might cause lab to not be identified. Consider adjusting the spelling in the spreadsheet!)\n")
    for region, value_ in submitting_labs.items():
        if region not in lab_collection:
            lab_collection[region] = {}
        for country in sorted(submitting_labs[region]):
            if country not in lab_collection[region]:
                lab_collection[region][country] = []

            s = country + ":\n"
            for lab in value_[country]:
                s += lab + ": "
                if country in lab_dictionary and lab.lower() in lab_dictionary[country]:
                    k = lab_dictionary[country][lab.lower()]
                    for l in lab_dictionary[country][lab.lower()].split(", "):
                        if l not in lab_collection[region][country]:
                            lab_collection[region][country].append(l)
                else:
                    k = "?"
                    lab_collection[region][country].append("???")
                if country == "United Kingdom":
                    k = strike(k) + " " + lab_UK

                s += bold(k) + "\n"
            print(s)

    print("----------------------------------------------\n")
    print("Originating labs (only printed if found in excel sheet):\n")
    for region, value__ in originating_labs.items():
        for country in originating_labs[region]:
            s = country + ":\n"
            for lab in value__[country]:
                if country in lab_dictionary and lab.lower() in lab_dictionary[country]:
                    s += lab
                    s += ": "
                    k = lab_dictionary[country][lab.lower()]
                    if country == "United Kingdom":
                        k = strike(k) + " " + lab_UK
                    s += bold(k)
                    for l in lab_dictionary[country][lab.lower()].split(", "):
                        if l not in lab_collection[region][country]:
                            lab_collection[region][country].append(l)
                    s += "\n"
            if s != country + ":\n":
                print(s)


    print("----------------------------------------------\n")
    print("Authors (only printed if found in excel sheet):\n")
    for region, value in authors.items():
        for country in value:
            s = country + ":\n"
            for author in authors[region][country]:
                if country in lab_dictionary and author.lower() in lab_dictionary[country]:
                    s += author
                    s += ": "
                    k = lab_dictionary[country][author.lower()]
                    if country == "United Kingdom":
                        k = strike(k) + " " + lab_UK
                    s += bold(k)
                    for a in lab_dictionary[country][author.lower()].split(", "):
                        if a not in lab_collection[region][country]:
                            lab_collection[region][country].append(a)
                    s += "\n"
            if s != country + ":\n":
                print(s)


    if (
        "Europe" in lab_collection
        and "United Kingdom" in lab_collection["Europe"]
    ):
        lab_collection["Europe"]["United Kingdom"] = [lab_UK]

    return lab_collection

def _extracted_from_collect_labs_13(region, arg1, country, arg3):
    if region not in arg1:
        arg1[region] = {}
    if country not in arg1[region]:
        arg1[region][country] = []
    if arg3 not in arg1[region][country]:
        arg1[region][country].append(arg3)




# Produce a concise overview of the new sequences containing only strain name, sampling and submission date sorted by
# country, useful for compiling tweets and collecting screenshots
def overview_with_dates(data, file_name):

    data_sorted = {}
    for id in data:
        strain = data[id]["strain"]
        submission_date = data[id]["date_submitted"]
        samlpe_date = data[id]["date"]
        country = data[id]["country"]

        if country not in data_sorted:
            data_sorted[country] = []
        data_sorted[country].append(strain + "\t" + samlpe_date + "\t" + submission_date)

    with open(file_name, "w") as myfile:
        myfile.write("strain\tsampling date\tsubmission date\n")
        for country, value in data_sorted.items():
            myfile.write(country + "\n")
            for s in value:
                myfile.write(s + "\n")
            myfile.write("\n")


# Produce output file containing all sequences from a certain region after a certain date
def filter_for_date_region(data, path_to_outputs, params):
    (region, month) = params

    special_strains = {}
    for id in data:
        date = data[id]["date"]
        if int(date[5:7]) >= month and region == data[id]["region"]:
            country = data[id]["country"]
            if country not in special_strains:
                special_strains[country] = {}
            if date[:7] not in special_strains[country]:
                special_strains[country][date[:7]] = 0
            special_strains[country][date[:7]] += 1

    with open(path_to_outputs + "special_check_" + region + "_" + str(month) + ".txt", "w") as myfile:
        myfile.write("New sequences from " + region + " after month " + str(month) + "\n\n")
        for country, value in special_strains.items():
            myfile.write(country + "\n")
            for date in sorted(special_strains[country]):
                myfile.write(date + ": " + str(value[date]) + "\n")
            myfile.write("\n")

def prepare_tweet(counts, total_lab_collection, lab_collection):

    links = {
        "Africa": "nextstrain.org/ncov/africa",
        "Asia": "nextstrain.org/ncov/asia",
        "Europe": "nextstrain.org/ncov/europe",
        "North America": "nextstrain.org/ncov/north-america",
        "Oceania": "nextstrain.org/ncov/oceania",
        "South America": "nextstrain.org/ncov/south-america"
    }

    starters = [
        ("Check out new sequences from ", " on "),
        ("New sequences from ", " can be found on "),
        ("You can see new sequences from ", " on "),
        ("You can find new sequences from ", " on "),
        ("New sequences from ", " can be seen on ")
    ]

    starters_split = [
        ("Check out new sequences from ", " below"),
        ("New sequences from ", " can be found below"),
        ("You can see new sequences from ", " below"),
        ("You can find new sequences from ", " below"),
        ("New sequences from ", " can be seen below")
    ]

    the = ["USA", "United Kingdom", "Democratic Republic of the Congo"]

    counts_country = {region: {country: sum(counts[country].values()) for country in total_lab_collection[region]} for region in total_lab_collection}
    total = sum(sum(counts_country[region].values()) for region in counts_country)

    start_tweet = "Thanks to #opendata sharing via @GISAID, we've updated nextstrain.org/ncov with " + str(
        total) + " new #COVID19 #SARSCoV2 sequences!"
    char_total = 230
    char_available = char_total - len("Check out the new sequences from on ") - len("(Thanks to )") - len("1/1")

    tweet_collection_full = {}
    tweet_collection_split = {}
    lengths = {}
    for region in lab_collection:
        countries_list = list(lab_collection[region].keys())
        length_prediction = [len(country) + len(", ".join(lab_collection[region][country])) for country in lab_collection[region]]
        if sum(length_prediction) > char_available:
            countries_extra = [] #extra large countries
            while length_prediction and max(length_prediction) > char_available:
                i = np.argmax(length_prediction)
                countries_extra.append([countries_list[i]])
                countries_list.pop(i)
                length_prediction.pop(i)

            if countries_list:
                countries = []

                while(sum(length_prediction) > char_available):
                    length_prediction_sum = np.cumsum(length_prediction)
                    k = np.argmax(length_prediction_sum > char_available)
                    countries.append(countries_list[:k])
                    countries_list = countries_list[k:]
                    length_prediction = length_prediction[k:]

                countries.append(countries_list)
                countries += countries_extra
            else:
                countries = countries_extra

            for i, countries_list in enumerate(countries, start=1):

                h = []
                for country in countries_list:
                    for l in lab_collection[region][country]:
                        if l not in h:
                            h.append(l)
                c = ["the " + country if country in the else country for country in countries_list]
                r = region
                if i > 1:
                    r += str(i)
                tweet_collection_split[r] = (c, h)
        else:
            h = []
            for country in lab_collection[region]:
                for l in lab_collection[region][country]:
                    if l not in h:
                        h.append(l)

            c = ["the " + country if country in the else country for country in countries_list]
            tweet_collection_full[region] = (c, h)
            lengths[region] = len(", ".join(c)) + len(", ".join(h)) + len(links.get(region, ""))

    tweet = [(start_tweet + "\n\n", "\n\n[pic_Global]")]
    while lengths:
        current_region = min(lengths, key=lengths.get)
        best_partner = ""
        current_length = lengths[current_region]
        for region, length in sorted(lengths.items(), key=lambda x: x[1]):
            if region == current_region:
                continue
            if current_length + length > char_available:
                break
            best_partner = region

        lengths.pop(current_region)
        c = tweet_collection_full[current_region][0]
        h = tweet_collection_full[current_region][1]
        p = "[pic_" + current_region.replace(" ", "") + "]"
        l = links.get(current_region, "")
        if best_partner != "":
            current_length += lengths[best_partner]
            lengths.pop(best_partner)
            c += tweet_collection_full[best_partner][0]
            h += tweet_collection_full[best_partner][1]
            l += " and " + links[best_partner]
            p += " " + "[pic_" + best_partner.replace(" ", "") + "]"

        c = ", ".join(c[:-1]) + " and " + c[-1] if len(c) > 1 else c[0]
        h = " ".join(h) if current_length > char_available else ", ".join(h)
        starter = random.choice(starters)
        s = starter[0] + c + starter[1] + l + ".\n\n"
        s += "(Thanks to " + h + ")\n\n"
        tweet.append((s, "\n\n" + p))

    for region in tweet_collection_split:
        c = tweet_collection_split[region][0]
        h = tweet_collection_split[region][1]
        p = "[pic_" + region.replace(" ", "") + "]"
        if region in links:
            starter = random.choice(starters)
            l = links[region]
        else:
            starter = random.choice(starters_split)
            l = ""
        c = ", ".join(c[:-1]) + " and " + c[-1] if len(c) > 1 else c[0]
        if len(", ".join(c)) + len(", ".join(h)) + len(l) > char_available:
            h = " ".join(h)
        else:
            h = ", ".join(h)

        s = starter[0] + c + starter[1] + l + ".\n\n"
        s += "(Thanks to " + h + ")\n\n"
        tweet.append((s, "\n\n" + p))

    with open(path_to_outputs + "tweet_resources.txt", "a") as out:
        out.write("===============================\n\n")
        for i, t in enumerate(tweet):
            (s, p) = t
            out.write(s + str(i+1) + "/" + str(len(tweet)) + p + "\n\n\n")

def prepare_tweet_new_format(counts, rare_labs):
    links = {
        "Africa": "nextstrain.org/ncov/africa",
        "Asia": "nextstrain.org/ncov/asia",
        "Europe": "nextstrain.org/ncov/europe",
        "North America": "nextstrain.org/ncov/north-america",
        "Oceania": "nextstrain.org/ncov/oceania",
        "South America": "nextstrain.org/ncov/south-america"
    }

    counts_country = {region: {country: sum(counts[country].values()) for country in lab_collection[region]} for region
                      in lab_collection}
    total = sum(sum(counts_country[region].values()) for region in counts_country)

    start_tweet = "Thanks to #opendata sharing by @GISAID, we've updated nextstrain.org/ncov with " + str(
        total) + " new #COVID19 #SARSCoV2 sequences!"
    char_total = 260
    char_available = char_total - len("Check out the new sequences from on ") - len("(Thanks to )") - len("1/1")
    char_available_first = char_available - len(start_tweet)








path_to_input = "scripts/developer_scripts/inputs_new_sequences/"
path_to_outputs = "scripts/developer_scripts/outputs_new_sequences/"
table_file_name = path_to_input + "Who to Tag in Nextstrain Update Posts COVID-19.xlsx"
today = str(datetime.datetime.now())[:10]

if __name__ == '__main__':
    data, list_of_strains = read_data(path_to_input)
    data = check_dates(data, today)
    data = check_flagged_properties(data)
    #plot_dates(data, path_to_outputs + "plots/")
    counts = print_counts(data)
    lab_collection = collect_labs(data, table_file_name)

    rare_countries, rare_labs = check_for_recency(counts, list_of_strains, lab_collection, "data/", table_file_name)

    # Special checks for individual user requirements (e.g. produce concise overview over strain names, provide all new
    # sequences from certain countries, etc.)
    overview_with_dates(data, path_to_outputs + "strains_overview.txt")
    prepare_tweet(counts, lab_collection, rare_labs)
    prepare_tweet_new_format(counts, rare_labs)
