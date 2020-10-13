import os
import datetime
import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()



def bold(s):
    return('\033[1m' + s + '\033[0m')

# Cut a string of the format "key: content" into a tuple (key, content)
def cut(s):
    key = s.split(":")[0]
    content = ":".join(s.split(":")[1:])[1:]
    return (key, content)

# Read all files within given directory that start with "metadata-changes" in a sorted manner. Contents are stored in a
# dictionary first by GISAID epi isl, then by info-type (e.g. date, division, host)
def read_data(path):

    data = {} # added and changed sequences

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
                        if id in data:
                            print("Attention, same sequence added two times! (" + id + ")")
                        data[id] = {}
                        for j in range(-2,24):
                            k = metadata_changes[i+j].strip()
                            (key, content) = cut(k)
                            if key != "gisaid_epi_isl":
                                data[id][key] = content

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

    return (data)

# Double check sample dates for invalid format or unrealistic / impossible date
def check_dates(data, today):

    invalid_sample_date = {}
    suspicious_sample_date = {}

    for id in list(data.keys()):
        date = data[id]["date"]
        strain = data[id]["strain"]

        if len(date) != len(today):
            invalid_sample_date[strain] = date
            data.pop(id)
            continue

        day = int(date[8:])
        month = int(date[5:7])
        year = int(date[:4])

        day_today = int(date[8:])
        month_today = int(date[5:7])
        year_today = int(date[:4])


        #check for past dates
        if year < 2019 or ((year) == 2019 and month < 12):
            invalid_sample_date[strain] = date
            data.pop(id)
            continue

        # Check for future dates
        if (year > year_today) or (year == year_today and month > month_today) or (year == year_today and month == month_today and day > day_today):
            invalid_sample_date[strain] = date
            data.pop(id)
            continue

        #Check for early dates
        if (year == 2020 and (month == 2 or month == 1)) or year == 2019:
            suspicious_sample_date[strain] = date

    print("\n----------------------------------------------\n")
    print("Invalid sample dates (please check whether all are automatically excluded):")
    for strain in invalid_sample_date:
        print(strain + ": " + invalid_sample_date[strain])

    print("\nEarly sample dates (might require double checking depending on country):")
    for strain in suspicious_sample_date:
        print(strain + ": " + suspicious_sample_date[strain])

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
        os.remove(path + f)

    for country in dates_by_country:
        dates = list(dates_by_country[country].keys())
        values = list(dates_by_country[country].values())
        plt.figure()
        plt.bar(dates, values)
        plt.title(country)
        plt.xticks(rotation=45, ha="right", size = 7)
        plt.savefig(path + country)


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
    for country in counts:
        sum_country = 0
        for division in counts[country]:
            sum_country += counts[country][division]
        sum_total += sum_country

    print("\n----------------------------------------------\n")
    print("Total counts: " + str(sum_total))

    for country in counts:
        s = country + ": "
        sum_country = 0
        for division in counts[country]:
            sum_country += counts[country][division]
        s = s + str(sum_country)
        if len(counts[country]) == 1:
            s = s + " (" + division + ")"
        else:
            s = s + " (" + ", ".join([str(counts[country][division]) + " " + division for division in counts[country]]) + ")"
        print(s)


# Collect all submitting and originating labs as well as authors and try to infer as many twitter handles as possible
# from a given excel file
def collect_labs(data, table_file_name):
    print("\n----------------------------------------------")
    submitting_labs = {}
    originating_labs = {}
    authors = {}
    for id in data:
        country = data[id]["country"]
        submitting_lab = data[id]["submitting_lab"]
        originating_lab = data[id]["originating_lab"]
        author = data[id]["authors"]

        if country not in submitting_labs:
            submitting_labs[country] = []
        if submitting_lab not in submitting_labs[country]:
            submitting_labs[country].append(submitting_lab)

        if country not in originating_labs:
            originating_labs[country] = []
        if originating_lab not in originating_labs[country] and originating_lab != submitting_lab:
            originating_labs[country].append(originating_lab)

        if country not in authors:
            authors[country] = []
        if author not in authors[country]:
            authors[country].append(author)

    excel_table = pd.read_excel(table_file_name, index_col=0, skiprows=1)
    lab_dictionary = {}
    for country, row in excel_table.iterrows():
        description = row["Who"]
        handle = row["Who to tag"]
        if country not in lab_dictionary:
            lab_dictionary[country] = {}
        if description in lab_dictionary[country]:
            print("Warning: lab description is found two times in excel table in same country (" + country + ", " + description + ")" )
        lab_dictionary[country][description] = handle


    print("\nSubmitting labs:\n(Note: small differences in spelling might cause lab to not be identified. Consider adjusting the spelling in the spreadsheet!)\n")
    for country in submitting_labs:
        s = country + ":\n"
        for lab in submitting_labs[country]:
            s += lab + ": "
            if country in lab_dictionary and lab in lab_dictionary[country]:
                s += bold(lab_dictionary[country][lab])
            else:
                s += bold("?")
            s += "\n"
        print(s)

    print("----------------------------------------------\n")
    print("Originating labs (only printed if different from submitting lab):\n")
    for country in originating_labs:
        s = country + ":\n"
        for lab in originating_labs[country]:
            s += lab
            if country in lab_dictionary and lab in lab_dictionary[country]:
                s += ": " + bold(lab_dictionary[country][lab])
            s += "\n"
        print(s)

    print("----------------------------------------------\n")
    print("Authors:\n")
    for country in authors:
        s = country + ":\n"
        for author in authors[country]:
            s += author + "\n"
        print(s)


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
        for country in data_sorted:
            myfile.write(country + "\n")
            for s in data_sorted[country]:
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
        for country in special_strains:
            myfile.write(country + "\n")
            for date in sorted(special_strains[country]):
                myfile.write(date + ": " + str(special_strains[country][date]) + "\n")
            myfile.write("\n")




path_to_input = "scripts/developer_scripts/inputs_new_sequences/"
path_to_outputs = "scripts/developer_scripts/outputs_new_sequences/"
table_file_name = path_to_input + "Who to Tag in Nextstrain Update Posts COVID-19.xlsx"
today = str(datetime.datetime.now())[:10]

data = read_data(path_to_input)
data = check_dates(data, today)
plot_dates(data, path_to_outputs + "plots/")
print_counts(data)
collect_labs(data, table_file_name)

# Special checks for individual user requirements (e.g. produce concise overview over strain names, provide all new
# sequences from certain countries, etc.)
overview_with_dates(data, path_to_outputs + "strains_overview.txt")
filter_for_date_region(data, path_to_outputs, ("Europe", 6))
