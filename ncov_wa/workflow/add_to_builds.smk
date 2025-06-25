from dateutil import relativedelta

# Calculate dates
d = date.today()
two_m = d - relativedelta.relativedelta(months=2)
four_m = d - relativedelta.relativedelta(months=4)
six_m = d - relativedelta.relativedelta(months=6)
one_y = d - relativedelta.relativedelta(years=1)

# Set earliest_date & latest_date in builds
if "wa_two_mon" in config["builds"]:
    config["builds"]["wa_two_mon"]["earliest_date"]= two_m.strftime('%Y-%m-%d')
    config["builds"]["wa_two_mon"]["background_date"] = four_m.strftime('%Y-%m-%d')

if "wa_four_mon" in config["builds"]:
    config["builds"]["wa_four_mon"]["earliest_date"]= four_m.strftime('%Y-%m-%d')
    
if "wa_six_mon" in config["builds"]:
    config["builds"]["wa_six_mon"]["earliest_date"]= six_m.strftime('%Y-%m-%d')

if "wa_1y" in config["builds"]:
    config["builds"]["wa_1y"]["earliest_date"]= one_y.strftime('%Y-%m-%d')
