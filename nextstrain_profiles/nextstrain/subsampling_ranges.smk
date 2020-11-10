import datetime

# Set subsampling max date to today.
today = datetime.date.today()

# Set the earliest date to roughly 4 months ago (18 weeks).
early_late_cuotff = today - datetime.timedelta(weeks=18)

for build in config["subsampling"]:
    for scheme in config["subsampling"][build]:
        if "_early" in scheme:
            config["subsampling"][build][scheme]["max_date"] = f"--max-date {early_late_cuotff.strftime('%Y-%m-%d')}"
        if "_late" in scheme:
            config["subsampling"][build][scheme]["min_date"] = f"--min-date {early_late_cuotff.strftime('%Y-%m-%d')}"
