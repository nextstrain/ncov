import argparse
import pandas as pd
import requests
import xlrd


def export_imf(output):
    # Read data from the web
    # Data source web page: <https://www.imf.org/external/datamapper/LP@WEO/OEMDC/ADVEC/WEOWORLD>
    excel_url = "https://www.imf.org/external/datamapper//export/excel.php?indicator=LP"
    response = requests.get(excel_url)
    excel_contents = response.content
    workbook = xlrd.open_workbook(file_contents=excel_contents, ignore_workbook_corruption=True)
    df = pd.read_excel(workbook,
        sheet_name='LP',
        comment='©',  # Ignore last line "© IMF …"
        na_values='no data',
    )

    # Rename columns to match names in metadata
    column_name_map = {
        'Population (Millions of people)': 'country',
        2022: 'weight',
    }
    df = df.rename(columns=column_name_map)

    # Keep only the columns used above
    df = df[column_name_map.values()]

    # Set country as index
    df = df.set_index('country')

    # Rename countries to match values in metadata
    country_name_map = {
        "Brunei Darussalam": "Brunei",
        "China, People's Republic of": "China",
        "Hong Kong SAR": "Hong Kong",
        "Korea, Republic of": "South Korea",
        "Kyrgyz Republic": "Kyrgyzstan",
        "Lao P.D.R.": "Laos",
        "Macao SAR": "Macao",
        "Taiwan Province of China": "Taiwan",
        "West Bank and Gaza": "Palestine",
    }
    df = df.rename(index=country_name_map)

    # Remove rows without a weight (or country)
    df = df.dropna(how='any')

    # Add missing countries
    # Data sourced from <https://www.cia.gov/the-world-factbook/about/archives/2023/field/population/country-comparison/>
    # TODO: consider using this data source for everything - it doesn't seem to have any missing countries
    df.loc['Syria', 'weight'] = 22.933531

    # Export
    df.to_csv(output, index=True, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create population sizes file",
    )

    parser.add_argument('--output', type=str, metavar="FILE", required=True, help="Path to output population sizes file")
    args = parser.parse_args()

    export_imf(args.output)
