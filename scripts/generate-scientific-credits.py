import pandas as pd

data = pd.read_csv('../data/metadata.tsv', sep='\t')
credits = data.groupby('originating_lab')['strain'].apply(list).to_dict()

detailed_ofile = open('../data/detailed_credits.md', 'a')
ofile = open('../data/credits.md', 'a')

for institution in sorted(list(credits.keys())):
    if institution == 'unknown':
        continue

    ofile.write('* '+institution+'\n')

    strains = sorted(credits[institution])
    detailed_ofile.write('* '+institution+'\n')
    for s in strains:
        detailed_ofile.write('\t* '+s+'\n')
    detailed_ofile.write('\n')

detailed_ofile.close()
ofile.close()
