import argparse, json
from matplotlib import pyplot as plt
from datetime import datetime
from Bio import Phylo, SeqIO
import pandas as pd
import json
import numpy as np
from collections import defaultdict

def smooth(x, smoothing=None):
    # assumes no empty weeks, need to fix
    if smoothing is None:
        smoothing = np.exp(-np.arange(-5,5)**2/2)

    return np.convolve(x, smoothing, mode='same')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="calculate frequencies of specific mutations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-summary', type=str)
    parser.add_argument('--metadata', type=str)
    args = parser.parse_args()


    meta = pd.read_csv(args.metadata, sep='\t', index_col=0)
    mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0).fillna('')
    meta = pd.concat([meta, mutations['S']], axis=1)

    ref_date = datetime.strptime('2020-W1-1', "%Y-W%U-%d").toordinal()
    meta = meta.loc[meta.date.apply(lambda x:len(x)==10 and 'X' not in x),:]
    meta["CW"] = meta.date.apply(lambda x:int((datetime.strptime(x, '%Y-%m-%d').toordinal()-ref_date)//7 + 1))
    seq_counts= {x:len(y) for x,y in meta.groupby(['CW', 'country'])}

    # from https://raw.githubusercontent.com/samayo/country-json/master/src/country-by-population.json
    with open('data/country-by-population.json') as fh:
        pops = {x['country']: x['population'] for x in json.load(fh)}
    pops['Republic of the Congo'] = pops['Congo']
    pops["Côte d'Ivoire"] = pops['Ivory Coast']
    weights = {(cw, country): pops.get(country, 1e7)/seq_counts[(cw, country)]/1e6 for cw, country in seq_counts}
    meta['weight'] = [weights[(row['CW'], row['country'])] for s, row in meta.iterrows()]

    query =['E484K']
    query =['A222V']
    query =['N501Y']
    query =['E484K', 'N501Y']
    fig = plt.figure()
    plt.title(','.join(query))
    for region in ['Europe', 'Africa', 'North America', 'South America']:
        subset = meta.loc[meta.region==region]
        ind = np.array([all([x in muts for x in query]) for muts in subset.S.fillna('')])
        m = subset.loc[ind].groupby('CW').sum().weight
        denom = subset.groupby('CW').sum().weight
        max_week = int(max(denom.index))+1
        dates = pd.DataFrame([datetime.fromordinal(int(7*x + ref_date)) for x in range(max_week)], index=range(max_week))
        d = pd.concat([dates, m,denom], axis=1).fillna(0)
        d = d.loc[d.index>0]
        plt.plot(d.iloc[:,0], smooth(d.iloc[:,1])/smooth(d.iloc[:,2]), label=region)

    plt.legend()
    fig.autofmt_xdate()
