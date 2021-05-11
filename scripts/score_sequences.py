import argparse, json
from Bio import Phylo, SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict
from augur.distance import read_distance_map
# from binding_calculator import BindingCalculator
from scipy.optimize import minimize
from scipy.stats import linregress
from datetime import datetime

maps = {
"c1":"defaults/distance_maps/class1.json",
#"c1_dms":"defaults/distance_maps/class_1_dms.json",
"c1_dms_mean":"defaults/distance_maps/class_1_mean_dms.json",
"c1_dms_total":"defaults/distance_maps/class_1_total_dms.json",
"c2":"defaults/distance_maps/class2.json",
#"c2_dms":"defaults/distance_maps/class_2_dms.json",
"c2_dms_mean":"defaults/distance_maps/class_2_mean_dms.json",
"c2_dms_total":"defaults/distance_maps/class_2_total_dms.json",
"c3":"defaults/distance_maps/class3.json",
#"c3_dms":"defaults/distance_maps/class_3_dms.json",
"c3_dms_mean":"defaults/distance_maps/class_3_mean_dms.json",
"c3_dms_total":"defaults/distance_maps/class_3_total_dms.json",
#"c4_dms":"defaults/distance_maps/class_4_dms.json",
"c4_dms_mean":"defaults/distance_maps/class_4_mean_dms.json",
"c4_dms_total":"defaults/distance_maps/class_4_total_dms.json",
#"conv_serum_dms":"defaults/distance_maps/convalescent_serum_dms.json",
"conv_serum_dms_mean":"defaults/distance_maps/convalescent_serum_mean_dms.json",
"conv_serum_dms_total":"defaults/distance_maps/convalescent_serum_total_dms.json",
#"moderna_serum_dms":"defaults/distance_maps/Moderna_vaccine_serum_dms.json",
"moderna_serum_dms_mean":"defaults/distance_maps/Moderna_vaccine_serum_mean_dms.json",
"moderna_serum_dms_total":"defaults/distance_maps/Moderna_vaccine_serum_total_dms.json",
"S1":"defaults/distance_maps/S1.json",
"ace2":"defaults/distance_maps/ace2.json",
"NTD":"defaults/distance_maps/NTD.json",
"cleavage":"defaults/distance_maps/cleavage.json",
}

def ordinal_to_CW(x):
    return (x - datetime(2021,1,4).toordinal())//7

def CW_to_ordinal(x):
    return x*7 + datetime(2021,1,4).toordinal()

def logistic(t, s, tau):
    est = np.exp(s*(t-tau))
    return est/(1+est)

def logit(x):
    return np.log(x/(1-x))

def logit_inv(x):
    return np.exp(x)/(1+np.exp(x))

def fit_logistic(observations, all_time_points, method="Powell"):
    def cost(x, obs, t):
        st_obs = np.clip(x[0]*(obs-x[1]), -100,100)
        st_all = np.clip(x[0]*(t-x[1]), -100,100)
        return -np.sum(st_obs) + np.sum(np.log(1+np.exp(st_all)))

    def jac(x, obs, t):
        st_obs = np.clip(x[0]*(obs-x[1]), -100,100)
        st_all = np.clip(x[0]*(t-x[1]), -100,100)
        est = np.exp(st_all)
        denom = est/(1+est)
        return np.array([ -np.sum(obs-x[1]) + np.sum((t-x[1])*denom),
                          x[0]*(len(obs) - np.sum(denom))])

    y_obs, bins = np.histogram(observations, bins=4)
    y_all, b = np.histogram(all_time_points, bins=bins)

    logit = np.log((y_obs+1)/(y_all-y_obs+1))
    bc = 0.5*(bins[1:] + bins[:-1])
    res = linregress(bc, logit)
    sol = minimize(cost, [res.slope, np.clip(-res.intercept/res.slope, 737000,738500)],
                   args=(observations, all_time_points), jac=jac, method=method, tol=1e-5)
    return sol

def fit_group(data, G, verbose=False, method="BFGS"):
    obs = data.loc[G, "ordinal"]
    sol = fit_logistic(obs, data.ordinal, method=method)
    today = datetime.today().toordinal()
    if np.abs(sol['x'][0])>1:
        print(" ...failed")
        return None

    sig_s = np.sqrt(sol['hess_inv'][0,0])
    sig_tau = np.sqrt(sol['hess_inv'][1,1])

    if verbose:
        print(f": s={sol['x'][0]:1.2f}±{sig_s:1.3f}/day, freq={logistic(today, sol['x'][0], sol['x'][1]):1.2e}")
    return  {"s":sol['x'][0],
            "sigma_s":sig_s,
            "t50_float":sol['x'][1],
            "covariance":sol['hess_inv'],
            "sigma_t50": sig_tau,
            "t50":datetime.fromordinal(int(sol['x'][1])).strftime('%Y-%m-%d'),
            "current_freq": logistic(today, sol['x'][0], sol['x'][1]),
            "total_count": len(G)}

def fit_groups(data, region, min_date, method="BFGS"):
    data_subset = data.loc[(data.region==region)&(data.date>min_date)]
    by_spike_groups = data_subset.groupby(by='S1_mut_str')
    growth = {}
    for mut, G in by_spike_groups.groups.items():
        if len(G)<50:
            continue
        print("Fitting", mut, len(G), end='')
        res = fit_group(data_subset, G, method=method)
        if res:
            growth[mut] = res

    return growth

def plot_group(data, G, n_sigma=1.96):
    plt.figure()
    ax=plt.subplot(111)

    res = fit_group(data, G)
    total_count = data.groupby('week').count()['S1']
    sub_count = data.loc[G].groupby('week').count()['S1']
    freq = (sub_count/total_count).fillna(0)
    ordinals = np.array([CW_to_ordinal(x) for x in freq.index])
    ax.plot(freq.index, freq, 'o')

    if n_sigma:
        slope = res['s']
        t50 = res['t50_float']
        CW_points = np.linspace(freq.index[0], freq.index[-1], 100)
        x_vals = CW_to_ordinal(CW_points)
        y_vals = slope*(x_vals - t50)

        cov_matrix = np.copy(res['covariance'])
        dev = n_sigma*np.array([np.sqrt(cov_matrix.dot(np.array([x-t50, -slope])).dot(np.array([x-t50,-slope])))
                                for x in x_vals])
        ax.fill_between(CW_points, logit_inv(y_vals-dev), logit_inv(y_vals+dev), alpha=0.2)
        CI = f'±{n_sigma*res["sigma_s"]:1.3f}'
    else: CI=''

    ax.plot(freq.index, logistic(ordinals, res["s"], res['t50_float']), '-', label=f's={res["s"]:1.3f}{CI}/day')
    plt.legend()
    plt.ylabel('frequency')
    plt.xlabel('week')


def analyze_growth(data, min_date):
    growth_by_region = defaultdict(dict)
    for region in ['Europe', 'North America', 'South America', 'Africa', 'Asia']:
        growth = fit_groups(data, region=region, min_date=min_date)
        for muts in growth:
            growth_by_region[muts][region] = growth[muts]
    return growth_by_region

def nsp6_score(orf1a_muts):
    return int(3675 in [x[1] for x in orf1a_muts])

def RBD_score(data):
    return np.sum([data[f'c{i}_dms_mean'] for i in range(1,5)], axis=0)

def RBD_categorical_score(data, threshold=0.1):
    return np.sum([data[f'c{i}']>threshold for i in range(1,4)], axis=0)


def apply_map(spike_mutations, s_map):
    d_val = s_map['default']
    if type(list(s_map['map']['S'].values())[0])==dict:
        # amino acid specific map
        return spike_mutations.apply(lambda x:np.sum([s_map['map']['S'].get(m[1]-1, {}).get((m[0], m[2]), d_val)
                                                      for m in x]))
    else:
        return spike_mutations.apply(lambda x:np.sum([s_map['map']['S'].get(m[1]-1, d_val) for m in x if m[2] not in ['X']]))


def report_evoscore_variants(data, nMax = 100, minCount=5,
                            fields= ['RBD', 'RBD_cat', 'total', "conv_serum_dms_mean"]):
    by_spike_groups = data.groupby(by='S1_mut_str')
    grouped_scores = by_spike_groups.mean()

    # print the nMax sequences with the highest score for each map and the sum
    for a in fields:
        print("\n\n### ATTRIBUTE:",a)
        print(f"#{'mutations':<80}\tmin date\tmax date\tnumber\t{'PANGO'}\t{'score'}")
        for muts, row in grouped_scores.sort_values(by=a)[-nMax:].iterrows():
            G = data.loc[by_spike_groups.groups[muts]]
            if len(G)<minCount:
                continue
            print(f"{muts:<80}\t{G.date.min()}\t{G.date.max()}\t{len(G)}\t{G.pango_lineage.mode()[0]}\t{row[a]:1.2f}")


def report_growing_variants(data, growth_by_region, focal_region='Europe'):
    by_spike_groups = data.groupby(by='S1_mut_str')
    constellations_of_interest = []
    for m in growth_by_region:
        if (focal_region in growth_by_region[m]
            and growth_by_region[m][focal_region]['slope']>0.01
            and growth_by_region[m][focal_region]['total_count']>50):
            lineage = data.loc[by_spike_groups.groups[m]].pango_lineage.mode()
            total_score = data.loc[by_spike_groups.groups[m]].total.mean()
            report_str = ", ".join([f"{r}, s={growth_by_region[m][r]['slope']:1.2f} n={growth_by_region[m][r]['total_count']}" for r in growth_by_region[m]])
            constellations_of_interest.append((lineage[0], m,  report_str, total_score,
                                               growth_by_region[m][focal_region]['slope']))

    for c in sorted(constellations_of_interest, key=lambda x:x[-1]):
        print(f"{c[0]:<15}\t{c[-1]:1.3f}\t{c[-2]:1.2f}\t{c[1]}\t{c[2]}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="count distances based on distance maps with counting gaps as one event",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-summary', type=str)
    parser.add_argument('--metadata', type=str)
    parser.add_argument('--min-date', type=str)
    parser.add_argument('--min-sub-date', type=str)
    parser.add_argument('--country')
    parser.add_argument('--region')
    args = parser.parse_args()

    print("loading mutations")
    raw_mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0)[['ORF1a', 'S']].fillna('')
    mutations = pd.concat([raw_mutations[col].apply(lambda x: tuple([(y[0], int(y[1:-1]), y[-1]) for y in x.split(',')]) if x else tuple())
                           for col in raw_mutations.columns], axis=1)

    print("loading metadata")
    meta = pd.read_csv(args.metadata, sep='\t', index_col=0).fillna('')
    data = pd.concat([mutations, meta[["date", "date_submitted", "country", "region",
                                       "Nextstrain_clade", "pango_lineage"]]], axis=1).loc[mutations.index]

    good_sequences = data['S'].apply(lambda x:len(x)<30)
    data = data.loc[good_sequences]
    good_sequences = data['date'].astype(str).apply(lambda x:len(x)==10 and 'X' not in x)
    data = data.loc[good_sequences]

    if args.min_date:
        data = data.loc[data.date>args.min_date]
    if args.min_sub_date:
        data = data.loc[data.date_submitted>args.min_sub_date]
    if args.country:
        data = data.loc[data.country==args.country]
    if args.region:
        data = data.loc[data.region==args.region]

    data["ordinal"] = data['date'].apply(lambda x:datetime.strptime(x, '%Y-%m-%d').toordinal())
    data["week"] = data['ordinal'].apply(ordinal_to_CW)

    data["S1_mut_str"] = data.S.apply(lambda muts: ','.join([f'{x[0]}{x[1]}{x[2]}' for x in muts if x[1]<690]))

    data['nsp6'] = data['ORF1a'].apply(nsp6_score)
    for n, map_file in maps.items():
        print("processing",n)
        s_map = read_distance_map(map_file)
        data[n] = apply_map(data['S'], s_map)
        data[n] /= np.percentile(data[n], 99.99)

    data['RBD'] = RBD_score(data)
    data['RBD_cat'] = RBD_categorical_score(data, threshold=0.25)
    data['total'] = 2*data['RBD'] + 2*data['NTD'] + data['ace2'] + data['nsp6'] + data['cleavage']

    report_evoscore_variants(data)

    # growth_by_region = analyze_growth(data, min_date="2021-02-01")

    # report_growing_variants(data, growth_by_region, focal_region='Asia')
