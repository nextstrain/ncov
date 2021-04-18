import pandas as pd
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def nmf(V, k, max_iter=100, verbose=False):
    W = np.random.random(size=(V.shape[0],k))
    H = np.random.random(size=(k, V.shape[1]))

    for i in range(max_iter):
        H *= W.T.dot(V)/(W.T.dot(W).dot(H)+1e-10)
        W *= V.dot(H.T)/(W.dot(H).dot(H.T)+1e-10)
        if verbose:
            print(f"{i}: {np.sum((V - W.dot(H))**2)}")

    return W,H



data = pd.read_csv(r"https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_data.csv")

maps = {x:defaultdict(list) for x in data.condition_subtype.unique()}
sum_maps = {x:defaultdict(lambda: defaultdict(float)) for x in data.condition_subtype.unique()}
mean_maps = {x:defaultdict(lambda: defaultdict(float)) for x in data.condition_subtype.unique()}

for c, row in data.iterrows():
    maps[row.condition_subtype][(row.site, row.wildtype, row.mutation)].append(row.mut_escape)
    sum_maps[row.condition_subtype][row.site][row.condition] = row.normalized_site_total_escape
    mean_maps[row.condition_subtype][row.site][row.condition] = row.normalized_site_mean_escape

average_maps = {x:{} for x in maps}
for x in maps:
    for pos, wt, mut in maps[x]:
        val = np.mean(maps[x][(pos, wt, mut)])
        if val>0.05:
            if pos not in average_maps[x]:
                average_maps[x][pos]={}
            average_maps[x][pos][(wt, mut)]=val

for x in average_maps:
    out_json = {"default":0,
                "map":{"S":{ pos:[{'from':k[0], 'to':k[1], 'weight':y} for k,y in average_maps[x][pos].items()]
                            for pos in average_maps[x]}},
                "name":x}

    with open(f"defaults/distance_maps/{x.replace(' ', '_')}_dms.json", 'w') as fh:
        json.dump(out_json, fh, indent=2)

for x in sum_maps:
    out_json = {"default":0,
                "map":{"S":{}},
                "name":x}
    for pos in sum_maps[x]:
        out_json["map"]["S"][pos] = np.mean(list(sum_maps[x][pos].values()))

    with open(f"defaults/distance_maps/{x.replace(' ', '_')}_total_dms.json", 'w') as fh:
        json.dump(out_json, fh, indent=2)


for x in mean_maps: 
    out_json = {"default":0,
                "map":{"S":{}},
                "name":x}
    for pos in mean_maps[x]:
        out_json["map"]["S"][pos] = np.mean(list(mean_maps[x][pos].values()))

    with open(f"defaults/distance_maps/{x.replace(' ', '_')}_mean_dms.json", 'w') as fh:
        json.dump(out_json, fh, indent=2)

spectra = []
for i in [1,2,3,4]:
    spike = np.zeros(shape=(600))
    for p, v in mean_maps[f"class {i}"].items():
        spike[p] = np.mean(list(v.values()))
    spectra.append(spike)

class_spectra = np.array(spectra)

training = {}
types = average_maps.keys()
for name in types:
    tmp = []
    for c in data.loc[data.condition_subtype==name,'condition'].unique():
        spike = np.zeros(shape=(600))
        for p, v in mean_maps[name].items():
            spike[p] = v[c]
        tmp.append(spike)
    training[name] = np.array(tmp)

####
monoclonals = [k for k in training if 'class' in k]
rank = range(1,30)
for tsets in [monoclonals, ['convalescent serum'], ['Moderna vaccine serum'], 'all']:
    var_remaining = []
    train_set = np.concatenate(list(training.values()) if tsets=='all' else [training[k] for k in tsets])
    for k in rank:
        contributions, nmf_spectra =  nmf(train_set, k)
        pred = contributions.dot(nmf_spectra)
        var_remaining.append(np.sum((train_set - pred)**2)/len(train_set))
    plt.plot([0] + list(rank), [np.sum((train_set)**2)/len(train_set)] + var_remaining, 
            label = tsets if type(tsets)==str else ', '.join(tsets))

plt.legend()
plt.xlabel('rank of decomposition')
plt.ylabel('residual variance')
plt.yscale('log')
plt.savefig('NMF_residuals.png')

####
train_set = np.concatenate(list(training.values()))
fig, axs = plt.subplots(2,2, sharex=True, sharey=True, figsize=(12,6))
for ri, rank in enumerate([2,3,4,5]):
    ax = axs[ri//2, ri%2]
    contributions, nmf_spectra =  nmf(train_set, rank)
    relevant_sites = list(range(310,550))
    ax.plot(relevant_sites, nmf_spectra[:,relevant_sites].T)
    ax.set_title(f"rank {rank}")

    if ri%2==0:
        ax.set_ylabel('escape')
    if ri//2==1:
        ax.set_xlabel('position')

plt.savefig('NMF_spectra_by_rank.png')

rank = 5
contributions, nmf_spectra =  nmf(train_set, rank)
spectra = nmf_spectra
M = {}
for name in training:
    M[name] = {'experiments': data.loc[data.condition_subtype==name,'condition'].unique()}
    tmp = []
    for c in M[name]['experiments']:
        spike = np.zeros(shape=(600))
        for p, v in mean_maps[name].items():
            spike[p] = v[c]
        tmp.append(spike)

    M[name]['data'] = np.array(tmp).T
    partitions, res = np.linalg.lstsq(spectra.T, M[name]['data'])[:2]
    M[name]['residuals'] = res 
    M[name]['pred'] = spectra.T.dot(partitions) 
    M[name]['rel_residuals'] = res/np.sum(M[name]['data']**2, axis=0)
    M[name]['contributions'] = partitions 
    
    # plt.figure(name)
    # for i, ex in enumerate(M[name]['experiments'][:4]):
    #     plt.plot(range(300,550), M[name]['pred'][300:550, i], c=f'C{i%10}', ls='-')
    #     plt.plot(range(300,550), M[name]['data'][300:550, i], c=f'C{i%10}', ls='--')

for name in M:
    a = sns.clustermap(M[name]['contributions'].T, col_cluster=False, vmin=-0.1, vmax=1.5)
    ind = a.dendrogram_row.reordered_ind
    a.ax_heatmap.set_yticklabels([M[name]['experiments'][i] for i in ind], rotation=0)
    a.ax_heatmap.set_xticklabels([f'c {i+1}' for i in range(len(spectra))])
    plt.tight_layout()
    plt.title(name)
    plt.savefig(name.replace(' ', '_') + '_map.png')

    plt.figure()
    plt.plot(1-M[name]['rel_residuals'][ind], 'o-')
    plt.xticks(np.arange(len(ind)), [M[name]['experiments'][i] for i in ind], rotation=30, ha='right')
    plt.tick_params(labelsize=8)
    plt.title(name)
    plt.ylim([0,1])
    plt.ylabel("fraction of signal explained")
    plt.tight_layout()
    plt.savefig(name.replace(' ', '_') + '_explained.png')





