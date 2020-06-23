"""
Estimate TMRCA assuming a star topology and a poisson mutation process
"""
import argparse, math
import numpy as np
import datetime
import matplotlib.pyplot as plt
from Bio import Phylo
from treetime.utils import numeric_date
from augur.utils import read_node_data, read_metadata, get_numerical_dates

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate TMRCA assuming a star topology and a poisson mutation process",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--node-data", required=True, help="JSON with ancestral reconstruction")
    parser.add_argument("--metadata", required=True, help="JSON with ancestral reconstruction")
    parser.add_argument("--tree", required=True, help="newick tree")
    parser.add_argument("--output", required=True, help="figure file")
    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')

    metadata, columns = read_metadata(args.metadata)
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')
    node_data = read_node_data(args.node_data, args.tree)

    tips = {}
    for n in T.get_terminals():
        if type(dates[n.name])==list:
            continue
        tips[n.name] = {'numdate':dates[n.name],
                        'mutations': []}
        path = T.root.get_path(target = n)
        for c in path:
            tips[n.name]['mutations'].extend([x for x in node_data['nodes'][c.name]['muts'] if not (x[0] in ['N','-'] or x[-1] in ['N','-'])])

    tmrca = np.linspace(2019.7, np.min([x['numdate'] for x in tips.values()])-0.001,101)

    tsum = np.sum([np.mean(v['numdate']) for v in tips.values()])
    ntips = len(tips)
    L=29000
    for mu in [3e-4, 5e-4, 1e-3]:
        logp = -mu*(tsum - ntips*tmrca)*L
        for tip in tips.values():
            logp += len(tip['mutations'])*np.log(tip['numdate']-tmrca)
        p  = np.exp(logp)
        p /= p.sum()
        plt.plot(tmrca, p, label=f"rate={mu:1.1e} per site and year", lw=2)

    plt.title('TMRCA of 2019-nCov assuming a star tree\nand Poisson statistics of mutations', fontsize=16)
    plt.xlabel('TMRCA', fontsize=16)
    plt.ylabel('Probability density', fontsize=16)
    ticks = ['2019-10-01', '2019-11-01', '2019-12-01', '2020-01-01']
    plt.legend(loc=2, fontsize=12)
    plt.tick_params(labelsize=12)
    plt.xticks([numeric_date(datetime.datetime.strptime(x, '%Y-%m-%d')) for x in ticks], ticks)
    plt.savefig(args.output)
