"""
Estimate branching parameters
"""
import argparse
import numpy as np
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from treetime.utils import numeric_date, datestring_from_numeric

def mean_non_extinct(T,b,d):
    '''
    average number of individuals in a supercritical branching
    process with birth rate b and death rate d after time T.
    The birth rate b can be a numpy array, others need to be scalars
    '''
    s = b-d
    res = np.ones_like(b)/(1+d*T)
    ind = np.abs(s/d)>1e-6
    res[ind] = ((b*(np.exp(s*T)-1))/s)[ind]
    return res


def LH(n,T,b,d):
    '''
    return the probability of observing n cases in a branching process
    after time T with birth and death rate b,d, conditional
    on non-extinction
    '''
    s = b-d
    return np.exp(-n/mean_non_extinct(T,b,d))/mean_non_extinct(T,b,d)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate Tmrca assuming a star topology and a poisson mutation process",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--infectious-period", type=float, default=10, help="infections period in days")
    parser.add_argument("--start", type=str, default='2019-11-27', help="start of the outbreak")
    parser.add_argument("--population", nargs='+', type=int, default = [300,1000,3000,10000],
                         help="number of individuals now")
    parser.add_argument("--output", required=True, help="figure file for line graph")
    args = parser.parse_args()

    d = 365.0/args.infectious_period
    n_vals = np.array(args.population)

    T = numeric_date() - numeric_date(datetime.datetime.strptime(args.start, '%Y-%m-%d'))
    b_vals = d*np.linspace(0.5, 6, 111)

    weeks = int(np.round(T*365/7))
    inf_period = int(365/d)

    fs=16
    plt.figure()
    plt.title(f"Start: {weeks} weeks ago = {args.start}. Infectious period {inf_period} days")
    for n in n_vals:
        lh =  LH(n,T,b_vals,d)
        lh /= lh.sum()
        lh /= (b_vals[1]-b_vals[0])/d
        plt.plot(b_vals/d, lh, lw=3, label=f'n={n}')
    plt.ylabel('Probability density', fontsize=fs)
    plt.xlim([0.9,4])
    plt.xlabel('R0', fontsize=fs)
    plt.legend(fontsize=0.8*fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.tight_layout()
    plt.savefig(args.output)

    # plt.savefig(f"R0_{weeks}_{inf_period}.png")
