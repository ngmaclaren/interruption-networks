import math
import itertools

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import InterruptionAnalysis as ia

import random
random.seed(12345)
# For analysis, I want to compare the count of each kind of triad in the empirical networks (summed across all empirical networks) against the count of the same kind of triad in n_sims number of samples of randomized copies of each graph. 
print('Loading graph data...')
data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
a = 0.99

gIDs = pd.unique(data['gID'])
igraphs = {}
vgraphs = {}
bgraphs = {}
for gID in gIDs:
    igraphs[gID] = nx.read_gml(f'./data/networks-iss/{gID}.gml')
    ngraphs[gID] = nx.read_gml(f'./data/networks-nss/{gID}.gml')
    bgraphs[gID] = nx.read_gml(f'./data/networks-both/{gID}.gml')
    vgraphs[gID] = nx.read_gml(f'./data/votenets/{gID}.gml')
    
n_sims = 5000

#graphs = igraphs
#imgpath = './img/network-motifs-5000runs-iss.svg'
#graphs = ngraphs
#imgpath = './img/network-motifs-5000runs-nss.svg'
graphs = bgraphs
imgpath = './img/network-motifs-5000runs-both.svg'

print('Generating RE graphs...')
sims = {}
for sim in range(n_sims):
    if sim % 10 == 0:
        print(sim)
    batch = [ia.randomize_edges(g) for g in graphs.values()]
    sims[sim] = sum([pd.Series(nx.triadic_census(g)) for g in batch])

empirical_triads = sum([pd.Series(nx.triadic_census(g)) for g in graphs.values()])
simulated_triads = pd.DataFrame(sims).T

# now, count how many times the value in col(x) in simulated_triads is greater than the value in row(x) in empirical_triads
ps = pd.Series(np.nan, index = empirical_triads.index)
for triad in list(empirical_triads.index):
    dist = simulated_triads[triad]
    emp = empirical_triads.loc[triad]
    p = len([d for d in dist if d > emp])/len(dist)
    ps.loc[triad] = p

# significant triads have p < 0.01
# Bonferroni correction is α' = α/m, where α is the desired α value and m is the number of hypotheses
# but it doesn't necessarily matter because I just want to do the plot below, pulling the right count and p-value given the new data structure
pvals = ps.sort_values()
m = len(pvals)
α = 0.05
alphas = [α/(m - i) for i in range(m)]
compdata = pd.DataFrame({'p': pvals, 'comp': alphas})
compdata['result'] = compdata['p'] < compdata['comp']
compdata['count'] = empirical_triads
print(compdata)

plt.rcParams['font.size'] = 12
fig, axs = plt.subplots(4, 4, figsize = (10, 10))

pos = {'a': (1, 1), 'b': (6, 1),'c': (3.5, 5.333)}
# scaling...
newpos = {a: (b[0]/10, b[1]/10) for a, b in pos.items()}

# empirical counts are in empirical_triads, indexed by triad name
# p-values are in ps, indexed by triad name
triads = list(empirical_triads.index)
for ax, triad in zip(axs.flatten(), triads):
    g = nx.triad_graph(triad)
    nx.draw(g, pos = newpos, with_labels = False, ax = ax,
            node_size = 50, width = .5, node_color = 'k')
    ax.set_title(triad)

    ax.text(newpos['a'][0], newpos['a'][1] - 0.05, r'$n$ = ' + f'{empirical_triads.loc[triad]}',
            horizontalalignment = 'left', verticalalignment = 'top')
    ax.text(newpos['b'][0], newpos['b'][1] - 0.05, r'$p$ = ' + f'{ps.loc[triad]}',
            horizontalalignment = 'right', verticalalignment = 'top')
    ax.set_axis_off()
    ax.set_xlim((newpos['a'][0] - 0.05, newpos['b'][0] + 0.05))
    ax.set_ylim((newpos['a'][1] - 0.1, newpos['c'][1] + 0.05))
fig.tight_layout()
fig.savefig(imgpath)
fig.show()
plt.rcParams['font.size'] = 10
print('done')
