import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import InterruptionAnalysis as ia

readpath = './data/edgedir-sim'

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
votedata.set_index('pID', inplace = True)
surveydata = pd.read_csv('./data/speakingTime-data.csv', index_col = 0)
surveydata.set_index('pID', inplace = True)
a = 0.99

gIDs = pd.unique(data['gID'])
igraphs = {}
vgraphs = {}
bgraphs = {}
ngraphs = {}
for gID in gIDs:
    igraphs[gID] = nx.read_gml(f'./data/networks-iss/{gID}.gml')
    vgraphs[gID] = nx.read_gml(f'./data/votenets/{gID}.gml')
    bgraphs[gID] = nx.read_gml(f'./data/networks-both/{gID}.gml')
    ngraphs[gID] = nx.read_gml(f'./data/networks-nss/{gID}.gml')

number_of_runs = 100
edge_reverse_probs = [i/100 for i in range(1, 101, 1)]
subpaths = ['iss', 'nss', 'both']


results = []
print('reading sims...')
for graphs, subpath in zip([igraphs, ngraphs, bgraphs], subpaths):
    subreadpath = f'{readpath}/{subpath}'
    print(subpath)
        
    for g in list(graphs.items()):
        gidpath = f'{subreadpath}/{g[0]}'
        print(g[0])

        emp_leader = ia.find_leader(g[1], alpha = a, weight = 'weight')

        for prob in edge_reverse_probs:
            thisreadpath = f'{gidpath}/{int(prob * 100)}'
            random_leaders = []
            for run in range(number_of_runs):
                rg = nx.read_gml(f'{thisreadpath}/{int(run)}.gml')
                random_leader = ia.find_leader(rg, alpha = a, weight = 'weight')
                random_leaders.append(random_leader)

            try:
                prop = sum([rl == emp_leader for rl in random_leaders])/len(random_leaders)
            except:
                prop = np.nan
                
            result = {"cond": subpath, "gID": g[0], "prob": prob, "prop": prop}
            results.append(result)
results = pd.DataFrame(results)

print('reading sims complete')

summary = results.groupby(["cond", "prob"]).agg(
    mean_prop = ("prop", np.mean),
    med_prop = ("prop", np.median)
).reset_index()

plt.rcParams.update({'font.size': 16})
color0 = '#bfbfbf'
color1 = "#ee6a50"##"#b4eeb4"#'#698b69'
color2 = "#00688b"#'#8b7355'

fig, axs = plt.subplots(1, 3, figsize = (21, 7))
for cond, ax in zip(pd.unique(results["cond"]), axs.flatten()):
    for gID in gIDs:
        gr = results.loc[(results['gID'] == gID) & (results['cond'] == cond), ]
        ax.plot(gr['prob'], gr['prop'],
                c = color0, #'#586e75',
                linewidth = 1)#,
        #alpha = .5)
    ax.plot(summary.loc[summary['cond'] == cond, 'prob'],
            summary.loc[summary['cond'] == cond, 'mean_prop'],
            c = color1, linewidth = 3, label = 'Mean', zorder = 2.5)
    ax.plot(summary.loc[summary['cond'] == cond, 'prob'],
            summary.loc[summary['cond'] == cond, 'med_prop'],
            c = color2, linewidth = 3, label = 'Median')
    ax.set_xlabel('Pr(Reverse Edge)')

axs[0].set_ylabel('Pr(Leaders Match)')
axs[0].legend()
        
axs[0].text(0.005, 0.995, 'A', transform = axs[0].transAxes,# fontweight = "bold",
           horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')
axs[1].text(0.005, 0.995, 'B', transform = axs[1].transAxes,
           horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')
axs[2].text(0.005, 0.995, 'C', transform = axs[2].transAxes,
           horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')

fig.tight_layout()
#fig.savefig('./img/edge-randomization.svg')
fig.savefig('./img/edge-randomization.pdf')
fig.show()
