import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import InterruptionAnalysis as ia

savepath = './data/edgedir-sim'
if not os.path.isdir(savepath):
    os.mkdir(savepath)

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


#all_results = []
#all_summaries = []
print('running sims...')
#for graphs, ax in zip([igraphs, ngraphs, bgraphs], axs.flatten()):
for graphs, subpath in zip([igraphs, ngraphs, bgraphs], subpaths):
    subsavepath = f'{savepath}/{subpath}'
    print(subpath)
    if not os.path.isdir(subsavepath):
        os.mkdir(subsavepath)
        
    #results = []
    #no_diffs = {gID: {} for gID in gIDs}
    for g in list(graphs.items()):
        gidpath = f'{subsavepath}/{g[0]}'
        print(g[0])
        if not os.path.isdir(gidpath):
            os.mkdir(gidpath)

        for prob in edge_reverse_probs:
            for run in range(number_of_runs):
                rg = ia.randomize_edge_directions(g[1], p = prob)
                ns.write_gml(rg, f'{gidpath}/p{int(prob * 100)}r{int(run)}.gml', stringizer = str)



                
    for prob in edge_reverse_probs:
        print(prob)

        for g in list(graphs.items()):
            emp_leader = ia.find_leader(g[1], alpha = a, weight = 'weight')

            random_leaders = []
            no_diff = 0
            for run in range(number_of_runs):
                rg = ia.randomize_edge_directions(g[1], p = prob)
                #nx.write_gml(rg, f'{thissavepath}/{g[0]}{run}.gml', stringizer = str)
                if g[1].edges == rg.edges:
                    no_diff += 1
                random_leader = ia.find_leader(rg, alpha = a, weight = 'weight')
                random_leaders.append(random_leader)

            try:
                prop = sum([rl == emp_leader for rl in random_leaders])/len(random_leaders)
            except:
                prop = np.nan
            result = (g[0], prob, prop)
            no_diffs[g[0]][prob] = no_diff
            results.append(result)

    results = pd.DataFrame(results, columns = ['gID', 'prob', 'prop'])
    all_results.append(results)
    
    summary = results.groupby(['prob']).agg([np.mean, np.median])
    all_summaries.append(summary)

print('sim complete')

## Old version of plot
# plt.rcParams.update({'font.size': 14})
# fig, axs = plt.subplots(1, 3, figsize = (21, 7))
# for results, ax in zip(all_results, axs.flatten()):
#     for gID in gIDs:
#         gr = results[results['gID'] == gID]
#         ax.plot(gr['prob'], gr['prop'],
#                 c = '#bfbfbf', #'#586e75',
#                 linewidth = 1)#,
#                 #alpha = .5)
        
# for summary, ax in zip(all_summaries, axs.flatten()):
#     ax.plot(list(summary.index), summary['prop']['mean'], c = '#ee6a50', linewidth = 3, label = 'Mean')  # '#268bd2'
#     ax.plot(list(summary.index), summary['prop']['median'], c = '#5f9ea0', linewidth = 3, label = 'Median') # '#859900'
#     ax.legend()

#     #ax.axvline(x = .50, color = '#dc322f', alpha = .5)
#     #ax.axhline(y = summary['prop'].loc[.5, 'mean'], color = '#dc322f', alpha = .5)
# for ax in axs.flatten():
#     ax.set_xlabel('Pr(Reverse Edge)')
#     ax.set_ylabel('Pr(Leaders Match)')

# axs[0].text(0.01, 0.99, 'A', transform = axs[0].transAxes,
#            horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')
# axs[1].text(0.01, 0.99, 'B', transform = axs[1].transAxes,
#            horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')
# axs[2].text(0.01, 0.99, 'C', transform = axs[2].transAxes,
#            horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')



# fig.tight_layout()
# #fig.savefig('./img/edge-randomization.svg')
# fig.savefig('./img/edge-randomization.pdf')
# fig.show()
