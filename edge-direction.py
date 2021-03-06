# This is a new version for the LQ submission.
# Analyze both the ISS-only and ISS/NSS interruption networks
# place the resulting figures side by side

# I am not satisfied with this simulation.
# I think what I want to do is draw random graphs until I reach 100 without change (if so, move to the next group) or, if change within 100, then continue until I have a sample of 30 random graphs.
# Then, base the proportion of differences on the mean/median of that sample of 30
# The problem is, some of these small graphs can only be different in so many ways.
# It may just be better to leave it as it is, pointing out that these are small, dense graphs that are resistant to randomization. 

# Interesting point: can say that it is possible this simulation reveals hidden heterogeneity in the stability of the emergent leader, but that conjecture for a future study. 

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import InterruptionAnalysis as ia

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

# Commented because fixed in the file that writes the GMLs in the first place
# # here need to put back in the nodes that aren't in i- and bgraphs.
# for gID in gIDs:
#     for node in vgraphs[gID].nodes:
#         for othergraph in [igraphs[gID], bgraphs[gID]]:
#             if node not in othergraph.nodes:
#                 othergraph.add_node(node,
#                                     gender = vgraphs[gID].nodes[node]['gender'],
#                                     tst = vgraphs[gID].nodes[node]['tst'])

number_of_runs = 100
edge_reverse_probs = [i/100 for i in range(1, 101, 1)]

all_results = []
all_summaries = []
print('running sim...')
for graphs, ax in zip([igraphs, ngraphs, bgraphs], axs.flatten()): # bgraphs
    #graphs = bgraphs
    results = []
    no_diffs = {gID: {} for gID in gIDs}
    for prob in edge_reverse_probs:
        print(prob)

        for g in list(graphs.items()):
        #for g in [g for g in list(graphs.items()) if g[0] == 'FLM']:
            emp_leader = ia.find_leader(g[1], alpha = a, weight = 'weight')

            random_leaders = []
            no_diff = 0
            for _ in range(number_of_runs):
                rg = ia.randomize_edge_directions(g[1], p = prob)
                if g[1].edges == rg.edges:
                    no_diff += 1
                    #continue # then need to see what's different about these groups
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
    # can start some analysis here. I want the group IDs for those groups whose networks don't change to compare against the groups that don't respond to edge randomization.
    all_results.append(results)
    
    summary = results.groupby(['prob']).agg([np.mean, np.median])
    all_summaries.append(summary)
    # no_diffs = pd.DataFrame(no_diffs)
    # wccsize = {}
    # sccwcc = {}
    # for gID in gIDs:
    #     ig = igraphs[gID]
    #     wccsize[gID] = len(ia.get_wcc(ig))/len(ig)
    #     sccwcc[gID] = len(ia.get_scc(ig))/len(ia.get_wcc(ig))

    # nd_data = pd.DataFrame({'wccsize': wccsize, 'sccwcc': sccwcc, 'num_nd': no_diffs.sum()})

print('sim complete')

plt.rcParams.update({'font.size': 12})
fig, axs = plt.subplots(1, 3, figsize = (21, 7))
for results, ax in zip(all_results, axs.flatten()):
    for gID in gIDs:
        gr = results[results['gID'] == gID]
        ax.plot(gr['prob'], gr['prop'],
                c = '#586e75',
                linewidth = 1,
                alpha = .5)#1 - nx.density(graphs[gID]))
        
for summary, ax in zip(all_summaries, axs.flatten()):
    ax.plot(list(summary.index), summary['prop']['mean'], c = '#268bd2', linewidth = 3, label = 'Mean')
    ax.plot(list(summary.index), summary['prop']['median'], c = '#859900', linewidth = 3, label = 'Median')
    ax.legend()

    #ax.axvline(x = .50, color = '#dc322f', alpha = .5)
    #ax.axhline(y = summary['prop'].loc[.5, 'mean'], color = '#dc322f', alpha = .5)
for ax in axs.flatten():
    ax.set_xlabel('Pr(Reverse Edge)')
    ax.set_ylabel('Pr(Leaders Match)')

axs[0].text(0.01, 0.99, 'A', transform = axs[0].transAxes,
           horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')
axs[1].text(0.01, 0.99, 'B', transform = axs[1].transAxes,
           horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')
axs[2].text(0.01, 0.99, 'C', transform = axs[2].transAxes,
           horizontalalignment = 'left', verticalalignment = 'top', fontsize = 'x-large')



fig.tight_layout()
fig.savefig('./img/edge-randomization.svg')
fig.show()
#plt.savefig('./img/rand-dir-sim.png')
#plt.show()
print('done')
