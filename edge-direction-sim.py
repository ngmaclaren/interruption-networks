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


print('running sims...')
for graphs, subpath in zip([igraphs, ngraphs, bgraphs], subpaths):
    subsavepath = f'{savepath}/{subpath}'
    print(subpath)
    if not os.path.isdir(subsavepath):
        os.mkdir(subsavepath)
        
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
