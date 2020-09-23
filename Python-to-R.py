import numpy as np
import pandas as pd
import networkx as nx

import InterruptionAnalysis as ia

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
surveydata = pd.read_csv('./data/speakingTime-data.csv', index_col = 0)
#surveydata['gID'] = surveydata.index.to_series().map(lambda x: str(x)[:3])

gIDs = pd.unique(votedata['gID'])
votecols = ['voteFor_p1', 'voteFor_p2', 'voteFor_p3', 'voteFor_p4', 'voteFor_p5']

for gID in gIDs:
    pIDs = pd.unique(data[data['gID'] == gID]['pID'])
    C = ia.contact_sequence(data, pIDs)
    ig = ia.interruption_network(C, pIDs)
    tsts = data[data['pID'].isin(pIDs)].groupby(['pID']).agg({'dur': sum})
    for pID in pIDs:
        ig.nodes[pID]['gender'] = surveydata.loc[pID, 'Gender']
        ig.nodes[pID]['tst'] = int(tsts.loc[pID, 'dur'])
    nx.write_gml(ig, f'./data/networks/{gID}.gml')

for gID in gIDs:
    pIDs = pd.unique(votedata[votedata['gID'] == gID]['pID'])
    vg = ia.vote_network(votedata, pIDs, votecols)
    for pID in pIDs:
        vg.nodes[pID]['gender'] = surveydata.loc[pID, 'Gender']
    nx.write_gml(vg, f'./data/votenets/{gID}.gml')

print('done')
