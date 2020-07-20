import numpy as np
import pandas as pd
import networkx as nx

import InterruptionAnalysis as ia

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
a = 0.99

gID = 'XOA'
pIDs = pd.unique(votedata[votedata['gID'] == gID]['pID'])
votecols = ['voteFor_p1', 'voteFor_p2', 'voteFor_p3', 'voteFor_p4', 'voteFor_p5']

C = ia.contact_sequence(data, pIDs)
ig = ia.interruption_network(C, pIDs)
vg = ia.vote_network(votedata, pIDs, votecols)

ig_pr = nx.pagerank_numpy(ig, alpha = a, weight = 'weight')
vg_pr = nx.pagerank_numpy(vg, alpha = a)

print('done')
