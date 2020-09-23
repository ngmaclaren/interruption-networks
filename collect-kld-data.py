import numpy as np
import pandas as pd

import InterruptionAnalysis as ia

data = pd.read_csv('./data/timeseries.csv', index_col = 0)

gIDs = pd.unique(data['gID'])

print('group level kl')
group_level_kl = {}
for gID in gIDs:
    dat = data[data['gID'] == gID]
    dhvg = ia.directed_horizontal_visibility_graph(dat)
    kl = ia.kullback_leibler_divergence(dhvg)
    group_level_kl[gID] = kl

print('indiv level kl')
pIDs = pd.unique(data['pID'])
indiv_level_kl = {}
for pID in pIDs:
    dat = data[data['pID'] == pID]
    if len(dat) < 2:
        continue
    else:
        dhvg = ia.directed_horizontal_visibility_graph(dat)
        kl = ia.kullback_leibler_divergence(dhvg)
        indiv_level_kl[pID] = kl


group_kld = pd.Series(group_level_kl)
indiv_kld = pd.Series(indiv_level_kl)

gID_count = data.groupby(['gID']).size()
pID_count = data.groupby(['pID']).size()

group_data = pd.DataFrame({'D': group_kld, 'n': gID_count})
indiv_data = pd.DataFrame({'D': indiv_kld, 'n': pID_count})

group_data.reset_index(level = 0, inplace = True)
indiv_data.reset_index(level = 0, inplace = True)
group_data = group_data.rename(columns = {'index': 'gID'})
indiv_data = indiv_data.rename(columns = {'index': 'pID'})

group_data.to_csv('./data/group-kld.csv')
indiv_data.to_csv('./data/indiv-kld.csv')

print('done')
