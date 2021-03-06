import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import InterruptionAnalysis as ia

# Friday, February 12, 2021
# Don't report Big 5, simulation, institution... I kind of think here I probably shouldn't report anything that's in the speaking time article.
# Variable list should be: 
# - weighted, non-normalized in/out deg
# - PageRank as appropriate
# - the counts: ISS, NSS

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
surveydata = pd.read_csv('./data/all-surveys-calc-anon2.csv', index_col = 0)
a = 0.99

gIDs = pd.unique(data['gID'])
igraphs = {}
ngraphs = {}
bgraphs = {}
tgraphs = {}
vgraphs = {}
for gID in gIDs:
    igraphs[gID] = nx.read_gml(f'./data/networks-iss/{gID}.gml')
    ngraphs[gID] = nx.read_gml(f'./data/networks-nss/{gID}.gml')
    bgraphs[gID] = nx.read_gml(f'./data/networks-both/{gID}.gml')
    tgraphs[gID] = nx.read_gml(f'./data/turnnets/{gID}.gml')
    vgraphs[gID] = nx.read_gml(f'./data/votenets/{gID}.gml')

i_wid = {}
i_wod = {}
i_pr = {}
n_wid = {}
n_wod = {}
n_pr = {}
b_wid = {}
b_wod = {}
b_pr = {}
t_wid = {}
t_wod = {}
v_wid = {}
v_wod = {}
v_pr = {}

for gID in gIDs:
    ig = igraphs[gID]
    ng = ngraphs[gID]
    bg = bgraphs[gID]
    tg = tgraphs[gID]
    vg = vgraphs[gID]
    iprs = nx.pagerank_numpy(ig, weight = 'weight', alpha = a)
    nprs = nx.pagerank_numpy(ng, weight = 'weight', alpha = a)
    bprs = nx.pagerank_numpy(bg, weight = 'weight', alpha = a)
    vprs = nx.pagerank_numpy(vg, alpha = a)

    for pID in pd.unique(votedata[votedata['gID'] == gID]['pID']):
        i_wid[pID] = ig.in_degree(weight = 'weight')[pID]
        i_wod[pID] = ig.out_degree(weight = 'weight')[pID]
        i_pr[pID] = iprs[pID]
        n_wid[pID] = ng.in_degree(weight = 'weight')[pID]
        n_wod[pID] = ng.out_degree(weight = 'weight')[pID]
        n_pr[pID] = nprs[pID]
        b_wid[pID] = bg.in_degree(weight = 'weight')[pID]
        b_wod[pID] = bg.out_degree(weight = 'weight')[pID]
        b_pr[pID] = bprs[pID]
        t_wid[pID] = tg.in_degree(weight = 'weight')[pID]
        t_wod[pID] = tg.out_degree(weight = 'weight')[pID]
        v_wid[pID] = vg.in_degree()[pID]
        v_wod[pID] = vg.out_degree()[pID]
        v_pr[pID] = vprs[pID]

idata = pd.DataFrame({
    'i_wid': i_wid, 'i_wod': i_wod, 'i_pr': i_pr,
    'n_wid': n_wid, 'n_wod': n_wod, 'n_pr': n_pr,
    'b_wid': b_wid, 'b_wod': b_wod, 'b_pr': b_pr,
    't_wid': t_wid, 't_wod': t_wod,
    'v_wid': v_wid, 'v_wod': v_wod, 'v_pr': v_pr
    })
idata['gID'] = idata.index.map(lambda x: str(x)[:3])
idata.to_csv('./data/idata-networks.csv')

cols = list(idata.drop(columns = 'gID'))
print('i mean')
print(idata[cols].mean())
print('i sd')
print(idata[cols].std())
print('i min')
print(idata[cols].min())
print('i max')
print(idata[cols].max())
print('i ICC1')
for col in cols:
    X = idata[[col, 'gID']]
    model = f'{col} ~ gID'
    k = X.groupby('gID').count().mean()
    print(f'{col}: {ia.icc1(X, model, k)}')
icorr = idata[cols].corr('pearson')
icorr.to_csv('./data/icorr-networks.csv')

# group
# So, perhaps one table for group-level variables (size, % female, total TST, avg TST, total ISS, total NSS, avg ISS, avg NSS, perhaps also total number of turns and average number of turns) and another table for network statistics. Can do the table as a letter-sectioned column of smaller tables with mean, sd, min, max, corr within for each type of network.
size = {}
p_female = {}
t_tst = {}
a_tst = {}
t_iss = {}
a_iss = {}
t_nss = {}
a_nss = {}
t_turns = {}
a_turns = {}

for gID in gIDs:
    ig = igraphs[gID]
    ng = ngraphs[gID]
    bg = bgraphs[gID]
    tg = tgraphs[gID]
    vg = vgraphs[gID]

    size[gID] = vg.order()

    genders = list(nx.get_node_attributes(vg, 'gender').values())
    p_female[gID] = len([v for v in genders if v == 'female'])/len(genders)

    tsts = list(nx.get_node_attributes(ig, 'tst').values())
    t_tst[gID] = sum(tsts)/1000
    a_tst[gID] = np.mean(tsts)/1000

    isss = list(dict(ig.in_degree(weight = 'weight')).values())
    t_iss[gID] = sum(isss)
    a_iss[gID] = np.mean(isss)

    nsss = list(dict(ng.out_degree(weight = 'weight')).values())
    t_nss[gID] = sum(nsss)
    a_nss[gID] = np.mean(nsss)

    dat = data[data['gID'] == gID]
    t_count = dat.groupby('pID')['dur'].count()
    t_turns[gID] = t_count.sum()
    a_turns[gID] = t_count.mean()

gdata = pd.DataFrame({
    'size': size, 'p_female': p_female,
    't_tst': t_tst, 'a_tst': a_tst,
    't_iss': t_iss, 'a_iss': a_iss, 't_nss': t_nss, 'a_nss': a_nss,
    't_turns': t_turns, 'a_turns': a_turns
    })
gdata.to_csv('./data/gdata-nodes.csv')

print('g mean')
print(gdata.mean())
print('g sd')
print(gdata.std())
print('g min')
print(gdata.min())
print('g max')
print(gdata.max())

gcorr = gdata.corr('pearson')
gcorr.to_csv('./data/gcorr-nodes.csv')

# density, centralization, avg clustering, maybe avg shortest path length for ISS, NSS, Both, Turn, Vote.

graphs = [igraphs, ngraphs, bgraphs, tgraphs, vgraphs]
for G, name in zip(graphs, ['i', 'n', 'b', 't', 'v']):
    p_isol = {}
    density = {}
    centralization = {}
    avg_clust = {}
    avg_shortest_path = {}

    for gID in gIDs:
        g = G[gID]
        n_isol = len([x for x in list(nx.weakly_connected_components(g)) if len(x) == 1])
        p_isol[gID] = n_isol/g.order()

        density[gID] = nx.density(g)
        if G == tgraphs:
            centralization[gID] = ia.in_degree_centralization(g, weight = 'weight')
        elif G == vgraphs:
            centralization[gID] = ia.pagerank_centralization(g, alpha = a)
        else:
            centralization[gID] = ia.pagerank_centralization(g, alpha = a, weight = 'weight')
        avg_clust[gID] = nx.average_clustering(g)
        avg_shortest_path[gID] = nx.average_shortest_path_length(ia.get_wcc(g))

    ndata = pd.DataFrame({
        'p_isol': p_isol, 'density': density, 'centralization': centralization,
        'avg_clust': avg_clust, 'avg_shortest_path': avg_shortest_path
        })
    
    print(f'{name} means')
    print(ndata.mean())
    print(f'{name} sd')
    print(ndata.std())
    print(f'{name} min')
    print(ndata.min())
    print(f'{name} max')
    print(ndata.max())

    print(ndata.corr('pearson'))






    
# tst = {}
# gender = {}
# wid = {}
# wod = {}
# pr = {}
# vid = {}
# vod = {}
# vpr = {}
# bpr = {}
# for gID in gIDs:
#     ig = igraphs[gID]
#     vg = vgraphs[gID]
#     bg = bgraphs[gID]
#     prs = nx.pagerank_numpy(ig, weight = 'weight', alpha = a)
#     vprs = nx.pagerank_numpy(vg, alpha = a)
#     bprs = nx.pagerank_numpy(bg, weight = 'weight', alpha = a)
#     for pID in pd.unique(votedata[votedata['gID'] == gID]['pID']):
#         if pID in ig.nodes:
#             tst[pID] = nx.get_node_attributes(ig, 'tst')[pID]
#             wid[pID] = ig.in_degree(weight = 'weight')[pID]
#             wod[pID] = ig.out_degree(weight = 'weight')[pID]
#             pr[pID] = prs[pID]
#             bpr[pID] = bprs[pID]
#         else:
#             tst[pID] = 0
#             wid[pID] = 0
#             wod[pID] = 0
#             pr[pID] = 0
#         gender[pID] = surveydata.loc[pID, 'gender']
        
#         vid[pID] = vg.in_degree()[pID]
#         vod[pID] = vg.out_degree()[pID]
#         vpr[pID] = vprs[pID]
# idata = pd.DataFrame({'tst': tst, 'gender': gender, 'wid': wid, 'wod': wod, 'pr': pr, 'vid': vid, 'vod': vod, 'vpr': vpr, 'bpr': bpr})
# idata['tst'] = idata['tst']/1000
# idata['gender'] = idata['gender'].map({'female': 0, 'male': 1})
# idata['gID'] = idata.index.map(lambda x: str(x)[:3])
# idata.to_csv('./data/idata.csv')

# # size, density, wcc/size, scc/wcc, centralization, avg clustering, avg shortest path length
# size = {}
# density = {}
# vdensity = {}
# wccsize = {}
# sccwcc = {}
# cent = {}
# vcent = {}
# avgclust = {}
# vavgclust = {}
# #aspl = {}
# #vaspl = {}
# for gID in gIDs:
#     ig = igraphs[gID]
#     vg = vgraphs[gID]
#     size[gID] = len(ig)
#     density[gID] = nx.density(ig)
#     vdensity[gID] = nx.density(vg)
#     wccsize[gID] = len(ia.get_wcc(ig))/len(ig)
#     sccwcc[gID] = len(ia.get_scc(ig))/len(ia.get_wcc(ig))
#     cent[gID] = ia.pagerank_centralization(ig, alpha = a,  weight = 'weight')
#     vcent[gID] = ia.pagerank_centralization(vg, alpha = a)
#     avgclust[gID] = nx.average_clustering(ig, weight = 'weight')
#     vavgclust[gID] = nx.average_clustering(vg)
# gdata = pd.DataFrame({'size': size, 'density': density, 'vdensity': vdensity, 'wccsize': wccsize, 'sccwcc': sccwcc, 'cent': cent, 'vcent': vcent, 'avgclust': avgclust, 'vavgclust': vavgclust})
# gdata.to_csv('./data/gdata.csv')

# cols = list(idata.drop(columns = 'gID'))

# print(idata[cols].corr('spearman'))
# print(gdata.corr('spearman'))

# # Still need to do individual and group level means, medians, stdevs, ranges, and icc1s.
# print('i means')
# print(idata[cols].apply(np.mean))
# print('i medians')
# print(idata[cols].apply(np.median))
# print('i stdevs')
# print(idata[cols].apply(np.std))
# print('i range')
# print(idata[cols].apply(min))
# print(idata[cols].apply(max))
# print('ICC1')
# for col in cols:
#     X = idata[[col, 'gID']]
#     model = f'{col} ~ gID'
#     k = X.groupby('gID').count().mean()
#     print(f'{col}: {ia.icc1(X, model, k)}')

# print('g means')
# print(gdata.apply(np.mean))
# print('g medians')
# print(gdata.apply(np.median))
# print('g stdevs')
# print(gdata.apply(np.std))
# print('g range')
# print(gdata.apply(min))
# print(gdata.apply(max))




