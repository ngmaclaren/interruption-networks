import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import InterruptionAnalysis as ia

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
surveydata = pd.read_csv('./data/speakingTime-data.csv', index_col = 0)
surveydata.set_index("pID", inplace = True)
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

    
i_wid = {} # ISS count
i_wod = {}
i_pr = {}
n_wid = {}
n_wod = {} # NSS count
n_pr = {}
b_wid = {}
b_wod = {}
b_pr = {}
t_wid = {}
t_wod = {}
v_wid = {} # total number of votes
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

surveydata = pd.concat([surveydata, idata], axis = 1)
# for convenience
surveydata["ISS"] = idata["i_wid"]
surveydata["NSS"] = idata["n_wod"]
surveydata["gID"] = surveydata["Group_ID"]
# turns is not exactly equal to t_wid
surveydata["Turns"] = data.groupby("pID")["begin"].count()

surveydata["Turns"].fillna(0, inplace = True)

surveydata["d_male"] = surveydata["Gender"].map({"male": 1, "female": 0})
surveydata["d_english"] = surveydata["English_Second_Language"].map({"no": 0, "yes": 1})
surveydata["d_simulation"] = surveydata["Simulation"].map({"bct": 0, "cs": 1})
surveydata["d_institution"] = surveydata["Institution"].map({"S1": 0, "S2": 1})
surveydata["d_operator"] = surveydata["Participant_Is_Operator"].map({"no": 0, "yes": 1})

surveydata["Total_Speaking_Time"] /= 1000

# this keeps all the columns from the ORB data and adds node level data from each network
# save this data to CSV for 2SLS analysis in Stata
surveydata.to_csv("./data/speakingTime-data-extended.csv")

cols = ['ISS', 'NSS', 'Turns', 'Total_Speaking_Time', 'Planning_Phase_Vote_Total',
        'd_male', 'Age', 'Game_Knowledge_Quiz', 'd_english', 'd_operator',
        'Conscientiousness', 'Agreeableness', 'Neuroticism', 'Openness', 'Extraversion',
        'Group_Size', 'd_simulation', 'd_institution']

dat = surveydata.loc[:, cols]
dat.corr("pearson").to_csv("./data/icorr.csv")

summarystats = dat.agg([np.mean, np.median, np.std, min, max]).T
dat["gID"] = surveydata["gID"]
summarystats["icc1"] = np.nan
for col in cols:
    X = dat[[col, "gID"]]
    model = f'{col} ~ gID'
    k = X.groupby("gID").count().mean()
    summarystats.loc[col, "icc1"] = ia.icc1(X, model, k)
summarystats.to_csv("./data/i-sumstats.csv")

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
gsumstat = gdata.agg([np.mean, np.median, np.std, min, max]).T
gsumstat.to_csv("./data/g-sumstat.csv")

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
