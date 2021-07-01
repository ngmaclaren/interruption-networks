import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as Color
import InterruptionAnalysis as ia

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
numericcols = ['begin', 'end', 'dur', 'lat']
data[numericcols] /= 1000 # cannot be done in seconds
votedata = pd.read_csv('./data/vote-data.csv')
stdata = pd.read_csv('./data/speakingTime-data.csv', index_col = 0)
a = 0.99

gIDs = pd.unique(data['gID'])
inets = {gID: nx.read_gml(f'./data/networks-iss/{gID}.gml') for gID in gIDs}
bnets = {gID: nx.read_gml(f'./data/networks-both/{gID}.gml') for gID in gIDs}

# for ease of .loc use:
for df in [votedata, stdata]:
    df.set_index('pID', inplace = True)

tst_dfs = []
iss_dfs = []
nss_dfs = []
for gID in gIDs:
    dat = data[data['gID'] == gID]
    pIDs = pd.unique(dat['pID'])

    iss = ia.interruptive_simultaneous_speech(data, pIDs)
    sorter = iss.groupby(['i', 'begin'])['dur'].count()
    sorter = list(sorter[sorter > 1].index)
    for row in sorter:
        temp = iss[(iss['i'] == row[0]) & (iss['begin'] == row[1])]
        drop = temp[temp['dur'] != max(temp['dur'])]
        iss.drop(list(drop.index), inplace = True)

    nss = ia.non_interruptive_simultaneous_speech(data, pIDs)
    sorter = nss.groupby(['j', 'begin'])['dur'].count()
    sorter = list(sorter[sorter > 1].index)
    for row in sorter:
        temp = nss[(nss['j'] == row[0]) & (nss['begin'] == row[1])]
        drop = temp[temp['dur'] != max(temp['dur'])]
        nss.drop(list(drop.index), inplace = True)

    tsts = pd.DataFrame({'tst_sum': dat.groupby('pID')['dur'].sum(),
                         'tst_count': dat.groupby('pID')['dur'].count()})
    isss = pd.DataFrame({'iss_sum': iss.groupby('i')['dur'].sum(),
                         'iss_count': iss.groupby('i')['dur'].count()})
    nsss = pd.DataFrame({'nss_sum': nss.groupby('j')['dur'].sum(),
                         'nss_count': nss.groupby('j')['dur'].count()})
    
    tst_dfs.append(tsts)
    iss_dfs.append(isss)
    nss_dfs.append(nsss)
    
tst = pd.concat(tst_dfs)
iss = pd.concat(iss_dfs); iss.index.rename('pID', inplace = True)
nss = pd.concat(nss_dfs); nss.index.rename('pID', inplace = True)
votes = votedata[['gID', 'totalVotes_p', 'totalVotes_g']].copy()
votes.rename(columns = {'totalVotes_p': 'votes', 'totalVotes_g': 'votes_g'}, inplace = True)

df = pd.concat([votes, tst, iss, nss], axis = 1)
df['gender'] = stdata['Gender']

results = []
for gID in gIDs:
    pIDs = pd.unique(data[data['gID'] == gID]['pID'])
    iss_pr = nx.pagerank_numpy(inets[gID], weight = 'weight', alpha = a)
    both_pr = nx.pagerank_numpy(bnets[gID], weight = 'weight', alpha = a)
    rows = [{'pID': pID, 'iss_pr': iss_pr[pID], 'both_pr': both_pr[pID]#, 'adj_pr': adj_pr[pID]
             } for pID in pIDs]
    results.extend(rows)
prs = pd.DataFrame(results)
prs.set_index('pID', inplace = True)
df = pd.concat([df, prs], axis = 1)

grpmaxvotes = pd.DataFrame(df.groupby('gID')['votes'].max())
grpmaxvotes.reset_index(inplace = True)

df['grpmaxvotes'] = df['gID'].apply(lambda x: grpmaxvotes.loc[grpmaxvotes['gID'] == x, 'votes'].item())
df['both'] = df['iss_count'] + df['nss_count']
df['both_rate'] = df['both']/df['tst_sum']

cols = ['tst_sum', 'tst_count', 'iss_sum', 'iss_count', 'nss_sum', 'nss_count', 'iss_pr', 'both_pr']
df[cols] = df[cols].fillna(0)

df['iss_rate'] = df['iss_count']/df['tst_sum']

plt.rcParams['font.size'] = 14
fig, ax = plt.subplots(figsize = (7, 7))

color0 = Color.to_rgba('#bfbfbf')
color1 = Color.to_rgba("#ee6a50")
color2 = Color.to_rgba("#00688b")

ax.set_xlabel('TST')
ax.set_ylabel('ISS Rate')

ldats = []
odats = []
for i in grpmaxvotes.index:
    ldat = df[(df['gID'] == grpmaxvotes.loc[i, 'gID']) &
             (df['votes'] != grpmaxvotes.loc[i, 'votes'])]
    odat = df[(df['gID'] == grpmaxvotes.loc[i, 'gID']) &
             (df['votes'] == grpmaxvotes.loc[i, 'votes'])]
    ldats.append(ldat)
    odats.append(odat)
ldats = pd.concat(ldats)
odats = pd.concat(odats)

fc1 = color0[:3] + (.5, )
ec1 = color0
fc2 = color1[:3] + (.5, )
ec2 = color1

ax.scatter(ldats.tst_sum, ldats.iss_rate, marker = 'o',
              color = fc1, edgecolors = ec1, s = 100, linewidths = 2.5,
              label = "Group members with less than\nthe maximum number of votes")
ax.scatter(odats.tst_sum, odats.iss_rate, marker = 's',
              color = fc2, edgecolors = ec2, s = 100, linewidths = 2.5,
              zorder = 2.5,
              label = "Group members with the most votes\n(including ties) in each group")
ax.legend(fontsize = 12)

fig.savefig("./img/tst-iss-rate-fig.pdf")
fig.show()

fig.show()
