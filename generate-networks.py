import os
import numpy as np
import pandas as pd
import networkx as nx

import InterruptionAnalysis as ia

if not os.path.isdir('./data/networks-iss'):
    os.mkdir('./data/networks-iss')
if not os.path.isdir('./data/networks-nss'):
    os.mkdir('./data/networks-nss')
if not os.path.isdir('./data/networks-both'):
    os.mkdir('./data/networks-both')
if not os.path.isdir('./data/votenets'):
    os.mkdir('./data/votenets')
if not os.path.isdir('./data/turnnets'):
    os.mkdir('./data/turnnets')

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
surveydata = pd.read_csv('./data/speakingTime-data.csv', index_col = 0)
surveydata.set_index('pID', inplace = True)

gIDs = pd.unique(votedata['gID'])
votecols = ['voteFor_p1', 'voteFor_p2', 'voteFor_p3', 'voteFor_p4', 'voteFor_p5']

for gID in gIDs:
    print(gID)
    pIDs = pd.unique(votedata[votedata['gID'] == gID]['pID'])
    tsts = surveydata.loc[pIDs, 'Total_Speaking_Time']
    genders = surveydata.loc[pIDs, 'Gender']
    esls = surveydata.loc[pIDs, 'English_Second_Language']
    isops = surveydata.loc[pIDs, 'Participant_Is_Operator']
    ages = surveydata.loc[pIDs, 'Age']
    ints = surveydata.loc[pIDs, 'Intelligence']
    gks = surveydata.loc[pIDs, 'Game_Knowledge_Quiz']
    insts = surveydata.loc[pIDs, 'Institution']
    sims = surveydata.loc[pIDs, 'Simulation']
    conscs = surveydata.loc[pIDs, 'Conscientiousness']
    agrees = surveydata.loc[pIDs, 'Agreeableness']
    neurs = surveydata.loc[pIDs, 'Neuroticism']
    opens = surveydata.loc[pIDs, 'Openness']
    extras = surveydata.loc[pIDs, 'Extraversion']
    

    igi = ia.interruption_network_pandas(data, pIDs, use = 'iss')
    ign = ia.interruption_network_pandas(data, pIDs, use = 'nss')
    igb = ia.interruption_network_pandas(data, pIDs, use = 'both')
    vg = ia.vote_network(votedata, pIDs, votecols, self_loops = True)
    tg = ia.turn_based_network(data, pIDs)

    for pID in pIDs:
        for g in [igi, ign, igb, vg, tg]:
            g.nodes[pID]['gender'] = str(genders[pID])
            g.nodes[pID]['tst'] = int(tsts[pID])
            g.nodes[pID]['esl'] = str(esls[pID])
            g.nodes[pID]['isop'] = str(isops[pID])
            # g.nodes[pID]['age'] = ages[pID] # keeps np.nan
            try: # replaces np.nan with a string.
                g.nodes[pID]['age'] = int(ages[pID])
            except:
                g.nodes[pID]['age'] = str("NA")
            g.nodes[pID]['intel'] = float(ints[pID])
            g.nodes[pID]['gameknowl'] = int(gks[pID])
            g.nodes[pID]['inst'] = str(insts[pID])
            g.nodes[pID]['sim'] = str(sims[pID])
            g.nodes[pID]['consc'] = float(conscs[pID])
            g.nodes[pID]['agree'] = float(agrees[pID])
            g.nodes[pID]['neur'] = float(neurs[pID])
            g.nodes[pID]['open'] = float(opens[pID])
            g.nodes[pID]['extra'] = float(extras[pID])

    nx.write_gml(igi, f'./data/networks-iss/{gID}.gml')
    nx.write_gml(ign, f'./data/networks-nss/{gID}.gml')
    nx.write_gml(igb, f'./data/networks-both/{gID}.gml')
    nx.write_gml(vg, f'./data/votenets/{gID}.gml')
    nx.write_gml(tg, f'./data/turnnets/{gID}.gml')

print('done')
