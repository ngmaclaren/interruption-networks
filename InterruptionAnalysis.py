import random
import itertools
import collections
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import networkx as nx
import statsmodels.api as sm
from statsmodels.formula.api import ols

import scipy.stats as stats
from scipy.stats._continuous_distns import _distn_names
from scipy.special import rel_entr
import statsmodels.api as sm

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def collect_vocalization_sequence(datafile, loc = './data/diarizations/', ext = 'pl.eaf'):
    """
    Given the ELAN diarization files, this function returns a pd.DataFrame with two cols, 'start' and 'dur'. Two further cols, 'gID' and 'pID', identify the speakers.

    This is the primary function that turns a diarization sequence (in the ELAN format used in this study) into time series data that can be further analyzed.
    """
    tree = ET.parse(datafile)
    root = tree.getroot()

    time_order = root.find('TIME_ORDER')
    start_time = int(time_order[1].get('TIME_VALUE'))
    participants = root.findall('TIER')[2:]

    gID = datafile.replace(loc, '')
    gID = gID.replace(ext, '')
        
    data = []
    for participant in participants:
        pID = gID + participant.get('TIER_ID')[4:]

        for u in participant.iter('ALIGNABLE_ANNOTATION'):
            ts1 = u.get('TIME_SLOT_REF1')
            ts2 = u.get('TIME_SLOT_REF2')
            begin, end = None, None
            for ts in time_order.iter('TIME_SLOT'):
                if ts.get('TIME_SLOT_ID') == ts1:
                    begin = int(ts.get('TIME_VALUE'))
            for ts in time_order.iter('TIME_SLOT'):
                if ts.get('TIME_SLOT_ID') == ts2:
                    end = int(ts.get('TIME_VALUE'))
            begin = begin - start_time
            end = end - start_time
            dur = end - begin
            data.append( (gID, pID, begin, end, dur) )
        
        
    data = pd.DataFrame(data, columns = ['gID', 'pID', 'begin', 'end', 'dur'])
    
    return data

########################################
#### Time Series Analysis Functions ####
########################################

def lines_transform(data):
    data = data.sort_values(['begin'])
    lines = [[(i, 0), (i, j)] for i, j in zip(data.begin, data.dur)]
    return lines

def directed_horizontal_visibility_graph(data):
    """
    From Lacasa et al 2012. Data needs to have data['begin'] and data['dur']. Nodes will be named by the start time of the speaking event.
    """
    edges = []
    lines = lines_transform(data) #[[(i, 0), (i, j)] for i, j in zip(data.begin, data.dur)]
    
    for i in lines:
        ix = i[0][0]
        iy = i[1][1]
        for j in lines:
            jx = j[0][0]
            jy = j[1][1]
            if jx < ix:
                continue
            if jx == ix:
                continue
            ns = [n for n in lines if n[0][0] > ix and n[0][0] < jx]
            if ns:
                maxn = max(ns, key = lambda x: x[1][1])[1][1]
                if iy > maxn and jy > maxn:
                    edges.append((ix, jx))
            else:
                edges.append((ix, jx))

    g = nx.DiGraph()
    g.add_edges_from(edges)

    return g

def kullback_leibler_divergence(g):
    """
    Given a directed horizontal visibility graph (doesn't currently test for that...), return the Kullback-Leibler divergence between the out-degree probability distribution (forward in time) and the in-degree probability distribution (backward in time). This function is specifically for dealing with directed horizontal visibility graphs and doesn't generalize. Keeping the name for now for backwards compatibility. The more general function, based on SciPy's rel_entropy(), is get_kld().
    """
    in_degree_sequence = sorted([d for n, d in g.in_degree])
    in_degree_count = dict(collections.Counter(in_degree_sequence))
    in_deg, in_deg_count = zip(*in_degree_count.items())
    in_deg_prop = [i/sum(in_deg_count) for i in in_deg_count]
    
    out_degree_sequence = sorted([d for n, d in g.out_degree])
    out_degree_count = dict(collections.Counter(out_degree_sequence))
    out_deg, out_deg_count = zip(*out_degree_count.items())
    out_deg_prop = [i/sum(out_deg_count) for i in out_deg_count]
    
    len_diff = len(out_deg) - len(in_deg)
    zvec = list(np.zeros(abs(len_diff)))

    if len_diff > 0:
        in_deg_prop.extend(zvec)
    elif len_diff < 0:
        out_deg_prop.extend(zvec)
        
    D = sum([out_deg_prop[k] * np.log(out_deg_prop[k]/in_deg_prop[k]) for k in range(len(out_deg)) if in_deg_prop[k] > 0])

    return D

def dhvg_local_clustering(g, v, direction):
    """
    Ref: Donges et al (2013)
    """
    # this only works because the node labels are integers whose order makes sense
    if direction == "forward":
        nbr = sorted(list(g.successors(v)))
    elif direction == "backward":
        nbr = sorted(list(g.predecessors(v)))

    deg = len(nbr)

    if deg < 2:
        c = 0
        return c
    else:
        count = 0
        for j in nbr:
            jnbr = list(g.successors(j))
            for k in nbr:
                if k in jnbr:
                    count += 1

        denominator = comb(deg, 2)
        c = count/denominator
        return c


###################################
#### Markov Analysis Functions ####
###################################

def get_transition_matrix(data, x, ignore_simultaneous = False): #where ~x~ is a pID ####### should be ~r~ to be consistent
    """
    For a given data frame ~data~, return the Markov transition probability matrix associated with group member ~x~. ~data~ should have the columns 'pID', 'gID', 'begin', 'end', and 'dur'. 
    Let u_X be the number of distinct speaking events coded for participant X.
    Let A represent the event that participant X is not talking (state 'o' in the ABM) and B represent the event that participant X is talking (state 'e'). 
    We have a matrix P_X for each participant the values of which are represented by 
from\to  A      B
A        P(A|A) P(B|A)
B        P(A|B) P(B|B)
which translates to:

  AA AB  -->  [0, 0]  [0, 1]
  BA BB       [1, 0]  [1, 1]
    If ignore_simultaneous = True, remove all interruptive and non-interruptive simultaneous speech from consideration.
    """
    gID = list(data[data['pID'] == x]['gID'])[0]
    grpdata = data[data['gID'] == gID]
    X = data[data['pID'] == x]

    if ignore_simultaneous:
        pIDs = pd.unique(grpdata['pID'])
        interruptions = interruptive_simultaneous_speech(data, pIDs)
        interjections = non_interruptive_simultaneous_speech(data, pIDs)
        ξ = interruptions[interruptions['i'] == x]
        ζ = interjections[interjections['j'] == x]
        X = X[~(X['begin'].isin(ξ['begin'])) & ~(X['begin'].isin(ζ['begin']))]

    if X.empty:
        P = np.full((2, 2), np.nan)
        return P
                                        # These are actually correct, according to the P(B|A) notation
                                        # Downstream files use the translation to AA, AB, etc.
                                        # and the accompanying NumPy array slicing. 
    u_X = len(X)
    N_B = sum(X['dur'])
    N_A = max(grpdata['end']) - N_B
    P_AB = u_X/N_B
    P_BA = u_X/N_A
    P_BB = 1 - P_AB
    P_AA = 1 - P_BA

    P = np.array([[P_AA, P_BA],
                  [P_AB, P_BB]])
    return P

def get_threestate_transition_matrix(data, x):
    """
    For a given data frame ~data~, return the three-state Markov transition probability matrix associated with group member ~pID~ (~x~). ~data~ should have the columns 'pID', 'gID', 'start', 'end', and 'dur'. 
    Let A represent the event that no participant is talking (x is "listening" to silence), B represent the event that participant x is not talking but at least one participant x_i is talking (x is "listening" to one or more x_i talking), and C represent the event that partixipant x is talking (x is talking).
    We have a matrix P_X for each participant the values of which are represented by 
from\to  A      B      C
A        P(A|A) P(B|A) P(C|A)
B        P(A|B) P(B|B) P(C|B)
C        P(A|C) P(B|C) P(C|C)
    """
    gID = list(data[data['pID'] == x]['gID'])[0]
    grpdata = data[data['gID'] == gID].sort_values('begin', ignore_index = True)
    pIDs = pd.unique(grpdata['pID'])
    X = data[data['pID'] == x]
    nox = grpdata[-(grpdata['pID'].isin([x]))]

    interruptions = interruptive_simultaneous_speech(data, pIDs)
    interjections = non_interruptive_simultaneous_speech(data, pIDs)
    ξ = interruptions[interruptions['i'] == x]
    ζ = interjections[interjections['j'] == x]
    ξp = interruptions[interruptions['j'] == x] # someone interrupts me
    ζp = interjections[interjections['i'] == x] # someone interjects me

    nox_ξ = interruptions[-(interruptions['i'].isin([x]))]
    nox_ζ = interjections[-(interjections['j'].isin([x]))]

    T = max(grpdata['end'])

    u_x = len(X)
    u_nox = len(nox)
    u_ξ = len(ξ)
    u_ζ = len(ζ)
    u_ξp = len(ξp)
    u_ζp = len(ζp)
    u_nox_ξ = len(nox_ξ)
    u_nox_ζ = len(nox_ζ)
    N_C = sum(X['dur'])
    #N_B = sum(nox['dur'])
    N_B = 0
    for i in range(len(nox)): # i think this needs to be nox
        if i == 0:
            N_B += nox.loc[i, 'dur']
        else:
            before = nox[nox.index < i]
            if nox.loc[i, 'begin'] < max(before['end']):
                overlap = max(before['end']) - nox.loc[i, 'begin']
                N_B += nox.loc[i, 'dur'] - overlap
            else:
                N_B += nox.loc[i, 'dur']

    N_A = T - N_B - N_C

    # some of these are yielding numbers outside [0,1]
    #P_CA = (u_x - u_ξ - u_ζ)/N_A
    P_CA = len(
        pd.unique(X[
            -(X['begin'].isin(ξ['begin']) |
              X['begin'].isin(ζ['begin']))
        ]['begin']
                  ))/N_A
    #P_BA = (u_nox - u_nox_ξ - u_nox_ζ)/N_A
    P_BA = len(
        pd.unique(nox[
            -(nox['begin'].isin(nox_ξ['begin']) |
              nox['begin'].isin(nox_ζ['begin']))
        ]['begin']
                  ))/N_A
    P_AA = 1 - (P_CA + P_BA)

    #P_CB = (u_ξ + u_ζ)/N_B
    P_CB = (len(pd.unique(ξ['begin'])) +
            len(pd.unique(ζ[
                -(ζ['begin'].isin(ξ['begin']))
            ]['begin']))
            )/N_B
    #P_AB = (u_nox - u_nox_ξ - u_nox_ζ)/N_B # this is wrong --> (u_ξ + u_ζ + u_ξp + u_ζp)
    P_AB = len(
        pd.unique(nox[
            -(nox['end'].isin(ζ['end']))]['end']
                  ))/N_B
    P_BB = 1 - (P_CB + P_AB)

    #P_BC = u_ζ/N_C
    P_BC = (len(pd.unique(ζ['begin'])) + len(pd.unique(ξp['begin'])))/N_C # ζp['begin']
    #P_AC = (u_x - u_ζ)/N_C # this will miss any cases where I start speaking during someone's turn, interrupt them, they stop, but someone else is still speaking when I stop
    P_AC = len(X[-(X['end'].isin(ζ['end']))])/N_C
    P_CC = 1 - (P_BC + P_AC)

    P = np.array([
        [P_AA, P_BA, P_CA],
        [P_AB, P_BB, P_CB],
        [P_AC, P_BC, P_CC]])
    
    return P

#def otherstate_dependent_markov(data, pID): # this could probably be cleaned up a little now
def get_twolevel_transition_matrix(data, x):
    
    """
    For a given data frame ~data~, return the Markov transition probability matrices associated with group member ~x~. ~data~ should have the columns 'pID', 'gID', 'start', 'end', and 'dur'. In this case, there is a hard-coded pair of state sets: a 'perceptual' state of 'someone else is talking/not talking' and a 'self' state of 'talking/not talking'.
    """
    
    gID = list(data[data['pID'] == x]['gID'])[0]
    grpdata = data[data['gID'] == gID].sort_values('begin', ignore_index = True)
    pIDs = pd.unique(grpdata['pID']) 
    X = data[data['pID'] == x]
    nox = grpdata[-(grpdata['pID'].isin([x]))]#data[data['pID'] != pID]

    interruptions = interruptive_simultaneous_speech(data, pIDs)
    interjections = non_interruptive_simultaneous_speech(data, pIDs)
    ξ = interruptions[interruptions['i'] == x]
    ζ = interjections[interjections['j'] == x]
    ξp = interruptions[interruptions['j'] == x]
    ζp = interjections[interjections['i'] == x]#overlap = ξ['dur'].sum() + ζ['dur'].sum() + ξp['dur'].sum()
    nox_ξ = interruptions[-(interruptions['i'].isin([x]))]
    nox_ζ = interjections[-(interjections['j'].isin([x]))]

    T = max(data['end'])

    N_B = X['dur'].sum()
    N_A = T - N_B
    # these are wrong, below:
    #N_D = nox['dur'].sum()#- (interruptions[~interruptions.isin([pID])]['dur'].sum() + interjections[~interjections.isin([pID])]['dur'].sum())
    #N_C = T - N_D

    N_D = 0
    for i in range(len(grpdata)):
        if i == 0:
            N_D += grpdata.loc[i, 'dur']
        else:
            before = grpdata[grpdata.index < i]
            if grpdata.loc[i, 'begin'] < max(before['end']):
                overlap = max(before['end']) - grpdata.loc[i, 'begin']
                N_D += grpdata.loc[i, 'dur'] - overlap
            else:
                N_D += grpdata.loc[i, 'dur']
    N_C = T - N_D
    print(x)
    print(ξ)
    print(ζ)
    print([N_A, N_B, N_C, N_D])

    P_ABC = np.nan
    P_BBC = np.nan

    P_BAC = len(X[-(X['begin'].isin(ξ['begin']) | X['begin'].isin(ζ['begin']))]['begin']) / N_A
    P_AAC = 1 - P_BAC

    P_BAD = len(pd.unique(ξ['begin']))/(N_D - X['dur'].sum())
    P_AAD = 1 - P_BAD

    P_ABD = len(pd.unique(ζ['begin']))/(X['dur'].sum())
    P_BBD = 1 - P_ABD

    # Monday, November 23, 2020
    # there are negative numbers in BAC, ABC, BBC, and BBD
    # I believe this is because there are double counts in ξ and ζ
    # The trick I used above was to use pd.unique.
    # Study the threestate solution, and implement something similar that makes sense here.
    # P_BAC = (len(x) - (len(ξ) + len(ζ)))/(T - N_B - N_D)#(N_C - N_B)
    # P_AAC = 1 - P_BAC
    #  try:
    #      P_ABC = (len(x) - len(ξp))/(N_B - overlap)#(N_C + N_B - (ξ['dur'].sum() + ζ['dur'].sum()))
    #      P_BBC = 1 - P_ABC
    #  except:
    #      P_ABC = np.nan
    #      P_BBC = np.nan

    PC = np.array([
        [P_AAC, P_BAC],
        [P_ABC, P_BBC]
    ])

    # P_BAD = (len(ξ) + len(ζ))/(N_D - overlap)#(N_D - (ξ['dur'].sum() + ζ['dur'].sum()))
    # P_AAD = 1 - P_BAD
    # try:
    #     P_ABD = len(ξp)/(overlap)#(N_D + (ξ['dur'].sum() + ζ['dur'].sum()))
    #     P_BBD = 1 - P_ABD
    # except:
    #     P_ABD = np.nan
    #     P_BBD = np.nan

    PD = np.array([
        [P_AAD, P_BAD],
        [P_ABD, P_BBD]
    ])

    return [PC, PD]

# def g_n_dependent_markov(data, pID):
     # the idea here is going to be to make the calculations for otherstate_dependent_markov(), above, be recalculated for each other g_n. There will be lots of NA values. For this function and the one above, the except clause shouldn't assign P_ABD as np.nan, but rather as the corresponding value for no one talking, P_ABC---right? Or just shut off? No, I think P_ABC, because it's not that the interruption couldn't happen, it's that we didn't observe it.
#     pass
 
############################################
#### General Network Analysis Functions ####
############################################

def get_scc(g):
    """
    Given a graph, /g/, return the largest strongly connected component of that graph. 

    Function works on the assumption that there is only one relevant largest strongly connected component of the graph. Function will not work as expected in a larger graph with multiple strongly connected components, which may be of theoretical interest at some future time.
    """
    sccs = list(nx.strongly_connected_components(g))
    largest_scc = list(sorted(sccs, key = lambda x: len(x), reverse = True)[0])
    scc = g.subgraph(largest_scc) # may not always work if there are sometimes more than one scc
    return scc

def get_wcc(g):
    """
    Given a graph, /g/, return the largest weakly connected component of that graph. 

    Function works on the assumption that there is only one relevant largest weakly connected component of the graph. Function will not work as expected in a larger graph with multiple strongly connected components, which may be of theoretical interest at some future time.
    """
    wccs = list(nx.weakly_connected_components(g))
    largest_wcc = list(sorted(wccs, key = lambda x: len(x), reverse = True)[0])
    wcc = g.subgraph(largest_wcc)
    return wcc

# def indegree_centralization(g):
#     """
#     Find the in-degree centralization for a graph. Assumes a DiGraph with an edge attribute called 'weight'.
#     Following Sauer & Kauffeld, who use the standard formula.
#     """
#     idc = dict(g.in_degree(weight = 'weight')).values()
#     max_idc = max(idc)
#     numerator = sum([max_idc - idc[i] for i in idc])
#     denominator = 
#     C = 
#     return C
    

#################################################
#### Interruption Network Analysis Functions ####
#################################################

def find_leader(g, alpha = 0.85, weight = None):
    """
    Finds the 'leader' in a given graph, /g/, where the leader is defined as the node with the highest PageRank centrality in g.
    """
    pr = nx.pagerank_numpy(g, alpha = alpha, weight = weight)
    leader = max(pr, key = pr.get)
    return leader

def pagerank_centralization(g, alpha = 0.85, weight = None):
    """
    This function returns the centralization, currently hard-coded to PageRank, for a given graph according to the following formula:
    
    \sum (max(PR_k) - PR_{ik}) / same calculation for a star graph of the same size as k

    """
    h = nx.DiGraph()
    h.add_edges_from([(n, len(g)) for n in range(len(g))])
    pr_h = nx.pagerank_numpy(h, alpha = alpha)
    pr_h_max = list(sorted(pr_h.items(), key = lambda x: x[1], reverse = True))[0][1]

    pr = nx.pagerank_numpy(g, alpha = alpha, weight = weight)
    pr_max = list(sorted(pr.items(), key = lambda x: x[1], reverse = True))[0][1]

    H = sum([pr_h_max - pr_h[i] for i in pr_h.keys()])

    centralization = sum([pr_max - pr[i] for i in pr.keys()])/H

    return centralization

def in_degree_centralization(g, weight = None):
    """
    """
    h = nx.DiGraph()
    h.add_edges_from([(n, len(g)) for n in range(len(g))])
    idc_h = dict(h.in_degree())
    idc_h_max = max(idc_h.items(), key = lambda x: x[1])[1]

    idc = dict(g.in_degree(weight = weight))
    idc_max = max(idc.items(), key = lambda x: x[1])[1]

    H = sum([idc_h_max - idc_h[i] for i in idc_h.keys()])

    centralization = sum([idc_max - idc[i] for i in idc.keys()])/H

    return centralization
################################
#### Other Helper Functions ####
################################

def correlation_ratio(X, obs, group):
    """
    Within-groups eta^2. Given data (pd.DataFrame) X, where obs is the column of observations and group is the column of groups. obs and group should be str column names.
    https://en.wikipedia.org/wiki/Correlation_ratio
    """
    grouped = X[[obs, group]].groupby([group])

    group_n = grouped.count()
    group_means = grouped.mean()
    grand_mean = X[[obs]].mean()

    numerator = np.sum(group_n * (group_means - grand_mean)**2)
    denominator = np.sum((X[[obs]] - grand_mean)**2)

    eta2 = numerator/denominator

    return eta2

def cv(x):
    """
    Calculate the coefficient of variation for a 1d array. 

    Example of use:

    cv([leader_distance_tonext(graphs[sID]) for sID in graphs])
    # 0.8218378733935942
    cv([leader_distance_tomean(graphs[sID]) for sID in graphs])
    # 0.3429549155096544
    
    """
    m = np.mean(x)
    sd = np.std(x)
    cv = sd/m

    return cv

def icc1(X, model, k):
    """
    Given a pandas data frame that contains a vector of data and a vector of values representing a grouping variable, a model of the form given below, and an average group size, give Bliese's ICC1 (Bliese, 2000) for the data frame.

    Example of use:
    icc_data = {sID: list(all_pr[sID].values()) for sID in all_pr.keys()}
    
    all_pr = {sID: nx.pagerank(graphs[sID]) for sID in graphs.keys()}
    pr_data = pd.DataFrame(all_pr)
    pr_data['pr'] = pr_data.sum(axis = 1)
    
    def get_sID(x):
        sID = x[:4] if len(x) == 6 else x[:5]
        return sID
    
    pr_data['sID'] = pr_data.index.map(get_sID)

    pr_data = pr_data[['pr', 'sID']]

    X = pr_data
    model = 'pr ~ sID'
    k = X.groupby('sID').count().mean()
    
    print(icc1(X, model, k))
    """
    icc_ols = ols(model, data = X).fit()
    icc_anova = sm.stats.anova_lm(icc_ols)

    MSR = icc_anova.sum_sq[0]/icc_anova.df[0]
    MSE = icc_anova.sum_sq[1]/icc_anova.df[1]
    
    ICC1 = (MSR - MSE) / (MSR + k*MSE)

    return ICC1.item()

def get_kld(dat, ref, dist_name, bins = 75):
    """
    Find the Kullback-Leibler divergence from ~ref~ to ~dat~ using SciPy's rel_entr() function. A hypothesized distribution is required, as this function fits bot the ref and the data to the same probability distribution (fitted separately): ~dist_name~ should be one of SciPy's probability distributions (https://docs.scipy.org/doc/scipy/reference/stats.html).
    """
    dist = getattr(stats, dist_name)
    
    y, x = np.histogram(ref, bins = bins)
    x = (x + np.roll(x, -1))[:-1] / 2.0
    #y = y/np.sum(y)

    d_params = dist.fit(dat)
    d_args = d_params[:-2]
    d_loc = d_params[-2]
    d_scale = d_params[-1]
    dpdf = dist.pdf(x, loc = d_loc, scale = d_scale, *d_args)

    r_params = dist.fit(ref)
    r_args = r_params[:-2]
    r_loc = r_params[-2]
    r_scale = r_params[-1]
    rpdf = dist.pdf(x, loc = r_loc, scale = r_scale, *r_args)
       
    # dw = dist.fit(dat) # returns params
    # rw = dist.fit(ref)

    # dpdf = dist.pdf(x, c = dw[0], loc = dw[1], scale = dw[2])
    # rpdf = dist.pdf(x, c = rw[0], loc = rw[1], scale = rw[2])

    dy = dpdf/np.sum(dpdf)
    ry = rpdf/np.sum(rpdf)

    return sum(rel_entr(dy, ry))


#####################################
#### Network Generator Functions ####
#####################################

# Note for all interruption networks the 'j' and 'i' may seem backwards: they are based on the idea that 'i' interrupts 'j', so the edge in an interruption goes from 'j' to 'i'. "Unsuccessful" interruptions were labeled accordingly: the edge still goes from 'j' to 'i', but the action now also goes from 'j' to 'i'.
# TODO: add the disconnected... syntax to each network constructor.

def contact_sequence(data, pIDs):
    """
    Given a pd.DataFrame with variables ['gID', 'pID', 'begin', 'end'], generate a contact sequence of events (j, i, t) where j is the group member being interrupted, i is the group member doing the interrupting, and t is the start time of the interruptive speech. Ignores non-interruptive simultaneous speech.
    Ref: Holme & Saramäki (2012)
    """
    C = []

    data = data[data['pID'].isin(pIDs)]
    for i in range(len(data)):
        ego = data.iloc[i, :]
        interruptions = data[
            (data['end'] > ego['end']) &
            (data['begin'] < ego['end']) &
            (data['begin'] >= ego['begin'])
        ]
        
        if interruptions.empty:
            continue
        else:
            for j in range(len(interruptions)):
                alter = interruptions.iloc[j, :]
                c = (ego['pID'], alter['pID'], alter['begin'])
                C.append(c)
    C = sorted(C, key = lambda x: x[2])
    return C

def interruptive_simultaneous_speech(data, pIDs):
    """
    Given a pd.DataFrame with variables ['gID', 'pID', 'begin', 'end'], generate an interval sequnce of events (j, i, t, t') where j is the group member being interrupted, i is the group member doing the interrupting, t is the start time of the interruptive speech, and t' is the end time of the interruptive speech. Needs to be sorted by group first. 
    ? Returns a data frame.
    Refs: Feldstein & W___ (1987), Holme & Saramäki (2012). This algorithm accepts overlapping edges, which Holme & Saramäki do not consider.
    """
    L = []

    data = data[data['pID'].isin(pIDs)]
    for i in range(len(data)):
        ego = data.iloc[i, :]
        overlaps = data[
            ((data['end'] > ego['end']) &
             (data['begin'] < ego['end']) &
             (data['begin'] >= ego['begin']))
        ]

        if overlaps.empty:
            continue
        else:
            for j in range(len(overlaps)):
                alter = overlaps.iloc[j, :]
                # interval is from the start of the interruptive speech until the end of the interruptive speech
                # the start is alter['begin']
                # the end is ego['end']
                start = max([ego['begin'], alter['begin']])
                stop = min([ego['end'], alter['end']])
                l = (ego['pID'], alter['pID'], start, stop, stop - start) 
                L.append(l)
    L = sorted(L, key = lambda x: x[2])
    L = pd.DataFrame(L, columns = ['j', 'i', 'begin', 'end', 'dur'])
    return L

def non_interruptive_simultaneous_speech(data, pIDs):
    """
    Given a pd.DataFrame with variables ['gID', 'pID', 'begin', 'end'], generate an interval sequnce of events (j, i, t, t') where j is the group member being interjected, i is the group member doing the interjective, t is the start time of the non-interruptive speech, and t' is the end time of the non-interruptive speech. Needs to be sorted by group first.
    ? Returns a data frame.
    Refs: Feldstein & W___ (1987), Holme & Saramäki (2012). This algorithm accepts overlapping edges, which Holme & Saramäki do not consider.
    """
    L = []

    data = data[data['pID'].isin(pIDs)]
    for i in range(len(data)):
        ego = data.iloc[i, :]
        overlaps = data[
            # I fail to interrupt someone else
            ((data['begin'] <= ego['begin']) &
             (data['end'] > ego['end']))
        ]

        if overlaps.empty:
            continue
        else:
            for j in range(len(overlaps)):
                alter = overlaps.iloc[j, :]
                # interval is from the start of the interruptive speech until the end of the interruptive speech
                # the start is alter['begin']
                # the end is ego['end']
                start = max([ego['begin'], alter['begin']])
                stop = min([ego['end'], alter['end']])
                l = (ego['pID'], alter['pID'], start, stop, stop - start) 
                L.append(l)
    L = sorted(L, key = lambda x: x[2])
    L = pd.DataFrame(L, columns = ['j', 'i', 'begin', 'end', 'dur'])
    return L

    
def time_to_contact_network(C, pIDs):
    """
    Given a contact sequence, C, return the corresponding time-to-contact network in which nodes are group members, edges are the first time i interrupts j (edge goes from j to i), and edges are weighted by the time elapsing from the beginning of the observation period until the focal interruption.
    """
    time_to_contact_edges = []

    for pID in pIDs:
        contacts = [(j, i, t) for j, i, t in C if i == pID]
        if contacts:
            first_contact = min(contacts, key = lambda x:x[2])
            time_to_contact_edges.append(first_contact)
    cnet = nx.DiGraph()
    cnet.add_weighted_edges_from(time_to_contact_edges)
    return cnet

def interruption_network(C, pIDs):
    """
    Given a contact sequence, C, return the corresponding interruption network in which nodes are group members, edges are interruptions (j, i, c) if i interrupts j c times. 
    """
    interrupts = [(j, i, 0) for j, i, t in C]
    interrupt_set = set(interrupts)
    interrupt_count = [[j, i, w] for j, i, w in interrupt_set]
    for edge in interrupt_count:
        for interrupt in interrupts:
            if (interrupt[0] == edge[0]) & (interrupt[1] == edge[1]):
                edge[2] += 1
    interruption_edges = [(j, i, w) for j, i, w in interrupt_count]
    inet = nx.DiGraph()
    for pID in pIDs:
        inet.add_node(pID)
    inet.add_weighted_edges_from(interruption_edges)
    return inet

def interruption_network_pandas(data, pIDs, use = 'iss'):
    """
    Given a data frame with variables ['gID', 'pID', 'begin', 'end'], generate an interruption network. 
    Default is interruptions only. For a network based on non-interruptive simultaneous speech only, use 'nss'. Use 'both' to include both.
    Not yet implemented: retain the counts of 'iss' and 'nss' in the edge characteristics and track the total count of both in the node characteristics.
    """
    ξ = interruptive_simultaneous_speech(data, pIDs)
    ζ = non_interruptive_simultaneous_speech(data, pIDs)

    #print(len(ξ))
    #print(sum(list(dict(nx.get_edge_attributes(og, 'weight')).values())))

    iss = ξ.groupby(['j', 'i']).agg(
        weight = pd.NamedAgg(column = 'begin', aggfunc = "count"))
    iss.reset_index(inplace = True)

    nss = ζ.groupby(['j', 'i']).agg(
        weight = pd.NamedAgg(column = 'begin', aggfunc = "count"))
    nss.reset_index(inplace = True)

    if use == 'iss':
        g = nx.from_pandas_edgelist(iss, source = 'j', target = 'i', edge_attr = 'weight', create_using = nx.DiGraph)
        disconnected = [pID for pID in pIDs if pID not in list(g.nodes)]
        g.add_nodes_from(disconnected)
        return g
    
    elif use == 'nss':
        g = nx.from_pandas_edgelist(nss, source = 'j', target = 'i', edge_attr = 'weight', create_using = nx.DiGraph)
        disconnected = [pID for pID in pIDs if pID not in list(g.nodes)]
        g.add_nodes_from(disconnected)
        return g
    
    elif use == 'both':
        dat = pd.concat([iss, nss], ignore_index = True)
        subset = ['j', 'i']
        duplicates = dat[dat.duplicated(subset = subset)]
        dat.drop_duplicates(subset = subset, inplace = True)
        for d in duplicates.index:
            for o in dat.index:
                if duplicates.loc[d, 'j'] == dat.loc[o, 'j'] and duplicates.loc[d, 'i'] == dat.loc[o, 'i']:
                    dat.loc[o, 'weight'] += duplicates.loc[d, 'weight']

        g = nx.from_pandas_edgelist(dat, source = 'j', target = 'i', edge_attr = 'weight', create_using = nx.DiGraph)
        disconnected = [pID for pID in pIDs if pID not in list(g.nodes)]
        g.add_nodes_from(disconnected)
        return g
    
    else:
        print('Use one of "iss", "nss", or "both".')

def vote_network(data, pIDs, vote_cols, self_loops = False):
    """
    Given a set of participants, ~pIDs~, and a set of columns containing directed vote data, ~vote_cols~, return a directed "vote network," sometimes called a "leadership network."
    This function expects a pd.DataFrame, ~data~, with a column for participant IDs, a corresponding integer workstation ID ('wID'), and one or more columns of votes, the values of which are also wIDs corresponding to other pIDs.
    """
    dat = data[data['pID'].isin(pIDs)]

    dat = dat.set_index('pID', drop = False)
    
    vnet_edgelist = []
    vnet = nx.DiGraph()
    
    for pID in pIDs:
        vnet.add_node(pID)

    for pID in pIDs:
        for v in vote_cols:
            if pd.notna(dat.loc[pID, v]):
                i = pID
                j = dat[dat['wID'] == dat.loc[pID, v]]['pID'][0]
                if self_loops == False and i == j:
                    continue
                else:
                    new_edge = (str(i), str(j))
                vnet_edgelist.append(new_edge)
    
    vnet.add_edges_from(vnet_edgelist)

    return vnet

def turn_based_network(data, pIDs):
    """
    Given a pd.DataFrame with variables ['gID', 'pID', 'begin', 'end'], generate a "turn-based" (my term for Sauer & Kauffeld's proposed network representation of a conversation) network of the conversation.
    Ref: Sauer & Kauffeld (2013)
    """
    data = data[data['pID'].isin(pIDs)]
    #data.sort_values(by = 'begin', inplace = True, ignore_index = True)
    
    interruptions = interruptive_simultaneous_speech(data, pIDs)
    interjections = non_interruptive_simultaneous_speech(data, pIDs)
    # Sauer & Kauffeld drop simultaneous speech from consideration
    # this is the old code:
    # data = data[~(data['begin'].isin(interruptions['begin'])) &
    #             ~(data['begin'].isin(interjections['begin']))]
    # the below removes interjections, since they don't switch a turn, but keeps interruptions, since they do. How to do deal with ties? I think drop ties, or make the arrow go to the one who spoke longer. 
    data = data.loc[~data["begin"].isin(interjections["begin"]), ]
    # remove speech events that have the same begin and end times but different speakers. This seems to remove the fewest speaking events while still ignoring Dabbs & Ruback's "group turns" as Sauer & Kauffeld seem to have done. By removing interjections and duplicates, data are increasing in both begin and end after the sort_values call, so arrows will always go to the person who finished speaking next.
    data = data.loc[~(data.duplicated(subset = ["begin", "end"], keep = False)), ]
    data.reset_index(inplace = True, drop = True)
    data.sort_values(by = 'begin', inplace = True, ignore_index = True)

    towhom = [np.nan] * len(data)
    for i in data.index[:-1]:
        towhom[i] = data.loc[i + 1, 'pID']
    data['towhom'] = towhom

    wtw = np.full(shape = (len(pIDs), len(pIDs)), fill_value = 0)
    wtw = pd.DataFrame(wtw, columns = pIDs, index = pIDs)

    for i in data.index[:-1]:
        who = data.loc[i, 'pID']
        whom = data.loc[i, 'towhom']
        wtw.loc[who, whom] += 1

    # Successive speaking events from the same speaker are considered part of the same turn
    for i, j in zip(wtw.index, list(wtw)):
        if i == j:
            wtw.loc[i, j] = 0

    g = nx.from_pandas_adjacency(wtw, create_using = nx.DiGraph)

    return g

############################################
#### Random Network Generator Functions ####
############################################

def gnp_from_data(sizes, densities, directed = True):#, p_disconnect_node = None):
    """
    Given a set of graph sizes (number of nodes) and densities, generate a new gnp (Bernoulli/Erdos-Renyi) random graph with size selected from the given graph sizes. Density is estimated based on a linear model of the observed sizes and densities. Currently defaults to directed graph.
    """
    sizes = np.array(sizes)
    densities = np.array(densities)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x = sizes, y = densities)
    dhats = sizes*slope + intercept
    resid = densities - dhats
    resid_std = np.std(resid)

    size = random.choice(sizes)
    dhat = size*slope + intercept
    err = stats.norm(0, resid_std).rvs()
    d = dhat + err

    rg = nx.gnp_random_graph(size, d, directed = directed)

    return rg

def randomize_edges(g, p = 0.5):
    """
    Given graph g, with probability p = 0.5 either [(i, j), (k, l)] --> [(i, l), (k, j)] or [(i, k), (j, l)].
    
    Milo et al. (2002) Algorithm A (SI)
    Kivelä et al (2014) Randomize Edges (proximal source)
    """
    rg = g.copy()

    for z in range(g.number_of_edges()):
        conditions_met = False

        trial = 0
        while conditions_met == False:
            if trial == 1000:

               break
            trial += 1
            edges = list(rg.edges())
            e = edges[z]
            f = random.choice(edges)

            i, j = e[0], e[1]
            k, l = f[0], f[1]

            if random.random() < p:
                new_e = (i, l)
                new_f = (k, j)
            else:
                new_e = (i, k)
                new_f = (j, l)

            if i != l and k != j and i != k and j != l and new_e != new_f and new_e not in rg.edges() and new_f not in rg.edges():
                conditions_met = True

        if trial < 1000:
            rg.remove_edges_from([e, f])
            rg.add_edges_from([new_e, new_f])

        if rg.number_of_edges() < g.number_of_edges():
            print(f'fault at {z}')
            break

    return rg

def randomize_edge_directions(g, p = .5):
    """
    Given g, for each (i, j) in g with probability p reverse the direction of all edges between i and j. Weights remain associated with edges even if edges change direction.
    """
    """
    Hiroki's Algorithm:
    1. Create a list of all node pairs that are connected by at least one way (or both ways).
    
    2. For each pair in the list created above, independently decide whether you want to swap the directions of their edges (with, say, 50% probability).
    
    3. For those node pairs for which the edge direction reversal was decided in 2, swap the directions of the edges between them (i.e., i –> j becomes j –> i, and j –> i becomes i –> j; you can just remove the old edges and create new ones).
    """
    old_nodes = list(g.nodes)
    node_pairs = list(itertools.combinations(g.nodes, 2))
    connected_nodes = {pair: [] for pair in node_pairs}
    for pair in connected_nodes.keys():
        u = pair[0]
        v = pair[1]
        for a, b, c in g.edges(data = True):
            if (a, b) == (u, v) or (a, b) == (v, u):
                connected_nodes[pair].append((a, b, c))
                
    for pair in list(connected_nodes.keys()):
        if connected_nodes[pair] == []:
            del connected_nodes[pair]
    
    new_edges = []
    for pair in list(connected_nodes.keys()):
        if random.random() < p:
            for edge in connected_nodes[pair]:
                new = (edge[1], edge[0], edge[2])
                new_edges.append(new)
        else:
            for edge in connected_nodes[pair]:
                new_edges.append(edge)

    new_graph = nx.DiGraph()
    new_graph.add_nodes_from(old_nodes)
    new_graph.add_edges_from(new_edges)

    return new_graph

def gnp_from_graph(g, seed = None, directed = False):
    """
    Given graph g, generate a new gnp random graph using n and p from g.
    Defaults match nx.gnp_random_graph().
    """
    n = g.number_of_nodes()
    p = nx.density(g)
    rg = nx.gnp_random_graph(n = n, p = p, seed = seed, directed = directed)

    return rg

## Need a randomly permuted times (RP) generator

##################################
#### Motif Analysis Functions ####
##################################

def count_motifs(g, motif_size, only_connected = True): # only_connected not implemented?
    freqs = {}
    for bunch in itertools.combinations(g.nodes, motif_size):
        sg = g.subgraph(bunch)
        connected = nx.is_weakly_connected(sg) if sg.is_directed() == True else nx.is_connected(sg)

        if only_connected == True and connected == False:
            continue
        else:
            for motif in freqs:
                if nx.is_isomorphic(sg, motif):
                    freqs[motif] += 1
                    break
            else:
                freqs[sg] = 1
    return freqs

def collect_motif_counts(count_list):
    freqs = {}
    for d in count_list:
        for key in d.keys():
            for motif in freqs:
                if nx.is_isomorphic(key, motif):
                    freqs[motif] += 1
                    break
            else:
                freqs[key] = 1
    return freqs

##############################
#### Burstiness Functions ####
##############################

def bursty_coef(data, finite = True):
    """
    Return the burstiness parameter, B, from a series of inter-event times, τ (here, `data`).
    Ref is Karsai et al 2018, "Bursty Human Dynamics."
    This function will expect a time series of data consisting of inter-event times.
    τ: the time series of inter-event times
    σ: the standard deviation of τ
    τ_bar: 〈τ〉, the mean of τ
    r: σ/τ_bar, the coefficient of variation
    B: burstiness parameter, ≔ (r - 1)/(r + 1) = (σ - τ_bar)/(σ + τ_bar)
    Alt:
    n = number of events in τ
    B_n: burstiness parameter with a correction for finite sample size
         ((√(n + 1)*r) - √(n - 1))/(((√(n + 1) - 2)*r)   + √(n - 1))
    """
    r = np.std(data)/np.mean(data)
    n = len(data)

    if n < 3:
        B = np.nan
    elif finite:
        numerator = (np.sqrt(n + 1)*r) - np.sqrt(n - 1)
        denominator = ((np.sqrt(n + 1) - 2)*r) + np.sqrt(n - 1)
        B = numerator/denominator
        #B = ((np.sqrt(n + 1)*r) - np.sqrt(n - 1))/(((np.sqrt(n + 1) - 2)*r) + np.sqrt(n - 1))
    else:
        B = (r - 1)/(r + 1)

    return B

def memory_coef(data, m = 1):
    """
    Return the memory coefficient, M, from a series of inter-event times. Depends on m, the the "lag" in the memory, representing m-1 intermediate events.
    Ref is Karsai et al 2018, "Bursty Human Dynamics."
    M: the memory coefficient
    m: the "lag"
    τ, τ_bar, σ, and n as for bursty_coef().
    There is both a τ_1 and a τ_2, and corresponding values, separated by m
    """
    data = list(data)
    n = len(data)

    if n < m + 3:
        M = np.nan
        return M
    
    unlagged = data[:n-m-1]
    tau1bar = np.mean(unlagged)
    sigma1 = np.std(unlagged)
    lagged = data[m:-1]
    tau2bar = np.mean(lagged)
    sigma2 = np.std(lagged)

    try:
        norm = 1/(n - m - 1)
    except:
        M = np.nan
        return M
    
    summation = []
    for i in range(n - m - 1): # n - m - 1
        numerator = (unlagged[i] - tau1bar)*(lagged[i] - tau2bar)
        denominator = sigma1*sigma2
        summation.append(numerator/denominator)
        
    M = norm*sum(summation)
    
    return M

##########################
## ABM Helper Functions ##
##########################

def Y_to_X(Y, ns):#(Y, R, ns):
    """
    Take the simulation output, Y, and convert it into a pd.DataFrame, X, for further analysis.
    """
    X = []
    for n in ns:
        # gID = R
        # pID = n # pID = R + str(n)
        curr = 0 # current state
        last = 0 #last state
        dur = 0
        begin = 0

        #for t in range(len(Y[n])):
        for t in range(len(Y)):
            #curr = Y[n][t] # could this just change to Y[t, n]?
            curr = Y[t, n]
            if curr == 1:
                if last == 1:
                    dur += 1
                else:#elif: last == 0:
                    begin = t
            else:#elif curr == 0:
                if last == 1:
                    dur += 1
                    #row = (R, n, begin, dur) # why am I doing this? R must be a group marker
                    # better to keep the group and X as separate constructs
                    # can use a dict set up instead
                    row = (n, begin, dur)
                    X.append(row)
                    begin = 0
                    dur = 0
                else:#elif last == 0:
                    continue
            last = curr
    X = pd.DataFrame(X, columns = ['pID', 'begin', 'dur'])#['gID', 'pID', 'begin', 'dur'])

    X['end'] = X['begin'] + X['dur']
    X['lat'] = np.nan

    #_gID = list(X).index('gID')
    _pID = list(X).index('pID')
    _begin = list(X).index('begin')
    _end = list(X).index('end')
    _dur = list(X).index('dur')
    _lat = list(X).index('lat')

    for i in range(len(X.index)):
        if X.iloc[i, _pID] != X.iloc[i - 1, _pID] or i == 0:
            X.iloc[i, _lat] = X.iloc[i, _begin]
        else:
            X.iloc[i, _lat] = X.iloc[i, _begin] - X.iloc[i - 1, _begin]

    return X

#############################
## Visualization Functions ##
#############################

solarized = {# https://ethanschoonover.com/solarized/
    'base03': '#002b36', 'base02': '#073642', 'base01': '#586e75', 'base00': '#657b83',
    'base0': '#839496', 'base1': '#93a1a1', 'base2': '#eee8d5', 'base3': '#fdf6e3',
    'yellow': '#b58900', 'orange': '#cb4b16', 'red': '#dc322f', 'magenta': '#d33682',
    'violet': '#6c71c4', 'blue': '#268bd2', 'cyan': '#2aa198', 'green': '#859900'}
whiteboard = {# from the Emacs whiteboard theme
    "Black": "#000000", "Gray40": "#666666", "Gray50": "#7F7F7F", "Gray60": "#999999", "Gray75": "#BFBFBF",
    "Gainsboro": "#DCDCDC", "Whitesmoke": "#F5F5F5", "Green4": "#008B00", "DarkOliveGreen4": "#6E8B3D",
    "Burlywood4": "#8B7355", "DarkOrange3": "#CD6600", "Gold3": "#CDAD00", "Blue4": "#00008B",
    "SkyBlue1": "#00BFFF", "Red": "#FF0000"}   

def visualize_speaking_data(data, gID = None, pID = None, cond = None, rep = None, ax = None, colors = 'k'):
    """
    Generate a plot that uses horizontal bars to represent periods of time when group members (or agents) are speaking. Plotting set up like plt.subplots() and plt.show() should be done outside of the function.
    """
    if gID:
        figdat = data.loc[data["gID"] == gID, ]
        pIDs = pd.unique(figdat['pID'])
        if cond:
            figdat = data[(data['type'] == cond) &
                          (data['gID'] == gID)]
        if rep:
            figdat = figdat[figdat['rep'] == rep]
    else:
        figdat = data.loc[data["pID"] == pID, ]
        pIDs = [pID]
        

    b, c = .25, 0
    yticks = []
    for pID in pIDs:
       dat = figdat[figdat['pID'] == pID].sort_values(['begin'])
       lines = [[(i, b), (j, b)] for i, j in zip(dat['begin'], dat['end'])]
       lc = LineCollection(lines, linewidths = 10, colors = colors, label = pID) # colors[c]
       ax.add_collection(lc)
       
       yticks.append(b)
       b += 0.25
       if dat.empty:
           continue
       c += 1

    ax.set_xlim([0, max(data['end']) + 5])
    ax.set_ylim([0, b])
    ax.set_yticks(yticks)
    ax.set_yticklabels(pIDs)
    ax.set_xlabel('t (seconds)')
    #ax.set_ylabel('Group Member')
    if cond:
        ax.set_title(cond.capitalize())
