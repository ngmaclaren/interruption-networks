import random
import itertools
import collections
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import networkx as nx
import statsmodels.api as sm
from statsmodels.formula.api import ols

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
    From Lacasa et al 2012. Data needs to have data['begin'] and data['dur'].
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
    Given a directed horizontal visibility graph (doesn't currently test for that...), return the Kullback-Leibler divergence between the out-degree probability distribution (forward in time) and the in-degree probability distribution (backward in time)
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

###################################
#### Markov Analysis Functions ####
###################################

def get_transition_matrix(data, x): #where ~x~ is a pID
    """
    For a given data frame ~data~, return the Markov transition probability matrix associated with group member ~x~. ~data~ should have the columns 'pID', 'gID', 'start', 'end', and 'dur'. 
    Let u_X be the number of distinct speaking events coded for participant X.
    Let tst be the sum of 
    Let A represent the event that participant X is talking and B represent the event that participant X is not talking. 
    We have a matrix P_X for each participant the values of which are represented by 
from\to  A      B
A        P(A|A) P(B|A)
B        P(A|B) P(B|B)
    """
    gID = list(data[data['pID'] == x]['gID'])[0]
    grpdata = data[data['gID'] == gID]
    X = data[data['pID'] == x]

    u_X = len(X)
    N_A = sum(X['dur'])
    N_B = max(grpdata['end']) - N_A
    P_BA = u_X/N_A
    P_AB = u_X/N_B
    P_AA = 1 - P_BA
    P_BB = 1 - P_AB

    P = np.array([[P_AA, P_BA],
                  [P_AB, P_BB]])
    return P


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

    return ICC1

#####################################
#### Network Generator Functions ####
#####################################

def contact_sequence(data, pIDs):
    """
    Given a pd.DataFrame with variables ['gID', 'pID', 'begin', 'end'], generate a contact sequence of events (j, i, t) where j is the group member being interrupted, i is the group member doing the interrupting, and t is the start time of the interruptive speech.
    """
    C = []

    data = data[data['pID'].isin(pIDs)]
    for i in range(len(data)):
        ego = data.iloc[i, :]
        interruptions = data[(data['end'] > ego['end']) &
                              (data['begin'] < ego['end'])]
        if interruptions.empty:
            continue
        else:
            for j in range(len(interruptions)):
                alter = interruptions.iloc[j, :]
                c = (ego['pID'], alter['pID'], alter['begin'])
                C.append(c)
    C = sorted(C, key = lambda x: x[2])
    return C

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
        
def vote_network(data, pIDs, vote_cols):
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
                new_edge = (str(i), str(j))
                vnet_edgelist.append(new_edge)
    
    vnet.add_edges_from(vnet_edgelist)

    return vnet    
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
        # connected = nx.is_weakly_connected(sg) if sg.is_directed() == True else nx.is_connected(sg)

        # if only_connected == True and connected == False:
        #     continue
        # else:
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
