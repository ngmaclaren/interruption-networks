import numpy as np
import mpmath
import pandas as pd
import matplotlib.pyplot as plt
import random
import powerlaw

def step1(data):
    fit = powerlaw.Fit(dat, verbose = False)
    n = len(dat)
    xmin = fit.xmin
    alpha = fit.alpha
    tail = dat[dat >= xmin]
    ntail = len(tail)
    not_tail = dat[dat < xmin]

    return (fit, n, xmin, alpha, tail, ntail, not_tail)

def step2(fit, n, xmin, alpha, tail, ntail, not_tail, n_synth = 2500):
    emp_ks = fit.power_law.D
    
    ks_dist = []
    for _ in range(n_synth):
        theoretical = powerlaw.Power_Law(xmin = xmin, parameters = [alpha])
        theosim = theoretical.generate_random(1000)
        
        synth = []
        for _ in range(n):
            prob = ntail/n
            if random.random() <= prob:
                synth.append(random.choice(theosim))
            else:
                synth.append(random.choice(list(not_tail)))

        synth_fit = powerlaw.Fit(synth, verbose = False)
        synth_ks = synth_fit.power_law.D
        ks_dist.append(synth_ks)

    p = len([ks for ks in ks_dist if ks > emp_ks])/len(ks_dist)
    
    return p

def step3(data, fit, ID):
    comp_dists = ['lognormal', 'lognormal_positive', 'exponential', 'stretched_exponential',
              'truncated_power_law']
    comp = {}

    dists = {}
    for dist in comp_dists:
        R, p = fit.distribution_compare('power_law', dist)
        dists[f'{dist}_R'] = R
        dists[f'{dist}_p'] = p
    comp[ID] = dists

    comp = pd.DataFrame(comp).T
    comp['n'] = len(data)

    return comp


data = pd.read_csv('./data/timeseries.csv', index_col = 0)
group_kld = pd.read_csv('./data/group-kld.csv', index_col = 0)
indiv_kld = pd.read_csv('./data/indiv-kld.csv', index_col = 0)

gIDs = pd.unique(data['gID'])
test = random.sample(list(gIDs), 2)

comp_results = []
for gID in gIDs:
    print(gID)
    dat = data[data['gID'] == gID]['lat']
    fit, n, xmin, alpha, tail, ntail, not_tail = step1(dat)

    fit_p = step2(fit, n, xmin, alpha, tail, ntail, not_tail)

    comp = step3(dat, fit, gID)
    comp['a'] = alpha
    comp['xmin'] = xmin
    comp['fit_p'] = fit_p
    comp_results.append(comp)#[gID] = comp

comps = pd.concat(comp_results)

pIDs = pd.unique(data['pID'])

ind_comp_results = []
for pID in pIDs:
    print(pID)
    dat = data[data['pID'] == pID]['lat']
    if len(dat) < 10:
        print(f'{pID} has fewer than 10 data points. Continuing.')
        continue
    
    fit, n, xmin, alpha, tail, ntail, not_tail = step1(dat)

    fit_p = step2(fit, n, xmin, alpha, tail, ntail, not_tail)

    comp = step3(dat, fit, pID)
    comp['a'] = alpha
    comp['xmin'] = xmin
    comp['fit_p'] = fit_p
    ind_comp_results.append(comp)#[gID] = comp

ind_comps = pd.concat(ind_comp_results)


comps.to_csv('./data/group-powerlaw-results.csv')
ind_comps.to_csv('./data/indiv-powerlaw-results.csv')
print('done')

