import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import InterruptionAnalysis as ia

data = pd.read_csv('./data/timeseries.csv', index_col = 0)
votedata = pd.read_csv('./data/vote-data.csv')
a = 0.99

gID = 'YMY'
pIDs = pd.unique(votedata[votedata['gID'] == gID]['pID'])

votecols = ['voteFor_p1', 'voteFor_p2', 'voteFor_p3', 'voteFor_p4', 'voteFor_p5']
numeric_cols = ['begin', 'end', 'dur', 'lat']
data[numeric_cols] = data[numeric_cols].apply(lambda x: x/1000)

mc = "#00688b"
fc = "#cdad00"
colors = [mc, mc, mc, mc, fc, mc, fc, fc, fc, mc]

plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(figsize = [12, 4])

b, c = .25, 0
yticks = []
for pID, color in zip(pIDs, colors):
    dat = data[data['pID'] == pID]
    dat = dat.sort_values(['begin'])
    lines = [[(i, b), (j, b)] for i, j in zip(dat.begin, dat.end)]
    lc = LineCollection(lines, linewidths = 10, colors = color, label = pID)
    ax.add_collection(lc)

    yticks.append(b)
    b += .25
    if dat.empty:
        continue
    c += 1
ax.set_xlim([0, max(data['end']) + 5])
ax.set_ylim([0, b])
ax.set_yticks(yticks)
ax.set_yticklabels([pid[3:] for pid in pIDs])
ax.set_xlabel('t (sec)')
ax.set_ylabel('Participant')

fig.tight_layout()
#fig.savefig('./img/ymy-speechactivity.svg', transparent = True)
fig.show()
