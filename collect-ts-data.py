import glob
import numpy as np
import pandas as pd

import InterruptionAnalysis as ia

loc = './data/diarizations'
emp = 'orig'
ext = 'pl.eaf'

datafiles = glob.glob(f'{loc}/*{ext}')

data = [ia.collect_vocalization_sequence(df) for df in datafiles]
data = pd.concat(data)

data['lat'] = np.nan

_gID = list(data).index('gID')
_pID = list(data).index('pID')
_begin = list(data).index('begin')
_end = list(data).index('end')
_dur = list(data).index('dur')
_lat = list(data).index('lat')

for i in range(len(data.index)):
    if data.iloc[i, _pID] != data.iloc[i - 1, _pID]:
        data.iloc[i, _lat] = data.iloc[i, _begin]
    else:
        data.iloc[i, _lat] = data.iloc[i, _begin] - data.iloc[i - 1, _begin]

data.to_csv('./data/timeseries.csv')

print('done')
