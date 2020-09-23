# interruption-networks
Data handling, analysis, and simulation files for working with interruption networks. A work in progress.

---

The working data for this study can be downloaded from Binghamton University's [Open Repository @ Binghamton](https://orb.binghamton.edu/management_fac/2/). To work with these files, download the repository and place `InterruptionAnalysis.py` and any script files in the home directory. Paths in the script files in this repository assume a path `<home>/data/diarizations/`. `Python-to-R.py` depends on the file `timeseries.csv` produced by `collect-ts-data.py`.

## File Notes

- `InterruptionAnalysis.py` contains the functions needed for the data collection and analysis steps used in this project.
- `collect-ts-data.py` transforms the diarization files into a `pandas` DataFrame that contains the group ID ('gID') and participant ID ('pID'); beginning ('begin') and ending ('end') time, in milliseconds, for each vocalization; and two calculated columns for the duration of each vocaliation ('dur') and the "latency" ('lat') for each vocalization for inter-event time analysis.
- `interruption-network-example.py` produces example interruption and vote networks for analysis.

## Software

The original diarizations were done using [ELAN](https://archive.mpi.nl/tla/elan).

Analysis and simulations are done in Python 3, relying on the following packages: `numpy`, `pandas`, `matplotlib`, `statsmodels`, and `networkx`.