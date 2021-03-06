# interruption-networks
Data handling, analysis, and simulation files for working with interruption networks.

---

The working data for this study can be downloaded from Binghamton University's [Open Repository @ Binghamton](https://orb.binghamton.edu/management_fac/2/). To work with these files, download the repository and place `InterruptionAnalysis.py` and any script files in the home directory. Paths in the script files in this repository assume a path `<home>/data/diarizations/`. `Python-to-R.py` depends on the file `timeseries.csv` produced by `collect-ts-data.py`.

## File Notes

The files as written expect the following files and directories:

- `./data/diarizations/` with `{gID}.eaf` annotation files, where `gID` is the group ID. 
- `./data/vote-data.csv` containing leader emergence votes (who voted for whom).
- `./data/speakingTime-data.csv` containing the results of demographic, psychometric, and sociometric surveys. 

These files are available from Binghamton University's [Open Repository @ Binghamton](https://orb.binghamton.edu/management_fac/2/).

With the above files and directory structure, run `collect-ts-data.py` to produce `./data/timeseries.csv`, which is the base data structure for building and analyzing interruption networks. A row contains the group ID, individual id (`pID`), and the start and stop time of the speaking event, its duration, and "latency" (how long since the last speaking event for that `pID`).

From this point, you can run `generate-networks.py` to make a set of interruption networks, turn-based networks (see Sauer & Kauffeld, 2013), and vote networks (also called "leadership networks", see Carter et al., 2015). Analysis files depend on these networks.

<!-- - `InterruptionAnalysis.py` contains the functions needed for the data collection and analysis steps used in this project. -->
<!-- - `collect-ts-data.py` transforms the diarization files into a `pandas` DataFrame that contains the group ID ('gID') and participant ID ('pID'); beginning ('begin') and ending ('end') time, in milliseconds, for each vocalization; and two calculated columns for the duration of each vocaliation ('dur') and the "latency" ('lat') for each vocalization for inter-event time analysis. -->
<!-- - `interruption-network-example.py` produces example interruption and vote networks for analysis. -->

## Software

The original diarizations were done using [ELAN](https://archive.mpi.nl/tla/elan).

Analysis and simulations are done in Python 3, relying on the following packages: `numpy`, `pandas`, `matplotlib`, `statsmodels`, and `networkx`. ERGM analysis uses R 4.0.3 and the `statnet` and `tidyverse` package families, as well as `igraph` and `intergraph`. 
