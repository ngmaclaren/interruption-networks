# interruption-networks
Data handling, analysis, and simulation files for working with interruption networks. 

---

This code implements interruption networks as described in [Maclaren et al. (2021)](https://psyarxiv.com/m8y5n), turn-based networks as described in [Sauer & Kauffeld (2013)](https://www.tandfonline.com/doi/abs/10.1080/19312458.2012.760729), and networks of emergent leader votes (also called leadership networks, see [Carter et al. (2015)](https://psycnet.apa.org/record/2015-12426-001)). This code was developed  to work with a particular data set which can be downloaded from Binghamton University's [Open Repository @ Binghamton](https://orb.binghamton.edu/management_fac/2/). 

## Use

These files expect to be together in a working directory with a `./data/` subdirectory containing the files and directories from the ORB [repository](https://orb.binghamton.edu/management_fac/2/). `InterruptionAnalysis.py` is a module that contains functions used throughout this project (e.g., `import InterruptionAnalysis as ia`). 

With the above files and directory structure, run `collect-ts-data.py` to produce `./data/timeseries.csv`, which is the base data structure for building and analyzing interruption networks. A row contains the group ID, individual id (`pID`), and the start and stop time of the speaking event, its duration, and "latency" (how long since the last speaking event for that `pID`).

Once you have the diarizations stored in `timeseries.csv`, you can run `generate-networks.py` to create each type of network referred to in [Maclaren et al. (2021)](https://psyarxiv.com/m8y5n) for each group whose data is in the ORB repository. Other files conduct specific analyses or create specific figures (or parts ot them). 

## Known Issues

- The simulations in `edge-direction-sim.py` are very slow, and load even slower in `edge-direction-plot.py`. Any suggestions for speeding these steps up would be greatly appreciated. 

## Software

The original diarizations were done using [ELAN](https://archive.mpi.nl/tla/elan).

Analysis and simulations are done in Python 3, relying on the following packages: `numpy`, `pandas`, `matplotlib`, `statsmodels`, and `networkx`. ERGM analysis uses R 4.0.3 and the `statnet` package families, as well as `igraph` and `intergraph`. Stata analysis conducted in Stata 17. 
