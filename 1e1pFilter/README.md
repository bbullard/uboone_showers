# Low-energy excess analyzer
This module is a [LArSoft](http://www.larsoft.org) analyzer that builds a ROOT `TTree` with information from &nu;<sub>e</sub> candidates, reconstructed with the [Pandora framework](https://github.com/PandoraPFA).
It looks for neutrino PFParticles with at least one daughter shower and at least one daughter track.

There are two FCL files, one for data (`run_PandoraOnly_data.fcl`) and one for Monte Carlo (`run_PandoraOnly.fcl`).

## Requirements

- `uboonecode v06_26_01_10` (MCC8.6)
- `uboonecode v06_26_01_08` (MCC8.4)
- [XGBoost](http://xgboost.readthedocs.io/en/latest/). Follow the instructions in [DocDB 10685](https://microboone-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=10685).
