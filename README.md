# PLATCOV-SAP

Statistical analysis plan for the PLATCOV trial.
The PLATCOV trial is registered at clinicaltrials.gov number NCT05041907


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



## Overview

This github repo provides the statistical analysis plan (see PLATCOV_SAP.pdf) and the code used for the statistical analysis.

Each interim analysis is done by running the workflow given in *Full_Analysis.Rmd*. This does the following:

* Loads data and makes some summary data plots
* Does some QC analysis of the PCR data (compares standard curves and check all plates have approximately similar results for the control samples)
* Runs a series of Bayesian hierarchical linear regression models on the available data (these models are coded in *stan* provided in the folder *Stan_models*)
* Checks for convergence and compares model fits
* Displays results
* Some senstivity analysis

The underlying data are not made publicly available until publication of results.



## Software needed

The R packages needed are:

* *rstan* (interfaces with *stan*)
* *loo* (for model comparison)
* *censReg* (censored regression - sensitivity analysis)
* *RColorBrewer* and *tictoc* (plotting and timing)


Any questions or comments drop me a message at jwatowatson at gmail dot com


