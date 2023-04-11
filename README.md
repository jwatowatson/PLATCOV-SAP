# PLATCOV-SAP

Statistical analysis plan for the PLATCOV trial.
The PLATCOV trial is registered at clinicaltrials.gov number [NCT05041907](https://clinicaltrials.gov/ct2/show/NCT05041907)

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Publications and outputs

This git repo will change over time. For reproducibility of the individual results, each trial output has its own git repo.

Ivermectin data and results: github repo [PLATCOV-Ivermectin](https://github.com/jwatowatson/PLATCOV-Ivermectin)
 and [preprint](https://www.medrxiv.org/content/10.1101/2022.07.15.22277570v1).
 
Casirivimab/imdevimab and remdesivir (the parenteral treatment arms): [github repo](https://github.com/jwatowatson/PLATCOV-Remdesivir-Regeneron) and preprint not yet online

## Overview

This github repo provides the study protocol (master protocol) and the statistical analysis plan (both in folder Analysis_Plan and PLATCOV_SAP_*.pdf), and the generic code used for the statistical analysis of each study arm (Generic_intervention_analysis.qmd)

Each interim analysis is done by running the full workflow given in *Generic_intervention_analysis.qmd* applied to the data from each arm separately. Each analysis is encoded as a separate csv file with the intervention and the concurrent controls (either negative controls for the futility/success analyses, or both negative and positive controls for the non-inferiority analyses). 


The generic statistical analysis for each arm does the following:

* Loads data: the dataset is specified via the variable *intervention*. The reference arm also needs to be specified (eg *no study drug*);
* Checks patients in dataset against ITT data and produces a mITT variable (does the patient have at least 3 days where samples are per protocol?);
* Makes some summary data plots and tables;
* Sets up the model runs for stan along with run parameters (number of chains etc..);
* Model fitting is done separately via the code in *run_models_local.R* (to do on local machine) or *run_models.R* (to do on cluster). The stan models are provided in the folder *Stan_models*;
* Checks for convergence and compares model fits
* Displays results
* Some sensitivity analyses

The underlying data are not made publicly available until publication of results.

## Model structure

The stan models all treat the PCR data as left-censored. All viral load data are on the log base 10 viral copies per mL scale. We fit regression models with varying degrees of complexity:

* Base model M0: individual/site random effects for slope and intercept
* M1: add human RNaseP correction (more human cells taken up by swab should in theory indicate more virus)
* M2: covariate effects (age, number of vaccine doses, serology, time since symtom onset)
* M3: non-linear version of M1, whereby the virus can still be in a growth phase at the start of the follow-up and then decreases.


## Software needed

The R packages needed are:

* *rstan* (interfaces with *stan*)
* *loo* (for model comparison)
* *censReg* (censored regression - sensitivity analysis)
* *RColorBrewer* and *tictoc* (plotting and timing)


Any questions or comments drop me a message at jwatowatson at gmail dot com


