---
title: "PLATCOV Statistical Analysis Plan"
author: "James Watson"
date: "18 February, 2022"
output: 
  html_document: 
    toc: yes
    fig_caption: yes
    keep_md: yes
---



TODO list:

* Add QC analysis (comparison of standard curves)
* Up and down model to account for the fact that subjects are in an increasing phase
* Data need to be also on the log viral copies per ml
* How much variance explained by the RNaseP
* Check missing data (vaccine/serology/age)
* Subgroup analysis model for temporal shift
* Define Epochs for the trial to get drift effects
* Model comparison using loo package

## Preambule

This Statistical Analysis Plan (SAP) is written for the PLATCOV trial (registered at ClinicalTrials.gov: https://clinicaltrials.gov/ct2/show/NCT05041907).

Data preparation is done in a different R script called data_prep.R.
This Markdown script assumes that the data are saved in two .csv file in long format. 

The first file interim_control_dat.csv contains the PCR quality control data: CT values from the control samples from each plate (duplicate control samples of known viral densities over the range 100 to 10^7/mL viral copies per mL). This has the following headers:

- ID: unique control ID
- Plate: unique Plate ID for the PCR assay
- CT_NS: observed Ct value for the N/S gene
- log10_true_density: number of viral copies in control sample (log10 scale)

The second file interim_dat.csv contains the patient viral load data with the following column headers:

- ID: anonymised patient id code
- BARCODE: unique sample id code (in the following format: PLT-site-number-swab type-timepoint)
- Swab_ID: RTS or TSL (right versus left tonsil)
- Plate: unique Plate ID for the PCR assay (matching with plate identifiers in interim_control_dat.csv)
- Site: site at enrollment
- Trt: treatment allocation as written in CRF
- Time: time from randomisation
- Variant: variant of concern (using standard WHO terminology for the main lineages, reference will be the predominant variant in the dataset at the start of the study)
- Vaccinated: 0/1 (1: any number of doses and any manufacturer)
- Vaccine_type: manufacturer of vaccine dose 1 (AZ: astra zeneca; SPH: sinopharm; SV: sinovac; MD: moderna; PZ: pfizer)
- Antibody_test: 0/1 (positive/negative for SARS-CoV-2 antibody rapid test)
- Age: (years - has to be between 18-50)
- Sex: 0/1 (male: 1; female/other: 0)
- Symptom_onset: time since onset of symptoms (days)
- CT_NS: observed Ct value for the N/S gene
- CT_RNaseP: observed Ct value for the human RNase P gene
- log10_viral_load: log10 number of viral copies per mL (estimated from plat specific standard curve)


## Setup


```
##                _                           
## platform       x86_64-apple-darwin17.0     
## arch           x86_64                      
## os             darwin17.0                  
## system         x86_64, darwin17.0          
## status                                     
## major          4                           
## minor          0.2                         
## year           2020                        
## month          06                          
## day            22                          
## svn rev        78730                       
## language       R                           
## version.string R version 4.0.2 (2020-06-22)
## nickname       Taking Off Again
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS  10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] censReg_0.5-32       maxLik_1.5-2         miscTools_0.6-26    
## [4] RColorBrewer_1.1-2   rstan_2.21.2         ggplot2_3.3.5       
## [7] StanHeaders_2.21.0-7
## 
## loaded via a namespace (and not attached):
##  [1] bdsmatrix_1.3-4    Rcpp_1.0.7         lattice_0.20-44    prettyunits_1.1.1 
##  [5] ps_1.6.0           zoo_1.8-9          lmtest_0.9-38      assertthat_0.2.1  
##  [9] glmmML_1.1.1       digest_0.6.27      utf8_1.2.2         V8_3.4.2          
## [13] R6_2.5.1           stats4_4.0.2       evaluate_0.14      pillar_1.6.2      
## [17] Rdpack_2.1.3       rlang_0.4.11       curl_4.3.2         rstudioapi_0.13   
## [21] callr_3.7.0        jquerylib_0.1.4    collapse_1.7.6     rmarkdown_2.11    
## [25] stringr_1.4.0      loo_2.4.1          munsell_0.5.0      compiler_4.0.2    
## [29] xfun_0.26          pkgconfig_2.0.3    pkgbuild_1.2.0     htmltools_0.5.2   
## [33] tidyselect_1.1.1   tibble_3.1.4       gridExtra_2.3      codetools_0.2-18  
## [37] matrixStats_0.61.0 fansi_0.5.0        crayon_1.4.1       dplyr_1.0.7       
## [41] withr_2.4.2        rbibutils_2.2.7    MASS_7.3-54        grid_4.0.2        
## [45] nlme_3.1-153       jsonlite_1.7.2     gtable_0.3.0       lifecycle_1.0.0   
## [49] DBI_1.1.1          magrittr_2.0.1     scales_1.1.1       RcppParallel_5.1.4
## [53] cli_3.0.1          stringi_1.7.4      plm_2.6-0          bslib_0.3.0       
## [57] ellipsis_0.3.2     generics_0.1.0     vctrs_0.3.8        sandwich_3.0-1    
## [61] Formula_1.2-4      tools_4.0.2        glue_1.4.2         purrr_0.3.4       
## [65] processx_3.5.2     parallel_4.0.2     fastmap_1.1.0      yaml_2.2.1        
## [69] inline_0.3.19      colorspace_2.0-2   knitr_1.34         sass_0.4.0
```


## Load data and provide data summaries


```
## We have 720 PCR datapoints on 36 patients from 1 sites
```

```
##                Site
## Arm             th001
##   Favipiravir       8
##   Ivermectin        7
##   No study drug     7
##   Regeneron         8
##   Remdesivir        6
```

```
## Age distribution of patients:
```

```
##    Site Age.mean Age.min Age.max
## 1 th001       30      21      49
```

```
## Proportion male (%):
```

```
##    Site Sex
## 1 th001 100
```

```
## [1] "!! WARNING - sex data made up as not in current clinical database that was sent !!"
```

```
## Proportion vaccinated (%):
```

```
##    Site Vaccinated
## 1 th001         81
```

```
## Vaccine type by site:
```

```
##        
## Site    AZ None SPH SV
##   th001 14    7   9  6
```

```
## Number of days since symptom onset by site:
```

```
##        Sympom_Onset_days
## Site     1  2  3  4
##   th001 14 11 10  1
```

```
## Number of baseline viral load samples:
```

```
## 
## PLT-TH1-001 PLT-TH1-002 PLT-TH1-003 PLT-TH1-004 PLT-TH1-005 PLT-TH1-006 
##           4           4           4           4           4           4 
## PLT-TH1-007 PLT-TH1-008 PLT-TH1-009 PLT-TH1-010 PLT-TH1-011 PLT-TH1-012 
##           4           4           4           4           4           4 
## PLT-TH1-013 PLT-TH1-014 PLT-TH1-015 PLT-TH1-016 PLT-TH1-017 PLT-TH1-018 
##           4           4           4           4           4           4 
## PLT-TH1-019 PLT-TH1-020 PLT-TH1-021 PLT-TH1-022 PLT-TH1-023 PLT-TH1-024 
##           4           4           4           4           4           4 
## PLT-TH1-025 PLT-TH1-026 PLT-TH1-027 PLT-TH1-028 PLT-TH1-029 PLT-TH1-030 
##           4           4           4           4           4           4 
## PLT-TH1-031 PLT-TH1-032 PLT-TH1-033 PLT-TH1-034 PLT-TH1-042 PLT-TH1-045 
##           4           4           4           4           4           4
```

```
## [1] TRUE
```

```
## In the 36 patients, the geometric mean baseline (defined as samples taken within 6 hours of randomisation) viral load was 291643 copies per mL (95 percentile interval from 1300 to 65433435; range from 370 to 76230035)
```

![](SAP_files/figure-html/load_data-1.png)<!-- -->

## Data visualisation

### Quality Control for PCR


```
## SD of the intercept is 0.52
```

```
## Maximum difference in CT values for a tenfold change in viral load is 0.2
```

![](SAP_files/figure-html/PCR_QC-1.png)<!-- -->![](SAP_files/figure-html/PCR_QC-2.png)<!-- -->



### Overall data

![](SAP_files/figure-html/trt_data_plot-1.png)<!-- -->


### RNaseP

![](SAP_files/figure-html/RNaseP-1.png)<!-- -->

```
## the standard deviation of the RNaseP CT distribution is 1.26.
```

```
## 95% of the distribution is within 4.92 CT values difference
```


## Model fitting


Main models: Bayesian hierarchical models with and without covariate adjustment. All models adjust for human RNaseP quantification (scaled to have mean 0).

Covariates that we use in model 2:

* Vaccination (yes/no)
* Time since symptom onset (days)
* Variant (made up data)
* Serology rapid test (+/-)



### Specify priors


```
## $A0_prior_mean
## [1] 18
## 
## $A0_prior_sd
## [1] 5
## 
## $alpha_prior_mean
## [1] -2
## 
## $alpha_prior_sd
## [1] 2
## 
## $sigma_trt_effect
## [1] 0.5
```


### Prepare model

compile model in stan




make stan data set


```
## Total number of datapoints up until day 8 is 648
```

### Run models

run models


```
## Warning: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

```
## Warning: Examine the pairs() plot to diagnose sampling problems
```

```
## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## Running the chains for more iterations may help. See
## http://mc-stan.org/misc/warnings.html#bulk-ess
```

```
## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
## Running the chains for more iterations may help. See
## http://mc-stan.org/misc/warnings.html#tail-ess
```

### Model fits: summaries


```
##                      mean     se_mean         sd
## A0            15.59702468 0.022097510 0.60486781
## alpha         -1.55513326 0.008534694 0.23770647
## trt_effect[1] -0.12349436 0.008216898 0.21680978
## trt_effect[2] -0.18421896 0.008128133 0.23116343
## trt_effect[3]  0.03530532 0.008262219 0.22618872
## trt_effect[4] -0.28992107 0.008306994 0.22280697
## sigmaCT        2.76514261 0.005257483 0.15625396
## t_dof          4.24553065 0.030799762 0.88629362
## gamma_rnasep   0.61207607 0.005090706 0.14137562
## sigma_plate    0.46948595 0.004398391 0.12834999
## sigmasq_u[1]   3.10789550 0.014524483 0.41953654
## sigmasq_u[2]   0.46798625 0.002985289 0.08366125
```

```
##                      mean     se_mean         sd
## A0            15.48982757 0.020192251 0.55578496
## alpha         -1.48491231 0.008606630 0.24665477
## trt_effect[1] -0.08891000 0.008974062 0.23486123
## trt_effect[2] -0.27911215 0.008524743 0.24769692
## trt_effect[3]  0.02702020 0.008458358 0.24551768
## trt_effect[4] -0.30815815 0.008808810 0.23819128
## sigmaCT        2.82548131 0.005607236 0.15244860
## t_dof          4.31963425 0.032490185 0.89971666
## gamma_rnasep   0.06180949 0.033906488 0.95703557
## sigma_plate    0.47180664 0.004422641 0.13216581
## sigmasq_u[1]   3.02138770 0.014120854 0.40460707
## sigmasq_u[2]   0.51894263 0.003717863 0.09512628
```

```
##                      mean     se_mean         sd
## A0            18.31985034 0.047624662 1.18095757
## alpha         -1.50830118 0.024970259 0.52053121
## trt_effect[1] -0.10875017 0.009377879 0.21786808
## trt_effect[2] -0.16257121 0.009581008 0.23245511
## trt_effect[3]  0.09959289 0.015913977 0.25984214
## trt_effect[4] -0.25535884 0.015090239 0.23688704
## sigmaCT        2.77592102 0.005680938 0.15580236
## t_dof          4.25974008 0.031165005 0.87990725
## gamma_rnasep   0.61795092 0.006084076 0.15011465
## sigma_plate    0.46113897 0.004355924 0.13052445
## sigmasq_u[1]   2.78633554 0.013974384 0.37964871
## sigmasq_u[2]   0.49012577 0.003639807 0.09168617
```

![](SAP_files/figure-html/summary-1.png)<!-- -->![](SAP_files/figure-html/summary-2.png)<!-- -->![](SAP_files/figure-html/summary-3.png)<!-- -->


## Results

### Estimated treatment effects: no covariate adjustment

Posterior distributions over the treatment effects for the interventions. Red: no effect; blue: median inferred effect.


```
##   Regeneron  Ivermectin  Remdesivir Favipiravir 
##   0.8838266   0.8317537   1.0359359   0.7483226
```

```
## Probability that Regeneron increases the viral clearance slope by more than 5%: 0.21
```

```
## Probability that Ivermectin increases the viral clearance slope by more than 5%: 0.15
```

```
## Probability that Remdesivir increases the viral clearance slope by more than 5%: 0.48
```

```
## Probability that Favipiravir increases the viral clearance slope by more than 5%: 0.06
```

![](SAP_files/figure-html/treatment_effects-1.png)<!-- -->


### Estimated treatment effects: with covariate adjustment


Posterior distributions over the treatment effects for the interventions, after adjustment for the key covariates (. Red: no effect; blue: median inferred effect.


```
##   Regeneron  Ivermectin  Remdesivir Favipiravir 
##   0.8969545   0.8499556   1.1047211   0.7746385
```

```
## Probability that Regeneron increases the viral clearance slope by more than 5%: 0.26
```

```
## Probability that Ivermectin increases the viral clearance slope by more than 5%: 0.18
```

```
## Probability that Remdesivir increases the viral clearance slope by more than 5%: 0.58
```

```
## Probability that Favipiravir increases the viral clearance slope by more than 5%: 0.1
```

![](SAP_files/figure-html/treatment_effects_cov_model-1.png)<!-- -->

### Estimated covariate effects

#### RNaseP

Gamma is the parameter

![](SAP_files/figure-html/unnamed-chunk-3-1.png)<!-- -->



#### Other covariates

![](SAP_files/figure-html/cov_effects-1.png)<!-- -->

### Individual fits

Individual plots. Thick line no covariate adjustment; dashed line covariate adjustment.


![](SAP_files/figure-html/individ_fits-1.png)<!-- -->![](SAP_files/figure-html/individ_fits-2.png)<!-- -->![](SAP_files/figure-html/individ_fits-3.png)<!-- -->![](SAP_files/figure-html/individ_fits-4.png)<!-- -->


## Sensitivity analysis

### Left vs right tonsil

![](SAP_files/figure-html/left_versus_right-1.png)<!-- -->![](SAP_files/figure-html/left_versus_right-2.png)<!-- -->![](SAP_files/figure-html/left_versus_right-3.png)<!-- -->![](SAP_files/figure-html/left_versus_right-4.png)<!-- -->

### Basic slope - individually fit

Use left censored linear regression

![](SAP_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](SAP_files/figure-html/unnamed-chunk-4-2.png)<!-- -->
