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

* Up and down model to account for the fact that subjects are in an increasing phase
* How much variance explained by the RNaseP
* Check missing data (vaccine/serology/age)
* Subgroup analysis model for temporal shift
* Define Epochs for the trial to get drift effects
* Model comparison using loo package

## Preambule

This Statistical Analysis Plan (SAP) is written for the PLATCOV trial (registered at ClinicalTrials.gov: https://clinicaltrials.gov/ct2/show/NCT05041907).

Data preparation is done in a different R script called data_prep.R.
This Markdown script assumes that the data are saved in two .csv files in long format. 

The first file *interim_control_dat.csv* contains the PCR quality control data: CT values from the control samples from each plate (duplicate control samples of known viral densities over the range 100 to 10^7/mL viral copies per mL). This has the following headers:

- ID: unique control ID
- Plate: unique PCR Plate ID
- CT_NS: observed CT value for the N/S gene
- log10_true_density: number of viral copies in control sample (log10 scale)

The second file *interim_dat.csv* contains the patient viral load data with the following column headers:

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
- Antibody_test: 0/1 (negative/positive for SARS-CoV-2 antibody rapid test)
- Age: (years - has to be between 18-50)
- Sex: 0/1 (male: 1; female/other: 0)
- Symptom_onset: time since onset of symptoms (days)
- CT_NS: observed CT value for the N/S gene
- CT_RNaseP: observed CT value for the human RNase P gene
- log10_viral_load: log10 number of viral copies per mL (estimated from plate specific standard curve - 12 controls per plate)


## Computational setup


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


## Data summaries


```
## We have 720 PCR datapoints on 36 patients from 1 sites
```

```
##                Site
## Arm             th001
##   Favipiravir       8
##   Ivermectin        8
##   No study drug     6
##   Regeneron         7
##   Remdesivir        7
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

make stan data set


```
## Total number of datapoints up until day 8 is 648
```

### Run models


We fit a sequence of Bayesian hierarchical models.
To make sure there are no bugs in the code (all stan code is written specifically for this trial analysis), I fit the following sequence of five models of increasing complexity:

* Model 0: vanilla student-t regression with left censoring at 0 and with individual random effects for slope and intercept;
* Model 1: vanilla + RNaseP;
* Model 2: vanilla + RNaseP + batch adjustment (no control samples);
* Model 3: vanilla + RNaseP + full batch adjustment (decomposition of variance into intrinsic and biological);
* Model 4: Whole shebang with covariate adjustment

Covariates that we use in model 4:

* Vaccination (yes/no)
* Time since symptom onset (days)
* Variant (WHO variants of concern)
* Serology rapid test (+/-)



```
## We are running the 5 models with 4 chains and 10000 samples for each chain, discarding half for burn-in and thining every 20, thus giving a total of 1000 posterior samples per model.
```

### Model fits: summaries


```
## 
## *********************
## Summary of posterior distribution for model 0
##                      mean        2.5%      97.5%     n_eff      Rhat
## trt_effect[1]  0.44308082 -0.00180076  0.8714354  952.9113 0.9990750
## trt_effect[2] -0.06177754 -0.53406231  0.4216128 1153.0575 0.9997622
## trt_effect[3]  0.00387721 -0.49797350  0.4455998 1005.2521 0.9991905
## trt_effect[4] -0.10775682 -0.56861489  0.3584960 1027.8453 1.0050767
## A0            15.44159744 14.36535109 16.5062542  959.7849 1.0068817
## alpha         -1.23730018 -1.68386740 -0.8503413  952.5037 1.0003419
## sigmaCT        2.83569788  2.53706717  3.1593501  956.7500 0.9995531
## sigmasq_u[1]   3.03064614  2.29934598  3.9019438  817.1630 1.0001303
## sigmasq_u[2]   0.45507058  0.30244432  0.6446984 1095.1847 0.9973530
## t_dof          4.31540543  2.95864532  6.1661448  969.2608 0.9979292
## 
## *********************
## Summary of posterior distribution for model 1
##                      mean        2.5%      97.5%     n_eff      Rhat
## trt_effect[1]  0.42917194  0.01673384  0.8603816  995.7399 0.9975046
## trt_effect[2] -0.03843522 -0.45730226  0.3800871  999.7579 0.9987647
## trt_effect[3]  0.02806542 -0.40556477  0.4522204 1000.9441 1.0010875
## trt_effect[4] -0.05890801 -0.45317411  0.3474273 1037.8002 0.9981622
## A0            15.58898118 14.52600337 16.6883919 1033.3918 0.9994717
## alpha         -1.30147957 -1.74449176 -0.9494138 1007.9183 1.0017206
## sigmaCT        2.78150242  2.47311686  3.0918452 1012.8477 1.0017545
## sigmasq_u[1]   3.07908485  2.30728444  4.0063816 1086.7112 1.0001451
## sigmasq_u[2]   0.41496641  0.28025855  0.5966383  967.3092 1.0012380
## t_dof          4.32046403  2.91743469  6.3611105  939.3114 1.0035230
## gamma_rnasep   0.62717088  0.34836472  0.9158722 1004.8623 0.9988334
## 
## *********************
## Summary of posterior distribution for model 2
##                      mean        2.5%      97.5%     n_eff      Rhat
## trt_effect[1]  0.42717842  0.01325427  0.8329235 1033.6348 0.9976164
## trt_effect[2] -0.04233842 -0.47931333  0.3733308 1104.2782 0.9982885
## trt_effect[3]  0.03625824 -0.39428338  0.4366160  961.3626 0.9999458
## trt_effect[4] -0.06821906 -0.49206626  0.3463031  912.9963 0.9977148
## A0            15.62342740 14.43635898 16.8187411  934.8939 0.9990589
## alpha         -1.30225475 -1.76308069 -0.9383355  969.4963 0.9994353
## sigmaCT        2.77554889  2.47407162  3.0854830  992.4406 0.9985815
## sigmasq_u[1]   3.07158336  2.34828393  3.9114384  808.9888 0.9988343
## sigmasq_u[2]   0.41391123  0.26934591  0.5792730  982.8577 0.9981743
## t_dof          4.29059589  2.95069185  6.3221138  971.9300 1.0004445
## gamma_rnasep   0.61863745  0.35507570  0.8968539  875.5191 0.9979761
## sigma_plate    0.47771044  0.03794913  1.5720710  844.8719 1.0019391
## 
## *********************
## Summary of posterior distribution for model 3
##                      mean        2.5%      97.5%     n_eff      Rhat
## trt_effect[1]  0.43897238  0.01900602  0.8571180 1027.1734 0.9986244
## trt_effect[2] -0.02884385 -0.46616938  0.3975431  940.2678 0.9999576
## trt_effect[3]  0.03013027 -0.41635841  0.4418609  889.8027 0.9984391
## trt_effect[4] -0.05218531 -0.48637232  0.3379149  906.1705 0.9983845
## A0            15.63841917 14.43655044 16.7789143  877.6204 1.0065842
## alpha         -1.29168424 -1.73583478 -0.9451520  948.3694 0.9981208
## sigmaCT        2.74933070  2.42382888  3.0428063 1037.5583 0.9992601
## sigmasq_u[1]   3.11997295  2.42857146  4.0265604  965.7845 1.0009358
## sigmasq_u[2]   0.41180737  0.27319176  0.5880420  781.9592 1.0016029
## t_dof          4.26577404  2.90378570  6.2628431  915.8585 1.0005223
## gamma_rnasep   0.61592458  0.34441262  0.8851127 1005.4472 1.0001028
## sigma_plate    0.45940123  0.28040078  0.7487150  911.3476 1.0000764
## sc_intercept  -3.70483897 -4.03091746 -3.3758705  968.0614 1.0032535
## sc_slope       3.37529132  3.34239880  3.4102418 1049.2316 0.9988345
## sigma_sc       0.31605388  0.27551989  0.3614771  946.4619 0.9983001
## 
## *********************
## Summary of posterior distribution for model 4
##                          mean        2.5%      97.5%     n_eff      Rhat
## trt_effect[1]      0.45106576  0.05212062  0.8972870  966.2667 0.9990515
## trt_effect[2]     -0.01227150 -0.43283072  0.4065470  635.8972 1.0028169
## trt_effect[3]      0.03427789 -0.42785670  0.4649577 1099.6490 1.0000821
## trt_effect[4]     -0.09134146 -0.56302191  0.3739333  814.8499 0.9992616
## A0                18.29221027 15.69982438 20.5246773 1066.1878 0.9992499
## alpha             -1.18272087 -1.88033974 -0.6843873  993.8584 0.9973995
## sigmaCT            2.75853952  2.46085403  3.0556051 1089.5268 1.0004429
## sigmasq_u[1]       2.79436042  2.07918817  3.6419900 1096.6366 1.0011389
## sigmasq_u[2]       0.42997646  0.29431927  0.6150155 1002.3195 0.9995922
## t_dof              4.31815260  3.00732452  6.2716177 1055.3510 1.0022417
## gamma_rnasep       0.61211098  0.33784223  0.8938075  954.1244 0.9989681
## sigma_plate        0.45988426  0.28444188  0.7757865  992.8436 1.0010166
## beta_intercept[1] -0.44169795 -1.86244841  1.0631924  825.6571 0.9973102
## beta_intercept[2] -1.31544646 -2.25619151 -0.2028980 1042.9587 0.9998909
## beta_slope[1]      0.13615373 -0.19311431  0.4668086  937.1084 0.9970974
## beta_slope[2]      0.02480117 -0.18293651  0.2218254  986.9162 0.9991707
```

```
## 
## *********************
##  Traceplots for posterior distribution for model 0
```

![](SAP_files/figure-html/summary-1.png)<!-- -->

```
## 
## *********************
##  Traceplots for posterior distribution for model 1
```

![](SAP_files/figure-html/summary-2.png)<!-- -->

```
## 
## *********************
##  Traceplots for posterior distribution for model 2
```

![](SAP_files/figure-html/summary-3.png)<!-- -->

```
## 
## *********************
##  Traceplots for posterior distribution for model 3
```

![](SAP_files/figure-html/summary-4.png)<!-- -->

```
## 
## *********************
##  Traceplots for posterior distribution for model 4
```

![](SAP_files/figure-html/summary-5.png)<!-- -->


## Results

### Estimated treatment effects under the 5 models

Posterior distributions over the treatment effects for the interventions. Red: no effect; blue: median inferred effect.


```
## 
## *******************
## Mean estimated treatment effects (multiplicative):
```

```
##         Regeneron Ivermectin Remdesivir Favipiravir
## Model_0      1.56       0.94       1.00        0.90
## Model_1      1.54       0.96       1.03        0.94
## Model_2      1.53       0.96       1.04        0.93
## Model_3      1.55       0.97       1.03        0.95
## Model_4      1.57       0.99       1.03        0.91
```

```
## 
## *******************
## Probability of super-superiority:
```

```
##         Regeneron Ivermectin Remdesivir Favipiravir
## Model_0      95.5       32.6       42.5        24.1
## Model_1      96.1       35.1       47.2        31.0
## Model_2      95.8       34.4       47.7        28.9
## Model_3      96.6       37.1       48.1        33.2
## Model_4      97.5       38.9       47.7        25.9
```

![](SAP_files/figure-html/treatment_effects-1.png)<!-- -->


### Estimated covariate effects
#### RNaseP

Gamma is the parameter on the RNaseP concentration

![](SAP_files/figure-html/unnamed-chunk-3-1.png)<!-- -->



#### Other covariates


```
## ci_level: 0.8 (80% intervals)
```

```
## outer_level: 0.95 (95% intervals)
```

![](SAP_files/figure-html/cov_effects-1.png)<!-- -->

### Individual fits

Individual plots colored by model.

![](SAP_files/figure-html/individ_fits-1.png)<!-- -->![](SAP_files/figure-html/individ_fits-2.png)<!-- -->![](SAP_files/figure-html/individ_fits-3.png)<!-- -->![](SAP_files/figure-html/individ_fits-4.png)<!-- -->


## Sensitivity analysis

### Left vs right tonsil

![](SAP_files/figure-html/left_versus_right-1.png)<!-- -->![](SAP_files/figure-html/left_versus_right-2.png)<!-- -->![](SAP_files/figure-html/left_versus_right-3.png)<!-- -->![](SAP_files/figure-html/left_versus_right-4.png)<!-- -->

### Basic slope - individually fit

Use left censored linear regression

![](SAP_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](SAP_files/figure-html/unnamed-chunk-4-2.png)<!-- -->
