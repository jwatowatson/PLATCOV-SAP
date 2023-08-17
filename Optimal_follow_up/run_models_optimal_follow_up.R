# script for running all the stan models with all settings on the BMRC cluster

#args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
#i = as.numeric(args[6])
#print(paste0("job(i) = ", i)) # this will print out in the *.o file
############################################################################
## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)
library(tidyverse)
library(kableExtra)
library(finalfit)
library(RColorBrewer)
library(lubridate)
library(brms)

i <- 35
############################################################################
# Setting directories
user <- "Chang"   # Make it as args[]

# Main directory
if(user == "Chang"){
  mainDir <- "D:/PLATCOV-SAP"
}else{
  mainDir <- ""}

setwd(mainDir)
############################################################################
source('functions.R')
############################################################################
prep_study_time <- function(platcov_dat){
  platcov_dat = platcov_dat %>% group_by(ID, Timepoint_ID) %>%
    mutate(daily_VL = mean(log10_viral_load)) %>% 
    ungroup() %>% mutate(Sex = as.factor(ifelse(Sex==1,'Male','Female')),
                         Site = as.factor(Site),
                         Trt = factor(Trt, levels=c(ref_arm, trts)),
                         Vaccinated = as.factor(ifelse(N_dose>0,'Yes','No')),
                         Variant = as.factor(Variant),
                         trt_color = 
                           as.character(plyr::mapvalues(Trt,
                                                        from = names(trt_colors),
                                                        to = trt_colors)),
                         Study_time = as.numeric(difftime(Rand_date,min(Rand_date),units = 'weeks')),
                         Study_time = scale(Study_time) #normalise
    )
  return(as.data.frame(platcov_dat))
}

dat_prep_for_analysis <- function(platcov_dat, Dmax){
  platcov_dat_analysis <- 
    platcov_dat %>% ungroup() %>%
    filter(Time <= Dmax+1, mITT) %>%
    arrange(log10_viral_load==log10_cens_vl) %>%
    mutate(Variant = as.factor(Variant),
           Epoch = paste(month(Rand_date), year(Rand_date), sep = '_'),
           Site = as.factor(Site),
           RnaseP_scaled = scale(40 - CT_RNaseP,scale = F),
           Mean_age = mean(Age[!duplicated(ID)]),
           SD_age = sd(Age[!duplicated(ID)]),
           Age_scaled = (Age-Mean_age)/SD_age,
           Symptom_onset = ifelse(is.na(Symptom_onset),2,Symptom_onset)) 
  return(as.data.frame(platcov_dat_analysis))
}
############################################################################
covs_base = c('Variant','Site','Study_time')
covs_full=c(covs_base, 'Age_scaled','Symptom_onset')
add_epoch = T # if using non-concurrent controls
############################################################################
load('Optimal_follow_up/model_settings.RData')

Max_job = nrow(model_settings)
if(i > Max_job) stop('no model setting corresponding to job ID')

writeLines('Doing the following job:')
print(model_settings[i, ])

options(mc.cores = model_settings$Nchain[i])
stopifnot(model_settings$Nchain[i]>getDoParWorkers()) # check worker number assigned

mod = stan_model(file = as.character(model_settings$mod[i])) # compile 
############################################################################
platcov_dat_analysis <- data_list[[model_settings$data_ID[i]]]
#platcov_dat_analysis$Trt[platcov_dat_analysis$Trt == "Nirmatrelvir + Ritonavir"] <- "Nirmatrelvir"
Dmax <- model_settings$Dmax[i]

ref_arm <- model_settings$ref_arm[i]
trts <- model_settings$intervention[i]
trt_colors <- get_trt_colors()

platcov_dat_analysis <- platcov_dat_analysis[platcov_dat_analysis$Trt %in% c(ref_arm, trts),]
platcov_dat_analysis <- prep_study_time(platcov_dat_analysis)
platcov_dat_analysis <- dat_prep_for_analysis(platcov_dat_analysis, Dmax)


stan_input_job <- make_stan_inputs(input_data_fit = platcov_dat_analysis,
                   int_covs_base = c(covs_base,'Symptom_onset'),
                   int_covs_full = covs_full,
                   slope_covs_base = covs_base,
                   slope_covs_full = covs_full,
                   trt_frmla = formula('~ Trt'),
                   epoch = add_epoch,
                   Dmax = Dmax+1)

analysis_data_stan = stan_input_job$analysis_data_stan
analysis_data_stan$trt_mat = stan_input_job$Trt_matrix
analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)

x_intercept = stan_input_job$cov_matrices$X_int[[model_settings$cov_matrices[i]]]
if(ncol(x_intercept)==0) x_intercept = array(0, dim=c(nrow(x_intercept),1))
analysis_data_stan$x_intercept = x_intercept
analysis_data_stan$K_cov_intercept= ncol(x_intercept)


x_slope = stan_input_job$cov_matrices$X_slope[[model_settings$cov_matrices[i]]]
if(ncol(x_slope)==0) x_slope = array(0, dim=c(nrow(x_slope),1))
analysis_data_stan$x_slope = x_slope
analysis_data_stan$K_cov_slope=ncol(x_slope)

# sample posterior
out = sampling(mod, 
               data=c(analysis_data_stan,
                      all_priors[[model_settings$prior[i]]]),
               iter=model_settings$Niter[i],
               chain=model_settings$Nchain[i],
               thin=model_settings$Nthin[i],
               warmup=model_settings$Nwarmup[i],
               save_warmup = FALSE,
               seed=i,
               pars=c('L_Omega'), # we don't save this as it takes up lots of memory!
               include=FALSE)

save(out, file = paste0('output/model_fits_',i,'.RData'))# save output

writeLines('Finished job')










