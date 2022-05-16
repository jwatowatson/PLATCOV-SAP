# script for running all the stan models with all settings on the BMRC cluster

args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file

set.seed(seed = i)

## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)

load('Rout/model_run_setup.RData')

Max_job = nrow(model_settings)
if(i > Max_job) stop('no model setting corresponding to job ID')

writeLines('Doing the following job:')
print(model_settings[i, ])

options(mc.cores = model_settings$Nchain[i])
stopifnot(model_settings$Nchain[i]>getDoParWorkers()) # check worker number assigned

mod = stan_model(file = as.character(model_settings$mod[i])) # compile 

analysis_data_stan$trt_mat = Trt_mats[[model_settings$trt_mat[i]]]
analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)
out = sampling(mod, # sample posterior 
               data=c(analysis_data_stan, 
                      all_priors[[model_settings$prior[i]]]), 
               iter=model_settings$Niter[i], 
               chain=model_settings$Nchain[i],
               thin=model_settings$Nthin[i], 
               warmup=model_settings$Nwarmup[i],
               save_warmup = FALSE, 
               pars=c('L_Omega','theta_rand_id'),
               include=FALSE)

save(out, file = paste0('Rout/model_fits_',i,'.RData'))# save output

writeLines('Finished job')

