args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
job_i = as.numeric(args[6])
print(paste0("job(i) = ", job_i)) # this will print out in the *.o file

library(rstan)
library(matrixStats)
library(doParallel)
library(stringr)
source('sample_size_functions.R')
source('../priors.R')

# use the linear model fits for simplicity
load('Rout/model_fits1.RData')
load('Rout/model_settings.RData')
# get the individual slope estimates
thetas = rstan::extract(out); rm(out)
my_LOD = 1
ncores = 4
options(mc.cores = ncores)


# Set up all simulation settings
Ns = c(25,50,75,100)
Nsims = 50

day_plans <- c(paste(as.character(0:6),collapse = ','), 
               paste(as.character(c(0,2,4,6)),collapse = ','), 
               paste(as.character(c(0,3,6)),collapse = ','), 
               paste(as.character(c(0,6)),collapse = ',')) 
               
N_swabs_per_day <- c(2,4)

sim_settings = expand.grid(intervention = model_settings$intervention,
                           ref_arm = model_settings$ref_arm,
                           N=Ns,
                           sim_k = 1:Nsims,
                           day_plans = day_plans,
                           N_swabs_per_day = N_swabs_per_day,
                           trt_effects = c(0.8, 1, 1.2, 1.4, 1.6))

save(sim_settings, file = 'sim_settings_extended.RData')

### set up simulation for the settings job_i
print(sim_settings[job_i, ])

day_sampling <- as.numeric(unlist(str_split(sim_settings$day_plans[job_i], ",")))

t_design <- sort(rep(day_sampling,sim_settings$N_swabs_per_day[job_i]))

# simulate data
sim_vl = sim_individuals(thetas = thetas,
                         t_design = t_design,
                         N = sim_settings$N[job_i]*2,
                         Trt_arm = c(rep(1, sim_settings$N[job_i]), 
                                     rep(2,sim_settings$N[job_i])),
                         LOD = my_LOD,
                         f_sim = f_sim)

sim_vl$Trt = factor(sim_vl$Trt_arm,levels=1:2)

sim_vl$Censored = as.numeric(sim_vl$log10_viral_load==sim_vl$log10_cens_vl)
sim_vl = dplyr::arrange(sim_vl, Censored, ID)
analysis_data=make_stan_inputs(input_data_fit = sim_vl,
                               trt_frmla = as.formula('~Trt'),
                               Dmax = max(t_design)+1)

# fit model to simulated data
mod = stan_model(file = '../Stan_models/Linear_model_basic.stan') # compile
stan_out = sampling(mod,
                    data=c(analysis_data,
                           all_priors[['WIP']]),
                    iter=model_settings$Niter, #4000
                    chain=model_settings$Nchain, #4
                    thin=model_settings$Nthin, #8
                    warmup=model_settings$Nwarmup, #2000
                    save_warmup = FALSE,
                    pars=c('trt_effect'),
                    include=T, verbose=F)

out_trt_effect = as.data.frame(rstan::extract(stan_out, pars='trt_effect'))

f_name = paste0('sims_out/sim_out_extended_',job_i,'.csv')
print(f_name)
write.csv(out_trt_effect, file = f_name,row.names = F)
