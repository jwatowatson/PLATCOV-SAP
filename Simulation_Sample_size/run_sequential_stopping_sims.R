args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file


### Set up
library(rstan)
library(matrixStats)
library(doParallel)
rstan_options(auto_write = TRUE)
source('sample_size_functions.R')
source('priors.R')

# use the linear model fits for simplicity
load('Rout/model_fits1.RData')
# get the individual slope estimates
thetas = rstan::extract(out); rm(out)
my_LOD = 1
ncores = 4
options(mc.cores = ncores)

load('Rout/sequential_sim_params.RData')
if(i > nrow(sim_settings)){ stop() }

trt_effects = c(sim_settings$Trt_effect_neg_control[i],
                sim_settings$Trt_effect_pos_control[i],
                sim_settings$intervention_effect[i])
Trt_vector = c(rep(1, sim_settings$Nmin[i]),
               rep(2, sim_settings$Nmin[i]),
               rep(3, sim_settings$Nmin[i]))

day_sampling <- as.numeric(unlist(str_split(sim_settings$day_plans[i], ",")))
t_design <- sort(rep(day_sampling,sim_settings$N_swabs_per_day[i]))

current_data = sim_individuals(thetas = thetas,
                               t_design = t_design,
                               N = length(Trt_vector),
                               trt_effects = trt_effects,
                               Trt_arm = Trt_vector,
                               LOD = my_LOD,
                               f_sim = f_sim)
Ncurrent = sim_settings$Nmin[i]
mod = stan_model(file = 'Linear_model_basic.stan') # compile linear model
futility = F; N_futility_success=NA;
success = F;
non_inferiority = F; N_non_inferiority=NA;
inferiority = F; N_inferiority=NA;

stop_trial = F

while( !stop_trial & Ncurrent<=sim_settings$Nmax[i]){

  # comparison with negative control arm
  if(!futility & !success){
    writeLines(sprintf('Running comparison with negative control with %s patients per group...',Ncurrent))
    ind = current_data$Trt_arm %in% c(1,3)
    sim_vl = current_data[ind, ]
    sim_vl$Trt = factor(sim_vl$Trt_arm,levels=c(1,3))

    sim_vl$Censored = as.numeric(sim_vl$log10_viral_load==sim_vl$log10_cens_vl)
    sim_vl = dplyr::arrange(sim_vl, Censored, ID)
    analysis_data=make_stan_inputs(input_data_fit = sim_vl,
                                   trt_frmla = as.formula('~Trt'),
                                   Dmax = 8)
    stan_out = sampling(mod,
                        data=c(analysis_data,
                               all_priors[['WIP']]),
                        iter=2000,
                        chain=4,
                        thin=4,
                        warmup=1000,
                        save_warmup = FALSE,
                        pars=c('trt_effect'),
                        include=T, verbose=F, refresh= 0)
    trt_bounds=quantile(extract(stan_out, pars='trt_effect')$trt_effect, probs = c(0.1, 0.9))
    writeLines(sprintf('Trt effect relative to negative control estimated between %s%% and %s%%',
                       round(100*(exp(trt_bounds[1])-1),1),
                       round(100*(exp(trt_bounds[2])-1),1)))

    if(trt_bounds[1] > sim_settings$Futility_delta[i]){
      success = T
      writeLines(sprintf('Success with %s patients per group!',Ncurrent))
      N_futility_success=Ncurrent

    } else if (trt_bounds[2] < sim_settings$Futility_delta[i]){
      futility = T
      writeLines(sprintf('Futility with %s patients per group!',Ncurrent))
      N_futility_success=Ncurrent
    }
    
    trt_effs <- quantile(extract(stan_out, pars='trt_effect')$trt_effect, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
    
  }

  # comparison with positive control arm
  if(success & !non_inferiority & Ncurrent>=sim_settings$Nmin_NI){
    writeLines(sprintf('Running comparison with positive control with %s patients per group...',Ncurrent))

    ind = current_data$Trt_arm %in% c(2,3)
    sim_vl = current_data[ind, ]
    sim_vl$Trt = factor(sim_vl$Trt_arm,levels=c(2,3))

    sim_vl$Censored = as.numeric(sim_vl$log10_viral_load==sim_vl$log10_cens_vl)
    sim_vl = dplyr::arrange(sim_vl, Censored, ID)
    analysis_data=make_stan_inputs(input_data_fit = sim_vl,
                                   trt_frmla = as.formula('~Trt'),
                                   Dmax = 8)
    stan_out = sampling(mod,
                        data=c(analysis_data,
                               all_priors[['WIP']]),
                        iter=2000,
                        chain=4,
                        thin=4,
                        warmup=1000,
                        save_warmup = FALSE,
                        pars=c('trt_effect'),
                        include=T, verbose=F, refresh= 0)
    trt_bounds=quantile(extract(stan_out, pars='trt_effect')$trt_effect, probs = c(0.1, 0.9))
    writeLines(sprintf('Trt effect relative to positive control estimated between %s%% and %s%%',
                       round(100*(exp(trt_bounds[1])-1),1),
                       round(100*(exp(trt_bounds[2])-1),1)))

    if(trt_bounds[2] < -sim_settings$NI_delta[i]){
      inferiority = T
      writeLines(sprintf('Inferiority with %s patients per group!',Ncurrent))
      N_inferiority=Ncurrent
    } else if (trt_bounds[1] > -sim_settings$NI_delta[i]){
      non_inferiority = T
      writeLines(sprintf('Non-inferiority with %s patients per group!',Ncurrent))
      N_non_inferiority=Ncurrent
    }
  }
  
  
  if(any(c(futility, non_inferiority, inferiority))){
    stop_trial=T
  } else {
    writeLines(sprintf('Adding %s patients per group...', sim_settings$Nseq[i]))
    Trt_vector = c(rep(1, sim_settings$Nseq[i]),
                   rep(2, sim_settings$Nseq[i]),
                   rep(3, sim_settings$Nseq[i]))
    next_data = sim_individuals(thetas = thetas,
                                t_design = t_design,
                                N = length(Trt_vector),
                                trt_effects = trt_effects,
                                Trt_arm = Trt_vector,
                                LOD = my_LOD,
                                f_sim = f_sim)
    current_data = rbind(current_data, next_data)
    Ncurrent = Ncurrent+sim_settings$Nseq[i]
  }
}

out_sim = data.frame(t(unlist(c(sim_settings[i,,drop=F], stop_trial,
                                futility, success, non_inferiority, inferiority,
                                N_futility_success,N_inferiority, N_non_inferiority))), trt_effs)
colnames(out_sim) = c(colnames(sim_settings),
                      'stop_trial',
                      'futility', 'success', 'non_inferiority', 'inferiority',
                      'N_futility_success', 'N_inferiority', 'N_non_inferiority', 
                      'trt_l95', "trt_l90", "trt_med", "trt_u90", "trt_u95")

write.csv(x = out_sim, file = paste0('Sim_out/sim_sequential_',i,'.csv'),row.names = F)
