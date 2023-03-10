args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file

library(rstan)
library(matrixStats)
library(doParallel)
rstan_options(auto_write = TRUE)
source('sample_size_functions.R')
source('priors.R')

# use the linear model fits for simplicity
load('Rout/linear_fit.RData')
# get the individual slope estimates
thetas = extract(out); rm(out)
my_LOD = 1
ncores = 4
options(mc.cores = ncores)


## Set up all simulation settings
Trt_effect_pos_control = 1.6
Trt_effect_neg_control = 1

trt_intervention = seq(from=1.2, to=1.6, by=.1)
NI_delta = log(1.6)-log(1.4)

Ns = c(25,50,75,100)
Nsims = 100
sim_settings = expand.grid(N=Ns,
                           sim_k = 1:Nsims,
                           trt_effect_comp = trt_intervention,
                           trt_control = c(Trt_effect_pos_control,
                                           Trt_effect_neg_control))


save(sim_settings, file = 'sim_settings.RData')

### set up simulation for the settings i
print(sim_settings[i, ])

trt_effects = c(sim_settings$trt_control[i],
                sim_settings$trt_effect_comp[i])
Trt_vector = c(rep(1, sim_settings$N[i]),
               rep(2, sim_settings$N[i]))
sim_vl = sim_individuals(thetas = thetas,
                         t_design = rep(0:7,2),
                         N = length(Trt_vector),
                         trt_effects = trt_effects,
                         Trt_arm = Trt_vector,
                         LOD = my_LOD,
                         f_sim = f_sim)

sim_vl$Trt = factor(sim_vl$Trt_arm,levels=1:2)

sim_vl$Censored = as.numeric(sim_vl$log10_viral_load==sim_vl$log10_cens_vl)
sim_vl = dplyr::arrange(sim_vl, Censored, ID)
analysis_data=make_stan_inputs(input_data_fit = sim_vl,
                               trt_frmla = as.formula('~Trt'),
                               Dmax = 8)

# xx = aggregate(log10_viral_load~Time+Trt, data = sim_vl, median)
# ind = xx$Trt==1
# plot(xx$Time[ind], xx$log10_viral_load[ind],type='l',col='blue',lwd=3,ylim=c(1,6))
# lines(xx$Time[!ind], xx$log10_viral_load[!ind], col='red', lwd=3)

mod = stan_model(file = 'Linear_model_basic.stan') # compile

stan_out = sampling(mod,
                    data=c(analysis_data,
                           all_priors[['WIP']]),
                    iter=2000,
                    chain=4,
                    thin=4,
                    warmup=1000,
                    save_warmup = FALSE,
                    pars=c('trt_effect'),
                    include=T, verbose=F)

out_trt_effect = as.data.frame(rstan::extract(stan_out, pars='trt_effect'))

f_name = paste0('sims_out/sim_out',i,'.csv')
print(f_name)
write.csv(out_trt_effect, file = f_name,row.names = F)
