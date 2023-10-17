## Set up all simulation settings
Trt_effect_pos_control = 1.6 ## positive control
Trt_effect_neg_control = 1 ## negative control

trt_intervention = c(0.8, 1, 1.2, 1.4, 1.6)
NI_delta = c(0.1)
Futility_delta = log(1.2)

Nmax = 120
Nseq = 10
Nmin = 20
Nmin_NI = 40

Nsims = 50

day_plans <- c(paste(as.character(0:6),collapse = ','), 
               paste(as.character(c(0,2,4,6)),collapse = ','), 
               paste(as.character(c(0,3,6)),collapse = ','), 
               paste(as.character(c(0,6)),collapse = ',')) 

N_swabs_per_day <- c(2,4)

sim_settings = expand.grid(sim_k = 1:Nsims,
                           intervention_effect = trt_intervention,
                           Trt_effect_pos_control = Trt_effect_pos_control,
                           Trt_effect_neg_control = Trt_effect_neg_control,
                           NI_delta = NI_delta,
                           Futility_delta = Futility_delta,
                           Nmax = Nmax,
                           Nmin = Nmin,
                           Nseq = Nseq,
                           Nmin_NI = Nmin_NI,
                           day_plans = day_plans,
                           N_swabs_per_day = N_swabs_per_day)

save(sim_settings, file = 'Rout/sequential_sim_params.RData')


