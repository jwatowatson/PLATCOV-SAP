load("sim_settings_extended.RData")
load('Rout/model_fits1.RData')
####################################################################################
library(stringr)
library(ggplot2)
source('sample_size_functions.R')
source('../priors.R')
####################################################################################
job_i <- 5953

### set up simulation for the settings job_i
print(sim_settings[job_i, ])

day_sampling <- as.numeric(unlist(str_split(sim_settings$day_plans[job_i], ",")))

t_design <- sort(rep(day_sampling,sim_settings$N_swabs_per_day[job_i]))

# simulate data
#set.seed(job_i)
sim_vl = sim_individuals(thetas = thetas,
                         t_design = t_design,
                         N = sim_settings$N[job_i]*2,
                         Trt_arm = c(rep(1, sim_settings$N[job_i]), 
                                     rep(2,sim_settings$N[job_i])),
                         trt_effect = sim_settings$trt_effects[job_i],
                         LOD = my_LOD,
                         f_sim = f_sim)

sim_vl$Trt = factor(sim_vl$Trt_arm,levels=1:2)

sim_vl$Censored = as.numeric(sim_vl$log10_viral_load==sim_vl$log10_cens_vl)
sim_vl = dplyr::arrange(sim_vl, Censored, ID)

############################################################################################
sim_vl2 <- aggregate(log10_viral_load~ID+Time+Trt, 
                       data = sim_vl, FUN = median)
sim_vl3 <- aggregate(log10_viral_load~Time+Trt, 
                     data = sim_vl2, FUN = median)
############################################################################################
ggplot() +
  geom_point(data = sim_vl2, aes(x = Time, y = log10_viral_load, col = Trt), alpha = 0.4, size = 1.75, shape = 16) +
  geom_line(data = sim_vl2, aes(x =  Time, y = log10_viral_load, group = ID, col = Trt), 
            alpha = 0.5, linewidth = 0.5, linetype = 1) +
  scale_color_manual(values = c("#FA9494", "#87CEEB"), name = "") +
  ggnewscale::new_scale_color() +
  geom_point(data = sim_vl3, aes(x = Time, y = log10_viral_load, col = Trt), size = 3, shape = 17) +
  geom_line(data = sim_vl3, aes(x =  Time, y = log10_viral_load, group = Trt, col = Trt), linewidth = 1, linetype = 1) +
  scale_color_manual(values = c("#B31312", "#332FD0"), name = "") +
  theme_bw() +
  xlab("") +
  ylab("Log viral loads") + 
  theme(axis.title  = element_text(face = "bold"))
############################################################################################