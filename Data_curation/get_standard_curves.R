library(rstan)
library(brms)
library(readr)
options(mc.cores = parallel::detectCores())

control_dat=read_csv(file = "../Analysis_Data/interim_control_dat.csv")
sort(table(control_dat$`Lot no.`))
# random slope and intercept
conv_mod = brm(log10_true_density ~ 1 + CT*Lab + (1+CT|Plate),
                data = control_dat)
summary(conv_mod)
