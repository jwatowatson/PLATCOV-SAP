
library(rstan)

load('../Rout/model_fits_2.RData')


# get the individual slope estimates
thetas = extract(out)
N=length(thetas$gamma_rnasep)
j_sample = sample.int(n = N,size = 1)
source('functions.R')
f_sim(t_design = sort(rep(0:7, 2)),A0 = thetas$alpha_0[j]
