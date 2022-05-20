prior_params = list(alpha_0_prior_mean = 5,
                          alpha_0_prior_sd = 2,
                          beta_0_prior_mean = -.5,
                          beta_0_prior_sd = 1,
                          trt_effect_sd = .5,
                          sigma_logvl_mean = 1,
                          sigma_logvl_sd = 1,
                          gamma_rnasep_sd = 1,
                          slope_coefs_sd = 1,
                          intercept_coefs_sd = 1,
                          t_max_pop_mean = -3,
                          t_max_pop_sd = 3)

prior_params_sensitive = prior_params
for(i in grep('sd',x = names(prior_params))){
  prior_params_sensitive[[i]] = prior_params_sensitive[[i]]*10
}

all_priors = list(WIP=prior_params, NIP=prior_params_sensitive)