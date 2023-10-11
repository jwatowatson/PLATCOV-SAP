library(rstan)
load('Rout/model_settings.RData')
i=1
print(model_settings[i, ])

options(mc.cores = model_settings$Nchain[i])

mod = stan_model(file = '../Stan_models/Linear_model_RNaseP.stan') # compile

stan_input_job = stan_inputs[[model_settings$dataset[i]]]

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
               seed=i)
save(out, file = 'Rout/linear_fit.RData')
