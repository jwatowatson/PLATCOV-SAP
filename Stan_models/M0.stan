//   Copyright 2022 James Watson, Mahidol Oxford Tropical Medicine Research Unit
//
//   Licensed under a CC-BY License
/**
Model 0 with the following characteristics:
- analysis on the copies per ml scale (batch effect adjustement done without uncertainty propagation)
- no adjustment for RNaseP
- no covariate adjustment

**/



data {
  // Patient data
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  int<lower=1,upper=Ntot> ind_start[n_id];     // Starting index for each patient
  real obs_day[Ntot];                          // Time since randomisation for sample
  real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
  real log10_cens_vl[Ntot];                    // censoring value for censored observation
  real RNaseP[Ntot];                           // Scaled RNaseP CT values (mean 0)
  int<lower=1> K_trt;                          // Number of treatment arms
  matrix[Ntot,K_trt] trt_mat;                  // Trt matrix
  int<lower=0> K_cov;                          // number of columns in covariate design matrix
  matrix[Ntot,K_cov] x;                        // covariate design matrix
  int<lower=1> K_epoch;                        // number of trial epochs (time periods for temporal drift)
  int<lower=1,upper=K_epoch> epoch[Ntot];      // trial epochs

  // priors
  real alpha_0_prior_mean; // prior mean intercept
  real alpha_0_prior_sd;   // prior sd intercept

  real beta_0_prior_mean;  // prior mean slope
  real beta_0_prior_sd;    // prior sd slope

  real sigma_trt_effect;   // prior sd on treatment effect
}

transformed data {
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;     // correlation matrix for individ RE
  vector<lower=0>[2] sigmasq_u;        // variance of individual random effects
  vector<lower=0>[2] sigmasq_u2;       // variance of epoch random effects

  // Measurement error
  real<lower=0> sigma_logvl;

  // Population parameters
  real alpha_0;                        // population intercept
  real beta_0;                         // population slope
  vector[K_trt] trt_effect;            // Estimates of the treatment effect

  // Random effects
  vector[2] theta_rand_id[n_id];      // individual random effects vector
  vector[K_epoch-1] theta_epoch[2];   // intercept and slope random effects for the temporal epochs of the trial

  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[K_epoch-1] theta_epoch_prime[2];
  vector[Ntot] trt_slope;

  trt_slope = trt_mat * trt_effect;

  for(i in 1:2){
    theta_epoch_prime[i] = append_row(0, theta_epoch[1]);
  }

  // calculate predicted log viral load under the model parameters
  for(i in 1:Ntot){
    pred_log10_vl[i] =
    alpha_0 + theta_rand_id[id[i]][1] + theta_epoch_prime[1][epoch[i]] +
    beta_0*exp(trt_slope[i]+theta_rand_id[id[i]][2]+theta_epoch_prime[2][epoch[i]])*obs_day[i];
  }
}

model {
  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance
  sigma_logvl ~ normal(1.5,3) T[0,];

  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  sigmasq_u2 ~ normal(0,0.5);
  L_Omega ~ lkj_corr_cholesky(2); // covariance matrix - random effects for individs
  // individual random effects
  for(i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));
  // epoch random effects (independent)
  for(i in 1:2) {
    theta_epoch[i] ~ normal(0, sigmasq_u2[i]);
  }

  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd); // intercept
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);    // slope
  trt_effect ~ normal(0,sigma_trt_effect);               // Treatment effect

  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_logvl);

  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
}

generated quantities {
  real preds[Ntot]; // For plotting
  vector[Ntot] log_lik;
  vector[n_id] slope;

  for(i in 1:N_obs){
    preds[i] = pred_log10_vl[i];
    log_lik[i] = student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in (N_obs+1):Ntot){
    preds[i] = pred_log10_vl[i];
    log_lik[i] = student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in 1:n_id){
    int j =ind_start[i];
    slope[i] = beta_0*exp(trt_slope[j]+theta_rand_id[id[j]][2]+theta_epoch_prime[2][epoch[j]]);
  }
}
