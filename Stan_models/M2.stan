//   Copyright 2022 James Watson, Mahidol Oxford Tropical Medicine Research Unit
//
//   Licensed under a CC-BY License
/**
Model 2 with the following characteristics:
- analysis on the copies per ml scale (batch effect adjustement done without uncertainty propagation)
- Adjustment for RNaseP
- Covariate adjustment

**/



data {
  // Patient data
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  real obs_day[Ntot];                          // Time since randomisation for sample
  real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
  real log10_cens_vl[Ntot];                    // censoring value for censored observation
  real RNaseP[Ntot];                           // Scaled RNaseP CT values (mean 0)
  int<lower=1> K_trt;                          // Number of treatment arms
  int<lower=1,upper=K_trt+1> trt[Ntot];        // Trt index (negative control is 0)
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
  cholesky_factor_corr[2] L_Omega;     // correlation matrix
  vector<lower=0>[2] sigmasq_u;        // variance of random effects
  vector<lower=0>[2] sigmasq_u2;       // variance of epoch random effects

  // Measurement error
  real<lower=0> sigmaCT;

  // Population parameters
  real alpha_0;     // population intercept
  real beta_0;                        // population slope
  real gamma_rnasep;                  // Adjustment for RNaseP
  vector[K_trt] trt_effect;           // Estimates of the treatment effect
  vector[K_cov] slope_coefs;          // Covariate effects on slope
  vector[K_cov] intercept_coefs;      // Covariate effects on intercept

  // Random effects
  vector[2] theta_rand_id[n_id];      // individual random effects vector
  vector[K_epoch-1] theta_epoch[2];   // intercept and slope random effects for the temporal epochs of the trial

  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[K_trt+1] trt_effect_prime;
  vector[K_epoch-1] theta_epoch_prime[2];
  vector[Ntot] beta_cov;
  vector[Ntot] alpha_cov;

  for(i in 1:2){
    theta_epoch_prime[i] = append_row(0, theta_epoch[1]);
  }
  trt_effect_prime = append_row(0, trt_effect);

  // make individual covariate transform
  beta_cov = x*slope_coefs;
  alpha_cov = x*intercept_coefs;

  // calculate predicted log viral load under the model parameters
  for(i in 1:Ntot){
    pred_log10_vl[i] =
    alpha_0 + theta_rand_id[id[i]][1] + theta_epoch_prime[1][epoch[i]] + gamma_rnasep*RNaseP[i] + alpha_cov[i] +
    beta_0*exp(trt_effect_prime[trt[i]]+theta_rand_id[id[i]][2]+theta_epoch_prime[2][epoch[i]]+beta_cov[i])*obs_day[i];
  }
}

model {
  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance (models 0-2)
  sigmaCT ~ normal(3,3) T[0,];

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
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);
  trt_effect ~ normal(0,sigma_trt_effect);  // Treatment effect
  slope_coefs ~ normal(0,1);
  intercept_coefs ~ normal(0,1);

  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigmaCT);

  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigmaCT);
  }
}

generated quantities {
  real preds[Ntot]; // For plotting
  vector[Ntot] log_lik;

  for(i in 1:N_obs){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigmaCT);
  }
  for(i in (N_obs+1):Ntot){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigmaCT);
  }
}
