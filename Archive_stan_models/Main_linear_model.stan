data {
  // Patient PCR data
  int<lower=0> N;                          // Number of PCR data points
  int<lower=0> n_id;                       // Number of individuals
  int<lower=1,upper=n_id> id[N];           // Patient identifier for each PCR sample
  real obs_day[N];                         // Time since randomisation for sample
  real log_10_vl[N];                       // log base 10 viral load in copies per mL
  real log10_cens_vl[N];                   // censoring value for censored observation
  real RNaseP[N];                          // Scaled RNaseP CT values (mean 0)
  int<lower=1> K_trt;                      // Number of treatment arms
  int<lower=1,upper=K_trt+1> trt[N];       // Trt index (negative control is 0)

  // priors
  real alpha_0_prior_mean;
  real alpha_0_prior_sd;

  real beta_0_prior_mean;
  real beta_0_prior_sd;

  real sigma_trt_effect;
}

transformed data {
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;     // correlation matrix
  vector<lower=0>[2] sigmasq_u;        // variance of random effects

  // Measurement error
  real<lower=0> sigmaCT;

  // Population parameters
  real<lower=0,upper=40> alpha_0;     // population intercept
  real beta_0;                        // population slope
  real gamma_rnasep;                  // Adjustment for RNaseP

  // Random effects
  vector[2] theta_rand_id[n_id];      // individual random effects vector
  vector[K_trt] trt_effect;           // Estimates of the treatment effect

  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_log10_vl[N];
  vector[K_trt+1] trt_effect_prime;

  trt_effect_prime[1]=0; // no study drug arm
  for(i in 1:K_trt){
    trt_effect_prime[i+1]=trt_effect[i];
  }

  for(i in 1:N){
    pred_log10_vl[i] =
    alpha_0 + theta_rand_id[id[i]][1] + gamma_rnasep*RNaseP[i] +
    beta_0*exp(trt_effect_prime[trt[i]]+theta_rand_id[id[i]][2])*obs_day[i];
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
  L_Omega ~ lkj_corr_cholesky(2); // covariance matrix - random effects
  for (i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));

  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);

  // Treatment effect
  trt_effect ~ normal(0,sigma_trt_effect);

  //***** Likelihood *****
  // Main model specification:
  for(i in 1:N){
    if(log_10_vl[i]>log10_cens_vl[i]){
      target += student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigmaCT);
    } else {
      target += student_t_lcdf(log10_cens_vl | t_dof, pred_log10_vl[i], sigmaCT);
    }
  }

}

generated quantities {
  real pred_log10_vl_mean[N]; // For plotting
  vector[N] log_lik;

  for(i in 1:N){
    pred_log10_vl_mean[i] =
    alpha_0 + theta_rand_id[id[i]][1] +
    beta_0*exp(trt_effect_prime[trt[i]]+theta_rand_id[id[i]][2])*obs_day[i];

    if(log_10_vl[i]>log10_cens_vl[i]){
      log_lik[i] = student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl_mean[i], sigmaCT);
    } else {
      log_lik[i] = student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl_mean[i], sigmaCT);
    }
  }
}
