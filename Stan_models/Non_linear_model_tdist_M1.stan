data {
  // Patient PCR data
  int<lower=0> N;                          // Number of PCR data points
  int<lower=0> n_id;                       // Number of individuals
  int<lower=1,upper=n_id> id[N];           // Patient identifier for each PCR sample
  real<lower=0> obs_day[N];                // Time since randomisation for sample
  real<lower=0,upper=40> delta_CT[N];      // 40-CT value
  real RNaseP[N];                          // Scaled RNaseP CT values (mean 0)
  int<lower=1> K_trt;                      // Number of treatment arms
  int<lower=1,upper=K_trt+1> trt[N];       // Trt index (negative control is 0)
  int<lower=1> K_plate;                    // Number of plates (batch) used for the PCR data
  int<lower=1,upper=K_plate> id_plate[N];  // Plate identifier
  int<lower=0> K_cov;                      // number of columns in covariate design matrix
  matrix[N,K_cov] x;                       // covariate design matrix

  // PCR quality control data
  int<lower=1> N_control;                                    // Number of control PCR samples data
  real control_density[N_control];                           // Known viral density (log10 viral copies per ml)
  real<lower=0,upper=40> control_delta_CT[N_control];        // Observed 40-CT value for the controls
  int<lower=1,upper=K_plate> cont_id_plate[N_control];       // Plate identifier

  // priors
  real alpha_0_prior_mean;
  real alpha_0_prior_sd;

  real beta_0_prior_mean;
  real beta_0_prior_sd;

  real sigma_trt_effect;
}

transformed data {
  vector[4] zeros4;
  for(i in 1:4) zeros4[i] = 0;
}

parameters {
    // hyperparameters
  cholesky_factor_corr[4] L_Omega; // correlation matrix
  vector<lower=0>[4] sigmasq_u; // variance of random effects

  // Measurement error
  real<lower=0> sigmaCT;

  // Population parameters
  real<lower=0,upper=40> alpha_0;          // population intercept
  vector<lower=0>[2] beta_0;               // growth and decline rates
  real gamma_rnasep;                       // Adjustment for RNaseP

  // Random effects
  vector[4] theta_rand_id[n_id];      // individual random effects vector
  vector[K_trt] trt_effect;           // Estimates of the treatment effect

  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_CT[N];
  vector[K_trt+1] trt_effect_prime;

  trt_effect_prime[1]=0; // no study drug arm
  for(i in 1:K_trt){
    trt_effect_prime[i+1]=trt_effect[i];
  }

  for(i in 1:N){
    real intercept = alpha_0+theta_rand_id[id[i]][1];
    real a = beta_0[1]*exp(theta_rand_id[id[i]][2]);
    real b = beta_0[2]*exp(theta_rand_id[id[i]][3]+trt_effect_prime[trt[i]]);
    real tmax = theta_rand_id[id[i]][4];
    pred_CT[i] = intercept+gamma_rnasep*RNaseP[i]+log(a+b)-
    log_sum_exp(log(b)-a*(obs_day[i]-tmax),log(a)+b*(obs_day[i]-tmax));
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
  sigmasq_u[2] ~ normal(1,1);
  sigmasq_u[3] ~ exponential(1);
  sigmasq_u[4] ~ normal(1,1);

  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(sigmasq_u, L_Omega));

  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0[1] ~ normal(3,3);
  beta_0[2] ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);

  // Treatment effect
  trt_effect ~ normal(0,sigma_trt_effect);

  // Main model specification:
  for(i in 1:N){
    if(delta_CT[i]>0){
      target += student_t_lpdf(delta_CT[i] | t_dof, pred_CT[i], sigmaCT);
    } else {
      target += student_t_lpdf(0 | t_dof, pred_CT[i], sigmaCT);
    }
  }
}

generated quantities {
  real pred_CT_mean[N]; // For plotting
  vector[N] log_lik;

  for(i in 1:N){
    real intercept = alpha_0+theta_rand_id[id[i]][1];
    real a = beta_0[1]*exp(theta_rand_id[id[i]][2]);
    real b = beta_0[2]*exp(theta_rand_id[id[i]][3]+trt_effect_prime[trt[i]]);
    real tmax = theta_rand_id[id[i]][4];

    pred_CT_mean[i] = intercept+log(a+b)-
    log_sum_exp(log(b)-a*(obs_day[i]-tmax),log(a)+b*(obs_day[i]-tmax));

    if(delta_CT[i]>0){
      log_lik[i] = student_t_lpdf(delta_CT[i] | t_dof, pred_CT[i], sigmaCT);
    } else {
      log_lik[i] = student_t_lcdf(0 | t_dof, pred_CT[i], sigmaCT);
    }
  }
}

