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
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega; // correlation matrix
  vector<lower=0>[2] sigmasq_u; // variance of random effects

  // Measurement error
  real<lower=0> sigmaCT;

  // Population parameters
  real<lower=0,upper=40> alpha_0;     // population intercept
  real beta_0;                        // population slope
  real gamma_rnasep;                  // Adjustment for RNaseP
  real sc_intercept;                  // standard curve intercept
  real sc_slope;                      // standard curve slope
  real<lower=0> sigma_sc;             // PCR machine standard deviation
  vector[K_cov] slope_coefs;           // Covariate effects on slope
  vector[K_cov] intercept_coefs;       // Covariate effects on intercept

  // Random effects
  vector[2] theta_rand_id[n_id];      // individual random effects vector
  vector[K_trt] trt_effect;           // Estimates of the treatment effect
  vector[K_plate] a_plate;            // plate random effects for the intercept
  real<lower=0> sigma_plate;          // variance of plate effects

  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_CT[N];
  vector[K_trt+1] trt_effect_prime;
  real sigma_tot;
  vector[N] beta_cov;
  vector[N] alpha_cov;

  // measurement variance (sqrt)
  sigma_tot = sqrt(sigmaCT^2 + sigma_sc^2);

  // make individual covariate transform
  beta_cov = x*slope_coefs;
  alpha_cov = x*intercept_coefs;

  trt_effect_prime[1]=0; // no study drug arm
  for(i in 1:K_trt){
    trt_effect_prime[i+1]=trt_effect[i];
  }
  for(i in 1:N){
    pred_CT[i] =
    alpha_0 + theta_rand_id[id[i]][1] + gamma_rnasep*RNaseP[i] + a_plate[id_plate[i]] +
    alpha_cov[i] +
    beta_0*exp(trt_effect_prime[trt[i]]+theta_rand_id[id[i]][2]+beta_cov[i])*obs_day[i];
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
  a_plate ~ normal(0, sigma_plate);
  sigma_plate ~ exponential(1);

  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);
  sc_intercept ~ normal(-3,3);         // standard curve for control data
  sc_slope ~ normal(log2(10), 1);      // standard curve for control data
  sigma_sc ~ normal(0.5, 1);           // assay variance (sqrt) for control data
  slope_coefs ~ normal(0,1);
  intercept_coefs ~ normal(0,1);

  // Treatment effect
  trt_effect ~ normal(0,sigma_trt_effect);

  //***** Likelihood *****
  // Main model specification:
  for(i in 1:N){
    if(delta_CT[i]>0){
      target += student_t_lpdf(delta_CT[i] | t_dof, pred_CT[i], sigma_tot);
    } else {
      target += student_t_lcdf(0 | t_dof, pred_CT[i], sigma_tot);
    }
  }
  // control data:
  for(i in 1:N_control){
    control_delta_CT[i] ~ normal(sc_intercept+a_plate[cont_id_plate[i]] + sc_slope*control_density[i], sigma_sc);
  }
}

generated quantities {
  real pred_CT_mean[N]; // For plotting
  vector[N] log_lik;

  for(i in 1:N){
    pred_CT_mean[i] =
    alpha_0 + theta_rand_id[id[i]][1] + a_plate[id_plate[i]] + alpha_cov[i] +
    beta_0*exp(trt_effect_prime[trt[i]]+theta_rand_id[id[i]][2]+beta_cov[i])*obs_day[i];

    if(delta_CT[i]>0){
      log_lik[i] = student_t_lpdf(delta_CT[i] | t_dof, pred_CT[i], sigma_tot);
    } else {
      log_lik[i] = student_t_lcdf(0 | t_dof, pred_CT[i], sigma_tot);
    }
  }
}
