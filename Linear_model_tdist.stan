data {
  // Patient PCR data
  int<lower=0> N;                          // Number of concatenated data points
  int<lower=0> n_id;                       // Number of individuals
  int<lower=1,upper=n_id> id[N];           // Patient identifier
  real<lower=0> obs_day[N];                // Time of sample
  real<lower=0,upper=40> delta_CT[N];      // 40-CT value
  real RNaseP[N];                          // Scaled RNaseP CT values (Mean 0, SD 1)
  int<lower=1> K_trt;                      // number of treatment arms
  int<lower=1,upper=K_trt+1> trt[N];       // trt index
  int<lower=1> K_plate;                    // Number of plates used for the PCR data
  int<lower=1,upper=K_plate> id_plate[N];  // plate identifier

  // PCR quality control data
  int<lower=1> N_control;                  // Number of control PCR samples data
  real control_density[N_control];
  real control_delta_CT[N_control];
  int<lower=1,upper=K_plate> cont_id_plate[N_control];  // plate identifier

  // priors
  real A0_prior_mean;
  real A0_prior_sd;

  real alpha_prior_mean;
  real alpha_prior_sd;

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

  real<lower=0> sigmaCT;

  // plate random effects for the intercept
  vector[K_plate] a_plate;
  real<lower=0> sigma_plate; // variance of plate effects

  // standard curve parameters
  real sc_intercept;
  real sc_slope;
  real<lower=0> sigma_sc;

  real<lower=0,upper=40> A0;          // population intercept
  real alpha;                         // population slope
  vector[2] theta_rand_id[n_id];      // individual random effects vector
  real gamma_rnasep;                  // Adjustment for RNaseP
  vector[K_trt] trt_effect;           // Estimates of the treatment effect

  real<lower=0> t_dof;                // Degrees of freedom for the t-distribution error model

}

transformed parameters {
  real pred_CT[N];
  vector[K_trt+1] trt_effect_prime;

  trt_effect_prime[1]=0; // no study drug arm
  for(i in 1:K_trt){
    trt_effect_prime[i+1]=trt_effect[i];
  }
  for(i in 1:N){
    pred_CT[i] =
    A0 + theta_rand_id[id[i]][1] + a_plate[id_plate[i]] + gamma_rnasep*RNaseP[i] +
    alpha*exp(trt_effect_prime[trt[i]])*exp(theta_rand_id[id[i]][2])*obs_day[i];
  }
}

model {
  // error model
  t_dof ~ exponential(1);

  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  a_plate ~ normal(0, sigma_plate);
  sigma_plate ~ exponential(1);

  // measurement error
  sigmaCT ~ normal(3,3) T[0,];

  // standard curve
  sc_intercept ~ normal(-3,1);
  sc_slope ~ normal(3.3, 1);
  sigma_sc ~ exponential(1);

  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));

  gamma_rnasep ~ normal(0,1);

  A0 ~ normal(A0_prior_mean,A0_prior_sd);
  alpha ~ normal(alpha_prior_mean,alpha_prior_sd);

  trt_effect ~ normal(0,sigma_trt_effect);

  // Main model specification:
  for(i in 1:N){
    if(delta_CT[i]>0){
      target += student_t_lpdf(delta_CT[i] | t_dof, pred_CT[i], sigmaCT);
    } else {
      target += student_t_lcdf(0 | t_dof, pred_CT[i], sigmaCT);
    }
  }

  // standard curve data;
  for(i in 1:N_control){
    control_delta_CT[i] ~ normal(sc_intercept+a_plate[cont_id_plate[i]] + sc_slope*control_density[i], sigma_sc);
  }
}

generated quantities {
  real pred_CT_mean[N]; // without the RNaseP correction for plotting
  for(i in 1:N){
    pred_CT_mean[i] =
    A0 + theta_rand_id[id[i]][1] + a_plate[id_plate[i]] +
    alpha*exp(trt_effect_prime[trt[i]])*exp(theta_rand_id[id[i]][2])*obs_day[i];
  }
}
