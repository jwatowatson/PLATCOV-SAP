functions{
  real g_inverse(vector beta, real absorb_obs){//beta is log scale
  real log_x;
  log_x = -(1/exp(beta[4])) * log(exp(beta[2])/(absorb_obs - exp(beta[1])) - 1) + beta[3];
  return(exp(log_x));
  }
}

data {
  int<lower=0> N_controls;                                   // number of data points in the controls
  int<lower=0> K_plates;                                     // number of plates used
  vector<lower=0>[N_controls] conc_controls;                 // concentration values in controls
  vector<lower=0>[N_controls] absorb_controls;               // observed color intensity in controls
  int<lower=1,upper=K_plates> ind_plate_controls[N_controls];// index of plate for each control sample
  
  int N_sample;                                     // number of unique unknowns
  int<lower=0> N_obs;                               // number of observations
  vector<lower=0>[N_obs] dilution_factor;           // dilution for that sample
  int<lower=1,upper=N_sample> ind_sample[N_obs];    // index of sample giving unique sample ID
  vector<lower=0>[N_obs] absorb_obs;                // observed absorbance value
  int<lower=1,upper=K_plates> ind_plate_obs[N_obs]; // index of plate for each measurement
}

transformed data {
  int K_z = 4;
  vector[K_z] my_zeros;
  for(i in 1:K_z) my_zeros[i] = 0;
}

parameters {
  
  cholesky_factor_corr[K_z] L_Omega;     // correlation matrix for control RE across plates
  vector<lower=0>[K_z] sigmasq_u;        // variance of random effects
  
  // standard curve and error terms
  vector[4] beta;                        // standard curve parameters
  real<lower=0> alpha;                   // heterosckedasticty parameter
  real<lower=0> sigma;                   // measurement error sd term
  
  // Random effects
  vector[K_z] theta_rand[K_plates];       // individual random effects vector
  
  // unknown concentrations
  vector[N_sample] log_x_init;
}

transformed parameters {
  vector<lower=0>[N_controls] g_controls;
  vector<lower=0>[N_controls] tau_controls;
  vector<lower=0>[N_obs] g_unk;
  vector<lower=0>[N_obs] tau_unk;
  vector[4] beta_plate[K_plates];
  
  for(k in 1:K_plates){
    beta_plate[k]=beta+theta_rand[k];
  }
  
  // compute expected absorbance value for the controls
  for (i in 1:N_controls) {
    // expected absorbance value
    g_controls[i] = exp(beta_plate[ind_plate_controls[i]][1]) + 
    exp(beta_plate[ind_plate_controls[i]][2]) / 
    (1 + (conc_controls[i] / exp(beta_plate[ind_plate_controls[i]][3])) ^ (- exp(beta_plate[ind_plate_controls[i]][4])));
    // variance
    tau_controls[i] = (g_controls[i] ^ alpha) * (sigma);
  }
  
  // compute expected absorbance value for the unknowns
  for (i in 1:N_obs){
    // expected absorbance value
    g_unk[i] = exp(beta_plate[ind_plate_obs[i]][1]) + 
    exp(beta_plate[ind_plate_obs[i]][2]) / 
    (1 + ((exp(log_x_init[ind_sample[i]]) * dilution_factor[i]) / exp(beta_plate[ind_plate_obs[i]][3])) ^ (- exp(beta_plate[ind_plate_obs[i]][4])));
    // variance
    tau_unk[i] = (g_unk[i] ^ alpha) * (sigma);
  }
}

model {
  // Priors
  alpha ~ normal(1, .5);
  beta[1] ~ normal(-2, 1);
  beta[2] ~ normal(1, 1);
  beta[3] ~ normal(-0.5, 1);
  beta[4] ~ normal(.5, 1);
  
  sigmasq_u[1] ~ normal(.5, .5) T[0,]; 
  sigmasq_u[2] ~ normal(.05,.025) T[0,]; 
  sigmasq_u[3] ~ normal(.15, .1) T[0,]; 
  sigmasq_u[4] ~ normal(0.5, 0.25) T[0,]; 
  
  log_x_init ~ normal(3,2);
  
  sigma ~ normal(.1, .05);
  
  L_Omega ~ lkj_corr_cholesky(3); // covariance matrix - random effects across plates
  
  // individual random effects
  for(i in 1:K_plates) theta_rand[i] ~ multi_normal_cholesky(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));
  
  // Likelihood for the control data
  absorb_controls ~ normal(g_controls, tau_controls);
  
  // Likelihood for the unknown samples data
  absorb_obs ~ student_t(10, g_unk, tau_unk);
}
