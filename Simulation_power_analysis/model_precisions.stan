data {
  int<lower=0> N_sim;
  vector[N_sim] mu_hat;
  vector[N_sim] sigma_hat;
  vector[N_sim] mu;
  int<lower=0> condition_i[N_sim];
  int<lower=0> N_conditions;
}

parameters {
  real<lower=0> sigma_mu[N_conditions];
  real<lower=0> sigma_sigma[N_conditions];
  real mu_sigma[N_conditions];

}

model {
  
  for(i in 1:N_sim){
  mu_hat[i] ~ normal(mu[i], sigma_mu[condition_i[i]]);
  sigma_hat[i]  ~ normal(mu_sigma[condition_i[i]], sigma_sigma[condition_i[i]]);
  }
  
  
  //sigma_mu ~ exponential(1);
  //mu_sigma ~ exponential(1);
  //sigma_sigma ~ exponential(1);

}

