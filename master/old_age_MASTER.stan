data {
  # here we declare what data we are passing in from outside stan (e.g. from R)
  int N; # ages
  int n_years; # years of data
  
  # matrix of deaths
  int<lower=0> deaths[N, n_years];
  
  # matrix of exposures
  matrix<lower=0>[N, n_years] expos;
  
  # we also pass in stardardised year and age variables.
  # we could also easily do this within stan.
  vector[N] age_standardised;
  vector[n_years] years_standardised;
}

parameters {
  # here we declare the parameter we wish to sample
  real beta[4];
  real<lower=0> log_dispersion;
  row_vector[n_years-2] kk_raw;
  real<lower=0> sigma_k;
  
}

transformed parameters {
  # here we can write functions of parameters we want to keep track of.
  matrix[N, n_years] mu;
  row_vector[n_years] kk;
  real<lower=0> dispersion;
  
  dispersion = exp(log_dispersion);
  # construct k so that it sums to zero and has zero growth.
  kk[2] = 0;
  for (t in 3:n_years){
    kk[2] = kk[2] - kk_raw[t-2]*(t-1);
    kk[t] = kk_raw[t-2];
  }
  kk[1] = -sum(kk[2:n_years]);
  
    
  # construct the linear predictor
  # for each age
  for (i in 1:N){
    # for each year
    for (t in 1:n_years){
      mu[i,t] = inv_logit(beta[1] + age_standardised[i] * beta[2] + 
                           years_standardised[t] * beta[3] +
                           years_standardised[t] * age_standardised[i] * beta[4] + 
                           kk[t]);
      
    }
  }
}

model {
  # here we write the model

  # priors -------------------------------------------------------------
  
  beta ~ normal(0, 10);
  sigma_k ~ normal(0, 10);
  kk_raw ~ normal(0, sigma_k);
 
  
  # likelihood ---------------------------------------------------------
  
  for (i in 1:N){
    for (t in 1:n_years){
      # ignore if there is no exposure. 
      if (expos[i,t] != 0){
        deaths[i,t] ~ neg_binomial_2(mu[i, t] * expos[i, t], dispersion);
      }
    }
  }
}

generated quantities{
  # here we include quantities of interest we want to calculate and save 
  # but which do not affect the log-posterior.
  matrix[N,n_years] simulated_deaths;
  for (i in 1:N){
    for (t in 1:n_years){
       simulated_deaths[i,t] = poisson_rng(mu[i,t]*expos[i,t]);
    }
  }
}
