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
  real beta[3];
}

transformed parameters {
  # here we can write functions of parameters we want to keep track of.
  matrix[N, n_years] mu;
  # first construct the linear predictor
  # for each age
  for (i in 1:N){
    # for each year
    for (t in 1:n_years){
      mu[i,t] = inv_logit(beta[1] + age_standardised[i] * beta[2] + 
                           years_standardised[t] * beta[3]);
      
    }
  }
}

model {
  # here we write the model

  # priors -------------------------------------------------------------
  
  beta ~ normal(0, 10);
  
  # likelihood ---------------------------------------------------------
  
  for (i in 1:N){
    for (t in 1:n_years){
      # ignore if there is no exposure. 
      if (expos[i,t] != 0){
        # add log_poisson density p(D | mu*R) to the log-posterior
        deaths[i,t] ~ poisson(mu[i,t] * expos[i,t]);
        
        # the (commented) statement below does almost exactly the same thing
        # target += poisson_lpdf(deaths[i,t] | mu[i,t] * expos[i,t])
        # target is the name of stan variable in which lp(theta | y) is accumulated.
      }
    }
  }
  
}


