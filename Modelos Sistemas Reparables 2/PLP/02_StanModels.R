# Establecer un directorio temporal personalizado para evitar errores de espacio
Sys.setenv(TMPDIR = "/mnt/nfs/home/natalymr/tmp")
dir.create(Sys.getenv("TMPDIR"), recursive = TRUE, showWarnings = FALSE)

library(rstan)       # Para compilar y ejecutar modelos Stan
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores()-1)

# Stan Models ------------------------------------------------------------------
## gamma ----
# --- Modelo Stan con priors Gamma ---
modelo.stan.gamma <- "
functions {
  
  real lambda_PLP(real t, vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    
    
    
return mu_T * beta / time_trunc * pow(t, beta - 1);
}

  real Lambda_PLP(vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    
    
    
    return mu_T / time_trunc * pow(time_trunc, beta);
  }

  real log_like_PLP(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (i in 1:k) {
      for (j in start_idx[i]:end_idx[i]) {
        log_like += log(lambda_PLP(time_flat[j], param, time_trunc));
      }
    }
    
    log_like += - k * Lambda_PLP(param, time_trunc);
    
    return log_like;
  }
}

data {
  int<lower=1> k;                  // number of time vectors
  int<lower=1> len_times[k];       // lengths of each time vector
  vector[sum(len_times)] time_flat;  // flattened time vectors
  int<lower=1> start_idx[k];       // start indices for each subvector
  int<lower=1> end_idx[k];         // end indices for each subvector
  real<lower=0> time_trunc;        // time truncation point

  vector[2] beta0;
  vector[2] mu0;
}

parameters {
  real<lower=0> beta;               // beta parameter
  real<lower=0> mu_T;               // mu_T parameter
}

transformed parameters {
  vector[2] param = [beta, mu_T]';
}

model {
  // Prior distributions
  beta ~ gamma(beta0[1], beta0[2]);
  mu_T ~ gamma(mu0[1], mu0[2]);

  // Likelihood
  target += log_like_PLP(param,time_trunc, time_flat, start_idx, end_idx, k);
}
"

# --- Modelo Stan con parámetros correlacionados ---
modelo.stan.corr <- "
functions {
  

  real lambda_PLP(real t, vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    return mu_T * beta / time_trunc * pow(t, beta - 1);
  }

  
  real Lambda_PLP(vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    return mu_T / time_trunc * pow(time_trunc, beta);
  }

  real log_like_PLP(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;
    for (i in 1:k) {
      for (j in start_idx[i]:end_idx[i]) {
        log_like += log(lambda_PLP(time_flat[j], param, time_trunc));
      }
    }
    log_like += - k * Lambda_PLP(param, time_trunc);
    return log_like;
  }
}

data {
  int<lower=1> k;                  // number of time vectors
  int<lower=1> len_times[k];       // lengths of each time vector
  vector[sum(len_times)] time_flat;  // flattened time vectors
  int<lower=1> start_idx[k];       // start indices for each subvector
  int<lower=1> end_idx[k];         // end indices for each subvector
  real<lower=0> time_trunc;        // time truncation point

  real<lower=0> beta0;
  real<lower=0> mu0;
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_mu0;
  real<lower=0> eta;
}

parameters {
  vector[2] Z;
  cholesky_factor_corr[2] L;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_mu_T;
}

transformed parameters {
  vector[2] W;
  W = L * Z;
  
  real<lower=0> beta = exp(W[1] * sigma_beta + log(beta0));
  real<lower=0> mu_T = exp(W[2] * sigma_mu_T + log(mu0));

  vector[2] param = [beta, mu_T]';
}

model {
  // Prior distributions
  Z ~ normal(0, 1);
  L ~ lkj_corr_cholesky(eta);

  sigma_beta ~ uniform(0, sigma_beta0);
  sigma_mu_T ~ uniform(0, sigma_mu0);

  // Likelihood
  target += log_like_PLP(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"




