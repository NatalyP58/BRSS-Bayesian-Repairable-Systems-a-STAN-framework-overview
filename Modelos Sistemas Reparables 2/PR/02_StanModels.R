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

real indicadora(vector times, real t) {  #ojo los times deben estar ordenados para que funcione el break correctamente
  real max_val = 0;
  for (i in 1:num_elements(times)) {
    if (times[i] >= t) break;
    if (times[i] > max_val)
      max_val = times[i];
  }
  return max_val;
}

  real lambda_PR(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    
    real timess = indicadora(times, t);
    
    return mu_T * beta / (time_trunc) * pow(t - timess, beta-1);
  }

  real Lambda_PR(vector param, vector times, real time_trunc) {
     real beta = param[1];
     real mu_T = param[2];
  
     int N = num_elements(times);
     vector[N] delta_times;

  // Creamos las diferencias t_i - t_{i-1}  append_row(0, head(times, N - 1)); si (1,2,4) entonces (0,1,2)
     delta_times = times - append_row(0, head(times, N - 1));

  // Suma vectorizada de (t_i - t_{i-1})^β
     real sum_term = sum(pow(delta_times, beta));
 
  // Término de truncamiento: (T - t_n)^β
     real tail_term = pow(time_trunc - times[N], beta);

  return mu_T / time_trunc * (sum_term + tail_term);
}

  real log_like_PR(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {  
    int n_l = end_idx[l] - start_idx[l] + 1;  // Calculamos el tamaño de cada segmento 
  
     vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);  // Segmento de tiempos, extrae los k subcojuntos de tiempos
                                                                    
               // Suma de logaritmos de las tasas de fallos para cada evento dentro del segmento
             for (i in 1:n_l) {
                  log_like += log(lambda_PR(times_l[i], param, times_l, time_trunc));
  }

  // Resta la integral acumulada de la tasa de fallos para el segmento
  log_like -= Lambda_PR(param, times_l, time_trunc);
}
   return log_like; 
}
}
data {
  int<lower=1> k;                    // number of time vectors
  int<lower=1> len_times[k];         // lengths of each time vector
  vector[sum(len_times)] time_flat;  // flattened time vectors
  int<lower=1> start_idx[k];         // start indices for each subvector
  int<lower=1> end_idx[k];           // end indices for each subvector
  real<lower=0> time_trunc;          // time truncation point

  vector[2] beta0;
  vector[2] mu0;
}

parameters {
  real<lower=0> beta_raw;              // beta parameter
  real<lower=0> mu_T;                  // mu_T parameter
}

transformed parameters {
  real<lower=1> beta = beta_raw + 1;
  vector[2] param    = [beta, mu_T]';
}

model {
  // Prior distributions
  beta_raw ~ gamma(beta0[1], beta0[2]);
  mu_T     ~ gamma(mu0[1], mu0[2]);

  // Likelihood
  target += log_like_PR(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"

# --- Modelo Stan con parámetros correlacionados ---
modelo.stan.corr <- "
functions {

  real indicadora(vector times, real t) {  #ojo los times deben estar ordenados para que funcione el break correctamente
  real max_val = 0;
  for (i in 1:num_elements(times)) {
    if (times[i] >= t) break;
    if (times[i] > max_val)
      max_val = times[i];
  }
  return max_val;
}

  real lambda_PR(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    
    real timess = indicadora(times, t);
    
    return mu_T * beta / (time_trunc) * pow(t - timess, beta-1);
  }

  real Lambda_PR(vector param, vector times, real time_trunc) {
     real beta = param[1];
     real mu_T = param[2];
  
     int N = num_elements(times);
     vector[N] delta_times;

  // Creamos las diferencias t_i - t_{i-1}  append_row(0, head(times, N - 1)); si (1,2,4) entonces (0,1,2)
     delta_times = times - append_row(0, head(times, N - 1));

  // Suma vectorizada de (t_i - t_{i-1})^β
     real sum_term = sum(pow(delta_times, beta));
 
  // Término de truncamiento: (T - t_n)^β
     real tail_term = pow(time_trunc - times[N], beta);

  return mu_T / time_trunc * (sum_term + tail_term);
}

  real log_like_PR(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {  
    int n_l = end_idx[l] - start_idx[l] + 1;  // Calculamos el tamaño de cada segmento
    
     vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);  // Segmento de tiempos, extrae los k subcojuntosde tiempos
                                                                  
               // Suma de logaritmos de las tasas de fallos para cada evento dentro del segmento
             for (i in 1:n_l) {
                  log_like += log(lambda_PR(times_l[i], param, times_l, time_trunc));
  }

  // Resta la integral acumulada de la tasa de fallos para el segmento
  log_like -= Lambda_PR(param, times_l, time_trunc);
}
   return log_like; 
}
}

data {
  int<lower=1> k;                    // number of time vectors
  int<lower=1> len_times[k];         // lengths of each time vector
  vector[sum(len_times)] time_flat;  // flattened time vectors
  int<lower=1> start_idx[k];         // start indices for each subvector
  int<lower=1> end_idx[k];           // end indices for each subvector

  real<lower=0> time_trunc;          // time truncation point

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
  
  real<lower=1> beta = exp(W[1] * sigma_beta + log(beta0)) + 1;
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
  target += log_like_PR(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"





