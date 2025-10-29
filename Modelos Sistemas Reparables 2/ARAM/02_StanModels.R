# Establecer un directorio temporal personalizado para evitar errores de espacio
Sys.setenv(TMPDIR = "/mnt/nfs/home/natalymr/tmp")
dir.create(Sys.getenv("TMPDIR"), recursive = TRUE, showWarnings = FALSE)

library(rstan)       # Para compilar y ejecutar modelos Stan
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

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

real h(real theta, real beta) {
    return ((beta < 1 && theta >= 0) || (beta >= 1 && theta < 0)) ? log((1 + theta) / (1 - theta)) * 1e4 : theta;
  }

  int sign(real beta) {
    return (beta >= 1) ? 1 : -1;
  }

  real lambda_ARAM(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];
    
    real timess  = indicadora(times, t);
    real h_theta = h(theta, beta);
    real signo   = sign(beta);
    
    return (mu_T * beta / (time_trunc) ) * pow(t - signo * h_theta * timess, beta-1);
  }

  real Lambda_ARAM(vector param, vector times, real time_trunc) {
     real beta    = param[1];
     real mu_T    = param[2];
     real theta   = param[3];
     real h_theta = h(theta, beta);
     real signo   = sign(beta);

  int N = num_elements(times);
  vector[N] delta_times_1;
  vector[N] delta_times_2;
  vector[N] timesPrev;
  real      sum_term;
  real      tail_term;
  
  timesPrev     = append_row(0, head(times, N - 1));
  delta_times_1 = (times - signo * h_theta * timesPrev);
  delta_times_2 = (timesPrev - signo * h_theta *  timesPrev);

  
  sum_term  = sum(pow(delta_times_1, beta)./ (time_trunc)) - sum(pow(delta_times_2, beta)./ (time_trunc));
  tail_term = 1/(time_trunc) * (pow(time_trunc - signo * h_theta *  times[N], beta) - pow(times[N] - signo * h_theta *  times[N], beta));

  return mu_T * (sum_term + tail_term);
}
  real log_like_ARAM(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {  
    int n_l = end_idx[l] - start_idx[l] + 1;  // Calculamos el tamaño de cada segmento 
    
     vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);  // Segmento de tiempos, extrae los k subcojuntosde tiempso
               
               // Suma de logaritmos de las tasas de fallos para cada evento dentro del segmento
             for (i in 1:n_l) {
                  log_like += log(lambda_ARAM(times_l[i], param, times_l, time_trunc));
  }

  // Resta la integral acumulada de la tasa de fallos para el segmento
  log_like -= Lambda_ARAM(param, times_l, time_trunc);
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
  real<lower=0> beta;                  // beta parameter
  real<lower=0> mu_T;                  // mu_T parameter
  real<lower=-1, upper=1> theta;       // theta parameter
}

transformed parameters {
  vector[3] param    = [beta, mu_T, theta]';
}

model {
  // Prior distributions
  beta     ~ gamma(beta0[1], beta0[2]);
  mu_T     ~ gamma(mu0[1], mu0[2]);
  theta    ~ uniform(-1,1);

  // Likelihood
  target += log_like_ARAM(param, time_trunc, time_flat, start_idx, end_idx, k);
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

real h(real theta, real beta) {
    return ((beta < 1 && theta >= 0) || (beta >= 1 && theta < 0)) ? log((1 + theta) / (1 - theta)) * 1e4 : theta;
  }

  int sign(real beta) {
    return (beta >= 1) ? 1 : -1;
  }

  real lambda_ARAM(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];
    
    real timess  = indicadora(times, t);
    real h_theta = h(theta, beta);
    real signo   = sign(beta);
    
    return mu_T * beta / (time_trunc) * pow(t - signo * h_theta * timess, beta-1);
  }

  real Lambda_ARAM(vector param, vector times, real time_trunc) {
     real beta    = param[1];
     real mu_T    = param[2];
     real theta   = param[3];
     real h_theta = h(theta, beta);
     real signo   = sign(beta);

  int N = num_elements(times);
  vector[N] delta_times_1;
  vector[N] delta_times_2;
  vector[N] timesPrev;
  real      sum_term;
  real      tail_term;
  
  timesPrev     = append_row(0, head(times, N - 1));
  delta_times_1 = (times - signo * h_theta * timesPrev);
  delta_times_2 = (timesPrev - signo * h_theta *  timesPrev);

  
  sum_term  = sum(pow(delta_times_1, beta)./ (time_trunc)) - sum(pow(delta_times_2, beta)./ (time_trunc));
  tail_term = 1/(time_trunc) * (pow(time_trunc - signo * h_theta *  times[N], beta) - pow(times[N] - signo * h_theta *  times[N], beta));

  return mu_T * (sum_term + tail_term);
}

  real log_like_ARAM(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {  
    int n_l = end_idx[l] - start_idx[l] + 1;  // Calculamos el tamaño de cada segmento 
  
     vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);  // Segmento de tiempos, extrae los k subcojuntosde tiempso
     
               // Suma de logaritmos de las tasas de fallos para cada evento dentro del segmento
             for (i in 1:n_l) {
                  log_like += log(lambda_ARAM(times_l[i], param, times_l, time_trunc));
  }

  // Resta la integral acumulada de la tasa de fallos para el segmento
  log_like -= Lambda_ARAM(param, times_l, time_trunc);
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
  vector[3] Z;
  cholesky_factor_corr[3] L;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_mu_T;
}

transformed parameters {
  vector[3] W;
  W = L * Z;
  
  real<lower=1> beta            = exp(W[1] * sigma_beta + log(beta0));
  real<lower=0> mu_T            = exp(W[2] * sigma_mu_T + log(mu0));
  real<lower=-1, upper=1> theta = tanh(W[3]);

  vector[3] param = [beta, mu_T, theta]';
}

model {
  // Prior distributions
  Z ~ normal(0, 1);
  L ~ lkj_corr_cholesky(eta);

  sigma_beta ~ uniform(0, sigma_beta0);
  sigma_mu_T ~ uniform(0, sigma_mu0);

  // Likelihood
  target += log_like_ARAM(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"





