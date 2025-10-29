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
real indicadora(vector tau, real t) {
  real max_val = 0;
  for (i in 1:num_elements(tau)) {
    if (tau[i] >= t) break;
    if (tau[i] > max_val)
      max_val = tau[i];
  }
  return max_val;
}

real h(real theta, real beta) {
  return ((beta < 1 && theta >= 0) || (beta >= 1 && theta < 0)) ? log((1 + theta) / (1 - theta)) * 1e4 : theta;
}

real g(real theta) {
 return  theta;
}

int sign(real beta) {
  return (beta >= 1) ? 1 : -1;
}

real lambda_PIR(real t, vector param, vector tau, real time_trunc) {
  real beta  = param[1];
  real mu_T  = param[2];
  real theta = param[3];

  real taus    = indicadora(tau, t);
  real h_theta = h(theta, beta);
  real signo   = sign(beta);
  real g_theta = g(theta);

  return mu_T * beta / (time_trunc + g_theta * taus) * pow(t - signo * h_theta * taus, beta - 1);
}

real Lambda_PIR(vector param, vector tau, real time_trunc) {
  real beta    = param[1];
  real mu_T    = param[2];
  real theta   = param[3];
  real h_theta = h(theta, beta);
  real signo   = sign(beta);
  real g_theta = g(theta);

  int N = num_elements(tau);
  vector[N] tauPrev = append_row(0, head(tau, N - 1));
  vector[N] num1    = tau - signo * h_theta * tauPrev;
  vector[N] num2    = tauPrev - signo * h_theta * tauPrev;

  real sum_term  = sum( pow(num1, beta) ./ (time_trunc + g_theta * tauPrev) )
                 - sum( pow(num2, beta) ./ (time_trunc + g_theta * tauPrev) );

  real tail_term = ( pow(time_trunc - signo * h_theta * tau[N], beta)
                   - pow(tau[N]      - signo * h_theta * tau[N], beta) )
                   / (time_trunc + g_theta * tau[N]);

  return mu_T * (sum_term + tail_term);
}


real log_like_PIR(vector param, real time_trunc, vector time_flat, vector tau, int[] start_idx, int[] end_idx, int k) {
  real log_like = 0;
  for (l in 1:k) {
    int n_l = end_idx[l] - start_idx[l] + 1;
    vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);
    for (i in 1:n_l) {
      log_like += log(lambda_PIR(times_l[i], param, tau, time_trunc));
    }
    log_like -= Lambda_PIR(param, tau, time_trunc);
  }
  return log_like;
}
}

data {
  int<lower=1> k;
  int<lower=1> N;
  int<lower=1> len_times[k];
  vector[sum(len_times)] time_flat;
  int<lower=1> start_idx[k];
  int<lower=1> end_idx[k];
  real<lower=0> time_trunc;
  vector<lower=0>[N] tau;
  vector[2] beta0;
  vector[2] mu0;
}

parameters {
  real<lower=0> beta;
  real<lower=0> mu_T;
  real<lower=-1, upper=1> theta;
}

transformed parameters {
  vector[3] param = [beta, mu_T, theta]';
}

model {
  beta ~ gamma(beta0[1], beta0[2]);
  mu_T ~ gamma(mu0[1], mu0[2]);
  theta ~ uniform(-1, 1);
  target += log_like_PIR(param, time_trunc, time_flat, tau, start_idx, end_idx, k);
}
"

# --- Modelo Stan con parámetros correlacionados ---
modelo.stan.corr <- "
functions {
real indicadora(vector tau, real t) {
  real max_val = 0;
  for (i in 1:num_elements(tau)) {
    if (tau[i] >= t) break;
    if (tau[i] > max_val)
      max_val = tau[i];
  }
  return max_val;
}

real h(real theta, real beta) {
  return ((beta < 1 && theta >= 0) || (beta >= 1 && theta < 0)) ? log((1 + theta) / (1 - theta)) * 1e4 : theta;
}

real g(real theta) {
 return  theta;
}

int sign(real beta) {
  return (beta >= 1) ? 1 : -1;
}

real lambda_PIR(real t, vector param, vector tau, real time_trunc) {
  real beta  = param[1];
  real mu_T  = param[2];
  real theta = param[3];

  real taus    = indicadora(tau, t);
  real h_theta = h(theta, beta);
  real signo   = sign(beta);
  real g_theta = g(theta);

  return mu_T * beta / (time_trunc + g_theta * taus) * pow(t - signo * h_theta * taus, beta - 1);
}

real Lambda_PIR(vector param, vector tau, real time_trunc) {
  real beta    = param[1];
  real mu_T    = param[2];
  real theta   = param[3];
  real h_theta = h(theta, beta);
  real signo   = sign(beta);
  real g_theta = g(theta);

  int N = num_elements(tau);
  vector[N] tauPrev = append_row(0, head(tau, N - 1));
  vector[N] num1    = tau - signo * h_theta * tauPrev;
  vector[N] num2    = tauPrev - signo * h_theta * tauPrev;

  real sum_term  = sum( pow(num1, beta) ./ (time_trunc + g_theta * tauPrev) )
                 - sum( pow(num2, beta) ./ (time_trunc + g_theta * tauPrev) );

  real tail_term = ( pow(time_trunc - signo * h_theta * tau[N], beta)
                   - pow(tau[N]      - signo * h_theta * tau[N], beta) )
                   / (time_trunc + g_theta * tau[N]);

  return mu_T * (sum_term + tail_term);
}


real log_like_PIR(vector param, real time_trunc, vector time_flat, vector tau, int[] start_idx, int[] end_idx, int k) {
  real log_like = 0;
  for (l in 1:k) {
    int n_l = end_idx[l] - start_idx[l] + 1;
    vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);
    for (i in 1:n_l) {
      log_like += log(lambda_PIR(times_l[i], param, tau, time_trunc));
    }
    log_like -= Lambda_PIR(param, tau, time_trunc);
  }
  return log_like;
}
}

data {
  int<lower=1> k;
  int<lower=1> N;
  int<lower=1> len_times[k];
  vector[sum(len_times)] time_flat;
  int<lower=1> start_idx[k];
  int<lower=1> end_idx[k];
  real<lower=0> time_trunc;
  vector<lower=0>[N] tau;
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

  real<lower=0> beta = exp(W[1] * sigma_beta + log(beta0));
  real<lower=0> mu_T = exp(W[2] * sigma_mu_T + log(mu0));
  real<lower=-1, upper=1> theta = tanh(W[3]);

  vector[3] param = [beta, mu_T, theta]';
}

model {
  Z ~ normal(0, 1);
  L ~ lkj_corr_cholesky(eta);

  sigma_beta ~ uniform(0, sigma_beta0);
  sigma_mu_T ~ uniform(0, sigma_mu0);

  target += log_like_PIR(param, time_trunc, time_flat, tau, start_idx, end_idx, k);
}
"





