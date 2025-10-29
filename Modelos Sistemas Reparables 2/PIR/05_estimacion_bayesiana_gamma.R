# Librerías necesarias
library(rstan)
library(coda)

rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

bayesian_estimation_gamma <- function(time, time_trunc,  beta_r = NULL, mu_T_r = NULL, theta_r = NULL,
                                      tau, compiled_model = NULL, simulation = TRUE, seed = NULL) {
  # Si no se proporciona un modelo compilado, usar el global
  if (is.null(compiled_model)) {
    if (!exists("modelo.compilado.gamma")) {
      stop("modelo.compilado.gamma no encontrado en el entorno global")
    }
    compiled_model <- modelo.compilado.gamma
  }
  
  # Validación: si es simulación, deben darse los valores reales
  if (simulation && (is.null(beta_r) || is.null(mu_T_r) || is.null(theta_r))) {
    stop("Si simulation = TRUE debes proporcionar los valores reales: beta_r, mu_T_r y theta_r.")
  }
  
  
  time_flat <- unlist(time)
  len_times <- sapply(time, length)
  start_idx <- cumsum(c(1, len_times[-length(len_times)]))
  end_idx   <- cumsum(len_times)
  
  datos.stan <- list(
    k          = length(time),                   
    N          = length(tau),
    len_times  = array(len_times, dim = length(len_times)),                       
    time_flat  = time_flat,                       
    start_idx  = array(start_idx, dim = length(start_idx)),                       
    end_idx    = array(end_idx, dim = length(end_idx)),                         
    time_trunc = time_trunc, 
    tau        = array(tau, dim = length(tau)),
    beta0      = c(0.16, 0.04), 
    mu0        = c(2, 0.05)   
  )
  
  ajuste <- sampling(compiled_model, data = datos.stan,
                     chains = 1, iter = 3000, warmup = 1000,
                     seed = if (is.null(seed)) 123L else as.integer(seed),
                     control = list(adapt_delta = 0.95), refresh = 0)
  
  zscores <- abs(geweke.diag(As.mcmc.list(ajuste))[[1]]$z[1:3])
  
  if (all(!is.na(zscores)) && any(zscores > 1.96)) {
    return(list(success = FALSE, estimate = NULL, coverage = c(0, 0, 0), fit = ajuste))
  }
  
  samples  <- rstan::extract(ajuste)
  beta_est <- samples$beta
  mu_T_est <- samples$mu_T
  theta_est<- samples$theta
  
  estimate <- c(median(beta_est), median(mu_T_est), median(theta_est))
  
  if (simulation) {
    beta_hpd  <- HPDinterval(as.mcmc(matrix(beta_est, ncol = 1)), prob = 0.95)
    mu_T_hpd  <- HPDinterval(as.mcmc(matrix(mu_T_est, ncol = 1)), prob = 0.95)
    theta_hpd <- HPDinterval(as.mcmc(matrix(theta_est, ncol = 1)), prob = 0.95)
    
    coverage <- c(
      as.integer(beta_hpd[1] <= beta_r & beta_r <= beta_hpd[2]),
      as.integer(mu_T_hpd[1] <= mu_T_r & mu_T_r <= mu_T_hpd[2]),
      as.integer(theta_hpd[1] <= theta_r & theta_r <= theta_hpd[2])
    )
  } else {
    coverage <- NULL
  }
  
  return(list(success = TRUE, estimate = estimate, coverage = coverage, fit = ajuste))
}
