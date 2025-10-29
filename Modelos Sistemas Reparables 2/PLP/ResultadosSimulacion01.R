# General Functions only freq &  gamma------------------------------------------------------------
Bias <- function(est1, est2, real_par) {
  est1 <- sweep(est1, 2, real_par) 
  est2 <- sweep(est2, 2, real_par)
  
  return(rbind(colMeans(est1, na.rm = TRUE),
               colMeans(est2, na.rm = TRUE)))
}


MSE <- function(est1, est2, real_par){
  est1 <- (sweep(est1, 2, real_par))^2
  est2 <- (sweep(est2, 2, real_par))^2
  
  return(rbind(colMeans(est1, na.rm = TRUE),
               colMeans(est2, na.rm = TRUE)))
}

sim_df  <- function(sim, param){ # simulación almacenamiento de datos y dataframe
  simul <- data.frame(
    Metodo = c("freq", "gamma"),
    cbind(Bias(sim$freq_est, sim$bayes_est_gamma, param),
          MSE(sim$freq_est, sim$bayes_est_gamma, param))
  )
  colnames(simul) <- c("Model", "Bias_beta", "Bias_mu_T", "Bias_theta", "MSE_beta", "MSE_mu_T", "MSE_theta" )
  return(simul)
}
#------------------------------
results <-  readRDS("resultadosk101_simulacion01.rds") 
# Extraer estimaciones frecuentistas
freq_est <- do.call(rbind, lapply(results, function(res) res$freq_est))
# Extraer estimaciones bayesianas con prior gamma
bayes_est_gamma <- do.call(rbind, lapply(results, function(res) res$bayes_est_gamma))

# Extraer cobertura bayesiana gamma
coverages_gamma <- do.call(rbind, lapply(results, function(res) res$coverage_gamma))
coverages_gamma_summary <- colMeans(coverages_gamma, na.rm = TRUE)



# Extraer ajustes (objetos de tipo 'stanfit', etc.)
ajuste_gamma <- lapply(results, function(res) res$ajuste_gamma)


sim <- list(freq_est          = freq_est,
              bayes_est_gamma = bayes_est_gamma,
              coverage_gamma  = coverages_gamma)

# Calcular Bias y MSE
param_true <- c(2, 30, 0.5)
im_df <- sim_df(sim,param_true)
im_df
#intevalos de credibilidad
coverages_gamma_summary

# General Functions ------------------------------------------------------------

Bias <- function(est1, est2, est3, real_par) {
  est1 <- sweep(est1, 2, real_par) #freq
  est2 <- sweep(est2, 2, real_par) # gamma
  est3 <- sweep(est3, 2, real_par) # corr
  
  return(rbind(colMeans(est1, na.rm = TRUE),
               colMeans(est2, na.rm = TRUE),
               colMeans(est3, na.rm = TRUE)))
}


MSE <- function(est1, est2, est3, real_par){
  est1 <- (sweep(est1, 2, real_par))^2
  est2 <- (sweep(est2, 2, real_par))^2
  est3 <- (sweep(est3, 2, real_par))^2
  
  return(rbind(colMeans(est1, na.rm = TRUE),
               colMeans(est2, na.rm = TRUE),
               colMeans(est3, na.rm = TRUE)))
}

sim_df  <- function(sim, param){ # simulación almacenamiento de datos y dataframe
  simul <- data.frame(
    Metodo = c("freq", "gamma", "correl"),
    cbind(Bias(sim$freq_est, sim$bayes_est_gamma, sim$bayes_est_corr, param),
          MSE(sim$freq_est, sim$bayes_est_gamma, sim$bayes_est_corr, param))
  )
  colnames(simul) <- c("Model", "Bias_beta", "Bias_mu_T","Bias_theta", "MSE_beta", "MSE_mu_T", "MSE_theta")
  return(simul)
}



sim_df_coverage   <- function(sim){
  simul           <- data.frame(
    Metodo = c("gamma", "correl"),
    rbind(sim$coverages_gamma, sim$coverages_corr)
  )
  colnames(simul) <- c("Model", "HPD_beta", "HPD_mu_T")
  return(simul)
}
#-------------------------------------
results <-  readRDS("resultadosk101_simulacion01.rds") 
# Extraer estimaciones frecuentistas
freq_est <- do.call(rbind, lapply(results, function(res) res$freq_est))
# Extraer estimaciones bayesianas con prior gamma
bayes_est_gamma <- do.call(rbind, lapply(results, function(res) res$bayes_est_gamma))
# Extraer estimaciones bayesianas con prior corr
bayes_est_corr <- do.call(rbind, lapply(results, function(res) res$bayes_est_corr))

# Extraer cobertura bayesiana gamma
coverages_gamma <- do.call(rbind, lapply(results, function(res) res$coverage_gamma))
coverages_gamma_summary <- colMeans(coverages_gamma, na.rm = TRUE)

# Extraer cobertura bayesiana corr
coverages_corr <- do.call(rbind, lapply(results, function(res) res$coverage_corr))
coverages_corr_summary <- colMeans(coverages_corr, na.rm = TRUE)



# Extraer ajustes (objetos de tipo 'stanfit', etc.)
ajuste_gamma <- lapply(results, function(res) res$ajuste_gamma)
ajuste_corr  <- lapply(results, function(res) res$ajuste_corr)

sim <- list(freq_est        = freq_est,
            bayes_est_gamma = bayes_est_gamma,
            coverage_gamma  = coverages_gamma,
            bayes_est_corr  = bayes_est_corr,
            coverage_corr   = coverages_corr)

# Calcular Bias y MSE
param_true <- c(2, 30, 0.5)
im_df <- sim_df(sim,param_true)
im_df
#intevalos de credibilidad
coverages_gamma_summary

