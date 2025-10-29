# Librerías necesarias solo para esta parte
library(pracma)

gen_failures_PIR <- function(k, param, time_trunc, tau) {
  

  h_theta <- function(theta, beta) {
    if ((beta < 1 & theta >= 0) || (beta >= 1 & theta < 0)) {
      return(log((1 + theta) / (1 - theta)) * 1e4)
    } else {
      return(theta)
    }
  }
  
  h_theta <- Vectorize(h_theta, "theta")
  
  
  g_theta <- function(theta) {return(theta)}
  
  
  g_theta <- Vectorize(g_theta, "theta")
  
  signo_beta <- function(beta) {
    ifelse(beta >= 1, 1, -1)
  }
  

  k_t <- function(t, tau, time_trunc) {
    findInterval(t, c(0, tau, time_trunc), left.open = TRUE)
  }
  
  
  lambda_PIR_acum <- function(t, param, tau, time_trunc) {
    if (t <= 0) return(0)
    
    beta  <- param[1]
    mu_T  <- param[2]
    theta <- param[3]
    
    kt    <- k_t(t, tau, time_trunc)
    
    tau0  <- c(0, tau)[1:(kt - 1)]
    tau1  <- c(tau, time_trunc)[1:(kt - 1)]
    tau_k <- ifelse(kt - 1 > 0, c(tau, time_trunc)[kt - 1], 0)
    
    h_val  <- h_theta(theta, beta)
    g_val  <- g_theta(theta)
    s_beta <- signo_beta(beta)
    
    lambda1 <- if (kt - 1 > 0) {
      num   <- (tau1 - s_beta * h_val * tau0)^beta - (tau0 - s_beta * h_val * tau0)^beta
      mu_T * num / (time_trunc + g_val * tau0)
    } else 0
    
    lambda1[is.na(lambda1)] <- 0
    
    num2    <- (t - s_beta * h_val * tau_k)^beta - (tau_k - s_beta * h_val * tau_k)^beta
    lambda2 <- mu_T * num2 / (time_trunc + g_val * tau_k)
    lambda2[is.na(lambda2)] <- 0
    
    return(sum(lambda1) + lambda2)
  }
  lambda_PIR_acum <- Vectorize(lambda_PIR_acum, "t")
  
 
  objetivo_incremento <- function(x, param, t, tau, time_trunc, u) {
    abs(lambda_PIR_acum(t + x, param, tau, time_trunc) -
          lambda_PIR_acum(t, param, tau, time_trunc) + log(1 - u))
  }
  
  
  simular_una_trayectoria <- function(param, tau, time_trunc) {
    t_actual <- 0
    fallas   <- numeric(0)
    
    repeat {
      u        <- runif(1)
      objetivo <- -log(1 - u)
      
      max_Lambda <- lambda_PIR_acum(time_trunc, param, tau, time_trunc) -
        lambda_PIR_acum(t_actual, param, tau, time_trunc)
      
      if (objetivo > max_Lambda) break
      
      resultado <- optimize(
        objetivo_incremento,
        interval   = c(0, time_trunc - t_actual),
        param      = param,
        t          = t_actual,
        tau        = tau,
        time_trunc = time_trunc,
        u          = u
      )
      
      incremento  <- resultado$minimum
      t_siguiente <- t_actual + incremento
      
      if (t_siguiente > time_trunc) break
      
      fallas   <- c(fallas, t_siguiente)
      t_actual <- t_siguiente
    }
    
    return(fallas)
  }
  
  # Ejecutar simulación para k trayectorias
  replicate(k, simular_una_trayectoria(param, tau, time_trunc), simplify = FALSE)
}




indicadora <- function(vector, t) {
  vals     <- vector[vector < t]
  if (length(vals) == 0) return(0)  
  max(vals)  
}

h <- function(theta, beta){
  if ((beta < 1 & theta >= 0) | (beta >= 1 & theta < 0)){
    return(as.numeric(log((1 + theta) / (1 - theta)) * 10^4))}
  else {
    return(as.numeric(theta))
  }
}

h    <- Vectorize(h, vectorize.args = "theta")

g    <- function(theta) {return(theta)}

g    <- Vectorize(g, "theta")

sign <- function(beta){ifelse(beta >= 1, 1, -1)}

lambda_PIR <- function(t, param, tau, time_trunc){
  if (t <= 0) return(0)
  beta  <- param[1]
  mu_T  <- param[2]
  theta <- param[3]
  
  taus  <- indicadora(tau, t) # taus= t_i-1
  h_theta <- h(theta, beta)
  signo   <- sign(beta)
  g_theta <- g(theta)
  
  mu_T * beta / (time_trunc + g_theta * taus) * (t - signo*h_theta*taus)^(beta-1)
}

lambda_PIR <- Vectorize(lambda_PIR, vectorize.args = "t")

Lambda_PIR <- function(param, tau, time_trunc){  # aqui tau={t_1,....,tn} t_o=0 y t_n+1=Time_trunc=T
  if (length(tau) == 0 || max(tau) <= 0) return(0)
  beta  <- param[1]
  mu_T  <- param[2]
  theta <- param[3]
  
  tau0 <- c(0, tau[-length(tau)])
  tau1 <- tau
  taup <- tau[length(tau)]
  
  h_theta <- h(theta, beta)
  g_theta <- g(theta)
  signo   <- sign(beta)
  
 
  sum(mu_T / (time_trunc + g_theta * tau0) * ((tau1 - signo*h_theta* tau0)^beta-(tau0 - signo*h_theta* tau0)^beta))   +   mu_T * ((time_trunc -  signo*h_theta*taup)^beta-(taup -  signo*h_theta*taup)^beta) / (time_trunc + g_theta * taup)
}

log_like_PIR <- function(param, times, tau, time_trunc) {
  
  
  log_like_single <- function(times_i) {
    if (length(times_i) == 0) return(0) 
    log_lambda_sum <- sum(log(lambda_PIR(times_i, param, tau, time_trunc)))
    
    
    lambda_total   <- Lambda_PIR(param, tau, time_trunc)
    log_lambda_sum - lambda_total
  }
  
  sum(sapply(times, log_like_single)) 
  }

