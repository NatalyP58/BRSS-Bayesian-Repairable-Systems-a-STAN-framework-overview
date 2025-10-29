# Librerías necesarias solo para esta parte
library(pracma)

gen_failures_ARAM <- function(k, param, time_trunc) {
  

  h <- function(theta, beta){
    if ((beta < 1 & theta >= 0) | (beta >= 1 & theta < 0)){
      return(as.numeric(log((1 + theta) / (1 - theta)) * 10^4))}
    else {
      return(as.numeric(theta))
    }
  }
  h <- Vectorize(h, vectorize.args = "theta")
  
  sign <- function(beta){ifelse(beta >= 1, 1, -1)}
  
  
  # Función interna k_t
  k_t <- function(t, times, time_trunc) {
    findInterval(t, c(0, times, time_trunc), left.open = TRUE) 
  }
  
  
  
  lambda_ARAM_acum <- function(t, param, times, time_trunc) {
    if (t <= 0) return(0)
    beta  <- param[1]
    mu_T  <- param[2]
    theta <- param[3]
    
    kt      <- k_t(t, times, time_trunc)
    times0  <- c(0, times)[1:(kt - 1)]
    times1  <- c(times, time_trunc)[1:(kt - 1)]
    times_k <- ifelse(kt - 1 > 0, c(times, time_trunc)[kt - 1], 0)
    
    h_theta <- h(theta, beta)
    signo   <- sign(beta)
    
    if (kt - 1 > 0) {
      
      lambda1 <- (mu_T  / time_trunc) * ((times1 - signo*h_theta*times0)^beta-(times0 - signo*h_theta*times0)^beta) 
      lambda1[is.na(lambda1)] <- 0
    } else {
      lambda1 <- 0
    }
    
    lambda2 <- (mu_T  / time_trunc ) * ((t - signo*h_theta*times_k)^beta-(times_k - signo*h_theta*times_k)^beta)
    lambda2[is.na(lambda2)] <- 0
    
    return(sum(lambda1) + lambda2)
  }
  
  lambda_ARAM_acum <- Vectorize(lambda_ARAM_acum, vectorize.args = "t")
  
 
  objetivo_incremento <- function(x, param, t, times, time_trunc, u) {
                         abs(lambda_ARAM_acum(t + x, param, times, time_trunc) -
                         lambda_ARAM_acum(t, param, times, time_trunc) + log(1 - u))
  }
  
  
  simular_una_trayectoria <- function(param, times, time_trunc) {
    t_actual <- 0
    times   <- numeric(0)
    
    repeat {
      u        <- runif(1)
      objetivo <- -log(1 - u)
      
      max_Lambda <- lambda_ARAM_acum(time_trunc, param, times, time_trunc) -
                    lambda_ARAM_acum(t_actual, param, times, time_trunc)
      
      if (objetivo > max_Lambda) break
      
      resultado <- optimize(
        objetivo_incremento,
        interval   = c(0, time_trunc - t_actual),
        param      = param,
        t          = t_actual,
        times      = times,
        time_trunc = time_trunc,
        u          = u
      )
      
      incremento  <- resultado$minimum
      t_siguiente <- t_actual + incremento
      
      if (t_siguiente > time_trunc) break
      
      times    <- c(times, t_siguiente)
      t_actual <- t_siguiente
    }
    
    return(times)
  }
  
  # Ejecutar simulación para k trayectorias
  replicate(k, simular_una_trayectoria(param, times, time_trunc), simplify = FALSE)
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

sign <- function(beta){ifelse(beta >= 1, 1, -1)}

lambda_ARAM <- function(t, param, times, time_trunc){
  if (t <= 0) return(0)
  beta  <- param[1]
  mu_T  <- param[2]
  theta <- param[3]
  
  timess  <- indicadora(times, t) 
  h_theta <- h(theta, beta)
  signo   <- sign(beta)
  
  (mu_T * beta / (time_trunc) )* (t - signo*h_theta*timess)^(beta-1)
}

lambda_ARAM <- Vectorize(lambda_ARAM, vectorize.args = "t")

Lambda_ARAM <- function(param, times, time_trunc){  
  if (length(times) == 0 || max(times) <= 0) return(0)
  beta  <- param[1]
  mu_T  <- param[2]
  theta <- param[3]
  
  times0 <- c(0, times[-length(times)])
  times1 <- times
  timesp <- times[length(times)]
  
  h_theta <- h(theta, beta)
  signo   <- sign(beta)
  
  sum((mu_T / (time_trunc)) * ((times1 - signo*h_theta* times0)^beta-(times0 - signo*h_theta* times0)^beta))   +   (mu_T / (time_trunc)) * ((time_trunc -  signo*h_theta*timesp)^beta-(timesp -  signo*h_theta*timesp)^beta) 
}


log_like_ARAM<- function(param, times, time_trunc) {
  
  
  log_like_single <- function(times_i) {
    if (length(times_i) == 0) return(0) 
    log_lambda_sum <- sum(log(lambda_ARAM(times_i, param, times_i, time_trunc)))
    
    
    lambda_total   <- Lambda_ARAM(param, times_i, time_trunc)
    log_lambda_sum - lambda_total
  }
  
  sum(sapply(times, log_like_single)) 
  }

