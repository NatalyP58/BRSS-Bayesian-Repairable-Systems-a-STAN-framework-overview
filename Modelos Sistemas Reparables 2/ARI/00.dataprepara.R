# Librerías necesarias solo para esta parte
library(pracma)

gen_failures_ARI <- function(k, param, time_trunc) {
  
  
  k_t <- function(t, times, time_trunc) {
    findInterval(t, c(0, times, time_trunc), left.open = TRUE)  
  }
  
  
  lambda_ARI_acum <- function(t, param, times, time_trunc) {
    
    if (t <= 0) return(0)
    
    
    
    beta  <- param[1]
    mu_T  <- param[2]
    rho   <- param[3]
    
    kt      <- k_t(t, times, time_trunc)
    times0  <- c(0, times)[1:(kt - 1)]
    times1  <- c(times, time_trunc)[1:(kt - 1)]
    times_k <- ifelse(kt - 1 > 0, c(times, time_trunc)[kt - 1], 0)
    
    
    if (kt - 1 > 0) {
      
      lambda1 <- ( mu_T / (time_trunc) ) * ((times1^beta - times0^beta)- rho * beta * times0^(beta-1) * (times1-times0))
      lambda1[is.na(lambda1)] <- 0
    } else {
      lambda1 <- 0
    }
    
    lambda2 <- ( mu_T / (time_trunc) )  * ((t^beta - (times_k)^beta)- rho * beta * (times_k)^(beta-1)*(t-times_k))
    lambda2[is.na(lambda2)] <- 0
    
    
    
    
    return(sum(lambda1) + lambda2)
  }
  lambda_ARI_acum <- Vectorize(lambda_ARI_acum, "t")
  
  
  objetivo_incremento <- function(x, param, t, times, time_trunc, u) {
    abs(lambda_ARI_acum(t + x, param, times, time_trunc) -
          lambda_ARI_acum(t, param, times, time_trunc) + log(1 - u))
  }
  
  
  simular_una_trayectoria <- function(param, time_trunc) {
    t_actual <- 0
    fallas   <- numeric(0)
    
    repeat {
      times    <- fallas 
      u        <- runif(1)
      objetivo <- -log(1 - u)
      
      max_Lambda <- lambda_ARI_acum(time_trunc, param, times, time_trunc) -
        lambda_ARI_acum(t_actual, param, times, time_trunc)
      
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
      
      fallas   <- c(fallas, t_siguiente)
      t_actual <- t_siguiente
    }
    
    return(fallas)
  }
  
  # Ejecutar simulación para k trayectorias
  replicate(k, simular_una_trayectoria(param, time_trunc), simplify = FALSE)
}


indicadora <- function(vector, t) {
  vals     <- vector[vector < t]
  if (length(vals) == 0) return(0)  
  max(vals)  
}


lambda_ARI <- function(t, param, times, time_trunc){
  if (t <= 0) return(0)
  beta  <- param[1]
  mu_T  <- param[2]
  rho   <- param[3]
  
  timess  <- indicadora(times, t) # timess= t_i-1
  
  
  mu_T * beta / (time_trunc) * (t^(beta-1)- rho*timess^(beta-1))
}

lambda_ARI <- Vectorize(lambda_ARI, vectorize.args = "t")

Lambda_ARI <- function(param, times, time_trunc){  #  times={t_1,....,tn} t_o=0 y t_n+1=Time_trunc=T
  if (length(times) == 0 || max(times) <= 0) return(0)
  beta  <- param[1]
  mu_T  <- param[2]
  rho   <- param[3]
  
  times0 <- c(0, times[-length(times)])
  times1 <- times
  timesp <- times[length(times)]
  
 
  sum(mu_T / (time_trunc) * ((times1^beta - times0^beta)-rho * beta * times0^(beta-1) * (times1 - times0)))   +   mu_T/ (time_trunc) * ((time_trunc^beta -  timesp^beta)-rho * beta * timesp^(beta-1) * (time_trunc -  timesp)) 

}

log_like_ARI <- function(param, times, time_trunc) {
  
  
  log_like_single <- function(times_i) {
    if (length(times_i) == 0) return(0) 
    log_lambda_sum <- sum(log(lambda_ARI(times_i, param, times_i, time_trunc)))
    
    
    lambda_total   <- Lambda_ARI(param, times_i, time_trunc)
    log_lambda_sum - lambda_total
  }
  
  sum(sapply(times, log_like_single)) 
  }


