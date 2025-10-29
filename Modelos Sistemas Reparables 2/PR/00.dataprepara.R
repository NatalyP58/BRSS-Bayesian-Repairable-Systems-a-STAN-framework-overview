library(pracma)

gen_failures_PR <-function(k,param, time_trunc) {
  
  k_t <- function(t, times, time_trunc) {
    findInterval(t, c(0, times, time_trunc), left.open = TRUE) 
  }
  
  
  lambda_PR_acum <- function(t, param, times, time_trunc) {
    if (t <= 0) return(0)
    beta <- param[1]
    mu_T <- param[2]
    
    kt      <- k_t(t, times, time_trunc)
    times0  <- c(0, times)[1:(kt - 1)]
    times1  <- c(times, time_trunc)[1:(kt - 1)]
    times_k <- ifelse(kt - 1 > 0, c(times, time_trunc)[kt - 1], 0)
    
    if (kt - 1 > 0) {
      lambda1 <- mu_T * ((times1 - times0)^beta) / time_trunc
      lambda1[is.na(lambda1)] <- 0
    } else {
      lambda1 <- 0
    }
    
    lambda2 <- mu_T * ((t - times_k)^beta) / time_trunc
    lambda2[is.na(lambda2)] <- 0
    
    return(sum(lambda1) + lambda2)
  }
  
  lambda_PR_acum <- Vectorize(lambda_PR_acum, vectorize.args = "t")
  
  
  
  optimizar <- function(x, param, t, times, time_trunc, u) {
    abs(lambda_PR_acum(t + x, param, times, time_trunc) -
          lambda_PR_acum(t, param, times, time_trunc) + log(1 - u))
  }
  
  
  simular_una_trayectoria <- function(param, time_trunc) {
    t_actual <- 0
    times    <- numeric(0)
    
    repeat {
      u          <- runif(1)
      objetivo   <- -log(1 - u)
      
      max_Lambda <- lambda_PR_acum(time_trunc, param, times, time_trunc) -
                    lambda_PR_acum(t_actual, param, times, time_trunc)
      
      if (objetivo > max_Lambda) break
      
      resultado <- optimize(
        optimizar,
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
  replicate(k, simular_una_trayectoria(param, time_trunc), simplify = FALSE)
}


indicadora <- function(vector, t) {
  vals     <- vector[vector < t]
  if (length(vals) == 0) return(0)  
  max(vals)
}

lambda_PR <- function(t, param, times, time_trunc){
  if (t <= 0) return(0)
  beta <- param[1]
  mu_T <- param[2]
  
  timess <- indicadora(times, t) # timess= t_i-1
  
  mu_T * beta / (time_trunc) * (t - timess)^(beta-1)
}

lambda_PR <- Vectorize(lambda_PR, vectorize.args = "t")

Lambda_PR <-function(param, times, time_trunc){  # aqui times={t_1,....,tn} t_o=0 y t_n+1=Time_trunc=T
  if (length(times) == 0 || max(times) <= 0) return(0)
  beta <- param[1]
  mu_T <- param[2]
  
  times0 <- c(0, times[-length(times)])
  times1 <- times
  timesp <- times[length(times)]
  
  sum(mu_T / (time_trunc) * ((times1 - times0)^beta))   +   mu_T * ((time_trunc -  timesp)^beta) / (time_trunc)
}


log_like_PR <- function(param, times, time_trunc) {
  
  
  log_like_single <- function(times_i) {
    if (length(times_i) == 0) return(0) # Si no hay fallas, el aporte a la verosimilitud es 0
    
    
    log_lambda_sum <- sum(log(lambda_PR(times_i, param, times_i, time_trunc)))
    
    
    lambda_total   <- Lambda_PR(param, times_i, time_trunc)
    
    
    log_lambda_sum - lambda_total
  }
 
  sum(sapply(times, log_like_single))   
}



