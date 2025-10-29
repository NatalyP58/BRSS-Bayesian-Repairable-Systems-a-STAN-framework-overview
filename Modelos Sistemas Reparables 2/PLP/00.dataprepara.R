# Librerías necesPLPas solo para esta parte
library(pracma)

gen_failures_PLP <- function(k, param, time_trunc) {
  beta <- param[1]
  mu_T <- param[2]
  
  generate_system <- function() {
    t             <- numeric(0)
    current_time  <- 0
    
    while(TRUE) {   #Λ(t+x)−Λ(t)=−ln(1−u)===>x=[t^β +T/μ (−ln(1−u))]^(1/β) -t
      u         <- runif(1)
      increment <- (current_time^beta + (time_trunc/mu_T)*(-log(1 - u)))^(1/beta) - current_time
      
      if(current_time + increment > time_trunc) break
      
      current_time <- current_time + increment
      t            <- c(t, current_time)
    }
    t
  }
  
  replicate(k, generate_system(), simplify = FALSE)  
}

lambda_PLP <- function(t, param, time_trunc){
  beta <- param[1]
  mu_T <- param[2]
  
  return((mu_T * beta/ time_trunc) * t^(beta-1))  
}

lambda_PLP <- Vectorize(lambda_PLP, vectorize.args = "t")

Lambda_PLP <- function(param, time_trunc){ #Λ_PLP(T)
  beta <- param[1]
  mu_T <- param[2]
  
  
  
  return((mu_T / time_trunc) * time_trunc^beta)
}

log_like_PLP <- function(param, time_trunc, time){
  
  log_like_single <- function(single_time) {
    sum(log(lambda_PLP(single_time, param, time_trunc))) - Lambda_PLP(param, time_trunc)
  }
  
  return(sum(sapply(time, log_like_single)))
}



