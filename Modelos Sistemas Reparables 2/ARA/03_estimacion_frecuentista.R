# Librerías necesarias
library(optimx)
library(pracma)
#source("00.dataprepara.R")
# metodo frecuentista
frequentist_estimation <- function(time, time_trunc) {
  result <- try(optimx::optimr(par = c(10, 10, 0.8), fn = log_like_ARA,
                               method = "L-BFGS-B", time_trunc = time_trunc,
                               time = time, control = list(maximize = TRUE),
                               lower = c(1, 0, 0), upper = c(Inf, Inf, 1))$par,
                silent = TRUE)
  
  if (inherits(result, "try-error") || any(is.na(result))) {
    return(list(success = FALSE, estimate = NULL))
  }
  
  return(list(success = TRUE, estimate = as.numeric(result)))
}


