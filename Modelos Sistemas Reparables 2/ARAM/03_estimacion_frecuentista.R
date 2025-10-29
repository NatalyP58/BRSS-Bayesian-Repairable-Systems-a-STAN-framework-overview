# Librerías necesarias
library(optimx)
library(pracma)
#source("00.dataprepara.R")
# metodo frecuentista
frequentist_estimation <- function(time, time_trunc) {
  result <- try(optimx::optimr(par = c(2, 20, 0.8), fn = log_like_ARAM,
                               method = "L-BFGS-B", time_trunc = time_trunc,
                               time = time, control = list(maximize = TRUE),
                               lower = c(0, 0, -1), upper = c(Inf, Inf, 1))$par,
                silent = TRUE)
  
  if (inherits(result, "try-error") || any(is.na(result))) {
    return(list(success = FALSE, estimate = NULL))
  }
  
  return(list(success = TRUE, estimate = as.numeric(result)))
}


