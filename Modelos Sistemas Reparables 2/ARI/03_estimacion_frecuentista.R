# Librerías necesarias
library(optimx)
library(pracma)
#source("00.dataprepara.R")
# metodo frecuentista
frequentist_estimation <- function(time, time_trunc) {
  result <- try(optimx::optimr(par = c(1.5, 5, 0.2), fn = log_like_ARI,
                               method = "L-BFGS-B", time_trunc = time_trunc,
                               times = time, control = list(maximize = TRUE),
                               lower = c(0, 0, -1), upper = c(Inf, Inf, 1))$par,
                silent = TRUE)
  
  if (inherits(result, "try-error") || any(is.na(result))) {
    return(list(success = FALSE, estimate = NULL))
  }
  
  return(list(success = TRUE, estimate = as.numeric(result)))
}


