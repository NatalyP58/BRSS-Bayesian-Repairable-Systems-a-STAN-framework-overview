# Librerías necesarias
library(optimx)
library(pracma)


# metodo frecuentista
frequentist_estimation <- function(time, time_trunc) {
  result <- try(optimx::optimr(par = c(5, 1), fn = log_like_PR,
                               method = "L-BFGS-B", time_trunc = time_trunc,
                               times = time, control = list(maximize = TRUE),
                               lower = c(1, 0), upper = c(Inf, Inf))$par,
                silent = TRUE)
  
  if (inherits(result, "try-error") || any(is.na(result))) {
    return(list(success = FALSE, estimate = NULL))
  }
  
  return(list(success = TRUE, estimate = as.numeric(result)))
}

