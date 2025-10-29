library(tidyverse)
library(dplyr)
# ==========================
# Leitura de dados e preparação
# ==========================

data <- read.csv("data.csv", sep=";")



vehiculo_sel <- "A"  # cambiar por el vehículo que quieras

# 1. Filtrar solo el vehículo seleccionado
data_veh <- data %>%
  filter(VEICULO == vehiculo_sel)

# 2. Convertir la columna de inicio a fecha/hora
data_veh$INICIO_PARADA2 <- as.POSIXct(data_veh$INICIO_PARADA, format="%d/%m/%Y %H:%M")

# --- Definir ventana de observación ---
start_obs <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M")
end_obs   <- as.POSIXct("2017-03-31 23:59", format="%Y-%m-%d %H:%M")

# 3. Calcular tiempo de falla desde la primera falla del vehículo
data_veh <- data_veh %>%
  arrange(INICIO_PARADA2) %>%
  mutate(time_failure_hours =  as.numeric(difftime(INICIO_PARADA2, start_obs, units = "hours")))

# 4. Sacar los 4 problemas más frecuentes dentro de ese vehículo
top4 <- data_veh %>%
  group_by(PROBLEMA) %>%
  summarise(n = n(), .groups="drop") %>%
  arrange(desc(n)) %>%
  slice(2:5) %>%
  pull(PROBLEMA)

# 5. Crear lista de tiempos por problema (cada vector inicia en t0=0)
times <- lapply(top4, function(prob) {
  tiempos <- data_veh %>%
    filter(PROBLEMA == prob) %>%
    arrange(time_failure_hours) %>%
    pull(time_failure_hours)
  tiempos <- tiempos[tiempos != min(tiempos)]   # asegurar que empiece en 0
  return(tiempos)
})



# 6. Nombrar la lista con los problemas
names(times) <- top4

# 7. Ver el resultado
times


# Ajustar tiempos quitandos los dos desfases mas grandes-----


ajustar_offseason <- function(x, start_obs, end_obs) {
  # Definir rangos de off-season (dic a mar) para cada año dentro del período
  years <- seq(from = as.numeric(format(start_obs, "%Y")),
               to   = as.numeric(format(end_obs, "%Y")),
               by = 1)
  
  # Lista de periodos de off-season
  off_seasons <- lapply(years, function(y) {
    inicio <- as.POSIXct(paste0(y, "-12-31 00:00"), tz = "UTC")
    fin    <- as.POSIXct(paste0(y+1, "-03-31 00:00"), tz = "UTC")
    c(inicio, fin)
  })
  
  # Filtrar los que caen dentro de la ventana de observación
  off_seasons <- Filter(function(r) r[1] < end_obs & r[2] > start_obs, off_seasons)
  
  # Transformar a data.frame con horas relativas a start_obs
  off_df <- data.frame(
    start = as.numeric(difftime(sapply(off_seasons, `[[`, 1), start_obs, units="hours")),
    end   = as.numeric(difftime(sapply(off_seasons, `[[`, 2), start_obs, units="hours"))
  )
  
  
  # Ajustar x eliminando huecos de off-season
  x_adj <- x
  dur_gap <- 0
  for (i in seq_len(nrow(off_df))) {
    despues <- off_df$end[i] - dur_gap
    dur_gap <- off_df$end[i] - off_df$start[i]
    x_adj[x_adj >= despues] <- x_adj[x_adj >= despues] - dur_gap
  }
  
  
  # Quitar eventos que cayeron dentro del off-season
  for (i in seq_len(nrow(off_df))) {
    x_adj <- x_adj[!(x >= off_df$start[i] & x < off_df$end[i])]
  }
  
  return(x_adj)
}



times_offseason <- lapply(times, function(vec) {
  ajustar_offseason(vec, start_obs, end_obs)
})

times_offseason  

end_obs
time_trunc <-  as.numeric(difftime(end_obs, start_obs, units = "hours"))



#----------------------------Inferencia

# Configuración inicial
current_dir <- getwd()

# Crear directorio temporal seguro para Stan
safe_temp_dir <- file.path(current_dir, "stan_temp")
if (!dir.exists(safe_temp_dir)) {
  dir.create(safe_temp_dir, recursive = TRUE)
}

# Configurar variables de entorno críticas
Sys.setenv(TMPDIR = safe_temp_dir)
Sys.setenv(TEMP = safe_temp_dir)
Sys.setenv(TMP = safe_temp_dir)
options(tmpdir = safe_temp_dir)

# Cargar librerías principales
library(parallel)
suppressPackageStartupMessages({
  library(rstan)
  library(coda)
})

# Configuraciones esenciales para Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
options(timeout = 600)

# Cargar scripts auxiliares usando rutas absolutas
source_files <- c(
  "00.dataprepara.R",
  "02_StanModels.R",
  "03_estimacion_frecuentista.R",
  "04_estimacion_bayesiana_corr.R",
  "05_estimacion_bayesiana_gamma.R"
)

for (file in source_files) {
  source(file.path(current_dir, file))
}

# Precompilar el modelo principal de forma segura
modelo.compilado.principal <- tryCatch({
  rstan::stan_model(
    model_code = modelo.stan.gamma,
    auto_write = FALSE,  # Evita escritura en disco
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error compilando modelo en nodo principal:", e$message, "\n")
  NULL
})

modelo.compilado.principal.corr <- tryCatch({
  rstan::stan_model(
    model_code = modelo.stan.corr,
    auto_write = FALSE,  # Evita escritura en disco
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error compilando modelo en nodo corr principal:", e$message, "\n")
  NULL
})

# Función principal de simulación paralela
run_application_parallel_cluster <- function(time , time_trunc, SEED = NULL) {
  
  # Si no se especifica, generar una semilla aleatoria reproducible
  if (is.null(SEED)) {
    SEED <- as.integer(Sys.time()) %% .Machine$integer.max
    cat("No se especificó semilla. Se generó automáticamente: ", SEED, "\n")
  } else {
    cat("Usando semilla fija: ", SEED, "\n")
  }
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  # Usar menos núcleos para mayor estabilidad
  ncores <- min(3, detectCores() - 1)
  cl     <- makeCluster(ncores, type = "PSOCK", outfile = file.path(current_dir, "cluster_log1.txt"))
  
  # Asegura corrientes reproducibles en cada worker
  parallel::clusterSetRNGStream(cl, iseed = SEED)
  
  # Exportar entorno a los nodos
  clusterExport(cl, varlist = c(
    "frequentist_estimation", "bayesian_estimation_gamma", "bayesian_estimation_corr",
    "Lambda_ARA", "lambda_ARA", "log_like_ARA", 
    "modelo.stan.gamma", "current_dir", "safe_temp_dir",
    "modelo.compilado.principal", "modelo.compilado.principal.corr", "source_files",  "SEED"
  ), envir = environment())
  
  # Configuración de nodos
  clusterEvalQ(cl, {
    Sys.setenv(TMPDIR = safe_temp_dir, TEMP = safe_temp_dir, TMP = safe_temp_dir)
    options(tmpdir = safe_temp_dir, mc.cores = 1, timeout = 600)
    rstan::rstan_options(auto_write = TRUE)
    Sys.setenv(STAN_NUM_THREADS = 1)
    
    suppressPackageStartupMessages({
      library(rstan)
      library(coda)
    })
    
    # Compilar modelos si no están disponibles
    if (!is.null(modelo.compilado.principal)) {
      assign("modelo.compilado.gamma", modelo.compilado.principal, envir = .GlobalEnv)
    } else {
      assign("modelo.compilado.gamma", rstan::stan_model(model_code = modelo.stan.gamma), envir = .GlobalEnv)
    }
    
    if (!is.null(modelo.compilado.principal.corr)) {
      assign("modelo.compilado.corr", modelo.compilado.principal.corr, envir = .GlobalEnv)
    } else {
      assign("modelo.compilado.corr", rstan::stan_model(model_code = modelo.stan.corr), envir = .GlobalEnv)
    }
  })
  
  # Estimación frecuentista
  freq <- tryCatch({
    frequentist_estimation(time, time_trunc)
  }, error = function(e) {
    cat("Error en estimación frecuentista:", e$message, "\n")
    list(success = FALSE)
  })
  
  # Estimación bayesiana gamma
  bayes_gamma <- tryCatch({
    bayesian_estimation_gamma(
      time, time_trunc,  beta_r = NULL, mu_T_r = NULL, theta_r = NULL,
      compiled_model = modelo.compilado.principal, simulation = FALSE,
      seed = SEED)
  }, error = function(e) {
    cat("Error en estimación bayesiana gamma:", e$message, "\n")
    list(success = FALSE)
  })
  
  # Estimación bayesiana con correlación
  bayes_corr <- tryCatch({
    bayesian_estimation_corr(
      time, time_trunc,beta_r = NULL, mu_T_r = NULL, theta_r = NULL,
      compiled_model_corr = modelo.compilado.principal.corr, simulation = FALSE,
      seed = SEED
    )
  }, error = function(e) {
    cat("Error en estimación bayesiana corr:", e$message, "\n")
    list(success = FALSE)
  })

  stopCluster(cl)
  
  # Resultados finales
  results <- list(
    freq_est        = if (freq$success) freq$estimate else NA,
    bayes_est_gamma = if (bayes_gamma$success) bayes_gamma$estimate else NA,
    ajuste_gamma    = bayes_gamma$fit,
    coverage_gamma  = if (bayes_gamma$success) bayes_gamma$coverage else NA,
    bayes_est_corr  = if (bayes_corr$success) bayes_corr$estimate else NA,
    ajuste_corr     = bayes_corr$fit,
    coverage_corr   = if (bayes_corr$success) bayes_corr$coverage else NA,
    success_freq    = freq$success,
    success_gamma   = bayes_gamma$success,
    success_corr    = bayes_corr$success,
    seed_used       = SEED
  )
  
  return(results)
}

# Ejecución principal con manejo de errores
tryCatch({
  results <- run_application_parallel_cluster (
    time       = times_offseason,
    time_trunc = time_trunc, 
    SEED       = 20251025
  )
  
  # Guardar resultados
  saveRDS(results, file = file.path(current_dir, "resultados_application01.rds"))
  
  # Limpieza final
  unlink(safe_temp_dir, recursive = TRUE)
}, error = function(e) {
  cat("Error en ejecución principal:", e$message, "\n")
}, finally = {
  # Intentar limpiar incluso si hay error
  if (dir.exists(safe_temp_dir)) {
    unlink(safe_temp_dir, recursive = TRUE)
  }
})

# AIC y BIC

calc_AIC_BIC <- function(results, times, time_trunc) {
  
  # Número de parámetros estimados
  k <- length(results$freq_est) 
  
  # Número de observaciones (total de tiempos en todos los sistemas)
  n <- sum(sapply(times, length))
  
  # --- Frecuentista ---
  loglik_freq <- log_like_ARA(results$freq_est, time=times, time_trunc=time_trunc)
  AIC_freq <- 2 * k - 2 * loglik_freq
  BIC_freq <- k * log(n) - 2 * loglik_freq
  
  # --- Bayesiano (gamma) ---
  loglik_bayes <- log_like_ARA(results$bayes_est_gamma, time=times, time_trunc=time_trunc)
  AIC_bayes <- 2 * k - 2 * loglik_bayes
  BIC_bayes <- k * log(n) - 2 * loglik_bayes
  
  # # --- Bayesiano (corr) ---
  loglik_bayes_corr <- log_like_ARA(results$bayes_est_corr, time=times, time_trunc=time_trunc)
  AIC_bayes_corr    <- 2 * k - 2 * loglik_bayes_corr
  BIC_bayes_corr    <- k * log(n) - 2 * loglik_bayes_corr
  #  
  return(list(
    n = n,
    k = k,
    AIC_freq = AIC_freq,
    BIC_freq = BIC_freq,
    AIC_bayes = AIC_bayes,
    BIC_bayes = BIC_bayes,
    AIC_bayes_corr = AIC_bayes_corr,
    BIC_bayes_corr = BIC_bayes_corr
  ))
}




ic_vals <- calc_AIC_BIC(results, times=times_offseason,  time_trunc = time_trunc)
print(ic_vals)






#DIC

calc_DIC <- function(fit_stan, times, time_trunc) {
  
  # Extraer draws de los parámetros
  draws <- as.data.frame(fit_stan)
  param_names <- c("beta", "mu_T", "theta")
  draws <- draws[, param_names]
  
  # Deviance en cada muestra posterior
  dev_samples <- apply(draws, 1, function(par){
    -2 * log_like_ARA(par, time=times, time_trunc=time_trunc)
  })
  
  # Promedio de deviance
  D_bar <- mean(dev_samples)
  
  # Deviance en la media posterior
  par_mean <- colMeans(draws)
  D_hat <- -2 * log_like_ARA(par_mean, time=times, time_trunc=time_trunc)
  
  # Número efectivo de parámetros
  p_D <- D_bar - D_hat
  
  # DIC
  DIC <- D_bar + p_D
  
  return(list(D_bar = D_bar, D_hat = D_hat, p_D = p_D, DIC = DIC))
}


calc_DIC(results$ajuste_gamma,times=times_offseason,  time_trunc = time_trunc)

calc_DIC(results$ajuste_corr, times=times_offseason,  time_trunc = time_trunc)



