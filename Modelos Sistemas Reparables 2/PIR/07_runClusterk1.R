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
run_simulation_parallel_cluster <- function(m, k, param, time_trunc, tau, SEED = NULL) {
  
  # Si no se especifica, generar una semilla aleatoria reproducible
  if (is.null(SEED)) {
    SEED <- as.integer(Sys.time()) %% .Machine$integer.max
    cat("No se especificó semilla. Se generó automáticamente: ", SEED, "\n")
  } else {
    cat("Usando semilla fija: ", SEED, "\n")
  }
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  
  # Semilla por cada simulación i = 1..m
  seeds <- SEED + seq_len(m)
  
  # Usar menos núcleos para mayor estabilidad
  ncores <- min(3, detectCores() - 1)
  cl <- makeCluster(ncores, type = "PSOCK", outfile = file.path(current_dir, "cluster_log1.txt"))
  
  # Asegura corrientes reproducibles en cada worker
  parallel::clusterSetRNGStream(cl, iseed = SEED)
  
  # Exportar entorno a los nodos
  clusterExport(cl, varlist = c(
    "gen_failures_PIR", "frequentist_estimation", "bayesian_estimation_gamma",
    "Lambda_PIR", "lambda_PIR", "log_like_PIR", "modelo.stan.gamma", 
    "current_dir", "safe_temp_dir", "modelo.compilado.principal", "modelo.compilado.principal.corr", "source_files", "seeds", "SEED"), 
    envir = environment())
  
  # Configuración de nodos
  clusterEvalQ(cl, {
    # Configurar directorio temporal en nodos
    Sys.setenv(TMPDIR = safe_temp_dir)
    Sys.setenv(TEMP = safe_temp_dir)
    Sys.setenv(TMP = safe_temp_dir)
    options(tmpdir = safe_temp_dir)
    
    # Cargar librerías
    suppressPackageStartupMessages({
      library(rstan)
      library(coda)
    })
    
    # Configuraciones Stan
    options(mc.cores = 1)
    options(timeout = 600)
    rstan::rstan_options(auto_write = TRUE)
    Sys.setenv(STAN_NUM_THREADS = 1)
    
    # Usar modelo precompilado o compilar localmente
    if (!is.null(modelo.compilado.principal)) {
      assign("modelo.compilado.gamma", modelo.compilado.principal, envir = .GlobalEnv)
    } else {
      tryCatch({
        assign("modelo.compilado.gamma", rstan::stan_model(
          model_code = modelo.stan.gamma,
          auto_write = FALSE,  # Evita escritura en disco
          verbose = FALSE
        ), envir = .GlobalEnv)
      }, error = function(e) {
        cat("Error compilando modelo en nodo:", e$message, "\n")
        assign("modelo.compilado.gamma", NULL, envir = .GlobalEnv)
      })
    }
    # Usar modelo precompilado o compilar localmente
    if (!is.null(modelo.compilado.principal.corr)) {
      assign("modelo.compilado.corr", modelo.compilado.principal.corr, envir = .GlobalEnv)
    } else {
      tryCatch({
        assign("modelo.compilado.corr", rstan::stan_model(
          model_code = modelo.stan.corr,
          auto_write = FALSE,  # Evita escritura en disco
          verbose = FALSE
        ), envir = .GlobalEnv)
      }, error = function(e) {
        cat("Error compilando modelo corr en nodo:", e$message, "\n")
        assign("modelo.compilado.corr", NULL, envir = .GlobalEnv)
      })
    }
    
    # Cargar funciones
    source_files <- c(
      "00.dataprepara.R",
      "05_estimacion_bayesiana_gamma.R",
      "04_estimacion_bayesiana_corr.R"
    )
    for (file in source_files) {
      source(file.path(current_dir, file))
    }
  })
  
  # Función de simulación individual
  single_run <- function(i) {
    set.seed(seeds[i])
    
    cat("Iniciando simulacion ", i, "\n")
    beta_r  <- param[1]
    mu_T_r  <- param[2]
    theta_r <- param[3]
    
    # Generar datos
    time <- tryCatch({
      gen_failures_PIR(k, c(beta_r, mu_T_r, theta_r), time_trunc, tau)
    }, error = function(e) {
      cat("Error generando datos:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(time)) {
      return(list(
        success_freq  = FALSE,
        success_gamma = FALSE,
        success_corr  = FALSE
      ))
    }
    
    # Estimación frecuentista
    freq <- tryCatch({
      frequentist_estimation(time, time_trunc, tau)
    }, error = function(e) {
      cat("Error en estimación frecuentista:", e$message, "\n")
      list(success = FALSE)
    })
    
    # Estimación bayesiana
    bayes_gamma <- tryCatch({
      if (exists("modelo.compilado.gamma") && !is.null(modelo.compilado.gamma)) {
        bayesian_estimation_gamma(
          time, time_trunc, beta_r, mu_T_r, theta_r, tau,
          compiled_model = modelo.compilado.gamma,
          seed = seeds[i] 
        )
      } else {
        stop("Modelo compilado no disponible")
      }
    }, error = function(e) {
      cat("Error en estimación bayesiana:", e$message, "\n")
      list(success = FALSE)
    })
    
    # Estimación bayesiana correlacion
    bayes_corr <- tryCatch({
      if (exists("modelo.compilado.corr") && !is.null(modelo.compilado.corr)) {
        bayesian_estimation_corr(
          time, time_trunc, beta_r, mu_T_r, theta_r, tau,
          compiled_model_corr = modelo.compilado.corr,
          seed = seeds[i] 
        )
      } else {
        stop("Modelo compilado no disponible")
      }
    }, error = function(e) {
      cat("Error en estimación bayesiana correlación:", e$message, "\n")
      list(success = FALSE)
    })
    
    
    # Resultados
    list(
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
      seed_used       = seeds[i]
    )
  }
  
  # Ejecutar simulaciones
  results <- tryCatch({
    parLapplyLB(cl, 1:m, single_run)
  }, finally = {
    stopCluster(cl)
  })
  
  # Estadísticas finales
  conv_freq  <- sum(sapply(results, function(x) x$success_freq))
  conv_bayes <- sum(sapply(results, function(x) x$success_gamma))
  conv_bayes_corr <- sum(sapply(results, function(x) x$success_corr))
  
  cat("\n--- RESUMEN FINAL ---\n")
  cat("Simulaciones completadas:", m, "\n")
  cat("Exitosis frecuentistas:", conv_freq, "\n")
  cat("Exitos Bayesianos:", conv_bayes, "\n")
  cat("Exitos Bayesianos con correlacion:", conv_bayes_corr, "\n")
  
  return(results)
}

# Ejecución principal con manejo de errores
tryCatch({
  results <- run_simulation_parallel_cluster(
    m = 1, 
    k = 2, 
    param = c(2, 30, 0.5), 
    time_trunc = 40, 
    tau = c(15, 25), 
    SEED = 20251024
  )
  
  # Guardar resultados
  saveRDS(results, file = file.path(current_dir, "resultadosk5_simulacion01.rds"))

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

