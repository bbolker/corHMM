args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L) {
  stop("Usage: Rscript tests/benchmark/run_real.R DATASET [STAGE] [PROFILE]")
}

dataset_name <- args[1]
stage <- if (length(args) >= 2L) args[2] else "all"
profile_name <- if (length(args) >= 3L) args[3] else "full"

source(file.path("tests", "benchmark", "methods.R"))
profile <- benchmark_profiles[[profile_name]]
if (is.null(profile)) stop("Unknown profile: ", profile_name)

if (dataset_name == "primates") {
  dataset <- load_primates_dataset()
  output_dir <- file.path("results", "benchmark", "primates")

  methods <- switch(
    stage,
    ml = c("cor_orig", "cor_rtmb", "cor_formula"),
    bayes = c("rtmb_mh", "rtmb_hmc", "rphylo", "btw"),
    rj = "btw_rj",
    all = c("cor_orig", "cor_rtmb", "cor_formula", "rtmb_mh", "rtmb_hmc", "rphylo", "btw"),
    stop("Unknown primates stage: ", stage)
  )

  for (method in methods) {
    warnings <- character()
    start <- proc.time()[[3]]

    result <- tryCatch(
      withCallingHandlers({
        fit <- run_method(
          method,
          dataset,
          profile,
          compute_ci = method %in% c("cor_orig", "cor_rtmb", "cor_formula")
        )
        raw <- fit$raw
        fit$raw <- NULL

        raw_file <- file.path(output_dir, paste0(method, "_raw.rds"))
        atomic_save_rds(raw, raw_file)

        list(
          status = "ok",
          dataset = dataset_name,
          stage = stage,
          profile = profile_name,
          method = method,
          fit = fit,
          raw_path = raw_file,
          warnings = warnings,
          metadata = benchmark_metadata()
        )
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        list(
          status = "error",
          dataset = dataset_name,
          stage = stage,
          profile = profile_name,
          method = method,
          error = conditionMessage(e),
          warnings = warnings,
          metadata = benchmark_metadata()
        )
      }
    )

    result$total_seconds <- proc.time()[[3]] - start
    atomic_save_rds(result, file.path(output_dir, paste0(method, ".rds")))
    cat(method, ":", result$status, "\n")
  }

} else if (dataset_name == "rayfish") {
  dataset <- load_rayfish_dataset()
  required <- c("ag", "care", "spawning")
  if (!all(required %in% dataset$trait_names)) {
    stop("c_traits.rds must contain Species, ag, care, and spawning columns.")
  }

  dataset$data <- dataset$data[, c("Species", required), drop = FALSE]
  dataset$trait_names <- required
  dataset$ntraits <- 3L
  output_dir <- file.path("results", "benchmark", "rayfish")

  fish_formula <- list(
    ag ~ care * spawning,
    care ~ 1,
    spawning ~ 1
  )

  methods <- switch(
    stage,
    ml = "fish_formula_ml",
    bayes = "fish_formula_hmc",
    application = c("fish_formula_ml", "fish_formula_hmc"),
    all = c("fish_formula_ml", "fish_formula_hmc"),
    stop("Unknown rayfish stage: ", stage)
  )

  for (method in methods) {
    warnings <- character()
    start <- proc.time()[[3]]

    result <- tryCatch(
      withCallingHandlers({
        if (method == "fish_formula_ml") {
          fit <- run_formula_ml_model(dataset, fish_formula, compute_ci = TRUE)
          main_effect <- NULL
        } else {
          fit <- run_rtmb_hmc(dataset, profile, formula_model = fish_formula)

          coefficient_names <- colnames(fit$coefficient_chains[[1]])
          care_index <- match("ag_gain:care", coefficient_names)
          interaction_index <- match("ag_gain:care:spawning", coefficient_names)
          if (is.na(interaction_index)) {
            interaction_index <- match("ag_gain:spawning:care", coefficient_names)
          }
          if (is.na(care_index)) stop("Could not locate the accessory-gland care coefficient.")

          ratio_chains <- lapply(fit$coefficient_chains, function(x) {
            log_ratio <- x[, care_index]
            if (!is.na(interaction_index)) log_ratio <- log_ratio + 0.5 * x[, interaction_index]
            matrix(exp(log_ratio), ncol = 1, dimnames = list(NULL, "care_effect_rate_ratio"))
          })
          main_effect <- summarize_rate_chains(ratio_chains)
          fit$main_effect <- main_effect
        }

        raw <- fit$raw
        fit$raw <- NULL
        raw_file <- file.path(output_dir, paste0(method, "_raw.rds"))
        atomic_save_rds(raw, raw_file)

        list(
          status = "ok",
          dataset = dataset_name,
          stage = stage,
          profile = profile_name,
          method = method,
          fit = fit,
          main_effect = main_effect,
          raw_path = raw_file,
          warnings = warnings,
          metadata = benchmark_metadata()
        )
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        list(
          status = "error",
          dataset = dataset_name,
          stage = stage,
          profile = profile_name,
          method = method,
          error = conditionMessage(e),
          warnings = warnings,
          metadata = benchmark_metadata()
        )
      }
    )

    result$total_seconds <- proc.time()[[3]] - start
    atomic_save_rds(result, file.path(output_dir, paste0(method, ".rds")))
    cat(method, ":", result$status, "\n")
  }

} else {
  stop("Unknown real dataset: ", dataset_name)
}
