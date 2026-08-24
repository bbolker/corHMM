args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args)) args[1] else "full"

source(file.path("tests", "benchmark", "config.R"))
if (!profile %in% names(benchmark_profiles)) stop("Unknown benchmark profile.")

dataset_file <- file.path("tests", "benchmark", "manifests", "datasets_simulation.csv")
if (!file.exists(dataset_file)) stop("Run generate_simulation_data.R first.")

datasets <- utils::read.csv(dataset_file, stringsAsFactors = FALSE)

make_jobs <- function(data, methods, analysis, save_raw = FALSE) {
  rows <- vector("list", nrow(data) * length(methods))
  k <- 1L

  for (i in seq_len(nrow(data))) {
    lower <- log(1)
    upper <- log(100 * data$ntips[i])

    for (method in methods) {
      rows[[k]] <- data.frame(
        job_id = sprintf("%s_%s_%s", analysis, data$dataset_id[i], method),
        dataset_id = data$dataset_id[i],
        analysis = analysis,
        method = method,
        profile = profile,
        run_seed = data$dataset_seed[i],
        ntraits = data$ntraits[i],
        ntips = data$ntips[i],
        tree_path = data$tree_path[i],
        data_path = data$data_path[i],
        prior_lower = lower,
        prior_upper = upper,
        prior_mean = (lower + upper) / 2,
        prior_sd = (upper - lower) / 3,
        compute_ci = FALSE,
        save_raw = save_raw,
        output_path = benchmark_path(
          "results", "benchmark", "simulation", analysis,
          paste0(data$dataset_id[i], "__", method, ".rds")
        ),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}

ml <- make_jobs(
  datasets,
  methods = c("cor_orig", "cor_rtmb", "cor_formula"),
  analysis = "ml"
)

bayes_data <- datasets[datasets$ntraits == 2L, , drop = FALSE]

bayes <- make_jobs(
  bayes_data,
  methods = c("rtmb_mh", "rtmb_hmc", "rphylo", "btw"),
  analysis = "bayes"
)

rj <- make_jobs(
  bayes_data,
  methods = "btw_rj",
  analysis = "rj"
)

manifest_dir <- file.path("tests", "benchmark", "manifests")
utils::write.csv(ml, file.path(manifest_dir, paste0("jobs_ml_", profile, ".csv")), row.names = FALSE)
utils::write.csv(bayes, file.path(manifest_dir, paste0("jobs_bayes_", profile, ".csv")), row.names = FALSE)
utils::write.csv(rj, file.path(manifest_dir, paste0("jobs_rj_", profile, ".csv")), row.names = FALSE)

cat("ML jobs:", nrow(ml), "\n")
cat("Bayesian fixed-model jobs:", nrow(bayes), "\n")
cat("RJMCMC jobs:", nrow(rj), "\n")
