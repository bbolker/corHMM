args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript tests/benchmark/run_job.R MANIFEST ROW")
}

manifest_file <- args[1]
row_index <- as.integer(args[2])

source(file.path("tests", "benchmark", "methods.R"))

manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
if (row_index < 1L || row_index > nrow(manifest)) stop("Manifest row is out of range.")
job <- manifest[row_index, , drop = FALSE]
profile <- benchmark_profiles[[job$profile]]
if (is.null(profile)) stop("Unknown benchmark profile: ", job$profile)

dataset <- load_manifest_dataset(job)
dataset$prior <- list(
  lower = as.numeric(job$prior_lower),
  upper = as.numeric(job$prior_upper),
  mean = as.numeric(job$prior_mean),
  sd = as.numeric(job$prior_sd)
)
warnings <- character()
start_time <- proc.time()[[3]]

result <- tryCatch(
  withCallingHandlers({
    fit <- run_method(
      method = job$method,
      dataset = dataset,
      profile = profile,
      compute_ci = isTRUE(as.logical(job$compute_ci))
    )

    raw <- fit$raw
    fit$raw <- NULL

    if (isTRUE(as.logical(job$save_raw))) {
      raw_path <- sub("\\.rds$", "_raw.rds", job$output_path)
      atomic_save_rds(raw, raw_path)
    } else {
      raw_path <- NA_character_
    }

    list(
      status = "ok",
      spec = as.list(job),
      dataset = list(
        ntraits = dataset$ntraits,
        ntips = dataset$ntips,
        trait_names = dataset$trait_names,
        original_total_branch_length = dataset$original_total_branch_length,
        scale_factor = dataset$scale_factor,
        scaled_total_branch_length = dataset$scaled_total_branch_length
      ),
      fit = fit,
      raw_path = raw_path,
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
      spec = as.list(job),
      error = conditionMessage(e),
      warnings = warnings,
      metadata = benchmark_metadata()
    )
  }
)

result$total_seconds <- proc.time()[[3]] - start_time
atomic_save_rds(result, job$output_path)

cat("job:", job$job_id, "\n")
cat("status:", result$status, "\n")
cat("output:", job$output_path, "\n")
cat("total seconds:", result$total_seconds, "\n")
