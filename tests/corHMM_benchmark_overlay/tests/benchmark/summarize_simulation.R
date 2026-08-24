source(file.path("tests", "benchmark", "common.R"))

input_dir <- file.path("results", "benchmark", "simulation")
output_dir <- file.path("results", "benchmark", "summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("_raw\\.rds$", files)]
results <- lapply(files, readRDS)
ok <- vapply(results, function(x) identical(x$status, "ok"), logical(1))

if (any(!ok)) {
  errors <- do.call(rbind, lapply(results[!ok], function(x) {
    data.frame(
      job_id = x$spec$job_id,
      method = x$spec$method,
      error = x$error,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(errors, file.path(output_dir, "simulation_errors.csv"), row.names = FALSE)
}
results <- results[ok]

ml <- results[vapply(results, function(x) identical(x$spec$analysis, "ml"), logical(1))]
if (length(ml)) {
  ml_fits <- do.call(rbind, lapply(ml, function(x) {
    data.frame(
      dataset_id = x$spec$dataset_id,
      ntraits = x$spec$ntraits,
      ntips = x$spec$ntips,
      dataset_seed = x$spec$run_seed,
      method = x$spec$method,
      fit_seconds = x$fit$fit_seconds,
      optimize_seconds = x$fit$optimize_seconds,
      loglik = x$fit$loglik,
      nll = -x$fit$loglik,
      convergence = x$fit$convergence,
      iterations = x$fit$iterations,
      stringsAsFactors = FALSE
    )
  }))

  ml_rates <- do.call(rbind, lapply(ml, function(x) {
    data.frame(
      dataset_id = x$spec$dataset_id,
      ntraits = x$spec$ntraits,
      ntips = x$spec$ntips,
      dataset_seed = x$spec$run_seed,
      method = x$spec$method,
      parameter = names(x$fit$rates),
      estimate = as.numeric(x$fit$rates),
      stringsAsFactors = FALSE
    )
  }))

  grouped <- split(ml, vapply(ml, function(x) x$spec$dataset_id, character(1)))
  pairwise <- list()
  winners <- list()
  k <- 1L
  w <- 1L

  for (dataset_id in names(grouped)) {
    z <- grouped[[dataset_id]]
    methods <- vapply(z, function(x) x$spec$method, character(1))
    best_nll <- min(vapply(z, function(x) -x$fit$loglik, numeric(1)))
    winner <- methods[abs(vapply(z, function(x) -x$fit$loglik, numeric(1)) - best_nll) <= 1e-6]
    winners[[w]] <- data.frame(
      dataset_id = dataset_id,
      winner = paste(winner, collapse = ";"),
      tie = length(winner) > 1L,
      stringsAsFactors = FALSE
    )
    w <- w + 1L

    if (length(z) >= 2L) {
      for (i in seq_len(length(z) - 1L)) {
        for (j in (i + 1L):length(z)) {
          pairwise[[k]] <- data.frame(
            dataset_id = dataset_id,
            ntraits = z[[i]]$spec$ntraits,
            ntips = z[[i]]$spec$ntips,
            method_a = z[[i]]$spec$method,
            method_b = z[[j]]$spec$method,
            delta_loglik = z[[i]]$fit$loglik - z[[j]]$fit$loglik,
            abs_delta_loglik = abs(z[[i]]$fit$loglik - z[[j]]$fit$loglik),
            rate_rmse = rate_rmse(z[[i]]$fit$rates, z[[j]]$fit$rates),
            max_abs_rate_difference = max(abs(
              z[[i]]$fit$rates[names(z[[j]]$fit$rates)] - z[[j]]$fit$rates
            )),
            stringsAsFactors = FALSE
          )
          k <- k + 1L
        }
      }
    }
  }

  utils::write.csv(ml_fits, file.path(output_dir, "simulation_ml_fits.csv"), row.names = FALSE)
  utils::write.csv(ml_rates, file.path(output_dir, "simulation_ml_rates.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, pairwise), file.path(output_dir, "simulation_ml_pairwise.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, winners), file.path(output_dir, "simulation_ml_winners.csv"), row.names = FALSE)
}

bayes <- results[vapply(results, function(x) identical(x$spec$analysis, "bayes"), logical(1))]
if (length(bayes)) {
  rate_summary <- do.call(rbind, lapply(bayes, function(x) {
    cbind(
      dataset_id = x$spec$dataset_id,
      ntips = x$spec$ntips,
      dataset_seed = x$spec$run_seed,
      method = x$spec$method,
      x$fit$posterior_summary
    )
  }))
  rownames(rate_summary) <- NULL

  method_summary <- do.call(rbind, lapply(bayes, function(x) {
    s <- x$fit$posterior_summary
    d <- x$fit$diagnostics
    native <- native_output_summary(x$fit)
    data.frame(
      dataset_id = x$spec$dataset_id,
      ntips = x$spec$ntips,
      dataset_seed = x$spec$run_seed,
      method = x$spec$method,
      prior_mean = as.numeric(x$spec$prior_mean),
      prior_sd = as.numeric(x$spec$prior_sd),
      fit_seconds = x$fit$fit_seconds,
      optimize_seconds = if (!is.null(x$fit$optimize_seconds)) x$fit$optimize_seconds else NA_real_,
      sampling_seconds = if (!is.null(x$fit$sampling_seconds)) x$fit$sampling_seconds else NA_real_,
      max_rhat = safe_max(s$rhat),
      min_bulk_ess = safe_min(s$ess_bulk),
      min_tail_ess = safe_min(s$ess_tail),
      min_bulk_ess_per_second = safe_min(s$ess_bulk) / x$fit$fit_seconds,
      min_tail_ess_per_second = safe_min(s$ess_tail) / x$fit$fit_seconds,
      max_mcse_over_sd = safe_max(s$mcse_mean / s$sd),
      divergences = if (!is.null(d$divergences)) d$divergences else NA_real_,
      treedepth_hits = if (!is.null(d$treedepth_hits)) d$treedepth_hits else NA_real_,
      min_ebfmi = if (!is.null(d$min_ebfmi)) d$min_ebfmi else NA_real_,
      mean_acceptance = if (!is.null(d$post_warmup_acceptance_rate)) {
        mean(d$post_warmup_acceptance_rate, na.rm = TRUE)
      } else if (!is.null(d$mean_accept)) {
        d$mean_accept
      } else {
        NA_real_
      },
      native_quantity = native$native_quantity,
      native_mean = native$native_mean,
      native_sd = native$native_sd,
      native_min = native$native_min,
      native_max = native$native_max,
      stringsAsFactors = FALSE
    )
  }))

  grouped <- split(bayes, vapply(bayes, function(x) x$spec$dataset_id, character(1)))
  pairwise <- list()
  k <- 1L

  for (dataset_id in names(grouped)) {
    z <- grouped[[dataset_id]]
    if (length(z) < 2L) next

    for (i in seq_len(length(z) - 1L)) {
      for (j in (i + 1L):length(z)) {
        a <- z[[i]]$fit$posterior_summary
        b <- z[[j]]$fit$posterior_summary
        mean_a <- setNames(a$mean, a$parameter)
        mean_b <- setNames(b$mean, b$parameter)
        sd_a <- setNames(a$sd, a$parameter)
        sd_b <- setNames(b$sd, b$parameter)

        pairwise[[k]] <- data.frame(
          dataset_id = dataset_id,
          ntips = z[[i]]$spec$ntips,
          method_a = z[[i]]$spec$method,
          method_b = z[[j]]$spec$method,
          posterior_mean_rmse = rate_rmse(mean_a, mean_b),
          posterior_sd_rmse = rate_rmse(sd_a, sd_b),
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }

  utils::write.csv(rate_summary, file.path(output_dir, "simulation_bayes_rates.csv"), row.names = FALSE)
  utils::write.csv(method_summary, file.path(output_dir, "simulation_bayes_methods.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, pairwise), file.path(output_dir, "simulation_bayes_pairwise.csv"), row.names = FALSE)
}

rj <- results[vapply(results, function(x) identical(x$spec$analysis, "rj"), logical(1))]
if (length(rj)) {
  rj_summary <- do.call(rbind, lapply(rj, function(x) {
    data.frame(
      dataset_id = x$spec$dataset_id,
      ntips = x$spec$ntips,
      dataset_seed = x$spec$run_seed,
      method = x$spec$method,
      fit_seconds = x$fit$fit_seconds,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(rj_summary, file.path(output_dir, "simulation_rj_methods.csv"), row.names = FALSE)

  model_counts <- list()
  k <- 1L
  for (x in rj) {
    for (chain in seq_along(x$fit$rj)) {
      z <- x$fit$rj[[chain]]
      if (!ncol(z)) next
      z$chain <- chain
      z$dataset_id <- x$spec$dataset_id
      key_columns <- setdiff(names(z), "chain")
      key_columns <- setdiff(key_columns, "dataset_id")
      key <- do.call(paste, c(z[key_columns], sep = "|"))
      tab <- table(key)
      first <- match(names(tab), key)
      out <- z[first, c("dataset_id", key_columns), drop = FALSE]
      out$chain <- chain
      out$count <- as.integer(tab)
      out$probability <- out$count / sum(out$count)
      model_counts[[k]] <- out
      k <- k + 1L
    }
  }
  if (length(model_counts)) {
    utils::write.csv(
      do.call(rbind, model_counts),
      file.path(output_dir, "simulation_rj_model_counts.csv"),
      row.names = FALSE
    )
  }
}

cat("Simulation summaries written to", output_dir, "\n")
