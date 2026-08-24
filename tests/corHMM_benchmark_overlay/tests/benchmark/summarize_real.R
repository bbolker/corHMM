source(file.path("tests", "benchmark", "common.R"))

summarize_directory <- function(dataset_name) {
  input_dir <- file.path("results", "benchmark", dataset_name)
  output_dir <- file.path("results", "benchmark", "summary")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  files <- files[!grepl("_raw\\.rds$", files)]
  if (!length(files)) {
    message("No results found for ", dataset_name, ".")
    return(invisible(NULL))
  }

  all_results <- lapply(files, readRDS)
  ok <- vapply(all_results, function(z) identical(z$status, "ok"), logical(1))

  if (any(!ok)) {
    errors <- do.call(rbind, lapply(all_results[!ok], function(z) {
      data.frame(
        dataset = dataset_name,
        method = if (!is.null(z$method)) z$method else NA_character_,
        error = if (!is.null(z$error)) z$error else NA_character_,
        warnings = paste(z$warnings, collapse = " | "),
        stringsAsFactors = FALSE
      )
    }))
    utils::write.csv(
      errors,
      file.path(output_dir, paste0(dataset_name, "_errors.csv")),
      row.names = FALSE
    )
  }

  x <- all_results[ok]
  if (!length(x)) return(invisible(NULL))

  method_table <- do.call(rbind, lapply(x, function(z) {
    f <- z$fit
    s <- f$posterior_summary
    d <- f$diagnostics
    native <- native_output_summary(f)

    mean_acceptance <- NA_real_
    if (!is.null(d$post_warmup_acceptance_rate)) {
      mean_acceptance <- mean(d$post_warmup_acceptance_rate, na.rm = TRUE)
    } else if (!is.null(d$mean_accept)) {
      mean_acceptance <- d$mean_accept
    }

    data.frame(
      dataset = dataset_name,
      method = z$method,
      profile = z$profile,
      fit_seconds = if (!is.null(f$fit_seconds)) f$fit_seconds else NA_real_,
      optimize_seconds = if (!is.null(f$optimize_seconds)) f$optimize_seconds else NA_real_,
      sampling_seconds = if (!is.null(f$sampling_seconds)) f$sampling_seconds else NA_real_,
      total_seconds = if (!is.null(z$total_seconds)) z$total_seconds else NA_real_,
      loglik = if (!is.null(f$loglik)) {
        f$loglik
      } else if (!is.null(f$map_loglik)) {
        f$map_loglik
      } else {
        NA_real_
      },
      logpost = if (!is.null(f$map_logpost)) f$map_logpost else NA_real_,
      convergence = if (!is.null(f$convergence)) f$convergence else NA_integer_,
      iterations = if (!is.null(f$iterations)) f$iterations else NA_integer_,
      max_rhat = if (!is.null(s)) safe_max(s$rhat) else NA_real_,
      min_bulk_ess = if (!is.null(s)) safe_min(s$ess_bulk) else NA_real_,
      min_tail_ess = if (!is.null(s)) safe_min(s$ess_tail) else NA_real_,
      min_bulk_ess_per_second = if (!is.null(s) && !is.null(f$fit_seconds)) {
        safe_min(s$ess_bulk) / f$fit_seconds
      } else {
        NA_real_
      },
      min_tail_ess_per_second = if (!is.null(s) && !is.null(f$fit_seconds)) {
        safe_min(s$ess_tail) / f$fit_seconds
      } else {
        NA_real_
      },
      max_mcse_over_sd = if (!is.null(s)) safe_max(s$mcse_mean / s$sd) else NA_real_,
      divergences = if (!is.null(d$divergences)) d$divergences else NA_real_,
      treedepth_hits = if (!is.null(d$treedepth_hits)) d$treedepth_hits else NA_real_,
      min_ebfmi = if (!is.null(d$min_ebfmi)) d$min_ebfmi else NA_real_,
      mean_acceptance = mean_acceptance,
      mean_leapfrog = if (!is.null(d$mean_leapfrog)) d$mean_leapfrog else NA_real_,
      mean_stepsize = if (!is.null(d$mean_stepsize)) d$mean_stepsize else NA_real_,
      native_quantity = native$native_quantity,
      native_mean = native$native_mean,
      native_sd = native$native_sd,
      native_min = native$native_min,
      native_max = native$native_max,
      warnings = paste(z$warnings, collapse = " | "),
      raw_path = if (!is.null(z$raw_path)) z$raw_path else NA_character_,
      stringsAsFactors = FALSE
    )
  }))

  utils::write.csv(
    method_table,
    file.path(output_dir, paste0(dataset_name, "_methods.csv")),
    row.names = FALSE
  )

  rate_rows <- list()
  k <- 1L

  for (z in x) {
    f <- z$fit

    if (!is.null(f$rates)) {
      parameters <- names(f$rates)
      ci <- f$ci
      if (!is.null(ci)) {
        ci <- ci[match(parameters, ci$parameter), , drop = FALSE]
      }

      rate_rows[[k]] <- data.frame(
        dataset = dataset_name,
        method = z$method,
        parameter = parameters,
        estimate = as.numeric(f$rates),
        ml_lower = if (!is.null(ci)) ci$lower else NA_real_,
        ml_upper = if (!is.null(ci)) ci$upper else NA_real_,
        ci_status = if (!is.null(ci)) ci$status else NA_character_,
        posterior_mean = NA_real_,
        posterior_median = NA_real_,
        posterior_sd = NA_real_,
        q025 = NA_real_,
        q975 = NA_real_,
        rhat = NA_real_,
        ess_bulk = NA_real_,
        ess_tail = NA_real_,
        mcse_mean = NA_real_,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }

    if (!is.null(f$posterior_summary)) {
      s <- f$posterior_summary
      rate_rows[[k]] <- data.frame(
        dataset = dataset_name,
        method = z$method,
        parameter = s$parameter,
        estimate = NA_real_,
        ml_lower = NA_real_,
        ml_upper = NA_real_,
        ci_status = NA_character_,
        posterior_mean = s$mean,
        posterior_median = s$median,
        posterior_sd = s$sd,
        q025 = s$q025,
        q975 = s$q975,
        rhat = s$rhat,
        ess_bulk = s$ess_bulk,
        ess_tail = s$ess_tail,
        mcse_mean = s$mcse_mean,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  if (length(rate_rows)) {
    rates <- do.call(rbind, rate_rows)
    utils::write.csv(
      rates,
      file.path(output_dir, paste0(dataset_name, "_rates.csv")),
      row.names = FALSE
    )
  }

  ml <- x[vapply(x, function(z) !is.null(z$fit$rates), logical(1))]
  bayes <- x[vapply(
    x,
    function(z) !is.null(z$fit$posterior_summary) && !grepl("rj", z$method, ignore.case = TRUE),
    logical(1)
  )]
  pairwise <- list()
  p <- 1L

  if (length(ml) >= 2L) {
    for (i in seq_len(length(ml) - 1L)) {
      for (j in (i + 1L):length(ml)) {
        pairwise[[p]] <- data.frame(
          dataset = dataset_name,
          comparison = "ML",
          method_a = ml[[i]]$method,
          method_b = ml[[j]]$method,
          delta_loglik = ml[[i]]$fit$loglik - ml[[j]]$fit$loglik,
          rate_mean_rmse = rate_rmse(ml[[i]]$fit$rates, ml[[j]]$fit$rates),
          rate_sd_rmse = NA_real_,
          stringsAsFactors = FALSE
        )
        p <- p + 1L
      }
    }
  }

  if (length(bayes) >= 2L) {
    for (i in seq_len(length(bayes) - 1L)) {
      for (j in (i + 1L):length(bayes)) {
        a <- bayes[[i]]$fit$posterior_summary
        b <- bayes[[j]]$fit$posterior_summary
        mean_a <- setNames(a$mean, a$parameter)
        mean_b <- setNames(b$mean, b$parameter)
        sd_a <- setNames(a$sd, a$parameter)
        sd_b <- setNames(b$sd, b$parameter)

        pairwise[[p]] <- data.frame(
          dataset = dataset_name,
          comparison = "Bayesian",
          method_a = bayes[[i]]$method,
          method_b = bayes[[j]]$method,
          delta_loglik = NA_real_,
          rate_mean_rmse = rate_rmse(mean_a, mean_b),
          rate_sd_rmse = rate_rmse(sd_a, sd_b),
          stringsAsFactors = FALSE
        )
        p <- p + 1L
      }
    }
  }

  if (length(pairwise)) {
    utils::write.csv(
      do.call(rbind, pairwise),
      file.path(output_dir, paste0(dataset_name, "_pairwise.csv")),
      row.names = FALSE
    )
  }

  if (dataset_name == "rayfish") {
    effect <- x[vapply(x, function(z) identical(z$method, "fish_formula_hmc"), logical(1))]
    if (length(effect) && !is.null(effect[[1]]$main_effect)) {
      utils::write.csv(
        effect[[1]]$main_effect,
        file.path(output_dir, "rayfish_main_effect.csv"),
        row.names = FALSE
      )
    }
  }

  invisible(list(methods = method_table, results = x))
}

summarize_directory("primates")
if (dir.exists(file.path("results", "benchmark", "rayfish"))) {
  summarize_directory("rayfish")
}
cat("Real-data summaries written.\n")
