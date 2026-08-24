source(file.path("tests", "benchmark", "common.R"))

run_cor_ml <- function(dataset, method, compute_ci = FALSE) {
  edge_table <- edge_specification(dataset$trait_names)
  root_p <- rep(1 / 2^dataset$ntraits, 2^dataset$ntraits)
  formula_model <- full_formula_model(dataset$trait_names)

  if (method == "cor_formula") {
    q <- translate(formula_model, nstate = 2L)
    ip <- formula_start(q)
    model <- formula_model
    use_rtmb <- TRUE
  } else {
    ip <- rep(0.1, nrow(edge_table))
    model <- "ARD"
    use_rtmb <- identical(method, "cor_rtmb")
  }

  elapsed <- system.time({
    fit <- corHMM(
      phy = dataset$tree,
      data = dataset$data,
      rate.cat = 1,
      model = model,
      node.states = "none",
      root.p = root_p,
      ip = ip,
      collapse = FALSE,
      use_RTMB = use_rtmb,
      return.devfun = compute_ci,
      verbose = FALSE
    )
  })

  rates <- if (method == "cor_formula") {
    extract_formula_rates(fit, edge_table)
  } else {
    extract_classic_rates(fit, edge_table)
  }

  ci <- if (compute_ci) {
    wald_ci_corhmm(fit, edge_table, method)
  } else {
    NULL
  }

  list(
    fit_seconds = unname(elapsed["elapsed"]),
    optimize_seconds = unname(fit$opt.time["elapsed"]),
    loglik = fit$loglik,
    rates = rates,
    ci = ci,
    convergence = fit$convergence,
    iterations = fit$iterations,
    raw = sanitize_corhmm_fit(fit)
  )
}

wald_ci_corhmm <- function(fit, edge_table, method, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)

  result <- tryCatch({
    if (identical(method, "cor_formula")) {
      theta <- fit$solution[, 1]
      obj <- fit$devfun()
      H <- if (!is.null(obj$he)) obj$he(theta) else numDeriv::hessian(obj$fn, theta)
      V <- solve(H)
      X <- fit$q_prep$rate_X
      X <- X[match(transition_names(edge_table), rownames(X)), , drop = FALSE]
      eta <- drop(X %*% theta)
      V_eta <- X %*% V %*% t(X)
      se <- sqrt(diag(V_eta))
    } else {
      theta <- fit$args.list$p
      if (identical(method, "cor_rtmb")) {
        obj <- fit$devfun()
        H <- if (!is.null(obj$he)) obj$he(theta) else numDeriv::hessian(obj$fn, theta)
      } else {
        H <- numDeriv::hessian(function(x) fit$devfun(p = x), theta)
      }
      V <- solve(H)
      i <- match(edge_table$from_label, rownames(fit$index.mat))
      j <- match(edge_table$to_label, colnames(fit$index.mat))
      index <- fit$index.mat[cbind(i, j)]
      eta <- theta[index]
      se <- sqrt(diag(V)[index])
    }

    data.frame(
      parameter = transition_names(edge_table),
      estimate = exp(eta),
      lower = exp(eta - z * se),
      upper = exp(eta + z * se),
      status = "ok",
      row.names = NULL
    )
  }, error = function(e) {
    data.frame(
      parameter = transition_names(edge_table),
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      status = conditionMessage(e),
      row.names = NULL
    )
  })

  result
}

extract_stan_p <- function(stanfit) {
  x <- rstan::extract(
    stanfit,
    permuted = FALSE,
    inc_warmup = FALSE
  )

  if (length(dim(x)) != 3L) stop("Unexpected tmbstan parameter array.")
  index <- grep("^p\\[", dimnames(x)[[3]])
  if (!length(index)) stop("Could not locate the RTMB parameter vector in stanfit.")
  x[, , index, drop = FALSE]
}

extract_stan_lp <- function(stanfit) {
  x <- rstan::extract(
    stanfit,
    permuted = FALSE,
    inc_warmup = FALSE
  )
  index <- match("lp__", dimnames(x)[[3]])
  if (is.na(index)) stop("Could not locate lp__ in stanfit.")
  x[, , index, drop = FALSE]
}

run_rtmb_hmc <- function(dataset, profile, formula_model = NULL) {
  edge_table <- edge_specification(dataset$trait_names)
  prior <- make_rate_prior(dataset$tree, dataset$prior)
  root_p <- rep(1 / 2^dataset$ntraits, 2^dataset$ntraits)
  formula_fit <- !is.null(formula_model)

  if (formula_fit) {
    q <- translate(formula_model, nstate = 2L)
    ip <- formula_start(q)
    model <- formula_model
  } else {
    ip <- rep(0.1, nrow(edge_table))
    model <- "ARD"
  }

  elapsed <- system.time({
    fit <- corHMM(
      phy = dataset$tree,
      data = dataset$data,
      rate.cat = 1,
      model = model,
      node.states = "none",
      root.p = root_p,
      ip = ip,
      collapse = FALSE,
      use_RTMB = TRUE,
      prior = list(target = "rates"),
      prior.mean = prior$mean,
      prior.sd = prior$sd,
      mcmc = list(
        chains = profile$chains,
        cores = profile$chains,
        iter = profile$hmc_iter,
        warmup = profile$hmc_warmup,
        seed = 427,
        open_progress = FALSE,
        control = list(
          adapt_delta = profile$hmc_adapt_delta,
          max_treedepth = profile$hmc_max_treedepth
        )
      ),
      return.devfun = TRUE,
      verbose = FALSE
    )
  })

  p_array <- extract_stan_p(fit$mcmc)

  if (formula_fit) {
    X <- fit$q_prep$rate_X
    X <- X[match(transition_names(edge_table), rownames(X)), , drop = FALSE]
    rate_chains <- lapply(seq_len(profile$chains), function(chain) {
      beta <- matrix(p_array[, chain, , drop = FALSE], ncol = fit$q_prep$n_par)
      out <- exp(beta %*% t(X))
      colnames(out) <- transition_names(edge_table)
      out
    })
    coefficient_chains <- lapply(seq_len(profile$chains), function(chain) {
      out <- matrix(p_array[, chain, , drop = FALSE], ncol = fit$q_prep$n_par)
      colnames(out) <- fit$q_prep$par_names
      out
    })
  } else {
    i <- match(edge_table$from_label, rownames(fit$index.mat))
    j <- match(edge_table$to_label, colnames(fit$index.mat))
    index <- fit$index.mat[cbind(i, j)]
    rate_chains <- lapply(seq_len(profile$chains), function(chain) {
      out <- exp(matrix(p_array[, chain, index, drop = FALSE], ncol = length(index)))
      colnames(out) <- transition_names(edge_table)
      out
    })
    coefficient_chains <- NULL
  }

  rate_chains <- trim_rate_chains(rate_chains, profile$retained_draws)
  if (!is.null(coefficient_chains)) {
    coefficient_chains <- lapply(coefficient_chains, tail, n = profile$retained_draws)
  }

  lp <- extract_stan_lp(fit$mcmc)
  log_target <- lapply(
    seq_len(profile$chains),
    function(chain) tail(drop(lp[, chain, 1]), profile$retained_draws)
  )

  list(
    fit_seconds = unname(elapsed["elapsed"]),
    optimize_seconds = unname(fit$opt.time["elapsed"]),
    sampling_seconds = max(0, unname(elapsed["elapsed"] - fit$opt.time["elapsed"])),
    map_loglik = fit$loglik,
    map_logpost = fit$logpost,
    rate_chains = rate_chains,
    coefficient_chains = coefficient_chains,
    log_target = log_target,
    native_quantity_name = "lp__ (Stan log-posterior target)",
    posterior_summary = summarize_rate_chains(rate_chains),
    diagnostics = hmc_diagnostics(fit$mcmc, profile$hmc_max_treedepth),
    prior = prior,
    raw = sanitize_corhmm_fit(fit)
  )
}

run_ram_chain <- function(obj, start, seed, warmup, iterations, thin,
                          initial_scale, target = 0.234, gamma = 2 / 3) {
  set.seed(seed)
  d <- length(start)
  S <- diag(initial_scale, d)
  current <- as.numeric(start)
  current_lp <- -obj$fn(current)

  n_save <- floor(iterations / thin)
  draws <- matrix(NA_real_, n_save, d)
  log_target <- numeric(n_save)
  accepted <- 0L
  accepted_post <- 0L
  save_index <- 0L
  total <- warmup + iterations

  for (iteration in seq_len(total)) {
    u <- stats::rnorm(d)
    proposal <- current + drop(S %*% u)
    proposal_lp <- -obj$fn(proposal)

    acceptance_probability <- if (is.finite(proposal_lp)) {
      min(1, exp(proposal_lp - current_lp))
    } else {
      0
    }

    if (stats::runif(1) < acceptance_probability) {
      current <- proposal
      current_lp <- proposal_lp
      accepted <- accepted + 1L
      if (iteration > warmup) accepted_post <- accepted_post + 1L
    }

    if (iteration <= warmup) {
      S <- ramcmc::adapt_S(
        S,
        u,
        acceptance_probability,
        iteration,
        target = target,
        gamma = gamma
      )
    } else if ((iteration - warmup) %% thin == 0L) {
      save_index <- save_index + 1L
      draws[save_index, ] <- current
      log_target[save_index] <- current_lp
    }
  }

  colnames(draws) <- names(start)

  list(
    draws = draws,
    log_target = log_target,
    acceptance_rate = accepted / total,
    post_warmup_acceptance_rate = accepted_post / iterations,
    final_proposal_cholesky = S,
    final_proposal_covariance = S %*% t(S)
  )
}

run_rtmb_mh <- function(dataset, profile) {
  edge_table <- edge_specification(dataset$trait_names)
  prior <- make_rate_prior(dataset$tree, dataset$prior)
  root_p <- rep(1 / 2^dataset$ntraits, 2^dataset$ntraits)

  map_time <- system.time({
    map_fit <- corHMM(
      phy = dataset$tree,
      data = dataset$data,
      rate.cat = 1,
      model = "ARD",
      node.states = "none",
      root.p = root_p,
      ip = rep(0.1, nrow(edge_table)),
      collapse = FALSE,
      use_RTMB = TRUE,
      prior = list(target = "rates"),
      prior.mean = prior$mean,
      prior.sd = prior$sd,
      return.devfun = TRUE,
      verbose = FALSE
    )
  })

  obj <- map_fit$devfun()
  start <- map_fit$args.list$p
  names(start) <- names(obj$par)

  chain_time <- system.time({
    chains <- parallel_lapply(seq_len(profile$chains), function(i) {
      run_ram_chain(
        obj = obj,
        start = start,
        seed = profile$chain_seeds[i],
        warmup = profile$rtmb_mh_warmup,
        iterations = profile$rtmb_mh_iterations,
        thin = profile$rtmb_mh_thin,
        initial_scale = profile$rtmb_mh_initial_scale
      )
    }, cores = profile$chains)
  })

  i <- match(edge_table$from_label, rownames(map_fit$index.mat))
  j <- match(edge_table$to_label, colnames(map_fit$index.mat))
  index <- map_fit$index.mat[cbind(i, j)]

  rate_chains <- lapply(chains, function(chain) {
    out <- exp(chain$draws[, index, drop = FALSE])
    colnames(out) <- transition_names(edge_table)
    out
  })
  rate_chains <- trim_rate_chains(rate_chains, profile$retained_draws)
  log_target <- lapply(
    chains,
    function(chain) tail(chain$log_target, profile$retained_draws)
  )

  list(
    fit_seconds = unname(map_time["elapsed"] + chain_time["elapsed"]),
    optimize_seconds = unname(map_time["elapsed"]),
    sampling_seconds = unname(chain_time["elapsed"]),
    map_loglik = map_fit$loglik,
    map_logpost = map_fit$logpost,
    rate_chains = rate_chains,
    log_target = log_target,
    native_quantity_name = "RTMB log-posterior target",
    posterior_summary = summarize_rate_chains(rate_chains),
    diagnostics = list(
      acceptance_rate = vapply(chains, `[[`, numeric(1), "acceptance_rate"),
      post_warmup_acceptance_rate = vapply(
        chains,
        `[[`,
        numeric(1),
        "post_warmup_acceptance_rate"
      ),
      final_proposal_covariance = lapply(chains, `[[`, "final_proposal_covariance")
    ),
    prior = prior,
    raw = list(map_fit = sanitize_corhmm_fit(map_fit), chains = chains)
  )
}

rphylo_column_names <- function(edge_table) {
  paste0(
    "q(",
    gsub("\\|", "", edge_table$from_label),
    "->",
    gsub("\\|", "", edge_table$to_label),
    ")"
  )
}

run_rphylo_fixed <- function(dataset, profile) {
  if (dataset$ntraits != 2L) stop("rphylo dependent models require two binary traits.")

  edge_table <- edge_specification(dataset$trait_names)
  prior <- make_rate_prior(dataset$tree, dataset$prior)
  trees <- structure(list(dataset$tree), class = "multiPhylo")
  traits <- dataset$data[-1]
  rownames(traits) <- dataset$data$Species

  elapsed <- system.time({
    jobs <- parallel_lapply(seq_len(profile$chains), function(i) {
      chain_time <- system.time({
        fit <- rphylo::run_dep_model(
          trees = trees,
          charfile = traits,
          burnin = profile$rphylo_burnin,
          iterations = profile$rphylo_iterations,
          sample_freq = profile$rphylo_sample_freq,
          prior = list(
            type = "lognormal",
            location = prior$mean,
            scale = prior$sd
          ),
          revjump = FALSE,
          pis = rep(1, 4),
          seed = profile$chain_seeds[i],
          verbose = FALSE
        )
      })
      list(fit = fit, seconds = unname(chain_time["elapsed"]))
    }, cores = profile$chains)
  })

  raw_names <- rphylo_column_names(edge_table)
  rate_chains <- lapply(jobs, function(job) {
    z <- as.matrix(job$fit$rates)
    if (!all(raw_names %in% colnames(z))) {
      stop("Unexpected rphylo rate-column names.")
    }
    z <- z[, raw_names, drop = FALSE]
    colnames(z) <- transition_names(edge_table)
    z
  })
  rate_chains <- trim_rate_chains(rate_chains, profile$retained_draws)
  native_loglik <- lapply(
    jobs,
    function(job) tail(job$fit$log_lh, profile$retained_draws)
  )

  list(
    fit_seconds = unname(elapsed["elapsed"]),
    chain_seconds = vapply(jobs, `[[`, numeric(1), "seconds"),
    rate_chains = rate_chains,
    native_loglik = native_loglik,
    native_quantity_name = "rphylo log_lh",
    posterior_summary = summarize_rate_chains(rate_chains),
    diagnostics = lapply(jobs, function(job) {
      job$fit[grep("accept", names(job$fit), ignore.case = TRUE)]
    }),
    prior = prior,
    raw = lapply(jobs, `[[`, "fit")
  )
}

btw_rate_columns <- function(edge_table) {
  lookup <- c(
    "0|0->0|1" = "q12",
    "0|0->1|0" = "q13",
    "0|1->0|0" = "q21",
    "0|1->1|1" = "q24",
    "1|0->0|0" = "q31",
    "1|0->1|1" = "q34",
    "1|1->0|1" = "q42",
    "1|1->1|0" = "q43"
  )
  unname(lookup[transition_names(edge_table)])
}

run_btw_one <- function(data, tree, commands, tag, remove_files = TRUE) {
  bt_path <- Sys.getenv("BAYESTRAITS_HOME")
  bt_bin <- Sys.getenv("BAYESTRAITS_BIN", unset = "BayesTraitsV5")
  options(bt_path = bt_path, bt_bin = bt_bin)

  windows <- .Platform$OS.type == "windows"
  executable <- file.path(bt_path, if (windows) paste0(bt_bin, ".exe") else bt_bin)
  if (!file.exists(executable)) stop("BayesTraits executable was not found.")

  root <- Sys.getenv("SLURM_TMPDIR", unset = tempdir())
  run_dir <- file.path(root, paste0("btw_", tag, "_", Sys.getpid()))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  if (windows) {
    file.copy(
      file.path(bt_path, "BayesTraitsV3.exe"),
      file.path(run_dir, "BayesTraitsV3.exe"),
      overwrite = TRUE
    )
  }

  old <- setwd(run_dir)
  on.exit({
    setwd(old)
    unlink(run_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  btw::bayestraits(
    data,
    tree,
    commands,
    silent = TRUE,
    remove_files = remove_files
  )
}

run_btw_fixed <- function(dataset, profile, reversible_jump = FALSE) {
  if (dataset$ntraits != 2L) stop("BayesTraits dependent models require two binary traits.")

  edge_table <- edge_specification(dataset$trait_names)
  prior <- make_rate_prior(dataset$tree, dataset$prior)
  traits <- dataset$data
  traits[-1] <- lapply(traits[-1], as.character)

  elapsed <- system.time({
    jobs <- parallel_lapply(seq_len(profile$chains), function(i) {
      if (reversible_jump) {
        prior_command <- "RevJump exp 10"
        burnin <- profile$rj_burnin
        iterations <- profile$rj_iterations
        sample_freq <- profile$rj_sample_freq
      } else {
        prior_command <- paste(
          "PriorAll lognormal",
          format(prior$mean, digits = 16),
          format(prior$sd, digits = 16)
        )
        burnin <- profile$btw_burnin
        iterations <- profile$btw_iterations
        sample_freq <- profile$btw_sample_freq
      }

      commands <- c(
        "3", "2",
        paste("Iterations", iterations),
        paste("Burnin", burnin),
        paste("Sample", sample_freq),
        prior_command,
        paste("Seed", profile$chain_seeds[i])
      )

      chain_time <- system.time({
        fit <- run_btw_one(
          traits,
          dataset$tree,
          commands,
          tag = paste0(if (reversible_jump) "rj" else "fixed", "_", i)
        )
      })
      list(fit = fit, seconds = unname(chain_time["elapsed"]))
    }, cores = profile$chains)
  })

  columns <- btw_rate_columns(edge_table)
  rate_chains <- lapply(jobs, function(job) {
    z <- job$fit$Log$results
    if (!all(columns %in% names(z))) stop("Unexpected BayesTraits rate-column names.")
    out <- as.matrix(z[, columns, drop = FALSE])
    colnames(out) <- transition_names(edge_table)
    out
  })
  rate_chains <- trim_rate_chains(rate_chains, profile$retained_draws)
  native_loglik <- lapply(
    jobs,
    function(job) tail(job$fit$Log$results$Lh, profile$retained_draws)
  )

  list(
    fit_seconds = unname(elapsed["elapsed"]),
    chain_seconds = vapply(jobs, `[[`, numeric(1), "seconds"),
    rate_chains = rate_chains,
    native_loglik = native_loglik,
    native_quantity_name = "BayesTraits Lh",
    posterior_summary = summarize_rate_chains(rate_chains),
    diagnostics = list(
      schedule = lapply(jobs, function(job) job$fit$Schedule),
      options = lapply(jobs, function(job) job$fit$Log$options)
    ),
    rj = if (reversible_jump) {
      lapply(jobs, function(job) {
        z <- job$fit$Log$results
        z[, intersect(
          c("No Of Parameters", "No Of Zero", "Model string", "Dep / InDep"),
          names(z)
        ), drop = FALSE]
      })
    } else {
      NULL
    },
    prior = if (reversible_jump) list(type = "exponential", mean = 10) else prior,
    raw = lapply(jobs, `[[`, "fit")
  )
}

run_method <- function(method, dataset, profile, compute_ci = FALSE) {
  switch(
    method,
    cor_orig = run_cor_ml(dataset, "cor_orig", compute_ci),
    cor_rtmb = run_cor_ml(dataset, "cor_rtmb", compute_ci),
    cor_formula = run_cor_ml(dataset, "cor_formula", compute_ci),
    rtmb_mh = run_rtmb_mh(dataset, profile),
    rtmb_hmc = run_rtmb_hmc(dataset, profile),
    formula_hmc = run_rtmb_hmc(dataset, profile, full_formula_model(dataset$trait_names)),
    rphylo = run_rphylo_fixed(dataset, profile),
    btw = run_btw_fixed(dataset, profile, reversible_jump = FALSE),
    btw_rj = run_btw_fixed(dataset, profile, reversible_jump = TRUE),
    stop("Unknown benchmark method: ", method)
  )
}

run_formula_ml_model <- function(dataset, formula_model, compute_ci = TRUE) {
  edge_table <- edge_specification(dataset$trait_names)
  root_p <- rep(1 / 2^dataset$ntraits, 2^dataset$ntraits)
  q <- translate(formula_model, nstate = 2L)

  elapsed <- system.time({
    fit <- corHMM(
      phy = dataset$tree,
      data = dataset$data,
      rate.cat = 1,
      model = formula_model,
      node.states = "none",
      root.p = root_p,
      ip = formula_start(q),
      collapse = FALSE,
      use_RTMB = TRUE,
      return.devfun = compute_ci,
      verbose = FALSE
    )
  })

  list(
    fit_seconds = unname(elapsed["elapsed"]),
    optimize_seconds = unname(fit$opt.time["elapsed"]),
    loglik = fit$loglik,
    rates = extract_formula_rates(fit, edge_table),
    ci = if (compute_ci) wald_ci_corhmm(fit, edge_table, "cor_formula") else NULL,
    convergence = fit$convergence,
    iterations = fit$iterations,
    raw = sanitize_corhmm_fit(fit)
  )
}
