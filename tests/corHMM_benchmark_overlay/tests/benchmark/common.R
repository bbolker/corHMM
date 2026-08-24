source(file.path("tests", "benchmark", "config.R"))

if (!file.exists("DESCRIPTION")) {
  stop("Run benchmark scripts from the corHMM repository root.")
}

pkgload::load_all(".", quiet = TRUE)

scale_tree_total_length <- function(tree, total = 1) {
  current <- sum(tree$edge.length)
  if (!is.finite(current) || current <= 0) {
    stop("Tree branch lengths must have a positive finite total.")
  }
  scale_factor <- total / current
  tree$edge.length <- tree$edge.length * scale_factor
  tree <- ape::reorder.phylo(tree, "pruningwise")
  attr(tree, "original_total_branch_length") <- current
  attr(tree, "branch_length_scale_factor") <- scale_factor
  tree
}

prepare_trait_data <- function(tree, data) {
  data <- as.data.frame(data)
  names(data)[1] <- "Species"
  data$Species <- as.character(data$Species)
  data <- data[match(tree$tip.label, data$Species), , drop = FALSE]
  rownames(data) <- NULL

  if (anyNA(data$Species) || !identical(tree$tip.label, data$Species)) {
    stop("Tree tips and trait rows do not match.")
  }

  data
}

prepare_dataset <- function(tree, data) {
  data <- as.data.frame(data)
  names(data)[1] <- "Species"
  data$Species <- as.character(data$Species)

  tree <- ape::drop.tip(
    tree,
    setdiff(tree$tip.label, data$Species)
  )

  data <- data[match(tree$tip.label, data$Species), , drop = FALSE]
  keep <- complete.cases(data)

  if (!all(keep)) {
    tree <- ape::drop.tip(tree, data$Species[!keep])
    data <- data[keep, , drop = FALSE]
    data <- data[match(tree$tip.label, data$Species), , drop = FALSE]
  }

  tree <- scale_tree_total_length(tree, 1)
  data <- prepare_trait_data(tree, data)

  for (j in seq.int(2L, ncol(data))) {
    data[[j]] <- as.integer(as.character(data[[j]]))
  }

  values <- unlist(data[-1], use.names = FALSE)
  if (anyNA(values) || !all(values %in% 0:1)) {
    stop("Benchmark traits must be complete binary variables coded 0/1.")
  }

  list(
    tree = tree,
    data = data,
    trait_names = names(data)[-1],
    ntraits = ncol(data) - 1L,
    ntips = ape::Ntip(tree),
    original_total_branch_length = attr(tree, "original_total_branch_length"),
    scale_factor = attr(tree, "branch_length_scale_factor"),
    scaled_total_branch_length = sum(tree$edge.length)
  )
}

load_primates_dataset <- function() {
  e <- new.env(parent = baseenv())
  data_file <- file.path("data", "primates.rda")

  if (file.exists(data_file)) {
    loaded <- load(data_file, envir = e)
    if ("primates" %in% loaded) {
      x <- e$primates
    } else {
      x <- e[[loaded[1]]]
    }
  } else {
    utils::data("primates", package = "corHMM", envir = e)
    x <- e$primates
  }

  phy_index <- which(vapply(x, inherits, logical(1), "phylo"))[1]
  data_index <- which(vapply(x, is.data.frame, logical(1)))[1]

  phy <- x[[phy_index]]
  dat <- x[[data_index]][, 1:3, drop = FALSE]
  names(dat) <- c("Species", "T1", "T2")

  prepare_dataset(phy, dat)
}

load_rayfish_dataset <- function() {
  tree_file <- file.path("data", "c_tree.rds")
  data_file <- file.path("data", "c_traits.rds")

  if (!file.exists(tree_file) || !file.exists(data_file)) {
    stop("Place c_tree.rds and c_traits.rds in data/ before running rayfish analyses.")
  }

  prepare_dataset(
    readRDS(tree_file),
    readRDS(data_file)
  )
}

load_manifest_dataset <- function(row) {
  prepare_dataset(
    readRDS(as.character(row$tree_path)),
    readRDS(as.character(row$data_path))
  )
}

make_rate_prior <- function(tree, override = NULL) {
  if (!is.null(override)) {
    required <- c("lower", "upper", "mean", "sd")
    if (!all(required %in% names(override))) {
      stop("Prior override must contain lower, upper, mean, and sd.")
    }
    return(override[required])
  }

  lower <- log(1)
  upper <- log(100 * ape::Ntip(tree))

  list(
    lower = lower,
    upper = upper,
    mean = (lower + upper) / 2,
    sd = (upper - lower) / 3
  )
}

full_formula_model <- function(trait_names) {
  lapply(trait_names, function(trait) {
    others <- setdiff(trait_names, trait)
    stats::as.formula(
      paste(
        trait,
        "~",
        if (length(others)) paste(others, collapse = " * ") else "1"
      )
    )
  })
}

edge_specification <- function(trait_names) {
  q <- translate(full_formula_model(trait_names), nstate = 2L)
  q$edge_table[, c(
    "edge_id", "from", "to", "from_label", "to_label",
    "changed_trait", "direction"
  )]
}

transition_names <- function(edge_table) {
  paste(edge_table$from_label, edge_table$to_label, sep = "->")
}

formula_start <- function(q_prep, rate = 0.1) {
  out <- numeric(q_prep$n_par)
  out[grep(":\\(Intercept\\)$", q_prep$par_names)] <- log(rate)
  names(out) <- q_prep$par_names
  out
}

extract_classic_rates <- function(fit, edge_table) {
  i <- match(edge_table$from_label, rownames(fit$solution))
  j <- match(edge_table$to_label, colnames(fit$solution))

  if (anyNA(i) || anyNA(j)) {
    stop("Could not match corHMM state labels to the common edge table.")
  }

  setNames(
    as.numeric(fit$solution[cbind(i, j)]),
    transition_names(edge_table)
  )
}

extract_formula_rates <- function(fit, edge_table) {
  X <- fit$q_prep$rate_X
  X <- X[match(transition_names(edge_table), rownames(X)), , drop = FALSE]
  beta <- fit$solution[, 1]

  setNames(
    exp(drop(X %*% beta)),
    transition_names(edge_table)
  )
}

rate_rmse <- function(x, y) {
  if (is.null(names(x)) || is.null(names(y))) {
    stop("Both rate vectors must be named.")
  }

  if (!setequal(names(x), names(y))) {
    stop("Rate names do not match.")
  }

  x <- x[names(y)]
  sqrt(mean((unname(x) - unname(y))^2))
}

parallel_lapply <- function(X, FUN, cores) {
  if (.Platform$OS.type == "windows" || cores <= 1L) {
    lapply(X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = cores)
  }
}

trim_rate_chains <- function(chains, n = NULL) {
  chains <- lapply(chains, as.matrix)
  available <- min(vapply(chains, nrow, integer(1)))
  if (is.null(n)) n <- available
  n <- min(n, available)
  lapply(chains, tail, n = n)
}

rate_draws_array <- function(chains) {
  chains <- trim_rate_chains(chains)
  n <- nrow(chains[[1]])
  p <- ncol(chains[[1]])

  out <- array(
    NA_real_,
    dim = c(n, length(chains), p),
    dimnames = list(NULL, NULL, colnames(chains[[1]]))
  )

  for (i in seq_along(chains)) out[, i, ] <- chains[[i]]
  posterior::as_draws_array(out)
}

summarize_rate_chains <- function(chains) {
  draws <- rate_draws_array(chains)
  x <- posterior::as_draws_matrix(draws)
  parameters <- colnames(chains[[1]])
  x <- x[, parameters, drop = FALSE]

  basic <- data.frame(
    parameter = parameters,
    mean = colMeans(x),
    median = apply(x, 2, stats::median),
    sd = apply(x, 2, stats::sd),
    q025 = apply(x, 2, stats::quantile, probs = 0.025),
    q975 = apply(x, 2, stats::quantile, probs = 0.975),
    row.names = NULL
  )

  diagnostics <- as.data.frame(
    posterior::summarise_draws(
      draws,
      "rhat", "ess_bulk", "ess_tail", "mcse_mean"
    )
  )
  names(diagnostics)[1] <- "parameter"

  out <- merge(basic, diagnostics, by = "parameter", sort = FALSE)
  out[match(parameters, out$parameter), , drop = FALSE]
}

hmc_diagnostics <- function(stanfit, max_treedepth) {
  x <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)

  ebfmi <- vapply(x, function(z) {
    e <- z[, "energy__"]
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))

  list(
    divergences = sum(vapply(x, function(z) sum(z[, "divergent__"]), numeric(1))),
    treedepth_hits = sum(vapply(
      x,
      function(z) sum(z[, "treedepth__"] >= max_treedepth),
      numeric(1)
    )),
    ebfmi = ebfmi,
    min_ebfmi = min(ebfmi, na.rm = TRUE),
    mean_accept = mean(vapply(x, function(z) mean(z[, "accept_stat__"]), numeric(1))),
    mean_leapfrog = mean(vapply(x, function(z) mean(z[, "n_leapfrog__"]), numeric(1))),
    mean_stepsize = mean(vapply(x, function(z) mean(z[, "stepsize__"]), numeric(1)))
  )
}


sanitize_corhmm_fit <- function(fit) {
  if (inherits(fit, "corhmm") && is.function(fit$devfun)) {
    fit$devfun <- NA
  }
  fit
}

native_output_summary <- function(fit) {
  values <- NULL
  name <- NULL

  if (!is.null(fit$native_loglik)) {
    values <- unlist(fit$native_loglik, use.names = FALSE)
    name <- if (!is.null(fit$native_quantity_name)) fit$native_quantity_name else "native log-likelihood"
  } else if (!is.null(fit$log_target)) {
    values <- unlist(fit$log_target, use.names = FALSE)
    name <- if (!is.null(fit$native_quantity_name)) fit$native_quantity_name else "native log-posterior target"
  }

  if (is.null(values) || !length(values)) {
    return(data.frame(
      native_quantity = NA_character_,
      native_mean = NA_real_,
      native_sd = NA_real_,
      native_min = NA_real_,
      native_max = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  values <- values[is.finite(values)]
  if (!length(values)) {
    return(data.frame(
      native_quantity = name,
      native_mean = NA_real_,
      native_sd = NA_real_,
      native_min = NA_real_,
      native_max = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    native_quantity = name,
    native_mean = mean(values),
    native_sd = stats::sd(values),
    native_min = min(values),
    native_max = max(values),
    stringsAsFactors = FALSE
  )
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

benchmark_metadata <- function() {
  git_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
    error = function(e) NA_character_
  )

  packages <- c(
    "corHMM", "RTMB", "rstan", "tmbstan", "posterior",
    "ramcmc", "rphylo", "btw"
  )

  versions <- setNames(vapply(packages, function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  }, character(1)), packages)

  list(
    time = as.character(Sys.time()),
    host = Sys.info()[["nodename"]],
    R = R.version.string,
    package_versions = versions,
    git_commit = git_commit,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
    slurm_array_job_id = Sys.getenv("SLURM_ARRAY_JOB_ID"),
    slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID"),
    slurm_cpus = Sys.getenv("SLURM_CPUS_PER_TASK")
  )
}

atomic_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp")
  saveRDS(object, tmp)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("Could not move temporary result file into place.")
  invisible(file)
}
