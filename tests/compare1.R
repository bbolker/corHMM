library(ape)
library(rphylo)
library(btw)
library(coda)

pkgload::load_all(".", quiet = TRUE)

## Settings -----------------------------------------------------------------

n_chains <- 4L
chain_seeds <- 427:430
n_keep <- 2000L

mh_burnin <- 100000L
mh_sample <- 500L
mh_post_iter <- n_keep * mh_sample
bt_total_iter <- mh_burnin + mh_post_iter

hmc_iter <- 2L * n_keep
hmc_warmup <- hmc_iter / 2L
hmc_cores <- min(n_chains, parallel::detectCores())

## A common broad Lognormal prior on rates.
prior_mu <- log(0.1)
prior_sd <- 3
prior_range <- 3
lower_rate <- exp(prior_mu - prior_range * prior_sd)
upper_rate <- exp(prior_mu + prior_range * prior_sd)

## Formula slopes are log-rate contrasts. This matches their marginal SD
## under independent N(prior_mu, prior_sd^2) log-rate priors.
cor_prior <- list(
  range = prior_range,
  coef.sd = sqrt(2) * prior_sd
)

options(mc.cores = hmc_cores)
rstan::rstan_options(auto_write = TRUE)

out_dir <- file.path("results", "primates_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Data ---------------------------------------------------------------------

data("primates", package = "corHMM")

phy <- primates[[which(vapply(primates, inherits, logical(1), "phylo"))[1]]]
dat <- primates[[which(vapply(primates, is.data.frame, logical(1)))[1]]]
dat <- dat[, 1:3]
names(dat) <- c("species", "T1", "T2")

dat$species <- as.character(dat$species)
dat$T1 <- as.integer(as.character(dat$T1))
dat$T2 <- as.integer(as.character(dat$T2))
dat <- dat[complete.cases(dat), ]

phy <- drop.tip(phy, setdiff(phy$tip.label, dat$species))
dat <- dat[match(phy$tip.label, dat$species), ]

stopifnot(
  identical(dat$species, phy$tip.label),
  all(dat$T1 %in% 0:1),
  all(dat$T2 %in% 0:1)
)

cor_data <- dat

rphylo_trees <- structure(list(phy), class = "multiPhylo")
rphylo_traits <- data.frame(
  T1 = dat$T1,
  T2 = dat$T2,
  row.names = dat$species
)

btw_traits <- data.frame(
  species = dat$species,
  T1 = as.character(dat$T1),
  T2 = as.character(dat$T2)
)

## rphylo: four sequential fixed-dimension chains --------------------------

fit_rphylo <- vector("list", n_chains)
rphylo_seconds <- numeric(n_chains)

for (i in seq_len(n_chains)) {
  tm <- system.time({
    fit_rphylo[[i]] <- rphylo::run_dep_model(
      trees = rphylo_trees,
      charfile = rphylo_traits,
      burnin = mh_burnin,
      iterations = mh_post_iter,
      sample_freq = mh_sample,
      prior = list(
        type = "lognormal",
        location = prior_mu,
        scale = prior_sd
      ),
      revjump = FALSE,
      pis = rep(1, 4),
      seed = chain_seeds[i],
      verbose = FALSE
    )
  })

  rphylo_seconds[i] <- unname(tm["elapsed"])
  saveRDS(
    fit_rphylo[[i]],
    file.path(out_dir, paste0("rphylo_chain_", i, ".rds"))
  )
}

## BayesTraits/btw: four sequential single-core chains ---------------------

bt_path <- Sys.getenv("BAYESTRAITS_HOME")
bt_bin <- Sys.getenv("BAYESTRAITS_BIN", "BayesTraitsV5")
options(bt_path = bt_path, bt_bin = bt_bin)

bt_exe <- file.path(
  bt_path,
  if (.Platform$OS.type == "windows") paste0(bt_bin, ".exe") else bt_bin
)
stopifnot(file.exists(bt_exe))

run_btw <- function(seed) {
  old_wd <- setwd(bt_path)
  on.exit(setwd(old_wd))

  btw::bayestraits(
    btw_traits,
    phy,
    c(
      "3",
      "2",
      paste("Iterations", bt_total_iter),
      paste("Burnin", mh_burnin),
      paste("Sample", mh_sample),
      "Pis none",
      paste("PriorAll lognormal", prior_mu, prior_sd),
      paste("Seed", seed)
    ),
    silent = TRUE,
    remove_files = TRUE
  )
}

fit_btw <- vector("list", n_chains)
btw_seconds <- numeric(n_chains)

for (i in seq_len(n_chains)) {
  tm <- system.time({
    fit_btw[[i]] <- run_btw(chain_seeds[i])
  })

  btw_seconds[i] <- unname(tm["elapsed"])
  saveRDS(
    fit_btw[[i]],
    file.path(out_dir, paste0("btw_chain_", i, ".rds"))
  )
}

## corHMM RTMB-HMC: four chains in parallel --------------------------------

set.seed(chain_seeds[1])
tm <- system.time({
  fit_hmc <- corHMM(
    phy = phy,
    data = cor_data,
    rate.cat = 1,
    model = "ARD",
    node.states = "none",
    root.p = rep(0.25, 4),
    ip = rep(0.1, 8),
    collapse = FALSE,
    lower.bound = lower_rate,
    upper.bound = upper_rate,
    use_RTMB = TRUE,
    prior = cor_prior,
    mcmc = list(chains = n_chains, iter = hmc_iter),
    verbose = FALSE
  )
})
hmc_seconds <- unname(tm["elapsed"])
saveRDS(fit_hmc, file.path(out_dir, "corhmm_rtmb_hmc.rds"))

set.seed(chain_seeds[1])
tm <- system.time({
  fit_formula <- corHMM(
    phy = phy,
    data = cor_data,
    rate.cat = 1,
    model = list(T1 ~ T2, T2 ~ T1),
    node.states = "none",
    root.p = rep(0.25, 4),
    ip = rep(c(log(0.1), 0), 4),
    collapse = FALSE,
    lower.bound = lower_rate,
    upper.bound = upper_rate,
    use_RTMB = TRUE,
    prior = cor_prior,
    mcmc = list(chains = n_chains, iter = hmc_iter),
    verbose = FALSE
  )
})
formula_seconds <- unname(tm["elapsed"])
saveRDS(fit_formula, file.path(out_dir, "corhmm_formula_hmc.rds"))

## Put every posterior into the same eight-rate order -----------------------

rate_names <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
from <- c("0|0", "0|0", "0|1", "0|1", "1|0", "1|0", "1|1", "1|1")
to <- c("0|1", "1|0", "0|0", "1|1", "0|0", "1|1", "0|1", "1|0")
transition_keys <- paste(from, to, sep = "->")

rphylo_draws <- lapply(fit_rphylo, function(x) {
  out <- as.matrix(x$rates)
  colnames(out) <- rate_names
  out
})

btw_draws <- lapply(fit_btw, function(x) {
  as.matrix(x$Log$results[, rate_names, drop = FALSE])
})

classic_index <- fit_hmc$index.mat[cbind(
  match(from, rownames(fit_hmc$index.mat)),
  match(to, colnames(fit_hmc$index.mat))
)]

hmc_array <- rstan::extract(
  fit_hmc$mcmc,
  permuted = FALSE,
  inc_warmup = FALSE
)

hmc_draws <- lapply(seq_len(n_chains), function(chain) {
  out <- matrix(
    hmc_array[, chain, classic_index, drop = FALSE],
    ncol = length(rate_names)
  )
  colnames(out) <- rate_names
  exp(out)
})

q <- fit_formula$q_prep
formula_design <- matrix(
  0,
  nrow = length(rate_names),
  ncol = q$n_par,
  dimnames = list(transition_keys, q$par_names)
)

for (block in q$blocks) {
  rows <- match(
    paste(block$from_label, block$to_label, sep = "->"),
    transition_keys
  )
  formula_design[rows, block$par_index] <- block$X
}

formula_array <- rstan::extract(
  fit_formula$mcmc,
  permuted = FALSE,
  inc_warmup = FALSE
)

formula_draws <- lapply(seq_len(n_chains), function(chain) {
  beta <- matrix(
    formula_array[, chain, seq_len(q$n_par), drop = FALSE],
    ncol = q$n_par
  )
  out <- exp(beta %*% t(formula_design))
  colnames(out) <- rate_names
  out
})

draws <- list(
  rphylo = rphylo_draws,
  `BayesTraits/btw` = btw_draws,
  `corHMM RTMB-HMC` = hmc_draws,
  `corHMM formula RTMB-HMC` = formula_draws
)

## R-hat, ESS, and MCSE of each posterior mean ------------------------------

diagnose <- function(chains) {
  mlist <- coda::mcmc.list(lapply(chains, coda::mcmc))
  rhat <- coda::gelman.diag(
    mlist,
    multivariate = FALSE,
    autoburnin = FALSE
  )$psrf[, "Point est."]
  ess <- coda::effectiveSize(mlist)
  pooled <- do.call(rbind, chains)
  post_sd <- apply(pooled, 2, sd)

  data.frame(
    parameter = colnames(pooled),
    mean = colMeans(pooled),
    sd = post_sd,
    rhat = unname(rhat[colnames(pooled)]),
    ess = unname(ess[colnames(pooled)]),
    mcse_mean = unname(post_sd / sqrt(ess[colnames(pooled)]))
  )
}

diagnostics <- do.call(
  rbind,
  lapply(names(draws), function(method) {
    cbind(method = method, diagnose(draws[[method]]))
  })
)
rownames(diagnostics) <- NULL

## Pairwise posterior-mean RMSE and Gaussian KL in log-rate space -----------

gaussian_kl <- function(x, y) {
  x <- log(pmax(x, 1e-12))
  y <- log(pmax(y, 1e-12))
  p <- ncol(x)
  sx <- cov(x) + diag(1e-8, p)
  sy <- cov(y) + diag(1e-8, p)
  d <- colMeans(y) - colMeans(x)

  out <- 0.5 * (
    sum(diag(solve(sy, sx))) +
      drop(crossprod(d, solve(sy, d))) - p +
      as.numeric(determinant(sy, logarithm = TRUE)$modulus) -
      as.numeric(determinant(sx, logarithm = TRUE)$modulus)
  )
  max(0, out)
}

pooled <- lapply(draws, function(x) do.call(rbind, x))
methods <- names(pooled)

rmse_matrix <- kl_matrix <- matrix(
  0,
  length(methods),
  length(methods),
  dimnames = list(methods, methods)
)

for (i in seq_along(methods)) {
  for (j in seq_along(methods)) {
    rmse_matrix[i, j] <- sqrt(mean(
      (colMeans(pooled[[i]]) - colMeans(pooled[[j]]))^2
    ))
    kl_matrix[i, j] <- gaussian_kl(pooled[[i]], pooled[[j]])
  }
}

## Main comparison ----------------------------------------------------------

timing <- data.frame(
  method = methods,
  chains = n_chains,
  cores_used = c(1L, 1L, hmc_cores, hmc_cores),
  parallel_chains = c(1L, 1L, min(n_chains, hmc_cores), min(n_chains, hmc_cores)),
  total_iterations_per_chain = c(
    mh_burnin + mh_post_iter,
    bt_total_iter,
    hmc_iter,
    hmc_iter
  ),
  burnin_or_warmup = c(mh_burnin, mh_burnin, hmc_warmup, hmc_warmup),
  retained_per_chain = vapply(draws, function(x) nrow(x[[1]]), integer(1)),
  wall_seconds = c(
    sum(rphylo_seconds),
    sum(btw_seconds),
    hmc_seconds,
    formula_seconds
  )
)

comparison <- transform(
  timing,
  wall_minutes = wall_seconds / 60,
  rmse_to_btw = rmse_matrix[method, "BayesTraits/btw"],
  max_rhat = vapply(method, function(x) {
    max(diagnostics$rhat[diagnostics$method == x], na.rm = TRUE)
  }, numeric(1)),
  min_ess = vapply(method, function(x) {
    min(diagnostics$ess[diagnostics$method == x], na.rm = TRUE)
  }, numeric(1)),
  max_mcse_mean = vapply(method, function(x) {
    max(diagnostics$mcse_mean[diagnostics$method == x], na.rm = TRUE)
  }, numeric(1))
)

hardware <- data.frame(
  os = paste(Sys.info()[c("sysname", "release")], collapse = " "),
  machine = Sys.info()[["machine"]],
  logical_cores = parallel::detectCores(logical = TRUE),
  physical_cores = parallel::detectCores(logical = FALSE),
  hmc_cores = hmc_cores
)

write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)
write.csv(rmse_matrix, file.path(out_dir, "rmse_matrix.csv"))
write.csv(kl_matrix, file.path(out_dir, "kl_matrix.csv"))
write.csv(hardware, file.path(out_dir, "hardware.csv"), row.names = FALSE)

saveRDS(
  list(
    comparison = comparison,
    diagnostics = diagnostics,
    rmse_matrix = rmse_matrix,
    kl_matrix = kl_matrix,
    hardware = hardware,
    settings = list(
      chains = n_chains,
      chain_seeds = chain_seeds,
      target_retained_per_chain = n_keep,
      prior_mu = prior_mu,
      prior_sd = prior_sd,
      formula_prior_note = paste(
        "The formula parameterization matches intercept and contrast",
        "marginals but not the covariance of independent log-rate priors."
      )
    )
  ),
  file.path(out_dir, "comparison_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

print(comparison)
print(round(rmse_matrix, 4))
print(round(kl_matrix, 3))
