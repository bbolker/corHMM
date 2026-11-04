library(testthat)

## negative log prior of log-rates x, evaluated through an RTMB tape
prior_nll_at <- function(x, prior) {
  ps <- corHMM:::prior_setup(prior)
  obj <- RTMB::MakeADFun(function(p) corHMM:::prior_nll_rates(p$x, ps),
                         list(x = x), silent = TRUE)
  obj$fn(x)
}

## brute-force marginal densities of the BayesTraits hyperpriors (integrating
## over the log hyper-scale, where the integrand is smooth)
p_hyper_gamma_bf <- function(q, A, B) {
  inner <- function(a) sapply(a, function(aa) {
    g <- function(s) dgamma(q, aa, scale = exp(s)) * exp(s)
    lims <- sort(c(log(q) - 40, min(log(q), log(B)), log(B)))
    integrate(g, lims[1], lims[2], rel.tol = 1e-11, subdivisions = 2000)$value +
      if (lims[3] > lims[2])
        integrate(g, lims[2], lims[3], rel.tol = 1e-11, subdivisions = 2000)$value else 0
  })
  integrate(inner, 0, A, rel.tol = 1e-9, subdivisions = 1000)$value / (A * B)
}
p_hyper_exp_bf <- function(q, M) {
  integrate(function(m) dexp(q, 1 / m), 0, M, rel.tol = 1e-10)$value / M
}


test_that("rate-scale prior densities include the log-scale Jacobian", {
  x <- log(c(0.003, 0.4, 7))
  q <- exp(x)
  expect_equal(prior_nll_at(x, list(family = "gamma", shape = 1.5, scale = 2)),
               -sum(dgamma(q, 1.5, scale = 2, log = TRUE) + x))
  expect_equal(prior_nll_at(x, list(family = "exponential", mean = 10)),
               -sum(dexp(q, 1 / 10, log = TRUE) + x))
  expect_equal(prior_nll_at(x, list(family = "uniform", min = 0, max = 100)),
               -sum(dunif(q, 0, 100, log = TRUE) + x))
  expect_equal(prior_nll_at(x, list(family = "lognormal", lb = log(1e-3),
                                    ub = log(10), range = 3)),
               -sum(dnorm(x, (log(1e-3) + log(10)) / 2,
                          (log(10) - log(1e-3)) / 6, log = TRUE)))
  ## outside the uniform's limits the density is (effectively) zero
  expect_gt(prior_nll_at(log(c(1, 200)), list(family = "uniform", min = 0, max = 100)),
            1e100)
})


test_that("hyperprior marginals match numerical integration", {
  skip_on_cran()
  qs <- 10^c(-8, -3, -1, 0, 1, 2)
  lp_g <- -sapply(log(qs), prior_nll_at,
                  prior = list(family = "hyper_gamma", shape = c(0, 10),
                               scale = c(0, 10))) - log(qs)
  bf_g <- log(sapply(qs, p_hyper_gamma_bf, A = 10, B = 10))
  expect_equal(lp_g, bf_g, tolerance = 1e-6)

  lp_e <- -sapply(log(qs), prior_nll_at,
                  prior = list(family = "hyper_exp", mean = c(0, 10))) - log(qs)
  bf_e <- log(sapply(qs, p_hyper_exp_bf, M = 10))
  expect_equal(lp_e, bf_e, tolerance = 1e-6)

  ## a non-zero lower hyper-scale limit
  pr <- list(family = "hyper_exp", mean = c(1, 10))
  bf <- integrate(function(m) dexp(0.5, 1 / m), 1, 10, rel.tol = 1e-10)$value / 9
  expect_equal(-prior_nll_at(log(0.5), pr) - log(0.5), log(bf), tolerance = 1e-6)
})


## brute-force (regularized) half-horseshoe density of u with global scale 1
## and slab scale r, integrating over log lambda
hs_g_bf <- function(u, r = Inf) {
  f <- function(w) {
    l <- exp(w)
    s <- if (is.finite(r)) l * r / sqrt(r^2 + l^2) else l
    2 * dnorm(u, 0, s) * (2 / pi) / (1 + l^2) * l
  }
  br <- sort(unique(c(log(u) - 10, log(u) + 3, if (is.finite(r)) log(r) + c(-3, 3),
                      max(log(u), 0, if (is.finite(r)) log(r) else 0) + 60)))
  sum(sapply(seq_len(length(br) - 1), function(i) {
    integrate(f, br[i], br[i + 1], rel.tol = 1e-12, subdivisions = 5000)$value
  }))
}

test_that("horseshoe densities match numerical integration", {
  skip_on_cran()
  v <- c(-20, -3, 0, 2)
  for (r in c(Inf, 5)) {
    expect_equal(corHMM:::hs_log_g(v, log(r)), log(sapply(exp(v), hs_g_bf, r = r)),
                 tolerance = 1e-8)
  }
  ## fixed global scale: independent rates, p(q) = g(q / tau) / tau
  x <- log(c(1e-5, 0.02, 3))
  for (sl in c(Inf, 0.5)) {
    pr <- list(family = "horseshoe", scale_global = 0.01, slab_scale = sl, fix_global = TRUE)
    bf <- sum(log(sapply(exp(x) / 0.01, hs_g_bf, r = sl / 0.01) / 0.01) + x)
    expect_equal(-prior_nll_at(x, pr), bf, tolerance = 1e-6)
  }
  ## shared global scale ~ C+(0, tau0), integrated over the tabulated range
  ## of log tau (rate bounds 1e-9, 100, widened by 3)
  q <- c(0.001, 0.3)
  tau0 <- 0.01
  lim <- c(log(1e-9) - 3, log(100) + 3)
  f <- function(t) sapply(t, function(tt) {
    tau <- exp(tt)
    prod(sapply(q / tau, hs_g_bf) / tau) * (2 / pi) * tau * tau0 / (tau0^2 + tau^2)
  })
  br <- seq(lim[1], lim[2], length.out = 40)
  num <- sum(sapply(1:39, function(i) integrate(f, br[i], br[i + 1], rel.tol = 1e-10)$value))
  den <- integrate(function(t) (2 / pi) * exp(t) * tau0 / (tau0^2 + exp(2 * t)),
                   lim[1], lim[2], rel.tol = 1e-12)$value
  expect_equal(-prior_nll_at(log(q), list(family = "horseshoe", scale_global = tau0)) -
                 sum(log(q)), log(num / den), tolerance = 1e-6)
})


test_that("priors give the same posterior with and without the formula interface", {
  set.seed(20261031)
  phy <- ape::rtree(30)
  dat <- data.frame(sp = phy$tip.label, T1 = rbinom(30, 1, 0.5),
                    T2 = rbinom(30, 1, 0.5))
  ## states (T1|T2) 0|0, 1|0, 0|1, 1|1; rates 5, 7 = T2 gain, 2, 4 = T2 loss
  rm_full <- getStateMat4Dat(dat)$rate.mat
  rm_restr <- equateStateMatPars(rm_full, list(c(5, 7), c(2, 4)))
  fit <- function(...) {
    corHMM(phy, dat, rate.cat = 1, node.states = "none", root.p = rep(0.25, 4),
           use_RTMB = TRUE, verbose = FALSE, ...)
  }
  priors <- list(
    list(lb = log(1e-3), ub = log(10)),
    list(family = "gamma", shape = 1.5, scale = 2),
    list(family = "exponential", mean = 10),
    list(family = "hyper_gamma", shape = c(0, 10), scale = c(0, 10)),
    list(family = "hyper_exp", mean = c(0, 10)),
    list(family = "horseshoe", scale_global = 0.1),
    list(family = "horseshoe", scale_global = 0.1, slab_scale = 5, fix_global = TRUE))
  for (pr in priors) {
    expect_equal(fit(model = list(T1 ~ T2, T2 ~ T1), prior = pr)$logpost,
                 fit(rate.mat = rm_full, prior = pr, collapse = FALSE)$logpost,
                 tolerance = 1e-6)
    ## T2 ~ 1: each shared T2 rate gets one prior term, as in the native model
    expect_equal(fit(model = list(T1 ~ T2, T2 ~ 1), prior = pr)$logpost,
                 fit(rate.mat = rm_restr, prior = pr, collapse = FALSE)$logpost,
                 tolerance = 1e-6)
  }
})


test_that("native uniform priors bound the rates", {
  set.seed(20261031)
  phy <- ape::rtree(20)
  dat <- data.frame(sp = phy$tip.label, T1 = rbinom(20, 1, 0.5))
  fit <- corHMM(phy, dat, rate.cat = 1, model = "ARD", node.states = "none",
                root.p = c(0.5, 0.5), use_RTMB = TRUE, verbose = FALSE,
                prior = list(family = "uniform", max = 2))
  expect_true(all(fit$solution <= 2 + 1e-8, na.rm = TRUE))
  expect_true(is.finite(fit$logpost))
})


test_that("invalid prior specifications are rejected", {
  set.seed(20261031)
  phy <- ape::rtree(8)
  dat <- data.frame(sp = phy$tip.label, T1 = rep(0:1, 4))
  fit <- function(...) {
    corHMM(phy, dat, rate.cat = 1, node.states = "none", verbose = FALSE, ...)
  }
  expect_error(fit(model = T1 ~ 1, prior = list(family = "uniform", max = 100)),
               "not available for formula models")
  expect_error(fit(model = T1 ~ 1, prior = list(family = "gamma", shape = 1)),
               "prior\\$scale")
  expect_error(fit(model = T1 ~ 1, prior = list(family = "hyper_gamma",
                                                shape = c(10, 0), scale = c(0, 10))),
               "increasing pair")
  expect_error(fit(model = T1 ~ 1, prior = list(family = "exponential", mean = 10,
                                                target = "coefficients")),
               "requires `prior\\$family = 'lognormal'`")
  expect_error(fit(model = T1 ~ 1, prior = list(family = "horseshoe")),
               "prior\\$scale_global")
  expect_error(fit(model = T1 ~ 1, prior = list(family = "horseshoe", scale_global = 1,
                                                slab_scale = -1)),
               "slab_scale")
  expect_error(fit(model = T1 ~ 1, prior = list(family = "cauchy")))
})


test_that("HMC recovers rate-scale priors when the data are uninformative", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if(.Platform$OS.type == "windows")
  set.seed(20261031)
  phy <- ape::rtree(8)
  probs <- c(0.1, 0.5, 0.9)
  mc <- list(chains = 2, cores = 2, iter = 3000, warmup = 1000, seed = 1031,
             refresh = 0, open_progress = FALSE)

  ## formula model, all traits missing: posterior = gamma hyperprior
  dat_f <- data.frame(sp = phy$tip.label, T1 = NA_integer_)
  fit_f <- suppressWarnings(
    corHMM(phy, dat_f, rate.cat = 1, model = T1 ~ 1, node.states = "none",
           root.p = c(0.5, 0.5), verbose = FALSE, mcmc = mc,
           prior = list(family = "hyper_gamma", shape = c(0, 10), scale = c(0, 10))))
  n <- 2e5
  ref <- log(rgamma(n, runif(n, 0, 10), scale = runif(n, 0, 10)))
  draws <- as.matrix(fit_f$mcmc)[, 1]
  expect_equal(unname(quantile(draws, probs)), unname(quantile(ref, probs)),
               tolerance = 0.3, scale = 1)

  ## native model, all tips ambiguous: posterior = U(0, 100) in each rate
  dat_n <- data.frame(sp = phy$tip.label, T1 = rep("0&1", 8))
  rm <- matrix(c(NA, 2, 1, NA), 2, dimnames = list(c("0", "1"), c("0", "1")))
  fit_n <- suppressWarnings(
    corHMM(phy, dat_n, rate.cat = 1, rate.mat = rm, node.states = "none",
           root.p = c(0.5, 0.5), use_RTMB = TRUE, verbose = FALSE, collapse = FALSE,
           mcmc = mc, prior = list(family = "uniform", max = 100)))
  r <- exp(as.matrix(fit_n$mcmc)[, 1:2])
  expect_true(all(r < 100))
  expect_equal(unname(quantile(r[, 1], probs)), qunif(probs, 0, 100),
               tolerance = 5, scale = 1)
  expect_equal(sum(rstan::get_divergent_iterations(fit_n$mcmc)), 0)

  ## native model, uninformative tips: posterior = horseshoe with a shared
  ## global scale (truncated to [1e-9, 100])
  fit_h <- suppressWarnings(
    corHMM(phy, dat_n, rate.cat = 1, rate.mat = rm, node.states = "none",
           root.p = c(0.5, 0.5), use_RTMB = TRUE, verbose = FALSE, collapse = FALSE,
           mcmc = mc, prior = list(family = "horseshoe", scale_global = 0.1)))
  tau <- abs(0.1 * rcauchy(n))
  ref_h <- cbind(abs(rnorm(n, 0, tau * abs(rcauchy(n)))), abs(rnorm(n, 0, tau * abs(rcauchy(n)))))
  keep <- rowSums(ref_h > 1e-9 & ref_h < 100) == 2
  lr <- as.matrix(fit_h$mcmc)[, 1:2]
  for (j in 1:2) {
    expect_equal(unname(quantile(lr[, j], probs)), unname(quantile(log(ref_h[keep, j]), probs)),
                 tolerance = 0.5, scale = 1)
  }
})


test_that("use.bounds keeps native HMC within the optimization bounds", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if(.Platform$OS.type == "windows")
  set.seed(20261031)
  phy <- ape::rtree(8)
  dat_n <- data.frame(sp = phy$tip.label, T1 = rep("0&1", 8))
  rm <- matrix(c(NA, 2, 1, NA), 2, dimnames = list(c("0", "1"), c("0", "1")))
  fit <- function(...) {
    suppressWarnings(
      corHMM(phy, dat_n, rate.cat = 1, rate.mat = rm, node.states = "none",
             root.p = c(0.5, 0.5), use_RTMB = TRUE, verbose = FALSE, collapse = FALSE,
             upper.bound = 5, prior = list(family = "exponential", mean = 10),
             mcmc = list(chains = 1, iter = 1000, seed = 1031, refresh = 0,
                         open_progress = FALSE, ...)))
  }
  ## posterior = prior (Exp(mean 10)), so P(rate > 5) is about 0.6 without bounds
  r_bounded <- exp(as.matrix(fit()$mcmc)[, 1:2])
  r_free <- exp(as.matrix(fit(use.bounds = FALSE)$mcmc)[, 1:2])
  expect_true(all(r_bounded <= 5))
  expect_gt(mean(r_free > 5), 0.3)
})
