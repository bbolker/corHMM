library(testthat)

test_that("horseshoe_scale() reduces to Piironen & Vehtari's Eq. 3.13 with equal sigma", {
  sigma <- rep(0.01, 12)
  for (pz in c(0.25, 0.5, 0.9)) {
    p0 <- (1 - pz) * 12
    expect_equal(c(horseshoe_scale(NULL, pz, sigma = sigma)), p0 / (12 - p0) * 0.01,
                 tolerance = 1e-8)
  }
  ## rates with sigma = Inf are always shrunk, so they bound prop_zero below
  expect_error(horseshoe_scale(NULL, 0.1, sigma = c(1, 1, Inf, Inf)), "must exceed 0.5")
})

test_that("a fixed global scale from horseshoe_scale() gives the target prior proportion", {
  set.seed(20261105)
  sigma <- 10^runif(10, -4, -1)
  tau0 <- horseshoe_scale(NULL, 0.6, sigma = sigma)
  pz <- horseshoe_prop_zero(tau0, sigma, fix_global = TRUE, n = 2e4)
  expect_equal(mean(pz), 0.6, tolerance = 0.01)
})

test_that("horseshoe_sigma() gives likelihood scales on the rate scale", {
  set.seed(20261105)
  phy <- ape::rtree(60)
  dat <- data.frame(sp = phy$tip.label, T1 = rbinom(60, 1, 0.5))
  rm <- matrix(c(NA, 2, 1, NA), 2, dimnames = list(c("0", "1"), c("0", "1")))
  fit <- corHMM(phy, dat, rate.cat = 1, rate.mat = rm, node.states = "none",
                root.p = c(0.5, 0.5), use_RTMB = TRUE, verbose = FALSE, collapse = FALSE,
                return.devfun = TRUE)
  sig <- horseshoe_sigma(fit)
  expect_length(sig, 2)
  expect_true(all(sig > 0))
  ## for an interior estimate, about the delta-method standard error
  x <- fit$args.list$p
  se <- exp(x) * sqrt(diag(solve(fit$devfun()$he(x))))
  ok <- x > log(1e-6) & is.finite(sig)
  expect_true(any(ok))
  expect_equal(sig[ok], unname(se[ok]), tolerance = 0.5)
  expect_error(horseshoe_sigma(corHMM(phy, dat, rate.cat = 1, rate.mat = rm,
                                      node.states = "none", root.p = c(0.5, 0.5),
                                      use_RTMB = TRUE, verbose = FALSE, collapse = FALSE)),
               "return.devfun")
})

test_that("horseshoe_scale_prior() gives the target prior probability of a zero rate", {
  set.seed(20261105)
  phy <- ape::rtree(50)
  L <- sum(phy$edge.length)
  for (fg in c(TRUE, FALSE)) {
    tau0 <- horseshoe_scale_prior(0.4, phy, fix_global = fg)
    expect_equal(attr(tau0, "threshold"), 1 / L)
    sim <- horseshoe_prior_zero(tau0, 1, 1 / L, fix_global = fg, n = 2e5)
    expect_equal(mean(sim), 0.4, tolerance = 0.01)
  }
  ## a slab only trims the upper tail, so it raises tau0 a little
  expect_gt(c(horseshoe_scale_prior(0.4, phy, slab_scale = 5 / L)),
            c(horseshoe_scale_prior(0.4, phy)))
  ## the exact probability below matches direct integration
  f <- function(l) (1 - 2 * pnorm(-0.6 / l)) * (2 / pi) / (1 + l^2)
  bf <- integrate(f, 0, 1, rel.tol = 1e-12)$value + integrate(f, 1, Inf, rel.tol = 1e-12)$value
  expect_equal(corHMM:::hs_prob_below(log(0.6)), bf, tolerance = 1e-8)
  expect_error(horseshoe_scale_prior(0.4), "give `phy`")
})
