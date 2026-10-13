library(testthat)


make_formula_data <- function(n = 12, two_trait = TRUE) {
  set.seed(20)
  phy <- ape::rtree(n)

  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = n)
  )

  if (two_trait) {
    dat$T2 <- rep(c(0, 0, 1, 1), length.out = n)
  }

  list(phy = phy, dat = dat)
}


test_that("Bayesian options are integrated into corHMM", {
  expect_true(all(
    c("prior", "mcmc") %in% names(formals(corHMM))
  ))

  x <- make_formula_data(8, two_trait = FALSE)

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = T1 ~ 1,
      node.states = "none",
      mcmc = list(
        chains = 1,
        cores = 1,
        iter = 1000,
        warmup = 500
      ),
      verbose = FALSE
    ),
    "HMC requires `prior`"
  )

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = "ER",
      node.states = "none",
      prior = list(),
      verbose = FALSE
    ),
    "Priors currently require `use_RTMB = TRUE`"
  )
})


test_that("formula priors are included in the RTMB objective", {
  x <- make_formula_data(10)

  model <- list(
    T1 ~ T2,
    T2 ~ 1
  )

  q_prep <- translate(
    formula_list = model,
    nstate = 2
  )

  pars <- setNames(
    seq(-2, -1, length.out = q_prep$n_par),
    q_prep$par_names
  )

  prior <- list(
    lb = log(1e-4),
    ub = log(10),
    range = 3,
    coef.sd = 2,
    target = "coefficients"
  )

  root.p <- rep(1 / q_prep$d, q_prep$d)

  fit_lik <- corHMM(
    phy = x$phy,
    data = x$dat,
    rate.cat = 1,
    model = model,
    node.states = "none",
    root.p = root.p,
    p = pars,
    verbose = FALSE
  )

  fit_post <- corHMM(
    phy = x$phy,
    data = x$dat,
    rate.cat = 1,
    model = model,
    node.states = "none",
    root.p = root.p,
    p = pars,
    prior = prior,
    return.devfun = TRUE,
    verbose = FALSE
  )

  intercept <- grepl(
    ":\\(Intercept\\)$",
    q_prep$par_names
  )

  prior.mean <- (prior$lb + prior$ub) / 2
  prior.sd <- (prior$ub - prior$lb) / (2 * prior$range)

  logprior <- sum(
    dnorm(
      pars[intercept],
      prior.mean,
      prior.sd,
      log = TRUE
    )
  ) + sum(
    dnorm(
      pars[!intercept],
      0,
      prior$coef.sd,
      log = TRUE
    )
  )

  expect_equal(
    fit_post$loglik,
    fit_lik$loglik,
    tolerance = 1e-6
  )

  expect_equal(
    fit_post$logpost,
    fit_lik$loglik + logprior,
    tolerance = 1e-6
  )

  expect_true(is.na(fit_post$AIC))
  expect_true(is.na(fit_post$AICc))

  expect_equal(
    as.numeric(fit_post$solution[, 1]),
    as.numeric(pars)
  )

  obj <- fit_post$devfun()

  expect_true(is.function(obj$fn))
  expect_true(is.function(obj$gr))

  expect_equal(
    -obj$fn(pars),
    fit_post$logpost,
    tolerance = 1e-6
  )

  expect_true(all(is.finite(obj$gr(pars))))
})


test_that("parallel HMC works through corHMM", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if(.Platform$OS.type == "windows")

  x <- make_formula_data(8, two_trait = FALSE)

  fit <- suppressWarnings(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = T1 ~ 1,
      node.states = "none",
      root.p = c(0.5, 0.5),
      prior = list(
        lb = log(1e-4),
        ub = log(10),
        range = 3,
        coef.sd = 2
      ),
      mcmc = list(
        chains = 2,
        cores = 2,
        iter = 1000,
        warmup = 500,
        seed = 427,
        open_progress = FALSE,
        refresh = 0,
        control = list(
          adapt_delta = 0.8,
          max_treedepth = 6
        )
      ),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "corhmm")
  expect_true(methods::is(fit$mcmc, "stanfit"))

  sampler <- rstan::get_sampler_params(
    fit$mcmc,
    inc_warmup = FALSE
  )

  expect_length(sampler, 2)

  draws <- as.matrix(fit$mcmc)

  expect_true(is.matrix(draws))
  expect_equal(nrow(draws), 2 * (1000 - 500))
  expect_true(all(is.finite(draws)))
})

test_that("one-parameter HMC warns that init is ignored and still samples", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")

  x <- make_formula_data(8, two_trait = FALSE)

  expect_warning(
    fit <- corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = "ER",
      node.states = "none",
      root.p = c(0.5, 0.5),
      use_RTMB = TRUE,
      prior = list(),
      mcmc = list(chains = 1, iter = 100, seed = 101, refresh = 0),
      verbose = FALSE
    ),
    "ignores `init`"
  )

  expect_equal(nrow(as.matrix(fit$mcmc)), 50)
})
