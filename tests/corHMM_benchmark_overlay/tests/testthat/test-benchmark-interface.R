test_that("mcmcControl supplies dependent defaults", {
  x <- corHMM:::mcmcControl()
  expect_equal(x$chains, 4)
  expect_equal(x$cores, 4)
  expect_equal(x$iter, 2000)
  expect_equal(x$warmup, 1000)
  expect_equal(x$seed, 427)
  expect_false(x$open_progress)
  expect_equal(x$control, list(adapt_delta = 0.9, max_treedepth = 12))

  y <- corHMM:::mcmcControl(chains = 2, iter = 80)
  expect_equal(y$cores, 2)
  expect_equal(y$warmup, 40)
})

test_that("saturated two-trait formula stores a rate mapping", {
  q <- corHMM:::translate(list(T1 ~ T2, T2 ~ T1))
  expect_equal(q$n_par, 8)
  expect_equal(dim(q$rate_X), c(8, 8))

  beta <- numeric(q$n_par)
  beta[grep(":\\(Intercept\\)$", q$par_names)] <- log(0.1)
  expect_equal(as.numeric(exp(q$rate_X %*% beta)), rep(0.1, 8))
})

test_that("RTMB backend contains one objective construction and no REPORT", {
  symbols <- all.names(body(corHMM:::mkdev.corhmm_rtmb), functions = TRUE)
  expect_equal(sum(symbols == "MakeADFun"), 1)
  expect_false("REPORT" %in% symbols)
  expect_false("TMBdata" %in% symbols)
})

test_that("corHMM exposes explicit prior mean and SD", {
  f <- formals(corHMM)
  expect_true(all(c("prior.mean", "prior.sd") %in% names(f)))
})
