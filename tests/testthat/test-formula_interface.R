library(testthat)
devtools::load_all()
test_that("formula-model detector distinguishes character and formula models", {
  expect_false(is_formula_model("ARD"))
  expect_false(is_formula_model("SYM"))
  expect_false(is_formula_model("ER"))

  expect_true(is_formula_model(T1 ~ 1))
  expect_true(is_formula_model(list(T1 ~ 1, T2 ~ T1)))

  expect_false(is_formula_model(1))
  expect_false(is_formula_model(TRUE))
  expect_false(is_formula_model(list("ARD")))
})


test_that("corHMM rejects invalid model types", {
  set.seed(1)
  phy <- ape::rtree(8)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 8),
    T2 = rep(c(0, 0, 1, 1), length.out = 8)
  )

  expect_error(
    corHMM(
      phy = phy,
      data = dat,
      rate.cat = 1,
      model = 1,
      node.states = "none",
      verbose = FALSE
    ),
    "`model` must be either"
  )
})


test_that("formula model rejects rate.mat", {
  set.seed(1)
  phy <- ape::rtree(8)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 8),
    T2 = rep(c(0, 0, 1, 1), length.out = 8)
  )

  rate.mat <- matrix(1, 4, 4)
  diag(rate.mat) <- 0

  expect_error(
    corHMM(
      phy = phy,
      data = dat,
      rate.cat = 1,
      model = list(T1 ~ 1, T2 ~ T1),
      rate.mat = rate.mat,
      node.states = "none",
      verbose = FALSE
    ),
    "`rate.mat` cannot be used"
  )
})


test_that("formula model currently requires node.states = 'none'", {
  set.seed(1)
  phy <- ape::rtree(8)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 8),
    T2 = rep(c(0, 0, 1, 1), length.out = 8)
  )

  expect_error(
    corHMM(
      phy = phy,
      data = dat,
      rate.cat = 1,
      model = list(T1 ~ 1, T2 ~ T1),
      node.states = "marginal",
      verbose = FALSE
    ),
    "node.states = 'none'"
  )
})


test_that("formula model dispatches through the model argument", {
  set.seed(1)
  phy <- ape::rtree(12)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 12),
    T2 = rep(c(0, 0, 1, 1), length.out = 12)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = list(T1 ~ 1, T2 ~ T1),
    node.states = "none",
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_true(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))
  expect_true(is.matrix(fit$solution))
  expect_equal(fit$rate.cat, 1)
})

test_that("character model still uses old corHMM interface", {
  set.seed(1)
  phy <- ape::rtree(8)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 8)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = "ER",
    node.states = "none",
    ip = 0.1,
    nstarts = 0,
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_false(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))
})


test_that("single formula model dispatches correctly", {
  set.seed(2)
  phy <- ape::rtree(10)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 10)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = T1 ~ 1,
    node.states = "none",
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_true(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))
  expect_true(is.matrix(fit$solution))
})

test_that("formula model returns coefficient-style solution", {
  set.seed(4)
  phy <- ape::rtree(12)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 12),
    T2 = rep(c(0, 0, 1, 1), length.out = 12)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = list(T1 ~ 1, T2 ~ T1),
    node.states = "none",
    verbose = FALSE
  )

  expect_true(isTRUE(fit$formula_model))
  expect_true(is.matrix(fit$solution))
  expect_equal(ncol(fit$solution), 1)
  expect_equal(colnames(fit$solution), "estimate")

  expect_true(any(grepl("^T1_gain", rownames(fit$solution))))
  expect_true(any(grepl("^T1_loss", rownames(fit$solution))))
  expect_true(any(grepl("^T2_gain", rownames(fit$solution))))
  expect_true(any(grepl("^T2_loss", rownames(fit$solution))))
})
