library(testthat)

test_that(".parse_formula basics",
          expect_identical(.parse_formula_spec(T1 ~ 1),
                           list(trait = "T1", symmetric = FALSE,
                                gain_formula = T1 ~ 1, 
                                loss_formula = T1 ~ 1)))


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


## FIXME: repeat data generation much less
set.seed(1)
phy <- ape::rtree(8)
dat <- data.frame(
  sp = phy$tip.label,
  T1 = rep(c(0, 1), length.out = 8),
  T2 = rep(c(0, 0, 1, 1), length.out = 8)
)

test_that("formula model rejects rate.mat", {

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
  ## FIXME: is it important that this one is length 12 vs length 8 ?
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

test_that("character model still works with RTMB path", {
  set.seed(11)
  phy <- ape::rtree(10)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 10)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = "ER",
    node.states = "none",
    nstarts = 1,
    use_RTMB = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_false(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))
})

test_that("formula RTMB devfun has finite gradient", {
  set.seed(12)
  phy <- ape::rtree(10)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 10),
    T2 = rep(c(0, 0, 1, 1), length.out = 10)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = list(T1 ~ 1, T2 ~ T1),
    node.states = "none",
    verbose = FALSE
  )

  g <- fit$devfun$gr(fit$devfun$par)

  expect_true(all(is.finite(g)))
  expect_equal(length(g), length(fit$devfun$par))
})


test_that("symm formula shares gain/loss parameters", {
  q_prep <- Q_prep(
    mode = "formula",
    formula_list = list(T1 ~ symm(T2), T2 ~ 1),
    nstate = 2
  )

  expect_true(any(grepl("^T1_symm", q_prep$par_names)))

  expect_false(any(grepl("^T1_gain", q_prep$par_names)))
  expect_false(any(grepl("^T1_loss", q_prep$par_names)))

  expect_identical(
    q_prep$blocks$T1_gain$par_index,
    q_prep$blocks$T1_loss$par_index
  )

  expect_false(identical(
    q_prep$blocks$T2_gain$par_index,
    q_prep$blocks$T2_loss$par_index
  ))

  expect_equal(q_prep$n_par, 4)
})

test_that("symm formula model fits", {
  set.seed(20)
  phy <- ape::rtree(10)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(c(0, 1), length.out = 10),
    T2 = rep(c(0, 0, 1, 1), length.out = 10)
  )

  fit <- corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = list(T1 ~ symm(T2), T2 ~ 1),
    node.states = "none",
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_true(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))
})
