library(testthat)

make_formula_data <- function(n = 12, two_trait = TRUE) {
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


test_that("formula-model detector and corHMM model validation work", {
  expect_false(is_formula_model("ARD"))
  expect_false(is_formula_model("SYM"))
  expect_false(is_formula_model("ER"))
  expect_false(is_formula_model(1))
  expect_false(is_formula_model(TRUE))
  expect_false(is_formula_model(list("ARD")))

  expect_true(is_formula_model(T1 ~ 1))
  expect_true(is_formula_model(list(T1 ~ 1, T2 ~ T1)))

  x <- make_formula_data(8)

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = 1,
      node.states = "none",
      verbose = FALSE
    ),
    "`model` must be either"
  )
})


test_that("formula interface rejects unsupported corHMM options", {
  x <- make_formula_data(8)

  rate.mat <- matrix(1, 4, 4)
  diag(rate.mat) <- 0

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = list(T1 ~ 1, T2 ~ T1),
      rate.mat = rate.mat,
      node.states = "none",
      verbose = FALSE
    ),
    "`rate.mat` cannot be used"
  )

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = list(T1 ~ 1, T2 ~ T1),
      node.states = "marginal",
      verbose = FALSE
    ),
    "node.states = 'none'"
  )

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = list(T1 ~ 1, T2 ~ T1),
      node.states = "none",
      tip.fog = 0.01,
      verbose = FALSE
    ),
    "tip.fog"
  )

  expect_error(
    corHMM(
      phy = x$phy,
      data = x$dat,
      rate.cat = 1,
      model = list(T1 ~ 1, T2 ~ T1),
      node.states = "none",
      lewis.asc.bias = TRUE,
      verbose = FALSE
    ),
    "lewis.asc.bias"
  )
})


test_that("formula model fits and returns coefficient-style output", {
  x <- make_formula_data(12)

  fit <- corHMM(
    phy = x$phy,
    data = x$dat,
    rate.cat = 1,
    model = list(T1 ~ 1, T2 ~ T1),
    node.states = "none",
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_true(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))

  expect_true(is.matrix(fit$solution))
  expect_equal(ncol(fit$solution), 1)
  expect_equal(colnames(fit$solution), "estimate")

  expect_true(any(grepl("^T1_gain", rownames(fit$solution))))
  expect_true(any(grepl("^T1_loss", rownames(fit$solution))))
  expect_true(any(grepl("^T2_gain", rownames(fit$solution))))
  expect_true(any(grepl("^T2_loss", rownames(fit$solution))))

  g <- fit$devfun$gr(fit$devfun$par)
  expect_true(all(is.finite(g)))
  expect_equal(length(g), length(fit$devfun$par))
})


test_that("single formula and old character RTMB paths still work", {
  x1 <- make_formula_data(10, two_trait = FALSE)

  fit_formula <- corHMM(
    phy = x1$phy,
    data = x1$dat,
    rate.cat = 1,
    model = T1 ~ 1,
    node.states = "none",
    verbose = FALSE
  )

  expect_s3_class(fit_formula, "corhmm")
  expect_true(isTRUE(fit_formula$formula_model))
  expect_true(is.finite(fit_formula$loglik))

  x2 <- make_formula_data(10, two_trait = FALSE)

  fit_char <- corHMM(
    phy = x2$phy,
    data = x2$dat,
    rate.cat = 1,
    model = "ER",
    node.states = "none",
    nstarts = 1,
    use_RTMB = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit_char, "corhmm")
  expect_false(isTRUE(fit_char$formula_model))
  expect_true(is.finite(fit_char$loglik))
})


test_that("symm formulas share gain/loss parameters and fit", {
  q_prep <- Q_prep(
    formula_list = list(T1 ~ symm(T2), T2 ~ 1),
    nstate = 2
  )

  expect_true(inherits(q_prep$Q0, "sparseMatrix"))
  expect_equal(q_prep$n_par, 4)

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

  x <- make_formula_data(10)

  fit <- corHMM(
    phy = x$phy,
    data = x$dat,
    rate.cat = 1,
    model = list(T1 ~ symm(T2), T2 ~ 1),
    node.states = "none",
    verbose = FALSE
  )

  expect_s3_class(fit, "corhmm")
  expect_true(isTRUE(fit$formula_model))
  expect_true(is.finite(fit$loglik))
})