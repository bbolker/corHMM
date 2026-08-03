test_that("base and RTMB share one nlminb path", {
  r_dir <- normalizePath(
    testthat::test_path("..", "..", "R"),
    mustWork = TRUE
  )
  r_lines <- unlist(
    lapply(
      list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE),
      readLines,
      warn = FALSE
    ),
    use.names = FALSE
  )

  expect_equal(
    sum(grepl("^\\s*corHMM\\s*<-\\s*function\\b", r_lines)),
    1
  )
  expect_equal(
    sum(grepl("^\\s*mkdev\\.corhmm_rtmb\\s*<-\\s*function\\b", r_lines)),
    1
  )

  cor_symbols <- all.names(body(corHMM), functions = TRUE)
  rtmb_symbols <- all.names(
    body(corHMM:::mkdev.corhmm_rtmb),
    functions = TRUE
  )

  expect_equal(sum(cor_symbols == "mkdev.corhmm_rtmb"), 1)
  expect_equal(sum(cor_symbols == "nlminb"), 1)
  expect_false(any(c("nloptr", "NLOPT_LN_SBPLX") %in% cor_symbols))

  expect_equal(sum(rtmb_symbols == "MakeADFun"), 1)
  expect_equal(sum(rtmb_symbols == "REPORT"), 1)
  expect_false("TMBdata" %in% rtmb_symbols)
})


test_that("mcmcControl supplies and updates chains and iterations", {
  default <- mcmcControl()

  expect_identical(default, list(chains = 4, iter = 2000))

  control <- mcmcControl(chains = 2, iter = 80)

  expect_identical(control, list(chains = 2, iter = 80))

  control$chains <- 1
  control$iter <- 40

  expect_identical(control, list(chains = 1, iter = 40))
})


test_that("classic and RTMB optimization give the same likelihood", {
  skip_if_not_installed("RTMB")

  set.seed(20)
  phy <- ape::rtree(12)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(0:1, 6)
  )

  args <- list(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = "ARD",
    node.states = "none",
    root.p = c(0.5, 0.5),
    ip = 0.1,
    opts = list(
      eval.max = 1000,
      iter.max = 1000,
      rel.tol = 1e-10
    ),
    verbose = FALSE
  )

  fit_base <- do.call(corHMM, c(args, list(use_RTMB = FALSE)))
  fit_rtmb <- do.call(corHMM, c(args, list(use_RTMB = TRUE)))

  expect_s3_class(fit_base, "corhmm")
  expect_s3_class(fit_rtmb, "corhmm")
  expect_true(is.finite(fit_base$loglik))
  expect_true(is.finite(fit_rtmb$loglik))
  expect_equal(fit_rtmb$loglik, fit_base$loglik, tolerance = 1e-6)
})


test_that("formula and symmetric models still use RTMB", {
  skip_if_not_installed("RTMB")

  set.seed(21)
  phy <- ape::rtree(12)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(0:1, 6),
    T2 = rep(c(0, 0, 1, 1), 3)
  )

  fit <- corHMM(
    phy, dat, 1,
    model = list(T1 ~ T2, T2 ~ T1),
    node.states = "none",
    root.p = rep(0.25, 4),
    verbose = FALSE
  )

  q <- corHMM:::translate(list(T1 ~ symm(T2), T2 ~ 1))
  fit_symm <- corHMM(
    phy, dat, 1,
    model = list(T1 ~ symm(T2), T2 ~ 1),
    node.states = "none",
    root.p = rep(0.25, 4),
    p = rep(0, q$n_par),
    verbose = FALSE
  )

  expect_true(fit$formula_model)
  expect_true(is.finite(fit$loglik))
  expect_equal(nrow(fit$solution), 8)
  expect_identical(
    q$blocks$T1_gain$par_index,
    q$blocks$T1_loss$par_index
  )
  expect_equal(nrow(fit_symm$solution), 4)
})


test_that("priors still produce classic and formula MAP fits", {
  skip_if_not_installed("RTMB")

  set.seed(22)
  phy <- ape::rtree(10)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(0:1, 5)
  )
  prior <- list(range = 3, coef.sd = 2)

  fits <- list(
    corHMM(
      phy, dat, 1,
      model = "ARD",
      node.states = "none",
      root.p = c(0.5, 0.5),
      use_RTMB = TRUE,
      prior = prior,
      verbose = FALSE
    ),
    corHMM(
      phy, dat, 1,
      model = T1 ~ 1,
      node.states = "none",
      root.p = c(0.5, 0.5),
      prior = prior,
      verbose = FALSE
    )
  )

  for (fit in fits) {
    expect_true(is.finite(fit$loglik))
    expect_true(is.finite(fit$logpost))
    expect_true(is.na(fit$AIC))
  }
})


test_that("mcmcControl runs classic, formula, and multiple-chain HMC", {
  skip_on_cran()
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  set.seed(23)
  phy <- ape::rtree(8)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(0:1, 4)
  )
  prior <- list(range = 3, coef.sd = 2)

  classic_control <- mcmcControl(
    chains = 1,
    iter = 40
  )

  formula_control <- mcmcControl(
    chains = 2,
    iter = 40
  )

  fit_classic <- suppressWarnings(corHMM(
    phy, dat, 1,
    model = "ARD",
    node.states = "none",
    root.p = c(0.5, 0.5),
    use_RTMB = TRUE,
    prior = prior,
    mcmc = classic_control,
    verbose = FALSE
  ))

  fit_formula <- suppressWarnings(corHMM(
    phy, dat, 1,
    model = T1 ~ 1,
    node.states = "none",
    root.p = c(0.5, 0.5),
    prior = prior,
    mcmc = formula_control,
    verbose = FALSE
  ))

  expect_s4_class(fit_classic$mcmc, "stanfit")
  expect_s4_class(fit_formula$mcmc, "stanfit")
  expect_equal(fit_classic$mcmc@sim$chains, 1)
  expect_equal(fit_formula$mcmc@sim$chains, 2)
})


test_that("btw runs two seeded BayesTraits V5 chains", {
  skip_on_cran()
  skip_if_not_installed("btw")

  bt_path <- Sys.getenv("BAYESTRAITS_HOME")
  bt_bin <- Sys.getenv(
    "BAYESTRAITS_BIN",
    unset = "BayesTraitsV5"
  )

  skip_if(
    bt_path == "",
    "Set BAYESTRAITS_HOME to the BayesTraits directory."
  )

  bt_exe <- file.path(
    bt_path,
    if (.Platform$OS.type == "windows") {
      paste0(bt_bin, ".exe")
    } else {
      bt_bin
    }
  )

  skip_if(
    !file.exists(bt_exe),
    paste("BayesTraits executable not found:", bt_exe)
  )

  if (.Platform$OS.type == "windows") {
    skip_if(
      !file.exists(file.path(bt_path, "BayesTraitsV3.exe")),
      "BayesTraitsV3.exe compatibility alias not found."
    )
  }

  old_options <- options(
    bt_path = bt_path,
    bt_bin = bt_bin
  )
  on.exit(options(old_options), add = TRUE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(bt_path)

  set.seed(24)
  phy <- ape::rtree(12)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = as.character(rep(0:1, 6)),
    T2 = as.character(rep(c(0, 0, 1, 1), 3))
  )

  fits <- lapply(c(427, 428), function(seed) {
    btw::bayestraits(
      dat,
      phy,
      c(
        "3", "2",
        "Iterations 40",
        "Burnin 20",
        "Sample 10",
        "PriorAll exp 10",
        paste("Seed", seed)
      ),
      silent = TRUE,
      remove_files = TRUE
    )
  })

  results <- lapply(fits, function(x) x$Log$results)

  expect_length(fits, 2)
  expect_true(all(vapply(
    results,
    function(x) is.data.frame(x) && nrow(x) > 0,
    logical(1)
  )))
  expect_false(identical(results[[1]], results[[2]]))
})


test_that("RTMB optimization preserves names and reports likelihood", {
  skip_if_not_installed("RTMB")

  set.seed(25)
  phy <- ape::rtree(8)
  dat <- data.frame(
    sp = phy$tip.label,
    T1 = rep(0:1, 4)
  )

  fit <- corHMM(
    phy, dat, 1,
    model = "ARD",
    node.states = "none",
    root.p = c(0.5, 0.5),
    use_RTMB = TRUE,
    prior = list(range = 3, coef.sd = 2),
    return.devfun = TRUE,
    verbose = FALSE
  )
  obj <- fit$devfun()

  expect_identical(
    names(obj$env$last.par.best),
    names(obj$par)
  )
  expect_true(all(nzchar(names(obj$env$last.par.best))))
  expect_length(
    obj$report(obj$env$last.par.best)$likelihood_nll,
    1
  )
})