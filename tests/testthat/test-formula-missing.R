## formula models with missing trait values (NA or "?") treat them as
## unknown, like "?" in the native interface: compare the two objectives at
## the same rates
test_that("formula interface handles missing trait values", {
  skip_if_not_installed("RTMB")
  data(primates, package = "corHMM")
  phy <- primates$tree
  d <- data.frame(species = as.character(primates$trait$Genus_sp),
                  estrus = primates$trait$T1, multimale = primates$trait$T2)
  d$multimale[1] <- NA
  d$estrus[2] <- "?"
  fit_f <- corHMM(phy, d, rate.cat = 1, model = list(estrus ~ multimale, multimale ~ estrus),
                  node.states = "none", root.p = rep(1/4, 4), use_RTMB = TRUE,
                  return.devfun = TRUE, verbose = FALSE)
  d_n <- d
  d_n[] <- lapply(d_n, function(x) { x <- as.character(x); x[is.na(x)] <- "?"; x })
  fit_n <- corHMM(phy, d_n, rate.cat = 1, model = "ARD", node.states = "none",
                  root.p = rep(1/4, 4), use_RTMB = TRUE, collapse = FALSE,
                  return.devfun = TRUE, verbose = FALSE)
  obj_f <- fit_f$devfun()
  obj_n <- fit_n$devfun()
  ## rates at the formula fit, in the native model's parameter order
  X <- fit_f$q_prep$rate_X
  log_rates <- drop(X %*% fit_f$solution[, 1])
  ends <- do.call(rbind, strsplit(rownames(X), "->"))
  im <- fit_n$index.mat
  p_n <- numeric(max(im, na.rm = TRUE))
  p_n[im[ends]] <- log_rates
  expect_equal(obj_f$fn(fit_f$solution[, 1]), obj_n$fn(p_n), tolerance = 1e-8)
  ## a tip with no observed traits contributes nothing (likelihood 1)
  d2 <- d
  d2[3, c("estrus", "multimale")] <- NA
  expect_warning(corHMM(phy, d2, rate.cat = 1, model = list(estrus ~ multimale, multimale ~ estrus),
                        node.states = "none", root.p = rep(1/4, 4), use_RTMB = TRUE,
                        verbose = FALSE), "no observed traits")
})
