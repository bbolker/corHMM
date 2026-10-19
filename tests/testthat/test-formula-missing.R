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

test_that("formula tips with no observed traits leave the likelihood unchanged", {
  skip_if_not_installed("RTMB")
  data(primates, package = "corHMM")
  phy <- primates$tree
  d <- data.frame(species = as.character(primates$trait$Genus_sp),
                  estrus = primates$trait$T1, multimale = primates$trait$T2)
  model <- list(estrus ~ multimale, multimale ~ estrus)
  ## an all-missing tip's partial likelihood is all ones, and P(t) %*% 1 = 1,
  ## so dropping that tip (and merging the branches on either side of its
  ## parent) gives the same likelihood. Use a tip whose parent is not the root:
  ## dropping a child of the root moves the root
  root <- ape::Ntip(phy) + 1L
  tip_par <- phy$edge[match(seq_len(ape::Ntip(phy)), phy$edge[, 2]), 1]
  drop_tip <- phy$tip.label[which(tip_par != root)[1]]
  d_na <- d
  d_na[d_na$species == drop_tip, c("estrus", "multimale")] <- NA
  fit_formula <- function(phy, d, p = NULL) {
    corHMM(phy, d, rate.cat = 1, model = model, node.states = "none",
           root.p = rep(1/4, 4), use_RTMB = TRUE, p = p, verbose = FALSE)
  }
  fit_na <- expect_warning(fit_formula(phy, d_na), "no observed traits")
  fit_drop <- fit_formula(ape::drop.tip(phy, drop_tip),
                          d[d$species != drop_tip, ], p = fit_na$solution[, 1])
  expect_equal(fit_drop$loglik, fit_na$loglik, tolerance = 1e-8)
})

test_that("formula traits must be binary", {
  skip_if_not_installed("RTMB")
  data(primates, package = "corHMM")
  d <- data.frame(species = as.character(primates$trait$Genus_sp),
                  estrus = primates$trait$T1, multimale = primates$trait$T2)
  d$estrus[1] <- 2
  expect_error(corHMM(primates$tree, d, rate.cat = 1,
                      model = list(estrus ~ multimale, multimale ~ estrus),
                      node.states = "none", root.p = rep(1/4, 4),
                      use_RTMB = TRUE, verbose = FALSE),
               "must be binary variables")
})
