## https://github.com/bbolker/corHMM/blob/bolker_clean/misc/test_RTMB.R
remotes::install_github("bbolker/corHMM", ref = "bolker_clean")
## source("~/R/pkgs/corHMM/R/RTMB_devcorhmm.R")
library(corHMM)
library(rbenchmark)
data("primates")
phy <- reorder(primates$tree, "pruningwise")

## clunky, have to find the right place in corHMM (inside random.restarts())
##  to save necessary bits ...
## trace(corHMM, at=46, browser)
## fit1 <- corHMM(phy = phy, data = primates$trait, rate.cat = 1)
## save("model.set.final", "starts", file = "model.set.final.rda")

cfun  <- function(...) {
  invisible(capture.output(x <- corHMM(...)))
  x
}
parfun <- function(x) log(na.omit(c(x$solution)))


fit_orig <- cfun(phy = phy, data = primates$trait, rate.cat = 1)
print(fit_orig)

fit_RTMB <- cfun(phy = phy, data = primates$trait, rate.cat = 1, use_RTMB = TRUE)

stopifnot(all.equal(fit_orig$loglik, fit_RTMB$loglik, tolerance = 2e-8)) ## mean diff: 1.6e-8
stopifnot(all.equal(parfun(fit_RTMB), parfun(fit_orig),
                    tolerance = 2e-4))
## Mean relative difference: 0.000103487

## takes about 1-2 minutes
bb <- benchmark(
  cfun(phy = phy, data = primates$trait, rate.cat = 1),
  cfun(phy = phy, data = primates$trait, rate.cat = 1, use_RTMB = TRUE),
  replications = 20,
  columns = c("test", "elapsed", "relative"))

## a slightly unfair comparison because fit_orig() also does some model setup tasks.
## but the bulk of the time is in the optimization ...

##        test elapsed relative
## 2 fit_new()   0.804    1.000
## 1 fit_orig()  50.929   63.345


###
library(ape)
library(phangorn)
nstate <- 2; ntrait <- 2; ntaxa <- 20
set.seed(7)
m3 <- ape::rtree(ntaxa)
g1 <- reorder(m3, "pruningwise")

Q <- matrix(
  c(0, 0.017, 0.124, 0, 0.234,
    0, 0, 0.26, 0.052,
    0, 0, 0.108, 0,
    0.136, 0.107, 0),
  nrow = 4L,
  ncol = 4L,
  dimnames = list(c("0|0", "1|0", "0|1", "1|1"), c("0|0", "1|0", "0|1", "1|1"))
)

traitMatrix <- data.frame(
  nm = c(
    "t15", "t8", "t12", "t19", "t2", "t20", "t6", "t14", "t11", "t5", "t4", "t13",
    "t16", "t7", "t17", "t10", "t9", "t3", "t18", "t1"
  ),
  V1 = rep(c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L), c(2L, 6L, 1L, 1L, 6L, 1L, 1L, 1L, 1L)),
  V2 = rep(c(1L, 0L, 1L), c(3L, 4L, 13L)),
  row.names = c(
    "t15", "t8", "t12", "t19", "t2", "t20", "t6", "t14", "t11", "t5", "t4", "t13",
    "t16", "t7", "t17", "t10", "t9", "t3", "t18", "t1"
  )
)

invisible(capture.output(cc <- corHMM(g1, traitMatrix, rate.cat = 1)))
expect_equal(c(log(na.omit(c(cc$solution))))
             c(-1.96837728100874, 1.901883681058, -2.28196888984159, 1.1263333445407))
expect_equal(cc$loglik, -15.0219701261726)

## simulate traits
## s <- phangorn::simSeq(g1, l = 1, Q = Q,
##        type = "USER", levels = seq(nstate^ntrait),
##        rate = 1)

## ## convert from enumerated traits (0-7) to 3 binary digits
## ## https://stackoverflow.com/questions/6614283/converting-decimal-to-binary-in-r
## to_bin <- function(x, n = 2) {
##   intToBits(x) |> rev() |> as.integer() |> tail(n)
## }
## setname <- function(x) {cbind(nm = rownames(x), x)}

## traitMatrix <- sapply(s, function(x) to_bin(x-1)) |>
##   t() |>
##   as.data.frame() |>
##   setname()
