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

