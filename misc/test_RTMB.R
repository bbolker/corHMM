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

## fit the basic model in corHMM (suppress output)
fitfun_orig <- function() {
  capture.output(fit <- corHMM(phy = phy, data = primates$trait, rate.cat = 1))
  fit
}

fit_orig <- fitfun_orig()
print(fit_orig)

## load previously constructed objects (rate index matrix, likelihood matrix, etc.)
L <- load("misc/model.set.final.rda")
print(L)
## [1] "model.set.final" "starts"         
names(model.set.final)
## [1] "rate.cat"     "np"           "rate"         "index.matrix" "liks"        
## [6] "Q"
##  we need liks, Q, rate
fitfun_RTMB <- function(set.fog = FALSE, lb = 1e-9, ub = 100) {
  ## mkdev.corhmm_rtmb isn't exported yet ...
  ff <- with(model.set.final,
             corHMM:::mkdev.corhmm_rtmb(log(starts), phy, liks, Q,
                                        rate, root.p = "yang",
                                        rate.cat = rate.cat, order.test = FALSE,
                                        lewis.asc.bias = FALSE, set.fog = set.fog,
                                        fog.vec = numeric(0)))
  opt <- with(ff, nlminb(par, fn, gr, lower = rep(log(lb), length(par)),
                         upper = rep(log(ub), length(par))))
  opt
}

fit_RTMB <- fitfun_RTMB()

stopifnot(all.equal(-1*fit_orig$loglik, fit_RTMB$objective, tolerance = 2e-8)) ## mean diff: 1.6e-8
stopifnot(all.equal(fit_RTMB$par, log(na.omit(c(fit_orig$solution))),
                    check.attributes  = FALSE, tolerance = 2e-4))
## Mean relative difference: 0.000103487

## takes about 1-2 minutes
bb <- benchmark(fit_orig(), fit_new(), replications = 20,
          columns = c("test", "elapsed", "relative"))

## a slightly unfair comparison because fit_orig() also does some model setup tasks.
## but the bulk of the time is in the optimization ...

##        test elapsed relative
## 2 fit_new()   0.804    1.000
## 1 fit_orig()  50.929   63.345

