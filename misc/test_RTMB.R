remotes::install_github("bbolker/corHMM", ref = "bolker_clean")
library(corHMM)
library(rbenchmark)
data("primates")
phy <- reorder(primates$tree, "pruningwise")

## clunky, have to find the right place in corHMM (inside random.restarts())
##  to save necessary bits ...
## trace(corHMM, at=46, browser)
## fit1 <- corHMM(phy = phy, data = primates$trait, rate.cat = 1)
## save("model.set.final", "starts", file = "model.set.final.rda")

fit_old <- function() {
  capture.output(fit <- corHMM(phy = phy, data = primates$trait, rate.cat = 1))
  fit
}

fit1 <- fit_old()
print(fit1)

## load previously constructed objects (rate index matrix, likelihood matrix, etc.)
load("misc/model.set.final.rda")
fit_new <- function() {
  ## mkdev.corhmm_rtmb isn't exported yet ...
  ff <- with(model.set.final, corHMM:::mkdev.corhmm_rtmb(log(starts), phy, liks, Q,
                                          rate, root.p = "yang",
                                          rate.cat = 1, order.test = FALSE,
                                          lewis.asc.bias = FALSE, set.fog = FALSE, fog.vec = numeric(0)))
  opt <- with(ff, nlminb(par, fn, gr, lower = rep(log(1e-9), length(par)),
                         upper = rep(log(100), length(par))))

  opt

}

fit2 <- fit_new()

stopifnot(all.equal(-1*fit1$loglik, fit2$objective, tolerance = 2e-8)) ## mean diff: 1.6e-8
stopifnot(all.equal(fit2$par, log(na.omit(c(fit1$solution))),
                    check.attributes  = FALSE, tolerance = 2e-4))
## Mean relative difference: 0.000103487

## takes about 1-2 minutes
bb <- benchmark(fit_old(), fit_new(), replications = 20,
          columns = c("test", "elapsed", "relative"))

## a slightly unfair comparison because fit_old() also does some model setup tasks.
## but the bulk of the time is in the optimization ...

##        test elapsed relative
## 2 fit_new()   0.804    1.000
## 1 fit_old()  50.929   63.345

