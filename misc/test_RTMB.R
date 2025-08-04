devtools::load_all()
data("primates")
phy <- reorder(primates$tree, "pruningwise")
## clunky, have to find the right place in corHMM (inside random.restarts())
##  to save necessary bits ...
## trace(corHMM, at=46, browser)
fit1 <- corHMM(phy = phy, data = primates$trait, rate.cat = 1)
## save("model.set.final", "starts", file = "model.set.final.rda")
load("model.set.final.rda")


## options(error=recover)
debug(mkdev.corhmm_rtmb)
## BMB: have to make this work with root.p = "yang" now ...
ff <- with(model.set.final, mkdev.corhmm_rtmb(starts, phy, liks, Q,
                                              rate, root.p = "yang",
                                              rate.cat = 1, order.test = FALSE,
                                              lewis.asc.bias = FALSE, set.fog = FALSE, fog.vec = numeric(0)))
ff$fn()
ff$gr()

opt <- with(ff, nlminb(par, fn, gr, lower = rep(log(1e-9), length(par)),
                upper = rep(log(100), length(par))))

opt$par
opt$objective

## compare basic pruning example
## what is different about the Q-setting machinery?
hdir <- "~/students/haoyu/pruning/R/" 
source(file.path(hdir, "pruning_funs.R"))
source(file.path(hdir, "Q_template.R"))
source(file.path(hdir, "loglik.R"))

Q2 <- model.set.final$rate
Q2[Q2==5] <- 0

traitList <- multi_to_single(primates$trait[,-1])
traitList[traitList==4] <- 3 ## collapse
Phylodata <- list(Q_template = Q2,
                  tree = phy, trait_values = traitList,
                  traitMatrix = primates$trait)
f0 <- cmb(prune_nll, Phylodata)
p0 <- list(log_trans_rates = starts)
debug(f0)
f0(p0)
ff <- MakeADFun(f0,
                list(log_trans_rates = starts), silent = TRUE)
ff$fn()


Q0 <- matrix(nrow=3, ncol=3, byrow = TRUE, c(
 -1.022380,1.02238,0.000000,
  1.029088,-2.03769,1.008601,
  0.000000,1.01117,-1.011170))
MASS::Null(Q0)
