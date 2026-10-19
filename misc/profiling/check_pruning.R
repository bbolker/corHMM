## check mkdev.corhmm_rtmb(pruning = "auto"/"eigen") against "expm"
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
suppressMessages({devtools::load_all("~/R/pkgs/corHMM", quiet = TRUE); library(ape)})
source("setup.R", local = new.env())  ## for mk()
mk <- get("mk", envir = (function() { e <- new.env(); sys.source("setup.R", e); e })())
mkobj <- function(a, pruning, ...) with(a, mkdev.corhmm_rtmb(p, phy, liks = liks, Q = Q, rate = rate,
  root.p = root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE,
  pruning = pruning, ...))
cmp <- function(a, p, lab, ...) {
  oe <- mkobj(a, "expm", ...); oa <- mkobj(a, "auto", ...)
  d <- c(fn = abs(oe$fn(p) - oa$fn(p)), gr = max(abs(oe$gr(p) - oa$gr(p))))
  cat(sprintf("%-45s fn diff %.1e  gr diff %.1e  %s\n", lab, d[1], d[2],
              if (all(d < 1e-6)) "OK" else "MISMATCH"))
}
set.seed(1)
a <- mk(300, 2); np <- length(a$p)
cmp(a, a$p, "ARD 4 states, distinct rates")
cmp(a, rep(0, np), "ARD 4 states, all rates equal (fallback)")
a1 <- mk(300, 1); cmp(a1, rep(-1, 2), "2 states, equal rates")
a3 <- mk(300, 3); cmp(a3, a3$p, "ARD 6 states")
## fixed nodes
af <- a; af$phy$node.label <- sample(c(NA, 1, 2), Nnode(af$phy), replace = TRUE, prob = c(.8, .1, .1))
cmp(af, af$p, "fixed node labels")
## ER with 3 observed states (structural degeneracy)
set.seed(2); phy <- phytools::pbtree(n = 200, scale = 1)
x <- sample(c("a","b","c"), 200, replace = TRUE)
dat <- data.frame(sp = phy$tip.label, x = x)
ms <- rate.cat.set.corHMM.JDB(phy = phy, data = dat, rate.cat = 1, ntraits = 3, model = "ER",
                              rate.mat = NULL, collapse = TRUE)
aer <- list(p = log(0.7), phy = reorder(phy, "pruningwise"), liks = ms$liks, Q = ms$Q, rate = ms$rate,
            root.p = "yang", rate.cat = 1)
oa <- mkobj(aer, "auto"); cat("ER 3 states: auto picks expm tape directly:", is.null(environment(oa$fn)$obj_expm), "\n")
cmp(aer, log(0.7), "ER 3 states (structural)")
## full corHMM fits, old vs new
data(primates, package = "corHMM")
pd <- primates$trait[, c(1, 2)]
f_new <- corHMM(primates$tree, pd, rate.cat = 2, model = "ARD", use_RTMB = TRUE, nstarts = 0)
