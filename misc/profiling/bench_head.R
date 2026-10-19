## RTMB backend before the eigen change (commit f321b24; rtmb_head.R is
## `git show f321b24:R/RTMB_devcorhmm.R`): which() per node + Matrix::expm per edge
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1)
RhpcBLASctl::blas_set_num_threads(1)
suppressMessages({devtools::load_all("~/R/pkgs/corHMM", quiet = TRUE); library(ape)})
e <- new.env(); sys.source("setup.R", e); mk <- e$mk
h <- new.env(parent = asNamespace("corHMM")); sys.source("rtmb_head.R", h)
set.seed(1)
res <- list()
for (rc in 1:3) for (nt in c(200, 1000, 5000)) {
  a <- mk(nt, rc)
  P <- replicate(20, a$p + rnorm(length(a$p), 0, 0.1))
  tt <- system.time(ob <- with(a, h$mkdev.corhmm_rtmb(p, phy, liks = liks, Q = Q, rate = rate,
    root.p = root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE)))[3]
  ob$fn(P[, 1]); ob$gr(P[, 1])
  f <- system.time(for (i in 1:20) ob$fn(P[, i]))[3]/20*1000
  g <- system.time(for (i in 1:20) ob$gr(P[, i]))[3]/20*1000
  res[[length(res)+1]] <- data.frame(states = 2*rc, ntip = nt, pruning = "HEAD", tape_s = tt,
                                     fn_ms = f, gr_ms = g, nll = ob$fn(a$p))
  cat(sprintf("states=%d ntip=%d HEAD tape=%.2fs fn=%.2fms gr=%.2fms\n", 2*rc, nt, tt, f, g))
}
res <- rbind(do.call(rbind, res), readRDS("bench_pkg.rds"))
saveRDS(res, "bench_pkg.rds")
w <- reshape(res[, c("states","ntip","pruning","tape_s","fn_ms","gr_ms")], idvar = c("states","ntip"),
             timevar = "pruning", direction = "wide")
nll <- tapply(res$nll, list(res$states, res$ntip), \(x) diff(range(x)))
cat("max nll difference across methods:", max(nll), "\n")
print(w, row.names = FALSE, digits = 3)
