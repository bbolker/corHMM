## benchmark mkdev.corhmm_rtmb: pruning = "expm" (per-edge matrix exponential,
## the previous method, now with precomputed child edges) vs "auto" (eigen)
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1)
RhpcBLASctl::blas_set_num_threads(1)
suppressMessages({devtools::load_all("~/R/pkgs/corHMM", quiet = TRUE); library(ape)})
e <- new.env(); sys.source("setup.R", e); mk <- e$mk
mkobj <- function(a, pruning) with(a, mkdev.corhmm_rtmb(p, phy, liks = liks, Q = Q, rate = rate,
  root.p = root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE, pruning = pruning))
set.seed(1)
res <- list()
for (rc in 1:3) for (nt in c(200, 1000, 5000)) {
  a <- mk(nt, rc)
  P <- replicate(20, a$p + rnorm(length(a$p), 0, 0.1))
  for (m in c("expm", "auto")) {
    tt <- system.time(ob <- mkobj(a, m))[3]
    ob$fn(P[, 1]); ob$gr(P[, 1])  ## warm-up
    f <- system.time(for (i in 1:20) ob$fn(P[, i]))[3]/20*1000
    g <- system.time(for (i in 1:20) ob$gr(P[, i]))[3]/20*1000
    res[[length(res)+1]] <- data.frame(states = 2*rc, ntip = nt, pruning = m, tape_s = tt,
                                       fn_ms = f, gr_ms = g, nll = ob$fn(a$p))
    cat(sprintf("states=%d ntip=%d %s tape=%.2fs fn=%.2fms gr=%.2fms\n", 2*rc, nt, m, tt, f, g))
  }
}
res <- do.call(rbind, res)
saveRDS(res, "bench_pkg.rds")
