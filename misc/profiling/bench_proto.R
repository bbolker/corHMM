source("proto_eigen.R")
set.seed(1)
res <- list()
for (rc in c(1, 2, 3)) for (nt in c(200, 1000, 5000)) {
  a <- mk(nt, rc)
  P <- replicate(20, a$p + rnorm(length(a$p), 0, 0.1))
  for (m in c("expm", "eigen")) {
    tt <- system.time(ob <- mk_obj(a, m))[3]
    f <- system.time(for (i in 1:20) ob$fn(P[, i]))[3]/20*1000
    g <- system.time(for (i in 1:20) ob$gr(P[, i]))[3]/20*1000
    res[[length(res)+1]] <- data.frame(states = 2*rc, ntip = nt, method = m,
      tape_s = round(tt, 2), fn_ms = round(f, 2), gr_ms = round(g, 2), nll = ob$fn(a$p))
    cat(sprintf("states=%d ntip=%d %s tape=%.2f fn=%.2f gr=%.2f\n", 2*rc, nt, m, tt, f, g))
  }
}
saveRDS(do.call(rbind, res), "bench_proto.rds")
