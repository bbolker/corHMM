source("proto_level.R")
a <- mk(200, 2); np <- length(a$p)
o1 <- mk_obj(a, "expm"); o3 <- mk_obj_lev(a)
Qf <- function(p) { Q <- a$rate; Q[Q <= np] <- exp(p)[Q[Q <= np]]; Q[a$rate == np+1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q); Q }
set.seed(3); z <- rnorm(np)
for (eps in 10^-(1:10)) {
  p <- eps * z; ev <- eigen(Qf(p))$values
  gap <- min(dist(ev)[dist(ev) > 0]); g1 <- o1$gr(p); g3 <- o3$gr(p)
  rgap <- min(Mod(outer(ev, ev, "-"))[upper.tri(diag(4))]) / max(Mod(ev))
  cat(sprintf("eps=%.0e  rel eigen gap=%.1e  max |grad err|=%.2e  kappa(V)=%.1e\n", eps, rgap,
      max(abs(g1 - g3)), kappa(eigen(Qf(p))$vectors, exact = TRUE)))
}
