source("setup.R")
mkrtmb <- function(a) with(a, corHMM:::mkdev.corhmm_rtmb(p, phy, liks=liks, Q=Q, rate=rate, root.p=root.p,
     rate.cat=rate.cat, order.test=FALSE, lewis.asc.bias=FALSE, set.fog=FALSE, fog.vec=NULL, q_prep=NULL, prior=NULL))
for (nt in c(200,1000,5000)) {
  a <- mk(nt, 2); tt <- system.time(ob <- mkrtmb(a))[3]
  P <- replicate(20, a$p + rnorm(length(a$p), 0, 0.1))
  f <- system.time(for (i in 1:20) ob$fn(P[,i]))[3]/20*1000
  g <- system.time(for (i in 1:20) ob$gr(P[,i]))[3]/20*1000
  chk <- c(ob$fn(P[,1]), do.call(dev_repo, modifyList(a, list(p=P[,1]))))
  cat(sprintf("ntip=%d tape=%.2fs fn=%.2fms gr=%.2fms  check %.6f %.6f\n", nt, tt, f, g, chk[1], chk[2]))
}
