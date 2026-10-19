source("proto_level.R")
a <- mk(200, 2); p1 <- a$p + rnorm(length(a$p), 0, 0.1)
o1 <- mk_obj(a, "expm"); o3 <- mk_obj_lev(a)
print(c(o1$fn(p1), o3$fn(p1))); print(rbind(o1$gr(p1), o3$gr(p1)))
for (rc in c(1,2,3)) for (nt in c(1000, 5000)) {
  a <- mk(nt, rc); P <- replicate(20, a$p + rnorm(length(a$p), 0, 0.1))
  tt <- system.time(ob <- mk_obj_lev(a))[3]
  f <- system.time(for (i in 1:20) ob$fn(P[, i]))[3]/20*1000
  g <- system.time(for (i in 1:20) ob$gr(P[, i]))[3]/20*1000
  cat(sprintf("states=%d ntip=%d level tape=%.2f fn=%.2f gr=%.2f nlev=%d\n", 2*rc, nt, tt, f, g, length(prune_levels(a$phy))))
}
