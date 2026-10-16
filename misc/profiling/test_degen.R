source("proto_level.R")
chk <- function(a, p, lab) {
  o1 <- mk_obj(a, "expm"); o3 <- mk_obj_lev(a)
  ev <- eigen({Q <- a$rate; np <- length(p); Q[Q <= np] <- exp(p)[Q[Q <= np]]; Q[a$rate == np+1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q); Q})$values
  cat("\n==", lab, " eigenvalues:", format(ev, digits = 3), "\n")
  print(rbind(expm = c(o1$fn(p), o1$gr(p)), level = c(o3$fn(p), o3$gr(p))))
}
a <- mk(200, 1); chk(a, rep(0, length(a$p)), "2 states, equal rates")
a <- mk(200, 2); chk(a, rep(0, length(a$p)), "4 states (rate.cat=2), all rates equal")
a <- mk(200, 3); chk(a, rep(log(0.5), length(a$p)), "6 states, all equal")
a <- mk(200, 2); chk(a, rep(0, length(a$p)) + c(1e-6, 0), "4 states, nearly equal")
