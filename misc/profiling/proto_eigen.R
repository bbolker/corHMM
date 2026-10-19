## prototype: RTMB pruning with one (AD) eigendecomposition per evaluation
source("setup.R")
library(RTMB)
mk_obj <- function(a, method = c("expm", "eigen")) {
  method <- match.arg(method)
  phy <- a$phy; liks0 <- a$liks; rate <- a$rate
  nb.tip <- Ntip(phy); nb.node <- Nnode(phy)
  anc <- unique(phy$edge[,1])
  child_rows <- split(seq_len(nrow(phy$edge)), phy$edge[,1])[as.character(anc)]
  f <- function(pars) {
    "[<-" <- RTMB::ADoverload("[<-"); "c" <- RTMB::ADoverload("c")
    p <- pars$p; liks <- liks0
    np <- length(p); p_rates <- exp(p)
    Q <- rate; Q[Q <= np] <- p_rates[Q[Q <= np]]; Q[rate == np + 1] <- 0
    for (i in seq_len(nrow(Q))) Q[i, i] <- -sum(Q[i, ])
    comp <- numeric(nb.tip + nb.node)
    if (method == "eigen") {
      ev <- eigen(Q, symmetric = FALSE)
      V <- ev$vectors; Vi <- solve(V); lam <- ev$values
    }
    for (i in seq_len(nb.node)) {
      focal <- anc[i]; rows <- child_rows[[i]]; v <- 1
      for (r in rows) {
        x <- liks[phy$edge[r, 2], ]
        if (method == "expm") {
          v <- drop(as.matrix(v * Matrix::expm(Q * phy$edge.length[r]) %*% x))
        } else {
          v <- v * Re(drop(V %*% (exp(lam * phy$edge.length[r]) * drop(Vi %*% x))))
        }
      }
      comp[focal] <- sum(v); liks[focal, ] <- v / comp[focal]
    }
    root.p <- Matrix::expm(10000 * Q)[1, ]
    -(sum(log(comp[-(1:nb.tip)])) + log(sum(root.p * liks[nb.tip + 1, ])))
  }
  MakeADFun(f, list(p = a$p), silent = TRUE)
}
a <- mk(200, 2)
o1 <- mk_obj(a, "expm"); o2 <- mk_obj(a, "eigen")
p1 <- a$p + rnorm(length(a$p), 0, 0.1)
print(c(o1$fn(p1), o2$fn(p1), do.call(dev_repo, modifyList(a, list(p = p1)))))
print(rbind(o1$gr(p1), o2$gr(p1), numDeriv::grad(\(p) do.call(dev_repo, modifyList(a, list(p = p))), p1)))
