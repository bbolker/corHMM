## level-wise vectorized eigen pruning: all edges whose parent is at the same
## height (max child level + 1) are handled with one set of matrix operations
source("proto_eigen.R")
prune_levels <- function(phy) {
  nb.tip <- Ntip(phy); e <- phy$edge
  lev <- integer(nb.tip + Nnode(phy))
  for (r in rev(seq_len(nrow(e)))) {}  ## placeholder
  ## pruningwise order: children before parents
  for (r in seq_len(nrow(e))) lev[e[r,1]] <- max(lev[e[r,1]], lev[e[r,2]] + 1L)
  lapply(sort(unique(lev[-(1:nb.tip)])), function(L) {
    rows <- which(lev[e[,1]] == L)
    rows <- rows[order(e[rows,1])]
    par <- e[rows,1]
    nodes <- unique(par)
    slot <- ave(seq_along(par), par, FUN = seq_along)
    list(rows = rows, child = e[rows,2], t = phy$edge.length[rows], nodes = nodes,
         slot_idx = lapply(seq_len(max(slot)), function(k) {
           list(src = which(slot == k), dst = match(par[slot == k], nodes))
         }))
  })
}
mk_obj_lev <- function(a) {
  phy <- a$phy; liks0 <- a$liks; rate <- a$rate
  nb.tip <- Ntip(phy); nb.node <- Nnode(phy)
  levs <- prune_levels(phy)
  f <- function(pars) {
    "[<-" <- RTMB::ADoverload("[<-"); "c" <- RTMB::ADoverload("c")
    p <- pars$p; liks <- liks0
    np <- length(p); p_rates <- exp(p)
    Q <- rate; Q[Q <= np] <- p_rates[Q[Q <= np]]; Q[rate == np + 1] <- 0
    for (i in seq_len(nrow(Q))) Q[i, i] <- -sum(Q[i, ])
    comp <- numeric(nb.tip + nb.node)
    ev <- eigen(Q, symmetric = FALSE)
    V <- ev$vectors; Vi <- solve(V); lam <- ev$values
    ## real and imaginary parts, so the loop uses real AD arithmetic only
    tVr <- t(Re(V)); tVim <- t(Im(V)); tVir <- t(Re(Vi)); tViim <- t(Im(Vi))
    lr <- matrix(Re(lam), nrow = 1); li <- matrix(Im(lam), nrow = 1)
    for (L in levs) {
      X <- liks[L$child, , drop = FALSE]
      tt <- matrix(L$t, ncol = 1)
      ea <- exp(tt %*% lr); tb <- tt %*% li
      Er <- ea * cos(tb); Ei <- ea * sin(tb)
      Wr <- X %*% tVir; Wi <- X %*% tViim
      Y <- (Wr * Er - Wi * Ei) %*% tVr - (Wr * Ei + Wi * Er) %*% tVim
      s1 <- L$slot_idx[[1]]
      v <- Y[s1$src, , drop = FALSE]
      for (s in L$slot_idx[-1]) v[s$dst, ] <- v[s$dst, , drop = FALSE] * Y[s$src, , drop = FALSE]
      cs <- rowSums(v)
      comp[L$nodes] <- cs
      liks[L$nodes, ] <- v / cs
    }
    root.p <- Matrix::expm(10000 * Q)[1, ]
    -(sum(log(comp[-(1:nb.tip)])) + log(sum(root.p * liks[nb.tip + 1, ])))
  }
  MakeADFun(f, list(p = a$p), silent = TRUE)
}
