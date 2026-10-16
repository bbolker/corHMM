source("setup.R")
Rcpp::sourceCpp("prune.cpp")
buildQ <- function(a) { Q <- a$Q; Q[] <- c(exp(a$p), 0)[a$rate]; diag(Q) <- -rowSums(Q); Q }
finish <- function(comp, L, Q, ntip, root) {
  rp <- Null(Q); rp <- c(rp/sum(rp))
  -(sum(log(comp[-(1:ntip)])) + log(sum(rp * L[root,])))
}
## B: same R loop, but P(t) from one eigendecomposition instead of expm per edge
dev_eigen <- function(a) {
  Q <- buildQ(a); phy <- a$phy; L <- a$liks; nt <- Ntip(phy)
  ev <- eigen(Q); V <- ev$vectors; Vi <- solve(V); lam <- ev$values
  comp <- numeric(nrow(L))
  for (focal in unique(phy$edge[,1])) {
    rows <- which(phy$edge[,1] == focal); v <- 1
    for (r in rows) v <- v * Re(V %*% (exp(lam * phy$edge.length[r]) * (Vi %*% L[phy$edge[r,2],])))
    comp[focal] <- sum(v); L[focal,] <- v/comp[focal]
  }
  finish(comp, L, Q, nt, nt+1)
}
dev_cpp <- function(a, method) {
  Q <- buildQ(a); nt <- Ntip(a$phy)
  r <- prune_cpp(a$phy$edge, a$phy$edge.length, a$liks, Q, method)
  finish(r$comp, r$L, Q, nt, nt+1)
}
mkrtmb <- function(a) with(a, corHMM:::mkdev.corhmm_rtmb(p, phy, liks=liks, Q=Q, rate=rate, root.p=root.p,
     rate.cat=rate.cat, order.test=FALSE, lewis.asc.bias=FALSE, set.fog=FALSE, fog.vec=NULL, q_prep=NULL, prior=NULL))

res <- list()
for (rc in 1:3) for (nt in c(200, 1000, 5000)) {
  a <- mk(nt, rate.cat=rc)
  vals <- c(A=do.call(dev_repo, a), B=dev_eigen(a), C=dev_cpp(a,0), D=dev_cpp(a,1))
  ob <- tryCatch(mkrtmb(a), error=function(e) {message("RTMB: ", conditionMessage(e)); NULL})
  if (!is.null(ob)) vals["RTMB"] <- ob$fn(a$p)
  stopifnot(diff(range(vals)) < 1e-6 * abs(vals[1]))
  nrep <- if (nt >= 5000) 3 else 10
  ex <- list(A=quote(do.call(dev_repo, a)), B=quote(dev_eigen(a)), C=quote(dev_cpp(a,0)), D=quote(dev_cpp(a,1)))
  if (!is.null(ob)) ex <- c(ex, list(RTMB_fn=quote(ob$fn(a$p)), RTMB_gr=quote(ob$gr(a$p))))
  mb <- microbenchmark::microbenchmark(list=ex, times=nrep)
  s <- summary(mb, unit="ms")
  res[[length(res)+1]] <- data.frame(states=2*rc, ntip=nt, expr=s$expr, median_ms=round(s$median, 2))
  cat(sprintf("states=%d ntip=%d loglik=%.4f done\n", 2*rc, nt, vals[1]))
}
res <- do.call(rbind, res)
w <- reshape(res, idvar=c("states","ntip"), timevar="expr", direction="wide")
names(w) <- sub("median_ms.", "", names(w))
print(w, row.names=FALSE)
saveRDS(res, "bench_res.rds")

## profile A
a <- mk(1000, rate.cat=2)
Rprof("prof.out", interval=0.005)
for (i in 1:20) do.call(dev_repo, a)
Rprof(NULL)
sp <- summaryRprof("prof.out")
print(head(sp$by.self, 12)); print(head(sp$by.total, 12))
