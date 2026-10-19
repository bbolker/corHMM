Sys.setenv(OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1)
RhpcBLASctl::blas_set_num_threads(1)
suppressMessages({library(corHMM); library(ape); library(expm); library(MASS)})
## use the repo's current dev.corhmm (repo HEAD doesn't parse as a whole)
src <- suppressWarnings(readLines("~/R/pkgs/corHMM/R/corHMM.R"))
i0 <- grep("^dev.corhmm <- function", src); i1 <- i0 - 1 + which(src[i0:length(src)] == "}")[1]
dev_repo <- eval(parse(text = sub("^dev.corhmm <-", "function_tmp <-", src[i0:i1]))[[1]][[3]]); environment(dev_repo) <- asNamespace("corHMM")
mk <- function(ntip, rate.cat=2, seed=101) {
  set.seed(seed)
  phy <- phytools::pbtree(n=ntip, scale=1)
  Qs <- matrix(c(-1,1,1,-1),2,dimnames=list(0:1,0:1))*2
  x <- phytools::sim.Mk(phy, Qs, anc="0")
  dat <- data.frame(sp=phy$tip.label, x=as.character(x))
  ns <- asNamespace("corHMM")
  CorData <- ns$corProcessData(dat, collapse=TRUE)
  m <- ns$match.tree.data(phy, CorData$corData); phy <- m$phy
  ms <- ns$rate.cat.set.corHMM.JDB(phy=phy, data=dat, rate.cat=rate.cat, ntraits=2,
                                   model="ARD", rate.mat=NULL, collapse=TRUE)
  phy <- reorder(phy, "pruningwise")
  np <- ms$np
  list(p=log(seq(3, 0.5, length.out=np)), phy=phy, liks=ms$liks, Q=ms$Q, rate=ms$rate,
       root.p="yang", rate.cat=rate.cat, order.test=FALSE, lewis.asc.bias=FALSE)
}
