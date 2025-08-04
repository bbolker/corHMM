## utility function, copied from lme4
## unnamed elements of ... will be given names from deparse()
namedList <- function (...)  {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}


## version of dev.corhmm that does core likelihood computation with RTMB
##' @param param parameter vector (rates + fog parameters)
##' @param phy phylogenetic tree
##' @param liks
##' @param Q
##' @param rate
##' @param root.p
##' @param rate.cat
##' @param order.test
##' @param lewis.asc.bias
##' @param set.fog
##' @param fog.vec
##' @param build (logical) rebuild RTMB function?
##' @examples
##' data("primates")
##' phy <- multi2di(primates[[1]])
##' data <- primates[[2]]
##' phy <- reorder(primates$tree, "pruningwise")
##' ## corHMM(phy = phy, data = primates$trait, rate.cat = 1)
##' ## save("model.set.final", "starts", file = "model.set.final.rda")
##' load("model.set.final.rda")
##' with(model.set.final, mkdev.corhmm_rtmb(starts, phy, liks, Q, rate, root.p = "maddfitz", rate.cat = 1, order.test = FALSE, lewis.asc.bias = FALSE, set.fog = FALSE, fog.vec = numeric(0)))
##'      

library(RTMB)
mkdev.corhmm_rtmb <- function(p, phy, liks, Q, rate, root.p, rate.cat, order.test, lewis.asc.bias, set.fog=FALSE, fog.vec) {
  ## isolate computations that depend explicitly on p within the RTMB functions
  ## do as much computation as possible (i.e. computation that depends only on other components)
  ##   outside the core function; pass data to core function through tmb_data (not strictly necessary,
  ##   but helps with clarity/modularization)
  ## core function **cannot** contain any hard if() statements (e.g. various bits here that return
  ##  <very-large-number> if some 'bad' condition is met
  nb.node <- Nnode(phy)
  nb.tip <- Ntip(phy)
  TIPS <- seq.int(nb.tip)
  anc <- unique(phy$edge[,1])
  k.rates <- dim(Q)[2] / 2
  
  tmb_data <- namedList(nb.node, nb.tip, TIPS, anc, k.rates)

  prune_fun <- function(pars) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    "diag<-" <- RTMB::ADoverload("diag<-")
    AD <- RTMB::AD
    RTMB::getAll(pars, tmb_data)
    p <- exp(p)
    cp_root.p <- root.p
    comp <- numeric(nb.tip + nb.node)
    
    ## Obtain an object of all the unique ancestors
    ## if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
    
    ## Sets up the liks matrix to account for tip fog:
    if (set.fog) {
      fog_pars <- seq.int(length(unique(fog.vec)))
      ## BMB: will subsetting p in this way confuse RTMB? do we have to do it differently?
      ##  (try it and find out)
      tip.fog.tmp <- p[fog_pars]
      p <- p[-fog_pars]
      tip.fog <- numeric(length(fog.vec))
      tip.fog[] <- c(tip.fog.tmp, 0)[fog.vec]

      ## BMB: not sure why these two code blocks are separate. The only thing I can see that's different is the division by 'rate.cat'
      ##  in one place, which doesn't make a difference since rate.cat==1 in the second block anyway?

      ## BMB: this does depend on tip.fog, which depends on a subset of the parameters, so we do have to do it inside the core function
      
      if(rate.cat > 1){
        ## Error only applies to observed states, but need to replicate across the rate categories:
        ## tip.fog <- rep(tip.fog, rate.cat)
        for(tip.index in seq(nb.tip)) {
          ## Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
          num.zeros <- length(liks[tip.index,which(liks[tip.index,]==0)])
          if(num.zeros > 0){
            liks[tip.index,which(liks[tip.index,]==1)] <- 1 - (sum(tip.fog[which(liks[tip.index,]!=1)])/rate.cat)
            liks[tip.index,which(liks[tip.index,]==0)] <- tip.fog[which(liks[tip.index,]==0)]
          }
        }
      } else {
        for(tip.index in seq(nb.tip)) {
                                        #Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
          num.zeros <- length(liks[tip.index,which(liks[tip.index,]==0)])
          if(num.zeros > 0){
            liks[tip.index,which(liks[tip.index,]==1)] <- 1 - sum(tip.fog[which(liks[tip.index,]!=1)])
            liks[tip.index,which(liks[tip.index,]==0)] <- tip.fog[which(liks[tip.index,]==0)]
          }
        }
      }
    }

    ## prune_nll uses:
    ## Q[Q!=0] <- trans_rates[Q[Q!=0]] 
    ## diag(Q) <- -1*rowSums(Q)
    ## return(Q)

    ## ... what's different here??
    ## jumping through all kinds of hoops with rowSums()/diag()<-, still not sure ...

    np <- length(p)
    Q <- rate
    dimnames(Q) <- NULL
    Q[Q<=np] <- p[Q[Q<=np]]
    Q[rate==np+1] <- 0

    for (i in 1:nrow(Q)) {
      Q[i,i] <- -1*sum(Q[i,])
    }
    
    ## an error check for rates that are really weird

    ## have to use a different expm(), but I believe this is the same method as expm::expm(., method = "Ward77")
    test_mat <- Matrix::expm(Q)

    ## rowSums is problematic again (and, we can't use if() anyway
    ## if(any(round(rowSums(test_mat))> 2)) {
    ##  stop("bad result in expm()")
    ## }
    
    ## # if the q matrix has columns not estimated, remove them
    ## row2rm <- apply(rate, 1, function(x) all(x == max(rate)))
    ## col2rm <- apply(rate, 2, function(x) all(x == max(rate)))
    ## Q.root <- Q[!row2rm | !col2rm, !row2rm | !col2rm]
    
    if(is.character(root.p)){
      if(root.p == "yang"){
        root.test <- Null(Q)
        if(dim(root.test)[2]>1){
          stop("bad root test for 'yang' root probabilities")
        }
      }      
    }

    if(order.test) {
                                        # ensure that the rate classes have mean rates in a consistent order (A > B > C > n)
      StateOrderMat <- matrix(1, (dim(Q)/rate.cat)[1], (dim(Q)/rate.cat)[2])
      RateClassOrderMat <- matrix(0, rate.cat, rate.cat)
      diag(RateClAssordermat) <- 1:rate.cat
      OrderMat <- RateClassOrderMat %x% StateOrderMat
      Rate01 <- vector("numeric", rate.cat)
      for(i in 1:rate.cat){
        tmp <- Q[OrderMat == i]
        Rate01[i] <- tmp[tmp>=0][1]
      }
      OrderTest <- all.equal(Rate01, sort(Rate01, decreasing = TRUE))
      if(!OrderTest) {
        stop("bad order test")
      }
    }

    ## pruning algo itself!
    for (i in seq(from = 1, length.out = nb.node)) {
      
      ##the ancestral node at row i is called focal
      focal <- anc[i]
      ##Get descendant information of focal
      desRows <- which(phy$edge[,1]==focal)
      desNodes <- phy$edge[desRows,2]
      v <- 1
      ##Loops through all descendants of focal (how we deal with polytomies):
      ## BMB need drop(as.matrix(...)) here, sigh
      for (desIndex in seq_along(desRows)) {
        v <- drop(as.matrix(v*Matrix::expm(Q * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]))
      }
      
      ##Allows for fixed nodes based on user input tree.
      if(!is.null(phy$node.label)){
        if(!is.na(phy$node.label[focal - nb.tip])){
          fixer.tmp <- numeric(dim(Q)[2]/rate.cat)
          fixer.tmp[phy$node.label[focal - nb.tip]] <- 1
          fixer <- rep(fixer.tmp, rate.cat)
          v <- v * fixer
        }
      }
      
      ##Sum the likelihoods:
      comp[focal] <- sum(v)
      ##Divide each likelihood by the sum to obtain probabilities:
      liks[focal, ] <- v/comp[focal]
    }

    ##Specifies the root:
    root <- nb.tip + 1L
    
    ##If any of the logs have NAs restart search:
    ## if (is.na(sum(log(comp[-TIPS])))){ stop("NA values in logs") }

    ## this was NULL: just for constructing vector?
    equil.root <- numeric(ncol(Q))

    ## assuming 'positive' means 'not diagonal'
    QQ <- Q
    diag(QQ) <- NA
    for(i in 1:ncol(Q)){
      rowsum <- sum(Q[,i], na.rm=TRUE)
      colsum <- sum(Q[i,], na.rm=TRUE)
      equil.root[i] <- rowsum/(rowsum+colsum)
    }
    if (is.null(root.p)){
      flat.root = equil.root
      k.rates <- 1/length(which(!is.na(equil.root)))
      flat.root[!is.na(flat.root)] = k.rates
      flat.root[is.na(flat.root)] = 0
      loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
    }

    if (is.character(root.p)){
      ## root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
      ## i.e. the *stationary distribution*. Null(Q) is giving RTMB a headache so substitute
      ## a big matrix exponential ...
      if(root.p == "yang"){
        ## root.p <- Null(Q)
        ## root.p <- c(root.p/sum(root.p))
        root.p <- Matrix::expm(10000*Q)[1,]
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
        ## if(is.infinite(loglik)){
        ##   stop("infinite loglik in yang root prob computation")
        ## }
      } else {
        ## root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
        root.p = liks[root,] / sum(liks[root,])
        ## BMB: switch to logsumexp formulation?
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
      }
    } else {
      if(is.numeric(root.p[1])){
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
        ## if(is.infinite(loglik)){
        ##   stop("infinite log-likelihood")
        ## }
      }
    }

    ## root.p!==NULL will fix root probabilities based on user supplied vector:
    if(lewis.asc.bias) {
      p <- log(p)
      dummy.liks.vec <- getLewisLikelihood(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = cp_root.p, rate.cat = rate.cat)
      loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
    }
    return(loglik)
  }

  RTMB::MakeADFun(prune_fun, list(p=p), silent = TRUE)
  
}

