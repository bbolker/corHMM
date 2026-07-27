namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}

#' control parameters for MCMC runs
#' @param chains number of MCMC chains
#' @param iter MCMC iterations per chain
#' note that these are all passed to rstan::sampling() through tmbstan::tmbstan(),
#' with some modifications - see the relevant help pages
mcmcControl <- function(chains = 4, iter = 2000, ...) {
    namedList(chains, iter)
}
