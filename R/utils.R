namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}

#' control parameters for MCMC runs
#'
#' All elements are passed to \code{rstan::sampling()} through
#' \code{tmbstan::tmbstan()}; see those help pages.
#' @param chains number of MCMC chains
#' @param iter MCMC iterations per chain
#' @param ... other arguments to \code{tmbstan::tmbstan()}, e.g. \code{warmup},
#' \code{cores}, \code{seed}, \code{control}, \code{open_progress}
#' @export
mcmcControl <- function(chains = 4, iter = 2000, ...) {
    c(namedList(chains, iter), list(...))
}
