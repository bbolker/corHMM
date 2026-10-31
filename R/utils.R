namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}

#' control parameters for MCMC runs
#'
#' All elements except \code{use.bounds} are passed to \code{rstan::sampling()}
#' through \code{tmbstan::tmbstan()}; see those help pages.
#' @param chains number of MCMC chains
#' @param iter MCMC iterations per chain
#' @param use.bounds for native (non-formula) models, sample the log-rates
#' within the optimizer's bounds (see \code{\link{corHMM}})
#' @param ... other arguments to \code{tmbstan::tmbstan()}, e.g. \code{warmup},
#' \code{cores}, \code{seed}, \code{control}, \code{open_progress}
#' @export
mcmc.control <- function(chains = 4, iter = 2000, use.bounds = TRUE, ...) {
    c(namedList(chains, iter, use.bounds), list(...))
}
