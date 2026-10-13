\name{mcmcControl}
\alias{mcmcControl}
\title{Control settings for HMC sampling}
\description{Builds the list of settings used by \code{\link{corHMM}} when its \code{mcmc} argument is supplied. All elements are passed to \code{rstan::sampling} through \code{tmbstan::tmbstan}.}
\usage{
mcmcControl(chains = 4, iter = 2000, ...)
}
\arguments{
\item{chains}{number of MCMC chains.}
\item{iter}{number of iterations per chain, including warmup.}
\item{...}{other arguments to \code{tmbstan::tmbstan}, e.g. \code{warmup} (default \code{floor(iter/2)}), \code{cores}, \code{seed}, \code{refresh}, \code{open_progress}, or \code{control} (a list such as \code{list(adapt_delta = 0.9)}). Settings that are not given take the \code{rstan} defaults.}
}
\value{a named list of sampler settings.}
\seealso{\code{\link{corHMM}}}
\examples{
mcmcControl(chains = 2, iter = 1000, warmup = 500, seed = 101)
\donttest{
if (requireNamespace("tmbstan", quietly = TRUE)) {
  data(primates)
  phy <- multi2di(primates[[1]])
  data <- primates[[2]][, 1:2]
  fit <- corHMM(phy, data, rate.cat = 1, model = "ER", node.states = "none",
                use_RTMB = TRUE, prior = list(),
                mcmc = mcmcControl(chains = 2, iter = 500, seed = 101, refresh = 0))
  fit$mcmc
}
}
}
\keyword{models}
