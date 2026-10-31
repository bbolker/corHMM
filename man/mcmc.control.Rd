\name{mcmc.control}
\alias{mcmc.control}
\title{Control settings for HMC sampling}
\description{Builds the list of settings used by \code{\link{corHMM}} when its \code{mcmc} argument is supplied. All elements except \code{use.bounds} are passed to \code{rstan::sampling} through \code{tmbstan::tmbstan}.}
\usage{
mcmc.control(chains = 4, iter = 2000, use.bounds = TRUE, ...)
}
\arguments{
\item{chains}{number of MCMC chains.}
\item{iter}{number of iterations per chain, including warmup.}
\item{use.bounds}{for native (non-formula) models, sample the log-rates within the bounds used for optimization, \code{[log(lower.bound), log(upper.bound)]}, narrowed to a uniform prior's limits if there is one (using Stan's transformation to a bounded interval). The defaults \code{lower.bound = 1e-9} and \code{upper.bound = 100} are close to BayesTraits' default rate limits (\code{SetMinMaxRate}, 1e-32 to 100), which BayesTraits also applies when sampling. Ignored for formula models, whose optimization bounds apply to the coefficients (including contrasts) and would truncate their posteriors arbitrarily; they are sampled without bounds. Also ignored if \code{lower} or \code{upper} are given explicitly. Setting \code{use.bounds = FALSE} with a uniform prior gives unreliable sampling (with a warning).}
\item{...}{other arguments to \code{tmbstan::tmbstan}, e.g. \code{warmup} (default \code{floor(iter/2)}), \code{cores}, \code{seed}, \code{refresh}, \code{open_progress}, or \code{control} (a list such as \code{list(adapt_delta = 0.9)}). Settings that are not given take the \code{rstan} defaults.}
}
\value{a named list of sampler settings.}
\seealso{\code{\link{corHMM}}}
\examples{
mcmc.control(chains = 2, iter = 1000, warmup = 500, seed = 101)
\donttest{
if (requireNamespace("tmbstan", quietly = TRUE)) {
  data(primates)
  phy <- multi2di(primates[[1]])
  data <- primates[[2]][, 1:2]
  fit <- corHMM(phy, data, rate.cat = 1, model = "ER", node.states = "none",
                use_RTMB = TRUE, prior = list(),
                mcmc = mcmc.control(chains = 2, iter = 500, seed = 101, refresh = 0))
  fit$mcmc
}
}
}
\keyword{models}
