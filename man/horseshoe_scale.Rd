\name{horseshoe_scale}
\alias{horseshoe_scale}
\alias{horseshoe_sigma}
\alias{horseshoe_prop_zero}
\title{Global scale of a horseshoe prior from a prior proportion of zero rates}
\description{
Chooses \code{scale_global} (\eqn{\tau_0}) for \code{prior = list(family = "horseshoe")} in \code{\link{corHMM}} so that a given proportion of the rates is expected to be shrunk to (effectively) zero, following Piironen and Vehtari (2017, Sec. 3.3), and simulates the prior proportion of zero rates that a choice of \eqn{\tau_0} implies.
}
\usage{
horseshoe_scale(fit = NULL, prop_zero, sigma = NULL, lower = 1e-9, upper = 1e8)
horseshoe_sigma(fit, lower = 1e-9, upper = 1e8)
horseshoe_prop_zero(tau0, sigma, fix_global = FALSE, slab_scale = Inf, n = 1e4)
}
\arguments{
\item{fit}{a maximum-likelihood \code{corHMM} fit (no \code{prior}) of a native model (specified with \code{rate.mat}, not a formula), with \code{use_RTMB = TRUE} and \code{return.devfun = TRUE}. Not used by \code{horseshoe_scale} if \code{sigma} is given.}
\item{prop_zero}{the prior guess for the proportion of rates that are zero, strictly between 0 and 1.}
\item{sigma}{the likelihood scale of each rate (see Details); by default computed by \code{horseshoe_sigma(fit)}.}
\item{lower, upper}{the range of rates searched for the likelihood scales (use the \code{lower.bound} of the fit and at least its \code{upper.bound}).}
\item{tau0}{the global scale \eqn{\tau_0} (\code{scale_global}).}
\item{fix_global, slab_scale}{as in the \code{"horseshoe"} prior (see \code{\link{corHMM}}): if \code{fix_global = FALSE}, \eqn{\tau \sim C^+(0, \tau_0)}, otherwise \eqn{\tau = \tau_0}; \code{slab_scale} is the slab scale \eqn{c} (\code{Inf}: no slab).}
\item{n}{the number of simulated draws.}
}
\details{
With a Gaussian approximation to the likelihood of each rate \eqn{q_j}, with scale \eqn{\sigma_j} and the rates treated as uncorrelated, the shrinkage factor of rate \eqn{j} under the horseshoe is \eqn{\kappa_j = 1/(1 + \tau^2\lambda_j^2/\sigma_j^2)} (0: not shrunk; 1: shrunk to zero), and \eqn{E(\kappa_j \mid \tau) = 1/(1 + \tau/\sigma_j)}. The effective number of non-zero rates is \eqn{m_{eff} = \sum_j (1 - \kappa_j)}, so \eqn{E(m_{eff} \mid \tau) = \sum_j 1/(1 + \sigma_j/\tau)}. \code{horseshoe_scale} solves \eqn{E(m_{eff} \mid \tau_0) = p_0} with \eqn{p_0 = (1 - } \code{prop_zero}\eqn{) D} for \eqn{D} rates. With equal \eqn{\sigma_j = \sigma} this is Piironen and Vehtari's Eq. 3.13, \eqn{\tau_0 = \sigma p_0/(D - p_0)}, which is returned as the attribute \code{tau0_pooled} (with \eqn{\sigma} the geometric mean of the finite \eqn{\sigma_j}).

\code{horseshoe_sigma} measures \eqn{\sigma_j} on the rate scale, as the half-width of the interval over which the log-likelihood, as a function of rate \eqn{j} with the other rates at their maximum-likelihood values, is within 1/2 of its maximum. For a well-estimated rate this is about its standard error; for a rate whose estimate is at or near zero (where the Hessian of the log-likelihood is degenerate), it is the rate at which the data first distinguish it from zero. A rate whose log-likelihood never drops by 1/2 below \code{upper} gets \eqn{\sigma_j = \infty}: it is shrunk whatever \eqn{\tau}, so the proportion of such rates is a lower limit for \code{prop_zero}.

With \code{fix_global = TRUE}, the prior mean of the proportion of zero rates, \eqn{1 - m_{eff}/D}, is \code{prop_zero}. With \eqn{\tau \sim C^+(0, \tau_0)} (the default of the \code{"horseshoe"} prior) its median is close to \code{prop_zero}, but the distribution is wider; \code{horseshoe_prop_zero} returns draws from it (as in Piironen and Vehtari's Fig. 3). With a slab, \eqn{1 - \kappa_j} is multiplied by \eqn{1 - b_j}, \eqn{b_j = 1/(1 + c^2/\sigma_j^2)} (their Sec. 3.4).
}
\value{
\code{horseshoe_scale}: \eqn{\tau_0}, with attributes \code{sigma}, \code{p0} and \code{tau0_pooled}. \code{horseshoe_sigma}: the vector of \eqn{\sigma_j}, in the order of the model's parameters (the index matrix). \code{horseshoe_prop_zero}: \code{n} draws of the prior proportion of zero rates.
}
\references{
Piironen, J., and A. Vehtari. 2017. Sparsity information and regularization in the horseshoe and other shrinkage priors. Electronic Journal of Statistics 11:5018-5051.
}
\seealso{\code{\link{corHMM}}}
\examples{
## likelihood scales of 12 rates, some informed by the data and one not
sigma <- c(10^seq(-4, -2, length.out = 11), Inf)
tau0 <- horseshoe_scale(prop_zero = 0.5, sigma = sigma)
tau0
summary(horseshoe_prop_zero(tau0, sigma))
\dontrun{
fit <- corHMM(phy, data, rate.cat = 1, rate.mat = rm, use_RTMB = TRUE,
              return.devfun = TRUE)
tau0 <- horseshoe_scale(fit, prop_zero = 0.5)
fit_hs <- corHMM(phy, data, rate.cat = 1, rate.mat = rm, use_RTMB = TRUE,
                 prior = list(family = "horseshoe", scale_global = tau0))
}
}
\keyword{models}
