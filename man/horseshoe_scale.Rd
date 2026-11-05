\name{horseshoe_scale}
\alias{horseshoe_scale}
\alias{horseshoe_sigma}
\alias{horseshoe_prop_zero}
\alias{horseshoe_scale_prior}
\alias{horseshoe_prior_zero}
\title{Global scale of a horseshoe prior from a proportion of zero rates}
\description{
Choose \code{scale_global} (\eqn{\tau_0}) for \code{prior = list(family = "horseshoe")} in \code{\link{corHMM}} from a guess of the proportion of rates that are zero, and simulate the proportion of zero rates that a choice of \eqn{\tau_0} implies. Since no rate is exactly zero under the horseshoe, "zero" needs a definition. There are two:
\itemize{
  \item \code{horseshoe_scale}: a rate is zero if the prior shrinks it below what the data can distinguish from zero (the effective number of non-zero rates of Piironen and Vehtari 2017, Sec. 3.3). This depends on the data, through a maximum-likelihood fit.
  \item \code{horseshoe_scale_prior}: a rate is zero if it is below \code{threshold}, by default \eqn{1/L} with \eqn{L} the total branch length (less than about one expected change on the tree). This uses the prior alone, like the prior probability of a zero rate in reversible-jump MCMC.
}
}
\usage{
horseshoe_scale(fit = NULL, prop_zero, sigma = NULL, lower = 1e-9, upper = 1e8)
horseshoe_sigma(fit, lower = 1e-9, upper = 1e8)
horseshoe_prop_zero(tau0, sigma, fix_global = FALSE, slab_scale = Inf, n = 1e4)
horseshoe_scale_prior(prop_zero, phy = NULL, threshold = 1 / sum(phy$edge.length),
                      fix_global = FALSE, slab_scale = Inf)
horseshoe_prior_zero(tau0, n_rates, threshold, fix_global = FALSE, slab_scale = Inf,
                     n = 1e4)
}
\arguments{
\item{fit}{a maximum-likelihood \code{corHMM} fit (no \code{prior}) of a native model (specified with \code{rate.mat}, not a formula), with \code{use_RTMB = TRUE} and \code{return.devfun = TRUE}. Not used by \code{horseshoe_scale} if \code{sigma} is given.}
\item{prop_zero}{the guess for the proportion of rates that are zero (for \code{horseshoe_scale_prior}, the prior probability that a rate is below \code{threshold}), strictly between 0 and 1.}
\item{sigma}{the likelihood scale of each rate (see Details); by default computed by \code{horseshoe_sigma(fit)}.}
\item{lower, upper}{the range of rates searched for the likelihood scales (use the \code{lower.bound} of the fit and at least its \code{upper.bound}).}
\item{tau0}{the global scale \eqn{\tau_0} (\code{scale_global}).}
\item{fix_global, slab_scale}{as in the \code{"horseshoe"} prior (see \code{\link{corHMM}}): if \code{fix_global = FALSE}, \eqn{\tau \sim C^+(0, \tau_0)}, otherwise \eqn{\tau = \tau_0}; \code{slab_scale} is the slab scale \eqn{c} (\code{Inf}: no slab).}
\item{phy}{the tree, for the default \code{threshold}.}
\item{threshold}{the rate below which a rate counts as zero.}
\item{n_rates}{the number of rates sharing the global scale.}
\item{n}{the number of simulated draws.}
}
\details{
With a Gaussian approximation to the likelihood of each rate \eqn{q_j}, with scale \eqn{\sigma_j} and the rates treated as uncorrelated, the shrinkage factor of rate \eqn{j} under the horseshoe is \eqn{\kappa_j = 1/(1 + \tau^2\lambda_j^2/\sigma_j^2)} (0: not shrunk; 1: shrunk to zero), and \eqn{E(\kappa_j \mid \tau) = 1/(1 + \tau/\sigma_j)}. The effective number of non-zero rates is \eqn{m_{eff} = \sum_j (1 - \kappa_j)}, so \eqn{E(m_{eff} \mid \tau) = \sum_j 1/(1 + \sigma_j/\tau)}. \code{horseshoe_scale} solves \eqn{E(m_{eff} \mid \tau_0) = p_0} with \eqn{p_0 = (1 - } \code{prop_zero}\eqn{) D} for \eqn{D} rates. With equal \eqn{\sigma_j = \sigma} this is Piironen and Vehtari's Eq. 3.13, \eqn{\tau_0 = \sigma p_0/(D - p_0)}, which is returned as the attribute \code{tau0_pooled} (with \eqn{\sigma} the geometric mean of the finite \eqn{\sigma_j}).

\code{horseshoe_sigma} measures \eqn{\sigma_j} on the rate scale, as the half-width of the interval over which the log-likelihood, as a function of rate \eqn{j} with the other rates at their maximum-likelihood values, is within 1/2 of its maximum. For a well-estimated rate this is about its standard error; for a rate whose estimate is at or near zero (where the Hessian of the log-likelihood is degenerate), it is the rate at which the data first distinguish it from zero. A rate whose log-likelihood never drops by 1/2 below \code{upper} gets \eqn{\sigma_j = \infty}: it is shrunk whatever \eqn{\tau}, so the proportion of such rates is a lower limit for \code{prop_zero}.

With \code{fix_global = TRUE}, the prior mean of the proportion of zero rates, \eqn{1 - m_{eff}/D}, is \code{prop_zero}. With \eqn{\tau \sim C^+(0, \tau_0)} (the default of the \code{"horseshoe"} prior) its median is close to \code{prop_zero}, but the distribution is wider; \code{horseshoe_prop_zero} returns draws from it (as in Piironen and Vehtari's Fig. 3). With a slab, \eqn{1 - \kappa_j} is multiplied by \eqn{1 - b_j}, \eqn{b_j = 1/(1 + c^2/\sigma_j^2)} (their Sec. 3.4).

Under the horseshoe, \eqn{q = s|z|} with \eqn{z \sim N(0, 1)} and \eqn{s = \tau\lambda} (with a slab, \eqn{s = c\tau\lambda/\sqrt{c^2 + \tau^2\lambda^2}}), so the prior probability \eqn{P(q < \epsilon)} depends only on \eqn{\epsilon/\tau_0} (and \eqn{c/\tau_0}). \code{horseshoe_scale_prior} computes it by quadrature over \eqn{\lambda} (and \eqn{\tau}, unless \code{fix_global = TRUE}) and solves \eqn{P(q < \epsilon) =} \code{prop_zero} for \eqn{\tau_0}. With a fixed global scale the rates are independent a priori, so the proportion of zero rates among \eqn{D} rates is binomial with mean \code{prop_zero}; with \eqn{\tau \sim C^+(0, \tau_0)} the rates are scaled together and the proportion is much more variable (from near 0 to 1), which \code{horseshoe_prior_zero} shows. The binomial spread is closer to that of the reversible-jump prior in BayesTraits.
}
\value{
\code{horseshoe_scale}: \eqn{\tau_0}, with attributes \code{sigma}, \code{p0} and \code{tau0_pooled}. \code{horseshoe_sigma}: the vector of \eqn{\sigma_j}, in the order of the model's parameters (the index matrix). \code{horseshoe_prop_zero}, \code{horseshoe_prior_zero}: \code{n} draws of the prior proportion of zero rates. \code{horseshoe_scale_prior}: \eqn{\tau_0}, with attributes \code{threshold}, \code{prop_zero} and \code{fix_global}.
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
## prior alone: P(rate < 1/L) = 0.5, for a tree of total length L = 100
tau0 <- horseshoe_scale_prior(0.5, threshold = 1 / 100, fix_global = TRUE)
tau0
table(horseshoe_prior_zero(tau0, n_rates = 8, threshold = 1 / 100, fix_global = TRUE))
\dontrun{
fit <- corHMM(phy, data, rate.cat = 1, rate.mat = rm, use_RTMB = TRUE,
              return.devfun = TRUE)
tau0 <- horseshoe_scale(fit, prop_zero = 0.5)
fit_hs <- corHMM(phy, data, rate.cat = 1, rate.mat = rm, use_RTMB = TRUE,
                 prior = list(family = "horseshoe", scale_global = tau0))
}
}
\keyword{models}
