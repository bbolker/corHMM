## Choosing the global scale of the horseshoe prior from a prior guess of the
## proportion of rates that are (effectively) zero, following Piironen and
## Vehtari (2017, Sec. 3.3).
##
## With a Gaussian approximation to the likelihood of each rate (standard
## error sigma_j, rates uncorrelated), the shrinkage factor of rate j is
## kappa_j = 1 / (1 + tau^2 lambda_j^2 / sigma_j^2), and with lambda_j ~
## C+(0, 1), E(kappa_j | tau) = 1 / (1 + tau / sigma_j). The expected number
## of rates distinguishable from zero is m_eff = sum_j (1 - kappa_j), so
## E(m_eff | tau) = sum_j 1 / (1 + sigma_j / tau). Setting it to
## p0 = (1 - prop_zero) D gives tau0; with equal sigma_j this is their
## Eq. 3.13, tau0 = p0 / (D - p0) sigma.

horseshoe_scale <- function(fit = NULL, prop_zero, sigma = NULL, lower = 1e-9, upper = 1e8) {
    if (!is.numeric(prop_zero) || length(prop_zero) != 1 || !(prop_zero > 0 && prop_zero < 1)) {
        stop("`prop_zero` must be a number strictly between 0 and 1.", call. = FALSE)
    }
    if (is.null(sigma)) sigma <- horseshoe_sigma(fit, lower = lower, upper = upper)
    if (!is.numeric(sigma) || any(is.na(sigma)) || any(sigma <= 0)) {
        stop("`sigma` must be positive numbers (Inf allowed).", call. = FALSE)
    }
    D <- length(sigma)
    p0 <- (1 - prop_zero) * D
    ## rates with sigma = Inf can't be distinguished from zero at any tau
    D_fin <- sum(is.finite(sigma))
    if (p0 >= D_fin) {
        stop(sprintf(paste("only %d of the %d rates are informed by the data (finite sigma),",
                           "so `prop_zero` must exceed %.3g."), D_fin, D, 1 - D_fin / D),
             call. = FALSE)
    }
    s <- sigma[is.finite(sigma)]
    f <- function(log_tau) sum(1 / (1 + s / exp(log_tau))) - p0
    rng <- log(range(s)) + c(-40, 40)
    log_tau0 <- stats::uniroot(f, rng, tol = 1e-10)$root
    structure(exp(log_tau0), sigma = sigma, p0 = p0,
              tau0_pooled = p0 / (D - p0) * exp(mean(log(s))))
}

## Per-rate likelihood scale sigma_j on the rate scale: the half-width of the
## interval over which the log-likelihood, as a function of rate j alone
## (the other rates at their ML values), is within 1/2 of its maximum. For
## a well-estimated rate this is about its standard error; for a rate whose
## estimate is at (or near) zero, it is the rate at which the data first
## distinguish it from zero. Inf if the log-likelihood never drops by 1/2
## within `upper`.
horseshoe_sigma <- function(fit, lower = 1e-9, upper = 1e8) {
    if (!inherits(fit, "corhmm") || !is.function(fit$devfun)) {
        stop("`fit` must be a corHMM fit with `use_RTMB = TRUE, return.devfun = TRUE`.",
             call. = FALSE)
    }
    if (isTRUE(fit$formula_model)) {
        stop("`fit` must be a native (rate.mat) model: for a formula model, fit the ",
             "equivalent native model (e.g. with equateStateMatPars()).", call. = FALSE)
    }
    if (!is.null(fit$prior)) {
        stop("`fit` must be a maximum-likelihood fit (no `prior`).", call. = FALSE)
    }
    obj <- fit$devfun()
    x_hat <- unname(fit$args.list$p)
    nll_hat <- obj$fn(x_hat)
    ## log-likelihood drop (minus 1/2) at log-rate y for rate j
    drop <- function(y, j) {
        x <- x_hat
        x[j] <- y
        obj$fn(x) - nll_hat - 0.5
    }
    vapply(seq_along(x_hat), function(j) {
        xj <- x_hat[j]
        q_up <- if (drop(log(upper), j) < 0) Inf else
            exp(stats::uniroot(drop, c(xj, log(upper)), j = j, tol = 1e-10)$root)
        if (!is.finite(q_up)) return(Inf)
        if (xj > log(lower) + 1e-8 && drop(log(lower), j) > 0) {
            q_low <- exp(stats::uniroot(drop, c(log(lower), xj), j = j, tol = 1e-10)$root)
            (q_up - q_low) / 2
        } else {
            q_up - exp(xj)
        }
    }, numeric(1))
}

## Simulated prior distribution of the effective proportion of zero rates,
## 1 - m_eff / D = mean_j kappa_j, for a horseshoe with global scale tau0 and
## per-rate likelihood scales sigma (Piironen and Vehtari 2017, Fig. 3). With
## a slab, 1 - kappa_j is multiplied by 1 - b_j, b_j = 1 / (1 + c^2 / sigma_j^2)
## (their Sec. 3.4).
horseshoe_prop_zero <- function(tau0, sigma, fix_global = FALSE, slab_scale = Inf, n = 1e4) {
    D <- length(sigma)
    tau <- if (fix_global) rep(tau0, n) else abs(tau0 * stats::rcauchy(n))
    lambda <- matrix(abs(stats::rcauchy(n * D)), n, D)
    a2 <- (tau * lambda)^2 / rep(sigma^2, each = n)
    kappa <- 1 / (1 + a2)
    ## (no slab: b = 0; a rate with sigma = Inf is fully shrunk either way)
    b <- if (is.finite(slab_scale)) 1 / (1 + slab_scale^2 / sigma^2) else rep(0, D)
    rowMeans(1 - (1 - kappa) * rep(1 - b, each = n))
}

## Choosing the global scale from the prior alone: a rate counts as zero if
## it is below `threshold` (by default 1/L, L the total branch length: less
## than about one expected change on the tree), and tau0 is set so that the
## prior probability of that is prop_zero. Under the horseshoe
## q = s |z|, z ~ N(0, 1), with s = tau lambda (slab: s = c tau lambda /
## sqrt(c^2 + tau^2 lambda^2)), so P(q < eps) depends on eps / tau0 (and
## c / tau0) only, and is decreasing in tau0.
horseshoe_scale_prior <- function(prop_zero, phy = NULL, threshold = 1 / sum(phy$edge.length),
                                  fix_global = FALSE, slab_scale = Inf) {
    if (!is.numeric(prop_zero) || length(prop_zero) != 1 || !(prop_zero > 0 && prop_zero < 1)) {
        stop("`prop_zero` must be a number strictly between 0 and 1.", call. = FALSE)
    }
    if (missing(threshold) && is.null(phy)) {
        stop("give `phy` (for the default threshold, 1/total branch length) or `threshold`.",
             call. = FALSE)
    }
    if (!is.numeric(threshold) || length(threshold) != 1 || !(threshold > 0)) {
        stop("`threshold` must be a positive number.", call. = FALSE)
    }
    f <- function(log_tau0) {
        hs_prob_below(log(threshold) - log_tau0, log(slab_scale) - log_tau0, fix_global) - prop_zero
    }
    rng <- log(threshold) + c(-60, 60)
    log_tau0 <- stats::uniroot(f, rng, tol = 1e-10)$root
    structure(exp(log_tau0), threshold = threshold, prop_zero = prop_zero,
              fix_global = fix_global)
}

## P(q < u) for the (regularized) half-horseshoe with tau0 = 1 and slab scale
## r = exp(log_r), at log u: with tau fixed at 1, E_lambda[2 Phi(u / s) - 1];
## with tau ~ C+(0, 1), that averaged over tau. Both by Gauss-Legendre
## quadrature over log lambda (and log tau), where the half-Cauchy density
## (2 / pi) e^w / (1 + e^(2 w)) decays like e^-|w|, so +/- 40 suffices.
hs_prob_below <- function(log_u, log_r = Inf, fix_global = TRUE) {
    gl <- gauss_legendre_panels(seq(-40, 40, by = 1), 10)
    log_dens <- log(gl$w) + log(2 / pi) + gl$x - log_add(0, 2 * gl$x)
    ## P(q < u | tau = e^t) for a vector t
    given_tau <- function(t) {
        vapply(t, function(tt) {
            w <- gl$x + tt  ## log of the scale tau lambda before the slab
            log_s <- if (is.finite(log_r)) w + log_r - 0.5 * log_add(2 * log_r, 2 * w) else w
            sum(exp(log_dens) * (1 - 2 * stats::pnorm(-exp(log_u - log_s))))
        }, numeric(1))
    }
    if (fix_global) return(given_tau(0))
    sum(exp(log_dens) * given_tau(gl$x))
}

## Simulated prior proportion of rates below `threshold`, among n_rates rates
## sharing the global scale (fix_global = FALSE: tau ~ C+(0, tau0) per draw)
horseshoe_prior_zero <- function(tau0, n_rates, threshold, fix_global = FALSE,
                                 slab_scale = Inf, n = 1e4) {
    tau <- if (fix_global) rep(tau0, n) else abs(tau0 * stats::rcauchy(n))
    s <- tau * matrix(abs(stats::rcauchy(n * n_rates)), n, n_rates)
    if (is.finite(slab_scale)) s <- slab_scale * s / sqrt(slab_scale^2 + s^2)
    rowMeans(s * abs(matrix(stats::rnorm(n * n_rates), n, n_rates)) < threshold)
}
