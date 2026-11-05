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
