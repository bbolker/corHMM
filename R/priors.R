## Priors on transition rates for the RTMB objective.
##
## The model parameters are log-rates x = log(q) (or, for formula models,
## coefficients mapped linearly to log-rates). Every family except "lognormal"
## is specified on the rate scale, as in BayesTraits, and its density is
## transformed to the log scale (Jacobian: + x), so the objective is the
## negative log posterior density of the log-rates -- the density that HMC
## samples and whose mode corHMM reports.

prior_families <- c("lognormal", "gamma", "exponential", "uniform",
                     "hyper_gamma", "hyper_exp", "horseshoe")

## fill in defaults and check a user-supplied `prior` list (see ?corHMM)
check_prior <- function(prior, lower.bound, upper.bound, model_is_formula) {
    if (!is.list(prior)) stop("`prior` must be a list.", call. = FALSE)
    if (is.null(prior$family)) prior$family <- "lognormal"
    prior$family <- match.arg(prior$family, prior_families)
    if (is.null(prior$target)) prior$target <- "rates"
    prior$target <- match.arg(prior$target, c("rates", "coefficients"))
    if (!model_is_formula && prior$target != "rates") {
        stop("`prior$target = 'coefficients'` is only available for formula models.",
             call. = FALSE)
    }
    if (prior$family != "lognormal" && prior$target != "rates") {
        stop("`prior$target = 'coefficients'` requires `prior$family = 'lognormal'`.",
             call. = FALSE)
    }

    need <- function(nm, len, lower = 0, strict = TRUE) {
        v <- prior[[nm]]
        ok <- is.numeric(v) && length(v) == len && all(is.finite(v)) &&
            all(if (strict) v > lower else v >= lower)
        if (!ok) {
            stop(sprintf("`prior$%s` must be %s %s number%s for family \"%s\".",
                         nm, if (len == 1) "a" else paste(len, ""),
                         if (strict) "positive" else "non-negative",
                         if (len == 1) "" else "s", prior$family),
                 call. = FALSE)
        }
    }
    need_range <- function(nm) {
        need(nm, 2, strict = FALSE)
        if (prior[[nm]][2] <= prior[[nm]][1]) {
            stop(sprintf("`prior$%s` must be an increasing pair c(lower, upper).", nm),
                 call. = FALSE)
        }
    }

    switch(prior$family,
           lognormal = {
               if (is.null(prior$lb)) prior$lb <- log(lower.bound)
               if (is.null(prior$ub)) prior$ub <- log(upper.bound)
               if (is.null(prior$range)) prior$range <- 3
               if (is.null(prior$coef.sd)) prior$coef.sd <- 2
           },
           gamma = { need("shape", 1); need("scale", 1) },
           exponential = need("mean", 1),
           uniform = {
               ## its limits apply to linear combinations of formula
               ## coefficients, which can't be box constraints, and HMC
               ## against a hard barrier is badly biased
               if (model_is_formula) {
                   stop("uniform priors are not available for formula models: ",
                        "use the native interface (`rate.mat`, e.g. built with ",
                        "`equateStateMatPars()`), which samples within the ",
                        "prior's limits.", call. = FALSE)
               }
               if (is.null(prior$min)) prior$min <- 0
               need("min", 1, strict = FALSE)
               need("max", 1)
               if (prior$max <= prior$min) {
                   stop("`prior$max` must be greater than `prior$min`.", call. = FALSE)
               }
           },
           hyper_gamma = { need_range("shape"); need_range("scale") },
           hyper_exp = need_range("mean"),
           horseshoe = {
               need("scale_global", 1)
               if (is.null(prior$slab_scale)) prior$slab_scale <- Inf
               if (!is.numeric(prior$slab_scale) || length(prior$slab_scale) != 1 ||
                   is.na(prior$slab_scale) || prior$slab_scale <= 0) {
                   stop("`prior$slab_scale` must be a positive number (or Inf) ",
                        "for family \"horseshoe\".", call. = FALSE)
               }
               if (is.null(prior$fix_global)) prior$fix_global <- FALSE
               if (!isTRUE(prior$fix_global) && !isFALSE(prior$fix_global)) {
                   stop("`prior$fix_global` must be TRUE or FALSE.", call. = FALSE)
               }
               ## the marginal over the global scale is tabulated over the
               ## range of rates that the optimizer and HMC can reach
               prior$rate_bounds <- c(lower.bound, upper.bound)
           })
    prior
}

## everything the objective needs to evaluate the prior, computed once
prior_setup <- function(prior) {
    ps <- list(family = prior$family)
    switch(prior$family,
           lognormal = {
               ps$mean <- (prior$lb + prior$ub) / 2
               ps$sd <- (prior$ub - prior$lb) / (2 * prior$range)
           },
           gamma = {
               ps$shape <- prior$shape
               ps$scale <- prior$scale
           },
           exponential = {
               ps$family <- "gamma"
               ps$shape <- 1
               ps$scale <- prior$mean
           },
           uniform = {
               qmin <- if (is.null(prior$min)) 0 else prior$min
               ps$log_width <- log(prior$max - qmin)
               ps$log_max <- log(prior$max)
               ps$log_min <- if (qmin > 0) log(qmin) else -Inf
           },
           hyper_gamma = ,
           hyper_exp = {
               hpar <- if (prior$family == "hyper_gamma") prior$shape else NULL
               hscale <- if (prior$family == "hyper_gamma") prior$scale else prior$mean
               tab <- hyper_log_tail(hpar)
               ps$logH <- RTMB::splinefun(tab$L, tab$logH, method = "natural")
               ps$log_b <- log(hscale)
               ps$log_norm <- log(diff(hscale)) +
                   if (is.null(hpar)) 0 else log(diff(hpar))
           },
           horseshoe = {
               ps <- c(ps, horseshoe_setup(prior))
           })
    ps
}

## negative log prior density of log-rates x (AD-safe)
prior_nll_rates <- function(x, ps) {
    "c" <- RTMB::ADoverload("c")
    n <- length(x)
    switch(ps$family,
           lognormal = sum(0.5 * log(2 * pi * ps$sd^2) +
                           0.5 * ((x - ps$mean) / ps$sd)^2),
           ## -log dgamma(exp(x)) - x
           gamma = sum(exp(x) / ps$scale - ps$shape * x) +
               n * (lgamma(ps$shape) + ps$shape * log(ps$scale)),
           uniform = {
               ## corHMM() turns (min, max) into box bounds for both
               ## optimization and HMC; this term only guards evaluations
               ## outside them (e.g. user-supplied `p`). Taped comparisons are
               ## not available, and TMBad folds 0 * (anything) to 0, so a
               ## 0 * log(.) barrier vanishes from the tape; max() is taped as a
               ## conditional, so this term is exactly 0 (zero gradient) inside
               ## and huge outside.
               excess <- if (is.finite(ps$log_min)) {
                   max(c(0, x - ps$log_max, ps$log_min - x))
               } else {
                   max(c(0, x - ps$log_max))
               }
               n * ps$log_width - sum(x) + 1e300 * excess
           },
           hyper_gamma = ,
           hyper_exp = {
               ## marginal density (H(x - log s1) - H(x - log s0)) / norm (see
               ## hyper_log_tail()); the second term is 0 when s0 = 0
               lp <- ps$logH(x - ps$log_b[2])
               if (ps$log_b[1] > -Inf) {
                   lp <- RTMB::logspace_sub(lp, ps$logH(x - ps$log_b[1]))
               }
               n * ps$log_norm - sum(lp) - sum(x)
           },
           horseshoe = {
               ## log of sum_k w_k prod_i p(q_i | tau_k) over the nodes of
               ## the global scale (one node, weight 1, if it is fixed)
               terms <- lapply(seq_along(ps$log_tau), function(k) {
                   g <- ps$log_g[[min(k, length(ps$log_g))]]
                   ps$log_w[k] + sum(g(x - ps$log_tau[k])) - n * ps$log_tau[k]
               })
               -Reduce(RTMB::logspace_add, terms) - sum(x)
           })
}

## Horseshoe prior on the rates (Carvalho et al. 2010), optionally
## regularized with a slab (Piironen & Vehtari 2017):
##   q_i ~ N+(0, s_i^2),  s_i^2 = c^2 tau^2 lambda_i^2 / (c^2 + tau^2 lambda_i^2),
##   lambda_i ~ C+(0, 1),  tau ~ C+(0, scale_global) (or tau = scale_global),
## with slab scale c (Inf: the plain horseshoe, s_i = tau lambda_i).
## Given tau, p(q | tau) = g(q / tau; c / tau) / tau, where g(u; r) is the
## density of u with tau = 1 and slab scale r. log g is tabulated as a
## function of log u by quadrature over log lambda and interpolated by an
## AD-safe spline: one table for the plain horseshoe, one per node of tau
## with a slab. A shared tau is integrated out by Gauss-Legendre quadrature
## over log tau, so the rates' prior is their exact joint marginal (up to
## quadrature error) and no hyperparameters are estimated. The log tau nodes
## cover the reachable rates (rate_bounds, from lower.bound and upper.bound)
## widened by 3 units, and scale_global +/- 8 units; the C+ density is
## renormalized over that range.
horseshoe_setup <- function(prior) {
    lb <- log(if (is.null(prior$rate_bounds)) c(1e-9, 100) else prior$rate_bounds)
    tab <- horseshoe_tables(lb, log(prior$scale_global),
                            log(if (is.null(prior$slab_scale)) Inf else prior$slab_scale),
                            isTRUE(prior$fix_global))
    list(log_tau = tab$log_tau, log_w = tab$log_w,
         log_g = lapply(tab$tables, function(t) RTMB::splinefun(t$v, t$log_g, method = "natural")))
}

## nodes and weights of log tau, and the tables of log g (v, log_g) for
## horseshoe_setup(); memoised, since tables with a slab take seconds
horseshoe_tables <- local({
    cache <- new.env(parent = emptyenv())
    function(lb, log_tau0, log_c, fix_global) {
        key <- paste(c(lb, log_tau0, log_c, fix_global), collapse = "_")
        if (!is.null(cache[[key]])) return(cache[[key]])
        if (fix_global) {
            log_tau <- log_tau0
            log_w <- 0
        } else {
            lo <- min(lb[1] - 3, log_tau0 - 8)
            hi <- max(lb[2] + 3, log_tau0 + 8)
            gl <- gauss_legendre_panels(seq(lo, hi, length.out = ceiling(hi - lo) + 1), 8)
            log_tau <- gl$x
            ## C+(0, tau0) density of log tau: (2/pi) tau tau0 / (tau0^2 + tau^2)
            log_w <- log(gl$w) + log(2 / pi) + log_tau + log_tau0 -
                log_add(2 * log_tau0, 2 * log_tau)
            log_w <- log_w - log_sum_exp_cols(matrix(log_w))
        }
        ## each table covers log u = log q - log tau for the reachable rates,
        ## with 5 units to spare (a natural spline extrapolates linearly)
        tables <- if (is.finite(log_c)) {
            lapply(log_tau, function(t) {
                v <- seq(lb[1] - t - 5, lb[2] - t + 5, by = 0.05)
                list(v = v, log_g = hs_log_g(v, log_c - t))
            })
        } else {
            v <- seq(lb[1] - max(log_tau) - 5, lb[2] - min(log_tau) + 5, by = 0.05)
            list(list(v = v, log_g = hs_log_g(v)))
        }
        res <- list(log_tau = log_tau, log_w = log_w, tables = tables)
        cache[[key]] <- res
        res
    }
})

## log(exp(a) + exp(b)), elementwise (numeric)
log_add <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))

## log density g(u; r) of the half-horseshoe with global scale 1 and slab
## scale r = exp(log_r) (Inf: no slab) at u = exp(v), by composite
## Gauss-Legendre quadrature over w = log lambda. The integrand vanishes
## (super-exponentially) for lambda below about u e^-8 and decays like
## 1/lambda^2 (no slab) or 1/lambda (slab, once the scale reaches r) above
## max(u, 1, r), so the range is cut off there.
hs_log_g <- function(v, log_r = Inf) {
    slab <- is.finite(log_r)
    lo <- v - 8
    hi <- pmax(v, 0, if (slab) log_r else 0) + if (slab) 40 else 18
    gl <- gauss_legendre_panels(seq(0, 1, length.out = 33), 12)
    ## nodes (rows) for each v (columns), mapped from [0, 1] to [lo, hi]
    w <- outer(gl$x, hi - lo) + rep(lo, each = length(gl$x))
    log_dw <- log(outer(gl$w, hi - lo))
    log_sd <- if (slab) w + log_r - 0.5 * log_add(2 * log_r, 2 * w) else w
    z2 <- exp(2 * (rep(v, each = nrow(w)) - log_sd))
    ## 2 dnorm(u, 0, sd) * (2 / pi) / (1 + lambda^2) * lambda
    log_f <- log(2) - 0.5 * log(2 * pi) - log_sd - z2 / 2 + log(2 / pi) + w -
        log_add(0, 2 * w)
    log_sum_exp_cols(log_dw + log_f)
}

## Gauss-Legendre nodes and weights on [a, b] (Golub-Welsch)
gauss_legendre <- function(n, a = -1, b = 1) {
    i <- seq_len(n - 1)
    J <- matrix(0, n, n)
    J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- i / sqrt(4 * i^2 - 1)
    e <- eigen(J, symmetric = TRUE)
    o <- order(e$values)
    list(x = (b - a) / 2 * e$values[o] + (a + b) / 2,
         w = (b - a) * e$vectors[1, o]^2)
}

## composite Gauss-Legendre rule over the panels defined by `breaks`
gauss_legendre_panels <- function(breaks, n = 20) {
    g <- lapply(seq_len(length(breaks) - 1),
                function(i) gauss_legendre(n, breaks[i], breaks[i + 1]))
    list(x = unlist(lapply(g, `[[`, "x")), w = unlist(lapply(g, `[[`, "w")))
}

log_sum_exp_cols <- function(m) {
    mx <- apply(m, 2, max)
    mx + log(colSums(exp(m - rep(mx, each = nrow(m)))))
}

## Tabulated log H(L), the building block of the hyperprior marginals.
##
## With q a rate and w = log(q / s) for hyper-scale s:
##  - gamma prior, shape a ~ U(a0, a1), scale s ~ U(s0, s1):
##      p(q) = [H(log q - log s1) - H(log q - log s0)] / ((a1 - a0) (s1 - s0)),
##      H(L) = int_L^Inf h(w) dw,  h(w) = exp(-e^w - w) int_{a0}^{a1} exp(a w) / Gamma(a) da
##  - exponential prior, mean s ~ U(s0, s1): the same with h(w) = exp(-e^w),
##      so H(L) = E1(exp(L)).
## H depends only on the shape range, so the table is memoised.
## As L -> -Inf, log H is asymptotically linear in L, so a natural spline
## extrapolates it sensibly below the table.
hyper_log_tail <- local({
    cache <- new.env(parent = emptyenv())
    function(shape_range = NULL, L = seq(-100, 8, by = 0.05)) {
        key <- paste(c(shape_range, range(L), length(L)), collapse = "_")
        if (!is.null(cache[[key]])) return(cache[[key]])

        ## nodes within each L interval, plus a tail interval above max(L)
        edges <- c(L, max(L) + 2)
        gl <- gauss_legendre(10)
        nint <- length(edges) - 1
        half <- diff(edges) / 2
        mid <- (edges[-1] + edges[-length(edges)]) / 2
        w <- outer(gl$x, half) + rep(mid, each = length(gl$x))
        logwt <- log(outer(gl$w / 2, 2 * half))

        logh <- -exp(w)
        if (!is.null(shape_range)) {
            a0 <- shape_range[1]
            a1 <- shape_range[2]
            ## panels refined towards the lower end, where exp(a w) / Gamma(a)
            ## changes fastest in a when w is very negative
            br <- if (a0 == 0) c(0, a1 * 10^(-6:0)) else
                exp(seq(log(a0), log(a1), length.out = 8))
            ga <- gauss_legendre_panels(br, 20)
            logg <- log_sum_exp_cols(outer(log(ga$w) - lgamma(ga$x), rep(1, length(w))) +
                                     outer(ga$x, as.vector(w)))
            logh <- logh - w + matrix(logg, nrow(w))
        }
        ## log of the integral of h over each interval
        logI <- log_sum_exp_cols(logwt + logh)
        ## H(L_i) = sum of interval integrals from i upwards (log scale)
        logH <- numeric(nint)
        acc <- -Inf
        for (i in rev(seq_len(nint))) {
            m <- max(acc, logI[i])
            acc <- m + log(exp(acc - m) + exp(logI[i] - m))
            logH[i] <- acc
        }
        res <- list(L = L, logH = logH)
        cache[[key]] <- res
        res
    }
})
