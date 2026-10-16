translate <- function(formula_list, nstate = 2) {
  if (inherits(formula_list, "formula")) {
    formula_list <- list(formula_list)
  }

  specs <- vector("list", length(formula_list))

  for (i in seq_along(formula_list)) {
    f <- formula_list[[i]]
    lhs <- f[[2]]
    rhs <- f[[3]]
    trait <- deparse(lhs)

    if (is.call(rhs) && identical(rhs[[1]], quote(symm))) {
      rhs <- rhs[[2]]
      f2 <- stats::as.formula(bquote(.(lhs) ~ .(rhs)))

      specs[[i]] <- list(
        trait = trait,
        symmetric = TRUE,
        gain_formula = f2,
        loss_formula = f2
      )
    } else {
      if (is.call(rhs) && identical(rhs[[1]], quote(list))) {
        stop(
          "different gain/loss formulas are not currently supported.",
          call. = FALSE
        )
      }

      specs[[i]] <- list(
        trait = trait,
        symmetric = FALSE,
        gain_formula = f,
        loss_formula = f
      )
    }
  }

  trait_names <- vapply(specs, `[[`, character(1), "trait")
  names(formula_list) <- trait_names
  names(specs) <- trait_names

  stateList <- rep.int(as.integer(nstate), length(trait_names))
  names(stateList) <- trait_names

  traitMatrix <- do.call(
    expand.grid,
    lapply(stateList, function(x) 0:(x - 1))
  )
  traitMatrix$label <- do.call(
    paste,
    c(traitMatrix[trait_names], sep = "|")
  )

  ns <- nrow(traitMatrix)
  Xstate <- as.matrix(traitMatrix[, trait_names, drop = FALSE])

  idx_grid <- expand.grid(
    to = seq_len(ns),
    from = seq_len(ns)
  )
  idx_grid <- idx_grid[idx_grid$from != idx_grid$to, , drop = FALSE]

  diff_matrix <- Xstate[idx_grid$to, , drop = FALSE] -
    Xstate[idx_grid$from, , drop = FALSE]

  one_step <- which(rowSums(diff_matrix != 0) == 1L)
  edges <- idx_grid[one_step, , drop = FALSE]
  final_diffs <- diff_matrix[one_step, , drop = FALSE]
  change_col <- max.col(final_diffs != 0)

  edge_tab <- data.frame(
    edge_id = seq_along(one_step),
    from = edges$from,
    to = edges$to,
    from_label = traitMatrix$label[edges$from],
    to_label = traitMatrix$label[edges$to],
    changed_trait = trait_names[change_col],
    delta = final_diffs[cbind(seq_along(change_col), change_col)]
  )

  edge_tab$direction <- ifelse(edge_tab$delta > 0, "gain", "loss")
  edge_tab$focal_from <- Xstate[cbind(edges$from, change_col)]
  edge_tab$focal_to <- Xstate[cbind(edges$to, change_col)]

  edge_tab <- cbind(
    edge_tab,
    traitMatrix[edges$from, trait_names, drop = FALSE]
  )

  edge_tab <- edge_tab[order(edge_tab$to, edge_tab$from), , drop = FALSE]
  edge_tab$edge_id <- seq_len(nrow(edge_tab))

  blocks <- list()
  shared_symm <- list()
  par_counter <- 1L

  edge_groups <- split(
    edge_tab,
    paste(edge_tab$changed_trait, edge_tab$direction, sep = "_")
  )
  edge_groups <- edge_groups[sort(names(edge_groups))]

  for (block_name in names(edge_groups)) {
    dat <- edge_groups[[block_name]]
    tr <- as.character(dat$changed_trait[1])
    dir <- as.character(dat$direction[1])
    spec <- specs[[tr]]

    curr_formula <- if (dir == "gain") {
      spec$gain_formula
    } else {
      spec$loss_formula
    }

    rhs_terms <- delete.response(stats::terms(curr_formula))
    X <- stats::model.matrix(rhs_terms, data = dat)
    n_col <- ncol(X)

    if (spec$symmetric) {
      shared_name <- paste0(tr, "_symm")

      if (is.null(shared_symm[[shared_name]])) {
        par_index <- par_counter + seq_len(n_col) - 1L
        coef_names <- paste0(shared_name, ":", colnames(X))

        shared_symm[[shared_name]] <- list(
          par_index = par_index,
          coef_names = coef_names
        )

        par_counter <- par_counter + n_col
      } else {
        par_index <- shared_symm[[shared_name]]$par_index
        coef_names <- shared_symm[[shared_name]]$coef_names
      }
    } else {
      par_index <- par_counter + seq_len(n_col) - 1L
      coef_names <- paste0(block_name, ":", colnames(X))
      par_counter <- par_counter + n_col
    }

    blocks[[block_name]] <- list(
      block_name = block_name,
      trait = tr,
      direction = dir,
      formula = curr_formula,
      rhs_terms = rhs_terms,
      edge_id = dat$edge_id,
      from = dat$from,
      to = dat$to,
      from_label = dat$from_label,
      to_label = dat$to_label,
      X = X,
      n_par = n_col,
      coef_names = coef_names,
      par_index = par_index
    )
  }

  Q_indicator <- matrix(0L, ns, ns)
  Q_indicator[as.matrix(edge_tab[, c("from", "to")])] <- edge_tab$edge_id
  rownames(Q_indicator) <- colnames(Q_indicator) <- traitMatrix$label

  Q0 <- Matrix::Matrix(
    0,
    ns,
    ns,
    sparse = TRUE,
    dimnames = list(traitMatrix$label, traitMatrix$label)
  )

  n_par <- par_counter - 1L
  par_names <- character(n_par)

  for (b in blocks) {
    par_names[b$par_index] <- b$coef_names
  }

  rate_X <- matrix(
    0,
    nrow = nrow(edge_tab),
    ncol = n_par,
    dimnames = list(
      paste(edge_tab$from_label, edge_tab$to_label, sep = "->"),
      par_names
    )
  )

  for (b in blocks) {
    rate_X[b$edge_id, b$par_index] <- b$X
  }

  list(
    d = ns,
    trait_names = trait_names,
    stateList = stateList,
    formulas = formula_list,
    traitMatrix = traitMatrix,
    edge_table = edge_tab,
    blocks = blocks,
    Q_indicator = Q_indicator,
    Q0 = Q0,
    n_par = n_par,
    par_names = par_names,
    par_index = lapply(blocks, `[[`, "par_index"),
    rate_X = rate_X
  )
}


## Group the edges for level-wise pruning. A node's level is 1 + the maximum
## level of its children (tips are level 0), so all edges whose parent is at
## the same level can be processed together once the lower levels are done.
## For each level: the edge rows, child nodes, edge lengths and parent nodes,
## plus, for the k-th child of each parent (`slots`), which rows of the
## per-edge result go to which parent.
prune_levels <- function(phy) {
  nb.tip <- Ntip(phy)
  e <- phy$edge
  lev <- integer(nb.tip + Nnode(phy))
  for (r in reorder(phy, "postorder", index.only = TRUE)) {
    lev[e[r, 1]] <- max(lev[e[r, 1]], lev[e[r, 2]] + 1L)
  }
  lapply(sort(unique(lev[-seq_len(nb.tip)])), function(L) {
    rows <- which(lev[e[, 1]] == L)
    rows <- rows[order(e[rows, 1])]
    par <- e[rows, 1]
    nodes <- unique(par)
    slot <- stats::ave(seq_along(par), par, FUN = seq_along)
    list(
      child = e[rows, 2],
      t = phy$edge.length[rows],
      nodes = nodes,
      slots = lapply(seq_len(max(slot)), function(k) {
        list(src = which(slot == k), dst = match(par[slot == k], nodes))
      })
    )
  })
}

## TRUE if Q has (nearly) repeated eigenvalues or (nearly) defective
## eigenvectors, where derivatives through the eigendecomposition are unreliable
eigen_degenerate <- function(Q, gap_tol = 1e-6, kappa_tol = 1e8) {
  if (any(!is.finite(Q))) return(TRUE)
  ev <- eigen(Q, symmetric = FALSE)
  lam <- ev$values
  scale <- max(Mod(lam), .Machine$double.xmin)
  gaps <- Mod(outer(lam, lam, "-"))[upper.tri(diag(length(lam)))]
  any(!is.finite(lam)) ||
    min(gaps) / scale < gap_tol ||
    kappa(ev$vectors, exact = TRUE) > kappa_tol
}

## settings for the RTMB back end of corHMM() (argument RTMB_control):
## user-supplied values override the defaults; unknown names are an error
corHMM_RTMB_control <- function(control = list()) {
  defaults <- list(pruning = "auto")
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("`RTMB_control` must be a list.", call. = FALSE)
  bad <- setdiff(names(control), names(defaults))
  if (length(bad) || (length(control) && is.null(names(control)))) {
    stop("unknown `RTMB_control` element(s): ",
         paste(if (length(bad)) bad else "(unnamed)", collapse = ", "),
         "; valid names: ", paste(names(defaults), collapse = ", "),
         call. = FALSE)
  }
  control <- utils::modifyList(defaults, control)
  control$pruning <- match.arg(control$pruning, c("auto", "eigen", "expm"))
  control
}

## pruning = "eigen" uses one eigendecomposition of Q per evaluation, with
## all edges at the same tree level handled together; "expm" uses a matrix
## exponential per edge. "auto" (default) uses "expm" throughout if Q has
## repeated eigenvalues for any parameter values (e.g. ER models with > 2
## states). Otherwise it returns an object whose tape (used directly by
## tmbstan, sdreport, $he) is "expm", but whose $fn, $gr and $report use
## "eigen", falling back to "expm" where eigenvalues (nearly) coincide (e.g.
## equal starting values) or the eigen tape fails. So optimization gets the
## speedup and HMC stays safe; pruning = "eigen" opts in for HMC as well.
mkdev.corhmm_rtmb <- function(p, phy, liks, Q, rate, root.p, rate.cat,
                              order.test, lewis.asc.bias,
                              set.fog = FALSE, fog.vec = NULL,
                              pen.type = NULL, lambda = 1,
                              q_prep = NULL, prior = NULL,
                              pruning = c("auto", "eigen", "expm")) {
  pruning <- match.arg(pruning)
  formula_interface <- !is.null(q_prep)
  n_fog <- if (set.fog) length(unique(fog.vec)) else 0L

  names(p) <- c(
    if (n_fog) paste0("fog", seq_len(n_fog)),
    if (formula_interface) {
      q_prep$par_names
    } else {
      paste0("rate", seq_len(length(p) - n_fog))
    }
  )

  if (lewis.asc.bias) {
    stop("can't do lewis.asc.bias yet; recursive call to likelihood function ...")
  }

  if (order.test) {
    stop("can't do order.test (non-differentiable ...)")
  }

  if (set.fog) {
    warning("set.fog is untested in RTMB implementation", call. = FALSE)
  }

  if (formula_interface && set.fog) {
    warning("tip.fog with the formula interface is untested.", call. = FALSE)
  }

  if (formula_interface && !is.null(pen.type)) {
    warning("pen.type with the formula interface is untested.", call. = FALSE)
  }

  if (formula_interface && rate.cat != 1) {
    stop("formula models currently require `rate.cat = 1`.", call. = FALSE)
  }

  nb.node <- Nnode(phy)
  nb.tip <- Ntip(phy)
  TIPS <- seq.int(nb.tip)
  anc <- unique(phy$edge[,1])
  k.rates <- dim(Q)[2] / 2

  pen_indices <- if (!formula_interface && !is.null(pen.type) && rate.cat > 1) {
    rate_class_names <- paste0("R", 1:rate.cat)
    lapply(rate_class_names, function(rc) grep(rc, colnames(rate)))
  } else {
    NULL
  }

  nq <- nrow(Q)

  fog_idx <- if (set.fog) {
    lapply(seq_len(nb.tip), function(tip.index) {
      list(
        zeros = which(liks[tip.index, ] == 0),
        ones  = which(liks[tip.index, ] == 1),
        not1  = which(liks[tip.index, ] != 1)
      )
    })
  } else {
    NULL
  }

  has_prior <- !is.null(prior)
  prior_on_rates <- FALSE
  prior_on_coefficients <- FALSE
  prior.mean <- prior.sd <- coef.sd <- 0
  intercept_idx <- slope_idx <- integer(0)

  if (has_prior) {
    prior.mean <- (prior$lb + prior$ub) / 2
    prior.sd <- (prior$ub - prior$lb) / (2 * prior$range)
    coef.sd <- prior$coef.sd

    if (formula_interface) {
      prior_on_rates <- identical(prior$target, "rates")
      prior_on_coefficients <- identical(prior$target, "coefficients")

      if (prior_on_coefficients) {
        intercept_idx <- grep(":\\(Intercept\\)$", q_prep$par_names)
        slope_idx <- setdiff(seq_along(q_prep$par_names), intercept_idx)
      }
    }
  }

  ## tree topology, precomputed once: child edges of each internal node (for
  ## "expm"), edge groups by level (for "eigen"), and the fixed-node masks
  child_rows <- split(seq_len(nrow(phy$edge)), phy$edge[, 1])[as.character(anc)]
  level_groups <- prune_levels(phy)
  fixer_mat <- NULL
  if (!is.null(phy$node.label)) {
    fixer_mat <- matrix(1, nb.tip + nb.node, nq)
    for (focal in anc) {
      if (!is.na(phy$node.label[focal - nb.tip])) {
        fixer.tmp <- numeric(nq / rate.cat)
        fixer.tmp[phy$node.label[focal - nb.tip]] <- 1
        fixer_mat[focal, ] <- rep(fixer.tmp, rate.cat)
      }
    }
  }

  if (is.null(root.p)) {
    adjacency <- if (formula_interface) {
      q_prep$Q_indicator > 0
    } else {
      rate < max(rate)
    }
    diag(adjacency) <- FALSE
    active <- rowSums(adjacency) + colSums(adjacency) > 0
    root.p <- as.numeric(active) / sum(active)
  }

  model_data <- list(
    phy = phy,
    liks = liks,
    Q = Q,
    rate = rate,
    root.p = root.p,
    rate.cat = rate.cat,
    set.fog = set.fog,
    fog.vec = fog.vec,
    pen.type = pen.type,
    lambda = lambda,
    q_prep = q_prep,
    formula_interface = formula_interface,
    nb.node = nb.node,
    nb.tip = nb.tip,
    TIPS = TIPS,
    anc = anc,
    k.rates = k.rates,
    nq = nq,
    pen_indices = pen_indices,
    fog_idx = fog_idx,
    has_prior = has_prior,
    prior_on_rates = prior_on_rates,
    prior_on_coefficients = prior_on_coefficients,
    prior.mean = prior.mean,
    prior.sd = prior.sd,
    coef.sd = coef.sd,
    intercept_idx = intercept_idx,
    slope_idx = slope_idx,
    child_rows = child_rows,
    level_groups = level_groups,
    fixer_mat = fixer_mat
  )

  ## rate matrix from the model parameters (fog parameters removed); works
  ## for both AD and plain numeric parameters
  build_Q <- function(p_model) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    eta_rates <- NULL
    if (formula_interface) {
      eta_rates <- drop(q_prep$rate_X %*% p_model)
      Q <- as.matrix(q_prep$Q0)
      for (b in q_prep$blocks) {
        Q[cbind(b$from, b$to)] <-
          exp(drop(b$X %*% p_model[b$par_index]))
      }
      for (i in seq_len(nrow(Q))) {
        Q[i, i] <- -sum(Q[i, ])
      }

      p_rates <- exp(eta_rates)
      np <- length(p_rates)
    } else {
      p_rates <- exp(p_model)
      np <- length(p_rates)

      Q <- rate
      Q[Q <= np] <- p_rates[Q[Q <= np]]
      Q[rate == np + 1] <- 0

      for (i in seq_len(nrow(Q))) {
        Q[i, i] <- -sum(Q[i, ])
      }
    }
    list(Q = Q, eta_rates = eta_rates, p_rates = p_rates, np = np)
  }

  ## plain numeric Q at parameter vector p (for the degeneracy check)
  Q_at <- function(p) {
    p_model <- if (n_fog) p[-seq_len(n_fog)] else p
    build_Q(as.numeric(p_model))$Q
  }

  make_prune_nll <- function(method) function(pars) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    "diag<-" <- RTMB::ADoverload("diag<-")
    RTMB::getAll(pars, model_data, warn = FALSE)

    p_raw <- p
    comp <- numeric(nb.tip + nb.node)

    p_model <- p_raw

    if (set.fog) {
      fog_pars <- seq.int(length(unique(fog.vec)))

      ## fog parameters are still on log scale and must be exponentiated
      tip.fog.tmp <- exp(p_raw[fog_pars])

      ## remaining parameters are either old log-rates or formula coefficients
      p_model <- p_raw[-fog_pars]

      tip.fog <- numeric(length(fog.vec))
      tip.fog[] <- c(tip.fog.tmp, 0)[fog.vec]

      for (tip.index in seq_len(nb.tip)) {
        idx <- fog_idx[[tip.index]]
        if (length(idx$zeros) > 0) {
          if (rate.cat > 1) {
            liks[tip.index, idx$ones] <-
              1 - (sum(tip.fog[idx$not1]) / rate.cat)
          } else {
            liks[tip.index, idx$ones] <-
              1 - sum(tip.fog[idx$not1])
          }

          liks[tip.index, idx$zeros] <- tip.fog[idx$zeros]
        }
      }
    }

    qq <- build_Q(p_model)
    Q <- qq$Q
    eta_rates <- qq$eta_rates
    p_rates <- qq$p_rates
    np <- qq$np

    pen_score <- 0
    if (!is.null(pen.type)) {
      if (rate.cat == 1) {
        diag_vals <- numeric(nq)
        for (i in 1:nq) diag_vals[i] <- Q[i, i]

        if (pen.type == "l1") {
          pen_score <- mean(-diag_vals)
        } else if (pen.type == "l2") {
          pen_score <- mean(diag_vals^2)
        } else if (pen.type == "er") {
          if (np > 1) {
            pen_score <- mean((p_rates - mean(p_rates))^2)
          }
        }
      } else {
        pen_score_vec <- numeric(rate.cat)
        for (rc_i in seq_len(rate.cat)) {
          idx <- pen_indices[[rc_i]]
          rc_offdiag <- numeric(0)
          for (r in idx) {
            for (cc in idx) {
              if (r != cc) rc_offdiag <- c(rc_offdiag, Q[r, cc])
            }
          }
          if (pen.type == "l1") {
            pen_score_vec[rc_i] <- mean(rc_offdiag)
          } else if (pen.type == "l2") {
            pen_score_vec[rc_i] <- mean(rc_offdiag^2)
          } else if (pen.type == "er") {
            n_rc <- length(rc_offdiag)
            if (n_rc > 1) {
              pen_score_vec[rc_i] <- mean((rc_offdiag - mean(rc_offdiag))^2)
            }
          }
        }
        pen_score <- sum(pen_score_vec)
      }
    }

    if (method == "eigen") {
      ## P(t) x = V exp(Lambda t) V^{-1} x, from one eigendecomposition of Q.
      ## The eigensystem is complex in general (Q is not symmetric). RTMB's
      ## complex AD type (adcomplex, a pair of real advectors) lacks methods for
      ## matrix(), outer() etc., so carry the real and imaginary parts as plain
      ## advectors and use real AD arithmetic
      ev <- RTMB::eigen(Q, symmetric = FALSE)
      V <- ev$vectors
      Vi <- RTMB::solve(V)
      tVr <- t(Re(V)); tVim <- t(Im(V))
      tVir <- t(Re(Vi)); tViim <- t(Im(Vi))
      ## t() rather than matrix(), which drops the AD class
      lam_r <- t(Re(ev$values))
      lam_i <- t(Im(ev$values))
      ## One pass per tree level: each row of Y is P(t_e) x for one edge e,
      ## computed in the eigenbasis as P(t) x = V exp(Lambda t) V^{-1} x for
      ## all edges at the level at once. Rows are transposed (row vectors), so
      ## the code multiplies on the right by t(V^{-1}) and t(V), and the complex
      ## product is expanded into real and imaginary parts: with W = x V^{-1}
      ## and E = exp(Lambda t), only Re((W * E) V) is needed. Sibling rows are
      ## then multiplied into their parent's row (via L$slots), any fixed-node
      ## mask is applied, and the result is normalized.
      for (L in level_groups) {
        X <- liks[L$child, , drop = FALSE]
        tt <- matrix(L$t, ncol = 1)
        ea <- exp(tt %*% lam_r)
        tb <- tt %*% lam_i
        Er <- ea * cos(tb)
        Ei <- ea * sin(tb)
        Wr <- X %*% tVir
        Wi <- X %*% tViim
        ## one row per edge: (P(t) x)^T for that edge's child
        Y <- (Wr * Er - Wi * Ei) %*% tVr - (Wr * Ei + Wi * Er) %*% tVim
        s1 <- L$slots[[1]]
        v <- Y[s1$src, , drop = FALSE]
        for (s in L$slots[-1]) {
          v[s$dst, ] <- v[s$dst, , drop = FALSE] * Y[s$src, , drop = FALSE]
        }
        if (!is.null(fixer_mat)) v <- v * fixer_mat[L$nodes, , drop = FALSE]
        ## row sums as a matrix product: base rowSums() drops the AD class
        cs <- (v %*% rep(1, nq))[, 1]
        comp[L$nodes] <- cs
        liks[L$nodes, ] <- v / cs
      }
    } else {
      for (i in seq_len(nb.node)) {
        focal <- anc[i]
        desRows <- child_rows[[i]]
        desNodes <- phy$edge[desRows, 2]
        v <- 1
        for (desIndex in seq_along(desRows)) {
          v <- drop(as.matrix(v * Matrix::expm(Q * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]))
        }

        if (!is.null(fixer_mat)) v <- v * fixer_mat[focal, ]

        comp[focal] <- sum(v)
        liks[focal, ] <- v / comp[focal]
      }
    }

    root <- nb.tip + 1L

    if (is.character(root.p)) {
      if (root.p == "yang") {
        root.p <- Matrix::expm(10000 * Q)[1,]
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
      } else {
        root.p <- liks[root,] / sum(liks[root,])
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root,])))))
      }
    } else {
      if (is.numeric(root.p[1])) {
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root,])))))
      }
    }

    prior_nll <- p_model[1] * 0

    if (has_prior) {
      if (formula_interface && prior_on_rates) {
        x <- eta_rates
        prior_nll <- prior_nll + sum(
          0.5 * log(2 * pi * prior.sd^2) +
            0.5 * ((x - prior.mean) / prior.sd)^2
        )
      } else if (formula_interface && prior_on_coefficients) {
        if (length(intercept_idx) > 0) {
          x <- p_model[intercept_idx]
          prior_nll <- prior_nll + sum(
            0.5 * log(2 * pi * prior.sd^2) +
              0.5 * ((x - prior.mean) / prior.sd)^2
          )
        }

        if (length(slope_idx) > 0) {
          x <- p_model[slope_idx]
          prior_nll <- prior_nll + sum(
            0.5 * log(2 * pi * coef.sd^2) +
              0.5 * (x / coef.sd)^2
          )
        }
      } else {
        prior_nll <- prior_nll + sum(
          0.5 * log(2 * pi * prior.sd^2) +
            0.5 * ((p_model - prior.mean) / prior.sd)^2
        )
      }
    }

    likelihood_nll <- loglik + pen_score * lambda
    RTMB::REPORT(likelihood_nll)
    likelihood_nll + prior_nll
  }

  mk_obj <- function(method) {
    RTMB::MakeADFun(
      make_prune_nll(method),
      list(p = p),
      silent = TRUE
    )
  }

  if (pruning == "auto") {
    ## structural check: repeated eigenvalues at two arbitrary, irregular
    ## parameter vectors mean they are repeated everywhere (e.g. ER with > 2
    ## states). Deterministic, so the user's random-number stream is untouched
    k <- seq_along(p)
    p_test <- list(-1 - 0.3 * sqrt(k), -0.5 + 0.2 * log(k + 1.7) * (-1)^k)
    structural <- all(vapply(p_test, function(pp) eigen_degenerate(Q_at(pp)),
                             logical(1)))
    if (structural) pruning <- "expm"
  }
  if (pruning == "expm") return(mk_obj("expm"))
  if (pruning == "eigen") return(mk_obj("eigen"))

  ## "auto": the returned object is built on the per-edge expm tape, so
  ## anything that uses the tape directly (tmbstan, sdreport, $he, $env$f)
  ## gets the safe method. Only $fn, $gr and $report use the faster eigen
  ## tape, falling back to expm where Q is (nearly) degenerate or where the
  ## eigen tape fails. RTMB's eigen atomic normalizes eigenvectors by
  ## sqrt(v^T v), which is ~0 for some complex eigenvectors (e.g. when a
  ## cycle of transitions dominates Q); its reverse pass then calls solve()
  ## on a nearly singular matrix and throws an R error, which aborts a Stan
  ## chain rather than rejecting the proposal.
  obj <- mk_obj("expm")
  obj_eigen <- mk_obj("eigen")
  fn_expm <- obj$fn
  gr_expm <- obj$gr
  report_expm <- obj$report
  env <- obj$env

  ## evaluate with the eigen tape unless Q is degenerate there or the
  ## result is an error or non-finite; otherwise use the expm tape
  eval_safe <- function(f_eigen, f_expm, x, check_finite = TRUE, ...) {
    if (eigen_degenerate(Q_at(x))) return(f_expm(x, ...))
    ans <- tryCatch(f_eigen(x, ...), error = function(e) NULL)
    if (is.null(ans) || (check_finite && any(!is.finite(ans)))) {
      ans <- f_expm(x, ...)
    }
    ans
  }

  obj$fn <- function(x = obj$par, ...) {
    ans <- eval_safe(obj_eigen$fn, fn_expm, x, ...)
    ## keep the returned object's bookkeeping (as its own $fn would)
    env$last.par <- x
    if (is.finite(ans) && ans < env$value.best) {
      env$last.par.best <- x
      env$value.best <- ans
    }
    ans
  }
  obj$gr <- function(x = obj$par, ...) {
    eval_safe(obj_eigen$gr, gr_expm, x, ...)
  }
  obj$report <- function(par = env$last.par, ...) {
    eval_safe(obj_eigen$report, report_expm, par, check_finite = FALSE, ...)
  }
  obj
}