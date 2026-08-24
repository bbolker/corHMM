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


mkdev.corhmm_rtmb <- function(p, phy, liks, Q, rate, root.p, rate.cat,
                              order.test, lewis.asc.bias,
                              set.fog = FALSE, fog.vec = NULL,
                              pen.type = NULL, lambda = 1,
                              q_prep = NULL, prior = NULL) {
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
    prior.mean <- prior$mean
    prior.sd <- prior$sd
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
    slope_idx = slope_idx
  )

  prune_nll <- function(pars) {
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

    if (formula_interface) {
      eta_rates <- drop(q_prep$rate_X %*% p_model)
      p_rates <- exp(eta_rates)
      Q <- as.matrix(q_prep$Q0)
      Q[cbind(q_prep$edge_table$from, q_prep$edge_table$to)] <- p_rates
      for (i in seq_len(nrow(Q))) {
        Q[i, i] <- -sum(Q[i, ])
      }

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

    for (i in seq_len(nb.node)) {
      focal <- anc[i]
      desRows <- which(phy$edge[,1] == focal)
      desNodes <- phy$edge[desRows, 2]
      v <- 1
      for (desIndex in seq_along(desRows)) {
        v <- drop(as.matrix(v * Matrix::expm(Q * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]))
      }

      if (!is.null(phy$node.label)) {
        if (!is.na(phy$node.label[focal - nb.tip])) {
          fixer.tmp <- numeric(dim(Q)[2] / rate.cat)
          fixer.tmp[phy$node.label[focal - nb.tip]] <- 1
          fixer <- rep(fixer.tmp, rate.cat)
          v <- v * fixer
        }
      }

      comp[focal] <- sum(v)
      liks[focal, ] <- v / comp[focal]
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
    likelihood_nll + prior_nll
  }

  RTMB::MakeADFun(
    prune_nll,
    list(p = p),
    silent = TRUE
  )
}