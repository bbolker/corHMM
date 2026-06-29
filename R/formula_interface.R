is_formula_model <- function(model) {
  inherits(model, "formula") ||
    (
      is.list(model) &&
        length(model) > 0L &&
        all(vapply(model, inherits, logical(1), "formula"))
    )
}


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
        par_index <- seq(par_counter, length.out = n_col)
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
      par_index <- seq(par_counter, length.out = n_col)
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

  list(
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
    par_index = lapply(blocks, `[[`, "par_index")
  )
}


Q_prep <- function(formula_list, nstate = 2) {
  trans <- translate(formula_list = formula_list, nstate = nstate)

  list(
    d = nrow(trans$traitMatrix),
    Q0 = trans$Q0,
    Q_indicator = trans$Q_indicator,
    n_par = trans$n_par,
    par_names = trans$par_names,
    blocks = trans$blocks,
    trans = trans
  )
}


build_Q_formula <- function(q_par, q_prep) {
  "[<-" <- RTMB::ADoverload("[<-")

  Q <- as.matrix(q_prep$Q0)

  for (nm in names(q_prep$blocks)) {
    b <- q_prep$blocks[[nm]]

    beta_block <- q_par[b$par_index]
    eta_block <- drop(b$X %*% beta_block)
    rate_block <- exp(eta_block)

    Q[cbind(b$from, b$to)] <- rate_block
  }

  for (i in seq_len(nrow(Q))) {
    Q[i, i] <- -sum(Q[i, ])
  }

  Q
}


multi_to_single <- function(traits, n = NULL) {
  if (is.null(n)) {
    n <- apply(traits, 2, max) + 1
    warning("guessing number of states per trait from max() + 1")
  }

  x <- rev(cumprod(rev(n)))
  x <- c(x[-1], 1)

  rowSums(sweep(traits, MARGIN = 2, x, "*")) + 1
}


corHMM_formula <- function(phy, data, model, rate.cat, root.p,
                           lower.bound, upper.bound, call,
                           verbose = TRUE) {
  if (rate.cat != 1) {
    stop("formula models currently require `rate.cat = 1`.",
         call. = FALSE)
  }

  if (!inherits(data, "data.frame")) {
    data <- as.data.frame(data)
  }

  formula_list <- if (inherits(model, "formula")) list(model) else model

  matching <- match.tree.data(phy, data)
  phy <- ape::reorder.phylo(matching$phy, "pruningwise")
  data <- matching$data
  data <- data[match(phy$tip.label, data[[1]]), , drop = FALSE]

  q_prep <- Q_prep(
    formula_list = formula_list,
    nstate = 2L
  )

  trait_data <- data[, q_prep$trans$trait_names, drop = FALSE]
  trait_data[] <- lapply(trait_data, as.numeric)

  trait_values <- multi_to_single(
    traits = trait_data,
    n = rep(2L, length(q_prep$trans$trait_names))
  )

  nb.tip <- ape::Ntip(phy)
  nb.node <- ape::Nnode(phy)

  liks <- matrix(0, nrow = nb.tip + nb.node, ncol = q_prep$d)
  for (i in seq_len(nb.tip)) {
    liks[i, trait_values[i]] <- 1
  }

  pars <- setNames(
    rep(log(0.1), q_prep$n_par),
    q_prep$par_names
  )

  RTMB_obj <- mkdev.corhmm_rtmb(
    p = pars,
    phy = phy,
    liks = liks,
    Q = q_prep$Q0,
    rate = NULL,
    root.p = root.p,
    rate.cat = rate.cat,
    order.test = FALSE,
    lewis.asc.bias = FALSE,
    set.fog = FALSE,
    fog.vec = NULL,
    pen.type = NULL,
    lambda = 1,
    q_prep = q_prep
  )

  lower <- rep(log(lower.bound), length(pars))
  upper <- rep(log(upper.bound), length(pars))

  opt.time <- system.time(
    out <- stats::nlminb(
      start = pars,
      objective = RTMB_obj$fn,
      gradient = RTMB_obj$gr,
      lower = lower,
      upper = upper
    )
  )

  est.pars <- as.numeric(out$par)
  names(est.pars) <- names(pars)

  loglik <- -as.numeric(out$objective)

  k <- length(est.pars)
  ntip <- ape::Ntip(phy)

  AIC <- -2 * loglik + 2 * k
  AICc <- -2 * loglik + (2 * k * (ntip / (ntip - k - 1)))

  solution <- matrix(est.pars, ncol = 1)
  rownames(solution) <- names(est.pars)
  colnames(solution) <- "estimate"

  obj <- list(
    loglik = loglik,
    AIC = AIC,
    AICc = AICc,
    rate.cat = rate.cat,
    solution = solution,
    index.mat = q_prep$Q_indicator,
    data = data,
    data.legend = data,
    phy = phy,
    states = NA,
    tip.states = NA,
    states.info = NA,
    iterations = out$iterations,
    collapse = NA,
    root.p = root.p,
    args.list = list(
      phy = phy,
      data = data,
      model = model,
      q_prep = q_prep,
      trait_values = trait_values,
      pars = est.pars
    ),
    call = call,
    devfun = RTMB_obj,
    opt.time = opt.time,
    formula_model = TRUE,
    formula_list = formula_list
  )

  class(obj) <- "corhmm"
  obj
}
