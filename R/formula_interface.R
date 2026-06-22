is_formula_model <- function(model) {
  inherits(model, "formula") ||
    (
      is.list(model) &&
        length(model) > 0L &&
        all(vapply(model, inherits, logical(1), "formula"))
    )
}

.is_call <- function(x, name) {
  is.call(x) && identical(as.character(x[[1]]), name)
}

.deparse_rhs <- function(x) {
  paste(deparse(x), collapse = "")
}

.formula_from_rhs <- function(lhs, rhs) {
  stats::as.formula(
    paste(lhs, "~", .deparse_rhs(rhs)),
    env = parent.frame()
  )
}

.parse_formula_spec <- function(f) {
  lhs <- as.character(f[[2]])
  rhs <- f[[3]]

  if (.is_call(rhs, "symm")) {
    rhs2 <- rhs[[2]]

    return(list(
      trait = lhs,
      symmetric = TRUE,
      gain_formula = .formula_from_rhs(lhs, rhs2),
      loss_formula = .formula_from_rhs(lhs, rhs2)
    ))
  }

  if (.is_call(rhs, "list")) {
    stop(
      "Different gain/loss formulas are not currently supported.",
      call. = FALSE
    )
  }

  list(
    trait = lhs,
    symmetric = FALSE,
    gain_formula = f,
    loss_formula = f
  )
}

translate <- function(formula_list, nstate = 2) {
  specs <- lapply(formula_list, .parse_formula_spec)
  trait_names <- vapply(specs, `[[`, character(1), "trait")

  names(formula_list) <- trait_names
  names(specs) <- trait_names

  stateList <- rep.int(as.integer(nstate), length(trait_names))
  names(stateList) <- trait_names

  traitMatrix <- do.call(expand.grid, lapply(stateList, function(x) 0:(x - 1)))
  traitMatrix$label <- do.call(paste, c(traitMatrix, sep = "|"))

  edge_table <- function(traitM, trait_names) {
    ns <- nrow(traitM)

    Xstate <- as.matrix(traitM[, trait_names, drop = FALSE])

    idx_grid <- expand.grid(
      to = seq_len(ns),
      from = seq_len(ns)
    )
    idx_grid <- idx_grid[idx_grid$from != idx_grid$to, , drop = FALSE]

    diff_matrix <- Xstate[idx_grid$to, , drop = FALSE] -
      Xstate[idx_grid$from, , drop = FALSE]

    diff_count <- rowSums(diff_matrix != 0)

    valid_idx <- which(diff_count == 1L)
    edges <- idx_grid[valid_idx, , drop = FALSE]
    final_diffs <- diff_matrix[valid_idx, , drop = FALSE]

    change_col_idx <- max.col(final_diffs != 0)

    res <- data.frame(
      edge_id = seq_along(valid_idx),
      from = edges$from,
      to = edges$to,
      from_label = traitM$label[edges$from],
      to_label = traitM$label[edges$to],
      changed_trait = trait_names[change_col_idx],
      delta = final_diffs[cbind(seq_along(change_col_idx), change_col_idx)]
    )

    res$direction <- ifelse(res$delta > 0, "gain", "loss")
    res$focal_from <- Xstate[cbind(edges$from, change_col_idx)]
    res$focal_to <- Xstate[cbind(edges$to, change_col_idx)]

    res <- cbind(
      res,
      traitM[edges$from, trait_names, drop = FALSE]
    )

    res
  }
  edge_tab <- edge_table(traitMatrix, trait_names)
  edge_tab <- edge_tab[order(edge_tab$to, edge_tab$from), ]
  edge_tab$edge_id <- seq_len(nrow(edge_tab))

  blocks <- list()
  par_counter <- 1L
  shared_symm <- list()

  edge_groups <- split(
    edge_tab,
    list(edge_tab$changed_trait, edge_tab$direction)
  )
  edge_groups <- edge_groups[sort(names(edge_groups))]

  for (group_name in names(edge_groups)) {
    dat <- edge_groups[[group_name]]
    if (nrow(dat) == 0) next

    tr <- as.character(dat$changed_trait[1])
    dir <- as.character(dat$direction[1])
    spec <- specs[[tr]]

    block_name <- paste(tr, dir, sep = "_")

    curr_formula <- if (dir == "gain") {
      spec$gain_formula
    } else {
      spec$loss_formula
    }

    rhs_terms <- delete.response(stats::terms(curr_formula))
    X <- stats::model.matrix(rhs_terms, data = dat)
    n_col <- ncol(X)

    if (isTRUE(spec$symmetric)) {
      shared_name <- paste0(tr, "_symm")

      if (is.null(shared_symm[[shared_name]])) {
        par_index <- seq(par_counter, length.out = n_col)
        coef_names <- paste0(shared_name, ":", colnames(X))

        shared_symm[[shared_name]] <- list(
          par_index = par_index,
          coef_names = coef_names,
          colnames = colnames(X)
        )

        par_counter <- par_counter + n_col
      } else {
        if (!identical(shared_symm[[shared_name]]$colnames, colnames(X))) {
          stop(
            "gain and loss design matrices differ for symmetric formula.",
            call. = FALSE
          )
        }

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

  ns <- nrow(traitMatrix)
  Q_indicator <- matrix(0L, ns, ns)
  Q_indicator[as.matrix(edge_tab[, c("from", "to")])] <- edge_tab$edge_id
  rownames(Q_indicator) <- colnames(Q_indicator) <- traitMatrix$label

  Q0 <- Matrix::Matrix(0, ns, ns, sparse = TRUE, dimnames = list(traitMatrix$label, traitMatrix$label))
  ## TODO: uncomment this and see what breaks (what's the error?)
  ## Q0 <- matrix(0, ns, ns, dimnames = list(traitMatrix$label, traitMatrix$label))
  ## Q0 <- RTMB::AD(Q0)  ## convert base-R to RTMB/AD type

  n_par <- par_counter - 1L
  par_names <- character(n_par)

  for (nm in names(blocks)) {
    b <- blocks[[nm]]

    old_names <- par_names[b$par_index]
    new_names <- b$coef_names

    if (any(nzchar(old_names) & old_names != new_names)) {
      stop(
        "internal error: conflicting parameter names for shared parameter indices.",
        call. = FALSE
      )
    }

    par_names[b$par_index] <- new_names
  }

  if (any(!nzchar(par_names))) {
    stop(
      "internal error: some formula-model parameters were not named.",
      call. = FALSE
    )
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


Q_prep <- function(mode = "formula", formula_list = NULL,
                   nstate = 2, ntrait = NULL) {
  if (!identical(mode, "formula")) {
    stop("formula interface only supports `mode = 'formula'.",
         call. = FALSE)
  }

  trans <- translate(formula_list = formula_list, nstate = nstate)

  list(
    mode = "formula",
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

prune_nll <- function(pars, Phylodata) {
  if (!require("RTMB")) stop("install RTMB package")
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, Phylodata)

  ntips <- ape::Ntip(tree)
  if (length(trait_values) != ntips) {
    stop("length(trait_values) must equal the number of tips in tree.")
  }

  d <- q_prep$d
  Q <- build_Q(q_par, q_prep)

  liks <- matrix(NA_real_, nrow = ntips + tree$Nnode, ncol = d)

  ## initialize tip likelihoods
  liks[seq_len(ntips), ] <- 0
  for (i in seq_len(ntips)) {
    liks[i, trait_values[i]] <- 1
  }

  comp <- numeric(nrow(liks))
  anc <- unique(tree$edge[, 1])

  for (i in anc) {
    desRows <- which(tree$edge[, 1] == i)
    desNodes <- tree$edge[desRows, 2]

    v <- rep(1, d)

    for (j in seq_along(desRows)) {
      t <- tree$edge.length[desRows[j]]

      u <- drop(as.matrix(
        Matrix::expm(Q * t) %*% liks[desNodes[j], ]
      ))

      v <- drop(v * u)
    }
    comp[i] <- sum(v)
    liks[i, ] <- v / comp[i]
  }

  TIPS <- 1:ape::Ntip(tree)
  root <- ape::Ntip(tree) + 1
  root.p <- rep(1/d, d)
  neg_loglik <- -1*(sum(log(comp[-TIPS])) +log(sum(root.p * liks[root, ])))
  return(neg_loglik)
}



#' translate trait matrix (no species name!) to single trait
#' @param traits trait matrix (trait values, 0-indexed)
#' @param n number of states per trait
multi_to_single <- function(traits, n = NULL) {
  if (is.null(n)) {
      n <- apply(traits, 2, max)+1
      warning("guessing number of traits per state from max()+1")
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

  ntrait <- length(formula_list)
  nstate <- 2L

  q_prep <- Q_prep(
    mode = "formula",
    formula_list = formula_list,
    nstate = nstate,
    ntrait = ntrait
  )

  trait_data <- data[, -1, drop = FALSE]
  trait_data[] <- lapply(trait_data, as.numeric)

  trait_values <- multi_to_single(
    traits = trait_data,
    n = rep(nstate, ntrait)
  )

  d <- q_prep$d
  nb.tip <- ape::Ntip(phy)
  nb.node <- ape::Nnode(phy)

  liks <- matrix(0, nrow = nb.tip + nb.node, ncol = d)
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