library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(expm)
library(RTMB)
library(future)
library(future.apply)
library(parallel)

is_formula_model <- function(model) {
  inherits(model, "formula") ||
    (
      is.list(model) &&
        length(model) > 0L &&
        all(vapply(model, inherits, logical(1), "formula"))
    )
}

translate <- function(formula_list, nstate = 2) {
  trait_names <- vapply(formula_list, \(f) as.character(f[[2]]), character(1))
  names(formula_list) <- trait_names

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
  
  edge_groups <- split(edge_tab, list(edge_tab$changed_trait, edge_tab$direction))
  edge_groups <- edge_groups[sort(names(edge_groups))]

  for (group_name in names(edge_groups)) {
    dat <- edge_groups[[group_name]]
    if (nrow(dat) == 0) next
    
    tr <- as.character(dat$changed_trait[1])
    dir <- as.character(dat$direction[1])
    
    curr_formula <- formula_list[[tr]]
    rhs_terms <- delete.response(terms(curr_formula))
    
    X <- model.matrix(rhs_terms, data = dat)
    n_col <- ncol(X)
    
    block_name <- paste(tr, dir, sep = "_")
    par_index <- seq(par_counter, length.out = n_col)
    coef_names <- paste0(block_name, ":", colnames(X))

    blocks[[block_name]] <- list(block_name = block_name, trait = tr, direction = dir,
                                  formula = curr_formula, rhs_terms = rhs_terms, edge_id = dat$edge_id,
                                  from = dat$from, to = dat$to,
                                  from_label = dat$from_label, to_label   = dat$to_label,
                                  X = X, n_par = n_col, coef_names = coef_names, par_index = par_index
                                )
    par_counter <- par_counter + n_col
  }

  ns <- nrow(traitMatrix)
  Q_indicator <- matrix(0L, ns, ns)
  Q_indicator[as.matrix(edge_tab[, c("from", "to")])] <- edge_tab$edge_id
  rownames(Q_indicator) <- colnames(Q_indicator) <- traitMatrix$label

  ## Q0_sparse <- Matrix::Matrix(0, ns, ns, sparse = TRUE, dimnames = list(traitMatrix$label, traitMatrix$label))
  Q0 <- matrix(0, ns, ns, dimnames = list(traitMatrix$label, traitMatrix$label))
  ## TODO: uncomment this and see what breaks (what's the error?)
  ## Q0 <- matrix(0, ns, ns, dimnames = list(traitMatrix$label, traitMatrix$label))
  ## Q0 <- RTMB::AD(Q0)  ## convert base-R to RTMB/AD type

  list(trait_names = trait_names, stateList = stateList, formulas = formula_list, traitMatrix = traitMatrix,
    edge_table = edge_tab, blocks = blocks, Q_indicator = Q_indicator, Q0 = Q0, # Q0 = Q0_sparse,
    n_par = par_counter - 1L, par_names = unlist(lapply(blocks, `[[`, "coef_names"), use.names = FALSE), par_index = lapply(blocks, `[[`, "par_index")
  )
}

#https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html
cmb <- function(f, d) function(p) f(p, d)

Q_template <- function(n=2, k= 3, set_indices = TRUE) {
  if (length(n) == 1) {
    n <- rep(n, k)
  }
  all_states <- do.call(expand.grid, lapply(n, \(x) 0:(x-1)))
  ## dimnms to match corHMM standard
  dimnms <- apply(all_states, 1, \(x) paste(x, collapse = "|"))
  ns <- prod(n)
  m <- matrix(0, ns, ns)
  for (i in 1:ns) {
    ## exactly one state changes ...
    for (j in 1:ns) {
      m[i,j] <- as.numeric(sum(all_states[i,] != all_states[j,])== 1)
    }
  }
  dimnames(m) <- list(dimnms, dimnms)
  if (set_indices) {
    m[m!=0] <- seq_len(sum(m==1))
  }
  return(m)
}

Q_prep <- function(mode = c("rates", "formula"), formula_list = NULL, nstate = 2, ntrait = 3) {
  mode <- match.arg(mode)
  if (mode == "formula") {
    trans <- translate(formula_list = formula_list, nstate = nstate)

    return(list(mode = "formula", d = nrow(trans$traitMatrix),
                Q0 = trans$Q0, Q_indicator = trans$Q_indicator,
                n_par = trans$n_par, par_names = trans$par_names,
                blocks = trans$blocks, trans = trans
              ))
  }
  ## rates mode: build Q_indicator internally
  Q_indicator <- Q_template(n = nstate, k = ntrait, set_indices = TRUE)
  ns <- nrow(Q_indicator)
  Q0 <- matrix(0, ns, ns, dimnames = dimnames(Q_indicator))
  # Q0_sparse <- Matrix(0, ns, ns, sparse = TRUE, dimnames = dimnames(Q_indicator))
  list(mode = "rates", d = nrow(Q_indicator),
      Q_indicator = Q_indicator, Q0 = Q0, # Q0 = Q0_sparse,
      n_par = max(Q_indicator), par_names = paste0("log_rate_", seq_len(max(Q_indicator)))
  )
}

build_Q <- function(q_par, q_prep) {
  "[<-" <- RTMB::ADoverload("[<-")

  Q <- as.matrix(q_prep$Q0)

  if (q_prep$mode == "formula") {
    for (nm in names(q_prep$blocks)) {
      b <- q_prep$blocks[[nm]]

      beta_block <- q_par[b$par_index]
      eta_block  <- drop(b$X %*% beta_block)
      rate_block <- exp(eta_block)

      Q[cbind(b$from, b$to)] <- rate_block
    }

  } else if (q_prep$mode == "rates") {
    idx <- which(q_prep$Q_indicator != 0)
    map <- as.integer(q_prep$Q_indicator[idx])
    rates <- exp(q_par)

    Q[idx] <- rates[map]

  } else {
    stop("unknown q_prep mode", call. = FALSE)
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

gl_pairs <- function(q_prep) {
  if (q_prep$mode == "rates") {
    Qind <- q_prep$Q_indicator
    idx <- which(upper.tri(Qind) & Qind != 0, arr.ind = TRUE)

    pairs <- vector("list", 0)
    for (r in seq_len(nrow(idx))) {
      i <- idx[r, "row"]
      j <- idx[r, "col"]

      a <- Qind[i, j]
      b <- Qind[j, i]

      if (a != 0 && b != 0) {
        pairs[[length(pairs) + 1L]] <- c(max(a, b), min(a, b))
      }
    }
    return(pairs)
  } else if (q_prep$mode == "formula") {
    nm <- q_prep$par_names
    gain_idx <- grep("_gain:", nm, fixed = TRUE)

    pairs <- vector("list", 0)
    for (i in gain_idx) {
      partner <- sub("_gain:", "_loss:", nm[i], fixed = TRUE)
      j <- match(partner, nm)
      if (!is.na(j)) {
        pairs[[length(pairs) + 1L]] <- c(i, j)
      }
    }
    return(pairs)
  }
}

#' @param p parameters (log-hazard rates)
#' @param lb lower bound(s) for baseline priors
#' @param ub upper bound(s)
#' @param range width of Gaussian (+/- SD between mean and lower/upper bounds)
#' @param gainloss_pairs
#' @param lb_gainloss
#' @param ub_gainloss
#' @param range_gainloss number of SDs from center to lower/upper bounds
#' @param nllfun \emph{negative} log-likelihood function
#' @param negative return negative log posterior?
postfun <- function(pars, Phylodata,
                    ## add whatever arguments the RTMB pruning algorithm loglik function
                    ## needs (tree, trait data, etc.)
                    # p,
                    lb = log(1e-9), ub = log(1e2), range = 3,
                    # gainloss_pairs = NULL,
                    lb_gainloss = log(1e-3), ub_gainloss = log(1e3), range_gainloss = 3,
                    negative = FALSE
                    ) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, Phylodata)
  
  ## call RTMB pruning-algorithm code here to compute log-likelihood ...
  nll <- prune_nll(pars, Phylodata)
  loglik <- -1*nll

  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  p <- q_par
  logdnorm <- function(x, mu, sd) {
    -0.5 * log(2*pi*sd^2) - 0.5 * ((x - mu)/sd)^2
  }
  log.prior  <- sum(logdnorm(p, prior.mean, prior.sd))
  ## product of likelihood and prior -> sum of LL and log-prior
  res <- loglik + log.prior
  ## calculate gain/loss priors

  gl.prior.mean <- (lb_gainloss + ub_gainloss) / 2
  gl.prior.sd <- (ub_gainloss - lb_gainloss) / (2 * range_gainloss)
  ## vapply might not work in RTMB? replace with for loop?
  gl.values <- p[1] * 0 + numeric(length(gainloss_pairs))
  for (i in seq_len(length(gainloss_pairs))) {
    idx <- gainloss_pairs[[i]]
    gl.values[i] <- p[idx[1]] - p[idx[2]]
  }
  gl.log.prior <- sum(logdnorm(gl.values, gl.prior.mean, gl.prior.sd))

  res <- -1*(res + gl.log.prior)
  
  return(res)
}

pack_refit <- function(start, res) {
  start <- as.numeric(start)
  fitted <- as.numeric(res$par)

  names(start)  <- paste0("start",  seq_along(start))
  names(fitted) <- paste0("fitted", seq_along(fitted))

  data.frame(objective = res$objective, convergence = res$convergence, message = res$message,
             rbind(start),rbind(fitted), check.names = FALSE)
}

random_refit <- function(task, Phylodata, pars0, opt.args = NULL) {

  t0 <- proc.time()[["elapsed"]]

  ff <- RTMB::MakeADFun(
    func = cmb(postfun, Phylodata),
    parameters = list(q_par = pars0),
    silent = TRUE
  )

  t1 <- proc.time()[["elapsed"]]

  res <- with(ff, do.call(nlminb, c(list(task$start, fn, gr), opt.args)))

  t2 <- proc.time()[["elapsed"]]

  df <- pack_refit(task$start, res)

  df$refit_id <- task$id
  df$start_seed  <- task$seed
  df$pid <- Sys.getpid()
  df$nodename <- Sys.info()[["nodename"]]
  df$time_MakeADFun <- t1 - t0          # MakeADFun
  df$time_opt <- t2 - t1          # nlminb 
  df$time_total <- t2 - t0   

  df
}

#' will simulate random trait data if `traitMatrix` is missing
#' @param verbose print output?
#' @examples
#' set.seed(427)
#' g1 <- reorder(ape::rtree(20), "pruningwise")
#' log_trans_rates <- log(abs(rnorm(8)))
#' res <- postAD(tree = g1, trait = 2, state = 2, pars = log_trans_rates, multistart = 10, seed = 427)
#' objfun <- postAD(g1, 2, 2, log_trans_rates, multistart = 10, seed = 427, return_obj = TRUE)
#' objfun$fn(objfun$par)
## FIXME: trait, state should probably be 'ntrait', 'nstate'
postAD <- function(tree, state = 2, pars = NULL, traitMatrix = NULL, formula_list = NULL,
                    multistart = 1, parallel = FALSE, jitter.sd = 0.5, seed = 427, 
                    rng_misuse = c("warning","error","ignore"), opt.args = NULL, keep_all = FALSE,
                    return_obj = FALSE, return_obj_type = c("AD", "raw"), verbose = FALSE
                  ) {

  return_obj_type <- match.arg(return_obj_type)
  mode <- if (is.null(formula_list)) 'rates' else  'formula'
  rng_misuse <- match.arg(rng_misuse)
  t_total0 <- proc.time()[["elapsed"]]
  
  n_non_numeric <- sum(!vapply(traitMatrix, is.numeric, logical(1)))
  if (n_non_numeric != 1L) stop("traitMatrix must look like: one non-numeric column for Species, and all remaining columns must be numeric trait columns.")
  ntrait <- ncol(traitMatrix) - 1L
  q_prep <- Q_prep(mode = mode, formula_list = formula_list, nstate = state, ntrait = ntrait)
  set.seed(seed)
  if (is.null(pars)) pars <- setNames(log(abs(rnorm(q_prep$n_par))), q_prep$par_names)
  
  set.seed(seed)
  tasks <- lapply(seq_len(multistart - 1), function(i) {
    si <- seed + i
    set.seed(si)
    list(id = i, seed = si, start = as.numeric(pars) + rnorm(length(pars), 0, jitter.sd))
  })
  
  if (is.null(traitMatrix)) {
    repeat {
      set.seed(seed)
      s <- phangorn::simSeq(tree, l = 1, Q = q_prep$Q_indicator,
                            type = "USER", levels = seq(state^ntrait), rate = 1)
      if (nrow(unique(as.character(s))) == prod(state^ntrait)) break
    }
    s <- as.numeric(unlist(s))
  } else {
    nvec <- if (length(state) == 1) rep(state, ntrait) else state
    spnames <- NULL

    if (is.data.frame(traitMatrix) && ncol(traitMatrix) > 0) {
        first_col <- traitMatrix[[1]]

        if (is.character(first_col) || is.factor(first_col)) {
            spnames <- as.character(first_col)
            traitMatrix <- traitMatrix[, -1, drop = FALSE]
        }
    }
    traitMatrix <- as.data.frame(traitMatrix)
    traitMatrix[] <- lapply(traitMatrix, as.numeric)

    s <- multi_to_single(traitMatrix, nvec)
    if (!is.null(spnames)) {
      stopifnot(length(spnames) == ape::Ntip(tree))
      stopifnot(all(spnames == tree$tip.label))
      names(s) <- spnames
    }
  }
  gainloss_pairs <- gl_pairs(q_prep) 
  Phylodata <- list(q_prep = q_prep, tree = tree, trait_values = s, gainloss_pairs = gainloss_pairs)

  if (return_obj  && return_obj_type == "raw") return(list(fn = postfun, pars = list(q_par = pars), nll = function(par = pars) prune_nll(list(q_par = par), Phylodata), data = Phylodata))
    
  # baseline
  t5 <- proc.time()[["elapsed"]]
  ff0 <- RTMB::MakeADFun(func = cmb(postfun, Phylodata),
                          parameters = list(q_par = pars),
                          silent = TRUE
                        )
  if (return_obj) return(ff0)
  
  t6 <- proc.time()[["elapsed"]]
  res0 <- with(ff0, do.call(nlminb, c(list(ff0$par, ff0$fn, ff0$gr), opt.args)))
  t7 <- proc.time()[["elapsed"]]
  df0  <- pack_refit(ff0$par, res0)
  df0$refit_id <- 0L
  df0$start_seed <- seed
  df0$pid <- Sys.getpid()
  df0$nodename <- Sys.info()[["nodename"]]
  df0$time_MakeADFun <- t6 - t5
  df0$time_opt <- t7 - t6 #nlminb
  df0$time_total <- t7 - t5
  
  workers <- 1
  if (parallel && multistart > 1) {
    workers <- min(multistart - 1, parallel::detectCores() - 1)
    workers <- max(workers, 1)

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    if (workers > 1) {
      future::plan(future::cluster, workers = workers)
    } else {
      future::plan(future::sequential)
    }

    old_opt <- getOption("future.rng.onMisuse")
    on.exit(options(future.rng.onMisuse = old_opt), add = TRUE)
    options(future.rng.onMisuse = rng_misuse)
  }
    
  out_list <- if (multistart > 1) {
    future.apply::future_lapply(tasks, random_refit, Phylodata = Phylodata, pars0 = pars,
                                opt.args = opt.args, future.seed = TRUE
                                )
  } else list()
  
  result_frame <- do.call(rbind, c(list(df0), out_list))

  result_good <- subset(result_frame, convergence == 0)
  if (nrow(result_good) == 0) stop("No successful refits (convergence==0).")

  fit_cols <- grep("^fitted", names(result_good), value = TRUE)[seq_len(length(ff0$par))]
  best_row <- which.min(result_good$objective)
  pars.best <- as.numeric(result_good[best_row, fit_cols])
  obj.best <- ff0$fn(pars.best)
  gr.best  <- ff0$gr(pars.best)
  fn.cdf <- ecdf(result_good$objective)

  t_total1 <- proc.time()[["elapsed"]]
  time_total = t_total1 - t_total0

  out <- list(pars.best  = pars.best, obj.best = obj.best, gr.best = gr.best, Phylo = Phylodata,
              multistart = multistart, workers = workers, cores_detected = parallel::detectCores(),
              time_total = time_total
              )
  if (keep_all) out$result_frame <- result_frame
  if (verbose) {
     cat('loglik:\n'); print(obj.best)
     cat('pars:\n'); print(pars.best)
     cat('time:\n'); print(time_total)
     cat('loglik cdf:\n'); plot(fn.cdf)
  }
  out
}

corHMM_formula <- function(phy, data, model, rate.cat, root.p,
                           lower.bound, upper.bound, call,
                           verbose = TRUE) {
  if (rate.cat != 1) {
    stop("formula models currently require `rate.cat = 1`.", call. = FALSE)
  }

  if (!inherits(data, "data.frame")) {
    data <- as.data.frame(data)
  }

  formula_list <- if (inherits(model, "formula")) list(model) else model

  matching <- match.tree.data(phy, data)
  phy <- ape::reorder.phylo(matching$phy, "pruningwise")
  data <- matching$data
  data <- data[match(phy$tip.label, data[[1]]), , drop = FALSE]

  fit <- postAD(
    tree = phy,
    traitMatrix = data,
    state = 2,
    formula_list = formula_list,
    multistart = 1,
    parallel = FALSE,
    opt.args = list(
      lower = rep(log(lower.bound), 1000),
      upper = rep(log(upper.bound), 1000)
    ),
    verbose = FALSE
  )

  pars <- fit$pars.best
  names(pars) <- fit$Phylo$q_prep$par_names

  loglik <- -as.numeric(fit$obj.best)
  k <- length(pars)
  ntip <- ape::Ntip(phy)

  AIC <- -2 * loglik + 2 * k
  AICc <- -2 * loglik + (2 * k * (ntip / (ntip - k - 1)))

  solution <- matrix(pars, ncol = 1)
  rownames(solution) <- names(pars)
  colnames(solution) <- "estimate"

  obj <- list(
    loglik = loglik,
    AIC = AIC,
    AICc = AICc,
    rate.cat = rate.cat,
    solution = solution,
    index.mat = fit$Phylo$q_prep$Q_indicator,
    data = data,
    data.legend = data,
    phy = phy,
    states = NA,
    tip.states = NA,
    states.info = NA,
    iterations = NA,
    collapse = NA,
    root.p = root.p,
    args.list = list(
      phy = phy,
      data = data,
      model = model,
      pars = pars
    ),
    call = call,
    devfun = NA,
    opt.time = fit$time_total,
    formula_model = TRUE,
    formula_list = formula_list
  )

  class(obj) <- "corhmm"
  obj
}