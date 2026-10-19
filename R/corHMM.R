## used in devfun (taken from environment)
utils::globalVariables(c("liks", "Q"))


######################################################################################################################################
######################################################################################################################################
### corHMM -- Generalized hidden Markov Models
######################################################################################################################################
######################################################################################################################################

corHMM <- function(phy, data, rate.cat, rate.mat=NULL, model = "ARD",
                   node.states = "marginal", fixed.nodes=FALSE, p=NULL,
                   root.p="yang", tip.fog=NULL, ip=NULL, fog.ip = 0.01,
                   nstarts=0, n.cores=1, get.tip.states = FALSE,
                   lewis.asc.bias = FALSE, collapse=TRUE,
                   lower.bound = 1e-9, upper.bound = 100, opts=NULL,
                   return.devfun = FALSE, use_RTMB = FALSE,
                   prior = NULL, mcmc = NULL, RTMB_control = list(),
                   verbose=TRUE) {

    call <- match.call()
    RTMB_control <- corHMM_RTMB_control(RTMB_control)

    input_rate_mat_states <- NULL
    if(!is.null(rate.mat)) input_rate_mat_states <- colnames(rate.mat)

    if(is.null(node.states)){
        obj <- NULL
        obj$loglik <- NULL
        obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to corHMM command for parameter 'node.states': joint, marginal, scaled, or none.")
        return(obj)
    } else {
        valid.models <- c("joint", "marginal", "scaled", "none")
        if(!any(valid.models == node.states)){
            obj <- NULL
            obj$loglik <- NULL
            obj$diagnostic <- paste("'", node.states, "' is not valid for ancestral state reconstruction method.  Please pass one of the following to corHMM command for parameter 'node.states': joint, marginal, scaled, or none.", sep="")
            return(obj)
        }
        if(length(node.states) > 1){
            node.states <- "marginal"
            message("No model selected for 'node.states'. Will perform marginal ancestral state estimation.\n")
        }
    }

    if(fixed.nodes == FALSE && !is.null(phy$node.label)){
        phy$node.label <- NULL
        message("You specified 'fixed.nodes=FALSE' but included a phy object with node labels. These node labels have been removed.\n")
    }

    if(!inherits(data, "data.frame")){
        data <- as.data.frame(data)
        message("Input data is required to be of class 'data.frame' - converting now.")
    }

    if(!is.null(root.p) && !is.character(root.p)) root.p <- root.p/sum(root.p)

    input.data <- data
    model_is_formula <- inherits(model, "formula") ||
        (is.list(model) && length(model) > 0 &&
         all(vapply(model, inherits, logical(1), "formula")))

    if(!is.character(model) && !model_is_formula){
        stop("`model` must be either a character string or a formula/list of formulas.", call. = FALSE)
    }

    if(!is.null(mcmc) && !use_RTMB && !model_is_formula){
        stop("HMC requires `use_RTMB = TRUE`.", call. = FALSE)
    }
    if(!is.null(prior) && !use_RTMB && !model_is_formula){
        stop("Priors currently require `use_RTMB = TRUE`.", call. = FALSE)
    }
    if(!is.null(mcmc) && is.null(prior)){
        stop("HMC requires `prior`.", call. = FALSE)
    }

    if(!is.null(prior)){
        if(is.null(prior$lb)) prior$lb <- log(lower.bound)
        if(is.null(prior$ub)) prior$ub <- log(upper.bound)
        if(is.null(prior$range)) prior$range <- 3
        if(is.null(prior$coef.sd)) prior$coef.sd <- 2
        if(is.null(prior$target)) prior$target <- "rates"
        prior$target <- match.arg(prior$target, c("rates", "coefficients"))

        if(!model_is_formula && prior$target != "rates"){
            stop("`prior$target = 'coefficients'` is only available for formula models.", call. = FALSE)
        }
    }

    lb <- log(lower.bound)
    ub <- log(upper.bound)
    if(upper.bound < lower.bound) message("Your upper bound is smaller than your lower bound.\n")

    q_prep <- NULL
    formula_list <- NULL
    trait_values <- NULL
    fog.est <- NULL
    set.fog <- FALSE

    if(model_is_formula){
        if(!is.null(rate.mat)) stop("`rate.mat` cannot be used with formula models.", call. = FALSE)
        if(rate.cat != 1) stop("formula models currently require `rate.cat = 1`.", call. = FALSE)
        if(node.states != "none") stop("formula models currently require `node.states = 'none'.", call. = FALSE)
        if(!is.null(tip.fog)) stop("formula models do not currently support `tip.fog`.", call. = FALSE)
        if(lewis.asc.bias) stop("formula models do not currently support `lewis.asc.bias = TRUE`.", call. = FALSE)

        use_RTMB <- TRUE
        order.test <- FALSE
        formula_list <- if(inherits(model, "formula")) list(model) else model

        matching <- match.tree.data(phy, input.data)
        phy <- ape::reorder.phylo(matching$phy, "pruningwise")
        data <- matching$data
        data <- data[match(phy$tip.label, data[[1]]), , drop=FALSE]

        if(any(phy$edge.length <= .Machine$double.eps)){
            warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)), immediate. = TRUE)
            phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps)
        }

        q_prep <- translate(formula_list, nstate=2L)
        trait_data <- data[, q_prep$trait_names, drop=FALSE]
        ## missing values (NA or "?") are unknown: the tip can be in any
        ## joint state compatible with its observed traits
        trait_data[] <- lapply(trait_data, function(x) {
            x <- as.character(x)
            x[x %in% "?"] <- NA
            suppressWarnings(as.integer(x))
        })

        if(!all(unlist(trait_data, use.names=FALSE) %in% c(0:1, NA))){
            stop("formula traits must be binary variables coded 0/1 (NA or \"?\" for missing).",
                 call. = FALSE)
        }
        if(any(rowSums(!is.na(trait_data)) == 0)){
            warning("some tips have no observed traits.", call. = FALSE)
        }

        complete <- stats::complete.cases(trait_data)
        trait_values <- rep(NA_integer_, nrow(trait_data))
        trait_values[complete] <- unname(match(
            do.call(paste, c(trait_data[complete, , drop=FALSE], sep="|")),
            q_prep$traitMatrix$label
        ))

        nb.tip <- ape::Ntip(phy)
        nb.node <- ape::Nnode(phy)
        liks <- matrix(0, nb.tip + nb.node, nrow(q_prep$traitMatrix))
        state_traits <- as.matrix(q_prep$traitMatrix[, q_prep$trait_names, drop=FALSE])
        for(i in seq_len(nb.tip)){
            obs <- unlist(trait_data[i, ])
            known <- !is.na(obs)
            compatible <- rowSums(state_traits[, known, drop=FALSE] !=
                                  rep(obs[known], each=nrow(state_traits))) == 0
            liks[i, compatible] <- 1
        }

        model.set.final <- list(
            np=q_prep$n_par,
            liks=liks,
            Q=q_prep$Q0,
            rate=NULL,
            index.matrix=q_prep$Q_indicator,
            fog.vec=NULL
        )

        StateNames <- q_prep$traitMatrix$label
        ## (counts of tips with complete data only)
        counts <- table(factor(trait_values, levels=seq_along(StateNames)))
        data.legend <- data
        lower <- rep(lb, model.set.final$np)
        upper <- rep(ub, model.set.final$np)

        if (verbose) {
            message("State distribution in data:\n")
            message("States:", StateNames, "\n", sep="\t")
            message("Counts:", as.integer(counts), "\n", sep="\t")
        }
    } else {
        nCol <- dim(data)[2]

        CorData <- corProcessData(data, collapse=collapse)
        data.legend <- data <- CorData$corData

        if(length(grep("&", CorData$corData[,2])) > 0){
            non_and_chars <- as.numeric(CorData$corData[,2][-grep("&", CorData$corData[,2])])
            and_chars <- as.numeric(unlist(strsplit(CorData$corData[,2][grep("&", CorData$corData[,2])], "&")))
            nObs <- max(c(non_and_chars, and_chars))
        } else {
            nObs <- max(as.numeric(CorData$corData[,2]))
        }

        matching <- match.tree.data(phy, data)
        data <- matching$data
        phy <- matching$phy

        if(nlevels(as.factor(data[,1])) <= 1 & !is.null(p)){
            obj <- NULL
            obj$loglik <- NULL
            obj$diagnostic <- paste("Character is invariant. Analysis stopped.", sep="")
            return(obj)
        } else {
            lvls <- as.factor(data[,1])
            if(nlevels(as.factor(data[,1])) == 2 && length(which(lvls == "?"))){
                obj <- NULL
                obj$loglik <- NULL
                obj$diagnostic <- paste("Character is invariant. Analysis stopped.", sep="")
                return(obj)
            }
        }

        if(any(phy$edge.length <= .Machine$double.eps)){
            warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)), immediate. = TRUE)
            phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps)
        }

        data.sort <- data.frame(data[,2], data[,2], row.names=data[,1])
        data.sort <- data.sort[phy$tip.label,]

        counts <- table(data.sort[,1])
        levels <- levels(as.factor(data.sort[,1]))
        cols <- as.factor(data.sort[,1])
        k <- 2
        order.test <- TRUE

        nb.tip <- length(phy$tip.label)
        nb.node <- phy$Nnode

        model.set.final <- rate.cat.set.corHMM.JDB(
            phy=phy, data=input.data, rate.cat=rate.cat,
            ntraits=nObs, model=model, rate.mat=rate.mat,
            collapse=collapse, CorData=CorData, matching=matching
        )
        phy <- reorder(phy, "pruningwise")

        if(!is.null(rate.mat)){
            order.test <- FALSE
            rate.mat[rate.mat == 0] <- NA
            rate <- rate.mat
            model.set.final$np <- max(rate, na.rm=TRUE)
            rate[is.na(rate)] <- max(rate, na.rm=TRUE) + 1
            model.set.final$rate <- rate
            model.set.final$index.matrix <- rate.mat
            model.set.final$Q <- matrix(0, dim(rate.mat)[1], dim(rate.mat)[2])
            col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
            row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
            drop.states <- col.sums[which(col.sums == row.sums)]
            if(length(drop.states > 0)) model.set.final$liks[,drop.states] <- 0
        }

        if(collapse){
            StateNames <- gsub("_", "|", CorData$ObservedTraits)
        } else {
            StateNames <- gsub("_", "|", CorData$PossibleTraits)
        }

        print_counts <- rep(0, length(StateNames))
        if(length(grep("&", names(counts))) > 0) counts <- counts[-grep("&", names(counts))]

        if(verbose == TRUE){
            print_counts[as.numeric(names(counts))] <- counts
            message("State distribution in data:\n")
            message("States:", StateNames, "\n", sep="\t")
            message("Counts:", print_counts, "\n", sep="\t")
        }

        lower <- rep(lb, model.set.final$np)
        upper <- rep(ub, model.set.final$np)

        if(!is.null(tip.fog)){
            if(sum(tip.fog) < 1){
                if(length(tip.fog) == 1) tip.fog <- rep(tip.fog, length(StateNames))
                if(rate.cat > 1){
                    tip.fog <- rep(tip.fog, rate.cat)
                    for(tip.index in 1:Ntip(phy)){
                        num.zeros <- length(model.set.final$liks[tip.index, which(model.set.final$liks[tip.index,] == 0)])
                        if(num.zeros > 0){
                            model.set.final$liks[tip.index, which(model.set.final$liks[tip.index,] == 1)] <- 1 - (sum(tip.fog[which(model.set.final$liks[tip.index,] != 1)])/rate.cat)
                            model.set.final$liks[tip.index, which(model.set.final$liks[tip.index,] == 0)] <- tip.fog[which(model.set.final$liks[tip.index,] == 0)]
                            set.fog <- FALSE
                        }
                    }
                } else {
                    for(tip.index in 1:Ntip(phy)){
                        num.zeros <- length(model.set.final$liks[tip.index, which(model.set.final$liks[tip.index,] == 0)])
                        if(num.zeros > 0){
                            model.set.final$liks[tip.index, which(model.set.final$liks[tip.index,] == 1)] <- 1 - sum(tip.fog[which(model.set.final$liks[tip.index,] != 1)])
                            model.set.final$liks[tip.index, which(model.set.final$liks[tip.index,] == 0)] <- tip.fog[which(model.set.final$liks[tip.index,] == 0)]
                            set.fog <- FALSE
                        }
                    }
                }
                fog.est <- tip.fog
            } else {
                if(rate.cat > 1) tip.fog <- rep(tip.fog, rate.cat)
                model.set.final$fog.vec <- tip.fog
                set.fog <- TRUE
            }
        }
    }

    n.fog <- if(set.fog) length(unique(model.set.final$fog.vec)) else 0
    np_tot <- n.fog + model.set.final$np

    if(set.fog){
        lower <- c(rep(lb, n.fog), lower)
        upper <- c(rep(log(0.50), n.fog), upper)
    }

    if(use_RTMB){
        RTMB_obj <- mkdev.corhmm_rtmb(
            rep(0, np_tot), phy,
            liks=model.set.final$liks,
            Q=model.set.final$Q,
            rate=model.set.final$rate,
            root.p=root.p,
            rate.cat=rate.cat,
            order.test=FALSE,
            lewis.asc.bias=FALSE,
            set.fog=set.fog,
            fog.vec=model.set.final$fog.vec,
            q_prep=q_prep,
            prior=prior,
            pruning=RTMB_control$pruning
        )
    }

    if(is.null(opts)){
        opts <- list(
            eval.max=1000000,
            iter.max=1000000,
            rel.tol=sqrt(.Machine$double.eps)
        )
    }

    if(!is.null(p)){
        if(verbose == TRUE) message("Calculating likelihood from a set of fixed parameters", "\n")
        start <- if(model_is_formula) as.numeric(p) else log(as.numeric(p))
        if(set.fog) start <- c(rep(log(fog.ip), n.fog), start)
        if(use_RTMB) names(start) <- names(RTMB_obj$par)

        opt.time <- system.time({
            objective <- if(use_RTMB){
                RTMB_obj$fn(start)
            } else {
                dev.corhmm(
                    start, phy, model.set.final$liks, model.set.final$Q,
                    model.set.final$rate, root.p, rate.cat, order.test,
                    lewis.asc.bias, set.fog, model.set.final$fog.vec
                )
            }
        })

        out <- list(
            par=start,
            objective=objective,
            iterations=0,
            convergence=0,
            message="fixed parameters"
        )
    } else {
        nfit <- if(is.null(ip)) max(1, nstarts) else 1

        if(!model_is_formula && is.null(ip)){
            taxa.missing.data.drop <- which(is.na(data.sort[,1]))
            if(length(taxa.missing.data.drop) != 0){
                dat <- as.matrix(data.sort)
                dat.red <- dat[-taxa.missing.data.drop,]
                phy.red <- drop.tip(phy, taxa.missing.data.drop)
                dat.red <- phyDat(dat.red, type="USER", levels=levels)
                phy.tmp <- multi2di(phy.red)
                par.score <- parsimony(phy.tmp, dat.red, method="fitch")/2
            } else {
                dat <- phyDat(as.matrix(data.sort), type="USER", levels=levels)
                phy.tmp <- multi2di(phy)
                par.score <- parsimony(phy.tmp, dat, method="fitch")/2
            }
            mean.change <- par.score/sum(phy$edge.length)
        }

        fit_one <- function(i){
            if(model_is_formula){
                if(!is.null(ip)){
                    start <- rep(ip, length.out=model.set.final$np)
                } else if(i == 1){
                    start <- rep(log(0.1), model.set.final$np)
                } else {
                    start <- rnorm(model.set.final$np, log(0.1), 1)
                }
            } else {
                if(!is.null(ip)){
                    start <- rep(ip, length.out=model.set.final$np)
                } else if(mean.change == 0){
                    start <- rep(0.01 + exp(lb), model.set.final$np)
                } else {
                    start <- sort(
                        rexp(model.set.final$np, 1/mean.change),
                        decreasing=TRUE
                    )
                }
            }

            if(!model_is_formula){
                start[start < exp(lb)] <- exp(lb)
                start[start > exp(ub)] <- exp(lb)
                if(set.fog) start <- c(rep(fog.ip, n.fog), start)
                start <- log(start)
            }

            if(use_RTMB) names(start) <- names(RTMB_obj$par)

            objective <- if(use_RTMB){
                RTMB_obj$fn
            } else {
                function(p) dev.corhmm(
                    p, phy, model.set.final$liks, model.set.final$Q,
                    model.set.final$rate, root.p, rate.cat, order.test,
                    lewis.asc.bias, set.fog, model.set.final$fog.vec
                )
            }

            stats::nlminb(
                start=start,
                objective=objective,
                gradient=if(use_RTMB) RTMB_obj$gr else NULL,
                lower=lower,
                upper=upper,
                control=opts
            )
        }

        if(verbose == TRUE){
            message("Beginning optimization search -- performing ", nfit,
                    if(nfit == 1) " fit.\n" else " fits.\n")
        }

        opt.time <- system.time({
            if(n.cores > 1 && .Platform$OS.type != "windows"){
                fits <- parallel::mclapply(seq_len(nfit), fit_one, mc.cores=n.cores)
            } else {
                fits <- lapply(seq_len(nfit), fit_one)
            }
        })

        best.fit <- which.min(vapply(fits, function(x) x$objective, numeric(1)))
        out <- fits[[best.fit]]
    }

    if(model_is_formula){
        est.pars <- out$par
        names(est.pars) <- q_prep$par_names
    } else {
        est.pars <- exp(out$par)
        if(set.fog){
            fog.est <- est.pars[seq_len(n.fog)]
            est.pars <- est.pars[-seq_len(n.fog)]
        } else if(!is.numeric(tip.fog)){
            fog.est <- NULL
        }
    }

    if(is.null(prior)){
        loglik <- -as.numeric(out$objective)
        logpost <- NA_real_
    } else {
        loglik <- -as.numeric(RTMB_obj$report(out$par)$likelihood_nll)
        logpost <- -as.numeric(out$objective)
    }

    mcmc.fit <- NULL
    if(!is.null(mcmc)){
        if(!requireNamespace("tmbstan", quietly = TRUE)){
            stop("HMC requires the `tmbstan` package.", call. = FALSE)
        }
        mcmc <- do.call(mcmcControl, mcmc)

        ## start chains at the posterior mode, except for one-parameter
        ## models: tmbstan passes a length-1 init as a scalar where Stan
        ## expects a vector, and sampling fails ("no more scalars to read")
        if(length(out$par) == 1){
            warning("HMC for one-parameter models ignores `init` ",
                    if(is.null(mcmc$init)) "(posterior mode) " else "(user-supplied) ",
                    "and starts chains at random values (tmbstan limitation).",
                    call. = FALSE)
            mcmc$init <- NULL
        } else if(is.null(mcmc$init)){
            mcmc$init <- rep(list(unname(out$par)), mcmc$chains)
        }

        ## settings left unspecified must be omitted, not passed as NULL:
        ## rstan fails on e.g. seed = NULL and silently returns no samples
        mcmc.fit <- do.call(tmbstan::tmbstan, c(list(RTMB_obj), mcmc))
        if(length(mcmc.fit@sim$samples) == 0){
            stop("HMC sampling failed: no samples were returned.", call. = FALSE)
        }
    }

    if(model_is_formula){
        k <- length(est.pars)
        AIC <- if(is.null(prior)) -2*loglik + 2*k else NA_real_
        AICc <- if(is.null(prior)) -2*loglik + 2*k*(nb.tip/(nb.tip-k-1)) else NA_real_

        solution <- matrix(est.pars, ncol=1,
                           dimnames=list(names(est.pars), "estimate"))
        fun <- if(return.devfun) function(...) RTMB_obj else NA

        obj <- list(
            loglik=loglik,
            logpost=logpost,
            AIC=AIC,
            AICc=AICc,
            rate.cat=rate.cat,
            solution=solution,
            tip.fog.probs=NULL,
            index.mat=q_prep$Q_indicator,
            fog.vec=NULL,
            data=input.data,
            data.legend=data,
            phy=phy,
            states=NA,
            tip.states=NA,
            states.info=NA,
            iterations=out$iterations,
            convergence=out$convergence,
            message=out$message,
            collapse=collapse,
            root.p=root.p,
            args.list=list(
                phy=phy, data=data, model=model,
                q_prep=q_prep, trait_values=trait_values,
                p=est.pars
            ),
            call=call,
            devfun=fun,
            opt.time=opt.time,
            formula_model=TRUE,
            formula_list=formula_list,
            q_prep=q_prep,
            prior=prior,
            mcmc=mcmc.fit
        )
        class(obj) <- "corhmm"
        return(obj)
    }

    args.list <- list(
        phy=phy, liks=model.set.final$liks,
        Q=model.set.final$Q, rate=model.set.final$rate,
        root.p=root.p, rate.cat=rate.cat,
        order.test=order.test, lewis.asc.bias=lewis.asc.bias
    )

    if(node.states != "none" && verbose == TRUE){
        message("Finished. Inferring ancestral states using", node.states, "reconstruction.", "\n")
    }

    if(node.states == "marginal" || node.states == "scaled"){
        lik.anc <- ancRECON(
            phy, input.data, est.pars, rate.cat,
            rate.mat=rate.mat, method=node.states, ntraits=NULL,
            root.p=root.p, model=model,
            get.tip.states=get.tip.states, tip.fog=fog.est,
            collapse=collapse
        )
        phy$node.label <- apply(lik.anc$lik.anc.states, 1, which.max)
        tip.states <- lik.anc$lik.tip.states
        row.names(tip.states) <- phy$tip.label
    }
    if(node.states == "joint"){
        lik.anc <- ancRECON(
            phy, input.data, est.pars, rate.cat,
            rate.mat=rate.mat, method=node.states, ntraits=NULL,
            root.p=root.p, model=model,
            get.tip.states=get.tip.states, tip.fog=fog.est,
            collapse=collapse
        )
        phy$node.label <- lik.anc$lik.anc.states
        tip.states <- lik.anc$lik.tip.states
    }
    if(node.states == "none"){
        lik.anc <- list(lik.tip.states=NA, lik.anc.states=NA, info.anc.states=NA)
        phy$node.label <- NA
        tip.states <- NA
    }

    solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
    if(collapse){
        StateNames <- rep(gsub("_", "|", CorData$ObservedTraits), rate.cat)
        RCNames <- rep(paste("R", 1:rate.cat, sep=""), each=length(CorData$ObservedTraits))
    } else {
        StateNames <- rep(gsub("_", "|", CorData$PossibleTraits), rate.cat)
        RCNames <- rep(paste("R", 1:rate.cat, sep=""), each=length(CorData$PossibleTraits))
    }
    if(rate.cat > 1) StateNames <- paste(RCNames, StateNames)
    try({rownames(solution) <- colnames(solution) <- StateNames})

    if(set.fog == TRUE){
        tip.fog.probs <- numeric(length(model.set.final$fog.vec))
        tip.fog.probs[] <- c(fog.est, 0)[model.set.final$fog.vec]
        names(tip.fog.probs) <- StateNames
        model.set.final$np <- model.set.final$np + n.fog
    } else {
        tip.fog.probs <- NULL
    }

    if(is.null(prior)){
        AIC <- -2*loglik + 2*model.set.final$np
        AICc <- -2*loglik + 2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))
    } else {
        AIC <- AICc <- NA_real_
    }

    if(is.character(node.states) && (node.states == "marginal" || node.states == "scaled")){
        try({colnames(lik.anc$lik.anc.states) <- StateNames})
        try({colnames(tip.states) <- StateNames})
    }

    args.list <- c(list(p=log(est.pars)), args.list)

    if(loglik == -1e+06){
        warning("corHMM may have failed to optimize correctly, consider checking inputs and running again.", immediate. = TRUE)
    }

    if(!return.devfun){
        fun <- NA
    } else {
        fun <- function(...){
            L <- list(...)
            for(n in names(L)) assign(n, L[[n]])
            if(use_RTMB) return(RTMB_obj)
            dev.corhmm(p, phy, liks, Q, rate, root.p, rate.cat,
                       order.test, lewis.asc.bias)
        }
        environment(fun) <- list2env(list(
            p=log(est.pars), phy=phy,
            liks=model.set.final$liks,
            Q=model.set.final$Q,
            rate=model.set.final$rate,
            root.p=root.p, rate.cat=rate.cat,
            order.test=order.test,
            lewis.asc.bias=lewis.asc.bias
        ))
    }

    obj <- list(
        loglik=loglik,
        logpost=logpost,
        AIC=AIC,
        AICc=AICc,
        rate.cat=rate.cat,
        solution=solution,
        tip.fog.probs=tip.fog.probs,
        index.mat=model.set.final$index.matrix,
        fog.vec=model.set.final$fog.vec,
        data=input.data,
        data.legend=data.legend,
        phy=phy,
        states=lik.anc$lik.anc.states,
        tip.states=tip.states,
        states.info=lik.anc$info.anc.states,
        iterations=out$iterations,
        convergence=out$convergence,
        message=out$message,
        collapse=collapse,
        root.p=root.p,
        args.list=args.list,
        call=call,
        devfun=fun,
        opt.time=opt.time,
        formula_model=FALSE,
        prior=prior,
        mcmc=mcmc.fit
    )

    if(!is.null(input_rate_mat_states[1])){
        colnames(obj$solution) <- rownames(obj$solution) <- input_rate_mat_states
    }

    class(obj) <- "corhmm"
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################

dev.corhmm <- function(p, phy, liks, Q, rate, root.p, rate.cat, order.test, lewis.asc.bias, set.fog=FALSE, fog.vec) {
  
  p <- exp(p)
  cp_root.p <- root.p
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  #Obtain an object of all the unique ancestors
  edge1 <- phy$edge[,1]
  edge2 <- phy$edge[,2]
  edge.length <- phy$edge.length
  anc <- unique(edge1)
  #Map every edge to its ancestor in one pass instead of rescanning phy$edge per node
  desRowsList <- getDesRows(edge1, anc)
  k.rates <- dim(Q)[2] / 2
  if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)

  #Sets up the liks matrix to account for tip fog:
  if(set.fog == TRUE){
	  tip.fog.tmp <- p[1:length(unique(fog.vec))]
	  p <- p[-c(1:length(unique(fog.vec)))]
	  tip.fog <- numeric(length(fog.vec))
	  tip.fog[] <- c(tip.fog.tmp, 0)[fog.vec]
	  #Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
	  liks <- applyTipFog(liks, nb.tip, tip.fog, if(rate.cat > 1) rate.cat else 1)
  }

  Q[] <- c(p, 0)[rate]
  diag(Q) <- -rowSums(Q)

  # an error check for rates that are really weird
  test_mat <- expm(Q, method=c("Ward77"))
  if(any(round(rowSums(test_mat))> 2)){
    return(1000000)
  }
  # Q is fixed for the rest of this call, so decompose it once and reuse it for
  # every branch rather than exponentiating per edge (validated against test_mat).
  expm.funs <- makeExpmFuns(Q, reference = test_mat, ref.t = 1)
  Pv <- expm.funs$Pv

  # # if the q matrix has columns not estimated, remove them
  # row2rm <- apply(rate, 1, function(x) all(x == max(rate)))
  # col2rm <- apply(rate, 2, function(x) all(x == max(rate)))
  # Q.root <- Q[!row2rm | !col2rm, !row2rm | !col2rm]
  if(is.character(root.p)){
      if(root.p == "yang"){
          root.test <- Null(Q)
          if(dim(root.test)[2]>1){
              return(1000000)
          }
      }
  }

  if(order.test == TRUE){
      # ensure that the rate classes have mean rates in a consistent order (A > B > C > n)
      StateOrderMat <- matrix(1, (dim(Q)/rate.cat)[1], (dim(Q)/rate.cat)[2])
      RateClassOrderMat <- matrix(0, rate.cat, rate.cat)
      diag(RateClassOrderMat) <- 1:rate.cat
      OrderMat <- RateClassOrderMat %x% StateOrderMat
      Rate01 <- vector("numeric", rate.cat)
      for(i in 1:rate.cat){
          tmp <- Q[OrderMat == i]
          Rate01[i] <- tmp[tmp>=0][1]
      }
      OrderTest <- all.equal(Rate01, sort(Rate01, decreasing = TRUE))
      if(OrderTest != TRUE){
          return(1000000)
      }
  }

  ## FIXME: can we do this part in TMB/Rcpp?
  has.node.label <- !is.null(phy$node.label)
  for (i in seq(from = 1, length.out = nb.node)) {
      #the ancestral node at row i is called focal
      focal <- anc[i]
      #Get descendant information of focal
      desRows <- desRowsList[[i]]
      desNodes <- edge2[desRows]
      v <- 1
      #Loops through all descendants of focal (how we deal with polytomies):
      for (desIndex in seq_along(desRows)){
          v <- v*Pv(edge.length[desRows[desIndex]], liks[desNodes[desIndex],])
      }

      ##Allows for fixed nodes based on user input tree.
      if(has.node.label){
          if(!is.na(phy$node.label[focal - nb.tip])){
              fixer.tmp <- numeric(dim(Q)[2]/rate.cat)
              fixer.tmp[phy$node.label[focal - nb.tip]] <- 1
              fixer <- rep(fixer.tmp, rate.cat)
              v <- v * fixer
          }
      }

      #Sum the likelihoods:
      comp[focal] <- sum(v)
      #Divide each likelihood by the sum to obtain probabilities:
      liks[focal, ] <- v/comp[focal]
  }

  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  node.loglik <- sum(log(comp[-TIPS]))
  if (is.na(node.loglik)){return(1000000)}
  liks.root <- liks[root,]

  if (is.null(root.p)){
      equil.root <- getEquilRoot(Q)
      flat.root = equil.root
      k.rates <- 1/length(which(!is.na(equil.root)))
      flat.root[!is.na(flat.root)] = k.rates
      flat.root[is.na(flat.root)] = 0
      loglik<- -(node.loglik + log(sum(flat.root * liks.root)))
  }
  if(is.character(root.p)){
      # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
      if(root.p == "yang"){
          #root.test was already computed above from the same Q
          root.p <- c(root.test/sum(root.test))
          loglik <- -(node.loglik + log(sum(root.p * liks.root)))
          if(is.infinite(loglik)){
              return(1000000)
          }
      }else{
          # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
          root.p = liks.root / sum(liks.root)
          loglik <- -(node.loglik + log(sum(exp(log(root.p)+log(liks.root)))))
      }
  }else{
    if(is.numeric(root.p[1])){
      loglik <- -(node.loglik + log(sum(exp(log(root.p)+log(liks.root)))))
      if(is.infinite(loglik)){
        return(1000000)
      }
    }
  }
  # root.p!==NULL will fix root probabilities based on user supplied vector:
  if(lewis.asc.bias == TRUE){
    p <- log(p)
    dummy.liks.vec <- getLewisLikelihood(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = cp_root.p, rate.cat = rate.cat)
    loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
  }
  return(loglik)
}




######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################

getLewisLikelihood <- function(p, phy, liks, Q, rate, root.p, rate.cat){
  nTips <- length(phy$tip.label)
  nNodes <- length(unique(phy$edge[,1]))
  nStates <- dim(Q)[1]/rate.cat
  states_structure <- seq(from = 1, by = nStates, length.out = rate.cat)
  dummy.liks.vec <- vector("numeric", nStates)
  for(state_i in 1:nStates){
    lik_structure_i <- states_structure + (state_i - 1)
    liks[] <- 0
    liks[1:nTips, lik_structure_i] <- 1
    dummy.liks.vec[state_i] <- -dev.corhmm(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE)
  }
  return(dummy.liks.vec)
}

# JDB modified functions
rate.cat.set.corHMM.JDB <- function(phy,data,rate.cat, ntraits, model, rate.mat=NULL, collapse=TRUE, CorData=NULL, matching=NULL) {
    obj <- NULL
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    obj$rate.cat<-rate.cat
    # corProcessData / match.tree.data have usually already been run by the
    # caller on this same data; reuse those results when they are supplied.
    if(is.null(CorData)){
      CorData <- corProcessData(data, collapse = collapse)
    }
    if(is.null(rate.mat)){
      rate <- getStateMat4Dat(data, model, collapse = collapse, CorData = CorData)$rate.mat
      if(rate.cat > 1){
        StateMats <- vector("list", rate.cat)
        for(i in 1:rate.cat){
          StateMats[[i]] <- rate
        }
        rate <- getFullMat(StateMats)
      }
    }else{
      rate <- rate.mat
      nTraits <- dim(rate)[1]/rate.cat
    }
    nTraits <- dim(rate)[1]/rate.cat
    rate[rate == 0] <- NA
    index.matrix<-rate
    rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1

    data <- CorData$corData
    nObs <- length(CorData$ObservedTraits)

    if(is.null(matching)){
      matching <- match.tree.data(phy,data)
    }
    data <- matching$data

    # this is no longer needed since corProcessData will produce a dataset of a specific type every time
    # x <- as.numeric(data.sort[,1])
    # TIPS <- 1:nb.tip
    # if(min(x) !=0){
    #     x <- x - min(x)
    # }
    # this is being removed temporarily - this handles NA tip states
    # for(i in 1:nb.tip){
    #     if(is.na(x[i])){x[i]=2}
    # }
    
    tmp <- matrix(0, nb.tip + nb.node, nTraits)
    for(i in 1:nb.tip){
        focal_state <- matching$data[i,2]
        if(focal_state == "?"){
          tmp[i, ] <- 1
        }else{
          state_index <- as.numeric(unlist(strsplit(as.character(focal_state), "&")))
          tmp[i, state_index] <- 1
        }
    }
    liks <- matrix(rep(tmp, rate.cat), nb.tip + nb.node, nTraits*rate.cat)
    
    Q <- matrix(0, dim(rate)[1], dim(rate)[1])
	
    obj$np <- max(rate) - 1
    obj$rate <- rate
    obj$index.matrix <- index.matrix
    obj$liks <- liks
    obj$Q <- Q
    
    return(obj)
}


corProcessData <- function(data, rate.mat=NULL, collapse=FALSE){
  nCol <- dim(data)[2]
  LevelList <- StateMats <- vector("list", nCol-1)
  # detect the number of states in each column. & is treated as indicating polymorphism. ? is treated as unknown data.
  for(i in 2:nCol){
    # 4/2/25 Caetano: Fix bug that causes "States_i" to become NULL.
    data_tmp <- data[,i]
    if(!is.factor(data_tmp)){
      data_tmp <- as.factor(data_tmp)
    }
    States_i <- levels(data_tmp)
    if(any(States_i == "?")){
      States_i <- States_i[!States_i == "?"]
    }
    if(length(grep("&", States_i)) > 0){
      States_i <- unique(unlist(strsplit(States_i, "&")))
    }
    
    StateMats[[i-1]] <- getRateCatMat(length(States_i))
    LevelList[[i-1]] <- States_i
    # if(any(is.na(suppressWarnings(as.numeric(States_i))))){
    #   LevelList[[i-1]] <- sort(States_i)
    # }else{
    #   LevelList[[i-1]] <- States_i[sort(as.numeric(States_i), index.return=TRUE)$ix]
    # }
  }
  # identify the possible trait combinations
  TraitList <- expand.grid(LevelList, stringsAsFactors = FALSE)
  Traits <- do.call(paste, c(TraitList, sep = "_"))
  # convert each column into a numeric value associated with a member of the trait combinations. ? are associated with all values of that column, & indicates the combination of two or more
  # and symbolizes it can be any of the separated states
  col.patterns <- lapply(2:nCol, function(j)
    paste("(", gsub("&", "|", as.character(data[,j]), fixed = TRUE), ")", sep = ""))
  search.strings <- paste("^", do.call(paste, c(col.patterns, sep = "_")), "$", sep = "")
  # ? means it can be any of the states in that character
  search.strings <- gsub("(?)", ".*", search.strings, fixed=TRUE)
  # There are only ever a handful of distinct coding patterns no matter how many
  # taxa there are, so match each unique pattern against Traits once and reuse.
  uniq.strings <- unique(search.strings)
  string.index <- match(search.strings, uniq.strings)
  trait.hits <- lapply(uniq.strings, function(x) grep(x, Traits))
  observed.traits_index <- unlist(trait.hits, use.names = FALSE)
  ObservedTraits <- Traits[sort(unique(observed.traits_index))]
  if(collapse){
    obs.hits <- vapply(uniq.strings, function(x)
      paste(grep(x, ObservedTraits), collapse="&"), character(1), USE.NAMES = FALSE)
  }else{
    obs.hits <- vapply(trait.hits, paste, character(1), collapse = "&")
  }
  corData <- data.frame(sp = data[,1], d = obs.hits[string.index],
                        stringsAsFactors = FALSE)
  return(list(StateMats = StateMats,  PossibleTraits = Traits, ObservedTraits = ObservedTraits, corData = corData))
}


#' @export
#print.corhmm <- function(x,...){

#  ntips=Ntip(x$phy)
#  output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips, row.names="")
#  names(output)<-c("lnL","AIC","AICc","Rate.cat","ntax")
#  cat("\nFit\n")
#  print(output)
#  cat("\n")
#  UserStates <- colnames(x$solution)
#  ColNames <- paste0(colnames(x$data)[-1], collapse = "|")
  
#  cat("Legend\n")
#  print(ColNames)
#  print(UserStates)
#  cat("\n")
  
#  param.est<- x$solution
#  cat("Rates\n")
#  print(param.est)
#  cat("\n")

#	if(!is.null(x$tip.fog.probs)){
#               cat("Tip fog\n")
#		print(x$tip.fog.probs)
#		cat("\n")
#	}
	
#    if(any(x$eigval<0)){
#        index.matrix <- x$index.mat
        #If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
#        if (any(x$eigval<0)) {
#            warning("The objective function may be at a saddle point", "\n")
#        }
#    }
#    else{
#        cat("Arrived at a reliable solution","\n")
#    }
#}

#' @export
print.corhmm <- function(x,...){
    if (isTRUE(x$formula_model)) {
        ntips <- Ntip(x$phy)
        output <- data.frame(x$loglik, x$AIC, x$AICc, x$rate.cat, ntips, row.names = "")
        names(output) <- c("lnL", "AIC", "AICc", "Rate.cat", "ntax")
        cat("\nFit\n")
        print(output)
        cat("\nFormula coefficients\n")
        print(x$solution)
        cat("\n")
        return(invisible(x))
    }
	if (!is.null(x$phy_list)){
		ntips=Ntip(x$phy)
		if(is.null(x$BIC)) {
			output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips,x$ntrees, row.names="")
			names(output)<-c("lnL","AIC","AICc","Rate.cat","ntax","ntrees")
		} else {
			output<-data.frame(x$loglik,x$AIC,x$AICc,x$BIC,x$rate.cat,ntips,x$ntrees, row.names="")
			names(output)<-c("lnL","AIC","AICc","BIC","Rate.cat","ntax","ntrees")
		}
		cat("\nFit\n")
		print(output)
		cat("\n")
	} else {
		ntips=Ntip(x$phy)
		if(is.null(x$BIC)) {
			output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips, row.names="")
			names(output)<-c("lnL","AIC","AICc","Rate.cat","ntax")
		} else {
			output<-data.frame(x$loglik,x$AIC,x$AICc,x$BIC,x$rate.cat,ntips, row.names="")
			names(output)<-c("lnL","AIC","AICc","BIC","Rate.cat","ntax")
		}
		cat("\nFit\n")
		print(output)
		cat("\n")
	}

	UserStates <- gsub("_", "|", corProcessData(x$data)$PossibleTraits)
	ColNames <- paste0(colnames(x$data)[-1], collapse = "|")

	cat("Legend\n")
	print(ColNames)
	print(UserStates)
	cat("\n")
	
	param.est<- x$solution
	cat("Rates\n")
	print(param.est)
	cat("\n")

	if(!is.null(x$tip.fog.probs)){
		cat("Tip fog\n")
		print(x$tip.fog.probs)
		cat("\n")
	}
	
	if(any(x$eigval<0)){
		index.matrix <- x$index.mat
		#If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
		if (any(x$eigval<0)) {
			cat("The objective function may be at a saddle point", "\n")
		}
	}
	else{
		cat("Arrived at a reliable solution","\n")
	}

}


# function for calculating PIR according to gardner and organ (2021)
getPIR <- function(phy, data, collapse){
  # calculate the CI score
  CorData <- corProcessData(data, collapse = collapse)
  dat <- CorData$corData
  matching <- match.tree.data(phy,dat)
  dat <- matching$data
  phy <- matching$phy
  data.sort <- data.frame(dat[,2], dat[,2],row.names=dat[,1])
  data.sort <- data.sort[phy$tip.label,]
  data.sort <- as.matrix(data.sort)
  dat.phangorn <- phyDat(data.sort,type="USER", levels=1:max(as.numeric(dat[,2])))
  phy.tmp <- multi2di(phy)
  par.score <- parsimony(phy.tmp, dat.phangorn, method="fitch")/2
  ci.score <- CI(phy.tmp, dat.phangorn)

  # calcualte the NIR
  state_table <- table(dat[,2])
  levels <- 1:max(as.numeric(dat[,2]))
  # if all levels are observed, then we use the min from the observed states, otherwise the min is just 0
  if(dim(state_table) == length(levels)){
    nir.score <- (max(state_table) - min(state_table))/length(phy$tip.label)
  }else{
    nir.score <- (max(state_table) - 0)/length(phy$tip.label)
  }
  pir.score <- ci.score * nir.score
  return(pir.score)
}


######################################################################################################################################
######################################################################################################################################
### ORIGINAL CODE THAT IS NO LONGER USED
######################################################################################################################################
######################################################################################################################################

#Generalized ace() function that allows analysis to be carried out when there are polytomies:
dev.corhmm.ORIGINAL <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test) {
    p = exp(p)
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    TIPS <- 1:nb.tip
    comp <- numeric(nb.tip + nb.node)
    phy <- reorder(phy, "pruningwise")
    #Obtain an object of all the unique ancestors
    anc <- unique(phy$edge[,1])
    k.rates <- dim(Q)[2] / 2
    if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)

    Q[] <- c(p, 0)[rate]
    diag(Q) <- -rowSums(Q)
    for (i  in seq(from = 1, length.out = nb.node)) {
        #the ancestral node at row i is called focal
        focal <- anc[i]
        #Get descendant information of focal
        desRows<-which(phy$edge[,1]==focal)
        desNodes<-phy$edge[desRows,2]
        v <- 1
        #Loops through all descendants of focal (how we deal with polytomies):
        for (desIndex in sequence(length(desRows))){
            v<-v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
        }
        #Sum the likelihoods:
        comp[focal] <- sum(v)
        #Divide each likelihood by the sum to obtain probabilities:
        liks[focal, ] <- v/comp[focal]
    }
    #Temporary solution for ensuring an ordered Q with respect to the rate classes. If a simpler model is called this feature is automatically turned off:
    par.order<-NA
    if(k.rates == 2){
        try(par.order <- p[3] > p[8])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 3){
        try(par.order <- p[3] > p[9] | p[9] > p[14])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 4){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[20])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 5){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[26])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 6){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[27] | p[27] > p[32])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 7){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[27] | p[27] > p[33] | p[33] > p[38])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    #Specifies the root:
    root <- nb.tip + 1L
    #If any of the logs have NAs restart search:
    if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
    else{
        equil.root <- NULL
        for(i in 1:ncol(Q)){
            posrows <- which(Q[,i] >= 0)
            rowsum <- sum(Q[posrows,i])
            poscols <- which(Q[i,] >= 0)
            colsum <- sum(Q[i,poscols])
            equil.root <- c(equil.root,rowsum/(rowsum+colsum))
        }
        if (is.null(root.p)){
            flat.root = equil.root
            k.rates <- 1/length(which(!is.na(equil.root)))
            flat.root[!is.na(flat.root)] = k.rates
            flat.root[is.na(flat.root)] = 0
            loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
        }
        else{
            if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
                if(root.p == "yang"){
                    root.p <- Null(Q)
                    root.p <- c(root.p/sum(root.p))
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                    if(is.infinite(loglik)){
                        return(1000000)
                    }
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks[root,] / sum(liks[root,])
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                }
            }
            # root.p!==NULL will fix root probabilities based on user supplied vector:
            else{
                loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                if(is.infinite(loglik)){
                    return(1000000)
                }
            }
        }
    }
    loglik
}

rate.cat.set.corHMM <- function (phy, data.sort, rate.cat) {
    k = 2
    obj <- NULL
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    obj$rate.cat <- rate.cat
    rate <- rate.mat.maker(hrm = TRUE, rate.cat = rate.cat)
    index.matrix <- rate
    rate[is.na(rate)] <- max(rate, na.rm = TRUE) + 1
    x <- data.sort[, 1]
    TIPS <- 1:nb.tip
    for (i in 1:nb.tip) {
        if (is.na(x[i])) {
            x[i] = 2
        }
    }
    if (rate.cat == 1) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        TIPS <- 1:nb.tip
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, 1] = 1
            }
            if (x[i] == 1) {
                liks[i, 2] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:2] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 2) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:4] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 3) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:6] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 4) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:8] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 5) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7, 9)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8, 10)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:10] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 6) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7, 9, 11)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8, 10, 12)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:12] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 7) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7, 9, 11, 13)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8, 10, 12, 14)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:14] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    obj$np <- max(rate) - 1
    obj$rate <- rate
    obj$index.matrix <- index.matrix
    obj$liks <- liks
    obj$Q <- Q
    obj
}


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
