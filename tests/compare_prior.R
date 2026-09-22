library(ape)
library(rphylo)
library(btw)
library(coda)
library(rstan)
library(ramcmc)

pkgload::load_all(".", quiet = TRUE)
rstan::rstan_options(auto_write = TRUE)

## -------------------------------------------------------------------------
## Settings
## -------------------------------------------------------------------------

tree_sizes <- c(100L, 250L, 500L)
sim_seeds <- 900L                 # pilot; later use 900:903 or 900:909
prior_sds <- c(1.5, 3, 6)        # tighter, target, wider
prior_mu <- log(0.1)
prior_range <- 3

n_chains <- 4L
chain_seeds <- 427:430
n_keep <- 2000L

ram_burnin <- 5000L
ram_thin <- 10L
ram_iter <- ram_burnin + n_keep * ram_thin

hmc_iter <- 2L * n_keep
hmc_warmup <- hmc_iter / 2L
hmc_cores <- min(n_chains, parallel::detectCores())

mh_burnin <- 100000L
mh_sample <- 500L
mh_post_iter <- n_keep * mh_sample
bt_total_iter <- mh_burnin + mh_post_iter

options(mc.cores = hmc_cores)

out_dir <- file.path("results", "prior_stability")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bt_path <- Sys.getenv("BAYESTRAITS_HOME")
bt_bin <- Sys.getenv("BAYESTRAITS_BIN", "BayesTraitsV5")
options(bt_path = bt_path, bt_bin = bt_bin)

## -------------------------------------------------------------------------
## Known generating model: two binary traits, eight rates
## -------------------------------------------------------------------------

rate_names <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
true_rates <- c(q12=.08, q13=.18, q21=.06, q24=.25,
                q31=.14, q34=.07, q42=.20, q43=.10)

Q <- matrix(0, 4, 4)
Q[1,2] <- true_rates["q12"]; Q[1,3] <- true_rates["q13"]
Q[2,1] <- true_rates["q21"]; Q[2,4] <- true_rates["q24"]
Q[3,1] <- true_rates["q31"]; Q[3,4] <- true_rates["q34"]
Q[4,2] <- true_rates["q42"]; Q[4,3] <- true_rates["q43"]
diag(Q) <- -rowSums(Q)

state_T1 <- c(0, 0, 1, 1)
state_T2 <- c(0, 1, 0, 1)

summary_all <- list()
diag_all <- list()
pairwise_all <- list()
k_summary <- k_diag <- k_pair <- 0L

## -------------------------------------------------------------------------
## Main experiment
## -------------------------------------------------------------------------

for (n in tree_sizes) {
  for (sim_seed in sim_seeds) {

    ## Generate one tree/data set; reuse it for every prior
    set.seed(sim_seed)
    phy <- ape::rtree(n)

    depth <- ape::node.depth.edgelength(phy)[seq_len(n)]
    phy$edge.length <- phy$edge.length / mean(depth)

    set.seed(sim_seed + 100000L)
    joint <- ape::rTraitDisc(
      phy, model = Q, states = c("00","01","10","11"),
      root.value = sample(1:4, 1)
    )

    joint_id <- match(as.character(joint), c("00","01","10","11"))
    dat <- data.frame(species=phy$tip.label,
                      T1=state_T1[joint_id],
                      T2=state_T2[joint_id])

    rphylo_trees <- structure(list(phy), class="multiPhylo")
    rphylo_traits <- data.frame(T1=dat$T1, T2=dat$T2, row.names=dat$species)
    btw_traits <- data.frame(species=dat$species,
                             T1=as.character(dat$T1),
                             T2=as.character(dat$T2))

    for (prior_sd in prior_sds) {

      prior_id <- c("tighter","target","wider")[match(prior_sd, prior_sds)]
      prior_lb <- prior_mu - prior_range * prior_sd
      prior_ub <- prior_mu + prior_range * prior_sd

      cor_prior <- list(
        lb = prior_lb, ub = prior_ub, range = prior_range,
        coef.sd = sqrt(2) * prior_sd, target = "rates"
      )
      rphylo_prior <- list(
        type="lognormal", location=prior_mu, scale=prior_sd
      )
      btw_prior <- paste("PriorAll lognormal", prior_mu, prior_sd)

      cat("\n============================================================\n")
      cat("n =", n, "| seed =", sim_seed, "| prior =", prior_id,
          "| SD =", prior_sd, "\n")
      cat("============================================================\n")

      ## ------------------------------------------------------------------
      ## Build one RTMB posterior objective
      ## ------------------------------------------------------------------

      cat("Building RTMB objective ...\n")

      fit_obj <- corHMM(
        phy=phy, data=dat, rate.cat=1, model="ARD",
        node.states="none", root.p=rep(.25,4), ip=rep(.1,8),
        collapse=FALSE, lower.bound=1e-12, upper.bound=1e8,
        use_RTMB=TRUE, prior=cor_prior,
        return.devfun=TRUE, verbose=FALSE
      )

      obj <- fit_obj$devfun()
      start <- obj$env$last.par.best

      from <- c("0|0","0|0","0|1","0|1","1|0","1|0","1|1","1|1")
      to   <- c("0|1","1|0","0|0","1|1","0|0","1|1","0|1","1|0")
      rate_index <- fit_obj$index.mat[cbind(
        match(from, rownames(fit_obj$index.mat)),
        match(to, colnames(fit_obj$index.mat))
      )]

      ## ------------------------------------------------------------------
      ## 1. RTMB adaptive Metropolis-Hastings
      ## ------------------------------------------------------------------

      cat("RTMB-MH ...\n")

      ram_draws <- vector("list", n_chains)
      ram_seconds <- numeric(n_chains)
      ram_accept <- numeric(n_chains)

      for (ch in seq_len(n_chains)) {
        set.seed(chain_seeds[ch])

        tm <- system.time({
          x <- start + rnorm(length(start), 0, .05)
          lp <- -obj$fn(x)
          S <- diag(.1, length(start))
          keep <- matrix(NA_real_, n_keep, length(start))
          accepted <- 0L
          save_i <- 0L

          for (iter in seq_len(ram_iter)) {
            u <- rnorm(length(start))
            prop <- as.numeric(x + S %*% u)
            lp_prop <- -obj$fn(prop)
            alpha <- if (is.finite(lp_prop)) exp(min(0, lp_prop-lp)) else 0

            if (runif(1) < alpha) {
              x <- prop
              lp <- lp_prop
              if (iter > ram_burnin) accepted <- accepted + 1L
            }

            if (iter <= ram_burnin)
              S <- ramcmc::adapt_S(S, u, alpha, iter, target=.234)

            if (iter > ram_burnin && (iter-ram_burnin) %% ram_thin == 0) {
              save_i <- save_i + 1L
              keep[save_i,] <- x
            }
          }
        })

        ram_seconds[ch] <- unname(tm["elapsed"])
        ram_accept[ch] <- accepted / (n_keep * ram_thin)
        ram_draws[[ch]] <- exp(keep[, rate_index, drop=FALSE])
        colnames(ram_draws[[ch]]) <- rate_names
      }

      stopifnot(all(vapply(ram_draws, function(x) all(is.finite(x)), logical(1))))
      saveRDS(ram_draws, file.path(out_dir,
              sprintf("ram_n%d_seed%d_%s.rds", n, sim_seed, prior_id)))

      ## ------------------------------------------------------------------
      ## 2. RTMB-HMC: exactly the same RTMB objective
      ## ------------------------------------------------------------------

      cat("RTMB-HMC ...\n")

      tm <- system.time({
        fit_hmc <- tmbstan::tmbstan(
          obj, chains=n_chains, cores=hmc_cores,
          iter=hmc_iter, warmup=hmc_warmup,
          seed=chain_seeds[1],
          init=rep(list(start), n_chains),
          control=list(adapt_delta=.9, max_treedepth=12)
        )
      })
      hmc_seconds <- unname(tm["elapsed"])

      arr <- rstan::extract(fit_hmc, permuted=FALSE, inc_warmup=FALSE)
      hmc_draws <- lapply(seq_len(n_chains), function(ch) {
        x <- exp(matrix(arr[,ch,rate_index,drop=FALSE], ncol=8))
        colnames(x) <- rate_names
        x
      })

      sp <- rstan::get_sampler_params(fit_hmc, inc_warmup=FALSE)
      hmc_div <- sum(vapply(sp, function(x) sum(x[,"divergent__"]), numeric(1)))
      hmc_depth <- sum(vapply(sp, function(x) sum(x[,"treedepth__"] >= 12), numeric(1)))

      stopifnot(all(vapply(hmc_draws, function(x) all(is.finite(x)), logical(1))))
      saveRDS(fit_hmc, file.path(out_dir,
              sprintf("hmc_n%d_seed%d_%s.rds", n, sim_seed, prior_id)))

      ## ------------------------------------------------------------------
      ## 3. rphylo
      ## ------------------------------------------------------------------

      cat("rphylo ...\n")

      rphylo_draws <- vector("list", n_chains)
      rphylo_seconds <- numeric(n_chains)

      for (ch in seq_len(n_chains)) {
        tm <- system.time({
          fit_r <- rphylo::run_dep_model(
            trees=rphylo_trees, charfile=rphylo_traits,
            burnin=mh_burnin, iterations=mh_post_iter,
            sample_freq=mh_sample, prior=rphylo_prior,
            revjump=FALSE, pis=rep(1,4),
            seed=chain_seeds[ch], verbose=FALSE
          )
        })
        rphylo_seconds[ch] <- unname(tm["elapsed"])
        rphylo_draws[[ch]] <- as.matrix(fit_r$rates)
        colnames(rphylo_draws[[ch]]) <- rate_names
      }

      stopifnot(all(vapply(rphylo_draws, function(x) all(is.finite(x)), logical(1))))
      saveRDS(rphylo_draws, file.path(out_dir,
              sprintf("rphylo_n%d_seed%d_%s.rds", n, sim_seed, prior_id)))

      ## ------------------------------------------------------------------
      ## 4. BayesTraits fixed-model MH
      ## ------------------------------------------------------------------

      cat("BayesTraits ...\n")

      btw_draws <- vector("list", n_chains)
      btw_seconds <- numeric(n_chains)

      for (ch in seq_len(n_chains)) {
        old_wd <- setwd(bt_path)

        tm <- system.time({
          fit_bt <- btw::bayestraits(
            btw_traits, phy,
            c("3", "2",
              paste("Iterations", bt_total_iter),
              paste("Burnin", mh_burnin),
              paste("Sample", mh_sample),
              "Pis none", btw_prior,
              paste("Seed", chain_seeds[ch])),
            silent=TRUE, remove_files=TRUE
          )
        })

        setwd(old_wd)
        btw_seconds[ch] <- unname(tm["elapsed"])
        btw_draws[[ch]] <- as.matrix(
          fit_bt$Log$results[,rate_names,drop=FALSE]
        )
      }

      stopifnot(all(vapply(btw_draws, function(x) all(is.finite(x)), logical(1))))
      saveRDS(btw_draws, file.path(out_dir,
              sprintf("btw_n%d_seed%d_%s.rds", n, sim_seed, prior_id)))

      ## ------------------------------------------------------------------
      ## Diagnostics
      ## ------------------------------------------------------------------

      draws <- list(
        `RTMB-MH`=ram_draws,
        `RTMB-HMC`=hmc_draws,
        rphylo=rphylo_draws,
        BayesTraits=btw_draws
      )

      wall <- c(sum(ram_seconds), hmc_seconds,
                sum(rphylo_seconds), sum(btw_seconds))

      ## Sequential MH/rphylo/BT; four parallel HMC chains.
      ## On Slurm replace HMC estimate with TotalCPU when available.
      cpu <- c(sum(ram_seconds), hmc_seconds*hmc_cores,
               sum(rphylo_seconds), sum(btw_seconds))

      method_diag <- list()

      for (method in names(draws)) {
        ml <- coda::mcmc.list(lapply(draws[[method]], coda::mcmc))
        rhat <- coda::gelman.diag(
          ml, multivariate=FALSE, autoburnin=FALSE
        )$psrf[,"Point est."]
        ess <- coda::effectiveSize(ml)
        pooled <- do.call(rbind, draws[[method]])
        post_mean <- colMeans(pooled)

        k_diag <- k_diag + 1L
        diag_all[[k_diag]] <- data.frame(
          n=n, seed=sim_seed, prior=prior_id, prior_sd=prior_sd,
          method=method, parameter=rate_names,
          mean=post_mean, rhat=unname(rhat[rate_names]),
          ess=unname(ess[rate_names]),
          cpu_seconds_per_ess=cpu[match(method,names(draws))] /
                              unname(ess[rate_names])
        )

        k_summary <- k_summary + 1L
        summary_all[[k_summary]] <- data.frame(
          n=n, seed=sim_seed, prior=prior_id, prior_sd=prior_sd,
          method=method,
          wall_seconds=wall[match(method,names(draws))],
          cpu_seconds=cpu[match(method,names(draws))],
          rmse_truth=sqrt(mean((post_mean-true_rates)^2)),
          max_rhat=max(rhat, na.rm=TRUE),
          min_ess=min(ess, na.rm=TRUE),
          max_cpu_seconds_per_ess=max(
            cpu[match(method,names(draws))] / ess, na.rm=TRUE
          ),
          acceptance=if (method=="RTMB-MH") mean(ram_accept) else NA_real_,
          divergences=if (method=="RTMB-HMC") hmc_div else NA_integer_,
          treedepth_hits=if (method=="RTMB-HMC") hmc_depth else NA_integer_
        )

        method_diag[[method]] <- post_mean
      }

      ## Pairwise posterior-mean matching
      cmb <- combn(names(method_diag), 2)
      for (j in seq_len(ncol(cmb))) {
        a <- cmb[1,j]; b <- cmb[2,j]
        k_pair <- k_pair + 1L
        pairwise_all[[k_pair]] <- data.frame(
          n=n, seed=sim_seed, prior=prior_id, prior_sd=prior_sd,
          method1=a, method2=b,
          mean_rmse=sqrt(mean((method_diag[[a]]-method_diag[[b]])^2)),
          log_mean_rmse=sqrt(mean(
            (log(method_diag[[a]])-log(method_diag[[b]]))^2
          ))
        )
      }

      ## Save immediately after every complete condition
      write.csv(do.call(rbind, summary_all),
                file.path(out_dir, "summary.csv"), row.names=FALSE)
      write.csv(do.call(rbind, diag_all),
                file.path(out_dir, "diagnostics.csv"), row.names=FALSE)
      write.csv(do.call(rbind, pairwise_all),
                file.path(out_dir, "pairwise_matching.csv"), row.names=FALSE)

      cat("Finished:", n, sim_seed, prior_id, "\n")
    }
  }
}

## -------------------------------------------------------------------------
## Final compact summaries
## -------------------------------------------------------------------------

summary_df <- do.call(rbind, summary_all)
diag_df <- do.call(rbind, diag_all)
pairwise_df <- do.call(rbind, pairwise_all)

print(summary_df[,c(
  "n","seed","prior","method","wall_seconds","rmse_truth",
  "max_rhat","min_ess","max_cpu_seconds_per_ess",
  "acceptance","divergences","treedepth_hits"
)])

print(pairwise_df)