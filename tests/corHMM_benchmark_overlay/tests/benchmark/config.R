benchmark_design <- list(
  tree_sizes = c(50L, 100L, 150L, 234L),
  ntraits = c(2L, 3L),
  dataset_seeds = 900:999,
  root_state = "equal",
  tree_total_length = 1
)

benchmark_profiles <- list(
  smoke = list(
    chains = 2L,
    retained_draws = 100L,
    chain_seeds = 427:428,

    rtmb_mh_warmup = 200L,
    rtmb_mh_iterations = 1000L,
    rtmb_mh_thin = 5L,
    rtmb_mh_initial_scale = 0.1,

    hmc_iter = 200L,
    hmc_warmup = 100L,
    hmc_adapt_delta = 0.9,
    hmc_max_treedepth = 12L,

    rphylo_burnin = 100L,
    rphylo_iterations = 1000L,
    rphylo_sample_freq = 5L,

    btw_burnin = 1000L,
    btw_iterations = 11000L,
    btw_sample_freq = 50L,

    rj_burnin = 1000L,
    rj_iterations = 11000L,
    rj_sample_freq = 50L
  ),

  full = list(
    chains = 4L,
    retained_draws = 8000L,
    chain_seeds = 427:430,

    ## 8,000 retained samples per chain.
    rtmb_mh_warmup = 4000L,
    rtmb_mh_iterations = 80000L,
    rtmb_mh_thin = 10L,
    rtmb_mh_initial_scale = 0.1,

    hmc_iter = 10000L,
    hmc_warmup = 2000L,
    hmc_adapt_delta = 0.9,
    hmc_max_treedepth = 12L,

    rphylo_burnin = 100000L,
    rphylo_iterations = 800000L,
    rphylo_sample_freq = 100L,

    ## BayesTraits records at the burn-in boundary; post-processing keeps
    ## the last 8,000 draws from each chain.
    btw_burnin = 3000000L,
    btw_iterations = 11000000L,
    btw_sample_freq = 1000L,

    rj_burnin = 3000000L,
    rj_iterations = 11000000L,
    rj_sample_freq = 1000L
  )
)

benchmark_path <- function(...) {
  gsub("\\\\", "/", file.path(...))
}
