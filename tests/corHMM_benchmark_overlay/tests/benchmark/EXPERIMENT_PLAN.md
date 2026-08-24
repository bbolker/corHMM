# Experiment plan implemented by the benchmark scripts

## Simulation maximum likelihood

Data source: random binary traits on subtrees of the complete 234-tip
ray-finned fish tree.

Default factors:

```text
ntips: 50, 100, 150, 234
ntraits: 2, 3
seeds: 900:999
```

Methods:

```text
original corHMM
corHMM-RTMB
corHMM formula interface
```

Recorded quantities:

```text
wall time
optimizer time
optimized log-likelihood
transition rates
convergence code
optimizer iterations
peak memory from the Slurm resource log
```

Post-processing:

```text
pairwise log-likelihood differences
named-rate RMSE
maximum absolute rate difference
tie-aware NLL winner, tolerance 1e-6
```

## Simulation fixed-model Bayesian comparison

Only two-trait datasets are used because rphylo and BayesTraits Discrete are
pairwise binary-trait models.

Methods:

```text
RTMB robust adaptive Metropolis-Hastings
RTMB-HMC through tmbstan
rphylo fixed-model MCMC
BayesTraits fixed-model MCMC through bbolker/btw
```

Common prior:

```text
lower = log(1)
upper = log(100 * ntips)
mean = (lower + upper) / 2
sd = (upper - lower) / 3
```

Recorded quantities:

```text
rate draws
native likelihood or posterior-target output supplied by each method
posterior mean, median, SD, and 95% interval
R-hat
bulk ESS
tail ESS
MCSE
wall time
ESS per second
sampler-specific diagnostics
peak memory from the Slurm resource log
```

## Separate RJMCMC extension

BayesTraits RJMCMC is submitted and summarized separately. It is not included
in the fixed-model rate-RMSE or ESS comparison.

## Primates

The same ML and fixed-model Bayesian methods are applied once to the package
primates dataset. ML Wald intervals are calculated outside corHMM from the
fitted objective Hessian. Full raw fit objects are retained.

## Ray-finned fish application

The complete `c_tree.rds` and `c_traits.rds` data are fitted with:

```r
list(
  ag ~ care * spawning,
  care ~ 1,
  spawning ~ 1
)
```

The main reported quantity is the male-care effect on accessory-gland gain,
averaged across spawning states and expressed as a rate ratio. Full raw fit
objects and diagnostics are retained.
