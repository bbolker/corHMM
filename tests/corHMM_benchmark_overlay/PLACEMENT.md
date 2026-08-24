# Placement map

Extract the overlay archive directly into the corHMM repository root, or copy
these files while preserving their relative paths.

```text
R/corHMM.R                         -> replace repository R/corHMM.R
R/RTMB_devcorhmm.R                 -> replace repository R/RTMB_devcorhmm.R

tests/benchmark/config.R           -> add
tests/benchmark/common.R           -> add
tests/benchmark/simulation.R       -> add
tests/benchmark/generate_simulation_data.R -> add
tests/benchmark/build_manifests.R  -> add
tests/benchmark/methods.R          -> add
tests/benchmark/run_job.R          -> add
tests/benchmark/run_real.R         -> add
tests/benchmark/run_local_primates.R -> add
tests/benchmark/run_local_rayfinned.R -> add
tests/benchmark/summarize_simulation.R -> add
tests/benchmark/summarize_real.R   -> add
tests/benchmark/collect_resources.R -> add
tests/benchmark/install_benchmark_packages.R -> add
tests/benchmark/check_environment.R -> add
tests/benchmark/README.md          -> add
tests/benchmark/QUICKSTART_CN.md   -> add
tests/benchmark/EXPERIMENT_PLAN.md -> add
tests/benchmark/slurm/*.sh         -> add

tests/testthat/test-benchmark-interface.R -> add
results/benchmark/*/.gitkeep       -> add
```

Do not replace the existing package dataset:

```text
data/primates.rda
```

Before the ray-finned fish analysis and simulation generation, add:

```text
data/c_tree.rds
data/c_traits.rds
```

`c_traits.rds` must contain `Species`, `ag`, `care`, and `spawning`.
