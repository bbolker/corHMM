# corHMM benchmark workflow

This directory contains the simulation, local real-data, post-processing, and
Slurm scripts for the first benchmark round.

## 1. Copy the files into the repository

Copy the bundle contents into the root of the corHMM repository while
preserving paths.

Two package source files are replaced:

```text
R/corHMM.R
R/RTMB_devcorhmm.R
```

The benchmark scripts are added under:

```text
tests/benchmark/
```

A Chinese quick-start guide is provided at:

```text
tests/benchmark/QUICKSTART_CN.md
```

Generated simulation inputs are written under:

```text
tests/benchmark/data/sim/
```

Results are written under:

```text
results/benchmark/
```

## 2. Required real-data files

The primates data are loaded from the existing package dataset:

```text
data/primates.rda
```

For the ray-finned fish application, add exactly these files:

```text
data/c_tree.rds
data/c_traits.rds
```

`c_traits.rds` must contain the columns:

```text
Species
ag
care
spawning
```

The first column is treated as the species-name column. All three traits must
be complete binary variables coded 0/1.

## 3. Install and verify dependencies

From the repository root:

```r
Rscript tests/benchmark/install_benchmark_packages.R
Rscript tests/benchmark/check_environment.R
```

The BayesTraits wrapper must be the GitHub package `bbolker/btw`, not an
unrelated package with the same name.

Set these environment variables before running methods that use BayesTraits.
The executable name must not include `.exe`:

```text
BAYESTRAITS_HOME=<directory containing BayesTraitsV5>
BAYESTRAITS_BIN=BayesTraitsV5
```

On Windows the same directory must also contain the `BayesTraitsV3.exe`
compatibility copy required by the current `bbolker/btw` source.

## 4. Validate the modified package source

From the repository root:

```r
parse(file = "R/corHMM.R")
parse(file = "R/RTMB_devcorhmm.R")
devtools::load_all(reset = TRUE)

testthat::test_file(
  "tests/testthat/test-benchmark-interface.R",
  reporter = "summary"
)

testthat::test_file(
  "tests/testthat/testnew.R",
  reporter = "summary"
)
```

The revised `corHMM()` accepts explicit `prior.mean` and `prior.sd` arguments.
The formula interface stores `rate_X`, and the rate-scale prior is applied to
`rate_X %*% beta`. The RTMB backend does not use `REPORT()`.

## 5. Local primates analysis

### Fast workflow check

```bash
Rscript tests/benchmark/run_local_primates.R smoke
Rscript tests/benchmark/summarize_real.R
```

This runs:

```text
original corHMM ML
corHMM-RTMB ML
corHMM formula ML
RTMB adaptive MH
RTMB-HMC
rphylo fixed-model MCMC
BayesTraits fixed-model MCMC through btw
```

Results are written to:

```text
results/benchmark/primates/
results/benchmark/summary/primates_methods.csv
results/benchmark/summary/primates_rates.csv
```

### Full run

```bash
Rscript tests/benchmark/run_local_primates.R full
Rscript tests/benchmark/summarize_real.R
```

### Run only one stage

```bash
Rscript tests/benchmark/run_real.R primates ml full
Rscript tests/benchmark/run_real.R primates bayes full
```

### Separate reversible-jump analysis

```bash
Rscript tests/benchmark/run_real.R primates rj full
```

RJMCMC is stored separately and is not included in the fixed-model sampler
comparison.

## 6. Local ray-finned fish application

After adding `data/c_tree.rds` and `data/c_traits.rds`:

### Fast workflow check

```bash
Rscript tests/benchmark/run_local_rayfinned.R smoke
Rscript tests/benchmark/summarize_real.R
```

### Full run

```bash
Rscript tests/benchmark/run_local_rayfinned.R full
Rscript tests/benchmark/summarize_real.R
```

This fits:

```r
list(
  ag ~ care * spawning,
  care ~ 1,
  spawning ~ 1
)
```

with formula RTMB maximum likelihood and formula RTMB-HMC. The reported care
effect averages the log-rate contrast over the two spawning states:

```text
beta_care + 0.5 * beta_care:spawning
```

and exponentiates it to obtain the accessory-gland gain rate ratio.

Outputs are written to:

```text
results/benchmark/rayfish/
results/benchmark/summary/rayfish_methods.csv
results/benchmark/summary/rayfish_rates.csv
results/benchmark/summary/rayfish_main_effect.csv
```

## 7. Generate simulation datasets

The default design is defined in `tests/benchmark/config.R`:

```text
tree sizes: 50, 100, 150, 234
trait dimensions: 2 and 3
seeds: 900:999
```

Each subtree is scaled so that total branch length equals 1. Simulated binary
traits are regenerated until every joint state occurs at least once.

Generate all datasets:

```bash
Rscript tests/benchmark/generate_simulation_data.R
```

This creates:

```text
tests/benchmark/data/sim/
tests/benchmark/manifests/datasets_simulation.csv
```

Build smoke manifests:

```bash
Rscript tests/benchmark/build_manifests.R smoke
```

Build full manifests:

```bash
Rscript tests/benchmark/build_manifests.R full
```

The Bayesian and RJ manifests include only two-trait datasets because
`rphylo::run_dep_model()` and the BayesTraits Discrete dependent model are
pairwise binary-trait models. Three-trait simulations remain in the ML
benchmark.

## 8. Run one manifest row locally

This is useful before submitting an array:

```bash
Rscript tests/benchmark/run_job.R \
  tests/benchmark/manifests/jobs_ml_smoke.csv \
  1
```

For a Bayesian row:

```bash
Rscript tests/benchmark/run_job.R \
  tests/benchmark/manifests/jobs_bayes_smoke.csv \
  1
```

## 9. MSI Slurm submission

The scripts assume the repository is located at:

```text
$HOME/corHMM
```

and BayesTraits is located at:

```text
$HOME/software/BayesTraitsV5.0.3-Linux/BayesTraitsV5
```

The Slurm scripts load:

```text
r-rich/4.6.0_msi1.2
```

Submit smoke arrays first:

```bash
bash tests/benchmark/slurm/submit_ml.sh smoke
bash tests/benchmark/slurm/submit_bayes.sh smoke
bash tests/benchmark/slurm/submit_rj.sh smoke
```

After checking the output and error logs, submit the full arrays:

```bash
bash tests/benchmark/slurm/submit_ml.sh full
bash tests/benchmark/slurm/submit_bayes.sh full
bash tests/benchmark/slurm/submit_rj.sh full
```

Each Bayesian array task receives four CPUs and runs four chains within one
job. Every BayesTraits chain uses an independent temporary directory, so the
fixed file names generated by `btw` do not collide.

The submission wrappers split manifests into chunks of at most 1,000 array
tasks. The row offset is passed to `run_array.sh`, so manifests larger than the
cluster's usual single-array limit are still processed without changing the
manifest.

Resource information is written by `/usr/bin/time -v` to:

```text
logs/resource_<array-job>_<task>.txt
```

Convert the GNU time logs and their manifest-row metadata to one CSV with:

```bash
Rscript tests/benchmark/collect_resources.R
```

The default output is:

```text
results/benchmark/summary/slurm_resources.csv
```

Slurm accounting can additionally be queried with:

```bash
sacct -j <JOB_ID> \
  --format=JobID,State,Elapsed,TotalCPU,MaxRSS,ExitCode
```

## 10. Post-process simulation results

```bash
Rscript tests/benchmark/summarize_simulation.R
```

The main CSV files are written to:

```text
results/benchmark/summary/simulation_ml_fits.csv
results/benchmark/summary/simulation_ml_rates.csv
results/benchmark/summary/simulation_ml_pairwise.csv
results/benchmark/summary/simulation_ml_winners.csv
results/benchmark/summary/simulation_bayes_rates.csv
results/benchmark/summary/simulation_bayes_methods.csv
results/benchmark/summary/simulation_bayes_pairwise.csv
results/benchmark/summary/simulation_rj_methods.csv
results/benchmark/summary/simulation_rj_model_counts.csv
```

The NLL winner calculation treats absolute NLL differences no larger than
`1e-6` as ties. RMSE always uses named transitions and therefore works for both
8-rate two-trait models and 24-rate three-trait models.

## 11. Prior and tree scaling used in this first round

All fixed-model Bayesian methods use the current specification from
`analyses.md`:

```text
lower = log(1)
upper = log(100 * ntips)
mean = (lower + upper) / 2
sd = (upper - lower) / 3
```

These values are calculated per dataset and recorded in each job manifest.
They are independent of the wide numerical bounds used by the optimizer.

All trees are scaled to total branch length 1 before fitting. The original
branch length, scale factor, and scaled branch length are recorded for each
simulation dataset.

## 12. Local and HPC timing

Local Windows runs are intended to verify the complete inferential workflow.
The external four-chain methods run sequentially on Windows, whereas rstan can
parallelize HMC chains. Use the MSI four-CPU jobs for the wall-time values
reported in the paper.

## 13. Saved-object policy

Simulation jobs save standardized estimates, draws, likelihood/target output,
diagnostics, timing, warnings, and metadata. They do not save complete
`stanfit`, rphylo, or BTW fit objects.

Real-data jobs save both standardized results and complete raw fit objects.
Failed jobs still write an RDS file containing the job specification and error
message.
