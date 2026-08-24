# corHMM benchmark 快速使用说明

所有命令均从 corHMM 仓库根目录运行。

## 1. 文件放置

将压缩包内容直接解压到仓库根目录并保留相对路径：

```text
R/corHMM.R
R/RTMB_devcorhmm.R
tests/benchmark/*
tests/testthat/test-benchmark-interface.R
results/benchmark/*/.gitkeep
```

其中两个 `R/` 文件会替换现有源码。原有 `data/primates.rda` 不替换。

辐鳍鱼数据准备完成后，放置为：

```text
data/c_tree.rds
data/c_traits.rds
```

`c_traits.rds` 至少包含：

```text
Species
ag
care
spawning
```

三个性状必须是完整的 0/1 二元变量。

## 2. 本地环境变量

### Windows PowerShell

```powershell
$env:BAYESTRAITS_HOME = "E:/BayesTraitsV5.0.3-Win64"
$env:BAYESTRAITS_BIN = "BayesTraitsV5"
```

`BAYESTRAITS_BIN` 不写 `.exe`。当前 Windows 版 `bbolker/btw` 还要求同一目录存在 `BayesTraitsV3.exe` 兼容副本。

### Linux / MSI

```bash
export BAYESTRAITS_HOME="$HOME/software/BayesTraitsV5.0.3-Linux"
export BAYESTRAITS_BIN="BayesTraitsV5"
```

## 3. 安装依赖与检查源码

```bash
Rscript tests/benchmark/install_benchmark_packages.R
Rscript tests/benchmark/check_environment.R
```

随后在 R 中运行：

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

## 4. 本地运行 primates

### 先跑快速流程检查

```bash
Rscript tests/benchmark/run_local_primates.R smoke
Rscript tests/benchmark/summarize_real.R
```

该命令依次运行：

```text
original corHMM ML
corHMM-RTMB ML
corHMM formula ML
RTMB adaptive MH
RTMB-HMC
rphylo fixed-model MCMC
BayesTraits fixed-model MCMC through btw
```

结果位置：

```text
results/benchmark/primates/
results/benchmark/summary/primates_methods.csv
results/benchmark/summary/primates_rates.csv
results/benchmark/summary/primates_pairwise.csv
```

### 正式长链

```bash
Rscript tests/benchmark/run_local_primates.R full
Rscript tests/benchmark/summarize_real.R
```

### 分阶段运行

```bash
Rscript tests/benchmark/run_real.R primates ml full
Rscript tests/benchmark/run_real.R primates bayes full
```

### 单独运行 RJMCMC

```bash
Rscript tests/benchmark/run_real.R primates rj full
Rscript tests/benchmark/summarize_real.R
```

RJMCMC 不进入固定模型 sampler 的 RMSE/ESS 主比较。

## 5. 本地运行辐鳍鱼应用

确认两个文件存在：

```text
data/c_tree.rds
data/c_traits.rds
```

### 快速流程检查

```bash
Rscript tests/benchmark/run_local_rayfinned.R smoke
Rscript tests/benchmark/summarize_real.R
```

### 正式分析

```bash
Rscript tests/benchmark/run_local_rayfinned.R full
Rscript tests/benchmark/summarize_real.R
```

使用的模型为：

```r
list(
  ag ~ care * spawning,
  care ~ 1,
  spawning ~ 1
)
```

主要结果写入：

```text
results/benchmark/summary/rayfish_main_effect.csv
```

其余结果位置：

```text
results/benchmark/rayfish/
results/benchmark/summary/rayfish_methods.csv
results/benchmark/summary/rayfish_rates.csv
```

## 6. 本地运行 simulation 的一行任务

先生成模拟数据和 manifest：

```bash
Rscript tests/benchmark/generate_simulation_data.R
Rscript tests/benchmark/build_manifests.R smoke
```

测试第一个 ML 任务：

```bash
Rscript tests/benchmark/run_job.R \
  tests/benchmark/manifests/jobs_ml_smoke.csv \
  1
```

测试第一个 Bayesian 任务：

```bash
Rscript tests/benchmark/run_job.R \
  tests/benchmark/manifests/jobs_bayes_smoke.csv \
  1
```

汇总：

```bash
Rscript tests/benchmark/summarize_simulation.R
```

## 7. MSI 提交

```bash
Rscript tests/benchmark/generate_simulation_data.R
Rscript tests/benchmark/build_manifests.R smoke

bash tests/benchmark/slurm/submit_ml.sh smoke
bash tests/benchmark/slurm/submit_bayes.sh smoke
bash tests/benchmark/slurm/submit_rj.sh smoke
```

确认 smoke jobs 后：

```bash
Rscript tests/benchmark/build_manifests.R full

bash tests/benchmark/slurm/submit_ml.sh full
bash tests/benchmark/slurm/submit_bayes.sh full
bash tests/benchmark/slurm/submit_rj.sh full
```

汇总模型结果和资源使用：

```bash
Rscript tests/benchmark/summarize_simulation.R
Rscript tests/benchmark/collect_resources.R
```

资源 CSV 写入：

```text
results/benchmark/summary/slurm_resources.csv
```

## 8. 时间解释

本地 Windows 运行适合确认结果和完整流程。`rphylo`、BTW 和 RTMB-MH 的四条链在 Windows 上会顺序执行，而 HMC 可以由 rstan 并行执行。因此论文中的正式 wall-time 比较应使用 MSI 的统一四 CPU 配置，不应直接混用本地 Windows 时间。
