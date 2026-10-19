# Profiling `dev.corhmm` (2026-10-13)

## Setup

- Trees: Yule trees (`phytools::pbtree`, height 1).
- Data: a binary trait simulated under Mk.
- Model: ARD, `rate.cat` = 1/2/3 (2/4/6 states), `root.p = "yang"`, `order.test = FALSE`.
- BLAS limited to one thread.
- Arguments built with `corProcessData()` / `rate.cat.set.corHMM.JDB()`, phylogeny reordered `"pruningwise"`, as in `corHMM()`.
- All variants agree with `dev.corhmm` to within 1e-6 (relative).

Scripts are in `misc/profiling/`: `setup.R`, `bench.R`, `rtmb.R`, `prune.cpp` (this section). The RTMB eigendecomposition work (below) uses `bench_pkg.R`, `bench_head.R`, `check_pruning.R`, and prototypes `proto_*.R`, `test_*.R`.

## Timings (median ms per likelihood evaluation)

| states | tips | A: current | B: R loop + eigen | C: C++, `expmat` per edge | D: C++, one eigen | RTMB fn | RTMB gr |
|---|---|---|---|---|---|---|---|
| 2 | 200 | 62 | 13 | 3.6 | 0.9 | | |
| 2 | 1000 | 303 | 80 | 16 | 2.0 | | |
| 2 | 5000 | 2363 | 1125 | 79 | 9.3 | | |
| 4 | 200 | 56 | 11 | 4.5 | 1.0 | 4.8 | 20 |
| 4 | 1000 | 323 | 84 | 19 | 2.5 | 22 | 102 |
| 4 | 5000 | 2477 | 1244 | 106 | 11 | 107 | 473 |
| 6 | 1000 | 362 | 85 | 38 | 3.1 | | |
| 6 | 5000 | 2402 | 1092 | 178 | 14 | | |

RTMB was timed with different parameters on every call, because `fn` returns a cached value when the parameters haven't changed. Building the tape (`MakeADFun`) takes about 2 s at 1000 tips and 15 s at 5000 tips.

## Where the time goes (current code, 1000 tips, 4 states, `Rprof`)

- **`expm()`: 64% of total time, but only 17% self time.** Most of the cost is wrapper overhead (`match.arg`, `stopifnot`, `formals`, `Sys.getenv`) paid on each of the ~2000 per-edge calls (`R/corHMM.R:741`). There is one more call for the sanity check at `R/corHMM.R:696`.
- **`which(phy$edge[,1] == focal)`: 14%.** It scans all edges for every internal node, so the cost grows as O(n²). This explains why A and B scale worse than linearly from 1000 to 5000 tips.
- The rest is ordinary interpreter overhead in the node loop.

Since a single character is a single site pattern, there is little parallelism available inside one likelihood evaluation. BEAGLE-style parallelism across sites doesn't apply. The fixes below target per-call overhead and gradient cost.

## Fix 1 (cheap, pure R): precomputed topology + one eigendecomposition

**(a) Precompute the topology once, outside the objective function.** `phy` and `liks` don't change during optimization.

```r
anc <- unique(phy$edge[,1])
child_rows <- split(seq_len(nrow(phy$edge)), phy$edge[,1])[as.character(anc)]
## pass child_rows (or store it in the closure) and use
##   desRows <- child_rows[[i]]
## instead of which(phy$edge[,1] == focal)
```

**(b) Compute P(t) from a single eigendecomposition of Q per evaluation**, instead of one `expm()` call per edge:

```r
ev  <- eigen(Q)                       # complex in general (ARD Q is not symmetric)
V   <- ev$vectors; Vi <- solve(V); lam <- ev$values
## for edge r with child c:
##   P(t) %*% x = V %*% (exp(lam * t) * (Vi %*% x))
v <- v * Re(V %*% (exp(lam * phy$edge.length[r]) * (Vi %*% liks[c, ])))
```

You never need to form P(t) explicitly, because only P(t)·x is used.

**Caveat:** Q can be defective, i.e. have repeated eigenvalues without a full set of eigenvectors. This can happen with zero rates, dropped states, or rate classes with identical rates. In that case `V` is singular or ill-conditioned. Guard with something like `kappa(V) > 1e8` (or `tryCatch(solve(V))`) and fall back to `expm()` per edge for that evaluation.

You can likewise replace the `expm(Q)` check at line 696 with a check on the eigenvalues, e.g. real parts ≤ 0 within tolerance.

Expected speedup: about 4× at 200–1000 tips (column B), and more at large n once the O(n²) `which` is removed. Column B still includes it.

Going further, compiled code (Rcpp/RcppArmadillo) with a single eigendecomposition (column D) is 100–250× faster than the current code. This prototype covers only the basic model: it doesn't yet handle fixed nodes, tip fog, Lewis correction, or root options other than "yang".

## Fix 2 (RTMB path): faster gradients

RTMB's `fn` already runs at compiled speed (similar to column C). The gradient costs about 5× `fn`, which is typical for reverse-mode AD over a taped loop with an `expm` at every edge. The tape is also expensive to build at large n. Options, roughly in order of effort:

1. **Apply Fix 1 inside `prune_nll`** (`R/RTMB_devcorhmm.R`). Replace the per-edge `Matrix::expm(Q * t)` (line 457) with one eigendecomposition. This only works if RTMB can differentiate through the eigendecomposition, which it may not for non-symmetric Q. Check whether RTMB provides AD-able `eigen` for general matrices. If it doesn't, skip to option 2.

2. **Write the pruning as an RTMB atomic function with an analytic gradient.** The pre-order ("outside") pass of Ji et al. (2020, *MBE*) gives the derivative with respect to every branch length in one extra tree traversal. The derivative with respect to the Q entries follows from the derivative of the matrix exponential, i.e. the Fréchet derivative, which is cheap in the eigenbasis when Q is diagonalizable. The total gradient cost is then about 2–3× `fn`, versus the current ~5×, and no tape is built over the tree.

3. **Move the likelihood to a TMB C++ template.** The loop is then compiled and taped in C++, which avoids R-level taping cost and probably shortens tape construction substantially at 5000 tips. This means more code to maintain, and it moves away from the pure-R RTMB design.

## Implemented in the RTMB backend (2026-10-16)

`mkdev.corhmm_rtmb()` (`R/RTMB_devcorhmm.R`) now has `pruning = c("auto", "eigen", "expm")`.

- **Topology precomputed once:** each node's child edges (for `"expm"`) and the edges grouped by tree level (`prune_levels()`, for `"eigen"`). A node's level is 1 + the maximum level of its children (tips are level 0).
- **`"eigen"`:** one AD eigendecomposition of Q per evaluation (`RTMB::eigen(Q, symmetric = FALSE)`, complex in general), then P(t)·x = V exp(Λt) V⁻¹ x for all edges at a level in one set of matrix operations. Real and imaginary parts are carried separately so the loop uses real AD arithmetic. Done edge by edge instead of level by level, taping was ~7× *slower* than the old code (34 s at 5000 tips), because of the many small complex operations taped from R.
- **`"expm"`:** the previous per-edge `Matrix::expm`.
- **`"auto"` (default):** `"eigen"`, except where Q has repeated eigenvalues. There, derivatives through the eigendecomposition are wrong *and finite* (no NaN to flag them); away from exact ties, the gradient error is at roundoff level until eigenvalues agree to ~1 part in 10⁶ (`test_gap.R`). So:
  - eigenvalues repeated at two arbitrary parameter vectors (structural, e.g. ER with > 2 states) → `"expm"` throughout;
  - otherwise the returned object is built on the `"expm"` tape, which `tmbstan`, `sdreport`, `$he` and `$env$f` use directly. Only `$fn`, `$gr` and `$report` use the eigen tape, falling back to `"expm"` when the relative eigenvalue gap at the current parameters is < 1e-6 or κ(V) > 1e8 (e.g. equal starting values), or when the eigen tape throws an error or returns non-finite values. `$fn` keeps `last.par.best`/`value.best` up to date. Both tapes are built, so set-up costs about the old tape time plus the eigen tape time.

  An earlier version (commit `c6a524b`) gave `tmbstan` the eigen tape. An HMC chain (`benchmark/compare_prior.R`, prior sd 0.5, chain 4) then aborted early in warmup: Stan stepped to extreme log-rates where Q was numerically defective, and RTMB's complex `solve()` (`solve_complex_atomic`, which calls base `solve()` from R during tape evaluation) threw "system is computationally singular". An R error inside an atomic aborts the chain instead of giving a non-finite log density that Stan would reject. RTMB's eigen reverse pass also drops the terms coupling equal eigenvalues (`Delta <- Conj(Diff)/(Diff*Conj(Diff) + 1e-300)`, `diag(Delta) <- 0`), which explains the finite-but-wrong gradients at repeated eigenvalues. A minimal RTMB-only reproducer is drafted as an upstream issue. `pruning = "eigen"` still gives `tmbstan` the fast tape, at the user's risk.

Gotchas inside the package namespace (RTMB not attached): `eigen()`/`solve()` must be called as `RTMB::eigen()`/`RTMB::solve()`, and `matrix()` and `rowSums()` drop the AD class (use `t()` and `v %*% rep(1, n)` instead).

### Timings (`bench_pkg.R`, `bench_head.R`; same setup as above; ms per evaluation)

"old" is the code before this change (commit `f321b24`); "new" is `pruning = "auto"`. Negative log-likelihoods agree within 5e-13; gradients agree with `"expm"` within 1e-11 (`check_pruning.R`: distinct and equal rates, 6 states, fixed nodes, structural ER). Test suite passes.

| states | tips | tape (s): old → new | fn: old → new | gr: old → new |
|---|---|---|---|---|
| 2 | 200 | 0.19 → 0.56* | 1.8 → 0.95 (1.9×) | 8.5 → 2.6 (3.3×) |
| 2 | 1000 | 1.0 → 0.12 | 9.1 → 1.8 (5.2×) | 43 → 3.4 (13×) |
| 2 | 5000 | 5.4 → 0.40 | 46 → 4.4 (10×) | 210 → 9.7 (22×) |
| 4 | 200 | 0.14 → 0.07 | 2.5 → 1.7 (1.5×) | 11 → 2.4 (4.6×) |
| 4 | 1000 | 1.0 → 0.21 | 12 → 3.2 (3.9×) | 57 → 6.5 (8.7×) |
| 4 | 5000 | 7.2 → 1.1 | 61 → 11 (5.6×) | 328 → 27 (12×) |
| 6 | 200 | 0.16 → 0.10 | 3.9 → 2.1 (1.9×) | 18 → 3.7 (4.8×) |
| 6 | 1000 | 0.97 → 0.38 | 20 → 5.2 (3.8×) | 84 → 12 (7.2×) |
| 6 | 5000 | 9.6 → 2.2 | 95 → 22 (4.3×) | 407 → 53 (7.6×) |

\*First run in the session; includes one-off warm-up.

A gradient now costs ~2–2.5× a function evaluation (was ~5×), which covers most of what Fix 2 aimed for. Not yet checked: very unbalanced (caterpillar-like) trees, where nearly every edge is its own level and the level grouping helps less; timings of complete `corHMM()` fits.

## Parallelism

At the scale of a single likelihood evaluation, parallelism doesn't pay. Run in parallel across optimizer restarts, candidate models (`corHMMDredge`), or trees instead. When doing so, limit each worker to one BLAS/OpenMP thread to avoid oversubscription:

```r
RhpcBLASctl::blas_set_num_threads(1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1)
```
