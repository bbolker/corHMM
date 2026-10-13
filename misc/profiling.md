# Profiling `dev.corhmm` (2026-10-13)

## Setup

- Trees: Yule trees (`phytools::pbtree`, height 1).
- Data: a binary trait simulated under Mk.
- Model: ARD, `rate.cat` = 1/2/3 (2/4/6 states), `root.p = "yang"`, `order.test = FALSE`.
- BLAS limited to one thread.
- Arguments built with `corProcessData()` / `rate.cat.set.corHMM.JDB()`, phylogeny reordered `"pruningwise"`, as in `corHMM()`.
- All variants agree with `dev.corhmm` to within 1e-6 (relative).

Scripts (not in the repo): `setup.R`, `bench.R`, `rtmb.R`, `prune.cpp`.

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

## Parallelism

At the scale of a single likelihood evaluation, parallelism doesn't pay. Run in parallel across optimizer restarts, candidate models (`corHMMDredge`), or trees instead. When doing so, limit each worker to one BLAS/OpenMP thread to avoid oversubscription:

```r
RhpcBLASctl::blas_set_num_threads(1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1)
```
