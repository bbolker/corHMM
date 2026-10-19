## new

* we want everything (as much as possible) to run through a *single* code path
    * it's OK if the formula interface only works with RTMB
	* priors should ideally work with both base-corHMM and RTMB-corHMM, but it's OK if they only work with RTMB
	* `nlminb` should *only be called once*: don't duplicate any code you don't have to
	* HMC code can only work with RTMB backend (but should be able to work with or without formula interface)
* extend the formula interface to multi-state (>2) traits. Currently binary-only: `corHMM()` hard-codes `translate(nstate = 2L)` and rejects values other than 0/1/NA; `translate()` uses one `nstate` for all traits and codes states as `0:(nstate-1)`; transitions are grouped by trait and gain/loss direction only (so 0->1, 1->2, 0->2 would share a rate; `symm()` assumes one gain + one loss); predictor traits enter `model.matrix()` as numeric, not factors. Needs per-trait `nstate`, a decision on multistate parameterization (ordered/adjacent-only up/down formulas? per-pair formulas?), and factor predictors. Then drop/modify the "formula traits must be binary" test in `test-formula-missing.R`.
* irreversible model + `root.p = "yang"` + `use_RTMB = TRUE` fails in `nlminb` with "NA/NaN gradient evaluation" (e.g. `primates[[2]][, 1:2]`, `getStateMat4Dat()` rate matrix with `rm[2,1] <- 0`). Probably the stationary distribution putting zero weight on the non-absorbing state (log(0) in the gradient)? Not yet checked whether this predates the merge of 03a9323.


## old (porting formulas etc to corHMM)

* port `formula_interface.R` code to corHMM
* test `sym()` syntax (there are still `browser()` calls in the code!)
* ? allow `..rateclass` (or whatever) in formulas
* reverse-mapping functions? (i.e. going from Q parameters to beta, we can do this by inverting the X matrices with `solve()`) 
   * i.e. `Q = X %*% beta` → `beta = Xinv %*% Q`
* priors
* check benchmarking results: why are additive models doing better (BMB) ?

## old 

### checking vs corHMM

* `rate.cat` > 1 ?  (I think this should probably work, but I don't know if it actually does. Are there any examples in the documentation that use rate.cat > 1?)
  * go to one of the examples in misc/ratecat_ex_files.txt and try it with use_RTMB = TRUE and see if it breaks! and if not, if it gives (nearly) the same answer ...
  
* when are differences in loglik more likely/bigger? when fit is worse, is it resolved by multi-start? is there a pattern to which parameters differ?

* (small tweak) set things up so we can turn RTMB on and off with a global user option?

  change the default argument of `use_RTMB` to `use_RTMB = getOption("corHMM.use_RTMB", FALSE)`
  
  then if you set `option(corHMM.use_RTMB = TRUE)` then corHMM will use RTMB by default. This can be handy for testing, because you set the option and then run all the tests


* OpenMP parallelization? (Need to think about interaction with multi-start etc.) (see misc/openmp_ex.R ...)  [seems slightly lower priority]

* add code to corHMM that stores the log-likelihoods from *all* of the random starts? (This would be useful e.g. for figuring out whether multiple modes are really being found ...) [maybe BMB will do this?] 
    * if we had this we could take whatever summary we wanted of the distribution of log-likelihoods when we were summarizing results (min-max? standard dev? )

* `use.fog` ? test!

## miscellaneous

* `quiet` option
* return devfun etc etc.? (see main branch of BMB fork)
* priors? (tricky, they have to be built into the likelihood)
* construct contrasts?
* try HMC/`tmbstan` ... [tried it, it's slow without priors!]
