## checking vs corHMM

* `rate.cat` > 1 ?  (I think this should probably work, but I don't know if it actually does. Are there any examples in the documentation that use rate.cat > 1?)
  * go to one of the examples in misc/ratecat_ex_files.txt and try it with use_RTMB = TRUE and see if it breaks! and if not, if it gives (nearly) the same answer ...
  
* when are differences in loglik more likely/bigger? when fit is worse, is it resolved by multi-start? is there a pattern to which parameters differ?


* (small tweak) set things up so we can turn RTMB on and off with a global user option?

  change the default argument of `use_RTMB` to `use_RTMB = getOption("corHMM.use_RTMB", FALSE)`
  
  then if you set `option(corHMM.use_RTMB = TRUE)` then corHMM will use RTMB by default. This can be handy for testing, because you set the option and then run all the tests


* OpenMP parallelization? (Need to think about interaction with multi-start etc.) (see misc/openmp_ex.R ...)  [seems slightly lower priority]
* `use.fog` ? [seems hard]

## miscellaneous

* `quiet` option
* return devfun etc etc.? (see main branch of BMB fork)
* priors?
* construct contrasts?
* try HMC/`tmbstan` ...
