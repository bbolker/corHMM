cran <- c(
  "ape", "Matrix", "expm", "RTMB", "pkgload", "numDeriv",
  "posterior", "ramcmc", "rphylo", "remotes", "rstan", "tmbstan"
)

missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

btw_ok <- requireNamespace("btw", quietly = TRUE) &&
  identical(utils::packageDescription("btw")$RemoteUsername, "bbolker")

if (!btw_ok) {
  remotes::install_github("bbolker/btw", upgrade = "never", force = TRUE)
}

cat("Benchmark packages installed.\n")
