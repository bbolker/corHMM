packages <- c(
  "ape", "RTMB", "pkgload", "numDeriv", "posterior",
  "ramcmc", "rphylo", "rstan", "tmbstan", "btw"
)

status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
print(status)
if (!all(status)) stop("Install the missing benchmark packages.")

btw_description <- utils::packageDescription("btw")
if (!identical(btw_description$RemoteUsername, "bbolker")) {
  stop("The benchmark requires bbolker/btw from GitHub.")
}

bt_path <- Sys.getenv("BAYESTRAITS_HOME")
bt_bin <- Sys.getenv("BAYESTRAITS_BIN", unset = "BayesTraitsV5")
windows <- .Platform$OS.type == "windows"
executable <- file.path(bt_path, if (windows) paste0(bt_bin, ".exe") else bt_bin)

cat("R:", R.version.string, "\n")
cat("BayesTraits:", executable, "\n")
if (!file.exists(executable)) stop("BayesTraits executable was not found.")

pkgload::load_all(".", quiet = TRUE)
q <- translate(list(T1 ~ T2, T2 ~ T1))
stopifnot(q$n_par == 8L, identical(dim(q$rate_X), c(8L, 8L)))
cat("corHMM source and formula mapping loaded successfully.\n")
