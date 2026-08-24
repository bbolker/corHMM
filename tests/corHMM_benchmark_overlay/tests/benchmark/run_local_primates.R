args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args)) args[1] else "full"
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("tests/benchmark/run_real.R", "primates", "all", profile)
)
if (status != 0L) quit(status = status)
