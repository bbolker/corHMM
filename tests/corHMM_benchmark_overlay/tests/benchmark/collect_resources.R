args <- commandArgs(trailingOnly = TRUE)
logs_dir <- if (length(args) >= 1L) args[1] else "logs"
out_file <- if (length(args) >= 2L) {
  args[2]
} else {
  file.path("results", "benchmark", "summary", "slurm_resources.csv")
}

parse_meta <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- x[grepl("=", x, fixed = TRUE)]
  key <- sub("=.*$", "", x)
  value <- sub("^[^=]*=", "", x)
  setNames(as.list(value), key)
}

extract_value <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(pattern, "", hit[1]))
}

elapsed_seconds <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  fields <- as.numeric(strsplit(x, ":", fixed = TRUE)[[1]])
  if (anyNA(fields)) return(NA_real_)
  if (length(fields) == 3L) {
    fields[1] * 3600 + fields[2] * 60 + fields[3]
  } else if (length(fields) == 2L) {
    fields[1] * 60 + fields[2]
  } else {
    fields[1]
  }
}

meta_files <- list.files(
  logs_dir,
  pattern = "^resource_.*\\.meta$",
  full.names = TRUE
)

if (!length(meta_files)) {
  stop("No resource metadata files were found in ", logs_dir, ".")
}

rows <- lapply(meta_files, function(meta_file) {
  meta <- parse_meta(meta_file)
  time_file <- sub("\\.meta$", ".txt", meta_file)
  lines <- if (file.exists(time_file)) readLines(time_file, warn = FALSE) else character()

  elapsed <- extract_value(
    lines,
    "^\\s*Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\):\\s*"
  )

  data.frame(
    array_job_id = if (!is.null(meta$array_job_id)) meta$array_job_id else NA_character_,
    array_task_id = if (!is.null(meta$array_task_id)) meta$array_task_id else NA_character_,
    slurm_job_id = if (!is.null(meta$slurm_job_id)) meta$slurm_job_id else NA_character_,
    manifest = if (!is.null(meta$manifest)) meta$manifest else NA_character_,
    manifest_row = if (!is.null(meta$row)) as.integer(meta$row) else NA_integer_,
    host = if (!is.null(meta$host)) meta$host else NA_character_,
    user_seconds = as.numeric(extract_value(lines, "^\\s*User time \\(seconds\\):\\s*")),
    system_seconds = as.numeric(extract_value(lines, "^\\s*System time \\(seconds\\):\\s*")),
    elapsed = elapsed,
    elapsed_seconds = elapsed_seconds(elapsed),
    cpu_percent = extract_value(lines, "^\\s*Percent of CPU this job got:\\s*"),
    max_rss_kb = as.numeric(extract_value(
      lines,
      "^\\s*Maximum resident set size \\(kbytes\\):\\s*"
    )),
    file_system_inputs = as.numeric(extract_value(lines, "^\\s*File system inputs:\\s*")),
    file_system_outputs = as.numeric(extract_value(lines, "^\\s*File system outputs:\\s*")),
    exit_status = as.integer(extract_value(lines, "^\\s*Exit status:\\s*")),
    resource_file = time_file,
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, out_file, row.names = FALSE)
cat("Wrote", nrow(out), "resource records to", out_file, "\n")
