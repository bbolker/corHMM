source(file.path("tests", "benchmark", "config.R"))
source(file.path("tests", "benchmark", "simulation.R"))

if (!file.exists(file.path("data", "c_tree.rds")) ||
    !file.exists(file.path("data", "c_traits.rds"))) {
  stop("Place c_tree.rds and c_traits.rds in data/ before generating simulations.")
}

base_tree <- readRDS(file.path("data", "c_tree.rds"))
base_traits <- readRDS(file.path("data", "c_traits.rds"))
names(base_traits)[1] <- "Species"
base_traits$Species <- as.character(base_traits$Species)

rows <- list()
k <- 1L

for (ntrait in benchmark_design$ntraits) {
  for (ntips in benchmark_design$tree_sizes) {
    for (seed in benchmark_design$dataset_seeds) {
      set.seed(seed)

      sub <- sample_subtree(base_tree, base_traits, ntips)
      sim <- realfun(sub$tree, ntrait = ntrait)

      dataset_id <- sprintf("sim_t%d_n%03d_s%03d", ntrait, ntips, seed)
      out_dir <- file.path(
        "tests", "benchmark", "data", "sim",
        paste0("t", ntrait),
        sprintf("n%03d", ntips),
        sprintf("seed_%03d", seed)
      )
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      tree_path <- benchmark_path(out_dir, "tree.rds")
      data_path <- benchmark_path(out_dir, "traits.rds")
      saveRDS(sim$tree, tree_path)
      saveRDS(sim$data, data_path)

      state_counts <- table(do.call(paste, c(sim$data[-1], sep = "|")))

      rows[[k]] <- data.frame(
        dataset_id = dataset_id,
        dataset_type = "simulation",
        ntraits = ntrait,
        ntips = ntips,
        dataset_seed = seed,
        tree_path = tree_path,
        data_path = data_path,
        original_total_branch_length = sub$original_total_branch_length,
        tree_scale_factor = sub$scale_factor,
        scaled_total_branch_length = sub$scaled_total_branch_length,
        state_counts = paste(names(state_counts), as.integer(state_counts), sep = ":", collapse = ";"),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}

manifest <- do.call(rbind, rows)
out_file <- file.path("tests", "benchmark", "manifests", "datasets_simulation.csv")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, out_file, row.names = FALSE)
cat("Wrote", nrow(manifest), "datasets to", out_file, "\n")
