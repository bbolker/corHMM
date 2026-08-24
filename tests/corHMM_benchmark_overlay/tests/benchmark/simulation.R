sample_subtree <- function(tree, traits, ntaxa) {
  if (ntaxa > ape::Ntip(tree)) {
    stop("`ntaxa` exceeds the number of tree tips.")
  }

  names(traits)[1] <- "Species"
  traits$Species <- as.character(traits$Species)

  keep_tips <- sample(tree$tip.label, ntaxa)
  sub_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep_tips))

  original_total_branch_length <- sum(sub_tree$edge.length)
  sub_tree$edge.length <- sub_tree$edge.length / original_total_branch_length
  sub_tree <- ape::reorder.phylo(sub_tree, "pruningwise")

  sub_traits <- traits[
    match(sub_tree$tip.label, traits$Species),
    ,
    drop = FALSE
  ]
  rownames(sub_traits) <- NULL

  stopifnot(identical(sub_tree$tip.label, sub_traits$Species))

  list(
    tree = sub_tree,
    traits = sub_traits,
    original_total_branch_length = original_total_branch_length,
    scale_factor = 1 / original_total_branch_length,
    scaled_total_branch_length = sum(sub_tree$edge.length)
  )
}

realfun <- function(raw_tree, ntrait, nstate = 2L) {
  sample_trait <- function() {
    sample(0:(nstate - 1L), size = ape::Ntip(raw_tree), replace = TRUE)
  }

  repeat {
    d <- as.data.frame(replicate(ntrait, sample_trait(), simplify = FALSE))
    names(d) <- paste0("trait", seq_len(ntrait))
    if (nrow(unique(d)) == nstate^ntrait) break
  }

  d <- cbind(Species = raw_tree$tip.label, d)
  rownames(d) <- NULL

  list(tree = raw_tree, data = d)
}
