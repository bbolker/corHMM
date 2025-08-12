#' @param n (vector of) number of states (if n is a scalar, all traits have the same number of states)
#' @param k number of traits (only used if n is a scalar)
#' @examples
#' setup_Q_template(n = 3, k = 2)
#' setup_Q_template(n = 2, k = 3)
Q_template <- function(n=3, k= 1, set_indices = TRUE) {
  if (length(n) == 1) {
    n <- rep(n, k)
  }
  all_states <- do.call(expand.grid, lapply(n, \(x) 0:(x-1)))
  ## dimnms to match corHMM standard
  dimnms <- apply(all_states, 1, \(x) paste(x, collapse = "|"))
  ns <- prod(n)
  m <- matrix(0, ns, ns)
  for (i in 1:ns) {
    ## exactly one state changes ...
    for (j in 1:ns) {
      m[i,j] <- as.numeric(sum(all_states[i,] != all_states[j,])== 1)
    }
  }
  dimnames(m) <- list(dimnms, dimnms)
  if (set_indices) {
    m[m!=0] <- seq_len(sum(m==1))
  }
  return(m)
}


## convert from enumerated traits (0-7) to 3 binary digits
## https://stackoverflow.com/questions/6614283/converting-decimal-to-binary-in-r
to_bin <- function(x, n = 2) {
   intToBits(x) |> rev() |> as.integer() |> tail(n)
}

setname <- function(x) {cbind(nm = rownames(x), x)}
##' @param nstate number of states per trait (currently limited to 2)
##' @param ntrait number of traits
##' @param ntaxa number of taxa/phylogenty tips
##' @param seed random-number seed
##' @param meanrate mean transition rate
##' @examples
##' simfun(ntaxa = 8)
simfun <- function(nstate = 2, ntrait = 2, ntaxa = 200, seed = NULL,
                   meanrate = 1) {
  if (ntrait!=2) stop("oops, simSeq to trait matrix not implemented for ntrait>1")
  if (!is.null(seed)) set.seed(seed)
  require("ape")
  require("phangorn")
  phy <- ape::rtree(ntaxa)
  phy <- reorder(phy, "pruningwise")

  Q <- Q_template(n=nstate, k= ntrait)
  nrates <- sum(Q != 0)
  Q[Q!=0] <- rexp(nrates, rate = 1/meanrate)

  s <- phangorn::simSeq(phy, l = 1, Q = Q,
                        type = "USER", levels = seq(nstate^ntrait),
                        rate = 1)
  traitMatrix <- sapply(s, function(x) to_bin(x-1)) |>
    t() |>
    as.data.frame() |>
    setname()

  list(phy = phy, data = traitMatrix)
}

