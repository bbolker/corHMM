// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// edges in pruningwise order (1-based node ids); L has rows = nodes
// method 0: expmat per edge; method 1: one complex eigendecomposition
// [[Rcpp::export]]
Rcpp::List prune_cpp(const imat& edge, const vec& el, mat L, const mat& Q, int method) {
  int nb = L.n_rows, k = L.n_cols, ne = edge.n_rows;
  vec comp(nb, fill::zeros);
  std::vector<bool> init(nb, false), done(nb, false);
  cx_vec lam; cx_mat V, Vi;
  if (method == 1) { eig_gen(lam, V, Q); Vi = inv(V); }
  for (int e = 0; e < ne; e++) {
    int par = edge(e,0)-1, ch = edge(e,1)-1;
    rowvec x = L.row(ch);
    if (init[ch] && !done[ch]) { double s = accu(x); comp[ch] = s; x /= s; L.row(ch) = x; done[ch] = true; }
    vec y;
    if (method == 0) y = expmat(Q * el[e]) * x.t();
    else y = real(V * (exp(lam * el[e]) % (Vi * conv_to<cx_vec>::from(x.t()))));
    if (!init[par]) { L.row(par).ones(); init[par] = true; }
    L.row(par) %= y.t();
  }
  int root = edge(ne-1,0)-1;
  double s = accu(L.row(root)); comp[root] = s; L.row(root) /= s;
  return Rcpp::List::create(Rcpp::Named("comp")=comp, Rcpp::Named("L")=L);
}
