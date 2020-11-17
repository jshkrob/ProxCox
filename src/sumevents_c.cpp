#include <RcppArmadillo.h>
using namespace Rcpp;

//' helper function written in Rcpp/RcppArmadilllo
//' @param i current index
//' @param X feature matrix
//' @param Ftime matrix of values of the basis for certain time points
//' @param theta current parameters of partial Cox likelihood
//' @param n dimension
//' @param p number of variables
//' @param eventtimes indexes of individuals that had event of interest

// [[Rcpp::export]]
List sumevents_c(int i, arma::mat X, arma::mat Ftime, arma::mat theta, int n, int p, arma::vec eventtimes) {
    arma::mat xj = X.rows((i-1), (X.n_rows - 1));
    arma::uvec idx = arma::find(eventtimes == i);
    arma::mat Fi = Ftime.rows(idx);
    arma::mat rr = exp(xj * theta * Fi.t());
    double hi = 1 / arma::accu(rr);
    arma::mat xmean = hi * (xj.t()) * rr;
    List L = List::create(hi, ((-1 * xmean) * (Fi)));
    return L;
}
