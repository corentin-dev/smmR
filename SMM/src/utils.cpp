// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR

using namespace Rcpp;

//' Set the RNG Seed from within Rcpp
//' 
//' @param seed A \code{unsigned int} that is the seed one wishes to use.
//' @return A set RNG scope.
//' @examples
//' set.seed(10)
//' x <- rnorm(5, 0, 1)
//' setSeed(10)
//' y <- rnorm(5, 0, 1)
//' all.equal(x, y, check.attributes = FALSE)
//' 
//' @export
//' 
// [[Rcpp::export]]
void setSeed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function setSeedR = base_env["set.seed"];
  setSeedR(seed);
}
