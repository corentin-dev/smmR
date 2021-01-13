// [[Rcpp::depends(RcppArmadillo)]]
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

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





//' Return the semi-Markov chain given the processes J and T
//' 
//' @param J A vector giving the successively visited states.
//' @param T A vector giving the successive time points.
//' @return A vector giving the reconstructed semi-Markov chain.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::vec getChain(arma::vec J, arma::vec T) {
  
  arma::uword n = T.n_elem;
  arma::vec chain(T(n - 1) - 1, arma::fill::zeros);
  
  arma::uword k;
  arma::uword t = 1;
  arma::uword i = 0;
  
  while (i < n) {
    
    k = T(i) - t;
    chain.subvec(t - 1, t - 1 + k - 1).fill(J(i));
    t = T(i);
    i += 1;
    
  }
  
  return chain;
}





//' Return the processes Y, J, T, L and U given the sequences of states
//' 
//' @param sequences A list of sequences of states.
//' @return A list giving the processes.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List getProcesses(List sequences) {
  
  int nbseq = sequences.size();
  
  List Ym(nbseq);
  List Tm(nbseq);
  List Jm(nbseq);
  List Lm(nbseq);
  List Um(nbseq);
  
  for (int l = 0; l < nbseq; l++) {
    
    arma::vec sequence = sequences[l];
    
    // Process Y - Semi-Markov chain
    arma::vec Y = sequences[l];
    Ym[l] = Y;
    
    // Process T - Successive time points when state changes
    arma::uvec idx = arma::find(arma::diff(Y) != 0) + 1;
    arma::uvec T(idx.n_elem + 1, arma::fill::zeros);
    T.subvec(0, idx.n_elem - 1) = idx + 1;
    T(idx.n_elem) = Y.n_elem + 1;
    Tm[l] = T;
    
    // Process J - Successively visited states (coded with numbers)
    arma::vec J;
    J = arma::join_cols(Y.subvec(0, 0), Y.elem(idx));
    Jm[l] = J;
    
    // Process L - Sojourn time
    arma::uvec L(T.n_elem, arma::fill::zeros);
    
    L = T - arma::join_cols(arma::uvec(1, arma::fill::ones), idx + 1);
    Lm[l] = L;
    
    // Process U - Backward recurrence time
    arma::uword n = L.n_elem;
    arma::uword k = 0;
    arma::uvec U(Y.n_elem, arma::fill::zeros);
    
    for (arma::uword i = 0; i < n; i++) {
      
      U.subvec(k, k + L(i) - 1) = arma::regspace<arma::uvec>(0, L(i) - 1);
      k += L(i);
      
    }
    Um[l] = U;
    
  }
  
  return List::create(Named("Ym") = Ym,
                      _["Tm"] = Tm,
                      _["Jm"] = Jm,
                      _["Lm"] = Lm,
                      _["Um"] = Um);
  
}





//' Give the values of the counting processes Nij, Ni, Nj...
//' 
//' @param Jm A list of sequences of states.
//' @param Lm A list of sojourn time processes.
//' @param s Cardinal of the state space S.
//' @param kmax Maximum length of the sojourn times.
//' @return A list giving the counting processes.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List getCountingProcesses(List& Jm, List& Lm, arma::uword& s, arma::uword& kmax) {
  
  int nbseq = Jm.size();
  
  arma::vec Nstarti(s, arma::fill::zeros);
  arma::cube Nijk(s, s, kmax, arma::fill::zeros);
  arma::cube Nbijk(s, s, kmax, arma::fill::zeros);
  arma::mat Neik(s, kmax, arma::fill::zeros);
  
  arma::uword n;
  
  for (int l = 0; l < nbseq; l++) {
    
    arma::vec J = Jm[l];
    arma::vec L = Lm[l];
    n = J.n_elem;
    
    Nstarti(J(0)) += 1;
    Nbijk(J(0), J(1), L(0) - 1) += 1;
    Neik(J(n - 1), L(n - 1) - 1) += 1;
    
    for (arma::uword i = 1; i < n; i++) {
      
      Nijk(J(i - 1), J(i), L(i - 1) - 1) += 1;
      
    }
    
  }
  
  arma::mat Nij = sum(Nijk, 2);
  
  arma::mat Nbik = sum(Nbijk, 1);
  arma::mat Nbjk = sum(Nbijk, 0);
  arma::vec Nbk = sum(Nbjk, 0).t();
  
  arma::vec Nek = sum(Neik, 0).t();
  
  arma::mat Nebik = Neik + Nbik;
  arma::vec Nebk = Nek + Nbk;
  
  arma::vec Ni = sum(Nij, 1);
  arma::vec Nj = sum(Nij, 0).t();
  
  double N = accu(Ni);
  
  arma::mat Nik = sum(Nijk, 1);
  arma::mat Njk = sum(Nijk, 0);
  arma::vec Nk = sum(Nik, 0).t();
  
  return List::create(Named("Nijk") = Nijk,
                      _["Nij"] = Nij,
                      _["Nik"] = Nik,
                      _["Njk"] = Njk,
                      _["Ni"] = Ni,
                      _["Nj"] = Nj,
                      _["Nk"] = Nk,
                      _["N"] = N,
                      _["Nbijk"] = Nbijk,
                      _["Nbik"] = Nbik,
                      _["Nbjk"] = Nbjk,
                      _["Nbk"] = Nbk,
                      _["Neik"] = Neik,
                      _["Nek"] = Nek,
                      _["Nebik"] = Nebik,
                      _["Nebk"] = Nebk,
                      _["Nstarti"] = Nstarti);
  
}





//' Give the values of the counting processes (cf. article 
//' Exact MLE and asymptotic properties for nonparametric semi-Markov models)
//' 
//' @param Ym A list of sequences of states.
//' @param Um A list of backward recurrence time processes.
//' @param s Cardinal of the state space S.
//' @param kmax Maximum length of the sojourn times.
//' @return A list giving the counting processes.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::cube getCountingNiuj(List& Ym, List& Um, arma::uword& s, arma::uword& kmax) {
  
  int nbseq = Ym.size();
  
  arma::cube Niuj(s, kmax, s, arma::fill::zeros);
  
  arma::uword n;
  
  for (int l = 0; l < nbseq; l++) {
    
    arma::vec Y = Ym[l];
    arma::vec U = Um[l];
    n = Y.n_elem;
    
    for (arma::uword i = 1; i < n; i++) {
      
      Niuj(Y(i - 1), U(i - 1), Y(i)) += 1;
      
    }
    
  }
  
  return Niuj;
  
}





//' Give the values of the kernel q (cf. Proposition 3.1 article 
//' Exact MLE and asymptotic properties for nonparametric semi-Markov models)
//' 
//' @param p A cube representing the transition matrix of the chain (Y, U).
//' @return A cube giving the MLE of the kernel q.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::cube computeKernelNonParamEndcensoring(arma::cube& p) {
  
  arma::uword s = p.n_rows;
  arma::uword kmax = p.n_cols;
  
  arma::cube q(s, s, kmax + 1, arma::fill::zeros);
  
  for (arma::uword i = 0; i < s; i++) {
    
    arma::rowvec temp = p.subcube(arma::span(i), arma::span::all, arma::span(i));
    
    for (arma::uword j = 0; j < s; j++) {
      
      if (i != j) {
        
        // Case k = 1
        q(i, j, 1) = p(i, 0, j);
        
        // Case k >= 2
        for (arma::uword k = 2; k <= kmax; k++) {
          
          q(i, j, k) = p(i, k - 1, j) * prod(temp.subvec(0, k - 2));
          
        } 
        
      }
    }
  }
  
  return q;
  
}
