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





//' Generate random variable from a Discrete Weibull distribution of type 1
//' 
//' @param q The first parameter, \eqn{0 < q < 1}.
//' @param beta The second parameter, \eqn{\beta > 0}.
//' @return A random variable.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
double C_rdweibull(double& q, double& beta) {
  
  double u = R::runif(0, 1);
  return ceil(pow(log(1 - u) / log(q), 1 / beta));
  
}





//' Simulate parametric semi-Markov chain
//' 
//' @param seed A \code{unsigned int} that is the seed one wishes to use.
//' @param nsim A vector of integers specifying the length of the sequences.
//' @param init Vector of initial distribution of length s.
//' @param ptrans Matrix of transition probabilities of the embedded Markov chain.
//' @param distr A matrix giving the conditional sojourn time distributions for
//'   each transition from the current state \eqn{i} to the next state \eqn{j}.
//'   The values of the elements of the matrix \code{distr} are:
//'   \itemize{
//'     \item 0: if the sojourn time distribution is a uniform distribution;
//'     \item 1: if the sojourn time distribution is a geometric distribution;
//'     \item 2: if the sojourn time distribution is a poisson distribution;
//'     \item 3: if the sojourn time distribution is a discrete weibull 
//'       distribution;
//'     \item 4: if the sojourn time distribution is a negative binomial 
//'       distribution.
//'   }
//' @param param1 A matrix giving the first parameters of the sojourn time 
//'   distributions.
//' @param param2 A matrix giving the second parameters (if necessary, 
//'   otherwise the value is NA) of the sojourn time distributions.
//' @param censBeg A logical value indicating whether or not sequences are 
//'   censored at the beginning.
//' @param censEnd A logical value indicating whether or not sequences are 
//'   censored at the end.
//'   
//' @return A list of vectors representing the simulated semi-Markov chains.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List simulateParam(unsigned int& seed, arma::Col<arma::uword>& nsim, arma::vec& init, 
                 arma::mat& ptrans, arma::Mat<int>& distr, arma::mat& param1, arma::mat& param2, 
                 bool censBeg = false, bool censEnd = false) {
  
  setSeed(seed);
  
  arma::uword nbseq = nsim.n_elem;
  
  arma::uvec states = arma::regspace<arma::uvec>(0, init.n_elem - 1);
  
  arma::vec prob(init.n_elem, arma::fill::zeros);
  int distrib;
  double parameter1;
  double parameter2;
  double k = 0;
  
  arma::uword i;
  arma::uword t;
  
  List sequences(nbseq);
  
  for (arma::uword m = 0; m < nbseq; m++) {
    
    arma::vec J(nsim(m), arma::fill::zeros);
    arma::vec T(nsim(m), arma::fill::zeros);
    
    // Initial state ...
    J(0) = RcppArmadillo::sample(states, 1, false, init)(0);
    
    i = 0;
    t = 1;
    
    // ... and the other states
    while (t <= nsim(m)) {
      
      // Sample the following state ...
      prob = ptrans.row(J(i)).t();
      J(i + 1) = RcppArmadillo::sample(states, 1, false, prob)(0);
      
      // ... and sample the sojourn time according to the correct distribution
      distrib = distr(J(i), J(i + 1));
      
      parameter1 = param1(J(i), J(i + 1));
      parameter2 = param2(J(i), J(i + 1));
      
      switch (distrib) {
      case 0: // unif
        k = RcppArmadillo::sample(arma::regspace(1, (int)parameter1), 1, false)(0);
        break;
      case 1: // geom
        k = R::rgeom(parameter1) + 1;
        break;
      case 2: // pois
        k = R::rpois(parameter1) + 1;
        break;
      case 3: // dweibull
        k = C_rdweibull(parameter1, parameter2);
        break;
      case 4: // nbinom
        k = R::rnbinom(parameter1, parameter2) + 1;
        break;
      }
      
      T(i) = t + k;
      t = T(i);
      i += 1;
      
    }
    
    arma::vec sequence = getChain(J.subvec(0, i - 1), T.subvec(0, i - 1)) + 1;
    
    // Censoring sequences
    if (censBeg && censEnd) {// Censoring at the beginning and at the end
      
      arma::uword n = nsim(m);
      arma::uword l = t - n;
      arma::uword Nl = (arma::uword)std::floor(l / 2.0);
      
      if (Nl == 0) {
        sequences[m] = sequence.subvec(Nl, (t - 1 - Nl - 1));
      } else {
        sequences[m] = sequence.subvec(Nl - 1, (t - 1 - Nl - 1));
      }
      
    } else if (!censBeg && censEnd) {// Censoring at the end
      
      sequences[m] = sequence.subvec(0, nsim(m) - 1);
      
    } else if (censBeg && !censEnd) {
      
      arma::uword n = nsim(m);
      arma::uword l = t - n;
      
      sequences[m] = sequence.subvec(l - 1, (t - 1 - 1));
      
    } else {
      
      sequences[m] = sequence;
      
    }
    
  }
  
  return sequences;
  
}





//' Simulate nonparametric semi-Markov chain
//' 
//' @param seed A \code{unsigned int} that is the seed one wishes to use.
//' @param nsim A vector of integers specifying the length of the sequences.
//' @param init Vector of initial distribution of length s.
//' @param ptrans Matrix of transition probabilities of the embedded Markov chain.
//' @param distr A cube giving the conditional sojourn time distributions for
//'   each transition from the current state \eqn{i} to the next state \eqn{j}.
//' @param censBeg A logical value indicating whether or not sequences are 
//'   censored at the beginning.
//' @param censEnd A logical value indicating whether or not sequences are 
//'   censored at the end.
//'   
//' @return A list of vectors representing the simulated semi-Markov chains.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List simulateNonParam(unsigned int& seed, arma::Col<arma::uword>& nsim, arma::vec& init, 
                         arma::mat& ptrans, arma::cube& distr, bool censBeg = false, bool censEnd = false) {
  
  setSeed(seed);
  
  arma::uword nbseq = nsim.n_elem;
  arma::uword kmax = distr.n_slices;
  
  arma::uvec states = arma::regspace<arma::uvec>(0, init.n_elem - 1);
  arma::uvec range1kmax = arma::regspace<arma::uvec>(1, kmax);
  
  arma::vec prob(init.n_elem, arma::fill::zeros);
  double k = 0;
  
  arma::uword i;
  arma::uword t;
  
  List sequences(nbseq);
  
  for (arma::uword m = 0; m < nbseq; m++) {
    
    arma::vec J(nsim(m), arma::fill::zeros);
    arma::vec T(nsim(m), arma::fill::zeros);
    
    // Initial state ...
    J(0) = RcppArmadillo::sample(states, 1, false, init)(0);
    
    i = 0;
    t = 1;
    
    // ... and the other states
    while (t <= nsim(m)) {
      
      // Sample the following state ...
      prob = ptrans.row(J(i)).t();
      J(i + 1) = RcppArmadillo::sample(states, 1, false, prob)(0);
      
      // ... and sample the sojourn time according to the correct distribution
      k = RcppArmadillo::sample(range1kmax, 1, false, distr.tube(J(i), J(i + 1)))(0);
      
      T(i) = t + k;
      t = T(i);
      i += 1;
      
    }
    
    arma::vec sequence = getChain(J.subvec(0, i - 1), T.subvec(0, i - 1)) + 1;
    
    // Censoring sequences
    if (censBeg && censEnd) {// Censoring at the beginning and at the end
      
      arma::uword n = nsim(m);
      arma::uword l = t - n;
      arma::uword Nl = (arma::uword)std::floor(l / 2.0);
      
      if (Nl == 0) {
        sequences[m] = sequence.subvec(Nl, (t - 1 - Nl - 1));
      } else {
        sequences[m] = sequence.subvec(Nl - 1, (t - 1 - Nl - 1));
      }
      
    } else if (!censBeg && censEnd) {// Censoring at the end
      
      sequences[m] = sequence.subvec(0, nsim(m) - 1);
      
    } else if (censBeg && !censEnd) {
      
      arma::uword n = nsim(m);
      arma::uword l = t - n;
      
      sequences[m] = sequence.subvec(l - 1, (t - 1 - 1));
      
    } else {
      
      sequences[m] = sequence;
      
    }
    
  }
  
  return sequences;
  
}





//' Discrete-time convolution product of \eqn{f} and \eqn{g} 
//'   (See definition 2.2 p. 20)
//' 
//' @param f A vector giving the values of the function \eqn{f} for each 
//'   \eqn{k \in \mathbb{N}}.
//' @param g A vector giving the values of the function \eqn{g} for each 
//'   \eqn{k \in \mathbb{N}}.
//'   
//' @return A vector giving the values of the discrete-time convolution of 
//'   \eqn{f} and \eqn{g} for each \eqn{k \in \mathbb{N}}.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::vec convolution(arma::vec& f, arma::vec& g) {
  
  arma::uword k = f.n_elem - 1;
  
  arma::vec conv(k + 1, arma::fill::zeros);
  
  for (arma::uword l = 0; l <= k; l++) {
    conv(l) = sum(reverse(f.rows(0, l)) % g.rows(0, l));
  }
  
  return conv;
  
}





//' Discrete-time matrix convolution product 
//'   (See definition 3.5 p. 48)
//' 
//' @param A A cube of dimension \eqn{(S, S, k + 1)}.
//' @param B A cube of dimension \eqn{(S, S, k + 1)}.
//'   
//' @return A cube of dimension \eqn{(S, S, k + 1)} giving the discrete-time 
//'   matrix convolution product for each \eqn{k \in \mathbb{N}}.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::cube matrixConvolution(arma::cube& A, arma::cube& B) {
  
  // ###########################################################
  // ###########################################################
  // A and B must be of dimension (S, S, k + 1) where:
  //   - S represents the cardinal of the state space E;
  //   - k represents the time horizon;
  // 
  // Return: Cube C which is the matrix convolution product 
  //  A * B for each m, m = 0,..., k
  // ###########################################################
  // ###########################################################
  
  arma::uword k = A.n_slices - 1;
  
  arma::cube C(arma::size(A), arma::fill::zeros);
  
  for (arma::uword m = 0; m <= k; m++) {
    for (arma::uword l = 0; l <= m; l++) {
      C.slice(m) += A.slice(m - l) * B.slice(l);
    }
  }
  
  return C;
}





//' Compute the variance of the estimator of the transition function P 
//'   (See equation (4.29), p.91)
//'   
//'   @details The decomposition of the variance is as follows: 
//'     \deqn{\sigma_{P}^{2}(i, j, k)) = \sum_{m = 1}^{s} \mu_{mm} \left\{ \underbrace{\sum_{r = 1}^{s} \underbrace{\left[\delta_{mj}\Psi_{ij} - \underbrace{(1 - H_{j}) * \psi_{im} * \psi_{rj}}_{\text{part11}} \right]^2}_{\text{part12}} * \ q_{mr}(k)}_{\text{part1}} - \left[ \underbrace{\delta_{mj} \psi_{ij} * H_{m}(k)}_{\text{part22}} - \underbrace{\sum_{r = 1}^{s} (1 - H_{j}) * \psi_{im} * \psi_{rj} * q_{mr}}_{\text{part21}} \right]^{2}(k) \right\}}
//'   
//' @return A cube of dimension \eqn{(S, S, k + 1)} giving the values of the 
//'   variances for each transition from \eqn{i} to \eqn{j} and each time 
//'   horizon \eqn{k \in \mathbb{N}}.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::cube varP(arma::vec& mu, arma::cube& q, arma::cube& psi, arma::cube& Psi, arma::cube& H) {
  
  arma::uword s = psi.n_rows;
  arma::uword k = psi.n_slices - 1;
  
  arma::vec convolpsi(k + 1, arma::fill::zeros);
  arma::vec part11(k + 1, arma::fill::zeros); // (1 - H_{j}) * \psi_{im} * \psi_{rj}
  arma::vec part12(k + 1, arma::fill::zeros); // \left[\delta_{mj}\Psi_{ij} - (1 - H_{j}) * \psi_{im} * \psi_{rj}\right]^2
  arma::vec part1(k + 1, arma::fill::zeros); // \sum_{r = 1}^{s} \left[\delta_{mj}\Psi_{ij} - (1 - H_{j}) * \psi_{im} * \psi_{rj}\right]^2 * q_{mr}(k)
  
  arma::vec part21(k + 1, arma::fill::zeros); // \sum_{r = 1}^{s} (1 - H_{j}) * \psi_{im} * \psi_{rj} * q_{mr}
  arma::vec part22(k + 1, arma::fill::zeros); // \delta_{mj} \psi_{ij} * H_{m}(k)
  
  arma::vec part1bis(k + 1, arma::fill::zeros);
  arma::vec part2bis(k + 1, arma::fill::zeros);
  
  
  arma::cube sigma2_i_j_k(size(psi), arma::fill::zeros);
  
  
  arma::vec bar_H_j(k + 1, arma::fill::zeros);
  arma::vec Psi_i_j(k + 1, arma::fill::zeros);
  arma::vec psi_i_j(k + 1, arma::fill::zeros);
  arma::vec psi_i_m(k + 1, arma::fill::zeros);
  arma::vec psi_r_j(k + 1, arma::fill::zeros);
  arma::vec q_m_r(k + 1, arma::fill::zeros);
  arma::vec H_m(k + 1, arma::fill::zeros);
  int delta_m_j;
  
  for (arma::uword i = 0; i < s; i++) {
    
    for (arma::uword j = 0; j < s; j++) {
      
      bar_H_j = 1 - H.tube(j, j);
      Psi_i_j = Psi.tube(i, j);
      psi_i_j = psi.tube(i, j);
      
      part1bis.zeros();
      part2bis.zeros();
      
      for (arma::uword m = 0; m < s; m++) {
        
        psi_i_m = psi.tube(i, m);
        H_m = H.tube(m, m);
        
        if (m == j) {
          delta_m_j = 1;
        } else {
          delta_m_j = 0;
        }
        
        part1.zeros();
        part21.zeros();
        
        for (arma::uword r = 0; r < s; r++) {
          
          psi_r_j = psi.tube(r, j);
          q_m_r = q.tube(m, r);
          
          convolpsi = convolution(psi_i_m, psi_r_j);
          
          part11 = convolution(bar_H_j, convolpsi);
          part12 = arma::square(delta_m_j * Psi_i_j - part11);
          part1 += convolution(part12, q_m_r);
          
          part21 += convolution(part11, q_m_r);
          
        }
        
        part22 = delta_m_j * convolution(psi_i_j, H_m);
        
        part1bis += mu(m) * part1;
        part2bis += mu(m) * arma::square(part22 - part21);
        
      }
      
      sigma2_i_j_k.tube(i, j) = part1bis - part2bis;
      
    }
  }
  
  return sigma2_i_j_k;
}
