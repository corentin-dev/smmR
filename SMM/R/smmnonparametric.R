#' Non-parametric semi-Markov model specification
#'
#' @description Creates a non-parametric model specification for a semi-Markov model.
#'
#' @details This function creates a semi-Markov model object in the 
#' non-parametric case, taking into account the type of sojourn time and the 
#' censoring described in references. The non-parametric specification concerns 
#' sojourn time distributions defined by the user.
#'
#' The difference between the Markov model and the semi-Markov model concerns 
#' the modelisation of the sojourn time. With a Markov chain, the sojourn time 
#' distribution is modeled by a Geometric distribution (in discrete time). 
#' With a semi-Markov chain, the sojourn time can be any arbitrary distribution.
#' 
#' We define :
#'  \itemize{
#'    \item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in states} of the embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Y_1 = i)}, \eqn{i \in 1, 2, \dots, s};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in states,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )}, f is specified by the argument "param" in the parametric case and by "distr" in the non-parametric case.
#'  }
#'
#' In this package we can choose differents types of sojourn time.
#' Four options are available for the sojourn times:
#' \itemize{
#'   \item depending on the present state and on the next state (\eqn{f_{ij}});
#'   \item depending only on the present state (\eqn{f_{i}});
#'   \item depending only on the next state (\eqn{f_{j}});
#'   \item depending neither on the current, nor on the next state (\eqn{f}).
#' }
#' 
#' Let define kmax the maximum length of the sojourn times.
#' If  `type.sojourn = "fij"`, `distr` is an array of size SxSxKmax.
#' If `type.sojourn = "fi"` or `"fj"`, `distr` must be a matrix of size SxKmax.
#' If `type.sojourn = "f"`, `distr` must be a vector of length kmax.
#' 
#' If the sequence is censored at the beginning and at the end, `cens.beg` 
#' must be equal to `TRUE` and `cens.end` must be equal to `TRUE` too. 
#' All the sequences must be censored in the same way.
#'
#' @param states Vector of state space of length s.
#' @param init Vector of initial distribution of length s.
#' @param ptrans Matrix of transition probabilities of the embedded Markov chain 
#'   \eqn{J=(J_m)_{m}} of size sxs.
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param distr
#'   \itemize{
#'     \item Array of size SxSxKmax if `type.sojourn = "fij"`;
#'     \item Matrix of size SxKmax if `type.sojourn = "fi"` or `"fj"`;
#'     \item Vector of length kmax if the `type.sojourn = "f"`.
#'   }
#'   kmax is the maximum length of the sojourn times.
#' @param cens.beg Optional. A logical value indicating whether or not 
#'   sequences are censored at the beginning.
#' @param cens.end Optional. A logical value indicating whether or not 
#'   sequences are censored at the end.
#' @return Returns an object of class [smmnonparametric][smmnonparametric].
#' 
#' 
#' @seealso [simulate], [fitsemimarkovmodel], [smmparametric]
#' @export
#'
#' @examples 
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' 
#' # Creation of the initial distribution
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' 
#' # Creation of transition matrix
#' pij <- matrix(c(0, 0.2, 0.5, 0.3, 
#'                 0.2, 0, 0.3, 0.5, 
#'                 0.3, 0.5, 0, 0.2, 
#'                 0.4, 0.2, 0.4, 0), 
#'               ncol = s, byrow = TRUE)
#' 
#' # Creation of a matrix corresponding to the 
#' # conditional sojourn time distributions
#' kmax <- 6
#' nparam.matrix <- matrix(c(0.2, 0.1, 0.3, 0.2, 
#'                           0.2, 0, 0.4, 0.2, 
#'                           0.1, 0, 0.2, 0.1, 
#'                           0.5, 0.3, 0.15, 0.05, 
#'                           0, 0, 0.3, 0.2, 
#'                           0.1, 0.2, 0.2, 0), 
#'                         nrow = s, ncol = kmax, byrow = TRUE)
#' 
#' smm2 <- smmnonparametric(states = states, init = vect.init, ptrans = pij, 
#'                          type.sojourn = "fj", distr = nparam.matrix)
#' 
#' smm2
smmnonparametric <- function(states, init, ptrans, type.sojourn = c("fij", "fi", "fj", "f"), 
                             distr, cens.beg = FALSE, cens.end = FALSE) {
  
  #############################
  # Checking parameter states
  #############################
  
  s <- length(states)
  
  if (!(is.vector(states) && (length(unique(states)) == s))) {
    stop("The state space states is not a vector of unique elements")
  }
  
  #############################
  # Checking parameter init
  #############################
  
  if (!(is.vector(init) && (length(init) == s))) {
    stop("init is not a vector of length s")
  }
  
  if (!(all(init >= 0) && all(init <= 1))) {
    stop("Probabilities in init must be between [0, 1]")
  }
  
  if (!(sum(init) == 1)) {
    stop("The sum of init is not equal to one")
  }
  
  #############################
  # Checking parameter ptrans
  #############################
  
  if (!is.matrix(ptrans)) {
    stop("ptrans is not a matrix")
  }
  
  if (!(all(ptrans >= 0) && all(ptrans <= 1))) {
    stop("Probabilities in ptrans must be between [0, 1]")
  }
  
  if (!((dim(ptrans)[1] == s) && (dim(ptrans)[2] == s))) {
    stop("The size of the matrix ptrans must be equal to sxs")
  }
  
  if (!all(apply(ptrans, 1, sum) == 1)) {
    stop("ptrans is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }
  
  if (!all(diag(ptrans) == 0)) {
    stop("All the diagonal elements of ptrans must be equal to 0 since transitions to the same state are not allowed")
  }
  
  #############################
  # Checking parameter type.sojourn
  #############################
  
  type.sojourn <- match.arg(type.sojourn)
  
  #############################
  # Checking parameter distr
  #############################
  
  if (type.sojourn == "fij" && !(is.array(distr) && !is.matrix(distr))) {
    stop("distr must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !is.matrix(distr)) {
    stop("distr must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !is.vector(distr)) {
    stop("distr must be a vector of length kmax since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !(dim(distr)[1] == s && dim(distr)[2] == s)) {
    stop("distr must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !(dim(distr)[1] == s)) {
    stop("distr must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (!all(distr >= 0) | !all(distr <= 1)) {
    stop("Probabilities in distr must be between [0, 1]")
  }
  
  if (type.sojourn == "fij") {
    temp <- apply(distr, c(1, 2), sum)
    indexdiag <- seq(1, s * s, by = s + 1)
    
    if (!(all(diag(temp == 0)) && all(temp[-indexdiag] == 1))) {
      stop("distr is not a stochastic matrix")
    }
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !(all(apply(distr, 1, sum) == 1))) {
    stop("distr is not a stochastic matrix")
  }
  
  if (type.sojourn == "f" && !(sum(distr) == 1)) {
    stop("distr is not a stochastic matrix")
  }
  
  #############################
  # Checking parameter cens.beg and cens.end
  #############################
  
  if (!(is.logical(cens.beg) && is.logical(cens.end))) {
    stop("cens.beg and cens.end must be TRUE or FALSE")
  }
  
  if (type.sojourn == "fij") {
    kmax <- dim(distr)[3]
  } else if (type.sojourn %in% c("fi", "fj")) {
    kmax <- dim(distr)[2]
  } else {
    kmax <- length(distr)
  }
  
  ans <-
    list(
      states = states,
      s = s,
      kmax = kmax,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = distr,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(ans) <- c("smm", "smmnonparametric")
  
  return(ans)
}

# Function to check if an object is of class smmnonparametric
is.smmnonparametric <- function(x) {
  inherits(x, "smmnonparametric")
}

# Method to get the semi-Markov kernel q
.get.q.smmnonparametric <- function(x, kmax = x$kmax) {
  
  q <- array(data = 0, dim = c(x$s, x$s, kmax))
  
  if (x$type.sojourn == "fij") {
    q[, , 1:x$kmax] <- array(x$ptrans, c(x$s, x$s, x$kmax)) * x$distr
  } else if (x$type.sojourn == "fi") {
    q[, , 1:x$kmax] <- array(x$ptrans, c(x$s, x$s, x$kmax)) * aperm(array(x$distr, c(x$s, x$kmax, x$s)), c(1, 3, 2))
  } else if (x$type.sojourn == "fj") {
    q[, , 1:x$kmax] <- array(x$ptrans, c(x$s, x$s, x$kmax)) * aperm(array(x$distr, c(x$s, x$kmax, x$s)), c(3, 1, 2))
  } else if (x$type.sojourn == "f") {
    q[, , 1:x$kmax] <- array(x$ptrans, c(x$s, x$s, x$kmax)) * aperm(array(x$distr, c(x$kmax, x$s, x$s)), c(2, 3, 1))
  }
  
  return(q)
  
}

#' Loglikelihood
#'
#' @description Computation of the loglikelihood for a semi-Markov model
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @param states Vector of state space (of length s).
#' @return A vector giving the value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik.smmnonparametric <- function(x, sequences, states) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!is.list(sequences)) {
    stop("The parameter sequences should be a list")
  }
  
  if (!all(unique(unlist(sequences)) %in% states)) {
    stop("Some states in the list of observed sequences sequences are not in the state space states")
  }
  
  s <- length(states)
  
  #############################
  # Checking smm parameter
  #############################
  
  if ((x$s != s)) {
    stop("The size of the matrix ptrans must be equal to sxs with s = length(states)")
  }
  
  if (!all.equal(states, x$states)) {
    stop("The state space of the estimated SMM smm is different from the given state states")
  }
  
  
  sequences <- processes(sequences = sequences, states = states)
  kmax <- sequences$kmax
  
  if (!(kmax == x$kmax)) {
    stop("kmax of the given sequences is different from the kmax of the estimated SMM model")  
  }
  
  type.sojourn <- x$type.sojourn
  cens.end <- x$cens.end
  
  #############################
  # Let's compute the loglikelihood
  #############################
  
  init <- x$init # Initial distribution
  Nstarti <- sequences$counting$Nstarti
  maskNstarti <- Nstarti != 0 & init != 0
  
  if (!cens.end) {# No censoring
    
    pij <- x$ptrans # Transition matrix
    Nij <- sequences$counting$Nij
    maskNij <- Nij != 0 & pij != 0
    
    # Contribution of the initial distribution and 
    # the transition matrix to the loglikelihood
    loglik <- sum(Nstarti[maskNstarti] * log(init[maskNstarti])) +
      sum(Nij[maskNij] * log(pij[maskNij]))
    
    
    # Contribution of the sojourn time distribution
    if (type.sojourn == "fij") {
      
      Nijk <- sequences$counting$Nijk
      maskNijk <- Nijk != 0 & x$distr != 0
      
      loglik <- loglik + sum(Nijk[maskNijk] * log(x$distr[maskNijk]))
      
    } else if (type.sojourn == "fi") {
      
      Nik <- sequences$counting$Nik
      maskNik <- Nik != 0 & x$distr != 0
      
      loglik <- loglik + sum(Nik[maskNik] * log(x$distr[maskNik]))
      
    } else if (type.sojourn == "fj") {
      
      Njk <- sequences$counting$Njk
      maskNjk <- Njk != 0 & x$distr != 0
      
      loglik <- loglik + sum(Njk[maskNjk] * log(x$distr[maskNjk]))
      
    } else {
      
      Nk <- sequences$counting$Nk
      maskNk <- Nk != 0 & x$distr != 0
      
      loglik <- loglik + sum(Nk[maskNk] * log(x$distr[maskNk]))
      
    }
    
  } else {# Censoring
    
    s <- sequences$s
    Y <- sequences$Y
    U <- sequences$U
    
    
    # Computation of Niujv (couple Markov chain (Y, U))
    Niujv <- .getCountingNiujv(Y, U, s, kmax)
    Niu <- apply(Niujv, c(1, 2), sum)
    
    phat <- Niujv / array(Niu, c(s, kmax, s, kmax))
    phat[is.na(phat)] <- 0
    
    # Computation of q
    q <- .computeKernelNonParamEndcensoring(phat)
    
    piujv <-
      apply(
        X = phat,
        MARGIN = c(1, 2, 3),
        FUN = function(x)
          ifelse(length(x[x != 0]) > 0, x[x != 0], 0)
      )
    
    Niub <-
      apply(
        X = Niujv,
        MARGIN = c(1, 2, 3),
        FUN = function(x)
          ifelse(length(x[x != 0]) > 0, x[x != 0], 0)
      )
    
    maskNiub <- Niub != 0 & piujv != 0
    
    loglik <- sum(Nstarti * log(init)) +
      sum(Niub[maskNiub] * log(piujv[maskNiub]))
    
  }
  
  return(loglik)
}

# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.smmnonparametric <- function(x) {
  
  s <- x$s
  kmax <- x$kmax
  type.sojourn <- x$type.sojourn
  
  if (type.sojourn == "fij") {
    kpar <- s * (s - 2) * (kmax - 1)
  } else if (type.sojourn == "fi") {
    kpar <- s * (kmax - 1)
  } else if (type.sojourn == "fj") {
    kpar <- s * (kmax - 1)
  } else {
    kpar <- kmax - 1
  }
  
  return(kpar)
}

#' Akaike Information Criterion (AIC)
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @param states Vector of state space (of length s).
#' @return A numeric value giving the value of the AIC.
#' 
#' 
#' @export
#'
aic.smmnonparametric <- function(x, sequences, states) {
  
  loglik <- loglik(x, sequences, states)
  sequences <- processes(sequences = sequences, states = states)
  
  s <- x$s
  kmax <- sequences$kmax
  
  kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * kpar
  
  return(aic)
  
}

#' Bayesian Information Criterion (BIC)
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC criterion must be computed.
#' @param states Vector of state space (of length s).
#' @return A numeric value giving the value of the BIC.
#' 
#' 
#' @export
#'
bic.smmnonparametric <- function(x, sequences, states) {
  
  loglik <- loglik(x, sequences, states)
  sequences <- processes(sequences = sequences, states = states)
  
  s <- x$s
  kmax <- sequences$kmax
  
  kpar <- .getKpar(x)
  
  n <- sum(sapply(sequences$Ym, length))
  
  bic <- -2 * loglik + log(n) * kpar
  
  return(bic)
  
}