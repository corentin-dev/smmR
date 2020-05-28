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
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in E} of the embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Y_1 = i)}, \eqn{i \in 1, 2, \dots, S};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in E,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )}, f is specified by the argument "param" in the parametric case and by "laws" in the non-parametric case.
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
#' Let define Kmax the maximum length of the sojourn times.
#' If  `type.sojourn = "fij"`, `laws` is an array of size SxSxKmax.
#' If `type.sojourn = "fi"` or `"fj"`, `laws` must be a matrix of size SxKmax.
#' If `type.sojourn = "f"`, `laws` must be a vector of length Kmax.
#' 
#' If the sequence is censored at the beginning and at the end, `cens.beg` 
#' must be equal to `TRUE` and `cens.end` must be equal to `TRUE` too. 
#' All the sequences must be censored in the same way.
#'
#' @param E Vector of state space of length S.
#' @param init Vector of initial distribution of length S.
#' @param ptrans Matrix of transition probabilities of the embedded Markov chain 
#'   \eqn{J=(J_m)_{m}} of size SxS.
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param laws
#'   \itemize{
#'     \item Array of size SxSxKmax if `type.sojourn = "fij"`;
#'     \item Matrix of size SxKmax if `type.sojourn = "fi"` or `"fj"`;
#'     \item Vector of length Kmax if the `type.sojourn = "f"`.
#'   }
#'   Kmax is the maximum length of the sojourn times.
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
#' E <- c("a", "c", "g", "t")
#' S <- length(E)
#' 
#' # Creation of the initial distribution
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' 
#' # Creation of transition matrix
#' pij <- matrix(c(0, 0.2, 0.5, 0.3, 
#'                 0.2, 0, 0.3, 0.5, 
#'                 0.3, 0.5, 0, 0.2, 
#'                 0.4, 0.2, 0.4, 0), 
#'               ncol = S, byrow = TRUE)
#' 
#' # Creation of a matrix corresponding to the 
#' # conditional sojourn time distributions
#' Kmax <- 6
#' nparam.matrix <- matrix(c(0.2, 0.1, 0.3, 0.2, 
#'                           0.2, 0, 0.4, 0.2, 
#'                           0.1, 0, 0.2, 0.1, 
#'                           0.5, 0.3, 0.15, 0.05, 
#'                           0, 0, 0.3, 0.2, 
#'                           0.1, 0.2, 0.2, 0), 
#'                         nrow = S, ncol = Kmax, byrow = TRUE)
#' 
#' smm2 <- smmnonparametric(E = E, init = vect.init, ptrans = pij, 
#'                          type.sojourn = "fj", laws = nparam.matrix)
#' 
#' smm2
smmnonparametric <- function(E, init, ptrans, type.sojourn = c("fij", "fi", "fj", "f"), 
                             laws, cens.beg = FALSE, cens.end = FALSE) {
  
  #############################
  # Checking parameter E
  #############################
  
  S <- length(E)
  
  if (!(is.vector(E) && (length(unique(E)) == S))) {
    stop("The state space E is not a vector of unique elements")
  }
  
  #############################
  # Checking parameter init
  #############################
  
  if (!(is.vector(init) && (length(init) == S))) {
    stop("init is not a vector of length S")
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
  
  if (!((dim(ptrans)[1] == S) && (dim(ptrans)[2] == S))) {
    stop("The size of the matrix ptrans must be equal to SxS")
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
  # Checking parameter laws
  #############################
  
  if (type.sojourn == "fij" && !(is.array(laws) && !is.matrix(laws))) {
    stop("laws must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !is.matrix(laws)) {
    stop("laws must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !is.vector(laws)) {
    stop("laws must be a vector of length Kmax since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !(dim(laws)[1] == S && dim(laws)[2] == S)) {
    stop("laws must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !(dim(laws)[1] == S)) {
    stop("laws must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (!all(laws >= 0) | !all(laws <= 1)) {
    stop("Probabilities in laws must be between [0, 1]")
  }
  
  if (type.sojourn == "fij") {
    temp <- apply(laws, c(1, 2), sum)
    indexdiag <- seq(1, S * S, by = S + 1)
    
    if (!(all(diag(temp == 0)) && all(temp[-indexdiag] == 1))) {
      stop("laws is not a stochastic matrix")
    }
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !(all(apply(laws, 1, sum) == 1))) {
    stop("laws is not a stochastic matrix")
  }
  
  if (type.sojourn == "f" && !(sum(laws) == 1)) {
    stop("laws is not a stochastic matrix")
  }
  
  #############################
  # Checking parameter cens.beg and cens.end
  #############################
  
  if (!(is.logical(cens.beg) && is.logical(cens.end))) {
    stop("cens.beg and cens.end must be TRUE or FALSE")
  }
  
  if (type.sojourn == "fij") {
    Kmax <- dim(laws)[3]
  } else if (type.sojourn %in% c("fi", "fj")) {
    Kmax <- dim(laws)[2]
  } else {
    Kmax <- length(laws)
  }
  
  ans <-
    list(
      E = E,
      S = S,
      Kmax = Kmax,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      laws = laws,
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
.get.q.smmnonparametric <- function(x, Kmax = x$Kmax) {
  
  q <- array(data = 0, dim = c(x$S, x$S, Kmax))
  
  if (x$type.sojourn == "fij") {
    q[, , 1:x$Kmax] <- array(x$ptrans, c(x$S, x$S, x$Kmax)) * x$laws
  } else if (x$type.sojourn == "fi") {
    q[, , 1:x$Kmax] <- array(x$ptrans, c(x$S, x$S, x$Kmax)) * aperm(array(x$laws, c(x$S, x$Kmax, x$S)), c(1, 3, 2))
  } else if (x$type.sojourn == "fj") {
    q[, , 1:x$Kmax] <- array(x$ptrans, c(x$S, x$S, x$Kmax)) * aperm(array(x$laws, c(x$S, x$Kmax, x$S)), c(3, 1, 2))
  } else if (x$type.sojourn == "f") {
    q[, , 1:x$Kmax] <- array(x$ptrans, c(x$S, x$S, x$Kmax)) * aperm(array(x$laws, c(x$Kmax, x$S, x$S)), c(2, 3, 1))
  }
  
  return(q)
  
}

#' Loglikelihood
#'
#' @description Computation of the loglikelihood for a semi-Markov model
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param seq A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik.smmnonparametric <- function(x, seq, E) {
  
  #############################
  # Checking parameters seq and E
  #############################
  
  if (!is.list(seq)) {
    stop("The parameter seq should be a list")
  }
  
  if (!all(unique(unlist(seq)) %in% E)) {
    stop("Some states in the list of observed sequences seq are not in the state space E")
  }
  
  S <- length(E)
  
  #############################
  # Checking smm parameter
  #############################
  
  if ((x$S != S)) {
    stop("The size of the matrix ptrans must be equal to SxS with S = length(E)")
  }
  
  if (!all.equal(E, x$E)) {
    stop("The state space of the estimated SMM smm is different from the given state E")
  }
  
  
  seq <- sequences(seq = seq, E = E)
  Kmax <- seq$Kmax
  
  if (!(Kmax == x$Kmax)) {
    stop("Kmax of the given sequences is different from the Kmax of the estimated SMM model")  
  }
  
  type.sojourn <- x$type.sojourn
  cens.end <- x$cens.end
  
  #############################
  # Let's compute the loglikelihood
  #############################
  
  init <- x$init # Initial distribution
  Nstarti <- seq$counting$Nstarti
  maskNstarti <- Nstarti != 0 & init != 0
  
  if (!cens.end) {# No censoring
    
    pij <- x$ptrans # Transition matrix
    Nij <- seq$counting$Nij
    maskNij <- Nij != 0 & pij != 0
    
    # Contribution of the initial distribution and 
    # the transition matrix to the loglikelihood
    loglik <- sum(Nstarti[maskNstarti] * log(init[maskNstarti])) +
      sum(Nij[maskNij] * log(pij[maskNij]))
    
    
    # Contribution of the sojourn time distribution
    if (type.sojourn == "fij") {
      
      Nijk <- seq$counting$Nijk
      maskNijk <- Nijk != 0 & x$laws != 0
      
      loglik <- loglik + sum(Nijk[maskNijk] * log(x$laws[maskNijk]))
      
    } else if (type.sojourn == "fi") {
      
      Nik <- seq$counting$Nik
      maskNik <- Nik != 0 & x$laws != 0
      
      loglik <- loglik + sum(Nik[maskNik] * log(x$laws[maskNik]))
      
    } else if (type.sojourn == "fj") {
      
      Njk <- seq$counting$Njk
      maskNjk <- Njk != 0 & x$laws != 0
      
      loglik <- loglik + sum(Njk[maskNjk] * log(x$laws[maskNjk]))
      
    } else {
      
      Nk <- seq$counting$Nk
      maskNk <- Nk != 0 & x$laws != 0
      
      loglik <- loglik + sum(Nk[maskNk] * log(x$laws[maskNk]))
      
    }
    
  } else {# Censoring
    
    S <- seq$S
    Y <- seq$Y
    U <- seq$U
    
    
    # Computation of Niujv (couple Markov chain (Y, U))
    Niujv <- .getCountingNiujv(Y, U, S, Kmax)
    Niu <- apply(Niujv, c(1, 2), sum)
    
    phat <- Niujv / array(Niu, c(S, Kmax, S, Kmax))
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
  
  S <- x$S
  Kmax <- x$Kmax
  type.sojourn <- x$type.sojourn
  
  if (type.sojourn == "fij") {
    Kpar <- S * (S - 2) * (Kmax - 1)
  } else if (type.sojourn == "fi") {
    Kpar <- S * (Kmax - 1)
  } else if (type.sojourn == "fj") {
    Kpar <- S * (Kmax - 1)
  } else {
    Kpar <- Kmax - 1
  }
  
  return(Kpar)
}

#' Akaike Information Criterion (AIC)
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param seq A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A numeric value giving the value of the AIC.
#' 
#' 
#' @export
#'
aic.smmnonparametric <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  seq <- sequences(seq = seq, E = E)
  
  S <- x$S
  Kmax <- seq$Kmax
  
  Kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * Kpar
  
  return(aic)
  
}

#' Bayesian Information Criterion (BIC)
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param seq A list of vectors representing the sequences for which the 
#'   BIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A numeric value giving the value of the BIC.
#' 
#' 
#' @export
#'
bic.smmnonparametric <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  seq <- sequences(seq = seq, E = E)
  
  S <- x$S
  Kmax <- seq$Kmax
  
  Kpar <- .getKpar(x)
  
  n <- sum(sapply(seq$Ym, length))
  
  bic <- -2 * loglik + log(n) * Kpar
  
  return(bic)
  
}