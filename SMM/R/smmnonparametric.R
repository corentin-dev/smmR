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
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !is.matrix(laws)) {
    stop("laws must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !is.vector(laws)) {
    stop("laws must be a vector of length Kmax since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !(dim(laws)[1] == S && dim(laws)[2] == S)) {
    stop("laws must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !(dim(laws)[1] == S)) {
    stop("laws must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (!all(laws >= 0) || !all(laws <= 1)) {
    stop("Probabilities in laws must be between [0, 1]")
  }
  
  
  temp <- apply(laws, c(1, 2), sum)
  indexdiag <- seq(1, S * S, by = S + 1)
  if (type.sojourn == "fij" && !(all(diag(temp == 0)) && all(temp[-indexdiag] == 1))) {
    stop("laws is not a stochastic matrix")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !(all(apply(laws, 1, sum) == 1))) {
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

is.smmnonparametric <- function(x) {
  inherits(x, "smmnonparametric")
}

.get.f.smmnonparametric <- function(x, Kmax) {
  
  f <- x$laws
  type.sojourn <- x$type.sojourn
  S <- x$S
  
  if (type.sojourn == "fij") {
    
    fijk <- f
    
  } else if (type.sojourn == "fi") {
    
    f <- rep(as.vector(t(f)), each = S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fk <- array(as.vector(t(fmat)), c(S, S, Kmax))
    fijk <- apply(X = fk, MARGIN =  c(1, 3), FUN =  t)
    
  } else if (type.sojourn == "fj") {
    
    f <- rep(as.vector(t(f)), each = S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(S, S, Kmax))
    
  } else {
    
    f <- rep(f, each = S * S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(S, S, Kmax))
    
  }
  
  return(fijk)
}

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