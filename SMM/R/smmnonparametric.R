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
    stop("The dimension of the matrix ptrans must be equal to (s, s)")
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
    stop("distr must be an array of dimension (s, s, kmax) since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !is.matrix(distr)) {
    stop("distr must be a matrix of dimension (s, kmax) since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !is.vector(distr)) {
    stop("distr must be a vector of length kmax since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !(dim(distr)[1] == s && dim(distr)[2] == s)) {
    stop("distr must be an array of dimension (s, s, kmax) since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !(dim(distr)[1] == s)) {
    stop("distr must be a matrix of dimension (s, kmax) since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (!all(distr >= 0) | !all(distr <= 1)) {
    stop("Probabilities in distr must be between [0, 1]")
  }
  
  if (type.sojourn == "fij") {
    temp <- apply(distr, c(1, 2), sum)
    indexdiag <- seq(1, s * s, by = s + 1)
    
    if (!(all(diag(temp == 0)) && (temp[-indexdiag] == 0 || temp[-indexdiag]  == 1))) {
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
  
  # Add names to the attributes init, ptrans, distr and param for readability
  colnames(ptrans) <- words(length = 1, alphabet = states)
  row.names(ptrans) <- colnames(ptrans)
  names(init) <- colnames(ptrans)
  
  if (is.array(distr) & !(is.matrix(distr))) {
    dimnames(distr) <- rep(list(colnames(ptrans)), 2)
  } else if (is.matrix(distr)) {
    row.names(distr) <- colnames(ptrans)
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

#' Method to get the semi-Markov kernel \eqn{q}
#'
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}.
#' 
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return An array giving the value of \eqn{q_{ij}(k)} at each time between 0 
#'   and `k` if `var = FALSE`. If `var = TRUE`, a list containing the 
#'   following components:
#'   \itemize{
#'    \item{x: }{an array giving the value of \eqn{q_{ij}(k)} at each time 
#'      between 0 and `k`;}
#'    \item{sigma2: }{an array giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(i, j, k)}.}
#'  }
#'
#' @export
#' 
getKernel.smmnonparametric <- function(x, k, var = FALSE, klim = 10000) {
  
  q <- array(data = 0, dim = c(x$s, x$s, k + 1))
  
  if (k <= x$kmax) {
    end <- k
  } else {
    end <- x$kmax
  }
   
  if (x$type.sojourn == "fij") {
    q[, , 2:(end + 1)] <- array(x$ptrans, c(x$s, x$s, end)) * x$distr[, , 1:end]
  } else if (x$type.sojourn == "fi") {
    q[, , 2:(end + 1)] <- array(x$ptrans, c(x$s, x$s, end)) * aperm(array(x$distr[, 1:end], c(x$s, end, x$s)), c(1, 3, 2))
  } else if (x$type.sojourn == "fj") {
    q[, , 2:(end + 1)] <- array(x$ptrans, c(x$s, x$s, end)) * aperm(array(x$distr[, 1:end], c(x$s, end, x$s)), c(3, 1, 2))
  } else if (x$type.sojourn == "f") {
    q[, , 2:(end + 1)] <- array(x$ptrans, c(x$s, x$s, end)) * aperm(array(x$distr[1:end], c(end, x$s, x$s)), c(2, 3, 1))
  }
  
  if (var) {
    
    mu <- .get.mu(x = x, klim = klim)
    sigma2 <- array(data = mu, dim = c(x$s, x$s, k + 1)) * q * (1 - q)
    
    return(list(x = q, sigma2 = sigma2))
    
  } else {
    
    return(q)
    
  }
  
}

#' Loglikelihood
#'
#' @description Computation of the loglikelihood for a semi-Markov model
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @return Value of the loglikelihood.
#' 
#' 
#' @export
#'
loglik.smmnonparametric <- function(x, sequences) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) && all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter sequences should be a list of vectors")
  }
  
  if (!all(unique(unlist(sequences)) %in% x$states)) {
    stop("Some states in the list of observed sequences sequences are not in the state space given by the model x")
  }
  
  
  sequences <- processes(sequences = sequences, states = x$states, verbose = FALSE)
  kmax <- sequences$kmax
  
  if (!(kmax == x$kmax)) {
    stop("kmax of the given sequences is different from the kmax of the estimated SMM model")  
  }
  
  type.sojourn <- x$type.sojourn
  cens.end <- x$cens.end
  
  #############################
  # Let's compute the log-likelihood
  #############################
  
  init <- x$init # Initial distribution
  Nstarti <- sequences$counting$Nstarti
  maskNstarti <- Nstarti != 0 & init != 0
  
  if (!cens.end) {# No censoring
    
    pij <- x$ptrans # Transition matrix
    Nij <- sequences$counting$Nij
    maskNij <- Nij != 0 & pij != 0
    
    # Contribution of the initial distribution and 
    # the transition matrix to the log-likelihood
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
    Y <- lapply(Y, function(x) x - 1)
    Niujv <- getCountingNiuj(Y, U, s, kmax)
    Niu <- apply(Niujv, c(1, 2), sum)
    
    phat <- Niujv / array(Niu, c(s, kmax, s, kmax))
    phat[is.na(phat)] <- 0
    
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
#' @return Value of the AIC.
#' 
#' 
#' @export
#'
aic.smmnonparametric <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
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
#' @return Value of the BIC.
#' 
#' 
#' @export
#'
bic.smmnonparametric <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  n <- sum(unlist(lapply(sequences, length)))
  
  bic <- -2 * loglik + log(n) * kpar
  
  return(bic)
  
}

#' Plot function for an object of class smmnonparametric
#'
#' @description Displays the densities for the conditional sojourn time 
#'   distributions depending on the current state `i` and on the next state 
#'   `j`.
#'
#' @param x An object of class [smmnonparametric][smmnonparametric].
#' @param i An integer giving the current state in the following cases: 
#'   `type.sojourn = "fij"` or `type.sojourn = "fi"`, otherwise, `i` is
#'   ignored.
#' @param j An integer giving the next state in the following cases: 
#'   `type.sojourn = "fij"` or `type.sojourn = "fj"`, otherwise, `j` is
#'   ignored.
#' @param klim An integer giving the limit value for which the density will be 
#'   plotted. If `klim` is `NULL`, then quantile or order 0.95 is used.
#' @param ... Arguments passed to plot.
#' 
#' @export
#'
#' @import graphics
plot.smmnonparametric <- function(x, i = 1, j = 1, klim = NULL, ...) {
  
  #############################
  # Checking parameters i and j
  #############################
  
  if (x$type.sojourn != "f") {
    
    if (x$type.sojourn == "fi") {
      
      if (!((i > 0) & (i <= dim(x$distr)[1]) & ((i %% 1) == 0))) {
        stop(paste0("i must be an integer between 1 and ", dim(x$distr)[1]))
      }
      
    } else if (x$type.sojourn == "fj") {
      
      if (!((j > 0) & (j <= dim(x$distr)[1]) & ((j %% 1) == 0))) {
        stop(paste0("j must be an integer between 1 and ", dim(x$distr)[1]))
      }
      
    } else {
      
      if (i == j) {
        stop(paste0("the conditional distribution for the couple (i = ", i, ", j = ", j, ") doesn't exist.")) 
      }
      
      if (!((i > 0) & (i <= dim(x$distr)[1]) & ((i %% 1) == 0))) {
        stop(paste0("i must be an integer between 1 and ", dim(x$distr)[1]))
      }
      
      if (!((j > 0) & (j <= dim(x$distr)[2]) & ((j %% 1) == 0))) {
        stop(paste0("j must be an integer between 1 and ", dim(x$distr)[2]))
      }
      
    }
    
  }
  
  #############################
  # Checking parameter klim
  #############################
  
  if (!is.null(klim)) {
    if (!((klim > 0) && ((klim %% 1) == 0))) {
      stop("klim must be a strictly positive integer")
    }
  }
  
  
  if (x$type.sojourn == "fij") {
    f <- x$distr[i, j, ]
    ylab <- bquote(f["i = " ~ .(i) ~ ", j = " ~ .(j)])
    main <- paste0("Sojourn time density function for the \n current state i = ", i, " and the next state j = ", j)
  } else if (x$type.sojourn == "fi") {
    f <- x$distr[i, ]
    ylab <- bquote(f["i = " ~ .(i)])
    main <- paste0("Sojourn time density function for the current state i = ", i)
  } else if (x$type.sojourn == "fj") {
    f <- x$distr[j, ]
    ylab <- bquote(f["i = " ~ .(j)])
    main <- paste0("Sojourn time density function for the next state j = ", j)
  } else {
    f <- x$distr
    ylab <- "f"
    main <- paste0("Sojourn time density function")
  }
  
  # Compute the quantile of order alpha if klim is NULL
  alpha <- 0.95
  if (is.null(klim)) {
    cdf <- cumsum(f)
    klim <- ifelse(is.na(which(cdf >= alpha)[1]), x$klim, which(cdf >= alpha)[1])
  }
  
  plot.default(x = 1:klim, y = f[1:klim], xlab = "k", ylab = ylab, ...)
  title(main = main)
}
