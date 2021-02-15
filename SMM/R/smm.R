#' Function to check if an object is of class `smm`
#' 
#' @description `is.smm` returns `TRUE` if `x` is an object of class `smm`.
#' 
#' @param x An arbitrary R object.
#' 
#' @export
#' 
is.smm <- function(x) {
  inherits(x, "smm")
}


#' Method to get the semi-Markov kernel \eqn{q_{Y}}
#' 
#' @description Computes the semi-Markov kernel \eqn{q_{Y}(k)}
#'   (See proposition 5.1 p.106).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @return An array giving the value of \eqn{q_{Y}(k)} at each time between 0 
#'   and `k`.
#' 
#' @noRd
#' 
.get.qy <- function(x, k, upstates = x$upstates) {
  
  u <- length(upstates)
  
  q <- getKernel(x = x, k = k)
  
  q11 <- q[which(x$states %in% upstates), which(x$states %in% upstates), , drop = FALSE]
  q12 <- q[which(x$states %in% upstates), which(!(x$states %in% upstates)), , drop = FALSE]
  colq12 <- apply(X = q12, MARGIN = c(1, 3), sum)
  
  qy <-
    array(data = sapply(
      X = 1:(k + 1),
      FUN = function(l)
        rbind(cbind(q11[, , l], colq12[, l]), 0)
    ),
    dim = c(u + 1, u + 1, k + 1))
  
  return(qy)
  
}


#' Method to compute the value of \eqn{H}
#' 
#' @description Method to compute the value of \eqn{H} (See equation (3.4) p.46).
#' 
#' @param q An array giving the values of the kernel for a giving time horizon 
#'   \eqn{[0, \dots, k]} (This kernel `q` is the output of the method `getKernel` 
#'   or `.get.qy`).
#' @return An array giving the value of \eqn{H(k)} at each time between 0 
#'   and `k`.
#' 
#' @noRd
#' 
.get.H <- function(q) {
  
  k <- dim(q)[3] - 1
  
  hik <- apply(X = q, MARGIN = c(1, 3), sum)
  Hik <- t(apply(X = hik, MARGIN = 1, cumsum))
  
  H <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1))
  
  for (j in 1:(k + 1)) {
    H[, , j] <- diag(Hik[, j])
  }
  
  return(H)
}


#' Method to compute the value of \eqn{\psi}
#' 
#' @description Method to compute the value of \eqn{\psi}
#'   (See equation (3.16) p.53).
#' 
#' @param q An array giving the values of the kernel for a giving time horizon 
#'   \eqn{[0, \dots, k]} (This kernel `q` is the output of the method `getKernel` 
#'   or `.get.qy`).
#' @return An array giving the value of \eqn{\psi(k)} at each time between 0 
#'   and `k`.
#' 
#' @noRd
#' 
.get.psi <- function(q) {
  
  k <- dim(q)[3] - 1
  
  psi <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1)) # (S, S, k + 1)
  psi[, , 1] <- diag(x = 1, nrow = nrow(q), ncol = ncol(q)) # k = 0
  
  for (j in 1:k) {
    
    psi[, , j + 1] <-
      -Reduce('+', lapply(
        X = 0:(j - 1),
        FUN = function(l)
          psi[, , l + 1] %*% (-q[, , j - l + 1])
      ))
  }
  
  return(psi)
  
}


#' Method to compute the value of \eqn{P}
#' 
#' @description Method to compute the value of \eqn{P} 
#'   (See equation (3.33) p.59).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return An array giving the value of \eqn{P_{i,j}(k)} at each time between 0
#'   and `k` if `var = FALSE`. If `var = TRUE`, a list containing the 
#'   following components:
#'   \itemize{
#'    \item{x: }{an array giving the value of \eqn{P_{ij}(k)} at each time 
#'      between 0 and `k`;}
#'    \item{sigma2: }{an array giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{P}^{2}(i, j, k)}.}
#'  }
#'  
#' @noRd
#' 
.get.P <- function(x, k, states = x$states, var = FALSE, klim = 10000) {
  
  ###########################################################
  # Compute P, the transition function
  ###########################################################
  
  q <- getKernel(x = x, k = k)
  q11 <- q[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  
  psi <- .get.psi(q = q11)
  
  H <- .get.H(q = q)
  H1 <- H[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  B <- array(data = diag(nrow(H1)), dim = dim(H1)) - H1
  
  p <- matrixConvolution(psi, B)
  
  
  ###########################################################
  # Compute the variance (See equation (4.29), p.91)
  ###########################################################
  
  if (var) {
    
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    mu <- meanRecurrenceTimes(x = x, klim = klim)
    
    sigma2 <- varP(mu, q, psi, Psi, H)
    
    return(list(x = p, sigma2 = sigma2))
    
  } else {
    
    return(p)
    
  }
  
}


#' Method to compute the value of \eqn{P_{Y}}
#' 
#' @description Method to compute the value of \eqn{P_{Y}}
#'   (See Proposition 5.1 p.105-106).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @return An array giving the value of \eqn{P_{Y}(k)} at each time between 0
#'   and `k`.
#' 
#' @noRd
#' 
.get.Py <- function(x, k, upstates = x$states) {
  
  ###########################################################
  # Compute P, the transition function
  ###########################################################
  
  qy <- .get.qy(x = x, k = k, upstates = upstates)
  
  psi <- .get.psi(qy)
  
  H <- .get.H(q = qy)
  B <- array(data = diag(nrow(H)), dim = dim(H)) - H
  
  p <- matrixConvolution(psi, B)
  
  return(p)
  
}


#' @export
meanSojournTimes.smm <- function(x, states = x$states, klim = 10000) {
  
  q <- getKernel(x = x, k = klim)
  H1 <- .get.H(q)[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  
  if (dim(H1)[1] != 1) {
    m1 <- apply(1 - apply(H1, 3, diag), 1, sum)
  } else {
    m1 <- sum(1 - H1)
  }
  
  return(m1)
}


#' @export
meanRecurrenceTimes.smm <- function(x, klim = 10000) {
  
  nu <- .stationaryDistribution(ptrans = x$ptrans)
  m <- meanSojournTimes(x = x, klim = klim)
  mu <- sum(nu * m) / nu
  
  return(mu)
  
}


#' @export
reliability.smm <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("'k' must be a positive integer")
  }
  
  if ((!((k >= 0) && ((k %% 1) == 0)))) {
    stop("'k' must be a positive integer")
  }
  
  #############################
  # Checking parameters upstates
  #############################
  
  if (!(is.vector(upstates) && (length(unique(upstates)) == length(upstates)))) {
    stop("The subset of state space 'upstates' is not a vector of unique elements")
  }
  
  if (!all(upstates %in% x$states)) {
    stop("Every element of 'upstates' must be in 'states' ('upstates' is a subet of 'states')")
  }
  
  
  ###########################################################
  # Compute R, the reliability
  ###########################################################
  
  # alpha1 <- x$init[which(x$states %in% upstates)]
  # P11 <- .get.P(x = x, k = k, states = upstates)
  # 
  # reliab <-
  #   apply(X = P11, MARGIN = 3, function(x)
  #     rowSums(t(alpha1) %*% x))
  
  alpha1 <- x$init[which(x$states %in% upstates)]
  u <- length(upstates)
  Py <- .get.Py(x = x, k = k, upstates = upstates)
  
  reliab <-
    apply(X = Py, MARGIN = 3, function(x)
      rowSums(t(c(alpha1, 0)) %*% x[, 1:u, drop = FALSE]))
  
  
  ###########################################################
  # Compute the variance (See equation (5.29), p.116)
  ###########################################################
  
  if (var) {
    
    q <- .get.qy(x = x, k =  k, upstates = upstates)
    Q <- aperm(apply(q, c(1, 2), cumsum), c(2, 3, 1))
    psi <- .get.psi(q = q)
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    H <- .get.H(q)
    mu1 <- meanRecurrenceTimes(x = x, klim = klim)[which(x$states %in% upstates)]
    
    sigma2 <- as.numeric(varR(alpha1, mu1, q, psi, Psi, H, Q))
    
    return(list(x = reliab, sigma2 = sigma2))
    
  } else {
    
    return(reliab)
    
  }

}


#' @export
maintainability.smm <- function(x, k, downstates = x$states, var = FALSE, klim = 10000) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("'k' must be a positive integer")
  }
  
  if ((!((k >= 0) && ((k %% 1) == 0)))) {
    stop("'k' must be a positive integer")
  }
  
  #############################
  # Checking parameters upstates
  #############################
  
  if (!(is.vector(downstates) && (length(unique(downstates)) == length(downstates)))) {
    stop("The subset of state space 'downstates' is not a vector of unique elements")
  }
  
  if (!all(downstates %in% x$states)) {
    stop("Every element of 'downstates' must be in 'states' ('dowstates' is a subset of 'states')")
  }
  
  reliab <- reliability(x = x, k = k, upstates = downstates, var = var, klim = klim)
  
  if (var) {
    
    return(list(x = 1 - reliab$x, sigma2 = reliab$sigma2))
    
  } else {
    
    return(1 - reliab)
    
  }
  
}


#' @export
availability.smm <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("'k' must be a positive integer")
  }
  
  if ((!((k >= 0) && ((k %% 1) == 0)))) {
    stop("'k' must be a positive integer")
  }
  
  #############################
  # Checking parameters upstates
  #############################
  
  if (!(is.vector(upstates) && (length(unique(upstates)) == length(upstates)))) {
    stop("The subset of state space 'upstates' is not a vector of unique elements")
  }
  
  if (!all(upstates %in% x$states)) {
    stop("Every element of 'upstates' must be in 'states' ('upstates' is a subet of 'states')")
  }
  
  
  ###########################################################
  # Compute A, the availability
  ###########################################################
  
  P <- .get.P(x = x, k = k)
  
  avail <-
    apply(X = P[which(x$states %in% upstates), which(x$states %in% upstates), , drop = FALSE],
          MARGIN = 3, function(y)
            rowSums(t(x$init[which(x$states %in% upstates)]) %*% y))
  
  ###########################################################
  # Compute the variance (See equation (5.34), p.118)
  ###########################################################
  
  if (var) {
    
    indices_u <- which(x$states %in% upstates) - 1
    
    alpha <- x$init
    
    q <- getKernel(x = x, k = k)
    Q <- aperm(apply(q, c(1, 2), cumsum), c(2, 3, 1))
    
    psi <- .get.psi(q = q)
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    
    H <- .get.H(q = q)
    mu <- meanRecurrenceTimes(x = x, klim = klim)
    
    sigma2 <- as.numeric(varA(indices_u, alpha, mu, q, psi, Psi, H, Q))
    
    return(list(x = avail, sigma2 = sigma2))
  
  } else {
    
    return(avail)
    
  }
  
}


#' @export
failureRateBMP.smm <- function(x, k, upstates = x$states, var = FALSE, epsilon = 1e-3, klim = 10000) {
  
  ###########################################################
  # Compute \lambda, the BMP-failure rate
  ###########################################################
  
  reliab <- reliability(x = x, k = k, upstates = upstates)
  reliab[reliab < epsilon] <- 0
  
  lbda <- rep.int(0, k + 1)
  lbda[1] <- 1 - reliab[1] # k = 0
  
  for (j in 2:(k + 1)) {
    lbda[j] <- ifelse(reliab[j - 1] != 0, 1 - reliab[j] / reliab[j - 1], 0)
  }
  
  ###########################################################
  # Compute the variance (See equation (5.35), p.119)
  ###########################################################
  
  if (var) {
    
    alpha1 <- x$init[which(x$states %in% upstates)]
    mu1 <- meanRecurrenceTimes(x = x, klim = klim)[which(x$states %in% upstates)]
    
    qy <- .get.qy(x = x, k =  k, upstates = upstates)
    Q <- aperm(apply(qy, c(1, 2), cumsum), c(2, 3, 1))
    
    psi <- .get.psi(q = qy)
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    
    H <- .get.H(qy)
    
    sigma2 <- as.numeric(varBMP(reliab, alpha1, mu1, qy, psi, Psi, H, Q))
    
    return(list(x = lbda, sigma2 = sigma2))
    
  } else {
    
    return(lbda)
    
  }
  
}


#' @export
failureRateRG <- function(x, k, upstates = x$states, var = FALSE, epsilon = 1e-3, klim = 10000) {
  
  lbda <- failureRateBMP(x = x, k = k, upstates = upstates, var = var, epsilon = epsilon, klim = klim)
  
  if (var) {
    
    return(list(x = -log(1 - lbda$x), sigma2 = (1 / (1 - lbda$x) ^ 2) * lbda$sigma2))
    
  } else {
    
    return(-log(1 - lbda))
    
  }
  
}


#' @export
mttf.smm <- function(x, upstates = x$states, klim = 10000, var = FALSE) {
  
  p11 <- x$ptrans[which(x$states %in% upstates), which(x$states %in% upstates), drop = FALSE]
  
  m1 <- meanSojournTimes(x = x, states = upstates, klim = klim)
  
  if (length(m1) != 1) {
    mttf <- as.vector(solve(diag(nrow(p11)) - p11) %*% m1)
  } else {
    mttf <- as.vector(solve(diag(nrow(p11)) - p11) * m1)
  }
  names(mttf) <- upstates
  
  if (var) {
    
    indices_u <- which(x$states %in% upstates) - 1
    indices_d <- which(!(x$states %in% upstates)) - 1
    
    m <- meanSojournTimes(x = x, klim = klim)
    mu <- meanRecurrenceTimes(x = x, klim = klim)
    
    q <- getKernel(x = x, k = klim)
    
    sigma2 <- as.numeric(varMTTF(indices_u, indices_d, m, mu, x$ptrans, q))
    
    return(list(x = mttf, sigma2 = sigma2))
    
  } else {
    
    return(mttf)
    
  }
  
}


#' @export
mttr.smm <- function(x, downstates = x$states, klim = 10000, var = FALSE) {
  
  return(mttf(x = x, upstates = downstates, klim = klim, var = var))
  
}


#' Plot function for an object of class smm
#' 
#' @description Displays the densities for the conditional sojourn time 
#'   distributions depending on the current state `i` and on the next state 
#'   `j`.
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of 
#'   class [smmparametric] or [smmnonparametric]).
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
#' 
plot.smm <- function(x, i = 1, j = 1, klim = NULL, ...) {
  NextMethod(x)
}


#' Simulates semi-Markov chains
#' 
#' @description Simulates sequences from a semi-Markov model.
#' 
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#' 
#' @param object An object inheriting from the S3 class `smm` (an object of 
#'   class [smmparametric] or [smmnonparametric]).
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed `seed` for the random number generator.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' @seealso [smmparametric], [smmnonparametric], [fitsemimarkovmodel]
#' 
#' @export
#' 
simulate.smm <- function(object, nsim = 1, seed = NULL, ...) {
  NextMethod(object)
}
