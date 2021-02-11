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

#' Method to get the sojourn time distribution f
#'
#' @description Computes the conditional sojourn time distribution \eqn{f(k)}, 
#'   \eqn{f_{i}(k)}, \eqn{f_{j}(k)} or \eqn{f_{ij}(k)}.
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @return A vector, matrix or array giving the value of \eqn{f} at each time 
#'   between 0 and `k`.
#'
#' @noRd
#' 
.get.f <- function(x, k) {
  UseMethod(".get.f", x)
}

#' Method to get the semi-Markov kernel \eqn{q}
#'
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}.
#' 
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
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
getKernel <- function(x, k, var = FALSE, klim = 10000) {
  UseMethod("getKernel", x)
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

#' Method to get the number of parameters of the semi-Markov chain
#'
#' @description Method to get the number of parameters of the semi-Markov 
#'   chain. This method is useful for the computation of criteria such as AIC 
#'   and BIC.
#' 
#' @param x An object of class `smm`.
#' @return A positive integer giving the number of parameters.
#'
#' @noRd
#' 
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}

#' Method to get the mean recurrence times \eqn{\mu}
#'
#' @description Method to get the mean recurrence times \eqn{\mu}.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}.
#'   
#'   We are interested in investigating the mean recurrence times of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given 
#'   at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}.
#'   
#'   Let \eqn{S = (S_{n})_{n \in N}} denote the successive time points when 
#'   state changes in \eqn{(Z_{n})_{n \in N}} occur and let also 
#'   \eqn{J = (J_{n})_{n \in N}} denote the successively visited states at 
#'   these time points.
#'   
#'   The mean recurrence of an arbitrary state \eqn{j \in E} is given by:
#'   
#'   \deqn{\mu_{jj} = \frac{\sum_{i \in E} \nu(i) m_{i}}{\nu(j)}}
#'   
#'   where \eqn{m_{i}} is the mean sojourn time in state \eqn{i \in E} 
#'   (see [meanSojournTimes] function for the computation).
#' 
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function).
#' @return A vector giving the mean recurrence time 
#'   \eqn{(\mu_{i})_{i \in [1, \dots, s]}}.
#'
#' @export
#' 
meanRecurrenceTimes <- function(x, klim = 10000) {
  
  nu <- .stationaryDistribution(ptrans = x$ptrans)
  m <- meanSojournTimes(x = x, klim = klim)
  mu <- sum(nu * m) / nu
  
  return(mu)
  
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

#' Reliability Function
#'
#' @description Consider a system \eqn{S_{ystem}} starting to function at time 
#'   \eqn{k = 0}. The reliability of \eqn{S_{ystem}} at time \eqn{k \in N} is 
#'   the probability that the system has functioned without failure in the 
#'   period \eqn{[0, k]}.
#'
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the reliability theory of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#' 
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'  \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'
#'  The reliability at time \eqn{k \in N} of a discrete-time semi-Markov system is
#'  
#'  \deqn{R(k) := P(T_D > k) = P(Zn \in U,n = 0,\dots,k)}
#'  
#'  which can be rewritten as follows:
#'  
#'  \deqn{R(k) = \sum_{i \in U} P(Z_0 = i) P(T_D > k | Z_0 = i) = \sum_{i \in U} \alpha_i P(T_D > k | Z_0 = i)}
#'
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   reliability should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the reliability
#'   for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, a list
#'   containing the following components:
#'   \itemize{
#'    \item{x: }{a vector of length \eqn{k + 1} giving the values of the 
#'      reliability for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{R}^{2}(k)}.}
#'  }
#'
#' @export
#'
reliability <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
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

#' Maintainability Function
#'
#' @description For a reparable system \eqn{S_{ystem}} for which the failure 
#'   occurs at time \eqn{k = 0}, its maintainability at time \eqn{k \in N} is 
#'   the probability that the system is repaired up to time \eqn{k}, given that
#'   it has failed at time \eqn{k = 0}.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the maintainability theory of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to fail at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#' 
#'   Thus, we take \eqn{(\alpha_{i} := P(Z_{0} = i))_{i \in U} = 0} and we 
#'   denote by \eqn{T_U} the first hitting time of subset \eqn{U}, called the 
#'   duration of repair or repair time, that is,
#' 
#'   \deqn{T_U := \textrm{inf}\{ n \in N;\ Z_n \in U\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The maintainability at time \eqn{k \in N} of a discrete-time semi-Markov 
#'   system is 
#'   
#'   \deqn{M(k) = P(T_U \leq k) = 1 - P(T_{U} \geq k) = 1 - P(Z_{n} \in D,\ n = 0,\dots,k).}
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   maintainability should be computed.
#' @param downstates Vector giving the subset of non-operational states \eqn{D}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the maintainability
#'   for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, a list
#'   containing the following components:
#'   \itemize{
#'    \item{x: }{a vector of length \eqn{k + 1} giving the values of the 
#'      maintainability for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{M}^{2}(k)}.}
#'  }
#'
#' @export
#'
maintainability <- function(x, k, downstates = x$states, var = FALSE, klim = 10000) {
  
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

#' Availability Function
#'
#' @description The pointwise (or instantaneous) availability of a system 
#'   \eqn{S_{ystem}} at time \eqn{k \in N} is the probability that the system 
#'   is operational at time \eqn{k} (independently of the fact that the system 
#'   has failed or not in \eqn{[0, k)}).
#'
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the availability theory of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given 
#'   at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, 
#'   for a certain \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in 
#'   operating mode \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a 
#'   certain \eqn{j \in D}, means that the system is not operational at time 
#'   \eqn{k} due to the mode of failure \eqn{j} or that the system is under the
#'   repairing mode \eqn{j}.
#'   
#'   The pointwise (or instantaneous) availability of a system \eqn{S_{ystem}} 
#'   at time \eqn{k \in N} is the probability that the system is operational 
#'   at time \eqn{k} (independently of the fact that the system has failed or 
#'   not in \eqn{[0, k)}).
#'   
#'   Thus, the pointwise availability of a semi-Markov system at time 
#'   \eqn{k \in N} is
#'   
#'   \deqn{A(k) = P(Z_k \in U) = \sum_{i \in E} \alpha_i A_i(k),}
#'   
#'   where we have denoted by \eqn{A_i(k)} the conditional availability of the 
#'   system at time \eqn{k \in N}, given that it starts in state \eqn{i \in E},
#'   
#'   \deqn{A_i(k) = P(Z_k \in U | Z_0 = i).}
#'
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   availability should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Optional. A logical value indicating whether or not the variance 
#'   of the estimator should be computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the availability
#'   for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, a list
#'   containing the following components:
#'   \itemize{
#'    \item{x: }{a vector of length \eqn{k + 1} giving the values of the 
#'      availability for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{A}^{2}(k)}.}
#'  }
#'
#' @export
#'
availability <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
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

#' BMP-Failure Rate Function
#'
#' @description Consider a system \eqn{S_{ystem}} starting to work at time 
#'   \eqn{k = 0}. The BMP-failure rate at time \eqn{k \in N} is the conditional 
#'   probability that the failure of the system occurs at time \eqn{k}, given 
#'   that the system has worked until time \eqn{k - 1}.
#'
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the failure rate theory of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   The BMP-failure rate at time \eqn{k \in N} is the conditional probability 
#'   that the failure of the system occurs at time \eqn{k}, given that the 
#'   system has worked until time \eqn{k - 1}.
#'
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'   \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   For a discrete-time semi-Markov system, the failure rate at time 
#'   \eqn{k \geq 1} has the expression:
#'   
#'   \deqn{\lambda(k) := P(T_{D} = k | T_{D} \geq k)}
#'   
#'   We can rewrite it as follows :
#'   
#'   \deqn{\lambda(k) = 1 - \frac{R(k)}{R(k - 1)},\ \textrm{if } R(k - 1) \neq 0;\ \lambda(k) = 0, \textrm{otherwise}}
#'   
#'   The failure rate at time \eqn{k = 0} is defined by \eqn{\lambda(0) := 1 - R(0)},
#'   with \eqn{R} being the reliability function (see [reliability][reliability] 
#'   function).
#'   
#'   The calculation of the reliability \eqn{R} involves the computation of 
#'   many convolutions. It implies that the computation error, may be higher 
#'   (in value) than the "true" reliability itself for reliability close to 0.
#'   In order to avoid inconsistent values of the BMP-failure rate, we use the 
#'   following formula:
#'   
#'   \deqn{\lambda(k) = 1 - \frac{R(k)}{R(k - 1)},\ \textrm{if } R(k - 1) \geq \epsilon;\ \lambda(k) = 0, \textrm{otherwise}}
#'   
#'   with \eqn{\epsilon}, the threshold, the parameter `epsilon` in the 
#'   function `failureRateBMP`.
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   BMP-failure rate should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param epsilon Value of the reliability above which the latter is supposed 
#'   to be 0 because of computation errors (see Details).
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the BMP-failure 
#'   rate for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, 
#'   a list containing the following components:
#'   \itemize{
#'    \item{x: }{a vector of length \eqn{k + 1} giving the values of the 
#'      BMP-failure rate for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{\lambda}^{2}(k)}.}
#'  }
#'
#' @export
#'
failureRateBMP <- function(x, k, upstates = x$states, var = FALSE, epsilon = 1e-3, klim = 10000) {
  
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

#' RG-Failure Rate Function
#'
#' @description Discrete-time adapted failure rate, proposed by D. Roy and 
#'   R. Gupta. Classification of discrete lives. Microelectronics Reliability, 
#'   32(10):1459--1473, 1992.
#'   We call it the RG-failure rate and denote it by \eqn{r(k),\ k \in N}.
#' 
#' @details Expressing \eqn{r(k)} in terms of the reliability \eqn{R} we obtain 
#'   that the RG-failure rate function for a discrete-time system is given by:
#'   
#'   \deqn{r(k) = - \ln \frac{R(k)}{R(k - 1)},\ \textrm{if } k \geq 1;\ r(k) = - \ln R(0),\ \textrm{if } k = 0}
#'   
#'   for \eqn{R(k) \neq 0}. If \eqn{R(k) = 0}, we set \eqn{r(k) := 0}.
#'   
#'   Note that the RG-failure rate is related to the BMP-failure rate 
#'   (see [failureRateBMP] function) by:
#'   
#'   \deqn{r(k) = - \ln (1 - \lambda(k)),\ k \in N}
#'   
#'   The computation of the RG-failure rate is based on the [failureRateBMP] 
#'   function (See [failureRateBMP] for details about the parameter `epsilon`).
#' 
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   RG-failure rate should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param epsilon Value of the reliability above which the latter is supposed 
#'   to be 0 because of computation errors (see Details).
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the RG-failure 
#'   rate for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, 
#'   a list containing the following components:
#'   \itemize{
#'    \item{x: }{a vector of length \eqn{k + 1} giving the values of the 
#'      RG-failure rate for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{r}^{2}(k)}.}
#'  }
#'
#' @export
#'
failureRateRG <- function(x, k, upstates = x$states, var = FALSE, epsilon = 1e-3, klim = 10000) {
  
  lbda <- failureRateBMP(x = x, k = k, upstates = upstates, var = var, epsilon = epsilon, klim = klim)
  
  if (var) {
    
    return(list(x = -log(1 - lbda$x), sigma2 = (1 / (1 - lbda$x) ^ 2) * lbda$sigma2))
    
  } else {
    
    return(-log(1 - lbda))
    
  }
  
}

#' Mean Sojourn Times Function
#'
#' @description The mean sojourn time is the mean time spent in each state.
#'
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}.
#'   
#'   We are interested in investigating the mean sojourn times of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given 
#'   at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}.
#'   
#'   Let \eqn{S = (S_{n})_{n \in N}} denote the successive time points when 
#'   state changes in \eqn{(Z_{n})_{n \in N}} occur and let also 
#'   \eqn{J = (J_{n})_{n \in N}} denote the successively visited states at 
#'   these time points.
#'   
#'   The mean sojourn times vector is defined as follows:
#'   
#'   \deqn{m_{i} = E[S_{1} | Z_{0} = j] = \sum_{k \geq 0} (1 - P(Z_{n + 1} \leq k | J_{n} = j)) = \sum_{k \geq 0} (1 - H_{j}(k)),\ i \in E}
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function).
#' @return A vector of length \eqn{\textrm{card}(E)} giving the values of the 
#'   mean sojourn times for each state \eqn{i \in E}.
#' 
#' @export
#'
meanSojournTimes <- function(x, states = x$states, klim = 10000) {
  
  q <- getKernel(x = x, k = klim)
  H1 <- .get.H(q)[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  
  if (dim(H1)[1] != 1) {
    m1 <- apply(1 - apply(H1, 3, diag), 1, sum)
  } else {
    m1 <- sum(1 - H1)
  }
  
  return(m1)
}

#' Mean Time To Failure (MTTF) Function
#'
#' @description Consider a system \eqn{S_{ystem}} starting to work at time 
#'   \eqn{k = 0}. The mean time to failure (MTTF) is defined as the mean 
#'   lifetime.
#'
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the mean time to failure theory of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'   \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The mean time to failure (MTTF) is defined as the mean lifetime, i.e., the
#'   expectation of the hitting time to down set \eqn{D},
#'   
#'   \deqn{MTTF = E[T_{D}]}
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @param var Optional. A logical value indicating whether or not the variance 
#'   of the estimator should be computed.
#' @return If `var = FALSE`, a vector of length \eqn{\textrm{card}(U)} giving 
#'   the values of the mean time to failure for each state \eqn{i \in U}. If 
#'   `var = TRUE`, a list containing the following components:
#'   \itemize{
#'    \item{x: }{a vector of length \eqn{\textrm{card}(U) = s_{1}} giving 
#'      the values of the mean time to failure for each state \eqn{i \in U}.}
#'    \item{sigma2: }{the variances of the estimator for each estimation of the 
#'      mean time to failure \eqn{\sigma_{MTTF_{i}}^{2}(k)};}
#'  }
#' 
#' @export
#'
mttf <- function(x, upstates = x$states, klim = 10000, var = FALSE) {
  
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

#' Mean Time To Repair (MTTR) Function
#'
#' @description Consider a system \eqn{S_{ystem}} starting to fail at time 
#'   \eqn{k = 0}. The mean time to repair (MTTR) is defined as the mean of the 
#'   repair duration.
#'
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the mean time to repair theory of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to fail at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   Let \eqn{T_U} denote the first passage time in subset \eqn{U}, called the 
#'   duration of repair or repair time, i.e.,
#'   
#'   \deqn{T_U := \textrm{inf}\{ n \in N;\ Z_n \in U\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The mean time to repair (MTTR) is defined as the mean of the repair 
#'   duration, i.e., the expectation of the hitting time to up set \eqn{U},
#'   
#'   \deqn{MTTR = E[T_{U}]}
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param downstates Vector giving the subset of non-operational states \eqn{D}.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @param var Optional. A logical value indicating whether or not the variance 
#'   of the estimator should be computed.
#' @return If `var = FALSE`, a vector of length \eqn{\textrm{card}(D)} giving 
#'   the values of the mean time to repair for each state \eqn{i \in D}. If 
#'   `var = TRUE`, a list containing the following components:
#'   \itemize{
#'  \item{x: }{a vector of length \eqn{\textrm{card}(D)} giving the values
#'      of the mean time to repair for each state \eqn{i \in D}.}
#'    \item{sigma2: }{the variances of the estimator for each estimation of the 
#'      mean time to repair \eqn{\sigma_{MTTR_{i}}^{2}(k)};}
#'  }
#' 
#' @export
#'
mttr <- function(x, downstates = x$states, klim = 10000, var = FALSE) {
  
  return(mttf(x = x, upstates = downstates, klim = klim, var = var))
  
}
