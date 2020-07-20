# Function to check if an object is of class smm
is.smm <- function(x) {
  inherits(x, "smm")
}

# Method to get the sojourn time distribution f
.get.f <- function(x, ...) {
  UseMethod(".get.f", x)
}

# Method to get the semi-Markov kernel q
.get.q <- function(x, ...) {
  UseMethod(".get.q", x)
}

# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}

# Method to compute the value of psi (estimator p. 53 (3.16))
.get.psi <- function(x, k, states = x$states) {
  
  q <- .get.q(x, k)[which(x$states %in% states), which(x$states %in% states), ,  drop = FALSE]
  
  psi <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1)) # (S, S, k + 1)
  psi[, , 1] <- diag(x = 1, nrow = nrow(q), ncol = ncol(q)) # k = 0
  
  for (j in 1:k) {
    
    psi[, , j + 1] <-
      -Reduce('+', lapply(
        X = 0:(j - 1),
        FUN = function(l)
          psi[, , l + 1] %*% (-q[, , j - l])
      ))
  }
  
  return(psi)
  
}

# Method to compute the value of H (definition p. 46 (Definition 3.4))
.get.H <- function(x, k) {
  
  q <- .get.q(x, k)
  
  hik <- apply(X = q, MARGIN = c(1, 3), sum)
  Hik <- t(apply(X = hik, MARGIN = 1, cumsum))
  
  H <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1))
  
  for (j in 2:(k + 1)) {
    H[, , j] <- diag(Hik[, j - 1])
  }
  
  return(H)
}

# Method to compute the value of P (estimator p. 59 (3.33))
.get.P <- function(x, k, states = x$states) {
  
  psi <- .get.psi(x, k, states = states)
  
  H1 <- .get.H(x, k)[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  B <- array(data = diag(nrow(H1)), dim = dim(H1)) - H1
  
  p <- array(data = 0, dim = c(nrow(psi), ncol(psi), k + 1)) # (S, S, k + 1)
  p[, , 1] <- diag(nrow(psi)) # k = 0
  
  for (j in 1:k) {
    p[, , j + 1] <- .matrixConvolve(psi[, , 1:(j + 1), drop = FALSE], B[, , 1:(j + 1), drop = FALSE])
  }
  
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
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   reliablity should be computed.
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in U}} of 
#'   length \eqn{\textrm{card}(U) = s1} (see Details).
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @return A vector of length \eqn{k + 1} giving the values of the reliability
#'   for the period \eqn{[0,\dots,k]}.
#' 
#' @export
#'
reliability <- function(x, k, alpha = x$init, upstates = x$states) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("k must be a positive integer")
  }
  
  if ((!((k >= 0) && ((k %% 1) == 0)))) {
    stop("k must be a positive integer")
  }
  
  #############################
  # Checking parameters upstates
  #############################
  
  if (!(is.vector(upstates) && (length(unique(upstates)) == length(upstates)))) {
    stop("The subset of state space upstates is not a vector of unique elements")
  }
  
  if (!all(upstates %in% x$states)) {
    stop("Every element of upstates must be in states (U is a subet of E)")
  }
  
  #############################
  # Checking parameters alpha
  #############################
  
  if (!(is.vector(alpha) && (length(alpha) == length(upstates)))) {
    stop("alpha is not a vector of the same length as U")
  }
  
  if (!(all(alpha >= 0) && all(alpha <= 1))) {
    stop("Probabilities in alpha must be between [0, 1]")
  }
  
  P11 <- .get.P(x, k, states = upstates)
  
  reliab <-
    apply(X = P11, MARGIN = 3, function(x)
      rowSums(t(alpha) %*% x))
  
  return(reliab)
  
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
#'   Thus, we take \eqn{(\alpha_i)_{i \in U} = 0} and we denote by \eqn{T_U} 
#'   the first hitting time of subset \eqn{U}, called the duration of repair or
#'   repair time, that is,
#' 
#'   \deqn{T_U := \textrm{inf}\{ n \in N;\ Z_n \in U\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The maintainability at time \eqn{k \in N} of a discrete-time semi-Markov 
#'   system is 
#'   
#'   \deqn{M(k) = P(T_U \leq k) = 1 - P(T_{U} \geq k) = 1 - P(Z_{n} \in D,\ n = 0,\dots,k).}
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   maintainability should be computed.
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in U}} of 
#'   length \eqn{\textrm{card}(U) = s1} (see Details).
#' @param downstates Vector giving the subset of non-operational states \eqn{D}.
#' @return A vector of length \eqn{k + 1} giving the values of the maintainability
#'   for the period \eqn{[0,\dots,k]}.
#' 
#' @export
#'
maintainability <- function(x, k, alpha = x$init, downstates = x$states) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("k must be a positive integer")
  }
  
  if ((!((k >= 0) && ((k %% 1) == 0)))) {
    stop("k must be a positive integer")
  }
  
  #############################
  # Checking parameters upstates
  #############################
  
  if (!(is.vector(downstates) && (length(unique(downstates)) == length(downstates)))) {
    stop("The subset of state space downstates is not a vector of unique elements")
  }
  
  if (!all(downstates %in% x$states)) {
    stop("Every element of downstates must be in states (U is a subet of E)")
  }
  
  #############################
  # Checking parameters alpha
  #############################
  
  if (!(is.vector(alpha) && (length(alpha) == length(downstates)))) {
    stop("alpha is not a vector of the same length as U")
  }
  
  if (!(all(alpha >= 0) && all(alpha <= 1))) {
    stop("Probabilities in alpha must be between [0, 1]")
  }
  
  P22 <- .get.P(x, k, states = downstates)
  
  maintain <-
    1 - apply(X = P22, MARGIN = 3, function(x)
      rowSums(t(alpha) %*% x))
  
  return(maintain)
  
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
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   availability should be computed.
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in E}} of 
#'   length \eqn{\textrm{card}(E) = s}.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @return A vector of length \eqn{k + 1} giving the values of the availability
#'   for the period \eqn{[0,\dots,k]}.
#' 
#' @export
#'
availability <- function(x, k, alpha = x$init, upstates = x$states) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("k must be a positive integer")
  }
  
  if ((!((k >= 0) && ((k %% 1) == 0)))) {
    stop("k must be a positive integer")
  }
  
  #############################
  # Checking parameters upstates
  #############################
  
  if (!(is.vector(upstates) && (length(unique(upstates)) == length(upstates)))) {
    stop("The subset of state space upstates is not a vector of unique elements")
  }
  
  if (!all(upstates %in% x$states)) {
    stop("Every element of upstates must be in states (U is a subet of E)")
  }
  
  #############################
  # Checking parameters alpha
  #############################
  
  if (!(is.vector(alpha) && (length(alpha) == length(x$states)))) {
    stop("alpha is not a vector of the same length as states")
  }
  
  if (!(all(alpha >= 0) && all(alpha <= 1))) {
    stop("Probabilities in alpha must be between [0, 1]")
  }
  
  if (!(sum(alpha) == 1)) {
    stop("The sum of alpha is not equal to one")
  }
  
  
  P <- .get.P(x, k)
  
  avail <-
    apply(X = P[which(x$states %in% upstates), which(x$states %in% upstates),],
          MARGIN = 3, function(y)
            rowSums(t(alpha[which(x$states %in% upstates)]) %*% y))
  
  return(avail)
  
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
#'   The failure rate at time \eqn{k = 0} is defined by \eqn{\lambda(0) := 1 - R(0)},
#'   with \eqn{R} being the reliability function.
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   BMP-Failure Rate should be computed.
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in U}} of 
#'   length \eqn{\textrm{card}(U) = s1}.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @return A vector of length \eqn{k + 1} giving the values of the BMP-Failure 
#'   Rate for the period \eqn{[0,\dots,k]}.
#' 
#' @export
#'
failureRateBMP <- function(x, k, alpha = x$init, upstates = x$states) {
  
  reliab <- reliability(x = x, k = k, alpha = alpha, upstates = upstates)
  
  lbda <- rep.int(0, k + 1)
  lbda[1] <- 1 - reliab[1] # k = 0
  
  for (j in 2:(k + 1)) {
    lbda[j] <- ifelse(reliab[j - 1] != 0, 1 - reliab[j] / reliab[j - 1], 0)
  }
  
  return(lbda)
  
}

#' RG-Failure Rate Function
#'
#' @description Discrete-time adapted failure rate, proposed by Roy and Gupta (1992).
#'   We call it the RG-failure rate and denote it by \eqn{r(k),\ k \in N}.
#' 
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   RG-Failure Rate should be computed.
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in U}} of 
#'   length \eqn{\textrm{card}(U) = s1}.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @return A vector of length \eqn{k + 1} giving the values of the RG-Failure 
#'   Rate for the period \eqn{[0,\dots,k]}.
#' 
#' @export
#'
failureRateRG <- function(x, k, alpha = x$init, upstates = x$states) {
  
  reliab <- reliability(x = x, k = k, alpha = alpha, upstates = upstates)
  
  r <- rep.int(0, k + 1)
  r[1] <- ifelse(reliab[1] != 0, -log(reliab[1]), 0) # k = 0
  
  for (j in 2:(k + 1)) {
    r[j] <- ifelse(reliab[j - 1] != 0, -log(reliab[j] / reliab[j - 1]), 0)
  }
  
  return(r)
  
}

#' Mean Sojourn Times Function
#'
#' @description Consider a system \eqn{S_{ystem}} starting to work at time 
#'   \eqn{k = 0}. The mean time to failure (MTTF) is defined as the mean 
#'   lifetime.
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
#'   state changes in \eqn{(Z_{n})_{n \in N}} occur.
#'   
#'   The mean sojourn times vector is defined as follows:
#'   
#'   \deqn{m_{i} = E[S_{1} | Z_{0} = j] = \sum_{n \geq 0} (1 - H_{j}(n)),\ i \in E}
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @param k Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector.
#' @return A vector of length \eqn{\textrm{card}(E)} giving the values of the 
#'   mean sojourn times for each state \eqn{i \in E}.
#' 
#' @export
#'
meanSojournTimes <- function(x, states = x$states, k = 10000) {
  
  H1 <- .get.H(x = x, k = k)[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  
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
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in U}} of 
#'   length \eqn{\textrm{card}(U) = s1}.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param k Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function).
#' @return A vector of length \eqn{\textrm{card}(U)} giving the values of the 
#'   mean time to failure for each state \eqn{i \in U}.
#' 
#' @export
#'
mttf <- function(x, alpha = x$init, upstates = x$states, k = 10000) {
  
  p11 <- x$ptrans[which(x$states %in% upstates), which(x$states %in% upstates), drop = FALSE]
  
  m1 <- meanSojournTimes(x = x, states = upstates, k = k)
  
  if (length(m1) != 1) {
    mttf <- t(alpha) %*% solve(diag(nrow(p11)) - p11) %*% diag(m1)
  } else {
    mttf <- t(alpha) %*% solve(diag(nrow(p11)) - p11) * m1
  }
  
  return(mttf)
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
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param alpha Vector of initial distribution \eqn{(\alpha_i)_{i \in D}} of 
#'   length \eqn{\textrm{card}(D) = s - s1}.
#' @param downstates Vector giving the subset of non-operational states \eqn{D}.
#' @param k Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function).
#' @return A vector of length \eqn{\textrm{card}(D)} giving the values of the 
#'   mean time to repair for each state \eqn{i \in D}.
#' 
#' @export
#'
mttr <- function(x, alpha = x$init, downstates = x$states, k = 10000) {
  
  return(mttf(x = x, alpha = alpha, upstates = downstates, k = k))
  
}
