#' Function to check if an object is of class `smm`
#'
#' @description `is.smm` returns `TRUE` if `x` is an object of class `smm`.
#' 
#' @param x An arbitrary R object.
#'
#' @noRd
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
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}
#'   (Theorem 4.2 p.82).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return An array giving the value of \eqn{q_{ij}(k)} at each time between 0 
#'   and `k` if `var = FALSE`. If `var = TRUE`, a list containing the 
#'   following components:
#'   \itemize{
#'    \item{q: }{an array giving the value of \eqn{q_{ij}(k)} at each time 
#'      between 0 and `k`;}
#'    \item{sigma2: }{an array giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(i, j, k)}.}
#'  }
#'
#' @noRd
.get.q <- function(x, k, var = FALSE, klim = 10000) {
  UseMethod(".get.q", x)
}

#' Method to get the semi-Markov kernel \eqn{q_{Y}}
#'
#' @description Computes the semi-Markov kernel \eqn{q_{Y}(k)}
#'   (Proposition 5.1 p.106).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @return An array giving the value of \eqn{q_{Y}(k)} at each time between 0 
#'   and `k`.
#'
#' @noRd
.get.qy <- function(x, k, upstates = x$upstates) {
  
  u <- length(upstates)
  
  q <- .get.q(x = x, k = k)
  
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
#' @description Method to get the mean recurrence times \eqn{\mu}
#'   (Proposition 3.6 p.57).
#' 
#' @param x An object of class `smm`.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean sojourn times vector 
#'   \eqn{(m_i)_{i \in [1,\dots,S]}}.
#' @return A vector giving the mean recurrence time 
#'   \eqn{(\mu_{i})_{i \in [1,\dots,S]}}.
#'
#' @noRd
.get.mu <- function(x, klim = 10000) {
  
  nu <- .stationaryDistribution(ptrans = x$ptrans)
  m <- meanSojournTimes(x = x, klim = klim)
  mu <- sum(nu * m) / nu
  
  return(mu)
  
}

# # Method to compute the value of psi (estimator p.53 (3.16))
# #' @export
# .get.psi <- function(x, k, states = x$states) {
# 
#   q <- .get.q(x = x, k = k)[which(x$states %in% states), which(x$states %in% states), ,  drop = FALSE]
# 
#   psi <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1)) # (S, S, k + 1)
#   psi[, , 1] <- diag(x = 1, nrow = nrow(q), ncol = ncol(q)) # k = 0
# 
#   for (j in 1:k) {
# 
#     psi[, , j + 1] <-
#       -Reduce('+', lapply(
#         X = 0:(j - 1),
#         FUN = function(l)
#           psi[, , l + 1] %*% (-q[, , j - l + 1])
#       ))
#   }
# 
#   return(psi)
# 
# }

# # Method to compute the value of psi (estimator p.53 (3.16))
# #' @export
# .get.psiy <- function(x, k, downstates = x$states) {
# 
#   q <- .get.qy(x = x, k = k, downstates = downstates)
# 
#   psi <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1)) # (S, S, k + 1)
#   psi[, , 1] <- diag(x = 1, nrow = nrow(q), ncol = ncol(q)) # k = 0
# 
#   for (j in 1:k) {
# 
#     psi[, , j + 1] <-
#       -Reduce('+', lapply(
#         X = 0:(j - 1),
#         FUN = function(l)
#           psi[, , l + 1] %*% (-q[, , j - l + 1])
#       ))
#   }
# 
#   return(psi)
# 
# }

#' Method to compute the value of \eqn{\psi}
#'
#' @description Method to compute the value of \eqn{\psi}
#'   (estimator p.53 (3.16)).
#' 
#' @param q An array giving the values of the kernel for a giving time horizon 
#'   \eqn{[0, \dots, k]} (This kernel `q` is the output of the method `.get.q` 
#'   or `.get.qy`).
#' @return An array giving the value of \eqn{\psi(k)} at each time between 0 
#'   and `k`.
#'
#' @noRd
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

# # Method to compute the value of H (definition p.46 (Definition 3.4))
# #' @export
# .get.H <- function(x, k) {
#   
#   q <- .get.q(x = x, k = k)
#   
#   hik <- apply(X = q, MARGIN = c(1, 3), sum)
#   Hik <- t(apply(X = hik, MARGIN = 1, cumsum))
#   
#   H <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1))
#   
#   for (j in 1:(k + 1)) {
#     H[, , j] <- diag(Hik[, j])
#   }
#   
#   return(H)
# }

# # Method to compute the value of H (definition p.46 (Definition 3.4))
# #' @export
# .get.Hy <- function(x, k, downstates = x$states) {
#   
#   q <- .get.qy(x = x, k = k, downstates = downstates)
#   
#   hik <- apply(X = q, MARGIN = c(1, 3), sum)
#   Hik <- t(apply(X = hik, MARGIN = 1, cumsum))
#   
#   H <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1))
#   
#   for (j in 1:(k + 1)) {
#     H[, , j] <- diag(Hik[, j])
#   }
#   
#   return(H)
# }

#' Method to compute the value of \eqn{H}
#'
#' @description Method to compute the value of \eqn{H} (Definition 3.4 p.46).
#' 
#' @param q An array giving the values of the kernel for a giving time horizon 
#'   \eqn{[0, \dots, k]} (This kernel `q` is the output of the method `.get.q` 
#'   or `.get.qy`).
#' @return An array giving the value of \eqn{H(k)} at each time between 0 
#'   and `k`.
#'
#' @noRd
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
#' @description Method to compute the value of \eqn{P} (estimator p.59 (3.33)).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return An array giving the value of \eqn{P_{i,j}(k)} at each time between 0
#'   and `k` if `var = FALSE`. If `var = TRUE`, a list containing the 
#'   following components:
#'   \itemize{
#'    \item{p: }{an array giving the value of \eqn{P_{ij}(k)} at each time 
#'      between 0 and `k`;}
#'    \item{sigma2: }{an array giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{P}^{2}(i, j, k)}.}
#'  }
#'
#' @noRd
.get.P <- function(x, k, states = x$states, var = FALSE, klim = 10000) {
  
  ###########################################################
  # Compute P, the transition function
  ###########################################################
  
  q <- .get.q(x = x, k = k)
  q11 <- q[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  
  psi <- .get.psi(q = q11)
  
  H <- .get.H(q = q)
  H1 <- H[which(x$states %in% states), which(x$states %in% states), , drop = FALSE]
  B <- array(data = diag(nrow(H1)), dim = dim(H1)) - H1
  
  p <- array(data = 0, dim = c(nrow(psi), ncol(psi), k + 1)) # (S, S, k + 1)
  p[, , 1] <- diag(nrow(psi)) # k = 0
  
  for (t in 1:k) {
    p[, , t + 1] <- .matrixConvolve(psi[, , 1:(t + 1), drop = FALSE], B[, , 1:(t + 1), drop = FALSE])
  }
  
  
  ###########################################################
  # Compute the variance (equation (4.29), p.91)
  # The decomposition of the variance is as follows:
  #
  # \sigma_{P}^{2}(i, j, k)) = \sum_{m = 1}^{s} \mu_{mm} \left\{ \sum_{r = 1}^{s} \underbrace{\underbrace{\left[\delta_{mj}\Psi_{ij} - \underbrace{(1 - H_{j}) * \psi_{im} * \psi_{rj}}_{\text{part11}} \right]^2}_{\text{part12}} * \ q_{mr}(k)}_{\text{part1}} - \left[ \underbrace{\delta_{mj} \psi_{ij} * H_{m}(k)}_{\text{part22}} - \sum_{r = 1}^{s} \underbrace{(1 - H_{j}) * \psi_{im} * \psi_{rj} * q_{mr}}_{\text{part21}} \right]^{2}(k) \right\}
  #
  ###########################################################
  
  if (var) {
    
    u <- length(states)
    indices_u <- which(x$states %in% states)
    
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    mu <- .get.mu(x = x, klim = klim)[indices_u]
    
    
    convolpsi <- array(data = 0, dim = c(u, u, u, u, k + 1)) # (i, j, m, r, k)
    
    part11 <- array(data = 0, dim = c(u, u, u, u, k + 1)) # (1 - H_{j}) * \psi_{im} * \psi_{rj}
    part12 <- array(data = 0, dim = c(u, u, u, u, k + 1)) # [\delta_{mj}\Psi_{ij} - (1 - H_{j}) * \psi_{im} * \psi_{rj}]^2
    part1 <- array(data = 0, dim = c(u, u, u, u, k + 1)) # [\delta_{mj}\Psi_{ij} - (1 - H_{j}) * \psi_{im} * \psi_{rj}]^2 * q_{mr}(k)
    
    part21 <- array(data = 0, dim = c(u, u, u, u, k + 1)) # (1 - H_{j}) * \psi_{im} * \psi_{rj} * q_{mr}
    part22 <- array(data = 0, dim = c(u, u, u, k + 1)) # (i, j, m, k) \delta_{mj} \psi_{ij} * H_{m}(k)
    
    for (t in 0:k) {
      
      for (i in 1:u) {
        
        for (j in 1:u) {
          
          convolpsi[i, j, , , t + 1] <-
            outer(
              X = 1:u,
              Y = 1:u,
              FUN = function(m, r)
                Reduce('+', lapply(
                  X = 0:t,
                  FUN = function(h)
                    psi[i, m, h + 1] * psi[r, j, t - h + 1]
                ))
            )
          
          part11[i, j, , , t + 1] <- apply(X = convolpsi[i, j, , , 1:(t + 1)], MARGIN = c(1, 2), FUN = function(elt) Reduce('+', elt * (1 - H[j, j, (t + 1):1])))
          
          delta_mj <- ifelse(1:u == j, 1, 0)
          
          part12[i, j, , , t + 1] <- (delta_mj * Psi[i, j, t + 1] - part11[i, j, , , t + 1]) ^ 2
          
          part1[i, j, , , t + 1] <- apply(X = part12[i, j, , , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum)
          
          part21[i, j, , , t + 1] <- apply(X = part11[i, j, , , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum)
          
          part22[i, j, which(delta_mj != 0), t + 1] <- Reduce('+', psi[i, j, 1:(t + 1)] * H[delta_mj != 0, delta_mj != 0, (t + 1):1])
          
        }
        
      }
    
    }
      
    sigma2 <- apply(X = apply(part1, c(1, 2, 3, 5), sum) - 
                      (part22 - apply(part21, c(1, 2, 3, 5), sum)) ^ 2, MARGIN = c(1, 2, 4), function(elt) sum(elt * mu))
   
    return(list(p = p, sigma2 = sigma2))
     
  } else {
    
    return(p)
      
  }
  
}

#' Method to compute the value of \eqn{P_{Y}}
#'
#' @description Method to compute the value of \eqn{P_{Y}}
#'   (Proposition 5.1 p.105-106).
#' 
#' @param x An object of class `smm`.
#' @param k A positive integer giving the time horizon.
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @return An array giving the value of \eqn{P_{Y}(k)} at each time between 0
#'   and `k`.
#'
#' @noRd
.get.Py <- function(x, k, upstates = x$states) {
  
  ###########################################################
  # Compute P, the transition function
  ###########################################################
  
  qy <- .get.qy(x = x, k = k, upstates = upstates)
  
  psi <- .get.psi(qy)
  
  H <- .get.H(q = qy)
  B <- array(data = diag(nrow(H)), dim = dim(H)) - H
  
  p <- array(data = 0, dim = c(nrow(psi), ncol(psi), k + 1)) # (U + 1, U + 1, k + 1)
  p[, , 1] <- diag(nrow(psi)) # k = 0
  
  for (t in 1:k) {
    p[, , t + 1] <- .matrixConvolve(psi[, , 1:(t + 1), drop = FALSE], B[, , 1:(t + 1), drop = FALSE])
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
#'   reliability should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the reliability
#'   for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, a list
#'   containing the following components:
#'   \itemize{
#'    \item{reliab: }{a vector of length \eqn{k + 1} giving the values of the 
#'      reliability for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(k)}.}
#'  }
#'
#' @export
#'
reliability <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
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
  # Compute the variance (equation (5.29), p.116)
  # 
  # Be careful: 
  # 
  # In the formula (5.29), we use q_{Y} (Proposition 5.1 p.105-106) instead 
  # of q, and every others quantities such as \psi, \Psi,\dots derive from q_{Y}
  # 
  # 
  # The decomposition of the variance is as follows:
  # 
  # \sigma_{R}^{2}(k) = \sum_{i = 1}^{s} \mu_{ii} \left\{ \sum_{j = 1}^{s} \underbrace{\underbrace{\left[ D^{U}_{ij} - \mathbb{1}_{i \in U} \sum_{t \in U} \alpha_{t} \Psi_{ti} \right]^{2}}_{\text{part11}} * q_{ij}(k)}_{\text{part1}} - \left[ \sum_{j = 1}^{s} \underbrace{\left( \underbrace{D^{U}_{ij} * q_{ij}}_{\text{part22}} - \mathbb{1}_{i \in U} \sum_{t \in U} \alpha_{t} \underbrace{\psi_{ti} * Q_{ij}}_{\text{part21}} \right)}_{\text{part2}} \right]^{2}(k) \right\}
  # 
  # D^{U}_{ij} := \underbrace{\sum_{n \in U} \sum_{r \in U} \underbrace{\alpha_{n} \psi_{ni} * \psi_{jr} * (\text{I} - diag(\text{Q.1}))_{rr}}_{partduij}}_{duij}
  # 
  ###########################################################
  
  if (var) {
 
    u <- length(upstates)
    
    q <- .get.qy(x = x, k =  k, upstates = upstates)
    Q <- aperm(apply(q, c(1, 2), cumsum), c(2, 3, 1))
    
    psi <- .get.psi(q = q)
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    
    H <- .get.H(q)
    
    mu <- .get.mu(x = x, klim = klim)


    indicator <- matrix(data = 1:(u + 1) %in% 1:u, nrow = u + 1, ncol = u + 1) # indicator matrix \mathbb{1}_{i \in U}
        
    convolpsi <- array(data = 0, dim = c(u + 1, u + 1, u + 1, u + 1, k + 1)) # (i, j, n, r, k)
    
    partduij <- array(data = 0, dim = c(u + 1, u + 1, u, u, k + 1)) # # \alpha_{n} \psi_{ni} * \psi_{jr} * (\text{I} - diag(\text{Q.1}))_{rr}
    duij <- array(data = 0, dim = c(u + 1, u + 1, k + 1))
    
    part11 <- array(data = 0, dim = c(u + 1, u + 1, k + 1))
    part1 <- array(data = 0, dim = c(u + 1, u + 1, k + 1))
    
    part21 <- array(data = 0, dim = c(u, u + 1, u + 1, k + 1)) # (t, i, j, k)
    part22 <- array(data = 0, dim = c(u + 1, u + 1, k + 1))
    
    part2 <- array(data = 0, dim = c(u + 1, u + 1, k + 1))
    
    for (t in 0:k) {
      
      for (i in 1:(u + 1)) {
        
        for (j in 1:(u + 1)) {
          
          convolpsi[i, j, , , t + 1] <- outer(X = 1:(u + 1), Y = 1:(u + 1), FUN = function(n, r)
            Reduce('+', lapply(
              X = 0:t,
              FUN = function(h) psi[n, i, h + 1] * psi[j, r, t - h + 1]))
          )
          
          for (r in 1:u) {
            for (n in 1:u) {
              partduij[i, j, n, r, t + 1] <- sum(convolpsi[i, j, n, r, 1:(t + 1)] * (1 - H[r, r, (t + 1):1]))
              partduij[i, j, n, r, t + 1] <- partduij[i, j, n, r, t + 1] * alpha1[n]
            }
          }
          
          duij[, , t + 1] <- apply(partduij[, , , , t + 1], c(1, 2), sum)
          
          
          for (r in 1:u) {
            part21[r, i, j, t + 1] <- sum(psi[r, i, 1:(t + 1)] * Q[i, j, (t + 1):1])
          }
          
        }
      }
      
      part11[, , t + 1] <- (duij[, , t + 1] - indicator * (t(Psi[, , t + 1]) %*% c(alpha1, 0) %*% t(rep.int(times = u + 1, x = 1)))) ^ 2
      part1[, , t + 1] <- apply(X = part11[, , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum)
      
      part22[, , t + 1] <- apply(X = duij[, , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum)
      # part2[, , t + 1,] <- part22[, , t + 1] - indicator * apply(part21[, , , t + 1], c(2, 3), function(elt) elt %*% alpha1)
      part2[, , t + 1] <- part22[, , t + 1] - indicator * drop(apply(part21[, , , t + 1, drop = FALSE], c(2, 3, 4), function(elt) elt %*% alpha1))
      
    }
    
    sigma2 <-
      apply(apply(part1[1:u, , , drop = FALSE], c(1, 3), sum) - (apply(part2[1:u, , , drop = FALSE], c(1, 3), sum)) ^ 2, 2, function(elt)
        sum(elt * mu[which(x$states %in% upstates)]))
    
    return(list(reliab = reliab, sigma2 = sigma2))
    
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
#' @param downstates Vector giving the subset of non-operational states \eqn{D}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the maintainability
#'   for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, a list
#'   containing the following components:
#'   \itemize{
#'    \item{reliab: }{a vector of length \eqn{k + 1} giving the values of the 
#'      maintainability for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(k)}.}
#'  }
#'
#' @export
#'
maintainability <- function(x, k, downstates = x$states, var = FALSE, klim = 10000) {
  
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
  
  tmp <- reliability(x = x, k = k, upstates = downstates, var = var, klim = klim)
  
  if (var) {
    
    maintain <- list(maintain = 1 - tmp$reliab, sigma2 = tmp$sigma2)
    
  } else {
    
    maintain <- 1 - tmp
    
  }
  
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
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Optional. A logical value indicating whether or not the variance 
#'   of the estimator should be computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean recurrence times vector for the asymptotic variance.
#' @return A vector of length \eqn{k + 1} giving the values of the availability
#'   for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, a list
#'   containing the following components:
#'   \itemize{
#'    \item{avail: }{a vector of length \eqn{k + 1} giving the values of the 
#'      availability for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(k)}.}
#'  }
#'
#' @export
#'
availability <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
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
  
  
  ###########################################################
  # Compute A, the availability
  ###########################################################
  
  P <- .get.P(x = x, k = k)
  
  avail <-
    apply(X = P[which(x$states %in% upstates), which(x$states %in% upstates), , drop = FALSE],
          MARGIN = 3, function(y)
            rowSums(t(x$init[which(x$states %in% upstates)]) %*% y))
  
  ###########################################################
  # Compute the variance (Theorem 5.2. equation (5.34), p.118)
  # The decomposition of the variance is as follows:
  #
  # \sigma_{A}^{2}(k) = \sum_{i = 1}^{s} \mu_{ii} \left\{ \sum_{j = 1}^{s} \underbrace{\underbrace{\left[ D_{ij} - \mathbb{1}_{i \in U} \sum_{t = 1}^{s} \alpha_{t} \Psi_{ti} \right]^{2}}_{\text{part11}} * q_{ij}(k)}_{\text{part1}} - \left[ \sum_{j = 1}^{s} \underbrace{\left( \underbrace{D_{ij} * q_{ij}}_{\text{part22}} - \mathbb{1}_{i \in U} \sum_{t = 1}^{s} \alpha_{t} \underbrace{\psi_{ti} * Q_{ij}}_{\text{part21}} \right)}_{\text{part2}} \right]^{2}(k) \right\}
  #
  # D_{ij} := \underbrace{\sum_{n = 1}^{s} \sum_{r \in U} \underbrace{\alpha_{n} \psi_{ni} * \psi_{jr} * (\text{I} - diag(\text{Q.1}))_{rr}}_{partdij}}_{dij}
  # 
  ###########################################################
  
  if (var) {
  
    q <- .get.q(x = x, k = k)
    Q <- aperm(apply(q, c(1, 2), cumsum), c(2, 3, 1))
    
    psi <- .get.psi(q = q)
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    
    H <- .get.H(q = q)
    mu <- .get.mu(x = x, klim = klim)
    
    u <- length(upstates)
    indices_u <- which(x$states %in% upstates)
    
    
    indicator <- matrix(data = 1:x$s %in% indices_u, nrow = x$s, ncol = x$s) # indicator matrix \mathbb{1}_{i \in U}
    
    convolpsi <- array(data = 0, dim = c(x$s, x$s, x$s, u, k + 1)) # (i, j, n, r, k)
    
    partdij <- array(data = 0, dim = c(x$s, x$s, x$s, u, k + 1)) # \alpha_{n} \psi_{ni} * \psi_{jr} * (\text{I} - diag(\text{Q.1}))_{rr}
    dij <- array(data = 0, dim = c(x$s, x$s, k + 1))
    
    part11 <- array(data = 0, dim = c(x$s, x$s, k + 1))
    part1 <- array(data = 0, dim = c(x$s, x$s, k + 1))
    
    part21 <- array(data = 0, dim = c(x$s, x$s, x$s, k + 1)) # (t, i, j, k)
    part22 <- array(data = 0, dim = c(x$s, x$s, k + 1))
    
    part2 <- array(data = 0, dim = c(x$s, x$s, k + 1))
    
    for (t in 0:k) {
      
      for (i in 1:x$s) {
        
        for (j in 1:x$s) {
          
          convolpsi[i, j, , , t + 1] <- outer(X = 1:x$s, Y = 1:u, FUN = function(n, r)
            Reduce('+', lapply(
              X = 0:t,
              FUN = function(h) psi[n, i, h + 1] * psi[j, indices_u[r], t - h + 1]))
          )
          
          for (r in 1:u) {
            for (n in 1:x$s) {
              
              partdij[i, j, n, r, t + 1] <- sum(convolpsi[i, j, n, r, 1:(t + 1)] * (1 - H[indices_u[r], indices_u[r], (t + 1):1]))
              partdij[i, j, n, r, t + 1] <- partdij[i, j, n, r, t + 1] * x$init[n]
            }
          }
          
          dij[, , t + 1] <- apply(partdij[, , , , t + 1], c(1, 2), sum)
          
          for (r in 1:x$s) {
            part21[r, i, j, t + 1] <- sum(psi[r, i, 1:(t + 1)] * Q[i, j, (t + 1):1])
          }
          
        }
      }
      
      part11[, , t + 1] <- (dij[, , t + 1] - indicator * (t(Psi[, , t + 1]) %*% x$init %*% t(rep.int(times = x$s, x = 1)))) ^ 2
      part1[, , t + 1] <- apply(X = part11[, , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum)
      
      if (t > 0) {
        part22[, , t + 1] <- apply(X = dij[, , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum)
      }
      
      part2[, , t + 1] <- part22[, , t + 1] - indicator * apply(part21[, , , t + 1], c(2, 3), function(elt) elt %*% x$init)
      
    }
    
    sigma2 <-
      apply(apply(part1, c(1, 3), sum) - (apply(part2, c(1, 3), sum)) ^ 2, c(2), function(elt)
        sum(elt * mu))
    
    return(list(avail = avail, sigma2 = sigma2))
  
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
#'   The failure rate at time \eqn{k = 0} is defined by \eqn{\lambda(0) := 1 - R(0)},
#'   with \eqn{R} being the reliability function.
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   BMP-Failure Rate should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the BMP-Failure 
#'   Rate for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, 
#'   a list containing the following components:
#'   \itemize{
#'    \item{lbda: }{a vector of length \eqn{k + 1} giving the values of the 
#'      BMP-Failure Rate for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(k)}.}
#'  }
#'
#' @export
#'
failureRateBMP <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
  ###########################################################
  # Compute \lambda, the BMP-failure rate
  ###########################################################
  
  reliab <- reliability(x = x, k = k, upstates = upstates)
  
  lbda <- rep.int(0, k + 1)
  lbda[1] <- 1 - reliab[1] # k = 0
  
  for (j in 2:(k + 1)) {
    lbda[j] <- ifelse(reliab[j - 1] != 0, 1 - reliab[j] / reliab[j - 1], 0)
  }
  
  ###########################################################
  # Compute the variance (Theorem 5.4. equation (5.35), p.119)
  # 
  # Be careful: 
  # 
  # In the formula (5.35), we use q_{Y} (Proposition 5.1 p.105-106) instead 
  # of q, and every others quantities such as \psi, \Psi,\dots derive from q_{Y}
  # 
  # 
  # The decomposition of the variance is as follows:
  # 
  # 
  # 
  # D^{U}_{ij} := \underbrace{\sum_{n \in U} \sum_{r \in U} \underbrace{\alpha_{n} \psi_{ni} * \psi_{jr} * (\text{I} - diag(\text{Q.1}))_{rr}}_{partduij}}_{duij}
  # 
  ###########################################################
  
  if (var) {
    
    alpha1 <- x$init[which(x$states %in% upstates)]
    u <- length(upstates)
    
    q <- .get.qy(x = x, k =  k, upstates = upstates)
    Q <- aperm(apply(q, c(1, 2), cumsum), c(2, 3, 1))
    
    psi <- .get.psi(q = q)
    Psi <- aperm(a = apply(X = psi, MARGIN = c(1, 2), cumsum), perm = c(2, 3, 1))
    
    H <- .get.H(q)
    
    mu <- .get.mu(x = x, klim = klim)
    
    
    indicator <- matrix(data = 1:(u + 1) %in% 1:u, nrow = u + 1, ncol = u + 1) # indicator matrix \mathbb{1}_{i \in U}
    
    convolpsi <- array(data = 0, dim = c(u + 1, u + 1, u + 1, u + 1, k + 1)) # (i, j, n, r, k)
    
    partduij <- array(data = 0, dim = c(u + 1, u + 1, u, u, k + 1)) # # \alpha_{n} \psi_{ni} * \psi_{jr} * (\text{I} - diag(\text{Q.1}))_{rr}
    duij <- array(data = 0, dim = c(u + 1, u + 1, k + 1))
    
    part11 <- array(data = 0, dim = c(u + 1, u + 1, k + 1)) # \left[ D^{U}_{ij} - \mathbb{1}_{i \in U} \sum_{t = 1}^{s} \alpha_{t} \Psi_{ti} \right]^{2}
    part1 <- matrix(data = 0, nrow = u + 1, ncol = k + 1) # \sum_{j = 1}^{s} \left[ D^{U}_{ij} - \mathbb{1}_{i \in U} \sum_{t = 1}^{s} \alpha_{t} \Psi_{ti} \right]^{2} * q_{ij}(k - 1)
    
    part2 <- matrix(data = 0, nrow = u + 1, ncol = k + 1) # \sum_{j = 1}^{s} \left[ D^{U}_{ij} - \mathbb{1}_{i \in U} \sum_{t = 1}^{s} \alpha_{t} \Psi_{ti} \right]^{2} * q_{ij}(k)

    # T
    T <- matrix(data = 0, nrow = u + 1, ncol = k + 1)
    partT11 <- array(data = 0, dim = c(u, u + 1, u + 1, k + 1)) # (t, i, j, k) \psi_{ti} * Q_{ij}(k - 1)
    partT22 <- array(data = 0, dim = c(u, u + 1, u + 1, k + 1)) # (t, i, j, k) \psi_{ti} * Q_{ij}(k)
    partT1 <- array(data = 0, dim = c(u + 1, u + 1, k + 1)) # \mathbb{1}_{i \in U} \sum_{t \in U} \psi_{ti} * Q_{ij}(k - 1)
    partT2 <- array(data = 0, dim = c(u + 1, u + 1, k + 1)) # \mathbb{1}_{i \in U} \sum_{t \in U} \psi_{ti} * Q_{ij}(k)
    
    
    part31 <- array(data = 0, dim = c(u + 1, u + 1, k + 1)) # \left [  \mathbb{1}_{i \in U} D^{U}_{ij} \sum_{t \in U} \alpha_{t} \Psi^{+}_{ti} + \mathbb{1}_{i \in U} (D^{U}_{ij})^{+} \sum_{t \in U} \alpha_{t} \Psi_{ti} - (D^{U}_{ij})^{+}D^{U}_{ij} - \mathbb{1}_{i \in U} \left( \sum_{t \in U} \alpha_{t} \Psi_{ti} \right ) \left( \sum_{t \in U} \alpha_{t} \Psi^{+}_{ti} \right ) \right ]
    part3 <- matrix(data = 0, nrow = u + 1, ncol = k + 1) # \sum_{j = 1}^{s} \left [  \mathbb{1}_{i \in U} D^{U}_{ij} \sum_{t \in U} \alpha_{t} \Psi^{+}_{ti} + \mathbb{1}_{i \in U} (D^{U}_{ij})^{+} \sum_{t \in U} \alpha_{t} \Psi_{ti} - (D^{U}_{ij})^{+}D^{U}_{ij} - \mathbb{1}_{i \in U} \left( \sum_{t \in U} \alpha_{t} \Psi_{ti} \right ) \left( \sum_{t \in U} \alpha_{t} \Psi^{+}_{ti} \right ) \right ] * q_{ij}(k - 1)
    
    sigma2_1 <- rep.int(x = 0, times = k + 1)
    sigma2 <- rep.int(x = 0, times = k + 1)
    
    
    for (t in 0:k) {
      
      for (i in 1:(u + 1)) {
        
        for (j in 1:(u + 1)) {
          
          convolpsi[i, j, , , t + 1] <- outer(X = 1:(1 + u), Y = 1:(1 + u), FUN = function(n, r)
            Reduce('+', lapply(
              X = 0:t,
              FUN = function(h) psi[n, i, h + 1] * psi[j, r, t - h + 1]))
          )
          
          for (r in 1:u) {
            for (n in 1:u) {
              partduij[i, j, n, r, t + 1] <- sum(convolpsi[i, j, n, r, 1:(t + 1)] * (1 - H[r, r, (t + 1):1]))
              partduij[i, j, n, r, t + 1] <- partduij[i, j, n, r, t + 1] * alpha1[n]
            }
          }
          
          for (r in 1:u) {
            if (t >= 1) {
              partT11[r, i, j, t + 1] <- sum(psi[r, i, 1:t] * Q[i, j, t:1])
            }
            partT22[r, i, j, t + 1] <- sum(psi[r, i, 1:(t + 1)] * Q[i, j, (t + 1):1])
          }
        }
      }
      
      duij[, , t + 1] <- apply(partduij[, , , , t + 1], c(1, 2), sum)
      
      part11[, , t + 1] <- (duij[, , t + 1] - indicator * (t(Psi[, , t + 1]) %*% c(alpha1, 0) %*% t(rep.int(times = u + 1, x = 1)))) ^ 2
      
      if (t >= 1) {
        part1[, t + 1] <- apply(apply(X = part11[, , 1:t] * q[, , t:1], MARGIN = c(1, 2), sum), 1, sum)
      }
      
      part2[, t + 1] <- apply(apply(X = part11[, , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum), 1, sum)
      
      # partT1[, , t + 1] <- indicator * apply(partT11[, , , t + 1], c(2, 3), function(elt) elt %*% alpha1)
      # partT2[, , t + 1] <- indicator * apply(partT22[, , , t + 1], c(2, 3), function(elt) elt %*% alpha1)
      partT1[, , t + 1] <- indicator * drop(apply(partT11[, , , t + 1, drop = FALSE], c(2, 3, 4), function(elt) elt %*% alpha1))
      partT2[, , t + 1] <- indicator * drop(apply(partT22[, , , t + 1, drop = FALSE], c(2, 3, 4), function(elt) elt %*% alpha1))
      
      
      
      if (t >= 1) {
        T[, t + 1] <- apply(reliab[t + 1] * apply(X = duij[, , 1:t] * q[, , t:1], MARGIN = c(1, 2), sum) -
                              reliab[t] * apply(X = duij[, , 1:(t + 1)] * q[, , (t + 1):1], MARGIN = c(1, 2), sum) -
                              reliab[t + 1] * partT1[, , t + 1] + reliab[t] * partT2[, , t + 1], 1, sum)
      }
      
      if (t >= 1) {
        
        part31[, , t] <- indicator * duij[, , t] * (t(Psi[, , t + 1]) %*% c(alpha1, 0) %*% t(rep.int(times = u + 1, x = 1))) +
          indicator * duij[, , t + 1] * (t(Psi[, , t]) %*% c(alpha1, 0) %*% t(rep.int(times = u + 1, x = 1))) -
          duij[, , t + 1] * duij[, , t] -
          indicator * (t(Psi[, , t]) %*% c(alpha1, 0) %*% t(rep.int(times = u + 1, x = 1))) * (t(Psi[, , t + 1]) %*% c(alpha1, 0) %*% t(rep.int(times = u + 1, x = 1)))
        
      }
      
      if (t >= 1) {
        part3[ , t + 1] <- apply(apply(X = part31[, , 1:t] * q[, , t:1], MARGIN = c(1, 2), sum), 1, sum)
      }
      
      if (t > 1) {
        sigma2_1[t + 1] <- (reliab[t + 1] ^ 2 * part1[, t + 1] + reliab[t] ^ 2 * part2[, t + 1] - T[, t + 1] ^ 2 + 2 * reliab[t] * reliab[t + 1] * part3[, t + 1])[1:u, drop = FALSE] %*% mu[which(x$states %in% upstates)]
        sigma2[t + 1] <- (1 / reliab[t] ^ 4) * sigma2_1[t + 1]
      }
      
    }
    
    return(list(lbda = lbda, sigma2 = sigma2))
    
  } else {
    
    return(lbda)
    
  }
  
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
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return A vector of length \eqn{k + 1} giving the values of the RG-Failure 
#'   Rate for the period \eqn{[0,\dots,k]} if `var = FALSE`. If `var = TRUE`, 
#'   a list containing the following components:
#'   \itemize{
#'    \item{lbda: }{a vector of length \eqn{k + 1} giving the values of the 
#'      RG-Failure Rate for the period \eqn{[0,\dots,k]};}
#'    \item{sigma2: }{a vector giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(k)}.}
#'  }
#'
#' @export
#'
failureRateRG <- function(x, k, upstates = x$states, var = FALSE, klim = 10000) {
  
  # reliab <- reliability(x = x, k = k, upstates = upstates)
  # 
  # r <- rep.int(0, k + 1)
  # r[1] <- ifelse(reliab[1] != 0, -log(reliab[1]), 0) # k = 0
  # 
  # for (j in 2:(k + 1)) {
  #   r[j] <- ifelse(reliab[j - 1] != 0, -log(reliab[j] / reliab[j - 1]), 0)
  # }
  
  lbda <- failureRateBMP(x = x, k = k, upstates = upstates, var = var, klim = klim)
  
  if (var) {
    
    r <- list(r = -log(1 - lbda$lbda), sigma2 = (1 / (1 - lbda$lbda) ^ 2) * lbda$sigma2)
    
  } else {
    
    r <- -log(1 - lbda)
    
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
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector.
#' @return A vector of length \eqn{\textrm{card}(E)} giving the values of the 
#'   mean sojourn times for each state \eqn{i \in E}.
#' 
#' @export
#'
meanSojournTimes <- function(x, states = x$states, klim = 10000) {
  
  q <- .get.q(x = x, k = klim)
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
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function).
#' @param var Optional. A logical value indicating whether or not the variance 
#'   of the estimator should be computed.
#' @return If `var = FALSE`, a vector of length \eqn{\textrm{card}(U)} giving 
#'   the values of the mean time to failure for each state \eqn{i \in U}. If 
#'   `var = TRUE`, a list containing the following components:
#'   \itemize{
#'    \item{mttf: }{a vector of length \eqn{\textrm{card}(U)} giving the values
#'      of the mean time to failure for each state \eqn{i \in U}.}
#'    \item{sigma2: }{the variance of the estimator for each estimation of the 
#'      mean time to failure;}
#'  }
#' 
#' @export
#'
mttf <- function(x, upstates = x$states, klim = 10000, var = FALSE) {
  
  p11 <- x$ptrans[which(x$states %in% upstates), which(x$states %in% upstates), drop = FALSE]
  
  m1 <- meanSojournTimes(x = x, states = upstates, klim = klim)
  
  if (length(m1) != 1) {
    # mttf <- as.numeric(t(x$init) %*% solve(diag(nrow(p11)) - p11) %*% m1)
    mttf <- as.vector(solve(diag(nrow(p11)) - p11) %*% m1)
  } else {
    # mttf <- as.numeric(t(x$init) %*% solve(diag(nrow(p11)) - p11) * m1)
    mttf <- as.vector(solve(diag(nrow(p11)) - p11) * m1)
  }
  names(mttf) <- upstates
  
  if (var) {
    
    q <- .get.q(x = x, k = klim)
    q11 <- q[which(x$states %in% upstates), which(x$states %in% upstates), 2:(klim + 1), drop = FALSE]
    mu <- .get.mu(x = x, klim = klim)
    
    Qml <- apply(q11, 2, function(elt) rowSums(elt * t(sapply(m1, function(x) 1:klim - x))))
    a <- solve(diag(nrow(p11)) - p11)
    etal <- t(a) %*% m1
    etildal <- p11 %*% etal
    
    hik <- t(apply(X = q, MARGIN = c(1, 3), sum))[2:(klim + 1), which(x$states %in% upstates)]
    sigma2_m1 <- colSums(sapply(m1, function(x) (1:klim - x) ^ 2) * hik)
    
    sigma2 <- as.vector(a %*% a %*% (mu[which(x$states %in% upstates)] * (sigma2_m1 + apply(p11 * sapply(etal, function(x) x - etildal) ^ 2, 1, sum) + 2 * Qml %*% etal)))
    names(sigma2) <- upstates
    
    return(list(mttf = mttf, sigma2 = sigma2))
    
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
#'   [smmparametric][smmparametric] or [smmnonparametric][smmnonparametric]).
#' @param downstates Vector giving the subset of non-operational states \eqn{D}.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function).
#' @param var Optional. A logical value indicating whether or not the variance 
#'   of the estimator should be computed.
#' @return If `var = FALSE`, a vector of length \eqn{\textrm{card}(D)} giving 
#'   the values of the mean time to repair for each state \eqn{i \in D}. If 
#'   `var = TRUE`, a list containing the following components:
#'   \itemize{
#'    \item{mttf: }{a vector of length \eqn{\textrm{card}(D)} giving the values
#'      of the mean time to repair for each state \eqn{i \in D}.}
#'    \item{sigma2: }{the variance of the estimator for each estimation of the 
#'      mean time to failure;}
#'  }
#' 
#' @export
#'
mttr <- function(x, downstates = x$states, klim = 10000, var = FALSE) {
  
  return(mttf(x = x, upstates = downstates, klim = klim, var = var))
  
}
