## __________________________________________________________
## Functions to generate one random variable according to a specific distribution
## __________________________________________________________

.runif <- function(b, void) {
  return(sample(x = 1:b, size = 1))
}

.rgeom <- function(prob, void) {
  return(rgeom(n = 1, prob = prob) + 1)
}

.rpois <- function(lambda, void) {
  return(rpois(n = 1, lambda = lambda) + 1)
}

.rnbinom <- function(size, prob) {
  return(rnbinom(n = 1, size = size, prob = prob) + 1)
}

.rdweibull <- function(q, beta) {
  return(rdweibull(n = 1, q = q, beta = beta, zero = FALSE))
}

## __________________________________________________________
## Functions giving the value of the densities
## __________________________________________________________

.dunif <- function(x, b, void) {
  return(sapply(x, function(k) ifelse(k <= b, 1 / b, 0)))
}

.dgeom <- function(x, prob, void) {
  return(dgeom(x = x - 1, prob = prob))
}

.dpois <- function(x, lambda, void) {
  return(dpois(x = x - 1, lambda = lambda))
}

.dnbinom <- function(x, size, prob) {
  return(dnbinom(x = x - 1, size = size, prob = prob))
}

.ddweibull <- function(x, q, beta) {
  return(ddweibull(x = x, q = q, beta = beta, zero = FALSE))
}

## __________________________________________________________
## Functions giving the quantiles of the densities
## __________________________________________________________

.qunif <- function(p, b, void) {
  return(ceiling((b - 1) * p + 1))
}

.qgeom <- function(p, prob, void) {
  return(qgeom(p = p, prob = prob) + 1)
}

.qpois <- function(p, lambda, void) {
  return(qpois(p = p, lambda = lambda) + 1)
}

.qnbinom <- function(p, size, prob) {
  return(qnbinom(p = p, size = size, prob = prob) + 1)
}

.qdweibull <- function(p, q, beta) {
  return(qdweibull(p = p, q = q, beta = beta, zero = FALSE))
}

## __________________________________________________________
## .getSeq: Returns the character sequence given the processes J and T
## __________________________________________________________
.getSeq <- function(Jm, Tm) {
  
  sequences <- c()
  i <- 1
  n <- length(Tm)
  t <- 1
  
  while (i <= n) {
    k <- Tm[i] - t
    sequences <- c(sequences, rep(Jm[i], k))
    t <- Tm[i]
    i <- i + 1
  }
  
  return(sequences)
}

## __________________________________________________________
## .stationaryDistribution
## __________________________________________________________
.stationaryDistribution <- function(ptrans) {
  
  m <- dim(ptrans)[1] # Number of states
  
  A <- t(ptrans) - diag(1, m, m)
  A[m, ] <- 1
  b <- c(rep(0, (m - 1)), 1)
  statlaw <- solve(A, b)
  
  return(statlaw)
}

## __________________________________________________________
## .limitDistribution
## __________________________________________________________
.limitDistribution <- function(q = q, ptrans = ptrans) {
  
  kmax <- dim(q)[3]
  s <- dim(ptrans)[1]
  
  fik <- apply(q, c(1, 3), sum)
  mi <- apply(fik, 1, function(x) sum((1:kmax) * x))
  
  statlaw <- .stationaryDistribution(ptrans)
  
  out <- statlaw * mi / sum(statlaw * mi)
  
  return(out)
}

## __________________________________________________________
## .productProb: Useful to compute the initial distribution when init.estim == "prod"
## __________________________________________________________
.productProb <- function(length = 2, prob) {
  if (length == 1) {
    return(prob)
  } else {
    return(kronecker(prob, .productProb(length - 1, prob)))
  }
}

## __________________________________________________________
## Discrete-time matrix convolution product (definition 3.5 p. 48)
## __________________________________________________________
.matrixConvolve <- function(A, B) {
  
  ###########################################################
  ###########################################################
  # A and B must be of dimension (S, S, k + 1) where:
  #   - S represents the cardinal of the state space E;
  #   - k represents the time horizon;
  # 
  # Return: Matrix C which is the matrix convolution product A * B
  ###########################################################
  ###########################################################
  
  k <- dim(A)[3] - 1 # A[, , 1] represents k = 0
  
  C <- matrix(data = NA, nrow = nrow(A), ncol(A))
  
  C <-
    Reduce('+', lapply(
      X = 0:k,
      FUN = function(l)
        A[, , k - l + 1] %*% B[, , l + 1]
    ))
  
  return(C)
  
}

## __________________________________________________________
## Normalize transition matrix
## __________________________________________________________
.normalizePtrans <- function(p) {
  
  p <- p / apply(p, 1, sum)
  
  s <- nrow(p)
  ind <- 1:s
  for (i in 1:s) {
    col <- sample(x = ind[-i], size = 1)
    p[i, col] <- 1 - sum(p[i, -col])
  }
  
  return(p)
  
}
