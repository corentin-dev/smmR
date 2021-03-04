## __________________________________________________________
## Functions checking the parameters of the specified distributions
## __________________________________________________________

checkParameter <- function(distr, param) {
  
  out <- NULL
  
  if (distr == "unif") {
    
    if (is.na(param[1]) | !(is.na(param[2]))) {
      out <- "For Uniform distributions, only the first parameter must be specified"
    }
    
    if (!((param[1] > 0) & ((param[1] %% 1) == 0))) {
      out <- "For Uniform distributions, the value of the parameter must be a positive integer"
    }
    
  } else if (distr == "geom") {
    
    if (is.na(param[1]) | !(is.na(param[2]))) {
      out <- "For Geometric distributions, only the first parameter must be specified"
    }
    
    if (param[1] <= 0 | param[1] >= 1) {
      out <- "For Geometric distributions, the value of the parameter must be 
           between ]0, 1[ (the parameter is the probability of success)"
    }
    
  } else if (distr == "pois") {
    
    if (is.na(param[1]) | !(is.na(param[2]))) {
      out <- "For Poisson distributions, only the first parameter must be specified"
    }
    
    if (param[1] <= 0) {
      out <- "For Poisson distributions, the parameter must be a positive number"
    }
    
  } else if (distr == "dweibull") {
    
    if (anyNA(param)) {
      out <- "For Discrete Weibull distributions, two parameters must be specified"
    }
    
    if (param[1] <= 0 | param[1] >= 1) {
      out <- "For Discrete Weibull distributions, the value of the first parameter must be between ]0, 1["
    }
    
    if (param[2] <= 0) {
      out <- "For Discrete Weibull distributions, the second parameter must be a positive number"
    }
    
  } else if (distr == "nbinom") {
    
    if (anyNA(param)) {
      out <- "For Negative Binomial distributions, two parameters must be specified"
    }
    
    if (param[1] <= 0) {
      out <- "For Negative Binomial distributions, the first parameter must be a 
           positive number (parameter of overdispersion)"
    }
    
    if (param[2] <= 0 | param[2] >= 1) {
      out <- "For Negative Binomial distributions, the value of the second parameter must 
           be between ]0, 1[ (the parameter is the probability of success)"
    }
    
  }
  
  return(out)
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
## .stationaryDistribution
## __________________________________________________________
.stationaryDistribution <- function(ptrans) {
  
  m <- dim(ptrans)[1] # Number of states
  
  A <- t(ptrans) - diag(1, m, m)
  A[m, ] <- 1
  b <- c(rep(0, (m - 1)), 1)
  statdistr <- solve(A, b)
  
  return(statdistr)
}

## __________________________________________________________
## .limitDistribution
## __________________________________________________________
.limitDistribution <- function(q = q, ptrans = ptrans) {
  
  kmax <- dim(q)[3]
  
  fik <- apply(q, c(1, 3), sum)
  mi <- apply(fik, 1, function(x) sum((1:kmax) * x))
  
  statdistr <- .stationaryDistribution(ptrans)
  
  out <- statdistr * mi / sum(statdistr * mi)
  
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
## Normalize transition matrix
## __________________________________________________________
.normalizePtrans <- function(p) {
  
  p <- p / apply(p, 1, sum)
  
  # s <- nrow(p)
  # ind <- 1:s
  # for (i in 1:s) {
  #   col <- sample(x = ind[-i], size = 1)
  #   p[i, col] <- 1 - sum(p[i, -col])
  # }
  
  return(p)
  
}

## __________________________________________________________
## Functions used to estimate the stationary distribution of 
## a Markov chain higher than 1
## __________________________________________________________
.addZeros <- function(vector, alphaSize, n, index) {
  modulo <- index %% n # Get line modulo
  
  if (modulo == 0) {
    modulo <- n # Set modulo to |A|^order-1 if equal to 0
  }
  
  return(c(rep(0, alphaSize * (modulo - 1)), vector, rep(0, alphaSize * (n - modulo)))) # Add zeros to vector
}

.blockMatrix <- function(ptrans) {
  
  rwnms <- rownames(ptrans) # Get rownames
  alphaSize <- ncol(ptrans) # Get states size |ptrans|
  n <- nrow(ptrans) / alphaSize # Get |ptrans|^order-1
  
  # Overlap states in matrix
  block <- sapply(1:nrow(ptrans), function(i, n, alphaSize, ptrans){
    .addZeros(vector = ptrans[i,], n = n, alphaSize = alphaSize, index = i)
  }, n = n, alphaSize = alphaSize, ptrans = ptrans)
  
  # Set dim names
  rownames(block) <- rwnms
  colnames(block) <- rwnms
  
  return(t(block))
}
