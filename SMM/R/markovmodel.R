markovmodel <- function(E, init, ptrans, k = 1) {

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
  # Checking parameter k
  #############################
  
  if (!((k > 0) && ((k %% 1) == 0))) {
    stop("k must be a strictly positive integer")
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
  
  if (!((dim(ptrans)[1] == S ^ k) && (dim(ptrans)[2] == S))) {
    stop("The size of the matrix ptrans must be equal to SxS")
  }
  
  if (!all(apply(ptrans, 1, sum) == 1)) {
    stop("ptrans is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }


  ans <- list(E = E, init = init, ptrans = ptrans, k = k)
  
  class(ans) <- "markovmodel"
  
  return(ans)
}

is.markovmodel <- function(x) {
  inherits(x, "markovmodel")
}

.stationary.distribution.markovmodel <- function(x) {
  
  m <- dim(x$ptrans)[1] # Number of states
  
  A <- t(x$ptrans) - diag(1, m, m)
  A[m, ] <- 1
  b <- c(rep(0, (m - 1)), 1)
  statlaw <- solve(A, b)
  
  return(statlaw)
}

AIC.fittedmarkovmodel <- function(object, ..., k = 2) {
  
  S <- length(object$estimate$E)
  
  # Kpar: number of parameters of the model
  Kpar <- (S - 1) * S ^ object$estimate$k
  
  vecAIC <- -2 * object$logliks + k * Kpar
  
  return(vecAIC)
  
}

BIC.fittedmarkovmodel <- function(object, ...) {
  
  S <- length(object$estimate$E)
  
  # Kpar: number of parameters of the model
  Kpar <- (S - 1) * S ^ object$estimate$k
  
  n <- sapply(object$seq, length)
  
  vecBIC <- -2 * object$logliks + log(n) * Kpar
  
  return(vecBIC)
}