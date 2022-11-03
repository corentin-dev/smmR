.fit.nonparam.nocensoring <- function(processes, type.sojourn = c("fij", "fi", "fj", "f"), init.estim = init.estim, cens.beg = cens.beg) {
  
  s <- processes$s
  states <- processes$states
  kmax <- processes$kmax
  counting <- processes$counting
  
  
  Nij <- counting$Nij
  Ni <- counting$Ni
  Nj <- counting$Nj
  N <- counting$N
  Nk <- counting$Nk
  Nik <- counting$Nik
  Njk <- counting$Njk
  Nijk <- counting$Nijk
  Nstart <- counting$Nstarti
  
  # Estimation of the transition matrix
  p <- Nij / matrix(data = Ni, nrow = s, ncol = s)
  p[is.na(p) | (p < 0)] <- 0
  
  # Renormalize p for potential numerical issues
  p <- .normalizePtrans(p)
  
  # Estimation of the sojourn time distribution and the kernel
  if (type.sojourn == "fij") {
    
    f <- Nijk / array(Nij, c(s, s, kmax))
    f[is.na(f) | (f < 0)] <- 0
    f <- f / array(apply(f, c(1, 2), sum), c(s, s, kmax)) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
    q <- Nijk / array(data = Ni, dim = dim(Nijk))
    
  } else if (type.sojourn == "fi") {
    
    f <- Nik / matrix(data = Ni, nrow = s, ncol = kmax)
    f[is.na(f) | (f < 0)] <- 0
    f <- f / apply(f, 1, sum) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
    # q <- aperm(a = array(data = Nik, dim = c(s, kmax, s)), perm = c(1, 3, 2)) * array(data = Nij, dim = c(s, s, kmax)) / 
    #   array(data = Ni ^ 2, dim = c(s, s, kmax))
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(s, kmax, s)), c(1, 3, 2))
    
  } else if (type.sojourn == "fj") {
    
    f <- Njk / matrix(data = Nj, nrow = s, ncol = kmax)
    f[is.na(f) | (f < 0)] <- 0
    f <- f / apply(f, 1, sum) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
    # q <- aperm(a = array(data = Njk, dim = c(s, kmax, s)), perm = c(3, 1, 2)) * array(data = Nij, dim = c(s, s, kmax)) / 
    #   (aperm(a = array(data = Nj, dim = c(s, s, kmax)), perm = c(2, 1, 3)) * array(data = Ni, dim = c(s, s, kmax)))
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(s, kmax, s)), c(3, 1, 2))
    
  } else if (type.sojourn == "f") {
    
    f <- Nk / N
    f[is.na(f) | (f < 0)] <- 0
    f <- f / sum(f) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
    # q <- aperm(a = array(data = Nk, dim = c(kmax, s, s)), perm = c(2, 3, 1)) * array(data = Nij, dim = c(s, s, kmax)) / 
    #   (N * array(data = Ni, dim = c(s, s, kmax)))
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(kmax, s, s)), c(2, 3, 1))
    
  }
  
  # Initial distribution
  if (is.vector(init.estim) & length(init.estim) == 1) {
    if (init.estim == "mle") {
      init <- Nstart / sum(Nstart)
    } else if (init.estim == "limit") {
      init <- .limitDistribution(q = q, ptrans = p)
    } else if (init.estim == "freq") {
      init <- counting$Ni / counting$N
    } else if (init.estim == "unif") {
      init <- rep.int(x = 1 / s, times = s)
    } else {
      stop("'init.estim' must be equal to \"mle\", \"limit\", \"freq\" or \"unif\".
           'init.estim' can also be a vector of length s for custom initial distribution")
    }
  } else {
    if (!(is.numeric(init.estim) & !anyNA(init.estim) & is.vector(init.estim) & length(init.estim) == s)) {
      stop("'init.estim' is not a numeric vector of length s")
    }
    
    if (!(all(init.estim >= 0) & all(init.estim <= 1))) {
      stop("Probabilities in 'init.estim' must be between [0, 1]")
    }
    
    if (!((sum(init.estim) >= 1 - sqrt(.Machine$double.eps)) | (sum(init.estim) <= 1 + sqrt(.Machine$double.eps)))) {
      stop("The sum of 'init.estim' is not equal to one")
    }
    
    init <- init.estim
  }
  
  init <- as.vector(init / sum(init))
  
  estimate <-
    smmnonparametric(
      states = states,
      init = init,
      ptrans = p,
      type.sojourn = type.sojourn,
      distr = f,
      cens.beg = cens.beg,
      cens.end = FALSE
    )
  
  return(estimate)
}
