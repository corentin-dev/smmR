.fit.nonparam.nocensoring <- function(processes, type.sojourn = c("fij", "fi", "fj", "f"), init.estim = init.estim, cens.beg = cens.beg) {
  
  s <- processes$s
  states <- processes$states
  L <- processes$L
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
  p[which(is.na(p))] <- 0
  
  # Estimation of the sojourn time distribution and the kernel
  if (type.sojourn == "fij") {
    
    f <- Nijk / array(Nij, c(s, s, kmax))
    f[which(is.na(f))] <- 0
    
    q <- Nijk / array(data = Ni, dim = dim(Nijk))
    
  } else if (type.sojourn == "fi") {
    
    f <- Nik / matrix(data = Ni, nrow = s, ncol = kmax)
    f[which(is.na(f))] <- 0
    
    # q <- aperm(a = array(data = Nik, dim = c(s, kmax, s)), perm = c(1, 3, 2)) * array(data = Nij, dim = c(s, s, kmax)) / 
    #   array(data = Ni ^ 2, dim = c(s, s, kmax))
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(s, kmax, s)), c(1, 3, 2))
    
  } else if (type.sojourn == "fj") {
    
    f <- Njk / matrix(data = Nj, nrow = s, ncol = kmax)
    f[which(is.na(f))] <- 0
    
    # q <- aperm(a = array(data = Njk, dim = c(s, kmax, s)), perm = c(3, 1, 2)) * array(data = Nij, dim = c(s, s, kmax)) / 
    #   (aperm(a = array(data = Nj, dim = c(s, s, kmax)), perm = c(2, 1, 3)) * array(data = Ni, dim = c(s, s, kmax)))
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(s, kmax, s)), c(3, 1, 2))
    
  } else if (type.sojourn == "f") {
    
    f <- Nk / N
    f[which(is.na(f))] <- 0
    
    # q <- aperm(a = array(data = Nk, dim = c(kmax, s, s)), perm = c(2, 3, 1)) * array(data = Nij, dim = c(s, s, kmax)) / 
    #   (N * array(data = Ni, dim = c(s, s, kmax)))
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(kmax, s, s)), c(2, 3, 1))
    
  }
  
  # Initial distribution
  if (init.estim == "mle") {
    init <- Nstart / sum(Nstart)
  } else {# init.estim == "stationary"
    init <- .limitDistribution(q = q, ptrans = p)
  }

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
