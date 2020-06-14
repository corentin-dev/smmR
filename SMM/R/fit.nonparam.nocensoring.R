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
  p <- Nij / tcrossprod(Ni, rep.int(1, s))
  p[which(is.na(p))] <- 0
  
  # # Estimation of the transition matrix
  # a <- Nijk / array(Ni, c(s, s, kmax))
  # a[which(is.na(a))] <- 0
  
  # Estimation of the sojourn time distribution and the kernel
  if (type.sojourn == "fij") {
    
    f <- Nijk / array(Nij, c(s, s, kmax))
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(s, s, kmax)) * f
    
  } else if (type.sojourn == "fi") {
    
    f <- Nik / Ni %*% t(rep.int(1, kmax))
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(s, kmax, s)), c(1, 3, 2))
    
  } else if (type.sojourn == "fj") {
    
    f <- Njk / Nj %*% t(rep.int(1, kmax))
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(s, s, kmax)) * aperm(array(f, c(s, kmax, s)), c(3, 1, 2))
    
  } else if (type.sojourn == "f") {
    
    f <- Nk / N
    f[which(is.na(f))] <- 0
    
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
