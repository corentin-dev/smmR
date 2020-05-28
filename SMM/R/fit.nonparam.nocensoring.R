.fit.nonparam.nocensoring <- function(seq, type.sojourn = c("fij", "fi", "fj", "f"), cens.beg = cens.beg) {
  
  S <- seq$S
  E <- seq$E
  L <- seq$L
  Kmax <- seq$Kmax
  counting <- seq$counting
  

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
  p <- Nij / tcrossprod(Ni, rep.int(1, S))
  p[which(is.na(p))] <- 0
  
  # # Estimation of the transition matrix
  # a <- Nijk / array(Ni, c(S, S, Kmax))
  # a[which(is.na(a))] <- 0
  
  # Estimation of the sojourn time distribution and the kernel
  if (type.sojourn == "fij") {
    
    f <- Nijk / array(Nij, c(S, S, Kmax))
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(S, S, Kmax)) * f
    
  } else if (type.sojourn == "fi") {
    
    f <- Nik / Ni %*% t(rep.int(1, Kmax))
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(S, S, Kmax)) * aperm(array(f, c(S, Kmax, S)), c(1, 3, 2))
    
  } else if (type.sojourn == "fj") {
    
    f <- Njk / Nj %*% t(rep.int(1, Kmax))
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(S, S, Kmax)) * aperm(array(f, c(S, Kmax, S)), c(3, 1, 2))
    
  } else if (type.sojourn == "f") {
    
    f <- Nk / N
    f[which(is.na(f))] <- 0
    
    q <- array(p, c(S, S, Kmax)) * aperm(array(f, c(Kmax, S, S)), c(2, 3, 1))
    
  }
  
  # Initial distribution
  if (L >= S * 10) {
    init <- Nstart / sum(Nstart)
  } else {# Computation of the limit distribution
    init <- .limitDistribution(q = q, ptrans = p)
  }
  
  estimate <-
    list(
      E = E,
      S = S,
      Kmax = Kmax,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = p,
      laws = f,
      cens.beg = cens.beg,
      cens.end = FALSE
    )
  
  class(estimate) <- c("smm", "smmnonparametric")
  
  return(estimate)
}
