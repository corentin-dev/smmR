.fit.nonparam.nocensoring <- function(seq, type.sojourn = c("fij", "fi", "fj", "f"), cens.beg = cens.beg) {
  
  S <- seq$S
  E <- seq$E
  nbSeq <- seq$nbSeq
  S <- seq$S
  Y <- seq$Y
  J <- seq$J
  T <- seq$T
  L <- seq$L
  U <- seq$U
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
  
  if (length(which(Ni == 0)) != 0) {# Presence of all states
    warning("Missing states")
  }
  
  indexdiag <- seq(1, S * S, by = S + 1)
  Nij.temp <- as.vector(Nij)[-indexdiag]
  if (length(which(Nij.temp == 0)) != 0) {# Presence of all transitions
    warning("All the transitions are not observed")
  }
  
  p <- Nij / tcrossprod(Ni, rep.int(1, S))
  p[which(is.na(p))] <- 0
  
  if (type.sojourn == "fij") {
    
    f <- Nijk / array(Nij, c(S, S, Kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- Nik / Ni %*% t(rep.int(1, Kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fj") {
    
    f <- Njk / Nj %*% t(rep.int(1, Kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "f") {
    
    f <- Nk / N
    f[which(is.na(f))] <- 0
    
  }
  
  # Initial distribution
  if (nbSeq >= S * 10) {
    init <- Nstart / sum(Nstart)
  } else {# Computation of the limit distribution
    q <- Nijk / array(Ni, c(S, S, Kmax))
    q[which(is.na(q))] <- 0
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
