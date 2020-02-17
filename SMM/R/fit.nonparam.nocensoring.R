.fit.nonparam.nocensoring <- function(seq, E, type.sojourn = c("fij", "fi", "fj", "f")) {
  
  S <- length(E)
  nbseq <- length(seq)
  
  Kmax <- 0
  J <- list()
  T <- list()
  L <- list()
  
  for (m in 1:nbseq) {
    
    processes <- .getprocesses(seq[[m]], E)
    
    J[[m]] <- processes$J
    T[[m]] <- processes$T
    L[[m]] <- T[[m]] - c(1, T[[m]][-length(T[[m]] - 1)]) # Sojourn time
    Kmax <- max(Kmax, max(L[[m]])) # Maximal sojourn time
  }
  
  # Get the counts
  res <- .count(J, L, S, Kmax)
  
  Nij <- res$Nij
  Ni <- res$Ni
  Nj <- res$Nj
  N <- res$N
  Nk <- res$Nk
  Nik <- res$Nik
  Njk <- res$Njk
  Nijk <- res$Nijk
  Nstart <- res$Nstarti
  
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
  
  q <- Nijk / array(Ni, c(S, S, Kmax))
  q[which(is.na(q))] <- 0
  
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
  
  # Inital distribution
  if (nbseq >= S * 10) {
    init <- Nstart / sum(Nstart)
  } else{# Computation of the limit distribution
    init <- .limit.distribution(q = q, ptrans = p)
  }
  
  estimate <-
    list(
      E = E,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = p,
      distr = "nonparametric",
      param = NULL,
      laws = f,
      cens.beg = FALSE,
      cens.end = FALSE,
      q = q
    )
  
  return(estimate)
}
