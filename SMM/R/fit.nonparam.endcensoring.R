.fit.nonparam.endcensoring <- function(seq, E, type.sojourn = c("fij", "fi", "fj", "f")) {
  
  S <- length(E)
  
  # All sequences
  nbseq <- length(seq)
  Kmaxstart <- 0
  Kmax <- 0
  
  J <- list()
  T <- list()
  L <- list()
  
  Y <- list()
  U <- list()
  
  for (m in 1:nbseq) {
    
    processes <- .getprocesses(seq[[m]], E)
    
    J[[m]] <- processes$J
    T[[m]] <- processes$T
    L[[m]] <- T[[m]] - c(1, T[[m]][-length(T[[m]] - 1)]) # Sojourn time
    Kmaxstart <- max(Kmaxstart, max(L[[m]]))
    
    Y[[m]] <- processes$Y
    U[[m]] <- processes$U
    Kmax <- max(Kmax, max(U[[m]]) + 1)
  }
  
  # Get the counts in first position for each state
  res <- .count(J, L, S, Kmaxstart)
  Nstart <- res$Nstarti
  
  # Computation of Niujv
  Niujv <- .count.Niujv(Y, U, S, Kmax)
  Niu <- apply(Niujv, c(1, 2), sum)
  
  phat <- Niujv / array(Niu, c(S, Kmax, S, Kmax))
  phat[is.na(phat)] <- 0
  
  # Computation of qs
  q <- .compute.q(phat)
  
  pij <- rowSums(q, dims = 2)
  
  if (type.sojourn == "fij") {
    
    f <- q / array(pij, c(S, S, Kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- apply(q, c(1, 3), sum)
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fj") {
    
    f <- apply(q, c(2,3),sum)
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "f") {
    
    f <- apply(q, 3, sum) / S
    f[which(is.na(f))] <- 0
    
  }
  
  # Initial distribution
  if (nbseq >= S * 10) {
    init <- Nstart / sum(Nstart)
  } else {# Computation of limit law
    init = .limit.distribution(q = q, ptrans = pij)
  }
  
  estimate <-
    list(
      E = E,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = pij,
      distr = "nonparametric",
      param = NULL,
      laws = f,
      cens.beg = FALSE,
      cens.end = FALSE,
      q = q
    )
  
}
