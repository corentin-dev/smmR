.fit.nonparam.endcensoring <- function(seq, E, type.sojourn = c("fij", "fi", "fj", "f"), cens.beg = cens.beg) {
  
  S <- seq$S
  E <- seq$E
  L <- seq$L
  Ym <- seq$Ym
  Um <- seq$Um
  Kmax <- seq$Kmax
  counting <- seq$counting
  Nstart <- counting$Nstarti

  
  # Computation of Niujv
  Niujv <- .getCountingNiujv(Ym, Um, S, Kmax)
  Niu <- apply(Niujv, c(1, 2), sum)
  
  phat <- Niujv / array(Niu, c(S, Kmax, S, Kmax))
  phat[is.na(phat)] <- 0
  
  # Computation of q
  q <- .computeKernelNonParamEndcensoring(phat)
  
  ptrans <- rowSums(q, dims = 2)
  
  if (type.sojourn == "fij") {
    
    f <- q / array(ptrans, c(S, S, Kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- apply(q, c(1, 3), sum)
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fj") {
    
    f <- apply(q, c(2, 3), sum) / apply(ptrans, 2, sum)
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "f") {
    
    f <- apply(q, 3, sum) / S
    f[which(is.na(f))] <- 0
    
  }
  
  # Initial distribution
  if (L >= S * 10) {
    init <- Nstart / sum(Nstart)
  } else {# Computation of the limit distribution
    init <- .limitDistribution(q = q, ptrans = ptrans)
  }
  
  estimate <-
    list(
      E = E,
      S = S,
      Kmax = Kmax,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      laws = f,
      cens.beg = cens.beg,
      cens.end = TRUE
    )
  
  class(estimate) <- c("smm", "smmnonparametric")
  
  return(estimate)
  
}
