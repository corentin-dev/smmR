.fit.nonparam.endcensoring <- function(processes, states, type.sojourn = c("fij", "fi", "fj", "f"), cens.beg = cens.beg) {
  
  s <- processes$s
  states <- processes$states
  L <- processes$L
  Ym <- processes$Ym
  Um <- processes$Um
  kmax <- processes$kmax
  counting <- processes$counting
  Nstart <- counting$Nstarti

  
  # Computation of Niujv
  Niujv <- .getCountingNiujv(Ym, Um, s, kmax)
  Niu <- apply(Niujv, c(1, 2), sum)
  
  phat <- Niujv / array(Niu, c(s, kmax, s, kmax))
  phat[is.na(phat)] <- 0
  
  # Computation of q
  q <- .computeKernelNonParamEndcensoring(phat)
  
  ptrans <- rowSums(q, dims = 2)
  
  if (type.sojourn == "fij") {
    
    f <- q / array(ptrans, c(s, s, kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- apply(q, c(1, 3), sum)
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fj") {
    
    f <- apply(q, c(2, 3), sum) / apply(ptrans, 2, sum)
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "f") {
    
    f <- apply(q, 3, sum) / s
    f[which(is.na(f))] <- 0
    
  }
  
  # Initial distribution
  if (L >= s * 10) {
    init <- Nstart / sum(Nstart)
  } else {# Computation of the limit distribution
    init <- .limitDistribution(q = q, ptrans = ptrans)
  }
  
  estimate <-
    list(
      states = states,
      s = s,
      kmax = kmax,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = f,
      cens.beg = cens.beg,
      cens.end = TRUE
    )
  
  class(estimate) <- c("smm", "smmnonparametric")
  
  return(estimate)
  
}
