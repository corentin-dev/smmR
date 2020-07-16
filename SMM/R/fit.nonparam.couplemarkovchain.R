.fit.nonparam.couplemarkovchain <- function(processes, states, type.sojourn = c("fij", "fi", "fj", "f"), init.estim = c("mle", "stationary"), cens.end) {
  
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
  
  # Renormalize ptrans for potential numerical issues
  ptrans <- ptrans / apply(ptrans, 1, sum)
  ptrans <- ptrans / apply(ptrans, 1, sum)
  
  if (type.sojourn == "fij") {
    
    f <- q / array(ptrans, c(s, s, kmax))
    f <- f / array(apply(f, c(1, 2), sum), c(s, s, kmax))
    f <- f / array(apply(f, c(1, 2), sum), c(s, s, kmax))
    f[which(is.na(f))] <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- apply(q, c(1, 3), sum)
    f[which(is.na(f))] <- 0
    f <- f / apply(f, 1, sum)
    f <- f / apply(f, 1, sum)
    
  } else if (type.sojourn == "fj") {
    
    f <- apply(q, c(2, 3), sum) / apply(ptrans, 2, sum)
    f[which(is.na(f))] <- 0
    f <- f / apply(f, 1, sum)
    f <- f / apply(f, 1, sum)
    
  } else if (type.sojourn == "f") {
    
    f <- apply(q, 3, sum) / s
    f[which(is.na(f))] <- 0
    f <- f / sum(f)
    f <- f / sum(f)
    
  }
  
  # Initial distribution
  if (init.estim == "mle") {
    init <- Nstart / sum(Nstart)
  } else {# init.estim == "stationary"
    init <- .limitDistribution(q = q, ptrans = ptrans)
  }
  
  estimate <-
    smmnonparametric(
      states = states,
      init = init,
      ptrans = ptrans,
      type.sojourn = type.sojourn,
      distr = f,
      cens.beg = FALSE,
      cens.end = cens.end
    )
  
  return(estimate)
  
}
