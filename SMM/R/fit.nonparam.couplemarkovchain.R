.fit.nonparam.couplemarkovchain <- function(processes, states, type.sojourn = c("fij", "fi", "fj", "f"), init.estim = c("mle", "stationary"), cens.end) {
  
  s <- processes$s
  states <- processes$states
  Ym <- processes$Ym
  Um <- processes$Um
  kmax <- processes$kmax
  counting <- processes$counting
  Nstart <- counting$Nstarti

  
  Ym <- lapply(Ym, function(x) x - 1)
  Niuj <- getCountingNiuj(Ym, Um, s, kmax)
  Niu <- apply(Niuj, c(1, 2), sum)
  
  phat <- Niuj / array(Niu, c(s, kmax, s))
  phat[is.na(phat)] <- 0
  
  q <- computeKernelNonParamEndcensoring(phat)
  q <- q[, , -1]
  
  ptrans <- rowSums(q, dims = 2)
  
  # Renormalize ptrans for potential numerical issues
  ptrans <- .normalizePtrans(ptrans)
  
  if (type.sojourn == "fij") {
    
    f <- q / array(ptrans, c(s, s, kmax))
    f[which(is.na(f))] <- 0
    f <- f / array(apply(f, c(1, 2), sum), c(s, s, kmax)) # Renormalize f
    f[which(is.na(f))] <- 0
    f[, , dim(f)[3]] <- 1 - apply(f[, , -dim(f)[3]], c(1, 2), sum) # Renormalize f
    diag(f[, , dim(f)[3]]) <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- apply(q, c(1, 3), sum)
    f[which(is.na(f))] <- 0
    f <- f / apply(f, 1, sum) # Renormalize f
    f[, ncol(f)] <- 1 - apply(f[, -ncol(f)], 1, sum) # Renormalize f
    
  } else if (type.sojourn == "fj") {
    
    f <- apply(q, c(2, 3), sum) / apply(ptrans, 2, sum)
    f[which(is.na(f))] <- 0
    f <- f / apply(f, 1, sum) # Renormalize f
    f[, ncol(f)] <- 1 - apply(f[, -ncol(f)], 1, sum) # Renormalize f
    
  } else if (type.sojourn == "f") {
    
    f <- apply(q, 3, sum) / s
    f[which(is.na(f))] <- 0
    f <- f / sum(f) # Renormalize f
    f[length(f)] <- 1 - sum(f[-length(f)]) # Renormalize f
    
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
