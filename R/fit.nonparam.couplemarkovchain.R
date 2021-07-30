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
  phat[is.na(phat) | (phat < 0)] <- 0
  
  q <- computeKernelNonParamEndcensoring(phat)
  q <- q[, , -1]
  
  ptrans <- rowSums(q, dims = 2)
  
  # Renormalize ptrans for potential numerical issues
  ptrans <- .normalizePtrans(ptrans)
  
  if (type.sojourn == "fij") {
    
    f <- q / array(ptrans, c(s, s, kmax))
    f[is.na(f) | (f < 0)] <- 0
    f <- f / array(apply(f, c(1, 2), sum), c(s, s, kmax)) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
  } else if (type.sojourn == "fi") {
    
    f <- apply(q, c(1, 3), sum)
    f[is.na(f) | (f < 0)] <- 0
    f <- f / apply(f, 1, sum) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
  } else if (type.sojourn == "fj") {
    
    f <- apply(q, c(2, 3), sum) / apply(ptrans, 2, sum)
    f[is.na(f) | (f < 0)] <- 0
    f <- f / apply(f, 1, sum) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
  } else if (type.sojourn == "f") {
    
    f <- apply(q, 3, sum) / s
    f[is.na(f) | (f < 0)] <- 0
    f <- f / sum(f) # Renormalize f
    f[is.na(f) | (f < 0)] <- 0
    
  }
  
  # Initial distribution
  if (is.vector(init.estim) & length(init.estim) == 1) {
    if (init.estim == "mle") {
      init <- Nstart / sum(Nstart)
    } else if (init.estim == "limit") {
      init <- .limitDistribution(q = q, ptrans = ptrans)
    } else if (init.estim == "freq") {
      init <- counting$Ni / counting$N
    } else if (init.estim == "unif") {
      init <- rep.int(x = 1 / s, times = s)
    } else {
      stop("'init.estim' must be equal to \"mle\", \"limit\", \"freq\" or \"unif\".
           'init.estim' can also be a vector of length s for custom initial distribution")
    }
  } else {
    if (!(is.numeric(init.estim) & !anyNA(init.estim) & is.vector(init.estim) & length(init.estim) == s)) {
      stop("'init.estim' is not a numeric vector of length s")
    }
    
    if (!(all(init.estim >= 0) & all(init.estim <= 1))) {
      stop("Probabilities in 'init.estim' must be between [0, 1]")
    }
    
    if (!((sum(init.estim) >= 1 - sqrt(.Machine$double.eps)) | (sum(init.estim) <= 1 + sqrt(.Machine$double.eps)))) {
      stop("The sum of 'init.estim' is not equal to one")
    }
    
    init <- init.estim
  }
  
  init <- as.vector(init / sum(init))
  
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
