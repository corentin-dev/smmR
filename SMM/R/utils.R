## __________________________________________________________
## Functions to generate one random variable according to a specific distribution
## __________________________________________________________

.runif <- function(b, void) {
  return(sample(x = 1:b, size = 1))
}

.rgeom <- function(prob, void) {
  return(rgeom(n = 1, prob = prob) + 1)
}

.rpois <- function(lambda, void) {
  return(rpois(n = 1, lambda = lambda) + 1)
}

.rnbinom <- function(size, prob) {
  return(rnbinom(n = 1, size = size, prob = prob) + 1)
}

.rdweibull <- function(q, beta) {
  return(rdweibull(n = 1, q = q, beta = beta, zero = FALSE))
}

## __________________________________________________________
## Functions giving the value of the densities
## __________________________________________________________

.dunif <- function(x, b, void) {
  return(sapply(x, function(k) ifelse(k <= b, 1 / b, 0)))
}

.dgeom <- function(x, prob, void) {
  return(dgeom(x = x - 1, prob = prob))
}

.dpois <- function(x, lambda, void) {
  return(dpois(x = x - 1, lambda = lambda))
}

.dnbinom <- function(x, size, prob) {
  return(dnbinom(x = x - 1, size = size, prob = prob))
}

.ddweibull <- function(x, q, beta) {
  return(ddweibull(x = x, q = q, beta = beta, zero = FALSE))
}

## __________________________________________________________
## Functions giving the quantiles of the densities
## __________________________________________________________

.qunif <- function(p, b, void) {
  return(ceiling((b - 1) * p + 1))
}

.qgeom <- function(p, prob, void) {
  return(qgeom(p = p, prob = prob) + 1)
}

.qpois <- function(p, lambda, void) {
  return(qpois(p = p, lambda = lambda) + 1)
}

.qnbinom <- function(p, size, prob) {
  return(qnbinom(p = p, size = size, prob = prob) + 1)
}

.qdweibull <- function(p, q, beta) {
  return(qdweibull(p = p, q = q, beta = beta, zero = FALSE))
}

## __________________________________________________________
## .getSeq: Returns the character sequence given the processes J and T
## __________________________________________________________
.getSeq <- function(Jm, Tm) {
  
  sequences <- c()
  i <- 1
  n <- length(Tm)
  t <- 1
  
  while (i <= n) {
    k <- Tm[i] - t
    sequences <- c(sequences, rep(Jm[i], k))
    t <- Tm[i]
    i <- i + 1
  }
  
  return(sequences)
}

## __________________________________________________________
## .getProcesses: Returns the processes Y, J, T, L and U given the sequence of state
## __________________________________________________________
.getProcesses <- function(sequences, states) {
  
  Ym <- sapply(sequences, function(x) which(states == x), USE.NAMES = FALSE)
  
  n1 <- length(sequences)
  i1 <- 1
  i2 <- 2
  Jm <- c()
  Tm <- c()
  
  while (i1 <= n1) {
    while (i2 <= n1 && (Ym[i2] == Ym[i1])) {
      i2 <- i2 + 1
    }
    Jm <- c(Jm, Ym[i1])
    Tm <- c(Tm, i2)
    i1 <- i2
    i2 <- i1 + 1
  }
  
  n2 <- length(Jm)
  Lm <- Tm - c(1, Tm[-length(Tm)])
  Um <- c()
  
  for (i in 1:n2) {
    for (j in 0:(Lm[i] - 1)) {
      Um <- c(Um, j)
    }
  }
  
  return(list(Ym = Ym, Jm = Jm, Tm = Tm, Lm = Lm, Um = Um))
  
}

## __________________________________________________________
## .getCountingProcesses: Gives the values of the counting processes
## __________________________________________________________
.getCountingProcesses <- function(Jm, Lm, s, kmax) {
  
  L <- length(Jm)
  
  Nstart <- c()
  Nstarti <- rep.int(0, s)
  Nstartil <- matrix(0, s, L)
  Nijkl <- array(0, c(s, s, kmax, L))
  
  Nbijkl <- array(0, c(s, s, kmax, L))
  Neikl <- array(0, c(s, kmax, L))
  
  for (l in 1:L) {
    
    Jl <- Jm[[l]]
    Ll <- Lm[[l]]
    n <- length(Jl)
    
    Nstart <- c(Nstart, Jl[1]) # Initial state
    Nbijkl[Jl[1], Jl[2], Ll[1], l] <- Nbijkl[Jl[1], Jl[2], Ll[1], l] + 1
    Neikl[Jl[n], Ll[n], l] <- Neikl[Jl[n], Ll[n], l] + 1
    
    for (i in 2:n) {
      Nijkl[Jl[i - 1], Jl[i], Ll[i - 1], l] <- Nijkl[Jl[i - 1], Jl[i], Ll[i - 1], l] + 1
    }
  }
  
  Nijk <- apply(Nijkl, c(1, 2, 3), sum)
  
  Nij <- apply(Nijkl, c(1, 2), sum)
  Nijl <- apply(Nijkl, c(1, 2, 4), sum)
  
  indexdiag <- seq(1, s * s, by = s + 1)
  Nij.temp <- as.vector(Nij)[-indexdiag]
  
  statesi <- row(Nij)[-indexdiag][which(Nij.temp == 0)]
  statesj <- col(Nij)[-indexdiag][which(Nij.temp == 0)]
  
  if ((length(statesi) != 0) | (length(statesj) != 0)) {
    message(
      "Some transitions from state i to state j are not observed.
            The following are ",
      paste0(sapply(1:length(statesi), function(x)
        paste0("(i=", statesi[x], " to j=", statesj[x], ")")), collapse = ", "),
      "."
    )
  }
  
  Nbijk <- apply(Nbijkl, c(1, 2, 3), sum)
  Nbjk <- apply(Nbijkl, c(2, 3), sum)
  Nbk <- apply(Nbjk, 2, sum)
  Nbik <- apply(Nbijkl, c(1, 3), sum)
  Neik <- apply(Neikl, c(1, 2), sum)
  Nek <- apply(Neikl, 2, sum)
  
  Nebik <- Nbik + Neik
  Nebk <- Nek + Nbk
  
  Nstarti <- as.vector(count(seq = Nstart, wordsize = 1, alphabet = 1:s))
  Nstartil <- sapply(Nstart, function(x) as.vector(count(seq = x, wordsize = 1, alphabet = 1:s)))
  
  Ni <- apply(Nij, 1, sum)
  
  if (length(which(Ni == 0)) != 0) {
    message("missing state space(s)")
  }
  
  Nj <- apply(Nij, 2, sum)
  
  Nj <- colSums(Nij)
  
  if (length(which(Nj == 0)) != 0) {
    message("missing state space(s)")
  }
  
  N <- sum(Nij)
  
  Nk <- apply(Nijkl, 3, sum)
  Nkl <- apply(Nijkl, c(3, 4), sum)
  
  Nik <- apply(Nijkl, c(1, 3), sum)
  Nikl <- apply(Nijkl, c(1, 3, 4), sum)
  
  Njk <- apply(Nijkl, c(2, 3), sum)
  Njkl <- apply(Nijkl, c(2, 3, 4), sum)
  
  return(
    list(
      Nijkl = Nijkl,
      Nijk = Nijk,
      Nijl = Nijl,
      Nikl = Nikl,
      Njkl = Njkl,
      Nij = Nij,
      Nkl = Nkl,
      Nik = Nik,
      Njk = Njk,
      Ni = Ni,
      Nj = Nj,
      Nk = Nk,
      N = N,
      Nbijkl = Nbijkl,
      Nbijk = Nbijk,
      Nbjk = Nbjk,
      Nbik = Nbik,
      Nbk = Nbk,
      Neik = Neik,
      Nek = Nek,
      Nebik = Nebik,
      Nebk = Nebk,
      Nstarti = Nstarti,
      Nstart = Nstartil
    )
  )
  
}

## __________________________________________________________
## .count.Niujv: Gives the values of the counting processes for the couple (Y, U)
## __________________________________________________________
.getCountingNiujv <- function(Ym, Um, s, kmax) {
  
  L <-  length(Ym)
  Niujv <- array(0, c(s, kmax, s, kmax))
  
  for (k in 1:L) {
    
    Yi <- Ym[[k]]
    Ui <- Um[[k]]
    n <- length(Yi)
    
    for (i in 2:n) {
      Niujv[Yi[i - 1], Ui[i - 1] + 1, Yi[i], Ui[i] + 1] <-
        Niujv[Yi[i - 1], Ui[i - 1] + 1, Yi[i], Ui[i] + 1] + 1
    }
    
  }
  
  return(Niujv)
  
}

## __________________________________________________________
## .computeKernelNonParamEndcensoring
## __________________________________________________________
.computeKernelNonParamEndcensoring <- function(p) {
  
  s <- dim(p)[1]
  kmax <- dim(p)[2]
  q <- array(0, c(s, s, kmax))
  
  for (i in 1:s) {
    for (j in 1:s) {
      for (dk in 1:kmax) {
        
        pi <- p[i, dk, j, 1]
        
        pr <- pi
        if (dk >= 2) {
          for (t in 0:(dk - 2)) {
            pr = pr * p[i, t + 1, i, t + 2]
          }
        }
        
        q[i, j, dk] <- pr
      }
    }
  }
  
  return(q)
  
}

## __________________________________________________________
## .stationaryDistribution
## __________________________________________________________
.stationaryDistribution <- function(ptrans) {
  
  m <- dim(ptrans)[1] # Number of states
  
  A <- t(ptrans) - diag(1, m, m)
  A[m, ] <- 1
  b <- c(rep(0, (m - 1)), 1)
  statlaw <- solve(A, b)
  
  return(statlaw)
}

## __________________________________________________________
## .limitDistribution
## __________________________________________________________
.limitDistribution <- function(q = q, ptrans = ptrans) {
  
  kmax <- dim(q)[3]
  s <- dim(ptrans)[1]
  
  fik <- apply(q, c(1, 3), sum)
  mi <- apply(fik, 1, function(x) sum((1:kmax) * x))
  
  statlaw <- .stationaryDistribution(ptrans)
  
  out <- statlaw * mi / sum(statlaw * mi)
  
  return(out)
}

## __________________________________________________________
## .productProb: Useful to compute the initial distribution when init.estim == "prod"
## __________________________________________________________
.productProb <- function(length = 2, prob) {
  if (length == 1) {
    return(prob)
  } else {
    return(kronecker(prob, .productProb(length - 1, prob)))
  }
}

## __________________________________________________________
## Discrete-time matrix convolution product (definition 3.5 p. 48)
## __________________________________________________________
.matrixConvolve <- function(A, B) {
  
  ###########################################################
  ###########################################################
  # A and B must be of dimension (S, S, k + 1) where:
  #   - S represents the cardinal of the state space E;
  #   - k represents the time horizon;
  # 
  # Return: Matrix C which is the matrix convolution product A * B
  ###########################################################
  ###########################################################
  
  k <- dim(A)[3] - 1 # A[, , 1] represents k = 0
  
  C <- matrix(data = NA, nrow = nrow(A), ncol(A))
  
  C <-
    Reduce('+', lapply(
      X = 0:k,
      FUN = function(l)
        A[, , k - l + 1] %*% B[, , l + 1]
    ))
  
  return(C)
  
}

## __________________________________________________________
## Normalize transition matrix
## __________________________________________________________
.normalizePtrans <- function(p) {
  
  p <- p / apply(p, 1, sum)
  
  s <- nrow(p)
  ind <- 1:s
  for (i in 1:s) {
    col <- sample(x = ind[-i], size = 1)
    p[i, col] <- 1 - sum(p[i, -col])
  }
  
  return(p)
  
}
