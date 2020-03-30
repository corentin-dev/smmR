## __________________________________________________________
## .getSeq: Returns the character sequence given the processes J and T
## __________________________________________________________
.getSeq <- function(J, T) {
  
  seq <- c()
  i <- 1
  n <- length(T)
  t <- 1
  
  while (i <= n) {
    k <- T[i] - t
    seq <- c(seq, rep(J[i], k))
    t <- T[i]
    i <- i + 1
  }
  
  return(seq)
}

## __________________________________________________________
## .getProcesses: Returns the processes Y, J, T, L and U given the sequence of state
## __________________________________________________________
.getProcesses <- function(seq, E) {
  
  Y <- sapply(seq, function(x) which(E == x), USE.NAMES = FALSE)
  
  n1 <- length(seq)
  i1 <- 1
  i2 <- 2
  J <- c()
  T <- c()
  
  while (i1 <= n1) {
    while (i2 <= n1 && (Y[i2] == Y[i1])) {
      i2 <- i2 + 1
    }
    J <- c(J, Y[i1])
    T <- c(T, i2)
    i1 <- i2
    i2 <- i1 + 1
  }
  
  n2 <- length(J)
  L <- T - c(1, T[-length(T)])
  U <- c()
  
  for (i in 1:n2) {
    for (j in 0:(L[i] - 1)) {
      U <- c(U, j)
    }
  }
  
  return(list(Y = Y, J = J, T = T, L = L, U = U))
  
}

## __________________________________________________________
## .getCountingProcesses: Gives the values of the counting processes
## __________________________________________________________
.getCountingProcesses <- function(J, L, S, Kmax) {
  
  nbseq <- length(J)
  
  Nstart <- c()
  Nstarti <- rep.int(0, S)
  Nstartil <- matrix(0, S, nbseq)
  Nijkl <- array(0, c(S, S, Kmax, nbseq))
  
  Nbijkl <- array(0, c(S, S, Kmax, nbseq))
  Neikl <- array(0, c(S, Kmax, nbseq))
  
  for (l in 1:nbseq) {
    
    Jl <- J[[l]]
    Ll <- L[[l]]
    M <- length(Jl)
    
    Nstart <- c(Nstart, Jl[1]) # Initial state
    Nbijkl[Jl[1], Jl[2], Ll[1], l] <- Nbijkl[Jl[1], Jl[2], Ll[1], l] + 1
    Neikl[Jl[M], Ll[M], l] <- Neikl[Jl[M], Ll[M], l] + 1
    
    for (i in 2:M) {
      Nijkl[Jl[i - 1], Jl[i], Ll[i - 1], l] <- Nijkl[Jl[i - 1], Jl[i], Ll[i - 1], l] + 1
    }
  }
  
  Nijk <- apply(Nijkl, c(1, 2, 3), sum)
  
  Nij <- apply(Nijkl, c(1, 2), sum)
  Nijl <- apply(Nijkl, c(1, 2, 4), sum)
  
  indexdiag <- seq(1, S * S, by = S + 1)
  Nij.temp <- as.vector(Nij)[-indexdiag]
  if (length(which(Nij.temp == 0)) != 0) {
    warning("Warning : missing transitions")
  }
  
  Nbijk <- apply(Nbijkl, c(1, 2, 3), sum)
  Nbjk <- apply(Nbijkl, c(2, 3), sum)
  Nbk <- apply(Nbjk, 2, sum)
  Nbik <- apply(Nbijkl, c(1, 3), sum)
  Neik <- apply(Neikl, c(1, 2), sum)
  Nek <- apply(Neikl, 2, sum)
  
  Nebik <- Nbik + Neik
  Nebk <- Nek + Nbk
  
  Nstarti <- as.vector(count(seq = Nstart, wordsize = 1, alphabet = 1:S))
  Nstartil <- sapply(Nstart, function(x) as.vector(count(seq = x, wordsize = 1, alphabet = 1:S)))
  
  Ni <- apply(Nij, 1, sum)
  
  if (length(which(Ni == 0)) != 0) {
    warning("Warning : missing state space(s)")
  }
  
  Nj <- apply(Nij, 2, sum)
  
  Nj <- colSums(Nij)
  if (length(which(Nj == 0)) != 0) {
    warning("Warning : missing state space(s)")
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
.getCountingNiujv <- function(Y, U, S, Kmax) {
  
  nbseq <-  length(Y)
  Niujv <- array(0, c(S, Kmax, S, Kmax))
  
  for (k in 1:nbseq) {
    
    Yi <- Y[[k]]
    Ui <- U[[k]]
    M <- length(Yi)
    
    for (i in 2:M) {
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
  
  S <- dim(p)[1]
  Kmax <- dim(p)[2]
  q <- array(0, c(S, S, Kmax))
  
  for (i in 1:S) {
    for (j in 1:S) {
      for (dk in 1:Kmax) {
        
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
  
  Kmax <- dim(q)[3]
  S <- dim(ptrans)[1]
  
  fik <- apply(q, c(1, 3), sum)
  mi <- apply(fik, 1, function(x) sum((1:Kmax) * x))
  
  statlaw <- .stationaryDistribution(ptrans)
  
  out <- statlaw * mi / sum(statlaw * mi)
  
  return(out)
}
