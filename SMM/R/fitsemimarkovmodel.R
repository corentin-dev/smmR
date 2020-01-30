fitsemimarkovmodel <- function(seq, E, type.sojourn = c("fij", "fi", "fj", "f"),
                   distr = "nonparametric", cens.end = FALSE, cens.beg = FALSE) {

  #############################
  # Checking parameters seq and E
  #############################

  if (!is.vector(seq)) {
    stop("The parameter seq should be a vector")
  }

  if (!all(unique(unlist(seq)) %in% E)) {
    stop("Some states in the list of observed sequences seq are not in the state space E")
  }

  #############################
  # Checking parameter type.sojourn
  #############################

  type.sojourn <- match.arg(type.sojourn)

  #############################
  # Checking parameters distr
  #############################
  
  S <- length(E) # State space size
  
  if (!(distr == "nonparametric")) {

    if (type.sojourn == "fij" && !(is.matrix(distr))) {
      stop("distr must be a matrix of size SxS since type.sojourn == \"fij\"")
    }

    if ((type.sojourn == "fi" || type.sojourn == "fj") && !is.vector(distr)) {
      stop("distr must be a vector of length S since type.sojourn == \"fi\" or \"fj\"")
    }

    if (type.sojourn == "f" && !(length(distr) == 1)) {
      stop("distr must be one element since type.sojourn == \"f\"")
    }


    if (type.sojourn == "fij" && !(dim(distr)[1] == S && dim(distr)[2] == S)) {
      stop("distr must be a matrix of size SxS since type.sojourn == \"fij\"")
    }

    if ((type.sojourn == "fi" || type.sojourn == "fj") && !(length(distr) == S)) {
      stop("distr must be a vector of length S since type.sojourn == \"fi\" or \"fj\"")
    }

    distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom")
    if (!all(distr %in% distrib.vec)) {
      stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "),
           ".\n Incorrect distribution(s) found in distr: ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
    }
  }



  if (length(distr) == 1 && distr == "nonparametric") {
    if (!cens.beg && !cens.end) {
      .fit.nonparam.nocensoring(seq = seq, E = E, type.sojourn = type.sojourn)
    } else if (cens.beg && !cens.end) {
      # .estimNP_CensureFin(file = file, seq, E, TypeSojournTime = TypeSojournTime)
      warning("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = TRUE, cens.end = FALSE")
    } else if (!cens.beg && cens.end) {
      .fit.nonparam.endcensoring(seq = seq, E = E, type.sojourn = type.sojourn)
    } else {
      warning("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = FALSE, cens.end = TRUE")
    }
  } else {
    # .estim.plusTraj(seq = seq, E, TypeSojournTime = TypeSojournTime, distr = distr, cens.end = cens.end, cens.beg = cens.beg)
    warning("To be implemented...")
  }
}

## __________________________________________________________
## .fit.nonparam.nocensoring(seq, E, type.sojourn)
## __________________________________________________________
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

## __________________________________________________________
## .fit.nonparam.endcensoring
## __________________________________________________________
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
