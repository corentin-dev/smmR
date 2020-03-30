smmparametric <- function(E, init, ptrans, type.sojourn = c("fij", "fi", "fj", "f"), 
                          distr, param, cens.beg = FALSE, cens.end = FALSE) {
  
  #############################
  # Checking parameter E
  #############################
  
  S <- length(E)
  
  if (!(is.vector(E) && (length(unique(E)) == S))) {
    stop("The state space E is not a vector of unique elements")
  }
  
  #############################
  # Checking parameter init
  #############################
  
  if (!(is.vector(init) && (length(init) == S))) {
    stop("init is not a vector of length S")
  }
  
  if (!(all(init >= 0) && all(init <= 1))) {
    stop("Probabilities in init must be between [0, 1]")
  }
  
  if (!(sum(init) == 1)) {
    stop("The sum of init is not equal to one")
  }
  
  #############################
  # Checking parameter ptrans
  #############################
  
  if (!is.matrix(ptrans)) {
    stop("ptrans is not a matrix")
  }
  
  if (!(all(ptrans >= 0) && all(ptrans <= 1))) {
    stop("Probabilities in ptrans must be between [0, 1]")
  }
  
  if (!((dim(ptrans)[1] == S) && (dim(ptrans)[2] == S))) {
    stop("The size of the matrix ptrans must be equal to SxS")
  }
  
  if (!all(apply(ptrans, 1, sum) == 1)) {
    stop("ptrans is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }
  
  if (!all(diag(ptrans) == 0)) {
    stop("All the diagonal elements of ptrans must be equal to 0 since transitions to the same state are not allowed")
  }
  
  #############################
  # Checking parameter type.sojourn
  #############################
  
  type.sojourn <- match.arg(type.sojourn)
  
  #############################
  # Checking parameters distr and param
  #############################
  
  if (type.sojourn == "fij" && !(is.matrix(distr) && (is.array(param) && !is.matrix(param)))) {
    stop("distr must be a matrix of size SxS and param must be an array of size SxSx2 since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !(is.vector(distr) && is.matrix(param))) {
    stop("distr must be a vector of length S and param must be a matrix of size Sx2 since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !((length(distr) == 1) && is.vector(param))) {
    stop("distr must be one element and param must be a vector of length 2 since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !((dim(distr)[1] == S && dim(distr)[2] == S) && (dim(param)[1] == S && dim(param)[2] == S))) {
    stop("distr must be a matrix of size SxS and param must be an array of size SxSx2 since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !((length(distr) == S) && (dim(param)[1] == S && dim(param)[2] == 2))) {
    stop("distr must be a vector of length S and param must be a matrix of size Sx2 since type.sojourn == \"fi\" or \"fj\"")
  }
  
  
  distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom", NA)
  if (!all(distr %in% distrib.vec)) {
    stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "), 
         ".\n Incorrect distribution(s) found in distr: ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
  }
  
  if (!all(param >= 0 || is.na(param))) {
    stop("Every element of param must be positive")
  }
  
  if (type.sojourn == "fij") {
    if (!(all(is.na(diag(distr))))) {
      stop("All the diagonal elements of distr must be equal to NA since transitions to the same state are not allowed")
    }
    
    if (!(all(is.na(diag(param[, , 1]))) && all(is.na(diag(param[, , 2]))))) {
      stop("All the diagonal elements of param must be equal to NA since transitions to the same state are not allowed")
    }
    
    if (!all(!is.na(distr[row(distr) != col(distr)]))) {
      stop("All the non-diagonal elements of distr must be specified. Found NAs values.")
    }
  }
  
  #############################
  # Checking parameter cens.beg and cens.end
  #############################
  
  if (!(is.logical(cens.beg) && is.logical(cens.end))) {
    stop("cens.beg and cens.end must be TRUE or FALSE")
  }
  
  
  ans <-
    list(
      E = E,
      S = S,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = distr,
      param = param,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(ans) <- c("smm", "smmparametric")
  
  return(ans)
}

is.smmparametric <- function(x) {
  inherits(x, "smmparametric")
}

.get.f.smmparametric <- function(x, Kmax) {
  
  S <- x$S
  
  if (x$type.sojourn == "fij") {
    param1 <- x$param[, , 1]
    param2 <- x$param[, , 2]
    f <- matrix(0, nrow = S * S, ncol = Kmax)
  } else if (x$type.sojourn == "fj") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = S, ncol = Kmax)
  } else if (x$type.sojourn == "fi") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = S, ncol = Kmax)
  } else {
    param1 <- x$param[1]
    param2 <- x$param[2]
    f <- matrix(0, nrow = 1, ncol = Kmax)
  }
  
  if ("dweibull" %in% x$distr) {
    indices <- which(x$distr == "dweibull")
    for (j in indices) {
      f[j, ] <- ddweibull(1:Kmax, q = param1[j], beta = param2[j], zero = FALSE)
    }
  }
  if ("geom" %in% x$distr) {
    indices <- which(x$distr == "geom")
    for (j in indices) {
      f[j, ] <- dgeom(0:(Kmax - 1), prob = param1[j])
    }
  }
  if ("nbinom" %in% x$distr) {
    indices <- which(x$distr == "nbinom")
    for (j in indices) {
      f[j, ] <- dnbinom(0:(Kmax - 1), size = param1[j], mu = param2[j])
    }
  }
  if ("pois" %in% x$distr) {
    indices <- which(x$distr == "pois")
    for (j in indices) {
      f[j, ] <- dpois(0:(Kmax - 1), lambda = param1[j])
    }
  }
  if ("unif" %in% x$distr) {
    indices <- which(x$distr == "unif")
    for (j in indices) {
      f[j, ] <- sapply(1:Kmax, function(k) ifelse(k <= x$param[j], 1 / x$param[j], 0))
    }
  }
  
  if (x$type.sojourn == "fij") {
    fijk <- f
  } else if (x$type.sojourn == "fi") {
    f <- rep(as.vector(t(f)), each = S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fk <- array(as.vector(t(fmat)), c(S, S, Kmax))
    fijk <- apply(X = fk, MARGIN =  c(1, 3), FUN =  t)
  } else if (x$type.sojourn == "fj") {
    f <- rep(as.vector(t(f)), each = S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(S, S, Kmax))
  } else {
    f <- rep(f, each = S * S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(S, S, Kmax))
  }
  
  return(fijk)
  
}

.get.q.smmparametric <- function(x, Kmax) {
  
  S <- x$S
  
  fijk <- .get.f(x, Kmax)
  q <- array(x$ptrans, c(S, S, Kmax)) * fijk
  
  return(q)
  
}

.getKpar.smmparametric <- function(x) {
  
  distr <- x$distr
  
  nbDweibull <- length(which(distr == "dweibull"))
  nbGeom <- length(which(distr == "geom"))
  nbNbinom <- length(which(distr == "nbinom"))
  nbPois <- length(which(distr == "pois"))
  nbUnif <- length(which(distr == "unif"))
  
  Kpar <- 2 * nbDweibull + nbGeom + 2 * nbNbinom + nbPois + nbUnif
  
  return(Kpar)
}