smmnonparametric <- function(E, init, ptrans, type.sojourn = c("fij", "fi", "fj", "f"), 
                             laws, cens.beg = FALSE, cens.end = FALSE) {
  
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
  
  #############################
  # Checking parameter type.sojourn
  #############################
  
  type.sojourn <- match.arg(type.sojourn)
  
  #############################
  # Checking parameter laws
  #############################
  
  if (type.sojourn == "fij" && !(is.array(laws) && !is.matrix(laws))) {
    stop("laws must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !is.matrix(laws)) {
    stop("laws must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !is.vector(laws)) {
    stop("laws must be a vector of length Kmax since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !(dim(laws)[1] == S && dim(laws)[2] == S)) {
    stop("laws must be an array of size SxSxKmax since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !(dim(laws)[1] == S)) {
    stop("laws must be a matrix of size SxKmax since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (!all(laws >= 0) || !all(laws <= 1)) {
    stop("Probabilities in laws must be between [0, 1]")
  }
  
  
  temp <- apply(laws, c(1, 2), sum)
  indexdiag <- seq(1, S * S, by = S + 1)
  if (type.sojourn == "fij" && !(all(diag(temp == 0)) && all(temp[-indexdiag] == 1))) {
    stop("laws is not a stochastic matrix")
  }
  
  if ((type.sojourn == "fi" || type.sojourn == "fj") && !(all(apply(laws, 1, sum) == 1))) {
    stop("laws is not a stochastic matrix")
  }
  
  if (type.sojourn == "f" && !(sum(laws) == 1)) {
    stop("laws is not a stochastic matrix")
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
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      laws = laws,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(ans) <- "smmnonparametric"
  
  return(ans)
}

is.smmnonparametric <- function(x) {
  inherits(x, "smmnonparametric")
}
