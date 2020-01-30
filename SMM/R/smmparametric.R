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
  
  distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom")
  if (!all(distr %in% distrib.vec)) {
    stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "), 
         ".\n Incorrect distribution(s) found in distr: ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
  }
  
  if (!all(param >= 0)) {
    stop("Every element of param must be positive")
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
      distr = distr,
      param = param,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(ans) <- "smmparametric"
  
  return(ans)
}

is.smmparametric <- function(x) {
  inherits(x, "smmparametric")
}
