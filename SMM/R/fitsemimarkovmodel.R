fitsemimarkovmodel <- function(seq, E, type.sojourn = c("fij", "fi", "fj", "f"),
                   distr = "nonparametric", cens.end = FALSE, cens.beg = FALSE) {

  #############################
  # Checking parameters seq and E
  #############################

  if (!is.list(seq)) {
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
  
  if (!(length(distr) == 1 && all(distr == "nonparametric"))) {

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

    distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom", NA)
    if (!all(distr %in% distrib.vec)) {
      stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "), 
           ".\n Incorrect distribution(s) found in distr: ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
    }
    
    if (type.sojourn == "fij") {
      if (!(all(is.na(diag(distr))))) {
        stop("All the diagonal elements of distr must be equal to NA since transitions to the same state are not allowed")
      }
      
      if (!all(!is.na(distr[row(distr) != col(distr)]))) {
        stop("All the non-diagonal elements of distr must be specified. Found NAs values.")
      }
    }
    
  }



  if (length(distr) == 1 && distr == "nonparametric") {
    if (!cens.beg && !cens.end) {
      .fit.nonparam.nocensoring(seq = seq, E = E, type.sojourn = type.sojourn)
    } else if (cens.beg && !cens.end) {
      warning("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = TRUE, cens.end = FALSE")
    } else if (!cens.beg && cens.end) {
      .fit.nonparam.endcensoring(seq = seq, E = E, type.sojourn = type.sojourn)
    } else {
      warning("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = FALSE, cens.end = TRUE")
    }
  } else {
    .fit.param(seq = seq, E = E, type.sojourn = type.sojourn, distr = distr, cens.end = cens.end, cens.beg = cens.beg)
  }
}
