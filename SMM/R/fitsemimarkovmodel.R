#' Estimation of a semi-Markov chain
#'
#' @description Estimation of a semi-Markov chain starting from one or several
#'   sequences. This estimation can be parametric or non-parametric, 
#'   non-censored, censored at the beginning and/or at the end of the 
#'   sequence, with one or several trajectories. Several parametric 
#'   distributions are considered (Uniform, Geometric, Poisson, Discrete 
#'   Weibull and Negative Binomial).
#'
#' @details This function estimates a semi-Markov model in parametric or 
#' non-parametric case, taking into account the type of sojourn time and the 
#' censoring described in references. The non-parametric estimation concerns 
#' sojourn time distributions defined by the user. For the parametric 
#' estimation, several discrete distributions are considered (see below).
#'
#' The difference between the Markov model and the semi-Markov model concerns 
#' the modelisation of the sojourn time. With a Markov chain, the sojourn time 
#' distribution is modeled by a Geometric distribution (in discrete time). 
#' With a semi-Markov chain, the sojourn time can be any arbitrary distribution.
#' In this package, the available distribution for a semi-Markov model are :
#'  \itemize{
#'    \item Uniform: \eqn{f(x) = \frac{1}{n}} for \eqn{1 \le x \le n};
#'    \item Geometric: \eqn{f(x) = \theta (1-\theta)^{x - 1}} for \eqn{x = 1, 2,\ldots,n}, \eqn{0 < \theta \le 1}, \eqn{\theta} is the probability of success;
#'    \item Poisson: \eqn{f(x) = \frac{\lambda^x exp(-\lambda)}{x!}} for \eqn{x = 0, 1, 2,\ldots,n}, with \eqn{\lambda > 0};
#'    \item Discrete Weibull of type 1: \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}}, \eqn{x = 1, 2,\ldots,n}, 
#'      with \eqn{q}, the first parameter and \eqn{\beta} the second parameter;
#'    \item Negative binomial: \eqn{f(x)=\frac{\Gamma(x+\alpha)}{\Gamma(\alpha) x!} p^{\alpha} (1 - p)^x}, 
#'      for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{\Gamma} is the Gamma function, 
#'      \eqn{\alpha} is the parameter of overdispersion and \eqn{p} is the 
#'      probability of success, \eqn{0 < p < 1};
#'    \item Non-parametric.
#'  }
#'  
#'  
#' We define :
#'  \itemize{
#'    \item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in E} of the 
#'      embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Y_1 = i)}, \eqn{i \in 1, 2, \dots, S};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in E,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )},
#'      f is specified by the argument `param` in the parametric case and by 
#'      `laws` in the non-parametric case.
#'  }
#'
#' The estimation of the transition matrix of the embedded Markov chain is 
#' \eqn{\widehat{p_{trans}}(i,j) = \frac{N_{ij}}{N_{i.}}}.
#'
#' The estimation of the initial law is the limit law if the number of sequences 
#' is less than 10xS, with S the length of the state space, else the estimation 
#' of the inital law is \eqn{\widehat{\mu_i} = \frac{N_i^l}{N^l}}, where 
#' \eqn{N_i^l} is the number of times the state i appears in all the sequences 
#' and \eqn{N^l} is the size of sequences.
#'
#' In the parametric case, the distribution of sojourn time is calculated with 
#' the estimated parameters.
#'
#' Note that \eqn{q_{ij}(k) = p_{trans}(i,j) \ f_{ij}(k)} in the general case 
#' (depending on the present state and on the next state). For particular cases, 
#' we replace \eqn{f_{ij}(k)} by \eqn{f_{i.}(k)} (depending on the present 
#' state \eqn{i}), \eqn{f_{.j}(k)} (depending on the next state \eqn{j}) and 
#' \eqn{f_{..}(k)} (depending neither on the present state nor on the next 
#' state).
#' 
#' In this package we can choose different types of sojourn time. 
#' Four options are available for the sojourn times:
#' \itemize{
#'   \item depending on the present state and on the next state (`fij`);
#'   \item depending only on the present state (`fi`);
#'   \item depending only on the next state (`fj`);
#'   \item depending neither on the current, nor on the next state (`f`).
#' }
#' 
#' If  `type.sojourn = "fij"`, `distr` is a matrix of size SxS (e.g., if the 
#' row 1 of the 2nd column is `"pois"`, that is to say we go from the first 
#' state to the second state following a Poisson distribution).
#' If `type.sojourn = "fi"` or `"fj"`, `distr` must be a vector (e.g., if the 
#' first element of vector is `"geom"`, that is to say we go from the first 
#' state to any state according to a Geometric distribution).
#' If `type.sojourn = "f"`, `distr` must be one of `"unif"`, `"geom"`, `"pois"`, 
#' `"dweibull"`, `"nbinom"` (e.g., if `distr` is equal to `"nbinom"`, that is 
#' to say that the sojourn times when going from any state to any state follows 
#' a Negative Binomial distribution).
#' For the non-parametric case, `distr` is equal to `"nonparametric"` whatever 
#' type of sojourn time given.
#' 
#' If the sequence is censored at the beginning and at the end, `cens.beg` 
#' must be equal to `TRUE` and `cens.end` must be equal to `TRUE` too. 
#' All the sequences must be censored in the same way.
#'
#' @param seq A list of vectors representing the sequences.
#' @param E Vector of state space (of length S).
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param distr By default `"nonparametric"` for the non-parametric estimation 
#'   case.
#'   
#'   If the parametric estimation case is desired, `distr` should be:
#'   \itemize{
#'     \item A matrix of distributions of size SxS if `type.sojourn = "fij"`;
#'     \item A vector of distributions of size S if `type.sojourn = "fi"` or `"fj"`;
#'     \item A distribution if `type.sojourn = "f"`.
#'   }
#'   
#'   The distributions to be used in `distr` must be one of `"unif"`, `"geom"`, 
#'   `"pois"`, `"dweibull"`, `"nbinom"`.
#' @param cens.beg Optional. A logical value indicating whether or not 
#'   sequences are censored at the beginning.
#' @param cens.end Optional. A logical value indicating whether or not 
#'   sequences are censored at the end.
#' @return Returns an object of class smm (a [smmnonparametric][smmnonparametric] 
#'   object if `distr = "nonparametric"`, a [smmparametric][smmparametric] 
#'   otherwise).
#' 
#' 
#' @seealso [smmnonparametric], [smmparametric], [simulate.smmparametric]
#' @export
#'
#' @examples 
#' E <- c("a", "c", "g", "t")
#' S <- length(E)
#' 
#' # Creation of the initial distribution
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' 
#' # Creation of transition matrix
#' pij <- matrix(c(0, 0.2, 0.5, 0.3, 
#'                 0.2, 0, 0.3, 0.5, 
#'                 0.3, 0.5, 0, 0.2, 
#'                 0.4, 0.2, 0.4, 0), 
#'               ncol = S, byrow = TRUE)
#' 
#' # Creation of the distribution matrix
#' 
#' distr.matrix <- matrix(c(NA, "pois", "geom", "nbinom", 
#'                          "geom", NA, "pois", "dweibull",
#'                          "pois", "pois", NA, "geom", 
#'                          "pois", "geom", "geom", NA), 
#'                        nrow = S, ncol = S, byrow = TRUE)
#' 
#' # Creation of an array containing the parameters
#' param1.matrix <- matrix(c(NA, 2, 0.4, 4, 
#'                           0.7, NA, 5, 0.6, 
#'                           2, 3, NA, 0.6, 
#'                           4, 0.3, 0.4, NA), 
#'                         nrow = S, ncol = S, byrow = TRUE)
#' 
#' param2.matrix <- matrix(c(NA, NA, NA, 0.6, 
#'                           NA, NA, NA, 0.8, 
#'                           NA, NA, NA, NA, 
#'                           NA, NA, NA, NA), 
#'                         nrow = S, ncol = S, byrow = TRUE)
#' 
#' param.array <- array(c(param1.matrix, param2.matrix), c(S, S, 2))
#' 
#' # Specify the semi-Markov model
#' smm1 <- smmparametric(E = E, init = vect.init, ptrans = pij, 
#'                       type.sojourn = "fij", distr = distr.matrix, 
#'                       param = param.array)
#' 
#' seq1 <- simulate(object = smm1, nsim = c(1000, 10000, 2000), seed = 100)
#' 
#' # Estimation of simulated sequences
#' est1 <- fitsemimarkovmodel(seq = seq1, E = E, type.sojourn = "fij", 
#'                            distr = distr.matrix)
#' est1
#' 
fitsemimarkovmodel <- function(seq, E, type.sojourn = c("fij", "fi", "fj", "f"),
                   distr = "nonparametric", cens.beg = FALSE, cens.end = FALSE) {

  #############################
  # Checking parameters seq and E
  #############################

  if (!is.list(seq)) {
    stop("The parameter seq should be a list")
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

  seq <- sequences(seq = seq, E = E)

  if (length(distr) == 1 && distr == "nonparametric") {
    if (!cens.beg && !cens.end) {
      .fit.nonparam.nocensoring(seq = seq, type.sojourn = type.sojourn, cens.beg = cens.beg)
    } else if (cens.beg && !cens.end) {
      warning("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = TRUE, cens.end = FALSE")
    } else if (!cens.beg && cens.end) {
      .fit.nonparam.endcensoring(seq = seq, E = E, type.sojourn = type.sojourn, cens.beg = cens.beg)
    } else {
      warning("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = TRUE, cens.end = TRUE")
    }
  } else {
    .fit.param(seq = seq, E = E, type.sojourn = type.sojourn, distr = distr, cens.end = cens.end, cens.beg = cens.beg)
  }
}
