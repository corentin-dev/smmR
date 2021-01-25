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
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in states} of the 
#'      embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Y_1 = i)}, \eqn{i \in 1, 2, \dots, s};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in states,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )},
#'      f is specified by the argument `param` in the parametric case and by 
#'      `laws` in the non-parametric case.
#'  }
#'
#' The estimation of the transition matrix of the embedded Markov chain is 
#' \eqn{\widehat{p_{trans}}(i,j) = \frac{N_{ij}}{N_{i.}}}.
#'
#' The estimation of the initial law is the limit law if the number of sequences 
#' is less than 10xS, with s the length of the state space, else the estimation 
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
#' If  `type.sojourn = "fij"`, `distr` is a matrix of size sxs (e.g., if the 
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
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param distr By default `"nonparametric"` for the non-parametric estimation 
#'   case.
#'   
#'   If the parametric estimation case is desired, `distr` should be:
#'   \itemize{
#'     \item A matrix of distributions of size sxs if `type.sojourn = "fij"`;
#'     \item A vector of distributions of size s if `type.sojourn = "fi"` or `"fj"`;
#'     \item A distribution if `type.sojourn = "f"`.
#'   }
#'   
#'   The distributions to be used in `distr` must be one of `"unif"`, `"geom"`, 
#'   `"pois"`, `"dweibull"`, `"nbinom"`.
#' @param init.estim Optional. Method used to estimate the initial distribution.
#'   If `init.estim = "mle"`, then the classical Maximum Likelihood Estimator 
#'   is used. If `init.estim = "stationary"`, then the initial distribution
#'   is replaced by the stationary distribution of the semi-Markov chains.
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
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' 
#' # Creation of the initial distribution
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' 
#' # Creation of transition matrix
#' pij <- matrix(c(0, 0.2, 0.5, 0.3, 
#'                 0.2, 0, 0.3, 0.5, 
#'                 0.3, 0.5, 0, 0.2, 
#'                 0.4, 0.2, 0.4, 0), 
#'               ncol = s, byrow = TRUE)
#' 
#' # Creation of the distribution matrix
#' 
#' distr.matrix <- matrix(c(NA, "pois", "geom", "nbinom", 
#'                          "geom", NA, "pois", "dweibull",
#'                          "pois", "pois", NA, "geom", 
#'                          "pois", "geom", "geom", NA), 
#'                        nrow = s, ncol = s, byrow = TRUE)
#' 
#' # Creation of an array containing the parameters
#' param1.matrix <- matrix(c(NA, 2, 0.4, 4, 
#'                           0.7, NA, 5, 0.6, 
#'                           2, 3, NA, 0.6, 
#'                           4, 0.3, 0.4, NA), 
#'                         nrow = s, ncol = s, byrow = TRUE)
#' 
#' param2.matrix <- matrix(c(NA, NA, NA, 0.6, 
#'                           NA, NA, NA, 0.8, 
#'                           NA, NA, NA, NA, 
#'                           NA, NA, NA, NA), 
#'                         nrow = s, ncol = s, byrow = TRUE)
#' 
#' param.array <- array(c(param1.matrix, param2.matrix), c(s, s, 2))
#' 
#' # Specify the semi-Markov model
#' smm1 <- smmparametric(states = states, init = vect.init, ptrans = pij, 
#'                       type.sojourn = "fij", distr = distr.matrix, 
#'                       param = param.array)
#' 
#' seq1 <- simulate(object = smm1, nsim = c(1000, 10000, 2000), seed = 100)
#' 
#' # Estimation of simulated sequences
#' est1 <- fitsemimarkovmodel(sequences = seq1, states = states, type.sojourn = "fij", 
#'                            distr = distr.matrix)
#' est1
#' 
fitsemimarkovmodel <-
  function(sequences, states, type.sojourn = c("fij", "fi", "fj", "f"), distr = "nonparametric",
           init.estim = c("mle", "stationary"), cens.beg = FALSE, cens.end = FALSE) {
    

  #############################
  # Checking parameters sequences and states
  #############################

  if (!is.list(sequences)) {
    stop("The parameter sequences should be a list")
  }

  if (!all(unique(unlist(sequences)) %in% states)) {
    stop("Some states in the list of observed sequences sequences are not in the state space states")
  }

  #############################
  # Checking parameter type.sojourn
  #############################

  type.sojourn <- match.arg(type.sojourn)

  #############################
  # Checking parameters distr
  #############################
  
  s <- length(states) # State space size
  
  if (!(length(distr) == 1 && all(distr == "nonparametric"))) {

    if (type.sojourn == "fij" && !(is.matrix(distr))) {
      stop("distr must be a matrix of dimension (s, s) since type.sojourn == \"fij\"")
    }

    if ((type.sojourn == "fi" || type.sojourn == "fj") && !is.vector(distr)) {
      stop("distr must be a vector of length s since type.sojourn == \"fi\" or \"fj\"")
    }

    if (type.sojourn == "f" && !(length(distr) == 1)) {
      stop("distr must be one element since type.sojourn == \"f\"")
    }


    if (type.sojourn == "fij" && !(dim(distr)[1] == s && dim(distr)[2] == s)) {
      stop("distr must be a matrix of dimension (s, s) since type.sojourn == \"fij\"")
    }

    if ((type.sojourn == "fi" || type.sojourn == "fj") && !(length(distr) == s)) {
      stop("distr must be a vector of length s since type.sojourn == \"fi\" or \"fj\"")
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
    }
    
  }
  
  #############################
  # Checking parameter init.estim
  #############################
  
  init.estim <- match.arg(init.estim)
  

  sequences <- processes(sequences = sequences, states = states)

  if (length(distr) == 1 && distr == "nonparametric") {
    
    if (!(cens.beg)) {
      
      if (sequences$L > 1) {# If more than one sequence, use the estimation based on a couple Markov chain
        estimate <- .fit.nonparam.couplemarkovchain(processes = sequences, states = states, type.sojourn = type.sojourn, init.estim = init.estim, cens.end = cens.end)
      } else {
        if (!cens.end) {
          estimate <- .fit.nonparam.nocensoring(processes = sequences, type.sojourn = type.sojourn, init.estim = init.estim, cens.beg = cens.beg)
        } else {
          estimate <- .fit.nonparam.couplemarkovchain(processes = sequences, states = states, type.sojourn = type.sojourn, init.estim = init.estim, cens.end = cens.end)
        }
      }
      
    } else {
     
      if (!cens.end) {
        stop("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = TRUE, cens.end = FALSE")
      } else {
        stop("fitsmm not implemented in the case distr = \"nonparametric\", cens.beg = TRUE, cens.end = TRUE")
      }
       
    }
  } else {
    estimate <- .fit.param(sequences = sequences, states = states, type.sojourn = type.sojourn, distr = distr, init.estim = init.estim, cens.end = cens.end, cens.beg = cens.beg)
  }
  
  if (any(estimate$init == 0)) {
    message("The probabilities of the initial states ",
            paste0(names(which(estimate$init == 0)), collapse = ", "),
            " are 0.")
  }
  
  return(estimate)
}
