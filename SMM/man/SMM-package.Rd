\name{SMM-package}
\alias{SMM-package}
\alias{SMM}
\docType{package}
\title{
SMM : Semi-Markov and Markov Models
}
\description{
This package performs parametric and non-parametric estimation and simulation for multi-state discrete-time semi-Markov processes. For the parametric estimation, several discrete distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial. The non-parametric estimation concerns the sojourn time distributions, where no assumptions are done on the shape of distributions. Moreover, the estimation can be done on the basis of one or several sample paths, with or without censoring at the beginning or/and at the end of the sample paths. Estimation and simulation of discrete-time k-th order Markov chains are also considered. 
}
\details{
This R package provides the following functions: estimSM, simulSM, LoglikelihoodSM, AIC_SM, BIC_SM, InitialLawSM, estimMk, simulMk, LoglikelihoodMK, AIC_Mk, BIC_Mk and InitialLawMk.
}

\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


\keyword{Markov models}
\keyword{Semi-Markov models }
\keyword{Simulation}
\keyword{Estimation}
\keyword{Censored data}

\seealso{
\code{\link[=simulMk]{simulMk}}
\code{\link[=estimMk]{estimMk}}
\code{\link[=simulSM]{simulSM}}
\code{\link[=estimSM]{estimSM}}
}
\examples{
##--------------------------------------##
#  An example for the simulMk function   #
#  see more examples >help(simulMk)      #
##--------------------------------------##
### Example 1 ###
# Second order model with the finite state space {a,c,g,t}
E <- c("a","c","g","t")
S = length(E)
init.distribution <- c(1/6,1/6,5/12,3/12)
k<-2
p <- matrix(0.25, nrow = S^k, ncol = S)

# We simulate 3 sequences of size 1000, 10000 and 2000 respectively
simulMk(E = E, nbSeq = 3, lengthSeq = c(1000, 10000, 2000), Ptrans = p, 
init = init.distribution, k = k)

##--------------------------------------##
#  An example for the estimMk function   #
#  see more examples >help(estimMk)      #
##--------------------------------------##
### Example 1 ###
# Second order model with the finite state space {a,c,g,t}
E <- c("a","c","g","t")
S = length(E)
init.distribution <- c(1/6,1/6,5/12,3/12)
k<-2
p <- matrix(0.25, nrow = S^k, ncol = S)

## simulation of 3 sequences with the simulMk function
seq3 = simulMk(E = E, nbSeq = 3, lengthSeq = c(1000, 10000, 2000), Ptrans = p, 
init = init.distribution, k = 2)

## estimation of simulated sequences 
res3 = estimMk(seq = seq3, E = E, k = 2)

##--------------------------------------##
#  An example for the simulSM function   #
#  see more examples >help(simulSM)      #
##--------------------------------------##
alphabet = c("a","c","g","t")
S = length(alphabet)
## creation of the initial distribution
vect.init = c(1/4,1/4,1/4,1/4)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), 
nrow = S, ncol = S, byrow = TRUE)

################################
## Parametric simulation of a trajectory (of length equal to 50) 
## where the sojourn times do not depend on the present state and on the next state
################################
## Simulation of a sequence of length 50
seq50 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 50, TypeSojournTime = "f", 
                init = vect.init, Ptrans = Pij, distr = "pois", param = 2)
                
##--------------------------------------##
#  An example for the simulSM function   #
#  see more examples >help(simulSM)      #
##--------------------------------------##
alphabet = c("a","c","g","t")
S = length(alphabet)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), 
nrow = S, ncol = S, byrow = TRUE)

Pij
#     [,1] [,2] [,3] [,4]
#[1,]  0.0  0.2  0.3  0.5
#[2,]  0.4  0.0  0.2  0.4
#[3,]  0.1  0.2  0.0  0.7
#[4,]  0.8  0.1  0.1  0.0

################################
## Parametric estimation of a trajectory (of length equal to 5000) 
## where the sojourn times do not depend on the present state and on the next state
################################
## Simulation of a sequence of length 5000
seq5000 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 5000, TypeSojournTime = "f", 
                init = c(1/4,1/4,1/4,1/4), Ptrans = Pij, distr = "pois", param = 2)
                
                
## Estimation of the simulated sequence                
estSeq5000 = estimSM(seq = seq5000, E = alphabet, TypeSojournTime = "f", 
                   distr = "pois", cens.end = 0, cens.beg = 0)
                   
}