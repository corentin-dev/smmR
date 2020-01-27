\name{AIC_SM}
\alias{AIC_SM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AIC (semi-Markov model)
}
\description{
AIC
}
\usage{
AIC_SM(seq, E, mu, Ptrans, distr = "NP", param = NULL, laws = NULL, TypeSojournTime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{
List of sequence(s)
}
  \item{E}{
Vector of state space of length S
}
  \item{mu}{
Vector of initial distribution
}
  \item{Ptrans}{
Matrix of transition probabilities of the embedded Markov chain \eqn{J=(J_m)_{m}}
}
    \item{distr}{
- "NP" for nonparametric case, laws have to be used, param is useless

- Matrix of distributions of size SxS if TypeSojournTime is equal to "fij";

- Vector of distributions of size S if TypeSojournTime is equal to "fi" or "fj";

- A distribution if TypeSojournTime is equal to "f".

The distributions to be used in distr must be one of "uniform", "geom", "pois", "dweibull", "nbinom".
}
  \item{param}{
- Useless if distr = "NP"

- Array of distribution parameters of size SxSx2 (2 corresponds to the maximal number of distribution parameters) if TypeSojournTime is equal to "fij";

- Matrix of distribution parameters of size Sx2 if TypeSojournTime is equal to "fi" or "fj";

- Vector of distribution parameters of length 2 if TypeSojournTime is equal to "f". 
}
  \item{laws}{
- Useless if distr \eqn{\neq} "NP"

- Array of size SxSxKmax if TypeSojournTime is equal to "fij";

- Matrix of size SxKmax if TypeSojournTime is equal to "fi" or "fj";

- Vector of length Kmax if the TypeSojournTime is equal to "f". 

Kmax is the maximum length of the sojourn times. 
}
  \item{TypeSojournTime}{
Character: "fij", "fi", "fj", "f" 
}
}
\details{
\eqn{AIC(M) = -2*log{L} + 2*M}, where L is the log-likelihood, M is the number of parameters of the model.
}
\value{
\item{AIC}{ List: value of AIC for each sequence}
}

\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{simulSM},
\link{estimMk},
\link{simulMk},
\link{estimSM},
\link{LoglikelihoodSM},
\link{LoglikelihoodMk},
\link{AIC_Mk}
}
\examples{
alphabet = c("a","c","g","t")
S = length(alphabet)
## creation of the initial distribution
vect.init = c(1/4,1/4,1/4,1/4)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), 
nrow = S, ncol = S, byrow = TRUE)

#Pij
#     [,1] [,2] [,3] [,4]
#[1,]  0.0  0.2  0.3  0.5
#[2,]  0.4  0.0  0.2  0.4
#[3,]  0.1  0.2  0.0  0.7
#[4,]  0.8  0.1  0.1  0.0


################################
## Parametric simulation of several trajectories (3 trajectories of length 1000, 10 000 
## and 2000 respectively)
## with sojourn times depending on the present state and on the next state
## the sojourn times are modelled by different distributions
################################
lengthSeq3 = c(1000, 10000, 2000)
## creation of the distribution matrix
distr.matrix = matrix(c("dweibull", "pois", "geom", "nbinom", "geom", "nbinom", 
"pois", "dweibull", "pois", "pois", "dweibull", "geom", "pois","geom", "geom", 
"nbinom"), nrow = S, ncol = S, byrow = TRUE)
## creation of an array containing the parameters
param1.matrix = matrix(c(0.6,2,0.4,4,0.7,2,5,0.6,2,3,0.6,0.6,4,0.3,0.4,4), 
nrow = S, ncol = S, byrow = TRUE)
param2.matrix = matrix(c(0.8,0,0,2,0,5,0,0.8,0,0,0.8,0,4,0,0,4), 
nrow = S, ncol = S, byrow = TRUE)
param.array = array(c(param1.matrix, param2.matrix), c(S,S,2))
### simulation of 3 sequences
seq3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, 
TypeSojournTime = "fij", init = vect.init, Ptrans = Pij, distr = distr.matrix,
 param = param.array, File.out = NULL)
                
#################################
## Computation of AIC
#################################
AIC_SM(seq = seq3, E = alphabet, mu = rep(1/4,4), Ptrans = Pij, distr = distr.matrix, 
      param = param.array, TypeSojournTime = "fij")

#[[1]]
#[1] 1566.418
#
#[[2]]
#[1] 15683.48
#
#[[3]]
#[1] 3146.728
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{AIC}
\keyword{Semi-Markov models}% __ONLY ONE__ keyword per line
