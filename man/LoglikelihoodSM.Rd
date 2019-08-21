\name{LoglikelihoodSM}
\alias{LoglikelihoodSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loglikelihood (semi-Markov model)
}
\description{
Computation of the loglikelihood starting from sequence(s), alphabet, initial distribution, transition matrix and type of sojourn times}
\usage{
## parametric case
LoglikelihoodSM(seq, E, mu, Ptrans, distr, param, laws = NULL, TypeSojournTime)
## non-parametric case
LoglikelihoodSM(seq, E, mu, Ptrans, distr, param = NULL, laws, TypeSojournTime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{
List of sequence(s)
}
  \item{E}{
Vector of state space
}
  \item{mu}{
Vector of initial distribution of length S
}
  \item{Ptrans}{
Matrix of transition probabilities of the embedded Markov chain \eqn{J=(J_m)_{m}} of size SxS
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
Character: "fij", "fi", "fj", "f" (for more explanations, see Details)
}
}
\details{
In this package we can choose differents types of sojourn time. Four options are available for the sojourn times:
\itemize{
\item depending on the present state and on the next state ("fij");
\item depending only on the present state ("fi");
\item depending only on the next state ("fj");
\item depending neither on the current, nor on the next state ("f").
}
}
\value{
\item{L}{Value of loglikelihood for each sequence}
\item{Kmax}{The maximal observed sojourn time}
}

\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{simulSM},
\link{estimMk},
\link{simulMk},
\link{estimSM}
}
\examples{
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
## Parametric estimation of a trajectory (of length equal to 5000), 
## where the sojourn times depend neither on the present state nor on the next state.
################################
## Simulation of a sequence of length 5000
seq5000 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 5000, TypeSojournTime = "f", 
                init = c(1/4,1/4,1/4,1/4), Ptrans = Pij, distr = "pois", param = 2)
                
#################################
## Computation of the loglikelihood
#################################
LoglikelihoodSM(seq = seq5000, E = alphabet, mu = rep(1/4,4), Ptrans = Pij, 
distr = "pois", param = 2, TypeSojournTime = "f")

#$L
#$L[[1]]
#[1] -1475.348
#
#
#$Kmax
#[1] 10

#------------------------------#
################################
## Non-parametric simulation of several trajectories (3 trajectories of length 1000, 
## 10 000 and 2000 respectively),
## where the sojourn times depend on the present state and on the next state.
################################
## creation of a matrix corresponding to the conditional sojourn time distributions
lengthSeq3 = c(1000, 10000, 2000)
Kmax = 4 
mat1 = matrix(c(0,0.5,0.4,0.6,0.3,0,0.5,0.4,0.7,0.2,0,0.3,0.4,0.6,0.2,0), 
nrow = S, ncol = S, byrow = TRUE)
mat2 = matrix(c(0,0.2,0.3,0.1,0.2,0,0.2,0.3,0.1,0.4,0,0.3,0.2,0.1,0.3,0), 
nrow = S, ncol = S, byrow = TRUE)
mat3 = matrix(c(0,0.1,0.3,0.1,0.3,0,0.1,0.2,0.1,0.2,0,0.3,0.3,0.3,0.4,0), 
nrow = S, ncol = S, byrow = TRUE)
mat4 = matrix(c(0,0.2,0,0.2,0.2,0,0.2,0.1,0.1,0.2,0,0.1,0.1,0,0.1,0), 
nrow = S, ncol = S, byrow = TRUE)
f <- array(c(mat1,mat2,mat3,mat4), c(S,S,Kmax))
### Simulation of 3 sequences
seqNP3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, 
TypeSojournTime = "fij", init = rep(1/4,4), Ptrans = Pij, laws = f, 
File.out = NULL)

#################################
## Computation of the loglikelihood
#################################
LoglikelihoodSM(seq = seqNP3, E = alphabet, mu = rep(1/4,4), Ptrans = Pij, laws = f, 
TypeSojournTime = "fij")

#$L
#$L[[1]]
#[1] -429.35
#
#$L[[2]]
#[1] -4214.521
#
#$L[[3]]
#[1] -818.6451
#
#
#$Kmax
#[1] 4
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Semi-Markov models}
\keyword{Loglikelihood}% __ONLY ONE__ keyword per line
