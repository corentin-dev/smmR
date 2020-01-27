\name{LoglikelihoodMk}
\alias{LoglikelihoodMk}
%- Also NEED an '\alias' for EACH other topic documented here.
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loglikelihood (Markov model)
}
\description{
Computation of the loglikelihood starting from sequence(s), alphabet, initial distribution, transition matrix}
\usage{
LoglikelihoodMk(seq, E, mu, Ptrans, k)
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
Vector of initial distribution
}
  \item{Ptrans}{
Matrix of transition probabilities
}
  \item{k}{
Order of the Markov chain
}
}

\value{
\item{L}{Value of loglikelihood for each sequence}
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
\link{LoglikelihoodSM}
}
\examples{
alphabet = c("a","c","g","t")
S = length(alphabet)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), 
nrow = S, ncol = S, byrow = TRUE)

#Pij
#     [,1] [,2] [,3] [,4]
#[1,]  0.0  0.2  0.3  0.5
#[2,]  0.4  0.0  0.2  0.4
#[3,]  0.1  0.2  0.0  0.7
#[4,]  0.8  0.1  0.1  0.0


## Simulation of two sequences of length 20 and 50 respectively
seq2 = simulMk(E = alphabet, nbSeq = 2, lengthSeq = c(20,50), Ptrans = Pij, 
init = rep(1/4,4), k = 1)
                
####################################
## Computation of the loglikelihood
####################################
LoglikelihoodMk(seq = seq2, E = alphabet, mu = rep(1/4,4), Ptrans = Pij, k = 1)

#$L
#$L[[1]]
#[1] -13.90161
#
#$L[[2]]
#[1] -39.58438

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Markov models}
\keyword{Loglikelihood}% __ONLY ONE__ keyword per line
