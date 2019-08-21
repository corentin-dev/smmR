\name{InitialLawSM}
\alias{InitialLawSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the initial law (semi-Markov model)
}
\description{
For one sequence, estimation of the initial law by computating the limit law of the semi-Markovian chain.

For several sequences, estimation of the initial law by computating the first state frequencies.
}
\usage{
InitialLawSM(E, seq, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{ Vector of state space
}
  \item{seq}{ List of sequence(s)
}
  \item{q}{ Array of size SxSxKmax with the semi-Markov kernel (S = length(E))}
}
\value{
\item{init}{ Vector of the initial distribution}
}

\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
StationaryLaw, estimSM, simulSM, estimMk, simulMk
}
\examples{
seq = list(c("a","c","c","g","t","a","a","a","a","g","c","t","t","t","g"))
res = estimSM(seq = seq, E = c("a","c","g","t"), distr = "NP")
q = res$q
p = res$Ptrans

InitialLawSM(E = c("a","c","g","t"), seq = seq, q = q)
}
\keyword{Semi-Markov models}
