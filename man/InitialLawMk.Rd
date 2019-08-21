\name{InitialLawMk}
\alias{InitialLawMk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the initial law (Markov model)
}
\description{
For order 1, estimation of initial law by computing the stationary law of the Markov chain.

For order greater than 1, estimation of initial law by computing the state frequencies.
}
\usage{
InitialLawMk(E, seq, Ptrans, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{ Vector of state space
}
  \item{seq}{ List of sequence(s)
}
  \item{Ptrans}{ Matrix of transition probabilities of size (S^(k))xS, with S = length(E) 
}
  \item{k}{ Order of the Markov chain}
}
\value{
\item{init}{ Vector of the initial law}
}



\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
InitialLawMk, estimSM, simulSM, estimMk, simulMk
}
\examples{
seq = list(c("a","c","c","g","t","a","a","a","a","g","c","t","t","t","g"))
res = estimMk(seq = seq, E = c("a","c","g","t"), k = 1)
p = res$Ptrans

InitialLawMk(E = c("a","c","g","t"), seq = seq, Ptrans = p, k = 1)
}

\keyword{Markov models}
