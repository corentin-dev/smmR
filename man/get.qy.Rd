% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgenerics.R
\name{get.qy}
\alias{get.qy}
\title{Method to get the semi-Markov kernel \eqn{q_{Y}}}
\usage{
get.qy(x, k, upstates = x$states)
}
\arguments{
\item{x}{An object of S3 class \code{smmfit} or \code{smm}.}

\item{k}{A positive integer giving the time horizon.}

\item{upstates}{Vector giving the subset of operational states \eqn{U}.}
}
\value{
An array giving the value of \eqn{q_{Y}(k)} at each time between 0
and \code{k}.
}
\description{
Computes the semi-Markov kernel \eqn{q_{Y}(k)}
(See proposition 5.1 p.106).
}
