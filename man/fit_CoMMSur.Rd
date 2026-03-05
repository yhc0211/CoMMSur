\name{fit_CoMMSur}
\alias{fit_CoMMSur}
\title{Fit CoMM-Surv and return mediator-level posterior configuration probabilities}
\usage{
fit_CoMMSur(X, M, time, event, mediator_names = NULL, chains = 4,
  iter = 2000, warmup = floor(iter/2), seed = 1, ...)
}
\arguments{
  \item{X}{Numeric vector of length N (exposure).}
  \item{M}{Numeric matrix of dimension N x K (mediators).}
  \item{time}{Numeric vector of length N (observed time).}
  \item{event}{Integer/logical vector of length N (1 = event, 0 = censored).}
  \item{mediator_names}{Optional character vector of length K. If NULL, colnames(M) are used; otherwise defaults to M1,...,MK.}
  \item{chains}{Number of Markov chains.}
  \item{iter}{Total iterations per chain (including warmup).}
  \item{warmup}{Warmup iterations per chain.}
  \item{seed}{Random seed.}
  \item{...}{Additional arguments passed to \code{rstan::sampling}.}
}
\value{
A list with \code{results} (data.frame with columns \code{mediator}, \code{p_null} (= 1 - p_alt), and posterior mean probabilities for each of the 9 joint components) and \code{fit} (the \code{rstan} fit object).
}
\description{
Fits a Bayesian survival mediation model with multiple mediators and a configuration-based joint spike-and-slab prior, and returns mediator-level posterior probabilities for the null/alternative hypothesis and each of the 9 joint configuration components.
}
\examples{
\dontrun{
# Suppose you have X, M, time, event already prepared:
fit <- fit_CoMMSur(X, M, time, event, iter = 1000, chains = 2)
fit$results
}
}
