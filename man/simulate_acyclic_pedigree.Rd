% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_acyclic_pedigree.R
\name{simulate_acyclic_pedigree}
\alias{simulate_acyclic_pedigree}
\title{Simulate an acyclic pedigree}
\usage{
simulate_acyclic_pedigree(Nm = 10, Nf = 10, end_time = 4)
}
\arguments{
\item{Nm}{number of male founders.}

\item{Nf}{number of female founders.}

\item{end_time}{the generation level of the last set of offspring to simulate.
The initial founders are in generation 1.}
}
\description{
This keeps track of the connected members of all the individual
and forbids mating between members of the same connected component.
}
