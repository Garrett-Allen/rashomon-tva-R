% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate.R
\name{initialize_sigma}
\alias{initialize_sigma}
\title{Initialize sigma matrix}
\usage{
initialize_sigma(M, R)
}
\arguments{
\item{M}{the number of arms}

\item{R}{If the number of levels per arm is constant, an integer of the
number of levels. If not, a vector giving the levels for each arm}
}
\value{
A sigma matrix with 1's in all valid cut locations and NAs if
we cannot cut at that location (e.g. the level for that arm does not exist).
}
\description{
Helper function used to initialize sigma
given M arms and a vector R (or a single number)
of levels for each arm
}
