% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count.R
\name{num_pools}
\alias{num_pools}
\title{Helper function used to find number of pools given arbitrary sigma and
number of levels in each arm.}
\usage{
num_pools(sigma, R = NA)
}
\arguments{
\item{sigma}{Partition matrix for a given pooling structure}

\item{R}{A list (or integer) of the number of levels in each arm.}
}
\value{
The number of pools given by the pooling structure sigma.
}
\description{
Helper function used to find number of pools given arbitrary sigma and
number of levels in each arm.
}
