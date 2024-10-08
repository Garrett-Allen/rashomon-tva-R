% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RashomonSet.R
\name{predict.RashomonSet}
\alias{predict.RashomonSet}
\title{Make predictions from a model in the RashomonSet}
\usage{
\method{predict}{RashomonSet}(rashomon_set, universal_labels, model_id = 1)
}
\arguments{
\item{rashomon_set}{A RashomonSet object.}

\item{universal_labels}{A vector of universal_labels that you want to
extract predictions for. This vector gives the unique policy id for each
new observation you want to predict for, as assignd by assign_universal_labels().
Make sure that the assigned labels are the same as those present in the data
when aggregate_rashomon_profiles() is called.}

\item{model_id}{Which model in the RashomonSet you want a prediction from. Defaults to 1. If you're calling
this on a RashomonSet from the output of make_rashomon_objs, leave this at the default, as there is one
pool dictionary per object that holds all of the mappings from unique policy ids to pool means.}
}
\description{
Allows you to get predictions from a model in the RashomonSet.
}
