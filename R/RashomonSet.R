#' Inserts new model and it's corresponding loss into a RashomonSet Object
#'
#' @param new_sigma A M X R-2 partition matrix
#' @param loss The loss of the model when using this pooling structure
#' @param num_pools Number of pools present in new_sigma
#' @param insert_profile Profile for the given sigma matrix.
insert_model <- function(new_sigma, loss, num_pools, profile){

  models <<- append(models, list(new_sigma))
  losses <<- c(losses,loss)
  pools <<- append(pools, num_pools)
  profiles <<- append(profiles, profile)
}

#' Sorts this rashomon set according to the model with the smallest loss.

sort <- function(){
  order <- order(losses)
  losses <<- losses[order]
  models <<- models[order]
  pools  <<- pools[order]
  profiles <<- profiles[order]
}


#' RashomonSet class. This can be used to both form the RashomonSet for a given
#' profile and a single member of the RashomonSet across all profiles.
#' @field models The M X R-2 partition matrices that give the pooling structure for
#' their profile
#' @field losses The losses for each of the models when evaluated on their profile
#' @field pools The number of pools in the each of the models when evaluated on their profile
#' @field profiles The profiles for each of the models.
#' @export RashomonSet
#' @exportClass RashomonSet
RashomonSet <- setRefClass("RashomonSet", fields = list(models = "list",
                                                        losses = "numeric",
                                                        pools = "list",
                                                        profiles = "list"),
                           methods = list(insert_model = insert_model,
                                          sort = sort))


