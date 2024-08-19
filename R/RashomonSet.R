#' @title Rashomon Set
#' @description This class stores members of the RashomonSet. This can be used to both form the RashomonSet for a given
#' profile and a single member of the RashomonSet across all profiles.
#' @param models The M X R-2 partition matrices that give the pooling structure for
#' their profile
#' @param losses The losses for each of the models when evaluated on their profile
#' @param num_pools The number of pools in the each of the models when evaluated on their profile
#' @param profiles The profiles for each of the models.
#' @param pool_dictionaries A list of collections::dict() objects that maps policy_ids
#' to pool_means for the given pooling structure in model.
#' @export
new_RashomonSet <- function(models = list(),
                            losses = numeric(),
                            num_pools = list(),
                            profiles = list(),
                            pool_dictionaries = list()){

  RashomonSet <- list(models = models,
                      losses = losses,
                      num_pools = num_pools,
                      profiles = profiles,
                      pool_dictionaries = pool_dictionaries)

  class(RashomonSet) <- "RashomonSet"

  RashomonSet

}

#' Inserts new model and it's corresponding loss into a RashomonSet Object
#' #TODO: Consider changing this to S3
#' @param obj A RashomonSet object
#' @param new_sigma A M X R-2 partition matrix
#' @param new_loss The loss of the model when using this pooling structure
#' @param new_num_pools Number of pools present in new_sigma
#' @param new_profile Profile for the given sigma matrix.
#' @param new_dict A dictionary that maps policy_ids to pool_means for the giving pooling structure in model.
#' @export
insert_model <- function(obj, new_sigma, new_loss, new_num_pools, new_profile, new_dict) {
  obj$models <- append(obj$models, new_sigma)
  obj$losses <- append(obj$losses, new_loss)
  obj$num_pools <- append(obj$num_pools, new_num_pools)
  obj$profiles <- append(obj$profiles, new_profile)
  obj$pool_dictionaries <- append(obj$pool_dictionaries, new_dict)

  obj
}

#' Sorts this rashomon set according to the model with the smallest loss.
#' @param obj A RashomonSet object
#' @export
sort_rashomon <- function(obj){
  order <- order(obj$losses)

  obj$models <- obj$models[order]
  obj$losses <- obj$losses[order]
  obj$num_pools <- obj$num_pools[order]
  obj$profiles <- obj$profiles[order]
  obj$pool_dictionaries <- obj$pool_dictionaries[order]

  obj
}

#' Combines all of the dictionaries of a RashomonSet object to form one dictionary.
#' @param obj A RashomonSet object
#' @param newdata Data that we want to predict given our model object. Must contain
#' column universal_labels from output of assign_universal_labels(), and these labels
#' must correspond to the same policies as in the data passed to aggregate_rashomon_profiles()
#' or find_rashomon_profile()
#' @export
combine_pool_dictionaries <- function(obj){

  total_dict <- list()

  for(dict in obj$pool_dictionaries){
    total_dict <- append(total_dict, dict$as_list())
  }

  collections::dict(items = unname(total_dict), keys = names(total_dict))
}

#' Gernerates predictions from a model in the RashomonSet
#' @param obj A RashomonSet object
#' @param newdata A dataframe. Must contain column universal_labels from output of assign_universal_labels(), and these labels
#' must correspond to the same policies as in the data passed to aggregate_rashomon_profiles()
#' or find_rashomon_profile()
#' @export
predict.RashomonSet <- function(obj, data){

  total_dict <- combine_pool_dictionaries(obj)

  pred <- c()
  for(i in 1:nrow(data)){
    pred[i] <- total_dict$get(as.character(data$universal_label[i]))
  }
  data$predictions <- pred
  data
}

