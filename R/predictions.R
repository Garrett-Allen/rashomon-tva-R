
#' @title Make output from aggregate_rashomon_profiles into a list of models in RashomonSet
#' @description Makes output from aggregate_rashomon_profiles usable in predict function,
#' by making each object in the list a model in the RashomonSet, with each possessing dictionary
#' that maps universal_labels from assign_universal_labels() to pool means.
#' @param rset A list object from the output of aggregate_rashomon_profiles()
#' @importFrom purrr list_flatten
#' @importFrom collections dict
#' @export
make_rashomon_objects <- function(rset) {

  rash_models <- list()
  length(rash_models) <- length(rset[[1]])
  count <- 0
  num_profiles <- length(rset[[2]])

  for (x in rset[[1]]) {
    count <- count + 1
    rashomon_i <- new_RashomonSet(
      models = list(),
      losses = numeric(),
      num_pools = list(),
      profiles = list(),
      pool_dictionaries = list()
    )
    models <- list()
    losses <- list()
    pools <- list()
    profiles <- list()
    pool_dicts <- list()

    length(models) <- num_profiles
    length(losses) <- num_profiles
    length(pools) <- num_profiles
    length(profiles) <- num_profiles
    length(pool_dicts) <- num_profiles
    for (i in 1:length(rset[[2]])) {
      model_i <- rset[[2]][[i]]$models[[x[[i]]]]
      loss_i <- rset[[2]][[i]]$losses[[x[[i]]]]
      pools_i <- rset[[2]][[i]]$num_pools[[x[[i]]]]
      profiles_i <- rset[[2]][[i]]$profiles[[x[[i]]]]
      pool_dictionary_i <- rset[[2]][[i]]$pool_dictionaries[[x[[i]]]]

      models[[i]] <- model_i
      losses[[i]] <- loss_i
      pools[[i]] <- pools_i
      profiles[[i]] <- profiles_i
      pool_dicts[[i]] <- pool_dictionary_i


    }

    flattened_dict <- purrr::list_flatten(lapply(pool_dicts, function(x) x$as_list()))

    pool_dicts <- list(collections::dict(flattened_dict))

    rashomon_i <- new_RashomonSet(models, losses, pools, profiles, pool_dicts)

    rash_models[[count]] = rashomon_i
  }
  rash_models
}



