#' Makes output from aggregate_rashomon_partitions cleaner
#' @export
make_rashomon_objects <- function(rset) {
  rash_models <- list()
  length(rash_models) <- length(rset[[1]])
  count <- 0
  for (x in rset[[1]]) {
    count <- count + 1
    rashomon_i <- new_RashomonSet(
      models = list(),
      losses = numeric(),
      num_pools = list(),
      profiles = list(),
      pool_dictionaries = list()
    )
    for (i in 1:length(rset[[2]])) {
      model_i <- rset[[2]][[i]]$models[[x[[i]]]]
      loss_i <- rset[[2]][[i]]$losses[[x[[i]]]]
      pools_i <- rset[[2]][[i]]$num_pools[[x[[i]]]]
      profiles_i <- rset[[2]][[i]]$profiles[[x[[i]]]]
      pool_dictionary_i <- rset[[2]][[i]]$pool_dictionaries[[x[[i]]]]

      rashomon_i = insert_model(rashomon_i, list(model_i), loss_i, pools_i, list(profiles_i), pool_dictionary_i)
    }

    rash_models[[count]] = rashomon_i
  }
  rash_models
}
