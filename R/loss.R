#' Finds policy labels given data with each observation labeled with its policy
#' TODO: Change all function calls from other packages to be explicit, and
#' change if necessary for speedups.
#' @param data Data.table containing the treatment assignments
#' @param value Name of value column observed for each observation, supplied by the user. Must be character.
#' @returns A Data.table of policy means, where each row is the mean for a given policy.
#' @import magrittr
#' @import dplyr
#' @export
policy_means <- function(data, value) {

  result <- data[, .(
    sum = sum(.SD[[1]], na.rm = TRUE),
    n = .N,
    mean = mean(.SD[[1]], na.rm = TRUE),
    policy_label = policy_label[1],
    universal_label = universal_label[1]
  ), by = policy_label, .SDcols = value]

  return(result)
}

#' Assigns pools to observations in data given a dictionary of pools.
#'
#' @param data Data.table containing the treatment assignments
#' @param pools_dict A collections::dict object from extract_pools() that
#' takes in a policy id and outputs the pool it is in.
#' @returns A list of the sums of the products of size k, where the (k+1)th element
#' denotes the kth sum. The 1st element is 1, for use in other functions for easy looping.
#' @importFrom collections dict
#' @export
pools_to_data <- function(data, pools_dict) {
  policy_label <- as.integer(data[,policy_label])

  len_data <- length(policy_label)

  pool_label <- numeric(len_data)

  for (i in 1:len_data) {
    pool_label[i] <- pools_dict$get(policy_label[i])
  }

  data$pool <- pool_label
  data
}


#' Finds pool means given data, value, and pool label
#' @param data Data.table containing the treatment assignments
#' @param value Name of value column observed for each observation, supplied by the user. Must be character.
#' @returns A collections::dict() object, where the key is an integer i corresponding
#' to the pool id and its value is the mean of that pool.
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @export
pool_means <- function(data, value){

  data[, mean_pool := mean(.SD[[1]], na.rm = TRUE), by = pool, .SDcols = value]

  univ_labels <- as.integer(data[,universal_label])
  mean_pools <- data[,mean_pool]

  pool_means_dict <- collections::dict(items = mean_pools, keys = univ_labels)

  pool_means_dict
}
#TODO I CHANGED HERE, universal_label


#' Maximally cuts in row i for use in the lookahead loss, up to column j.
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param i Row of the sigma matrix to cut at
#' @param j Column to cut up to
#' @returns A new partition matrix cut from row i to column j
#' @noRd
partition_sigma <- function(i, j, sigma) {
  new_sigma <- sigma
  new_sigma[i, j:ncol(new_sigma)] <- 0
  new_sigma[is.na(sigma)] <- NA

  new_sigma
}

#' Helper function to extract pools from a given sigma matrix
#'
#' @param policy_list A list of policies
#' @param sigma Row of the sigma matrix to cut at
#' @param lattice_edges List of edges between policies given a pooling structure.
#' @importFrom collections dict
#' @returns A collections:dict object of all the pools, where the key is integer
#' id of a policy and the value is the pool id.
#' @export
extract_pools <- function(policy_list, sigma, lattice_edges = NA) {
  if (!is.list(lattice_edges)) {
    lattice_relations <- lattice_edges(sigma, policy_list)
  } else {
    lattice_relations <- prune_edges(sigma, lattice_edges, policy_list)
  }


  pools <- connected_components(length(policy_list), lattice_relations)

  pools
}


#' MSE loss for data (y_i - pool_mean_i)^2 for i in pool i.
#'
#' @param data Data.table containing the treatment assignments
#' @param value Label of column containing value observed for each individual
#' @param M Data.table of policy means
#' @param sigma Partition matrix that gives pooling structure
#' @param policy_list List of policies, which if a list of vectors of all policies implied
#' by the data.
#' @param reg Regularization parameter that penalizes partitions with more pools
#' @param normalize Whether or not to normalize loss
#' @param lattice_edges Edges of pooling structure
#' @param return_dict Whether or not to return the dictionary object that gives the pool
#' means for each model in the RashomonSet. Defaults to True.
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @returns MSE given sigma pooling structure and data.
#' @export
compute_mse_loss <- function(data, value, M, sigma, policy_list, reg = 1, normalize = 0, lattice_edges = NA, return_dict = TRUE) {
  # Compute pools for new maximal split
  pool_dict <- extract_pools(policy_list, sigma, lattice_edges)

  # assigning pools to policy means so that we can compute pool means
  M_pool <- pools_to_data(M, pool_dict)


  # dictionary of pool means of type cc:dictionary()
  fixed_pool_means_dict <- pool_means(M_pool, "mean")

  # assign pool labels to data and extract pool labels
  universal_policy_labels <- data$universal_label

  n <- length(universal_policy_labels)
  # vector for storing pool mean for each observation
  pool_mean_data <- numeric(n)

  # assigning pool mean to each observation
  for (k in 1:n) {
    pool_mean_data[k] <- fixed_pool_means_dict$get(as.integer(universal_policy_labels[k]))
  }

  y <- dplyr::pull(data, value)

  mse <- (yardstick::rmse_vec(pool_mean_data, y, na_rm = TRUE))^2

  if (normalize > 0) {
    n_data_i <- sum(!is.na(y))
    mse <- (mse * n_data_i / normalize)
  }

  if(!return_dict){
    return(mse)
  }

  return(list(mse = mse,
              dictionary = fixed_pool_means_dict))
}

#' Computes penalization loss for data (regularization parameter * number of pools)
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param R A list (or integer) of the number of levels in each arm.
#' @param reg Regularization parameter that penalizes partitions with more pools
#' @returns Penalization loss
#' @export
compute_penalization_loss <- function(sigma, R, reg) {
  if (all(is.na(sigma))) {
    return(1)
  } else {
    return(num_pools(sigma, R) * reg)
  }
}

#' Look-ahead loss (B) to know when we need not continue with a given sigma.
#'
#' @param data User supplied Data.table that contains values
#' @param value Label of column containing value observed for each individual
#' @param i Row to split at in sigma
#' @param j Column to split from in sigma
#' @param M Data.table of policy means
#' @param sigma Partition matrix that gives pooling structure
#' @param policy_list List of policies, which if a list of vectors of all policies implied
#' by the data.
#' @param reg Regularization parameter that penalizes partitions with more pools
#' @param normalize Whether or not to normalize loss (unsure?)
#' @param lattice_edges Edges of pooling structure
#' @param R A list (or integer) of the number of levels in each arm.
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @export
#' @returns MSE given sigma pooling structure and data
compute_B <- function(data, value, i, j, M, sigma, policy_list, reg = 1, normalize = 0, lattice_edges = NA, R) {
  # Split maximally across row, starting at point i, j:
  sigma_max_split <- partition_sigma(i, j, sigma)

  mse <- compute_mse_loss(data, value, M, sigma_max_split, policy_list, reg = reg, normalize = normalize, lattice_edges, return_dict = FALSE)

  # least number of pools
  # least bad penalty for complexity
  sigma_max_split[i, (j + 1):ncol(sigma_max_split)] <- 1
  sigma_max_split[is.na(sigma)] <- NA

  reg_loss <- compute_penalization_loss(sigma_max_split, R, reg)

  B <- mse + reg_loss

  B
}
#' Loss given a specific pooling stucture (sigma)
#'
#' @param data User supplied Data.table that contains values
#' @param value Label of column containing value observed for each observation
#' @param M Data.table of policy means
#' @param sigma Partition matrix that gives pooling structure
#' @param policy_list List of policies, which if a list of vectors of all policies implied
#' by the data.
#' @param reg Regularization parameter that penalizes partitions with more pools
#' @param normalize Whether or not to normalize loss
#' @param lattice_edges Edges of pooling structure
#' @param R A list (or integer) of the number of levels in each arm.
#' @param return_dict Whether or not to return the dictionary object that gives the pool
#' means for each model in the RashomonSet. Defaults to True.
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @export
#' @returns Loss given pool for the data
compute_loss <- function(data, value, M, sigma, policy_list, reg = 1, normalize = 0, lattice_edges = NA, R, return_dict = TRUE) {
  mse <- compute_mse_loss(data, value, M, sigma, policy_list, reg = 1, normalize = normalize, lattice_edges, return_dict = TRUE)
  reg_loss <- compute_penalization_loss(sigma, R, reg)


  if(!return_dict){
    loss = mse + reg_loss
    return(loss)
  }

  loss = mse[[1]] + reg_loss

  return(list(loss = loss, dict = mse[[2]], num_pools = reg_loss / reg))
}
