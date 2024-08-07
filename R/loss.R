
#' Finds policy labels given data with each observation labeled with its policy
#' TODO: Change all function calls from other packages to be explicit, and
#' change if necessary for speedups.
#' @param data User data appended with policy labels
#' @param value Name of value column observed for each observation, supplied by the user.
#' @returns A dataframe of policy means, where each row is the mean for a given policy.
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @export
policy_means <- function(data, value){

  #convert to data.table
  data %>%
    group_by(policy_label) %>%
    summarize(sum = sum({{value}}, na.rm = TRUE),
              n = n(),
              mean = mean({{value}}),
              .groups = "drop",
              policy_label = policy_label[1])
}

#' Assigns pools to observations in data given a dictionary of pools.
#'
#' @param arr A list
#' @returns A list of the sums of the products of size k, where the (k+1)th element
#' denotes the kth sum. The 1st element is 1, for use in other functions for easy looping.
#' @importFrom collections dict
#' @export
pools_to_data <- function(data, pools_dict){

  policy_label <- as.integer(dplyr::pull(data,policy_label))

  len_data = length(policy_label)

  pool_label = numeric(len_data)

  for(i in 1:len_data){
    pool_label[i] = pools_dict$get(policy_label[i])
  }

  data$pool = pool_label
  data

}


#' Finds pool means given data, value, and pool label
#' TODO: Change all function calls from other packages to be explicit, and
#' change if necessary for speedups.
#' @param data User data appended with policy labels
#' @param value Name of value column observed for each observation, supplied by the user.
#' @param pool Name of pool column
#' @returns A collections::dict() object, where the key is an integer i corresponding
#' to the pool id and its value is the mean of that pool.
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @export
pool_means <- function(data, value, pool){

  means <- data %>%
    group_by({{pool}}) %>%
    summarize(mean_pool = mean({{value}}))

  pool_means_dict = collections::dict()
  for(i in 1:nrow(means)){
    pool_means_dict$set(as.integer(pull(means,{{pool}})[i]), pull(means,mean_pool)[i])
  }

  pool_means_dict
}


#' Maximally cuts in row i for use in the lookahead loss, up to column j.
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param i Row of the sigma matrix to cut at
#' @param j Column to cut up to
#' @returns A new partition matrix cut from row i to column j
partition_sigma <- function(i, j, sigma) {

  new_sigma = sigma
  new_sigma[i, j:ncol(new_sigma)] = 0
  new_sigma[is.na(sigma)] = NA

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
extract_pools <- function(policy_list, sigma, lattice_edges = NA){

  if(any(is.na(lattice_edges))){
    lattice_relations = lattice_edges(sigma, policy_list)
  }
  else{
    lattice_relations = prune_edges(sigma, lattice_edges, policy_list)
  }


  pools = connected_components(length(policy_list), lattice_relations)

  pools


}


#' MSE loss for data (y_i - pool_mean_i)^2 for i in pool i.
#'
#' @param data User supplied dataframe that contains values
#' @param value Label of column containing value observed for each individual
#' @param M Dataframe of policy means
#' @param sigma Partition matrix that gives pooling structure
#' @param policy_list List of policies, which if a list of vectors of all policies implied
#' by the data.
#' @param reg Regularization parameter that penalizes partitions with more pools
#' @param normalize Whether or not to normalize loss (unsure?)
#' @param lattice_edges Edges of pooling structure
#' @importFrom collections dict
#' @import magrittr
#' @import dplyr
#' @returns MSE given sigma pooling structure and data.
#' @export
compute_mse_loss <- function(data, value, M, sigma, policy_list, reg = 1, normalize = 0, lattice_edges = NA){

  #Compute pools for new maximal split
  pool_dict = extract_pools(policy_list, sigma, lattice_edges)

  # assigning pools to policy means so that we can compute pool means
  M_pool <- pools_to_data(M, pool_dict)

  #dictionary of pool means of type cc:dictionary()
  fixed_pool_means_dict = pool_means(M_pool, mean, pool)

  #assign pool labels to data and extract pool labels
  pools_data = pools_to_data(data,pool_dict)$pool

  #vector for storing pool mean for each observation
  pool_mean_data <- numeric(length(pools_data))

  #assigning pool mean to each observation
  for(k in 1:length(pools_data)){
    pool_mean_data[k] = fixed_pool_means_dict$get(as.integer(pools_data[k]))
  }

  y = dplyr::pull(data,{{value}})

  mse = (yardstick::rmse_vec(pool_mean_data, y))^2

  if(normalize > 0){
    mse = (mse * nrow(data) / normalize)
  }

  mse
}

#' Computes penalization loss for data (regularization param * number of pools)
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param R A list (or integer) of the number of levels in each arm.
#' @param reg Regularization parameter that penalizes partitions with more pools
#' @returns Penalization loss
#' @export
compute_penalization_loss <- function(sigma, R, reg){
  if(all(is.na(sigma))){
    return(1)
  }
  else{
    return(num_pools(sigma, R) * reg)
  }

}

#' Lookahead loss (B) to know when we need not continue with a given sigma.
#'
#' @param data User supplied dataframe that contains values
#' @param value Label of column containing value observed for each individual
#' @param i Row to split at in sigma
#' @param j Column to split from in sigma
#' @param M Dataframe of policy means
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
compute_B <- function(data, value, i,j, M, sigma, policy_list, reg = 1, normalize = 0, lattice_edges = NA, R){

  #Split maximally across row, starting at point i, j:
  sigma_max_split = partition_sigma(i,j,sigma)

  mse = compute_mse_loss(data, {{value}}, M, sigma_max_split, policy_list, reg = reg, normalize = normalize, lattice_edges)

  #least number of pools
  #least bad penalty for complexity
  sigma_max_split[i, (j+1):ncol(sigma_max_split)] = 1
  sigma_max_split[is.na(sigma)] = NA

  reg_loss = compute_penalization_loss(sigma_max_split, R, reg)

  B = mse + reg_loss

  B

}
#' Loss given a specific pooling stucture (sigma)
#'
#' @param data User supplied dataframe that contains values
#' @param value Label of column containing value observed for each observation
#' @param i Row to split at in sigma
#' @param j Column to split from in sigma
#' @param M Dataframe of policy means
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
#' @returns Loss given pool for the data
compute_loss <- function(data, value, M, sigma, policy_list, reg = 1, normalize = 0, lattice_edges = NA, R){

  mse = compute_mse_loss(data, {{value}}, M, sigma, policy_list, reg = 1, normalize = normalize, lattice_edges)
  reg_loss = compute_penalization_loss(sigma, R, reg)

  mse + reg_loss
}
