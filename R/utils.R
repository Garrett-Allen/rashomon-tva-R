#' @description Finds the policies present in a given profile.
#' @param policies A list of policies, typically given by generate_policies_from_data
#' @param profile The profile we wish to see which policies are active form
#' @param inactive The level in policies that denotes an arm is inactive.
#'
#' @returns A list of policies present in a given profile
#' @noRd
pol_in_prof <- function(policies, profile, inactive = 0){

  policy_profiles = lapply(policies, function(x) x != inactive)
  sapply(policy_profiles, function(x) all(x == profile))

}

#' @description Finds subset of a dataset corresponding to a given profile
#' @param data A dataframe output from assign_policy_label
#' @param policy_list A list of policies given by generate_policies_from_data
#' @param profile The profile we wish to see which policies are active form
#' @param inactive The level in policies that denotes an arm is inactive.
#'
#' @returns The subset of a data that corresponds to the given profile
#' @noRd
subset_prof <- function(data, policy_list, profile,inactive = 0){

  policies_in_profile = pol_in_prof(policy_list, profile, inactive)
  subset(data, policies_in_profile[data$policy_label])


}

#' @description Finds lower bound for a given profile
#' @param data Data prefiltered to be from a given profile
#' @param value The column name of the y values in data
#'
#'
#' @returns Bound on how low a model's loss can be in a given profile.
#' @noRd
find_profile_lower_bound <- function(data, value){

  n_k = nrow(data)
  data_mean <- data %>%
    group_by(policy_label) %>%
    mutate(mean = mean({{value}}))

  (yardstick::rmse_vec(pull(data_mean,{{value}}), data_mean$mean))^2 * n_k

}

#' @description Finds combinations of models such that the sum of their losses is less than
#' theta, the maximum number of pools is less than H, and each profile has a model
#' that gives the pooling structure for that profile.
#' @param rashomon_profiles A list of RashomonSet objects, where each entry of the list
#' contains a RashomonSet corresponding to a different profile
#' @param theta Threshold value to be present in the RashomonSet.
#' @param H The maximum number of pools that can be present across all profiles.
#' @param sorted Whether or not the RashomonSet objects have been sorted by the internal $sort method.
#' Defaults to FALSE
#'
#' @returns A list of feasible combinations of poolings such that the sum of the losses across
#' profiles is less than theta and has less than H total pools. The ith entry in
#' each list corresponds to the ith RashomonSet profile, and the value of
#' the ith entry represents the model in the ith RashomonSet profile.
#' @noRd
find_feasible_combinations <- function(rashomon_profiles, theta, H, sorted = FALSE){

  #sorting so we can pass into the feasible_sums function
  if(!sorted){
    for(rset in rashomon_profiles){
      rset$sort()
    }
  }

  #coalescing all losses into one list of lists.
  all_losses <- lapply(rashomon_profiles, function(x) x$losses)

  loss_combinations = find_feasible_sum_subsets(all_losses, theta)

  feasible_combinations = list()

  #filtering so that we only have combinations with number of pools smaller
  #than H
  for(i in 1:length(loss_combinations)){
    pools = 0
    comb = loss_combinations[[i]]
    for(j in 1:length(comb)){

      r_prof = rashomon_profiles[[j]]
      model_id = comb[[j]]

      if(all(is.na(r_prof$models[[model_id]]))){
        if(r_prof$losses[[model_id]] > 0){
          pools = pools + 1
        }
      }

      else{
        pools = pools + r_prof$pools[[model_id]]
      }


    }

    if(pools <= H){
      feasible_combinations = append(feasible_combinations, list(comb))
    }

  }

  feasible_combinations
}

#' @description Find all combinations of indices from the list of arrays S such the sum of those values
#' is less than theta.
#' @param S List of lists of numbers
#' @param theta Threshold
#' @returns List of combinations of indices, one from each array in S
#' such that their sum is less than or equal to theta
#' @noRd
find_feasible_sum_subsets <- function(S,theta){

  numsets = length(S)

  if(numsets == 0){
    return(list())
  }

  S1 = S[[1]]
  S1_feasible_ids = which(S1 <= theta)

  if(numsets == 1){
    feasible_combs = lapply(S1_feasible_ids, list)
    return(feasible_combs)
  }

  # Since the list is sorted, add the first elements of the remaining sets
  # Then use this to check feasibility of indices in S1_feasible_idx
  first_element_sum = 0

  for(x in S[2:numsets]){
    first_element_sum = first_element_sum + x[1]
  }

  feasible_combs = list()
  for(i in S1_feasible_ids){
    theta_i = theta - S1[i]

    if(first_element_sum > theta_i){
      next
    }

    subproblem_res_i = find_feasible_sum_subsets(S[2:numsets], theta_i)
    subproblem_res_i = lapply(subproblem_res_i, function(x) c(i,x))

    feasible_combs = append(feasible_combs, subproblem_res_i)
  }

  feasible_combs

}
