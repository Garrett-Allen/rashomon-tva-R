#' @title Finds the RashomonSet for a given profile.
#'
#' @param data A dataframe containing the column whose name you supply in value.
#' @param value The column name of the y values in data
#' @param arm_cols A character vector containing the names of the arm columns
#' @param M The number of arms
#' @param R An integer (or vector) denoting the number of levels in each arm. If
#' this value is an integer, it assumes each arm has the same number of levels.
#' @param H The maximum number of pools allowed to be specified by a partition matrix
#' @param reg The value of the regularization penalty when computing the loss for
#' a pooling structure. Defaults to 1, so that each additional pool in a model adds 1 to the loss.
#' @param profile Binary vector that denotes which arms are active.
#' @param policies A list of policies present in the data. This should be the output
#' of create_policies_from_data.
#' @param policy_means A dataframe that gives the mean value of value for each policy.
#' @param normalize Normalization factor. When this is 0, no normalization is applied. Should be set
#' to the size of the total dataset, so that each loss represents contribution to the overall
#' loss of the model across all profiles. Defaults to zero.
#' @param theta Threshold value to be present in the RashomonSet. Set to Inf if you want all poolings.
#' @param filtered Whether or not data is already filtered to be only the observations in the given
#' profile. Defaults to True.
#' @param inactive The level that denotes an inactive arm in data. Defaults to zero.
#' @returns A RashomonSet Object that gives all of the partition matrices in the
#' Rashomon Set for a given theta.
#' @import utils
#'
#' @export

find_rashomon_profile <- function(data,
                                  value,
                                  arm_cols = c(),
                                  M,
                                  R,
                                  H,
                                  reg = 1,
                                  profile,
                                  policies = c(),
                                  policy_means = c(),
                                  normalize = 0,
                                  theta,
                                  filtered = FALSE,
                                  inactive = 0
                                  ){

  if(is.null(data$universal_label)){
    warning("No universal label assigned; please run your data through assign_universal_label() to use pool_dictionaries in output")
    data <- assign_universal_label(data, arm_cols)
  }

  if (!filtered) {
    data <- subset_prof(data, policy_list = policies, profile = profile, inactive = inactive)
  }

  if (max(R) == 2) {
    # TODO: Verify this is correct thing to do
    sigma <- matrix(nrow = M, ncol = 1)
    y <- pull(data, value)
    mean <- mean(y)
    mse <- mean((y - mean)^2)

    if (normalize > 0) {
      mse <- mse * nrow(data) / normalize
    }

    num_pools <- 1
    loss <- mse + reg * num_pools

    rashomon_set <- new_RashomonSet(
      models = list(sigma),
      losses = c(loss),
      num_pools = list(num_pools),
      profiles = list(list(profile),
      pools_dictionaries = list(mean))
    )
    return(rashomon_set)
  }

  sigma <- initialize_sigma(M, R)
  hasse_edges <- lattice_edges(sigma, policies)

  # creating list for compatibility
  if (length(R) == 1) {
    R <- rep(R, M)
  }

  # defining list, hash tables, and queue
  rashomon_set <- new_RashomonSet(
    models = list(),
    losses = numeric(),
    num_pools = list(),
    profiles = list(),
    pool_dictionaries = list()
  )

  seen_sigmas <- hashtab()
  seen_sigmas_sub <- hashtab()
  queue_sigma <- collections::queue(items = NULL)
  # pushing initial problems onto queue
  for (i in 1:M) {
    if (!is.na(sigma[i, 1])) {
      queue_sigma$push(list(sigma, i, 0))
    }
  }
  count <- 0
  # main queue loop

  while (queue_sigma$size() > 0) {
    count <- count + 1
    sigma_list <- queue_sigma$pop()

    sigma <- sigma_list[[1]]
    i <- sigma_list[[2]]
    j <- sigma_list[[3]]
    sethash(seen_sigmas_sub, sigma_list, 1)

    if (num_pools(sigma, R) > H) {
      next
    }

    sigma_1 <- sigma
    sigma_0 <- sigma
    sigma_1[i, j] <- 1
    sigma_0[i, j] <- 0

    # adding all subproblem variants
    for (m in 1:M) {
      R_m <- R[m]
      j1 <- 1
      sigma_1_list <- list(sigma_1, m, j1)

      # have we seen this subproblem?
      seen <- gethash(seen_sigmas_sub, sigma_1_list, nomatch = 0)


      # if we have seen it, and our cutting index is less than R_m - 2, consider further cuts
      while (seen != 0 & j1 <= R_m - 2) {
        j1 <- j1 + 1
        sigma_1_list <- list(sigma_1, m, j1)
        seen <- gethash(seen_sigmas_sub, sigma_1_list, nomatch = 0)
      }

      # if we haven't seen it and cutting is at a vlaid range, add it as a subproblem
      if (j1 <= R_m - 2 & seen == 0) {
        queue_sigma$push(sigma_1_list)
        sethash(seen_sigmas_sub, sigma_1_list, 1)
      }

      j0 <- 1
      sigma_0_list <- list(sigma_0, m, j0)
      seen <- gethash(seen_sigmas_sub, sigma_0_list, nomatch = 0)

      while (seen != 0 & j0 <= R_m - 2) {
        j0 <- j0 + 1
        sigma_0_list <- list(sigma_0, m, j0)
        seen <- gethash(seen_sigmas_sub, sigma_0_list, nomatch = 0)
      }

      if (j0 <= R_m - 2 & seen == 0) {
        queue_sigma$push(sigma_0_list)
        sethash(seen_sigmas_sub, sigma_0_list, 1)
      }
    }

    # computing B to check whether more splits is possible

    # i had to add this condition... why does apara's code not need it?
    if (j != ncol(sigma)) {
      B <- compute_B(data,
                     value,
                     i,
                     j,
                     policy_means,
                     sigma,
                     policies,
                     reg = reg,
                     normalize = normalize,
                     lattice_edges = hasse_edges,
                     R
      )


      if (B > theta) {
        next
      }
    }

    # Check if unsplit pool satisfies Rashomon threshold
    if (gethash(seen_sigmas, sigma_1, nomatch = 0) == 0) {
      sethash(seen_sigmas, sigma_1, 1)

      Q <- compute_loss(data,
                        value,
                        policy_means,
                        sigma_1,
                        policies,
                        reg = reg,
                        normalize = normalize,
                        lattice_edges = hasse_edges,
                        R,
                        return_dict = TRUE
      )

      if (Q[[1]] <= theta) {
        n_pools <- Q[[3]]
        rashomon_set = insert_model(rashomon_set, list(sigma_1), Q[[1]], n_pools, list(profile), Q[[2]])
      }
    }
    # Check if split pool satisfies Rashomon threshold
    if (gethash(seen_sigmas, sigma_0, nomatch = 0) == 0 & num_pools(sigma_0, R) <= H) {
      sethash(seen_sigmas, sigma_0, 1)

      Q <- compute_loss(data,
                        value,
                        policy_means,
                        sigma_0,
                        policies,
                        reg = reg,
                        normalize = normalize,
                        lattice_edges = hasse_edges,
                        R,
                        return_dict = TRUE
      )

      if (Q[[1]] <= theta) {
        n_pools <- Q[[3]]
        rashomon_set = insert_model(rashomon_set, list(sigma_0), Q[[1]], n_pools, list(profile), Q[[2]])
      }
    }

    sigma_1_list <- list(sigma_1, i, j + 1)
    sigma_0_list <- list(sigma_0, i, j + 1)

    # add child problems
    if (j + 1 < R[i] - 2) {
      if (gethash(seen_sigmas_sub, sigma_1_list, nomatch = 0) == 0) {
        queue_sigma$push(list(sigma_1, i, j + 1))
      }
      if (gethash(seen_sigmas_sub, sigma_0_list, nomatch = 0) == 0) {
        queue_sigma$push(list(sigma_0, i, j + 1))
      }
    }
  }
  if(length(rashomon_set$models) == 0){
    return(NULL)
  }
  rashomon_set
}
#' @title Finds the RashomonSet of models for a given dataset.
#'
#' @param data A dataframe containing the column whose name you supply in value
#' @param value The column name of the y values in data
#' @param arm_cols A character vector containing the names of the arm columns
#' @param M The number of arms
#' @param R An integer (or vector) denoting the number of levels in each arm. If
#' this value is an integer, it assumes each arm has the same number of levels.
#' @param H The maximum number of pools allowed to be specified by a partition matrix
#' @param reg The value of the regularization penalty when computing the loss for
#' a pooling structure. Defaults to 1, so that each additional pool in a model adds 1 to the loss.
#' @param theta Threshold value to be present in the RashomonSet. Set to Inf if you want all poolings.
#' @param bruteforce Whether the RashomonSet should be found via bruteforce. Currently unimplemented,
#' defaults to FALSE.
#' @param inactive The level that denotes an inactive arm in data. Defaults to zero.
#' @returns A list of RashomonSet Objects that gives all of the partition matrices in the
#' Rashomon Set for a given theta.
#'
#' @export
aggregate_rashomon_profiles <- function(data,
                                        value,
                                        arm_cols,
                                        M,
                                        R,
                                        H,
                                        reg = 1,
                                        theta,
                                        bruteforce = FALSE,
                                        inactive = 0) {


  if(is.null(data$universal_label)){
    warning("No universal label assigned; please run your data through prep_data() to use pool_dictionaries in output")
  }

  num_profiles <- 2^M
  profiles <- expand.grid(replicate(M, 0:1, simplify = FALSE))
  # creating policy labels on data as well as list of all policies

  if (length(R) == 1) {
    R <- rep(R, M)
  }

  data_labeled <- prep_data(data, arm_cols = arm_cols, value = value, R = R, drop_unobserved_combinations = FALSE)
  policy_list <- create_policies_from_data(data_labeled, arm_cols)
  num_data <- nrow(data)
  num_ids <- nrow(data_labeled)
  data_labeled$id <- 1:num_ids

  # Maximum number of pools for a profile derived from maximum number of pools
  H_profile <- H - num_profiles + 1
  data_profile_ids <- rep(0, num_ids)

  # initialize storage objects
  eq_lb_profiles <- rep(0, num_profiles)
  rashomon_profiles <- rep(0, num_profiles)
  loss_object <- rep(0, num_profiles)
  D_profile <- list()

  # Assign profile ids to each data point
  control_univ_id = -1
  control_no_data = FALSE
  for (i in 1:num_profiles) {
    # assign profile labels and extract appropriate subset
    profile_i <- as.numeric(profiles[i, ])
    data_i <- subset_prof(data_labeled, policy_list, profile_i, inactive)

    # store profile for each observation and subset corresponding to that profile
    data_profile_ids[data_i$id] <- i
    D_profile[i] <- list(data_i$id)

    # if no policies correspond to that profile
    if(nrow(data_i) == 0 | all(is.na(pull(data_i, value)))) {
      if(i == 1){
        control_univ_id = data_i$universal_label[1]
        control_mean = NA
        control_no_data = TRUE
      }
      eq_lb_profiles[i] <- 0
      H_profile <- H_profile + 1
    }

    else {
      eq_lb_profiles[i] <- find_profile_lower_bound(data_i, value)
      if(is.na(eq_lb_profiles[i]) | is.nan(eq_lb_profiles[i])){
        eq_lb_profiles[i] = 0
      }
      if(i == 1){
        control_univ_id = data_i$universal_label[1]
        control_mean = mean(pull(data,value), na.rm = TRUE)
      }
    }
  }

  eq_lb_profiles <- eq_lb_profiles / num_data
  eq_lb_sum <- sum(eq_lb_profiles)

  if(control_no_data) {
    control_loss = 0
  }
  else{
    control_loss <- eq_lb_profiles[[1]] + reg
  }


  # deal with control separately
  control_dict = collections::dict(keys = as.integer(control_univ_id), items = control_mean)
  rashomon_profiles[1] <- list(new_RashomonSet(
    models = list(NA),
    losses = control_loss,
    num_pools = list(1),
    profiles = list(as.numeric(profiles[1, ])),
    pool_dictionaries = list(control_dict)
  ))

  for (i in 2:num_profiles) {
    # extracting relevant things to find rashomon set for this profile
    data_i <- data_labeled[unlist(D_profile[i]), ]

    if(nrow(data_i) == 0 | all(is.na(pull(data_i, value))) | all(is.nan(pull(data_i, value)))) {
      rashomon_profiles[i] <- list(new_RashomonSet(
        models = list(NA),
        losses = 0,
        num_pools = list(0),
        profiles = list(as.numeric(profiles[1, ])),
        pool_dictionaries = list(collections::dict(items = NA, keys = as.integer(-i)))
      ))
      next
    }

    profile_i <- as.numeric(profiles[i, ])
    M_i <- sum(profile_i)
    R_i <- R[as.logical(profile_i)]
    # find profile_lower_bound to find theta_k
    # lower_bound = compute_mse_loss(data,{{value}},)
    theta_k <- theta - (eq_lb_sum - eq_lb_profiles[i])

    data_i <- assign_policy_label(data_i, arm_cols)
    policy_list_i <- create_policies_from_data(data_i, arm_cols) # NEED TO CHANGE
    policy_list_i_masked <- lapply(policy_list_i, function(x) x[as.logical(profile_i)])
    means_i <- policy_means(data_i, value)


    rashomon_i <- find_rashomon_profile(data_i,
      value = value,
      arm_cols = arm_cols,
      M = M_i,
      R = R_i,
      H = H_profile,
      reg = reg,
      profile = profile_i,
      policies = policy_list_i_masked,
      policy_means = means_i,
      normalize = num_data,
      theta = theta_k,
      filtered = TRUE,
      inactive = 0
    )

    if(is.null(rashomon_i)){
      cat(paste0("No models in the RashomonSet for profile: ", profile_i, "\n"))
      return(list())
    }
    rashomon_i = sort_rashomon(rashomon_i)
    rashomon_profiles[i] <- list(rashomon_i)
  }

  if(length(rashomon_profiles) < num_profiles){
    print("No models in the rashomon set")
    return(list())
  }
  R_set <- find_feasible_combinations(rashomon_profiles, theta, H, sorted = TRUE)

  rset = list(R_set, rashomon_profiles)
  rset
}
