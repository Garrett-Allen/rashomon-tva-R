#' Finds sum of all the products of size k for each k < n, where n is the size
#' of the list. Useful for counting the number of pools.
#'
#' @param arr A list
#' @returns A list of the sums of the products of size k, where the (k+1)th element
#' denotes the kth sum. The 1st element is 1, for use in other functions for easy looping
#' @noRd
sum_product_k <- function(arr) {
  # Length of the input array
  n <- length(arr)

  # Initialize cache and k_sum arrays
  cache <- numeric(n + 1)
  k_sum <- numeric(n + 1)

  k_sum[1] <- 1

  # Case k = 1
  for (i in 1:n) {
    cache[i + 1] <- arr[i]
    k_sum[2] <- k_sum[2] + arr[i]
  }

  # Case k = 2 through n
  for (k in 3:(n + 1)) {
    prev_sum <- k_sum[k - 1]

    for (i in 1:n) {
      prev_sum <- prev_sum - cache[i + 1]
      cache[i + 1] <- arr[i] * prev_sum
    }

    k_sum[k] <- sum(cache)
  }

  return(k_sum)
}

#' Finds R, a list (or integer) of the number of levels in each arm.
#'
#' @param sigma Partition matrix for a given pooling structure
#' @returns A list (or integer) of the number of levels in each arm.
#' @export
find_R <- function(sigma) {
  R <- rowSums(sigma, na.rm = TRUE) + 2

  if (all(R[1] == R)) {
    return(R[1])
  }
  return(R)
}


#' Finds the number of pools in a sigma with the same number of levels for each arm.
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param R The number of total numbers of levels for each arm (in this function,
#' we are assuming R is equal for all arms)
#' @returns The number of pools given by this pooling structure
#' @noRd
num_pools_fixed_r <- function(sigma, R) {
  m <- nrow(sigma)
  z <- rowSums(sigma)
  z_sums <- sum_product_k(z)
  H <- 0
  for (i in 1:(m + 1)) {
    sign <- (-1)^(i - 1)
    H <- H + sign * z_sums[i] * ((R - 1))^(m - i + 1)
  }
  H
}

#' Finds the number of pools in a sigma with differing numbers of levels for each arm.
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param R Vector of the number of levels for each arm
#' @returns The number of pools given by this pooling structure
#' @noRd
num_pools_change_r <- function(sigma, R) {
  m <- nrow(sigma)
  R <- R - 1
  R_prod <- prod(R)

  z <- rowSums(sigma, na.rm = TRUE)
  indices <- 1:m
  z_combs <- powerset(indices)

  # accounting for null case (subset of size 0)
  H <- R_prod

  # due to R quirks, I split this up into subsets of size {1,..., n-1} and
  # subset of size n, but its just the same code repeated.

  for (k in 1:length(z_combs)) {
    subsets_k <- z_combs[[k]]

    if (!is.null(nrow(subsets_k))) {
      for (i in 1:ncol(subsets_k)) {
        subset_ik <- subsets_k[, i]

        sign <- (-1)**length(subset_ik)
        z_sum <- prod(z[as.vector(subset_ik)])
        splits <- R_prod / prod(R[as.vector(subset_ik)])

        H <- H + sign * z_sum * splits
      }
    } else {
      subset_ik <- subsets_k

      sign <- (-1)**length(subset_ik)
      z_sum <- prod(z[as.vector(subset_ik)])
      splits <- R_prod / prod(R[as.vector(subset_ik)])

      H <- H + sign * z_sum * splits
    }
  }

  H
}

#' Finds powerset of a given list (all the possible subsets of a list)
#'
#' @param arr A list
#' @returns A list of lists, where element list[[k]] is a matrix of all of
#' the subsets of size k, such that list[[k]][,i] is the ith subset of size k.
#' @importFrom RcppAlgos comboGeneral
#' @noRd
powerset <- function(arr) {
  subsets <- lapply(1:length(arr), function(r) {
    t(RcppAlgos::comboGeneral(arr, r))
  })
  subsets
}

#' Helper function used to find number of pools given arbitrary sigma and
#' number of levels in each arm.
#'
#' @param sigma Partition matrix for a given pooling structure
#' @param R A list (or integer) of the number of levels in each arm.
#' @returns The number of pools given by the pooling structure sigma.
#' @export
num_pools <- function(sigma, R = NA) {
  if (all(is.na(sigma))) {
    return(1)
  }

  if (all(is.na(R))) {
    R <- find_R(sigma)
  }

  if (length(R) == 1) {
    H <- num_pools_fixed_r(sigma, R)
  } else {
    H <- num_pools_change_r(sigma, R)
  }

  H
}
