#' Creates all factor combinations in a factorial design
#'
#' @param levels A number indicating the number of levels for each arm, or a list indicating
#' the number of levels for each arm
#' @param arms The number of arms, if levels is an integer.
#' @returns A list of all factor combinations (stored as numeric vectors) for a given number of arms and
#' levels per arm.
#' @export
create_policies <- function(levels, arms = 0){

  if(length(levels) == 1){
    arm1_levels = seq.int(1,levels)
    level_list = rep(list(arm1_levels),arms)
  }
  else{
    level_list = lapply(levels, function(x) seq.int(1,x))
  }

  df = expand.grid(level_list)
  asplit(as.matrix(df),1)

}

#'Calculates edgelist for Hasse diagram given sigma and a list of factor combinations (policies).
#'
#'
#' @param sigma M X R-1 matrix that gives where cuts occur (an entry if 1 is no cut, and 0 if there is)
#' @param policy_list A list of all factorial combinations given by create_policies
#' @returns Edge list of Hasse diagram, specifying which policies are connected in the same pool
#' @export
lattice_edges <- function(sigma, policy_list){

  num_policies = length(policy_list)
  edges = list()

  for(i in 1:(num_policies - 1)){
    policy_1 = policy_list[[i]]

    for(j in (i+1):num_policies){

      policy_2 = policy_list[[j]]
      diff_vec = abs(policy_1 - policy_2)
      diff = sum(diff_vec) #calculating difference between abs of policy

      if(diff == 1){ #can only be connected if diff == 1

        arm_diff = which(diff_vec == 1)
        min_dose = min(policy_1[arm_diff],policy_2[arm_diff])

        if(sigma[arm_diff, min_dose] == 1){
          edges[length(edges) + 1] = list(c(i,j)) #adds to edge list
        }
      }
    }
  }
  edges
}

#' Prunes adjacencies given by lattice_adjacencies if our pooling structure (sigma) changes
#' in order to reduce computation when exploring all sigmas
#' @param sigma M X R-1 matrix that gives where cuts occur (an entry if 1 is no cut, and 0 if there is)
#' @param edges Edge list of Hasse diagram, as given by the function lattice_adjacencies
#' @param policy_list A list of all factorial combinations given by create_policies

#' @export
prune_edges <- function(sigma, edges, policy_list){
  new_edges = list()

  for(x in edges){

    pol1 = policy_list[[x[1]]]
    pol2 = policy_list[[x[2]]]
    diff = pol1 - pol2
    arm_diff = which(diff != 0)
    min_dose = min(pol1[arm_diff],pol2[arm_diff])

    if(sigma[arm_diff, min_dose] == 1){
      new_edges[length(new_edges) + 1] = list(c(x[[1]],x[[2]]))
    }
  }
  new_edges
}


#' helper function for use in union_find data structure
merge_components <- function(parent, x){
  if(parent[x] == x){
    return(x)
  }
  return(merge_components(parent, parent[x]))
}

#' Returns the pooling structure given a sigma, in the form of a list that gives
#' the pool number for each factorial combination (policy id)
#' @param n The initial number of parent nodes, usually the total number of factorial combinations
#' @param edges The edgelist of the Hasse diagram given by lattice_edges or prune_edges
#' @returns A collections:dict() that gives the pool for each factorial combination (policy id). The key
#' for each element of the list is the index of its occurrence in policy_list (given by create policies),
#' and its value is the pool it is assigned to.
#' @export
#' @importFrom collections dict
connected_components <- function(n, edges,policy_list){
  parent = seq.int(from = 1, to = n)

  for(x in edges){
    comp1 = merge_components(parent, x[1])
    comp2 = merge_components(parent, x[2])
    parent[comp1] = comp2
  }

  # merge all components to
  for(i in 1:n){
    parent[i] = merge_components(parent, parent[i])
  }

  all_cc = collections::dict()
  for(i in 1:n){
    all_cc$set(as.numeric(unname(policy_list[[i]])),parent[i])
  }
  all_cc

}
