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

#calculates edgelist for Hasse diagram given
#sigma (which is the cuts) along with a list of policys.
#input: sigma, M x R-1 matrix that gives where cuts occur (1 if no cut)
#policy_list, a list of policies
#output: edge list of hasse diagram, specifying which policies are connected
#under pooling scheme implied by sigma.

#' @export
lattice_adjacencies <- function(sigma, policy_list){

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


#prunes adjacencies when given adjacency and new sigma, reduces recalculation
#same basic idea as adjacency_list

#' @export
prune_adjacencies <- function(sigma, edges, policy_list){
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


# make union find structure, so its easy to maintain pools

merge_components <- function(parent, x){
  if(parent[x] == x){
    return(x)
  }
  return(merge_components(parent, parent[x]))
}

#' @export
connected_components <- function(n, edges){
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
    all_cc$set(i,parent[i])
  }
  all_cc

}
