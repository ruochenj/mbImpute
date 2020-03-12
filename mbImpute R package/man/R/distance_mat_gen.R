#' distance_mat_gen
#'
#' @description This function aims to generate the distance matrix of taxa based on their path lengths in the phylogenetic tree.
#' @param edges A matrix of dimension N * 2 corresponding to the edge set of the phylogenetic tree (similar to the edge set for a graph).
#' @param n_taxa A scalar corresponding to number of taxa in the dataset.
#' @param tree_height A scalar corresponding to the height of the phylogenetic tree. Any number larger than the height of the phylogenetic tree will also work.
#' The initial value is set to be 50, which is usually enough as for a complete binary tree, height of 50 corresponds to 2^50 nodes.
#' @return A matrix of dimension n_taxa * n_taxa corresponding to the matrix D in mbImpute function.
#' @export
distance_mat_gen <- function(edges, n_taxa, tree_height = 50){
  k = tree_height
  m = n_taxa
  nd_mat <- matrix(rep(1, k*m), k, m)
  l <- rep(1,k)

  for(i in 1:n_taxa){
    print(i)
    l <- rep(1,tree_height+1)
    l[1] = i
    for(j in 2:(tree_height+1)){
      if(sum(edges[,2] %in% l[j-1]) != 0){
        l[j] = edges[edges[,2] %in% l[j-1], 1]
      }
      else{
        l[j] = NA
      }
    }
    nd_mat[,i] = l[2:(tree_height+1)]
  }


  d1_mat <- matrix(0, nrow = n_taxa, ncol = n_taxa)
  #records the position of 1:n_taxa in the edges set.
  taxa_vec <- match(1:n_taxa, edges[,2])
  #generate the distance matrix
  for(i in 1:n_taxa){
    for(j in 1:n_taxa){
      int_sc <- intersect(nd_mat[,i], nd_mat[,j])
      leni <- sum(!is.na(int_sc))
      len1 <- sum(!is.na(nd_mat[,i]))
      len2 <- sum(!is.na(nd_mat[,j]))
      d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
    }
  }
  diag(d1_mat) = 0
  #d1_mat denotes the distance for two taxa
  return(d1_mat)
}

