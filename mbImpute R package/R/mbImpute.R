#' mbImpute: an accurate and robust imputation method for microbiome data
#'
#' @description mbImpute aims to impute likely false zero counts or low counts for a normalized OTU table of dimension n * m. Where n is the number of biological samples and m is the number
#' of taxa. In order to achieve our goal, we will borrow information from covariate matrix of dimension n * p and taxa distance matrix of dimension m * m. The final result of the
#' mbImpute function will ouput an imputed matrix of dimension n * m, which is exactly the same dimension as the input.
#' @param condition A vector of length n indicating which group (e.g. treatment/control, CRC/control/Adenoma) the biology samples belong to.
#' @param otu_tab A matrix of dimension n * m corresponding to OTU count table.
#' @param meta_data A matrix of dimension n * p corresponding to the covariate matrix/meta data comes with the biological sample. Please note that the conditions/groups of the samples
#' should be excluded from the meta data.
#' @param D A matrix of dimension m * m corresponding to the distance matrix of taxa generated from an edge set. If an edge set is provided, you can use distance_mat function in our
#' package to generate the distance matrix D.
#' @param psi A scalar corresponding to the importance weight of the taxa distance matrix. This arguement is only useful if cross_valid is FALSE.
#' @param k A scalar corresponding to the number of nearest taxa in a phylogenetic tree we will use to impute a missing value. Theoretically, the larger k, the more accurate our imputation
#' will be (required that k <= m).
#' @param cross_valid A boolean indicating whether we will use cross validation to find the best psi. If cross_valid is TRUE, we use grid search and find the best psi in the set 0.5, 1, 1.5,
#' 2.
#' @param parallel A boolean indicating whether to use parallel or not. Default is FALSE.
#' @param ncores A scalar corresponding to the number of cores to use. It is only used when parallel = TRUE. Default is 1.
#' @param type_1 A boolean indicating which type of mbImpute to use. Details can be found in the Vignette.
#' @return An imputed and normalized OTU table at log scale. User can rescale it back to counts by using \code{floor(10^(otu_tab) - 1.01)}.
#' @export
mbImpute <- function(condition = NULL, otu_tab = NULL, meta_data = NULL, D = NULL, psi = 2, k = 10, cross_valid = FALSE, parallel = FALSE, ncores = 1, type_1 = T){
  #separate the samples of the dataset into different conditions
  #process the otu_table
  set.seed(1)
  otu_tab <- otu_tab / rowSums(otu_tab) * 10^6
  otu_tab <- log10(otu_tab + 1.01)
  if(is.null(otu_tab)){
    print("No OTU table as input")
    return(NULL)
  }
  if(is.null(condition)){
    condition = rep(1, dim(otu_tab)[1])
  }
  if(is.null(meta_data)){
    print("Meta data information unavailable")
    meta_data <- cbind(rnorm(dim(y_imp)[1]), rnorm(dim(y_imp)[1]))
  }else{
    cond <- as.numeric(as.factor(condition))
    dm <- data.matrix(meta_data)
    na_idx <- which(is.na(dm), arr.ind = TRUE)
    if(length(unique(na_idx[,2])) != 0){
      meta_data <- meta_data[,-unique(na_idx[,2])]
    }
    cond_check <- cor(meta_data, cond)
    if(length(which(cond_check == 1)) != 0){
      meta_data <- meta_data[,-which(cond_check == 1)]
    }
  }
  if(is.null(D)){
    print("Phylogenentic information unavailable")
    D = matrix(1, ncol = dim(y_imp)[2], nrow = dim(y_imp)[2])
  }

  cond_set <- unique(condition)
  for(cond in cond_set){
    idx <- which(condition %in% cond)
    print(paste("condition ", cond, " is imputed", sep = ""))
    otu_tab_cond <- otu_tab[idx,]
    if(is.vector(meta_data)){
      meta_data_cond <- meta_data[idx]
    }else{
      meta_data_cond <- meta_data[idx,]
    }
    if(!parallel){
      if(cross_valid){
        psi_vec = c(0.5, 1, 1.5, 2)
        min_mse = 1000
        for(psi in psi_vec){
          impute_val <- data_fit2(otu_tab_cond, meta_data_cond, D, psi, k = k, type_1 = type_1)
          if(impute_val$mse < min_mse){
            otu_tab[idx,] = impute_val$y_imp
            min_mse = impute_val$mse
          }
        }
      }else{
        impute_val <- data_fit2(otu_tab_cond, meta_data_cond, D, psi, k = k, type_1 = type_1)
        print("impute_val generated")
        print(dim(impute_val$y_imp))
        otu_tab[idx,] = impute_val$y_imp
        #rm(impute_val)
      }
    }else{
      if(cross_valid){
        psi_vec = c(0.5, 1, 1.5, 2)
        min_mse = 1000
        for(psi in psi_vec){
          impute_val <- data_fit2(otu_tab_cond, meta_data_cond, D, psi, k = k, parallel = parallel, ncores = ncores, type_1 = type_1)
          if(impute_val$mse < min_mse){
            otu_tab[idx,] = impute_val$y_imp
            min_mse = impute_val$mse
          }
        }
      }else{
        impute_val <- data_fit2(otu_tab_cond, meta_data_cond, D, psi, k = k, parallel = parallel, ncores = ncores, type_1 = type_1)
        print("impute_val generated")
        print(dim(impute_val$y_imp))
        otu_tab[idx,] = impute_val$y_imp
        #rm(impute_val)
      }
    }
  }
  return(otu_tab)
}
