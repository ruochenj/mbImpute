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
#' @param k A scalar corresponding to the number of nearest taxa in a phylogenetic tree we will use to impute a missing value. Theoretically, the larger k, the more accurate our imputation
#' will be (required that k <= m).
#' @param parallel A boolean indicating whether to use parallel or not. Default is FALSE.
#' @param ncores A scalar corresponding to the number of cores to use. It is only used when parallel = TRUE. Default is 1.
#' @param unnormalized A boolean indicating whether the data has been normalized. Details can be found in the Vignette.
#' @return A list three imputed OTU matrices.
#' imp_count_mat_lognorm: imputed normalized and log transformed matrix.
#' imp_count_mat_norm : imputed normalized count matrix with library size of each sample equal to 10^6.
#' imp_count_mat_origlibsize: imputed countmatrix at the original library size.
#' We recommend to use the first one imp_count_mat_lognorm. So that the log scaled counts follows approximately Normal for each taxon across samples.
#' @export
mbImpute <- function(condition = NULL, otu_tab = NULL, meta_data = NULL, D = NULL, k = 5, parallel = FALSE, ncores = 1, unnormalized = T){
  #separate the samples of the dataset into different conditions
  #process the otu_table
  set.seed(1)
  if(is.null(otu_tab)){
    print("No OTU table as input")
    return(NULL)
  }
  if(is.null(condition)){
    condition = rep(1, dim(otu_tab)[1])
  }
  if(is.null(meta_data)){
    print("Meta data information unavailable")
    meta_data <- cbind(rnorm(dim(otu_tab)[1]), rnorm(dim(otu_tab)[1]))
  }else{
    cond <- as.numeric(as.factor(condition))
    dm <- data.matrix(meta_data)
    na_idx <- which(is.na(dm), arr.ind = TRUE)
    if(length(unique(na_idx[,2])) != 0){
      meta_data <- meta_data[,-unique(na_idx[,2])]
    }
    meta_data <- apply(dm, 2, FUN = function(x){
      (x - min(x))/sd(x)
    })
    if(dim(meta_data)[2] == 0){
      print("too many NA's in Meta data")
      meta_data <- cbind(rnorm(dim(otu_tab)[1]), rnorm(dim(otu_tab)[1]))
    }
    cond_check <- cor(dm, cond)
    if(length(which(cond_check == 1)) != 0){
      meta_data <- meta_data[,-which(cond_check == 1)]
    }

  }
  if(is.null(D)){
    print("Phylogenentic information unavailable")
    D = matrix(2, ncol = dim(otu_tab)[2], nrow = dim(otu_tab)[2])
    diag(D) <- 1
  }
  if(unnormalized){
    scale <- rowSums(otu_tab) / (10^6)
    otu_tab <- otu_tab / rowSums(otu_tab) * 10^6
    otu_tab <- log10(otu_tab + 1.01)
  }
  cond_set <- unique(condition)
  for(cond in cond_set){
    idx <- which(condition %in% cond)
    print(paste("condition ", cond, " is imputing", sep = ""))
    otu_tab_cond <- otu_tab[idx,]
    if(is.vector(meta_data)){
      meta_data_cond <- meta_data[idx]
    }else{
      meta_data_cond <- meta_data[idx,]
    }
    otu_tab[idx,] <- data_fit2(otu_tab_cond, meta_data_cond, D, k = k, unnormalized = unnormalized)
  }
  print("Finished.")
  imp_count_mat_norm <- floor(10^(otu_tab) - 1.01)
  return(list(imp_count_mat_lognorm = otu_tab, imp_count_mat_norm = imp_count_mat_norm, imp_count_mat_origlibsize = floor(imp_count_mat_norm * scale)) )
}
