setwd("~/Dropbox/mbimpute/code/Real data Karlsson")


###### simulation ##########
# mbImpute coefficients
control_coeff <- readRDS("02072020_mbImpute/control_dat1_sim_add_filter_coef.rds")
T2D_coeff <- readRDS("02072020_mbImpute/T2D_dat2_sim_add_filter_coef.rds")

otu_real_data_T2D <- read.csv("otu_real_data_T2D.csv", row.names = "X")
otu_real_data_control <- read.csv("otu_real_data_control.csv", row.names = "X")

D <- read.csv("D.csv", row.names = "X")

meta_data_T2D <- read.csv("meta_data_T2D.csv", row.names = "X")
meta_data_control <- read.csv("meta_data_control.csv", row.names = "X")

set.seed(1234)
y_sim = otu_real_data_T2D
x = meta_data_T2D
k = 5
c1 = T2D_coeff
c1 = c1*5
c1[1] = T2D_coeff[1]
c1[length(c1)-6] = 0.2
c1[length(c1)-11] = -0.3

#loading used functions
design_mat_row_gen2 <- function(count_mat, covariate_mat, row_index, col_index, close_taxa){
  n = dim(count_mat)[1]
  m = dim(count_mat)[2]
  k = length(close_taxa[[1]])
  if(is.vector(covariate_mat)){
    p = 1
    i = row_index
    j = col_index

    close_taxa_set <- close_taxa[[j]]
    #generate a row including response and a row of design matrix.
    row_gen <- rep(0, m*k + (n-1) * n + n*p)

    row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
    row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
    row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i])
  }else{
    p = dim(covariate_mat)[2]
    i = row_index
    j = col_index

    close_taxa_set <- close_taxa[[j]]
    #generate a row including response and a row of design matrix.
    row_gen <- rep(0, m*k + (n-1) * n + n*p)

    row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
    row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
    row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i,])
  }
  return(row_gen)
}

#identify the low abundance taxa by binomial test
filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
  y = y_sim[,col_i]
  n = length(y)
  nz <- sum(y <= (log10(1.01) + 1e-6))
  pz = 1 - nz/n
  test = pz - 1.96 * sqrt(pz * (1-pz)/n)
  if(nz == n || test <= 0){
    return(0)
  }else{
    return(1)
  }
})

#keep the original matrix
#y_imp <- y_sim
# print(dim(y_imp))
# #keep the vectors that we neither impute nor borrow information
# remain_vec <- which(unlist(filter) == 0)
# y_rem <- y_sim[, remain_vec]
filter_vec <- which(unlist(filter) == 1)
y_sim = y_sim[, filter_vec]
D = D[filter_vec,filter_vec]
y_preserve <- y_sim
#apply the imputation method on the simulated matrix
m = dim(y_sim)[2]
n = dim(y_sim)[1]
# 53 * 193

#x[1:50,]
#X <- cbind(rep(1,n), x[1:n,])
X <- as.matrix(x)

#identifying the group needs imputation and group doesn't need imputation

#find k closest taxa
# k = 5
#generate close set
close_taxa <- list()
for(j in 1:m){
  #close_taxa[[j-1]] = which(D[,j] %in% sort(unique(D[,j]))[2:(k+1)])
  close_dist <- D[D[,j] %in% sort(D[,j])[2:(k+1)],j]
  close_taxa_vec = which(D[,j] %in% close_dist[close_dist != max(close_dist)])
  if(length(close_taxa_vec) < k){
    close_taxa_vec <- c(close_taxa_vec, which(D[,j] == max(close_dist))[1:(k-length(close_taxa_vec))])
  }
  close_taxa[[j]] = close_taxa_vec
}
#generate design matrix
p = dim(x)[2]
if(is.null(p)){
  p = 1
}
#generate a row including response and a row of design matrix.
row_length <- m * k + (n-1) * n + n*p
idx_set <- matrix(NA, m*n, 2)
for(i in 1:n){
  for(j in 1:m){
    idx_set[(i-1)*m+j, ] <- c(i, j)
  }
}


# for(itr in 1:5){
#   print(itr)
  #perform imputation on the rest

  design_mat_gen <- matrix(0, nrow = dim(idx_set)[1], ncol = row_length)
  for(i in 1:dim(idx_set)[1]){
    design_mat_gen[i,] <- design_mat_row_gen2(y_sim, x[1:n,], idx_set[i,1], idx_set[i,2], close_taxa)
  }
  design_mat_gen <- cbind(rep(1, dim(design_mat_gen)[1]), design_mat_gen)
  imputed_value <- design_mat_gen %*% c1
  impute_mat <- y_sim
  for(i in 1:dim(idx_set)[1]){
    impute_mat[idx_set[i,1], idx_set[i,2]] = max(imputed_value[i], log10(1.01))
  }
  print(sqrt(sum((impute_mat - y_sim)^2)))
  y_sim <- impute_mat
# }

y_sim <- y_sim - 1.5
write.csv(y_sim, "Karlsson_simulated.csv")

gamma_norm_mix <- function(y, X){
  loglik <- function(p, alpha, beta, cov_par, var1, X, y){
    n = length(y)
    lkval <- 0
    fgam <- dgamma(y, shape = alpha, rate = beta)
    for(i in 1:n){
      if(!is.vector(X)){
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i,] %*% cov_par, sd = sqrt(var1)))
      }else{
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i] * cov_par, sd = sqrt(var1)))
      }

    }
    return(lkval)
  }
  n = length(y)
  alpha_init <- 1
  beta_init <- 10
  p_init <- 0.5
  cov_par_init <- solve(t(X) %*% X) %*% t(X) %*% y
  var_init <- t(y - X %*% cov_par_init) %*% (y - X %*% cov_par_init) / n

  alpha_t <- alpha_init
  beta_t <- beta_init
  cov_par_t <- cov_par_init
  var_t <- var_init
  p_t <- p_init

  #update gamam param
  #Wei's Method
  ### root-finding equation
  fn = function(alpha, target){
    log(alpha) - digamma(alpha) - target
  }
  update_gmm_pars = function(x, wt){
    if(max(wt) > 0.00001){
      tp_s = sum(wt)
      tp_t = sum(wt * x)
      tp_u = sum(wt * log(x))
      tp_v = -tp_u / tp_s - log(tp_s / tp_t)
      if (tp_v <= 0){
        alpha = 20
      }else{
        alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
        if (alpha0 >= 20){alpha = 20
        }else{
          alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v,
                          extendInt = "yes")$root
        }
      }
      ## need to solve log(x) - digamma(x) = tp_v
      ## We use this approximation to compute the initial value
      beta = tp_s / tp_t * alpha
    }else{
      alpha = 0.001
      beta = 1000
    }

    return(c(alpha, beta))
  }
  #convergence criteria
  flag = TRUE
  maxitr = 300
  itr = 0
  while(flag){
    #E_step
    mean_t <- X %*% cov_par_t
    n = length(y)
    a_hat_t <- rep(0, n)
    dg_t <- dgamma(y, shape = alpha_t, rate = beta_t)
    for(i in 1:n){
      if(dg_t[i] == 0){
        a_hat_t[i] = 0
      }else{
        a_hat_t[i] <- p_t * dg_t[i]/(p_t * dg_t[i]+ (1-p_t)*dnorm(y[i], mean = mean_t[i], sd = sqrt(var_t)))
      }
    }
    #maximization
    #fit p
    p_t1 <- sum(a_hat_t)/n
    X_tilta <- sqrt(1-a_hat_t) * X
    y_tilta <- sqrt(1-a_hat_t) * y
    #fit normal
    out <- tryCatch(
      {
        # Just to highlight: if you want to use more than one
        # R expression in the "try" part then you'll have to
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression
        # in case the "try" part was completed successfully
        cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
        TRUE
      },
      error=function(cond) {
        FALSE
      }
    )
    if(!out){
      return(list("d" = y < log10(1.01) + 10^(-3) ))
    }
    cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
    var_t1 <- sum((1 - a_hat_t) * (y - X %*% cov_par_t)^2) / sum(1-a_hat_t)
    #fit gamma
    par_gamm <- update_gmm_pars(x = y, wt = a_hat_t)
    alpha_t1 = par_gamm[1]
    beta_t1 <- par_gamm[2]
    loglik1 <- loglik(p = p_t, alpha = alpha_t, beta = beta_t, cov_par = cov_par_t, var1 = var_t, X = X, y = y)
    loglik2 <- loglik(p = p_t1, alpha = alpha_t1, beta = beta_t1, cov_par = cov_par_t1, var1 = var_t1, X = X, y = y)
    if((abs(loglik1 - loglik2)) < 0.05 || itr > maxitr){
      flag = FALSE
    }else{
      alpha_t <- alpha_t1
      beta_t <- beta_t1
      cov_par_t <- cov_par_t1
      var_t <- var_t1
      p_t <- p_t1
      itr = itr + 1
    }
  }
  #fit a normal curve
  eta_hat <- cov_par_init
  omega_hat <- var_init
  #calculate new likelihood
  norm_log_lik = 0
  for(i in 1:n){
    if(!is.vector(X)){
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i,] %*% eta_hat, sd = sqrt(omega_hat)))
    }else{
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i] * eta_hat, sd = sqrt(omega_hat)))
    }
  }
  Dev = -2 * norm_log_lik - (-2 * loglik2)
  judgement = pchisq(Dev, df = 3, lower.tail = FALSE) < 0.05

  if(!judgement || (alpha_t/beta_t > 1)){
    p_t = 0
    a_hat_t <- rep(0,n)
  }
  return(list("p" = p_t, "alpha" = alpha_t, "beta" = beta_t, "cov_par" = cov_par_t, "var" = var_t, "d" = a_hat_t, "eta" = eta_hat, "omega" = omega_hat, "Deviance" = Dev))
}

####### learn from learning missing rate empirically ############

otable <- y_preserve
meta_tab <- meta_data_T2D
D

dim(otable)
# read in the gamma_norm_mix function
mean_record <- c()
percentage_record <- c()
D_vale_record <- c()
beta_record <- list()
for(j in 1:dim(otable)[2]){
  result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab))
  beta_record[[j]] = result$cov_par
  mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
  percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
  D_vale_record <- c(D_vale_record, result$d)
}
#filter <- which(is.nan(mean_record) | mean_record < log10(1.01) + 0.001)

# filter out the nan values.
plot(mean_record, percentage_record)

# build a map between the mean and percentage of missing
missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/3
  mean_up <- mean + win_len
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}
# missing_rate(1, mean_record, percentage_record)

col_mean <- colMeans(y_sim)
zero_rate <- unlist(lapply(col_mean, FUN = function(x){
  print(x)
  return(missing_rate(x, mean_record, percentage_record))
}))
zero_mat <- matrix(NA, nrow = n, ncol = m)
for(i in 1:m){
  zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
}
sim_tab_zi = y_sim * zero_mat
sim_tab_zi[sim_tab_zi < log10(1.01)+1e-6] = log10(1.01)

sqrt(sum((sim_tab_zi - y_preserve)^2))
sim_tab_zi_mbImpute <- sim_tab_zi

write.csv(sim_tab_zi, "simulated_zi_matrix_mbImpute.csv")

sim_tab_zi <- 10^(sim_tab_zi) - 1.01
sim_tab_zi_trans <- t(sim_tab_zi)
write.csv(sim_tab_zi_trans, "simulated_zi_matrix.csv")

y_sim_rec <- y_sim
sim_tab_zi_rec <- sim_tab_zi

write.csv(y_sim_rec, "truth.csv")

print(max(colMeans(y_sim)))
print(max(colMeans(sim_tab_zi)))
print(max(rowMeans(sim_tab_zi_trans)))
###### Library and functions ######
library(scImpute)
library(ggplot2)
library(gplots)
library(SAVER)
library(softImpute)

# alra
randomized.svd <- function(A,K, q, method = 'rsvd', mkl.seed = -1) {
  out <- setNames(vector("list", 3), c("u", "d", "v"))
  if (method == 'rsvd') {
    library(rsvd)
    out <- rsvd(A,K,q=q)
  }else if (method == 'rsvd-mkl') {
    library(fastRPCA)
    fastPCAOut <- fastPCA(A, k=K, its=q, l=(K+10), seed=mkl.seed)
    out$u <- fastPCAOut$U
    out$v <- fastPCAOut$V
    out$d <- diag(fastPCAOut$S)
  }else{
    stop('Method not recognized')
  }
  return(out)
}

normalize_data <- function (A) {
  #  Simple convenience function to library and log normalize a matrix

  totalUMIPerCell <- rowSums(A);
  if (any(totalUMIPerCell == 0)) {
    toRemove <- which(totalUMIPerCell == 0)
    A <- A[-toRemove,]
    totalUMIPerCell <- totalUMIPerCell[-toRemove]
    cat(sprintf("Removed %d cells which did not express any genes\n", length(toRemove)))
  }

  A_norm <- sweep(A, 1, totalUMIPerCell, '/');
  A_norm <- A_norm * 10E3
  A_norm <- log(A_norm +1);
}

choose_k <- function (A_norm,K=100, thresh=6, noise_start=80,q=2,use.mkl=F, mkl.seed =-1) {
  #  Heuristic for choosing rank k for the low rank approximation based on
  #  statistics of the spacings between consecutive singular values. Finds
  #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
  #  is significantly different than spacings in the tail of the singular values.
  #
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   K: Number of singular values to compute. Must be less than the
  #   smallest dimension of the matrix.
  #   thresh: Number of standard deviations away from the ``noise'' singular
  #   values which you consider to be signal
  #   noise_start : Index for which all smaller singular values are
  #   considered noise
  #   q : Number of additional power iterations
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  # Returns:
  #   A list with three items
  #       1) Chosen k
  #       2) P values of each possible k
  #       3) Singular values of the matrix A_norm

  if (K > min(dim(A_norm))) {
    stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
  }
  if (noise_start >K-5) {
    stop("There need to be at least 5 singular values considered noise.\n")
  }
  noise_svals <- noise_start:K
  if (!use.mkl) {
    rsvd_out <- randomized.svd(A_norm,K,q=q)
  }else {
    rsvd_out <- randomized.svd(A_norm,K,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
  diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]
  mu <- mean(diffs[noise_svals-1])
  sigma <- sd(diffs[noise_svals-1])
  num_of_sds <- (diffs-mu)/sigma
  k <- max (which(num_of_sds > thresh))
  return (list( k=k,num_of_sds = num_of_sds,d=rsvd_out$d))
}

alra <- function( A_norm, k=0,q=10, quantile.prob = 0.001, use.mkl = F, mkl.seed=-1) {
  # Computes the k-rank approximation to A_norm and adjusts it according to the
  # error distribution learned from the negative values.
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   k : the rank of the rank-k approximation. Set to 0 for automated choice of k.
  #   q : the number of additional power iterations in randomized SVD
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  #
  # Returns:
  #   A list with three items
  #       1) The rank k approximation of A_norm.
  #       2) The rank k approximation of A_norm, adaptively thresholded
  #       3) The rank k approximation of A_norm, adaptively thresholded and
  #       with the first two moments of the non-zero values matched to the
  #       first two moments of the non-zeros of A_norm. This is the completed
  #       matrix most people will want to work with
  # Example:
  #     result.completed <- adjusted_svd(A_norm,15)
  #     A_norm_rank15 <- result.completed[[1]]     # The low rank approximation for reference purposes...not suggested for matrix completion
  #     A_norm_rank15_cor <- result.completed[[3]] # The actual adjusted, completed matrix

  cat(sprintf("Read matrix with %d cells and %d genes\n", nrow(A_norm), ncol(A_norm)))
  if (class(A_norm) != 'matrix') {
    stop(sprintf("A_norm is of class %s, but it should be of class matrix. Did you forget to run as.matrix()?",class(A_norm)))
  }

  if (k ==0 ) {
    k_choice <- choose_k(A_norm)
    k <-  k_choice$k
    cat(sprintf("Chose k=%d\n",k))
  }

  cat("Getting nonzeros\n")
  originally_nonzero <- A_norm >0

  cat("Randomized SVD\n")
  if (!use.mkl) {
    fastDecomp_noc <- randomized.svd(A_norm,k,q=q)
  }else {
    fastDecomp_noc <- randomized.svd(A_norm,k,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
  A_norm_rank_k <- fastDecomp_noc$u[,1:k]%*%diag(fastDecomp_noc$d[1:k])%*% t(fastDecomp_noc$v[,1:k])


  cat(sprintf("Find the %f quantile of each gene\n", quantile.prob))
  #A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,min))
  A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,FUN=function(x) quantile(x,quantile.prob)))
  cat("Sweep\n")
  A_norm_rank_k_cor <- replace(A_norm_rank_k, A_norm_rank_k <= A_norm_rank_k_mins[col(A_norm_rank_k)], 0)


  sd_nonzero <- function(x) sd(x[!x == 0])
  sigma_1 <- apply(A_norm_rank_k_cor, 2, sd_nonzero)
  sigma_2 <- apply(A_norm, 2, sd_nonzero)
  mu_1 <- colSums(A_norm_rank_k_cor)/colSums(!!A_norm_rank_k_cor)
  mu_2 <- colSums(A_norm)/colSums(!!A_norm)

  toscale <- !is.na(sigma_1) & !is.na(sigma_2) & !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)

  cat(sprintf("Scaling all except for %d columns\n", sum(!toscale)))

  sigma_1_2 <- sigma_2/sigma_1
  toadd  <- -1*mu_1*sigma_2/sigma_1 + mu_2

  A_norm_rank_k_temp <- A_norm_rank_k_cor[,toscale]
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, sigma_1_2[toscale],FUN = "*")
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, toadd[toscale],FUN = "+")

  A_norm_rank_k_cor_sc <- A_norm_rank_k_cor
  A_norm_rank_k_cor_sc[,toscale] <- A_norm_rank_k_temp
  A_norm_rank_k_cor_sc[A_norm_rank_k_cor==0] = 0

  lt0 <- A_norm_rank_k_cor_sc  <0
  A_norm_rank_k_cor_sc[lt0] <- 0
  cat(sprintf("%.2f%% of the values became negative in the scaling process and were set to zero\n", 100*sum(lt0)/(nrow(A_norm)*ncol(A_norm))))

  A_norm_rank_k_cor_sc[originally_nonzero & A_norm_rank_k_cor_sc ==0] <- A_norm[originally_nonzero & A_norm_rank_k_cor_sc ==0]

  colnames(A_norm_rank_k_cor) <- colnames(A_norm)
  colnames(A_norm_rank_k_cor_sc) <- colnames(A_norm)
  colnames(A_norm_rank_k) <- colnames(A_norm)

  original_nz <- sum(A_norm>0)/(nrow(A_norm)*ncol(A_norm))
  completed_nz <- sum(A_norm_rank_k_cor_sc>0)/(nrow(A_norm)*ncol(A_norm))
  cat(sprintf("The matrix went from %.2f%% nonzero to %.2f%% nonzero\n", 100*original_nz, 100*completed_nz))

  list(A_norm_rank_k=A_norm_rank_k,A_norm_rank_k_cor =A_norm_rank_k_cor, A_norm_rank_k_cor_sc=A_norm_rank_k_cor_sc)
}

si_impute <- function(sim_tab){
  set.seed(12345)
  y_sim = sim_tab
  # add filter
  filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    n = length(y)
    nz <- sum(y <= (log10(1.01) + 1e-6))
    pz = 1 - nz/n
    test = pz - 1.96 * sqrt(pz * (1-pz)/n)
    if(nz == n || test <= 0){
      return(0)
    }else{
      return(1)
    }
  })
  y_imp <- y_sim
  #perform imputation on the rest
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]

  na_idx <- y_sim < log10(1.01) + 0.001
  y_sim[y_sim < log10(1.01) + 0.001] = NA

  # y_sim_cv <- unlist(y_sim)
  # y_sim_validate <- matrix(y_sim_cv, nrow(y_sim), ncol = ncol(y_sim))
  # identical(y_sim_validate, y_sim)
  y_sim_cv <- unlist(y_sim)
  na_intro <- sample(which(!is.na(y_sim_cv)), floor(sum(!is.na(y_sim_cv))/10))
  y_sim_cv_intro <- y_sim_cv
  y_sim_cv_intro[na_intro] = NA
  y_sim_cv_intro <- matrix(y_sim_cv_intro, nrow = nrow(y_sim), ncol = ncol(y_sim))

  j = 1
  se = 1e10
  for(i in 1:5){
    si_cv_1 <- softImpute(y_sim_cv_intro, rank.max = i, lambda = 0)
    y_imp_cv <- complete(y_sim_cv_intro, si_cv_1)
    y_sim_vali <- as.vector(y_imp_cv)
    se2 <- sum((y_sim_cv[na_intro] - y_sim_vali[na_intro])^2)
    print(se2)
    if(se2 < se){
      se = se2
      j = i
    }
  }
  print(j)
  si1 <- softImpute(y_sim_cv_intro, rank.max = 10, lambda = 0, trace.it = TRUE)
  impute_mat <- complete(y_sim_cv_intro, si1)

  y_imp[, filter_vec] = impute_mat

  return(y_imp)
}

####################

print(max(colMeans(y_sim)))
print(max(colMeans(sim_tab_zi)))
print(max(rowMeans(sim_tab_zi_trans)))

# y_sim_rec_new <- y_sim_rec

eval_omit <- c(26, 34, 42, 117, 123, 158, 191, 192)

real_data <- data.matrix(sim_tab_zi_trans)
dim(real_data)
saver_mat <- saver(real_data, ncores = 1, estimates.only = TRUE)
saver_mat_eval <- log10(saver_mat + 1.01)
saver_mat_eval <- t(saver_mat_eval)
for(j in 1:dim(saver_mat_eval)[2]){
  saver_mat_eval[,j] <- saver_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(saver_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - saver_mat_eval[,-eval_omit])^2))
# 135.5029
# trans 201.5581
# omit:194.1245
max(rowMeans(sim_tab_zi_trans))

sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(saver_mat_eval[,-eval_omit]))^2)
sum((apply(y_sim_rec[,-eval_omit], 2, var) - apply(saver_mat_eval[,-eval_omit], 2, var))^2)


scimpute(count_path = "simulated_zi_matrix.csv", Kcluster = 1, out_dir = "sim_imp")
scImpute_mat <- read.csv("sim_impscimpute_count.csv", row.names = "X")
dim(scImpute_mat)
scImpute_mat_eval <- log10(scImpute_mat + 1.01)
scImpute_mat_eval <- t(scImpute_mat_eval)
for(j in 1:dim(scImpute_mat_eval)[2]){
  scImpute_mat_eval[,j] <- scImpute_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(scImpute_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - scImpute_mat_eval[,-eval_omit])^2))
# 120.8523
# trans 205.3515
# 136.0469
# omit:125.4932
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(scImpute_mat_eval[,-eval_omit]))^2)
sum((apply(y_sim_rec[,-eval_omit], 2, var) - apply(scImpute_mat_eval[,-eval_omit], 2, var))^2)
# 125.4932
# 55.83098

T2D_mat <- normalize_data(sim_tab_zi)
k_chosen <- choose_k(T2D_mat, K = 49, noise_start = 44)$k
T2D_mat <- as.matrix(T2D_mat)
alra_mat_eval <- alra(T2D_mat, k = k_chosen)$A_norm_rank_k_cor_sc
for(j in 1:dim(alra_mat_eval)[2]){
  alra_mat_eval[,j] <- alra_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(alra_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - alra_mat_eval[,-eval_omit])^2))
# 209.4349
# omit:202.26
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(alra_mat_eval[,-eval_omit]))^2)
sum((apply(y_sim_rec[,-eval_omit], 2, var) - apply(alra_mat_eval[,-eval_omit], 2, var))^2)
# 611.5041
# 272.6174

Magic <- read.csv("simulation_MAGIC.csv")
for(j in 1:dim(Magic)[2]){
  Magic[,j] <- Magic[,j] * max(sim_tab_zi_mbImpute[,j]) / max(Magic[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - Magic[,-eval_omit])^2))
# 73.80223
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(Magic[,-eval_omit]))^2)
# 57.65293
sum((apply(y_sim_rec[,-eval_omit], 2, var) - apply(Magic[,-eval_omit], 2, var))^2)
# 9

# mbImpute
mbImpute_mat <- readRDS("sim_imp/T2D_RUNA_imputed2_mat.rds")
filter_vec_mbImpute <- readRDS("sim_imp/T2Dfilter_vec.rds")
coeff <- readRDS("sim_imp/T2D_dat2_sim_add_filter_coef.rds")
sqrt(sum((y_sim_rec[,-eval_omit] - mbImpute_mat[,-eval_omit])^2))
#sqrt(sum((y_sim_rec - mbImpute_mat)^2))
# 95.52454
# omit:78.414
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(mbImpute_mat[,-eval_omit]))^2)
# 57.65293
sum((apply(y_sim_rec[,-eval_omit], 2, var) - apply(mbImpute_mat[,-eval_omit], 2, var))^2)

# softImpute
softImpute_mat = si_impute(sim_tab_zi)
softImpute_mat[softImpute_mat < 0] = 0
softImpute_mat_eval <- log10(softImpute_mat+1.01)
for(j in 1:dim(softImpute_mat_eval)[2]){
  softImpute_mat_eval[,j] <- softImpute_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(softImpute_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - softImpute_mat_eval[,-eval_omit])^2))
# 113.7773

# zero_flated
sqrt(sum((y_sim_rec[,-eval_omit] - sim_tab_zi_mbImpute[,-eval_omit])^2))
# 197.8466


bt_mse <- matrix(ncol = 800, nrow = 7)
set.seed(1234)
for(i in 1:800){
  bt_idx <-  sample(1:185, 185, replace = TRUE)
  bt_mse[,i] <- c(sum((y_sim_rec[,-eval_omit][,bt_idx] - Magic[,-eval_omit][,bt_idx])^2),
  sum((y_sim_rec[,-eval_omit][,bt_idx] - mbImpute_mat[,-eval_omit][,bt_idx])^2),
  sum((y_sim_rec[,-eval_omit][,bt_idx] - scImpute_mat_eval[,-eval_omit][,bt_idx])^2),
  sum((y_sim_rec[,-eval_omit][,bt_idx] - alra_mat_eval[,-eval_omit][,bt_idx])^2),
  sum((y_sim_rec[,-eval_omit][,bt_idx] - saver_mat_eval[,-eval_omit][,bt_idx])^2),
  sum((y_sim_rec[,-eval_omit][,bt_idx] - softImpute_mat_eval[,-eval_omit][,bt_idx])^2),
  sum((y_sim_rec[,-eval_omit][,bt_idx] - sim_tab_zi_mbImpute[,-eval_omit][,bt_idx])^2)) / (185 * 53)
}
apply(bt_mse, 1, mean)
apply(bt_mse, 1, sd)


cor <- matrix(NA, nrow = 7, ncol = dim(mbImpute_mat[,-eval_omit])[2])
rownames(cor) <- c("Magic", "mbImpute", "scImpute", "ALRA", "saver", "softImpute", "ZI")
idx_eval <- which(!(1:193 %in% eval_omit))
for(i in 1:dim(mbImpute_mat[,-eval_omit])[2]){
  cor[1,i] <- cor(y_sim_rec[,idx_eval[i]], Magic[,idx_eval[i]])
  cor[2,i] <- cor(y_sim_rec[,idx_eval[i]], mbImpute_mat[,idx_eval[i]])
  cor[3,i] <- cor(y_sim_rec[,idx_eval[i]], scImpute_mat_eval[,idx_eval[i]])
  cor[4,i] <- cor(y_sim_rec[,idx_eval[i]], alra_mat_eval[,idx_eval[i]])
  cor[5,i] <- cor(y_sim_rec[,idx_eval[i]], saver_mat_eval[,idx_eval[i]])
  cor[6,i] <- cor(y_sim_rec[,idx_eval[i]], softImpute_mat_eval[,idx_eval[i]])
  cor[7,i] <- cor(y_sim_rec[,idx_eval[i]], sim_tab_zi_mbImpute[,idx_eval[i]])
}
which(cor[3,] > 0.95)
hist(cor[1,])
hist(cor[2,])
hist(cor[3,])
hist(cor[4,])
hist(cor[5,])
rowMeans(cor)

bt_cor <- matrix(ncol = 800, nrow = 7)
set.seed(1234)
for(i in 1:800){
  bt_idx <-  sample(1:185, 185, replace = TRUE)
  bt_cor[,i] <- rowMeans(cor[,bt_idx])
}
apply(bt_cor, 1, mean)
apply(bt_cor, 1, sd)


c(sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(saver_mat_eval[,-eval_omit]))^2),
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(scImpute_mat_eval[,-eval_omit]))^2),
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(alra_mat_eval[,-eval_omit]))^2),
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(Magic[,-eval_omit]))^2),
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(softImpute_mat_eval[,-eval_omit]))^2),
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(mbImpute_mat[,-eval_omit]))^2),
sum((colMeans(y_sim_rec[,-eval_omit]) - colMeans(sim_tab_zi_mbImpute[,-eval_omit]))^2))


# 2, 25, 32, 89, 95, 101,  103, 133

c(sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd))^2),
sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd))^2),
sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd))^2),
sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(Magic[,-eval_omit])/apply(Magic[,-eval_omit], 2, sd))^2),
sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd))^2),
sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(sim_tab_zi_mbImpute[,-eval_omit])/apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd))^2),
sum((colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd) - colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd))^2))

pdf("sim_imp/mean_divide_sd_comparison.pdf")
par(mfrow = c(3,3))
hist(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), col = "#fb9a99", main = "TRUTH", xlab = "mean/sd")
hist(colMeans(sim_tab_zi_mbImpute[,-eval_omit])/apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd), col = "#BC6C36", main = "Zero Inflated", xlab = "mean/sd")
hist(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), col = "#BC6C36", main = "TRUTH", xlab = "mean/sd")
hist(colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd), col = "#CE95CE", main = "ALRA", xlab = "mean/sd")
hist(colMeans(Magic[,-eval_omit])/apply(Magic[,-eval_omit], 2, sd), col = "#B1CEE3", main = "MAGIC", xlab = "mean/sd")
hist(colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd), col = "#b2df8a", main = "SAVER", xlab = "mean/sd")
hist(colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd), col = "#1f78b4", main = "scImpute", xlab = "mean/sd")
hist(colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd), col = "#6BAED6", main = "softImpute", xlab = "mean/sd")
hist(colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd), col = "#36A02C", main = "mbImpute", xlab = "mean/sd")
dev.off()

# ks distance
#install.packages("provenance")
library(provenance)
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(sim_tab_zi_mbImpute[,-eval_omit])/apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd))
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd))
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(Magic[,-eval_omit])/apply(Magic[,-eval_omit], 2, sd))
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd))
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd))
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd))
KS.diss(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd))

# Wasserstein distance
#install.packages("transport")
library(transport)
c(wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(sim_tab_zi_mbImpute[,-eval_omit])/apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd)),
wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd)),
wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(Magic[,-eval_omit])/apply(Magic[,-eval_omit], 2, sd)),
wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd)),
wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd)),
wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd)),
wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd)))



library(scales)
pdf("sim_imp/mean_sd_scatter.pdf")
par(mfrow = c(3,3))
plot(y = colMeans(y_sim_rec[,-eval_omit]), x = apply(y_sim_rec[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "Complete", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(sim_tab_zi_mbImpute[,-eval_omit]), x = apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "ZI", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(y_sim_rec[,-eval_omit]), x = apply(y_sim_rec[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "Complete", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(alra_mat_eval[,-eval_omit]), x = apply(alra_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "ALRA", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(Magic[,-eval_omit]), x = apply(Magic[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "MAGIC", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(saver_mat_eval[,-eval_omit]), x = apply(saver_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "SAVER", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(scImpute_mat_eval[,-eval_omit]), x = apply(scImpute_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "scImpute", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(softImpute_mat_eval[,-eval_omit]), x = apply(softImpute_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "softImpute", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(mbImpute_mat[,-eval_omit]), x = apply(mbImpute_mat[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "mbImpute", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
dev.off()

y_sim_mean_sd <- cbind(colMeans(y_sim_rec[,-eval_omit]), apply(y_sim_rec[,-eval_omit], 2, sd))
alra_mean_sd <- cbind(colMeans(alra_mat_eval[,-eval_omit]), apply(alra_mat_eval[,-eval_omit], 2, sd))
MAGIC_mean_sd <- cbind(colMeans(Magic[,-eval_omit]), apply(Magic[,-eval_omit], 2, sd))
saver_mean_sd <- cbind(colMeans(saver_mat_eval[,-eval_omit]), apply(saver_mat_eval[,-eval_omit], 2, sd))
scImpute_mean_sd <- cbind(colMeans(scImpute_mat_eval[,-eval_omit]), apply(scImpute_mat_eval[,-eval_omit], 2, sd))
softImpute_mean_sd <- cbind(colMeans(softImpute_mat_eval[,-eval_omit]), apply(softImpute_mat_eval[,-eval_omit], 2, sd))
mbImpute_mean_sd <- cbind(colMeans(mbImpute_mat[,-eval_omit]), apply(mbImpute_mat[,-eval_omit], 2, sd))
sim_tab_zi_mean_sd <- cbind(colMeans(sim_tab_zi_mbImpute[,-eval_omit]), apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd))

c(sum((y_sim_mean_sd - alra_mean_sd)^2) / (185 * 2),
sum((y_sim_mean_sd - MAGIC_mean_sd)^2) / (185 * 2),
sum((y_sim_mean_sd - saver_mean_sd)^2) / (185 * 2),
sum((y_sim_mean_sd - scImpute_mean_sd)^2) / (185 * 2),
sum((y_sim_mean_sd - softImpute_mean_sd)^2) / (185 * 2),
sum((y_sim_mean_sd - mbImpute_mean_sd)^2) / (185 * 2),
sum((y_sim_mean_sd - sim_tab_zi_mean_sd)^2) / (185 * 2))


plot(mbImpute_mean_sd[,2], mbImpute_mean_sd[,1])
plot(y_sim_mean_sd[,2], y_sim_mean_sd[,1])
plot(MAGIC_mean_sd[,2], MAGIC_mean_sd[,1])

c(wasserstein1d(cbind(colMeans(y_sim_rec[,-eval_omit]), apply(y_sim_rec[,-eval_omit], 2, sd))),
  wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd)),
  wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(Magic[,-eval_omit])/apply(Magic[,-eval_omit], 2, sd)),
  wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd)),
  wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd)),
  wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd)),
  wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd)))


df <- data.frame(y = colMeans(y_sim_rec[,-eval_omit]), x = apply(y_sim_rec[,-eval_omit], 2, sd))
ggplot(df, aes(x=x, y=y) ) +
  geom_bin2d() +
  theme_bw()

summary(colMeans(Magic[,-eval_omit]))
summary(colMeans(mbImpute_mat[,-eval_omit]))
summary(apply(Magic[,-eval_omit], 2, sd))
summary(apply(mbImpute_mat[,-eval_omit], 2, sd))
summary(apply(scImpute_mat_eval[,-eval_omit], 2, sd))



library(ggplot2)
pdf("sim_imp/z_score_comparison.pdf")
z_score <- matrix(NA, nrow = 5, ncol = dim(mbImpute_mat[,-eval_omit])[2])
rownames(z_score) <- c("Magic", "mbImpute", "scImpute", "ALRA", "saver")
for(i in 1:dim(mbImpute_mat[,-eval_omit])[2]){
  z_score_y_sim <- (y_sim_rec[,idx_eval[i]] - mean(y_sim_rec[,idx_eval[i]]))/sd(y_sim_rec[,idx_eval[i]])
  z_score_magic <- (Magic[,idx_eval[i]] - mean(Magic[,idx_eval[i]]))/sd(Magic[,idx_eval[i]])
  z_score_mbImpute <- (mbImpute_mat[,idx_eval[i]] - mean(mbImpute_mat[,idx_eval[i]]))/sd(mbImpute_mat[,idx_eval[i]])
  z_score_scImpute <- (scImpute_mat_eval[,idx_eval[i]] - mean(scImpute_mat_eval[,idx_eval[i]]))/sd(scImpute_mat_eval[,idx_eval[i]])
  z_score_alra <- (alra_mat_eval[,idx_eval[i]] - mean(alra_mat_eval[,idx_eval[i]]))/sd(alra_mat_eval[,idx_eval[i]])
  z_score_saver <- (saver_mat_eval[,idx_eval[i]] - mean(saver_mat_eval[,idx_eval[i]]))/sd(saver_mat_eval[,idx_eval[i]])
  z_score[,i] <- c( sum((z_score_y_sim - z_score_magic)^2),
                    sum((z_score_y_sim - z_score_mbImpute)^2),
                    sum((z_score_y_sim - z_score_scImpute)^2),
                    sum((z_score_y_sim - z_score_alra)^2),
                    sum((z_score_y_sim - z_score_saver)^2))
  par(mfrow = c(3,2))
  range <- max(abs(c(z_score_y_sim, z_score_magic, z_score_saver, z_score_scImpute, z_score_alra, z_score_mbImpute))) + 1
  hist(z_score_y_sim, col = "#BC6C36", xlim = c(-range, range), ylim = c(0, 50), main = paste("TRUTH_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_magic, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("Magic_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_saver, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("Saver_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_scImpute, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("scImpute_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_alra, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("alra_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_mbImpute, col = "#1D9E78", xlim = c(-range, range), ylim = c(0, 50), main = paste("mbImpute_", i, collapse = ""), breaks = seq(-range, range, by = 1))
}
dev.off()

library(ggplot2)
pdf("sim_imp/mean_sd_comparison.pdf")
z_score <- matrix(NA, nrow = 5, ncol = dim(mbImpute_mat[,-eval_omit])[2])
rownames(z_score) <- c("Magic", "mbImpute", "scImpute", "ALRA", "saver")
for(i in 1:dim(mbImpute_mat[,-eval_omit])[2]){
  z_score_y_sim <- (y_sim_rec[,idx_eval[i]] - mean(y_sim_rec[,idx_eval[i]]))/sd(y_sim_rec[,idx_eval[i]])
  z_score_magic <- (Magic[,idx_eval[i]] - mean(Magic[,idx_eval[i]]))/sd(Magic[,idx_eval[i]])
  z_score_mbImpute <- (mbImpute_mat[,idx_eval[i]] - mean(mbImpute_mat[,idx_eval[i]]))/sd(mbImpute_mat[,idx_eval[i]])
  z_score_scImpute <- (scImpute_mat_eval[,idx_eval[i]] - mean(scImpute_mat_eval[,idx_eval[i]]))/sd(scImpute_mat_eval[,idx_eval[i]])
  z_score_alra <- (alra_mat_eval[,idx_eval[i]] - mean(alra_mat_eval[,idx_eval[i]]))/sd(alra_mat_eval[,idx_eval[i]])
  z_score_saver <- (saver_mat_eval[,idx_eval[i]] - mean(saver_mat_eval[,idx_eval[i]]))/sd(saver_mat_eval[,idx_eval[i]])
  z_score[,i] <- c( sum((z_score_y_sim - z_score_magic)^2),
                    sum((z_score_y_sim - z_score_mbImpute)^2),
                    sum((z_score_y_sim - z_score_scImpute)^2),
                    sum((z_score_y_sim - z_score_alra)^2),
                    sum((z_score_y_sim - z_score_saver)^2))
  par(mfrow = c(3,2))
  range <- max(abs(c(z_score_y_sim, z_score_magic, z_score_saver, z_score_scImpute, z_score_alra, z_score_mbImpute))) + 1
  hist(z_score_y_sim, col = "#BC6C36", xlim = c(-range, range), ylim = c(0, 50), main = paste("TRUTH_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_magic, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("Magic_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_saver, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("Saver_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_scImpute, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("scImpute_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_alra, col = "#8294CD", xlim = c(-range, range), ylim = c(0, 50), main = paste("alra_", i, collapse = ""), breaks = seq(-range, range, by = 1))
  hist(z_score_mbImpute, col = "#1D9E78", xlim = c(-range, range), ylim = c(0, 50), main = paste("mbImpute_", i, collapse = ""), breaks = seq(-range, range, by = 1))
}
dev.off()

pdf("sim_imp/Simulation_results_comparison.pdf")
for(i in 1:193){
  print(i)
  par(mfrow = c(3,3))
  # df <- data.frame(cbind(y_sim[,i]))
  # colnames(df) <- "abundance"
  # ggplot(data= df, aes(x = abundance)) + geom_histogram(bins = 10, fill = "#1D9E78") +ggtitle("TRUTH") +theme(plot.title = element_text(hjust = 0.5))
  hist(y_sim[,i], col = "#BC6C36", xlim = c(0,9), ylim = c(0, 50), main = paste("TRUTH_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(log10(1.01 + sim_tab_zi[,i]), col = "#E4BA7C", xlim = c(0,9), ylim = c(0, 50), main = paste("ZI_", i, collapse = ""),breaks = seq(0, 7, by = 0.4))
  hist(saver_mat_eval[,i], col = "#8294CD", xlim = c(0,9), ylim = c(0, 50), main = paste("Saver_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(scImpute_mat_eval[,i], col = "#8294CD", xlim = c(0,9), ylim = c(0, 50), main = paste("scImpute_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(alra_mat_eval[,i], col = "#8294CD", xlim = c(0,9), ylim = c(0, 50), main = paste("ALRA_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(Magic[,i], col = "#8294CD", xlim = c(0,9), ylim = c(0, 50), main = paste("Magic_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(softImpute_mat_eval[,i], col = "#8294CD", xlim = c(0,9), ylim = c(0, 50), main = paste("softImpute_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(mbImpute_mat[,i], col = "#1D9E78", xlim = c(0,9), ylim = c(0, 50), main = paste("mbImpute_", i, collapse = ""), breaks = seq(0, 9, by = 0.4))
}
dev.off()


library(ggplot2)
pdf("sim_imp/Simulation_results_magic.pdf")
for(i in 1:193){
  print(i)
  par(mfrow = c(2,2))
  # df <- data.frame(cbind(y_sim[,i]))
  # colnames(df) <- "abundance"
  # ggplot(data= df, aes(x = abundance)) + geom_histogram(bins = 10, fill = "#1D9E78") +ggtitle("TRUTH") +theme(plot.title = element_text(hjust = 0.5))

  hist(y_sim[,i], col = "#BC6C36", xlim = c(0,9), ylim = c(0, 50), main = paste("TRUTH_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(log10(1.01 + sim_tab_zi[,i]), col = "#E4BA7C", xlim = c(0,9), ylim = c(0, 50), main = paste("ZI_", i, collapse = ""),breaks = seq(0, 7, by = 0.4))
  hist(Magic[,i], col = "#8294CD", xlim = c(0,9), ylim = c(0, 50), main = paste("Magic_", i, collapse = ""), breaks = seq(0, 7, by = 0.4))
  hist(mbImpute_mat[,i], col = "#1D9E78", xlim = c(0,9), ylim = c(0, 50), main = paste("mbImpute_", i, collapse = ""), breaks = seq(0, 9, by = 0.4))
}
dev.off()
# > sqrt(sum((y_sim - mbImpute_mat)^2))
# [1] 95.52454
# > sqrt(sum((y_sim - alra_mat_eval)^2))
# [1] 209
# > sqrt(sum((y_sim - scImpute_mat_eval)^2))
# [1] 136
# > sqrt(sum((y_sim - saver_mat_eval)^2))
# [1] 201

# Three variants of mbImpute

