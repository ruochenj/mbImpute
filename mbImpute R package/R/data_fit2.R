#' data_fit2
#'
#' @description data_fit2 is a function that aims to perform one imputation in a group of biological samples, it's a subfunction of mbImpute and intend not to be
#' used by the user.
#' @param y_sim A matrix of dimension n * m corresponding to the normalized OTU table, with n denoting the number of subjects(patients), and m denoting the number of taxa.
#' @param metadata A matrix of dimension n * p corresponding to the meta data matrix, with n denotes the number of subjects(patients), and p denotes the number of covariates (age/BMI/...).
#' @param D A matrix of dimension m * m corresponding to the taxa distance matrix.
#' @param k  A scalar corresponding to the number of nearest taxa in a phylogenetic tree we will use to impute a missing value. Theoretically, the larger k, the more accurate our imputation
#' will be (required that k <= m).
#' @param parallel A boolean indicating whether to use parallel or not. Default is FALSE.
#' @param ncores A scalar corresponding to the number of cores to use. It is only used when parallel = TRUE. Default is 1.
#' @import glmnet
#' @import doParallel
#' @import Matrix
#' @return y_imp: the imputed matrix
#'
data_fit2 <- function(y_sim, metadata, D, k, parallel = F, ncores = 1){
  #loading used functions
  gamma_norm_mix <- function(y, X){
    loglik <- function(p, alpha, beta, cov_par, var1, X, y){
      n = length(y)
      lkval <- 0
      fgam <- dgamma(y, shape = alpha, rate = beta)
      for(i in 1:n){
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i,] %*% cov_par, sd = sqrt(var1)))
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
    maxitr = 30
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
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i,] %*% eta_hat, sd = sqrt(omega_hat)))
    }
    Dev = -2 * norm_log_lik - (-2 * loglik2)
    judgement = pchisq(Dev, df = 3, lower.tail = FALSE) < 0.05

    if(!judgement || (alpha_t/beta_t > 1)){
      p_t = 0
      a_hat_t <- rep(0,n)
    }
    return(list("p" = p_t, "alpha" = alpha_t, "beta" = beta_t, "cov_par" = cov_par_t, "var" = var_t, "d" = a_hat_t, "eta" = eta_hat, "omega" = omega_hat, "Deviance" = Dev))
  }
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
      # row_gen <- rep(0, m*k + (n-1) * n + n*p)
      #
      # row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      # row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      # row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i])

      non_zero_idx <- c(((j-1)*k+1):(j*k), (m*k + (i-1)*(n-1)+1):(k*m + i*(n-1)),(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p) )
      non_zero_values <- c( as.numeric(count_mat[i,close_taxa[[j]]]), as.numeric(count_mat[-i,j]), as.numeric(covariate_mat[i]))
    }else{
      p = dim(covariate_mat)[2]
      i = row_index
      j = col_index

      close_taxa_set <- close_taxa[[j]]
      #generate a row including response and a row of design matrix.
      # row_gen <- rep(0, m*k + (n-1) * n + n*p)
      #
      # row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      # row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      # row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i,])

      non_zero_idx <- c(((j-1)*k+1):(j*k), (m*k + (i-1)*(n-1)+1):(k*m + i*(n-1)),(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p) )
      non_zero_values <- c( as.numeric(count_mat[i,close_taxa[[j]]]), as.numeric(count_mat[-i,j]), as.numeric(covariate_mat[i,]))
    }
    return(list("nz_idx" = non_zero_idx, "nz_val" = non_zero_values))
  }
  design_mat_row_gen2_imp <- function(count_mat, covariate_mat, row_index, col_index, close_taxa){
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
    # nz <- sum(y <= (log10(1.01) + 1e-6))
    # pz = 1 - nz/n
    # test = pz - 1.96 * sqrt(pz * (1-pz)/n)
    # if(nz == n || test <= 0){
    if(sum(y <= (log10(1.01) + 1e-6)) / n > 0.89){
      return(0)
    }else{
      return(1)
    }
  })

  #keep the original matrix
  y_imp <- y_sim
  # #keep the vectors that we neither impute nor borrow information
  # remain_vec <- which(unlist(filter) == 0)
  # y_rem <- y_sim[, remain_vec]

  #perform imputation on the rest
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]
  D = D[filter_vec,filter_vec]

  #apply the imputation method on the simulated matrix
  m = dim(y_sim)[2]
  n = dim(y_sim)[1]

  #x[1:50,]
  X <- cbind(rep(1,n), metadata[1:n,])
  X <- as.matrix(X)

  #identifying the group needs imputation and group doesn't need imputation

  result <- lapply(1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    return(gamma_norm_mix(y, X)$d)
  })

  hd_result <- result[1:6]
  confidence_set <- c(0,0)
  impute_set <- c(0,0)
  for(i in 1:length(result)){
    d_set <- result[[i]]
    cf_set <- which(d_set < 0.5)
    im_set <- which(d_set > 0.5)
    n_c = length(cf_set)
    confidence_set <- rbind(confidence_set, cbind(cf_set, rep(i, n_c)))
    impute_set <- rbind(impute_set, cbind(im_set, rep(i, n-n_c)))
  }

  #find k closest taxa
  # k = 10
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
  X <- X[,-1]
  p = dim(X)[2]
  if(is.null(p)){
    p = 1
  }
  #generate a row including response and a row of design matrix.
  row_length <- m * k + (n-1) * n + n*p

  # size 10000 mat
  size = 5000
  if(dim(confidence_set)[1] - 1 < size){
    remaining = (dim(confidence_set)[1] - 1) %% size
    design_mat_fit =  sparseMatrix(i = 1, j =1, x = 0, dims = c(remaining, row_length))
    track = 1:remaining
    for(i in 1:remaining){
      if(is.vector(X)){
        result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
      else{
        result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
    }
  }else{
    mat_num = ceiling( (dim(confidence_set)[1] - 1) / size)
    mat_list = list()
    if(!parallel){
      for(mat_new in 1:(mat_num-1)){
        print(mat_num-1-mat_new)
        design_mat_fit = sparseMatrix(i = 1, j =1, x = 0, dims = c(size, row_length))
        track = ((mat_new-1)*size+1):(mat_new*size)
        for(i in 1:size){
          if(is.vector(X)){
            result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
          else{
            result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
        }
        mat_list[[mat_new]] = design_mat_fit
      }
    }else{
      no_cores <- max(ncores, detectCores() - 1)
      registerDoParallel(cores=no_cores)
      cl <- makeCluster(no_cores, "FORK")
      f <- function(mat_new){
        design_mat_fit = sparseMatrix(i = 1, j =1, x = 0, dims = c(size, row_length))
        track = ((mat_new-1)*size+1):(mat_new*size)
        for(i in 1:size){
          if(is.vector(X)){
            result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
          else{
            result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
        }
        return(design_mat_fit)
      }
      mat_list <- parLapply(cl, 1:(mat_num-1), f)

    }
    remaining = (dim(confidence_set)[1] - 1) %% size
    design_mat_fit =  sparseMatrix(i = 1, j =1, x = 0, dims = c(remaining, row_length))
    track = ( (mat_num-1)*size+1):((mat_num - 1)*size+remaining)
    for(i in 1:remaining){
      if(is.vector(X)){
        result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
      else{
        result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
    }
    for(j in length(mat_list):1){
      design_mat_fit = rbind(mat_list[[j]], design_mat_fit)
    }
  }

  # y_sim <- as.matrix(y_sim@.Data)
  #generate covariate matrix for values need imputation
  # design_mat_impute <- sparseMatrix(i = 1, j =1, x = 0, dims = c(dim(impute_set)[1] - 1, row_length))
  # for(i in 2:dim(impute_set)[1]){
  #   if(i %% 10000 == 0){
  #     print(i)
  #   }
  #   if(is.vector(X)){
  #     result <- design_mat_row_gen2(y_sim, X[1:n], impute_set[i,1], impute_set[i,2], close_taxa)
  #     design_mat_impute[i-1,result$nz_idx] <- result$nz_val
  #   }
  #   else{
  #     result <- design_mat_row_gen2(y_sim, X[1:n,], impute_set[i,1], impute_set[i,2], close_taxa)
  #     design_mat_impute[i-1,result$nz_idx] <- result$nz_val
  #   }
  # }

  mse <- rep(5,4)
  psi <- c(0.1, 0.5, 1, 2)
  c1 <- list()
  for(i in 1:4){
    weights_pen <- D^psi[i]
    penalized_weights <- c()
    for(j in 1:m){
      penalized_weights <- c(penalized_weights, weights_pen[close_taxa[[j]],j])
    }
    #penalized_weights <- c(penalized_weights, rep(1, n*(n-1)), rep(c(0, rep(1, p-1)), n))
    penalized_weights <- c(penalized_weights, rep(1, n*(n-1)), rep(1, n*p))
    #print(length(penalized_weights))
    #print(dim(design_mat_fit))
    response <- y_sim[confidence_set]
    set.seed(1)
    print(dim(design_mat_fit))
    cv.result <- cv.glmnet(x = design_mat_fit, y = response, family = "gaussian", penalty.factor = penalized_weights, intercept = TRUE)
    c1[[i]] <- coef(cv.result, s = cv.result$lambda.min)
    mse[i] = min(cv.result$cvm)
    print(mse)
  }
  c1 <- c1[[which.min(mse)]]
  impute_mat = y_imp
  for(i in 2:dim(impute_set)[1]){
    if(i %% 10000 == 0){
      print(i)
    }
    if(is.vector(X)){
      result <- design_mat_row_gen2_imp(y_sim, X[1:n], impute_set[i,1], impute_set[i,2], close_taxa)
      impute_mat[impute_set[i,1], impute_set[i,2]] = min(max(y_sim[,impute_set[i,2]]), max(c(1, result) %*% c1, log10(1.01)))
    }else{
      result <- design_mat_row_gen2_imp(y_sim, X[1:n,], impute_set[i,1], impute_set[i,2], close_taxa)
      impute_mat[impute_set[i,1], impute_set[i,2]] <- min(max(y_sim[,impute_set[i,2]]), max(c(1, result) %*% c1, log10(1.01)))
    }
  }

  #print("the preserved values in zero inflated matrix is: ")
  #print(sum(y_sim == impute_mat)/(600*30))
  #print(dim(y_imp))
  y_imp = impute_mat
  #print(dim(y_imp))
  saveRDS(c1, file = "dat_sim_add_filter_coef.rds")
  saveRDS(y_imp, file = "imputed_mat_condition_as_covariate.rds")
  # print(max(impute_mat))
  # print(min(impute_mat))
  # print("the mse for psi ")
  # print(psi)
  # print("is")
  # print(mse)
  # print(sum(abs(impute_mat - y_sim) < 1e-6)/(m * n))
  # print(dim(y_imp))
  return(y_imp)
}
