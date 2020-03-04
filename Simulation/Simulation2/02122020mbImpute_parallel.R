library(glmnet)
library(parallel)

# parallel
param_set <- matrix(NA, ncol = 2, nrow = 6)
param_set[,1] = rep(c(0.1, 1, 2), 2)
param_set[,2] = c(rep("control", 3) , rep("T2D", 3))
# Calculate the number of cores
no_cores <- 6
# Initiate cluster
cl <- makeCluster(no_cores)
parLapply(cl, 1:3,
          function(x){
            library(glmnet)
            param_set <- matrix(NA, ncol = 2, nrow = 6)
            param_set[,1] = rep(c(0.1, 1, 2), 2)
            param_set[,2] = c(rep("control", 3) , rep("T2D", 3))
            #apply mbimpute
            data_fit2 <- function(y_sim, x, D, psi, k, cond){
              #loading used functions
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
              y_imp <- y_sim
              # print(dim(y_imp))
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
              #X <- cbind(rep(1,n), x[1:n,])
              X <- as.matrix(x)
              
              #identifying the group needs imputation and group doesn't need imputation
              
              result <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
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
              #print(dim(confidence_set)[1]/(dim(confidence_set)[1] + dim(impute_set)[1]))
              
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
              p = dim(x)[2]
              if(is.null(p)){
                p = 1
              }
              #generate a row including response and a row of design matrix.
              row_length <- m * k + (n-1) * n + n*p
              design_mat_fit <- matrix(0, nrow = dim(confidence_set)[1] - 1, ncol = row_length)
              for(i in 2:dim(confidence_set)[1]){
                if(is.vector(x)){
                  design_mat_fit[i-1,] <- design_mat_row_gen2(y_sim, x[1:n], confidence_set[i,1], confidence_set[i,2], close_taxa)
                }else{
                  design_mat_fit[i-1,] <- design_mat_row_gen2(y_sim, x[1:n,], confidence_set[i,1], confidence_set[i,2], close_taxa)
                }
              }
              design_mat_fit <- scale(design_mat_fit)
              print("design_mat generated")
              # print(dim(design_mat_fit))
              design_mat_fit[is.nan(design_mat_fit)] <- 0
              
              #generate covariate matrix for values need imputation
              design_mat_impute <- matrix(0, nrow = dim(impute_set)[1] - 1, ncol = row_length)
              for(i in 2:dim(impute_set)[1]){
                if(is.vector(x)){
                  design_mat_impute[i-1,] <- design_mat_row_gen2(y_sim, x[1:n], impute_set[i,1], impute_set[i,2], close_taxa)
                }
                else{
                  design_mat_impute[i-1,] <- design_mat_row_gen2(y_sim, x[1:n,], impute_set[i,1], impute_set[i,2], close_taxa)
                }
              }
              
              weights_pen <- D^psi
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
              c1 <- coef(cv.result, s = cv.result$lambda.min)
              mse = min(cv.result$cvm)
              rm(penalized_weights)
              #print(mse)
              
              # print("the imputed dimension is ")
              # print(dim(design_mat_impute))
              # print("the max of impute matrix is ")
              # print(max(design_mat_impute))
              # print("the length of c1 is ")
              # print(length(c1))
              
              #deal with design_mat_impute scaling
              design_mat_impute <- cbind(rep(1, dim(design_mat_impute)[1]), scale(design_mat_impute))
              #print(sum(is.nan(design_mat_impute))/(dim(design_mat_impute)[1] * dim(design_mat_impute)[2]))
              design_mat_impute[is.nan(design_mat_impute)] <- 0
              # print(max(design_mat_impute))
              # print(dim(design_mat_impute))
              
              #imputed_value <- cbind(rep(1, dim(design_mat_impute)[1]), scale(design_mat_impute)) %*% c1
              imputed_value <- design_mat_impute %*% c1
              ##make sure there's no NA##
              #print(head(imputed_value, 100))
              #print("nan in imputed value will be")
              # print(sum(is.nan(imputed_value)))
              impute_mat <- y_sim
              for(i in 2:dim(impute_set)[1]){
                impute_mat[impute_set[i,1], impute_set[i,2]] = max(imputed_value[i-1], log10(1.01))
              }
              #print("the preserved values in zero inflated matrix is: ")
              #print(sum(y_sim == impute_mat)/(600*30))
              #print(dim(y_imp))
              y_imp[, filter_vec] = impute_mat
              #print(dim(y_imp))
              saveRDS(c1, file = paste0(cond, "_dat", psi,"_sim_add_filter_coef.rds", collapse = ""))
              saveRDS(y_imp, file = paste0(cond, "_RUNA_imputed", psi, "_mat.rds", collapse = ""))
              saveRDS(filter_vec, file = paste0(cond, "filter_vec.rds", collapse = ""))
              # print(max(impute_mat))
              # print(min(impute_mat))
              #print("the mse for psi ")
              print(psi)
              #print("is")
              print(mse)
              #print(sum(abs(impute_mat - y_sim) < 1e-6)/(m * n))
              #print(dim(y_imp))
              # return(list(y_imp = y_imp, mse = mse))
              result <- paste("psi = ", psi, " & mse = ", mse, sep = "")
              write.csv(result, paste(cond, "_", psi,"_result.csv", sep = ""))
            }
            
            #########generation step##################
            otu_real_data_T2D <- read.csv("simulated_zi_matrix_mbImpute.csv", row.names = "X")
            
            D <- read.csv("D_sim.csv", row.names = "X")
            
            meta_data_T2D <- read.csv("simulated_meta_data_mbImpute.csv", row.names = "X")
            print(dim(otu_real_data_T2D))
            print(dim(meta_data_T2D))
            #########finding psi step###################
            # setwd("~/Dropbox/mbimpute/Simulation_reformulization/10042019_simulation2/sim_1_3_1")
            if(x <= 3){
              data_fit2(y_sim = otu_real_data_T2D, x = meta_data_T2D, D = D, psi = as.numeric(param_set[x,1]), k = 5, cond = "T2D") 
            }else{
              data_fit2(y_sim = otu_real_data_T2D, x = meta_data_T2D, D = D, psi = as.numeric(param_set[x,1]), k = 5, cond = param_set[x,2]) 
            }
          })
stopCluster(cl)

