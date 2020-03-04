# simulation 3
setwd("~/Dropbox/mbimpute/Simulation_reformulization/06172019")
library(geoCount)
library(MASS)

otable <- read.csv("otu_real_data.csv", row.names = "X")
meta_tab <- read.csv("meta_data.csv")
D <- read.csv("D.csv", row.names = "X")

#randomly choose 200 taxa from the 500
set.seed(1234)
taxa_idx <- sample(1:449, size = 200, replace = FALSE)
sample_idx <- sample(1:154, size = 50, replace = FALSE)

covariate_mat <- meta_tab[sample_idx,]
D <- D[taxa_idx, taxa_idx]
D_ori = D
# Now generate coefficients for covariates.
numeric_cov <- data.matrix(covariate_mat)

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
filter <- lapply(X = 1:ncol(otable), FUN = function(col_i){
  y = otable[,col_i]
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

filter_vec <- which(unlist(filter) == 1)
otable = otable[, filter_vec]

mean_record <- c()
percentage_record <- c()
D_vale_record <- c()
beta_record <- list()
for(j in 1:dim(otable)[2]){
  result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab))
  beta_record[[j]] = result$cov_par
  mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
  # percentage of false zeros
  percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
  D_vale_record <- c(D_vale_record, result$d)
}

set.seed(1234)
m = 200
n = 50

sample_grp <- 5
sample_sd <- 1
taxa_sd <- 1
samp_param <- 6
taxa_grp = 10
taxa_param = 12

error = 1

set.seed(1234)
#separate D into groups. 
taxa_sim_idx <- ceiling(runif(m, 0, taxa_grp))
D_taxa <- matrix(0, nrow = m, ncol = m)
for(i in 2:m){
  for(j in 1:(i-1)){
    D_taxa[i,j] = abs(taxa_sim_idx[i] - taxa_sim_idx[j]) * 2 + 2
  }
}
D_taxa = D_taxa + t(D_taxa)

#use scale_D <- D/8 since otherwise the value for sigma_mat will be too small.
# scale_D <- D/21
# sigma_mat <- exp(-scale_D^2)
sigma_mat1 <- U2Z(D_taxa, cov.par = c(1, taxa_param, 1.3))
# for one sample we will have:
mvrnorm(n = 1, mu = rep(0, m), Sigma = sigma_mat1)
# Now generate coefficients for covariates.
numeric_cov <- data.matrix(covariate_mat)

samp_sim_idx <- ceiling(runif(n, 0, sample_grp))
D_samp <- matrix(0, nrow = n, ncol = n)
for(i in 2:n){
  for(j in 1:(i-1)){
    D_samp[i,j] = abs(samp_sim_idx[i] - samp_sim_idx[j]) * 2 + 2
  }
}
D_samp = D_samp + t(D_samp)
sigma_mat2 <- U2Z(D_samp, cov.par = c(1, samp_param, 1.3))
# for one sample we will have:

y_full <- matrix(NA, nrow = n, ncol = m)
j = 1
i = 1
while(j <= 200){
  beta_j <- beta_record[[i]]
  if(!is.null(beta_j)){
    y_full[,j] = numeric_cov %*% beta_j
    i = i+1
    j = j+1
  }else{
    i = i+1
  }
}

y_full_covariate = y_full
y_full_taxa = y_full
y_full_sample = y_full

# Next, we add taxa dependent value, sample dependent value and random error
for(i in 1:n){
  y_full_taxa[i,] <- taxa_sd * mvrnorm(n = 1, mu = rep(3, m), Sigma = sigma_mat1) + rnorm(m, sd = error)
}

for(j in 1:m){
  y_full_sample[,j] <- sample_sd * mvrnorm(n = 1, mu = rep(3, n), Sigma = sigma_mat2)+ rnorm(n, sd = error)
}

y_full_covariate <- y_full + matrix(rnorm(m*n, sd = error), n,m)

# Next, we add taxa dependent value, sample dependent value and random error
for(i in 1:n){
  y_full[i,] <- y_full[i,] + taxa_sd * mvrnorm(n = 1, mu = rep(0, m), Sigma = sigma_mat1)
}
for(j in 1:m){
  y_full[,j] <- y_full[,j] + sample_sd * mvrnorm(n = 1, mu = rep(0, n), Sigma = sigma_mat2)
}
#y_full <- y_full + matrix(rnorm(n*m), nrow = n, ncol = m)
y_full <- y_full + matrix(rnorm(m*n, sd = error), n,m)


y_full_sample[y_full_sample < log10(1.01) + 0.001] = log10(1.01)
y_full_covariate[y_full_covariate < log10(1.01) + 0.001] = log10(1.01)
y_full_taxa[y_full_taxa < log10(1.01) + 0.001] = log10(1.01)
y_full[y_full < log10(1.01) + 0.001] = log10(1.01)


# build a map between the mean and percentage of missing
missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/3
  mean_up <- mean + win_len 
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}
# missing_rate(1, mean_record, percentage_record)

# introduce missing values
zi_intro <- function(y_sim){
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
  return(sim_tab_zi)
}

setwd("~/Dropbox/mbimpute/02162020_simulation3")
write.csv(D_ori, "D.csv")
write.csv(numeric_cov, "meta_data.csv")
write.csv(y_full_covariate, "covariate_full.csv")
write.csv(y_full_taxa, "taxa_full.csv")
write.csv(y_full_sample, "sample_full.csv")
write.csv(y_full, "comb_full.csv")

sim_tab_zi_covariate <- zi_intro(y_full_covariate)
sim_tab_zi_taxa <- zi_intro(y_full_taxa)
sim_tab_zi_sample <- zi_intro(y_full_sample)
sim_tab_zi_full <- zi_intro(y_full)

sum((sim_tab_zi_full - y_full)^2)/ (50 * 200)
sum((sim_tab_zi_taxa - y_full_taxa)^2)/ (50 * 200)
sum((sim_tab_zi_sample - y_full_sample)^2)/ (50 * 200)
sum((sim_tab_zi_covariate - y_full_covariate)^2)/ (50 * 200)

write.csv(sim_tab_zi_covariate, "covariate_zi.csv")
write.csv(sim_tab_zi_taxa, "taxa_zi.csv")
write.csv(sim_tab_zi_sample, "sample_zi.csv")
write.csv(sim_tab_zi_full, "comb_zi.csv")

# load the mbImpute data

# comb 1
imputed_comb <- readRDS("RData_mbImpute/comb_RUNA_imputed1_mat.rds")

# covariate 2
imputed_covariate <- readRDS("RData_mbImpute/Covariate_RUNA_imputed2_mat.rds")

# sample 2
imputed_sample <- readRDS("RData_mbImpute/sample_RUNA_imputed2_mat.rds")

# taxa 
imputed_taxa <- readRDS("RData_mbImpute/Taxa_RUNA_imputed2_mat.rds")

sum((imputed_comb - y_full)^2)/ (50 * 200)
sum((imputed_taxa - y_full_taxa)^2)/ (50 * 200)
sum((imputed_sample - y_full_sample)^2)/ (50 * 200)
sum((imputed_covariate - y_full_covariate)^2)/ (50 * 200)


bt_zi_MSE <- matrix(nrow = 4, ncol = 800)
bt_imp_MSE <- matrix(nrow = 4, ncol = 800)
for(i in 1:800){
  bt_idx <- sample(1:dim(sim_tab_zi_full)[2], dim(sim_tab_zi_full)[2], replace = TRUE)
  bt_zi_MSE[,i] = c(sum((sim_tab_zi_full[,bt_idx] - y_full[,bt_idx])^2)/ (50 * 200),
                    sum((sim_tab_zi_taxa[,bt_idx]  - y_full_taxa[,bt_idx] )^2)/ (50 * 200),
                    sum((sim_tab_zi_sample[,bt_idx]  - y_full_sample[,bt_idx] )^2)/ (50 * 200),
                    sum((sim_tab_zi_covariate[,bt_idx]  - y_full_covariate[,bt_idx] )^2)/ (50 * 200))
  bt_imp_MSE[,i] <- c(sum((imputed_comb[,bt_idx] - y_full[,bt_idx])^2)/ (50 * 200),
  sum((imputed_taxa[,bt_idx] - y_full_taxa[,bt_idx])^2)/ (50 * 200),
  sum((imputed_sample[,bt_idx] - y_full_sample[,bt_idx])^2)/ (50 * 200),
  sum((imputed_covariate[,bt_idx] - y_full_covariate[,bt_idx])^2)/ (50 * 200))
}
apply(bt_zi_MSE, 1, sd) / c(sum((sim_tab_zi_full - y_full)^2)/ (50 * 200),
sum((sim_tab_zi_taxa - y_full_taxa)^2)/ (50 * 200),
sum((sim_tab_zi_sample - y_full_sample)^2)/ (50 * 200),
sum((sim_tab_zi_covariate - y_full_covariate)^2)/ (50 * 200))
# 0.05242838 0.03383291 0.04457460 0.04958907

apply(bt_imp_MSE, 1, sd) / c(sum((imputed_comb - y_full)^2)/ (50 * 200),
                             sum((imputed_taxa - y_full_taxa)^2)/ (50 * 200),
                             sum((imputed_sample - y_full_sample)^2)/ (50 * 200),
                             sum((imputed_covariate - y_full_covariate)^2)/ (50 * 200))
#  0.06506816 0.09390213 0.09319147 0.09198824

df2 <- data.frame(cbind(c(sum((imputed_comb - y_full)^2)/ (50 * 200),
                          sum((imputed_taxa - y_full_taxa)^2)/ (50 * 200),
                          sum((imputed_sample - y_full_sample)^2)/ (50 * 200),
                          sum((imputed_covariate - y_full_covariate)^2)/ (50 * 200),
                          sum((sim_tab_zi_full - y_full)^2)/ (50 * 200),
                          sum((sim_tab_zi_taxa - y_full_taxa)^2)/ (50 * 200),
                          sum((sim_tab_zi_sample - y_full_sample)^2)/ (50 * 200),
                          sum((sim_tab_zi_covariate - y_full_covariate)^2)/ (50 * 200))),
                        c(rep("MBImpute", 4), rep("zi", 4)),
                  rep(c("comb", "taxa", "sample", "covariate"),2))
colnames(df2) <- c("mse", "method", "condition")
df2[,1] <- as.numeric(as.character(df2[,1]))
df2[,2] <- factor(df2[,2], levels = c("zi", "MBImpute")) 
df2[,3] <- factor(df2[,3], levels = c("covariate", "sample", "taxa", "comb")) 
pdf("hist.pdf")
ggplot(df2, aes(x = condition, y = mse, fill = method)) + geom_bar(stat = "identity", position=position_dodge())  + coord_cartesian(ylim=c(0,8))+
  scale_fill_manual("method", values = c("MBImpute" = "#1D9E78", "zi" = "#B2B2B2" )) + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()










