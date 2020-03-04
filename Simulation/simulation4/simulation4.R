# simulation 2

#data simulation and processing
#1. proportion of DE taxa: p
#2. Difference in DE taxa: delta_mu, delta can be negative!
#3. structure level: a. just two means: average mean, average mean + delta_mu
#                    b. just two means: average mean, average mean - delta_mu
#                    a,b are not working since we cant do t test based on that.
#                    c1. two means + random noise: average mean + noise, average mean - delta_mu + noise
#                    c2. two means + structured noise: average mean + noise, average mean - delta_mu + noise, noise = X\beta, X includes male/female
#                    c1 and c2 are simple models, d and e are complex models
#                    d. mean depends on covariates: simple case: three covariates male/female, old/young, treatment/control, 
#                       previous two with smaller coefficients: introduced parameters: beta_gen, beta2_age, beta_treat
#                    e. mean depends on covariates and taxa: simple case: three covariates male/female, old/young, treatment/control, 
#                       previous two with smaller coefficients: introduced parameters: beta_gen, beta2_age, beta_treat.
#                       Then beta closer to each other are more similar in constructed "phylogenetic tree".
#4. zero inflation introduced by zeller data: the coefficients learned from a linear model
#5. total number of samples: n
#6. total number of taxa: m
#7. 
# install.packages("tsne")

########### some functions ###################
test_generation_test_thresh <- function(control_imp_mat, crc_imp_mat, imp_method, test_thresh){
  sum(abs(control_imp_mat - sim_tab_zi_control) < 0.00001)/(300 * 60)
  sum(abs(control_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  sum(abs(crc_imp_mat - sim_tab_zi_CRC) < 0.00001)/(300 * 60)
  sum(abs(crc_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  
  control_imp_mat[1:6,1:6]
  imputed_mat <- sim_tab_zi
  imputed_mat[which(condition == 1), ] = as.matrix(crc_imp_mat)
  imputed_mat[which(condition == 0), ] = as.matrix(control_imp_mat)
  
  # sum(((sim_tab - imputed_mat)^2))
  
  #explore the tsne plot, PCA and heatmap
  pdf(file = paste0("plots/", imp_method, ".pdf", collapse = ""))
  # sim_tsne_out <- tsne(imputed_mat)
  # plot(sim_tsne_out, col = as.factor(condition))
  
  sim_rtsne_out <- Rtsne(imputed_mat)
  plot(sim_rtsne_out$Y, col = as.factor(condition))
  
  sim_pca_out <- prcomp(imputed_mat, center = TRUE)
  plot(sim_pca_out$x[,1:2], col = as.factor(condition))
  
  heatmap(as.matrix(imputed_mat), scale="row")
  
  distmat <- as.matrix(dist(as.matrix(imputed_mat)))
  diag(distmat) = max(distmat)
  heatmap(distmat)
  dev.off()
  
  #now check the wilcoxon test results
  less_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
  greater_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
  # less_identified_imp <- which(p.adjust(less_pval_imp, method = "fdr") < test_thresh)
  # greater_identified_imp <- which(p.adjust(greater_pval_imp, method = "fdr") < test_thresh)
  # 
  less_identified_imp <- which(less_pval_imp < test_thresh)
  greater_identified_imp <- which(greater_pval_imp < test_thresh)
  
  
  sum(less_identified_imp %in% less_identified_full)
  sum(greater_identified_imp %in% greater_identified_full)
  
  length(less_identified_imp)
  length(greater_identified_imp)
  
  length(less_identified_full)
  length(greater_identified_full)
  
  #make a .csv file that contains all the comparisons between t chosen taxa
  result <- matrix(NA, nrow = 3, ncol = 4)
  result[1,1] = length(less_identified_full)
  result[1,2] = length(greater_identified_full)
  result[2,1] = length(less_identified_zi)
  result[2,2] = length(greater_identified_zi)
  result[3,1] = length(less_identified_imp)
  result[3,2] = length(greater_identified_imp)
  result[2,3] = sum(less_identified_zi %in% less_identified_full)
  result[2,4] = sum(greater_identified_zi %in% greater_identified_full)
  result[3,3] = sum(less_identified_imp %in% less_identified_full)
  result[3,4] = sum(greater_identified_imp %in% greater_identified_full)
  # result[2,5] = sum(less_identified_zi %in% less_identified_imp)
  # result[3,5] = sum(greater_identified_zi %in% greater_identified_imp)
  
  colnames(result) <- c("# of taxa in CRC < control", "# of taxa in CRC > control", "overlap with full data in <", "overlap with full data in >")
  rownames(result) <- c("Full", "ZI", "Imputed")
  
  write.csv(result, file = paste0("results/", imp_method, ".csv", collapse = ""))
  return(list("precision" = (result[3,3] + result[3,4])/(result[3,1] + result[3,2]), "recall" = (result[3,3] + result[3,4])/(result[1,1] + result[1,2]) ))
}

test_generation_test_thresh_raw <- function(control_imp_mat, crc_imp_mat, imp_method, test_thresh){
  sum(abs(control_imp_mat - sim_tab_zi_control) < 0.00001)/(300 * 60)
  sum(abs(control_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  sum(abs(crc_imp_mat - sim_tab_zi_CRC) < 0.00001)/(300 * 60)
  sum(abs(crc_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  imputed_mat <- sim_tab_zi
  imputed_mat[which(condition == 1), ] = as.matrix(crc_imp_mat)
  imputed_mat[which(condition == 0), ] = as.matrix(control_imp_mat)
  
  sum(((sim_tab - imputed_mat)^2))
  
  #explore the tsne plot, PCA and heatmap
  pdf(file = paste0("plots/", imp_method, "_dm", delta_mu,"_avgmu", avg_mu,"_noise", noise, "_test_thresh", test_thresh,".pdf", collapse = ""))
  # sim_tsne_out <- tsne(imputed_mat)
  # plot(sim_tsne_out, col = as.factor(condition))
  
  sim_rtsne_out <- Rtsne(imputed_mat)
  plot(sim_rtsne_out$Y, col = as.factor(condition))
  
  sim_pca_out <- prcomp(imputed_mat, center = TRUE)
  plot(sim_pca_out$x[,1:2], col = as.factor(condition))
  
  heatmap(as.matrix(imputed_mat), scale="row")
  
  distmat <- as.matrix(dist(as.matrix(imputed_mat)))
  diag(distmat) = max(distmat)
  heatmap(distmat)
  dev.off()
  
  #now check the wilcoxon test results
  less_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "bonf")$p.value)
  greater_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "bonf")$p.value)
  less_identified_imp <- which(less_pval_imp <= test_thresh)
  greater_identified_imp <- which(greater_pval_imp <= test_thresh)
  
  sum(less_identified_imp %in% less_identified_full)
  sum(greater_identified_imp %in% greater_identified_full)
  
  length(less_identified_imp)
  length(greater_identified_imp)
  
  length(less_identified_full)
  length(greater_identified_full)
  
  #make a .csv file that contains all the comparisons between t chosen taxa
  result <- matrix(NA, nrow = 3, ncol = 4)
  result[1,1] = length(less_identified_full)
  result[1,2] = length(greater_identified_full)
  result[2,1] = length(less_identified_zi)
  result[2,2] = length(greater_identified_zi)
  result[3,1] = length(less_identified_imp)
  result[3,2] = length(greater_identified_imp)
  result[2,3] = sum(less_identified_zi %in% less_identified_full)
  result[2,4] = sum(greater_identified_zi %in% greater_identified_full)
  result[3,3] = sum(less_identified_imp %in% less_identified_full)
  result[3,4] = sum(greater_identified_imp %in% greater_identified_full)
  # result[2,5] = sum(less_identified_zi %in% less_identified_imp)
  # result[3,5] = sum(greater_identified_zi %in% greater_identified_imp)
  
  colnames(result) <- c("# of taxa in CRC < control", "# of taxa in CRC > control", "overlap with full data in <", "overlap with full data in >")
  rownames(result) <- c("Full", "ZI", "Imputed")
  
  # write.csv(result, file = paste0("results/summary_dm", delta_mu,"_avgmu", avg_mu,"_noise", noise, "_test_thresh", test_thresh, ".csv", collapse = ""))
  return(list("precision" = (result[3,3] + result[3,4])/(result[3,1] + result[3,2]), "recall" = (result[3,3] + result[3,4])/(result[1,1] + result[1,2]) ))
}

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
      return(list("d" = rep(1, n)))
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

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    lapply(colnames(data_comp), FUN = function(x){
      
    })
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }

ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

########### Data gen ################################################

setwd("~/Dropbox/mbimpute/Simulation_reformulization/10092019_simulation/10122019_simulation")
library(tsne)
library(Rtsne)
library(ggplot2)
library(exactRankTests)
library(nlme)
library(gplots)
library(pROC)

####### begin generation ###############

otable <- read.csv("otu_real_data.csv", row.names = "X")
meta_tab <- read.csv("meta_data.csv", row.names = "X")
D <- read.csv("D.csv", row.names = "X")

n = 120
m = 300
p_sample = 1/2
p_taxa = 1/3
#1 indicies with disease, 0 indicies without
condition <- rep(0, n)
delta_mu = 1
#this defines the mu for average
avg_mu = 3
noise = 2

#this is the first case
struc_level = "c1"
#indicies for disease
set.seed(1234)
samp_index <- 1:60
condition[samp_index] = 1
taxa_index <- 1:200
sim_tab <- matrix(avg_mu, nrow = n, ncol = m)
length(taxa_index)
sim_tab[samp_index, taxa_index[1:(p_taxa*m)]] = avg_mu + delta_mu
sim_tab[samp_index, taxa_index[(p_taxa*m+1):(2*p_taxa*m)]] = avg_mu - delta_mu
noise_tab <- matrix(rnorm(n*m, mean = 0, sd = noise), nrow = n, ncol = m)
sim_tab = sim_tab + noise_tab
sim_tab[sim_tab < min(otable)] = min(otable)

# heatcol <- c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177")
# heatcol <- c("#fec44f", "#fee391", "#ffffd4", "#d1e5f0", "#67a9cf", "#2166ac")

sim_pca_out_full <- prcomp(sim_tab, center = TRUE)
df1 <- data.frame(sim_pca_out_full$x[,1:2], condition)
df1$condition <- as.factor(condition)
ggplot(data = df1, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 10) +
  scale_color_manual("condition", values = c("0" = "#DA5E03", "1" = "#389E78"))


# make colors

distmat <- as.matrix(dist(as.matrix(sim_tab)))
diag(distmat) = mean(distmat)
heatmap.2(distmat, col = heatcol, dendrogram = 'none', scale="row", Colv = FALSE, Rowv = FALSE, trace = 'none')

#use t-test to identify truly DE genes in full data
less_pval_full <- apply(sim_tab, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_full <- apply(sim_tab, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
less_identified_full <- which( less_pval_full < 0.05)
greater_identified_full <- which( greater_pval_full < 0.05)

sim_tab_CRC <- sim_tab[condition == 1, ]
sim_tab_CON <- sim_tab[condition == 0, ]

###################

sim_tab_zi <- sim_tab

missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/5
  mean_up <- mean + win_len 
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}

####### learn from learning missing rate empirically ############

otable <- read.csv("otu_real_data.csv", row.names = "X")
meta_data <- read.csv("meta_data.csv", row.names = "X")
D1 <- read.csv("D.csv", row.names = "X")

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
X <- meta_data[,-2]

dim(otable)
# read in the gamma_norm_mix function
mean_record <- c()
percentage_record <- c()
D_vale_record <- c()
beta_record <- list()
for(j in 1:dim(otable)[2]){
  result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab[,-2]))
  beta_record[[j]] = result$cov_par
  # if(sum(result$d < 0.5) != 0){
  mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
  # }else{
  #   mean_record <- c(mean_record, mean(otable[which(result$d > 0.5),j]))
  # }
  # percentage of false zeros
  percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
  D_vale_record <- c(D_vale_record, result$d)
}
#filter <- which(is.nan(mean_record) | mean_record < log10(1.01) + 0.001)

hist(D_vale_record)

#total zeros
sum(otable < log10(1.01) + 0.0001)/(dim(otable)[1] * dim(otable)[2])
#0.77951
#total false zero or low counts
sum(D_vale_record == 0)/length(D_vale_record)
#0.4320713

# filter out the nan values.
plot(mean_record, percentage_record)


# filter2 <- which(is.nan(mean_record))
# mean_non_nan <- mean_record[-filter2]
# percentage_non_nan <- percentage_record[-filter2]
# 
# plot(mean_non_nan, percentage_non_nan)
# 
# plot(log(mean_non_nan), percentage_non_nan)
# abline(b = -0.81445, a = 1.18119)
# 
# summary(lm(percentage_non_nan ~sqrt(mean_non_nan)))
# 
# model <- lm(percentage_non_nan ~ log(mean_non_nan), weights=1/(percentage_non_nan+0.01))
# summary(model)
# 
# 
# 
# 
# rand_range <- range(mean_non_nan[which(percentage_non_nan > 0.6 & percentage_non_nan <= 1)])
# linear_idx <- which(percentage_non_nan <= 0.6)
# sm_linear <- summary(lm(percentage_non_nan[linear_idx] ~ mean_non_nan[linear_idx]))
# sm_linear
# sm_linear$coefficients
# # > sm_linear$coefficients
# # Estimate Std. Error   t value     Pr(>|t|)
# # (Intercept)               0.8045082 0.05887122 13.665560 3.015587e-24
# # mean_non_nan[linear_idx] -0.1786938 0.02035298 -8.779737 6.204072e-14
# sm_linear <- summary(lm(log(percentage_non_nan[linear_idx]+ 0.01) ~ mean_non_nan[linear_idx]))
# sm_linear
# sm_linear$coefficients
# # > sm_linear$coefficients
# #                           Estimate Std. Error    t value     Pr(>|t|)
# # (Intercept)               1.394903  0.2920612   4.776063 6.411534e-06
# # mean_non_nan[linear_idx] -1.022493  0.1009715 -10.126548 7.955536e-17
# 
# # the percentage between the rand part and the linear part.
# length(which(percentage_non_nan > 0.6 & percentage_non_nan <= 1)) / length(which(percentage_non_nan <= 0.6))
# #1.918367
# 

# Now use the non parametric way
rm_idx <- which(!is.nan(mean_record))
mean_record <- mean_record[rm_idx]
percentage_record <- percentage_record[rm_idx]

plot(mean_record, percentage_record)

# build a map between the mean and percentage of missing
missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/10
  mean_up <- mean + win_len 
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}
# missing_rate(1, mean_record, percentage_record)

col_mean <- colMeans(sim_tab)
zero_rate <- unlist(lapply(col_mean, FUN = function(x){
  return(missing_rate(x, mean_record, percentage_record))
}))
zero_mat <- matrix(NA, nrow = n, ncol = m)
for(i in 1:m){
  zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
}
sim_tab_zi = sim_tab * zero_mat
zero_mat <- matrix(NA, nrow = n, ncol = m)
for(i in 1:m){
  zero_mat[,i] = rbinom(n,1, runif(1, min = 0.6, max = 1))
}
sim_tab_zi = sim_tab_zi * zero_mat
sum(sim_tab_zi == 0) / (300 * 120)
sim_tab_zi[sim_tab_zi < log10(1.01)] = log10(1.01)

sim_tab_zi_CRC = sim_tab_zi[condition == 1,]
sim_tab_zi_control = sim_tab_zi[condition == 0,]

filter1 <- lapply(X = 1:ncol(sim_tab_zi_CRC), FUN = function(col_i){
  y = sim_tab_zi_CRC[,col_i]
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

filter_vec_CRC <- which(unlist(filter1) == 0)

filter <- lapply(X = 1:ncol(sim_tab_zi_control), FUN = function(col_i){
  y = sim_tab_zi_control[,col_i]
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
filter_vec_control <- which(unlist(filter1) == 0)

filter_vec <- intersect(filter_vec_control, filter_vec_CRC)

sim_tab <- sim_tab[, -filter_vec]
#use t-test to identify truly DE genes in full data
less_pval_full <- apply(sim_tab, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_full <- apply(sim_tab, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
less_identified_full <- which( less_pval_full < 0.05)
greater_identified_full <- which( greater_pval_full < 0.05)

sim_tab_CRC <- sim_tab[condition == 1, ]
sim_tab_CON <- sim_tab[condition == 0, ]


sim_tab_zi <- sim_tab_zi[, -filter_vec]

sim_tab_zi_CRC = sim_tab_zi[condition == 1,]
sim_tab_zi_control = sim_tab_zi[condition == 0,]
which(less_identified_full %in% filter_vec)
which(greater_identified_full %in% filter_vec)
length(less_identified_full)
length(greater_identified_full)
#now check the wilcoxon test results
less_pval_zi <- apply(sim_tab_zi, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_zi <- apply(sim_tab_zi, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
# less_identified_zi <- which( p.adjust(less_pval_zi, method = "fdr") < 0.1)
# greater_identified_zi <- which(  p.adjust(greater_pval_zi, method = "fdr") < 0.1)
less_identified_zi <- which( less_pval_zi  < 0.1)
greater_identified_zi <- which(  greater_pval_zi < 0.1)


ground_truth_less <- rep(1, 252)
ground_truth_more <- rep(1, 252)
ground_truth_less[less_identified_full] = 0
ground_truth_more[greater_identified_full] = 0
ground_truth <- c(ground_truth_less, ground_truth_more)

label_score <- c( less_pval_zi, greater_pval_zi)
label_score[is.na(label_score)] = 0
pROC_result_zi <- roc( ground_truth, label_score , direction = "<")
auc_zi <-  auc(pROC_result_zi)


zi_raw_TP <- sum(less_identified_zi %in% less_identified_full) + sum(greater_identified_zi %in% greater_identified_full)
zi_raw_precison_recall <- c(zi_raw_TP/(length(less_identified_zi) + length(greater_identified_zi)), 
                            zi_raw_TP/(length(less_identified_full) + length(greater_identified_full)))

# plots gen

sim_pca_out_zi <- prcomp(sim_tab_zi, center = TRUE)
df1 <- data.frame(sim_pca_out_zi$x[,1:2], condition)
df1$condition <- as.factor(condition)
ggplot(data = df1, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 10) +
  scale_color_manual("condition", values = c("0" = "#DA5E03", "1" = "#389E78"))


# make colors

distmat_zi <- as.matrix(dist(as.matrix(sim_tab_zi)))
diag(distmat_zi) = mean(distmat_zi)
heatmap.2(distmat_zi, col = heatcol, dendrogram = 'none', scale="none", Colv = FALSE, Rowv = FALSE, trace = 'none')


######### mbImpute ############

#run mbimpute on this data. 
library(glmnet)

#perform mbimpute on sim_tab_CRC and sim_tab_CON seperately
dim(sim_tab_CRC)
dim(sim_tab_CON)

D <- matrix(0, ncol = 300, nrow = 300)
close1 <- taxa_index[1: (p_taxa * m)]
close2 <- taxa_index[(p_taxa * m + 1):(2 * p_taxa * m)]
close3 <- 1:m
close3 <- which(!close3  %in% taxa_index)
for(i in 1:3){
  for(j in 1:3){
    if(i == j){
      tastr = eval(parse(text = paste0("close", i, collapse = "")))
      D[tastr, tastr] = 2
    }else{
      tastr1 = eval(parse(text = paste0("close", i, collapse = "")))
      tastr2 = eval(parse(text = paste0("close", j, collapse = "")))
      D[tastr1, tastr2] = 16
    }
  }
}

psi = c(0.1, 0.5, 1, 1.5, 2)
diag(D) = 0
D <- D[-filter_vec, -filter_vec]

sum(sim_tab_zi_CRC == log10(1.01))/(dim(sim_tab_CRC)[1] * dim(sim_tab_CRC)[2])
sum(sim_tab_zi_control == log10(1.01))/(dim(sim_tab_CON)[1] * dim(sim_tab_CON)[2])


write.csv(sim_tab_zi_CRC, "sim_tab_zi_CRC.csv")
write.csv(sim_tab_zi_control, "sim_tab_zi_CON.csv")
write.csv(D, "sim_D1.csv")

# Get the RData file

CON_imp_mat <- readRDS("control_RUNA_imputed0.1_mat.rds")
CRC_imp_mat <- readRDS("T2D_RUNA_imputed2_mat.rds")

imputed_mat <- sim_tab_zi
imputed_mat[which(condition == 1), ] = as.matrix(CRC_imp_mat)
imputed_mat[which(condition == 0), ] = as.matrix(CON_imp_mat)
sum(((sim_tab - imputed_mat)^2)) / (dim(sim_tab)[1] * dim(sim_tab)[2])

sim_pca_out_imp <- prcomp(imputed_mat, center = TRUE)
df1 <- data.frame(sim_pca_out_imp$x[,1:2], condition)
df1$condition <- as.factor(condition)
ggplot(data = df1, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 10) +
  scale_color_manual("condition", values = c("0" = "#DA5E03", "1" = "#389E78"))


# make colors

distmat_imp <- as.matrix(dist(as.matrix(imputed_mat)))
diag(distmat_imp) = mean(distmat_imp)
heatmap.2(distmat_imp, col = heatcol, dendrogram = 'none', scale="none", Colv = FALSE, Rowv = FALSE, trace = 'none')

#now check the wilcoxon test results
less_pval_mbImpute <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_mbImpute <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
# less_identified_imp <- which(p.adjust(less_pval_imp, method = "fdr") < test_thresh)
# greater_identified_imp <- which(p.adjust(greater_pval_imp, method = "fdr") < test_thresh)
# 
less_identified_mbImpute <- which(less_pval_mbImpute < 0.1)
greater_identified_mbImpute <- which(greater_pval_mbImpute < 0.1)

sum(less_identified_mbImpute %in% less_identified_full)
sum(greater_identified_mbImpute %in% greater_identified_full)

mbImpute_precison_recall <- c( (sum(less_identified_mbImpute %in% less_identified_full) +
                              sum(greater_identified_mbImpute %in% greater_identified_full)) / (length(less_identified_mbImpute) + length(greater_identified_mbImpute)),
                           (sum(less_identified_mbImpute %in% less_identified_full) +
                              sum(greater_identified_mbImpute %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))


label_score <- c( less_pval_mbImpute, greater_pval_mbImpute)
label_score[is.na(label_score)] = 0
pROC_result_mbImpute <- roc( ground_truth, label_score , direction = "<")
auc_mbImpute <-  auc(pROC_result_mbImpute)

mbImpute_precion_recall <- test_generation_test_thresh(CON_imp_mat, CRC_imp_mat, "mbImpute", 0.1)

############ softImpute ###########
library(softImpute)

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
  
  # j = 1
  # se = 1e10
  # for(i in 1:5){
  #   si_cv_1 <- softImpute(y_sim_cv_intro, rank.max = i, lambda = 0)
  #   y_imp_cv <- complete(y_sim_cv_intro, si_cv_1)
  #   y_sim_vali <- as.vector(y_imp_cv)
  #   se2 <- sum((y_sim_cv[na_intro] - y_sim_vali[na_intro])^2)
  #   print(se2)
  #   if(se2 < se){
  #     se = se2
  #     j = i
  #   }
  # }
  # print(j)
  si1 <- softImpute(y_sim_cv_intro, rank.max = 10, lambda = 0, trace.it = TRUE)
  impute_mat <- complete(y_sim_cv_intro, si1)
  
  y_imp[, filter_vec] = impute_mat
  
  return(y_imp)
}

si_imp_CRC = si_impute(sim_tab_zi_CRC)
si_imp_CON = si_impute(sim_tab_zi_control)

control_imp_mat[1:6,1:6]
imputed_mat <- sim_tab_zi
imputed_mat[which(condition == 1), ] = as.matrix(si_imp_CRC)
imputed_mat[which(condition == 0), ] = as.matrix(si_imp_CON)
sum(((sim_tab - imputed_mat)^2)) / (dim(sim_tab)[1] * dim(sim_tab)[2])

#now check the wilcoxon test results
less_pval_softImpute <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_softImpute <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
# less_identified_imp <- which(p.adjust(less_pval_imp, method = "fdr") < test_thresh)
# greater_identified_imp <- which(p.adjust(greater_pval_imp, method = "fdr") < test_thresh)
# 
less_identified_softImpute <- which(less_pval_softImpute < 0.1)
greater_identified_softImpute <- which(greater_pval_softImpute < 0.1)

sum(less_identified_softImpute %in% less_identified_full)
sum(greater_identified_softImpute %in% greater_identified_full)

softImpute_precison_recall <- c( (sum(less_identified_softImpute %in% less_identified_full) +
                                    sum(greater_identified_softImpute %in% greater_identified_full)) / (length(less_identified_softImpute) + length(greater_identified_softImpute)),
                                 (sum(less_identified_softImpute %in% less_identified_full) +
                                    sum(greater_identified_softImpute %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))

label_score <- c( less_pval_softImpute, greater_pval_softImpute)
label_score[is.na(label_score)] = 0
pROC_result_softImpute <- roc( ground_truth, label_score , direction = "<")
auc_softImpute <-  auc(pROC_result_softImpute)


softImpute_precision_recall <- test_generation_test_thresh(si_imp_CON, si_imp_CRC, "softImpute", 0.1)

############ t test #############
t_test_precision_recall <- test_generation_test_thresh(sim_tab_zi_control, sim_tab_zi_CRC, "t test", 0.1)

############ ANCOM ##############

Vardat <- cbind(rnorm(120,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")
Vardat <- cbind(1:120, Vardat)
colnames(Vardat)[1] <- "Sample.ID"
OTUdat <- 10^sim_tab_zi - 1.01
taxanames <- unlist( lapply(1:dim(sim_tab_zi)[2], FUN = function(x){
  paste("taxa", x, sep = "")
}))
colnames(OTUdat) <- taxanames
OTUdat <- cbind(1:120, OTUdat)
colnames(OTUdat)[1] = "Sample.ID"

Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)

comparison_test=ANCOM.main(OTUdat,
                           Vardat,
                           adjusted=F,
                           repeated=F,
                           main.var="condition",
                           adj.formula=NULL,
                           repeat.var=NULL,
                           multcorr=2,
                           sig=0.1,
                           prev.cut=0.90,
                           longitudinal = F)

ANCOM_detected <- comparison_test$W.taxa$otu.names[comparison_test$W.taxa$detected_0.6]
ancom_detected <- unlist(lapply(strsplit(as.character(ANCOM_detected), split = ""), function(x){
  paste0(x[5:length(x)], collapse = "")
}))

# 0.05
#ancom_detected <- c(20, 92, 40, 207, 121, 245, 116, 135, 286, 93, 290, 71, 118)
sum(ancom_detected %in% less_identified_full) # 1
sum(ancom_detected %in% less_identified_zi) # 1
sum(ancom_detected %in% greater_identified_full) # 12
sum(ancom_detected %in% greater_identified_zi)
length(ancom_detected) # 13

sum(sim_tab_zi < log10(1.01) + 0.0001) / (dim(sim_tab)[1] * dim(sim_tab)[2])
ancom_precison_recall <- c( (sum(ancom_detected %in% less_identified_full) +
                               sum(ancom_detected %in% greater_identified_full)) / (length(ancom_detected)),
                            (sum(ancom_detected %in% less_identified_full) +
                               sum(ancom_detected %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))

# ancom_pr <- matrix(nrow = 9, ncol = 2)
# sig <- seq(0, 1, 0.1)[-c(1, 11)]
# for(i in 1:9){
#   print(sig[i])
#   comparison_test=ANCOM.main(OTUdat,
#                              Vardat,
#                              adjusted=F,
#                              repeated=F,
#                              main.var="condition",
#                              adj.formula=NULL,
#                              repeat.var=NULL,
#                              multcorr=2,
#                              sig=sig[i],
#                              prev.cut=0.90,
#                              longitudinal = F)
#   ANCOM_detected <- comparison_test$W.taxa$otu.names[comparison_test$W.taxa$detected_0.6]
#   ancom_detected <- unlist(lapply(strsplit(as.character(ANCOM_detected), split = ""), function(x){
#     paste0(x[5:length(x)], collapse = "")
#   })) 
#   sum(sim_tab_zi < log10(1.01) + 0.0001) / (dim(sim_tab)[1] * dim(sim_tab)[2])
#   ancom_pr[i,] <- c( (sum(ancom_detected %in% less_identified_full) +
#                         sum(ancom_detected %in% greater_identified_full)) / (length(ancom_detected)),
#                      (sum(ancom_detected %in% less_identified_full) +
#                         sum(ancom_detected %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))
# }
# ancom_pr_ori <- ancom_pr

Vardat <- cbind(rnorm(120,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")

OTUdat <- floor(10^sim_tab_zi - 1.01)
taxanames <- unlist( lapply(1:252, FUN = function(x){
  paste("taxa", x, sep = "")
}))
colnames(OTUdat) <- taxanames

Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)
df <- data.frame(cbind(OTUdat, condition))

ANCOM_pval <- apply(OTUdat, 2, FUN = function(x){
  if(sum(x == 0)/ length(x) > 0.9){
    1
  }else{
    summary(aov(x ~ condition))[[1]][1,5] 
  }
})

ancom_detected <- which(ANCOM_pval < 0.1)

c( (sum(ancom_detected %in% less_identified_full) +
      sum(ancom_detected %in% greater_identified_full)) / (length(ancom_detected)),
   (sum(ancom_detected %in% less_identified_full) +
      sum(ancom_detected %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))

ground_truth1 <- rep(1, 252)
ground_truth1[c(less_identified_full, greater_identified_full)] = 0
label_score <- c(ANCOM_pval)
label_score[is.na(label_score)] = 0
pROC_result_ANCOM <- roc( ground_truth1, label_score , direction = "<")
auc_ANCOM <-  auc(pROC_result_ANCOM)


############ ZINB ################

Vardat <- cbind(rnorm(120,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")

OTUdat <- floor(10^sim_tab_zi - 1.01)
taxanames <- unlist( lapply(1:252, FUN = function(x){
  paste("taxa", x, sep = "")
}))
colnames(OTUdat) <- taxanames

Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)

library(pscl)

df_taxa <- c()
p_value_zinb <- c()
for(i in 1:252){
  print(i)
  # #tryCatch(
  #   sr <- summary(ml <- zeroinfl(OTUdat[,i] ~ condition, dist = "negbin"))
  #   p_val <- sr$coefficients$count[2,4]
  # )
  
  
  out <- tryCatch({
    sr <- summary(ml <- zeroinfl(OTUdat[,i] ~ condition, dist = "negbin"))
  }, error = function(cond) {
    message("error")
    return(TRUE)
  },
  warning = function(cond) {
    message("warning")
    return(FALSE)
  }, 
  finally =  {
    p_val = NA
  }
  )
  if(!(isTRUE(out))){
    if(sr$coefficients$count[2,1] > 0){
      p_val <- sr$coefficients$count[2,4]
    }else{
      p_val <- 1
    }
    
  }
  p_value_zinb <- c(p_value_zinb, p_val)
}
p_value_zinb[is.na(p_value_zinb)] <- 1
p_value_zinb_greater <- p_value_zinb
ZINB_ident_greater <- which(p_value_zinb < 0.1)
#ZINB_ident_greater <- which(p_value_zinb < 0.1)
sum(ZINB_ident_greater %in% greater_identified_full)

df_taxa <- c()
p_value_zinb <- c()
for(i in 1:252){
  print(i)
  # #tryCatch(
  #   sr <- summary(ml <- zeroinfl(OTUdat[,i] ~ condition, dist = "negbin"))
  #   p_val <- sr$coefficients$count[2,4]
  # )
  
  
  out <- tryCatch({
    sr <- summary(ml <- zeroinfl(OTUdat[,i] ~ condition, dist = "negbin"))
  }, error = function(cond) {
    message("error")
    return(TRUE)
  },
  warning = function(cond) {
    message("warning")
    return(FALSE)
  }, 
  finally =  {
    p_val = NA
  }
  )
  if(!(isTRUE(out))){
    if(sr$coefficients$count[2,1] < 0){
      p_val <- sr$coefficients$count[2,4]
    }else{
      p_val <- 1
    }
    
  }
  p_value_zinb <- c(p_value_zinb, p_val)
}
ZINB_ident_less <- which(p_value_zinb < 0.1)
p_value_zinb_less <- p_value_zinb
p_value_zinb_less[is.na(p_value_zinb_less)] <- 1
#ZINB_ident_less <- which(p_value_zinb < 0.1)
sum(ZINB_ident_less %in% less_identified_full)

label_score <- c(p_value_zinb_less, p_value_zinb_greater)
label_score[is.na(label_score)] = 0
pROC_result_zinb <- roc( ground_truth, label_score , direction = "<")
auc_zinb <-  auc(pROC_result_zinb)

ZINB_precison_recall <- c( (sum(ZINB_ident_less %in% less_identified_full) +
                           sum(ZINB_ident_greater %in% greater_identified_full)) / (length(ZINB_ident_less) + length(ZINB_ident_greater)),
                           (sum(ZINB_ident_less %in% less_identified_full) +
                              sum(ZINB_ident_greater %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))

############ Generate Test results ###########

# precison , recall, F value
precison_recall_comparison <- matrix(nrow = 6, ncol = 3)
rownames(precison_recall_comparison) <- c("ZI_raw", "mbImpute", "softImpute", "ANCOM", "ZINB", "t test")
colnames(precison_recall_comparison) <- c("Precision", "Recall", "F score")
precison_recall_comparison[1, 1:2] = zi_raw_precison_recall
precison_recall_comparison[2, 1:2] = unlist(mbImpute_precion_recall )
precison_recall_comparison[3, 1:2] = unlist(softImpute_precision_recall)
precison_recall_comparison[4, 1:2] = ancom_precison_recall
precison_recall_comparison[5, 1:2] = ZINB_precison_recall
precison_recall_comparison[6, 1:2] = unlist(t_test_precision_recall)

precison_recall_comparison[,3] =2 * precison_recall_comparison[,1] * precison_recall_comparison[,2] / (precison_recall_comparison[,1] + precison_recall_comparison[,2])
write.csv(precison_recall_comparison, "precision_recall_comparison.csv")



roc_curve_plot <- function(crc_imp_mat, control_imp_mat, imp_rate = 0.1, method){
  
  control_imp_mat[1:6,1:6]
  imputed_mat <- sim_tab_zi
  imputed_mat[which(condition == 1), ] = as.matrix(crc_imp_mat)
  imputed_mat[which(condition == 0), ] = as.matrix(control_imp_mat)
  
  # sum(((sim_tab - imputed_mat)^2))
  
  #now check the wilcoxon test results
  less_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
  greater_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
  label_score <- c( less_pval_imp, greater_pval_imp)
  label_score[is.na(label_score)] = 0
  ground_truth_less <- rep(1, 252)
  ground_truth_more <- rep(1, 252)
  ground_truth_less[less_identified_full] = 0
  ground_truth_more[greater_identified_full] = 0
  ground_truth <- c(ground_truth_less, ground_truth_more)
  pROC_result <- roc( ground_truth, label_score , direction = "<")
  gp <- ggroc(pROC_result)
  return(list(auc = auc(pROC_result), plot = gp))
}

auc_summary <- matrix(NA, ncol = 2, nrow = 1)
colnames(auc_summary) <- c("mbImpute auc", "softImpute auc")

# CON_imp_mat <- readRDS("control_imputed1_mat.rds")
# CRC_imp_mat <- readRDS("T2D_imputed2_mat.rds")
# si_imp_CON
# si_imp_CRC

# ZI_raw   mbImpute softImpute      ANCOM       ZINB 
# 0.02877698 0.49549550 0.62893082 0.17105263 0.51228070 


#########################################
library(pROC)
obj <- roc_curve_plot(crc_imp_mat = CRC_imp_mat, control_imp_mat = CON_imp_mat, imp_rate = 0.05, method = "mbImpute_0.05_")
pdf("mbImpute_0.05_roc_curve.pdf")
obj$plot
dev.off()
auc_summary[1,1] = obj$auc

obj <- roc_curve_plot(crc_imp_mat = si_imp_CRC, control_imp_mat = si_imp_CON, imp_rate = 0.05, method = "softImpute_0.05_")
pdf("softImpute_0.05_roc_curve.pdf")
obj$plot
dev.off()
auc_summary[1,2] = obj$auc

dim(obj$plot$data)
plot(obj$plot$data[,2] ~ (1-obj$plot$data[,3]))





# Bootstrap on the precision recall and F1 score
bt_precision <- matrix(nrow = 6, ncol = 800)
bt_Recall <- matrix(nrow = 6, ncol = 800)
bt_roc_mbImpute <- list()
bt_roc_softImpute <- list()
bt_roc_zinb <- list()
bt_roc_Wilcoxon <- list()
bt_auc <- matrix(nrow = 6, ncol = 800)
bt_curve_mbImpute <- matrix(nrow = 10, ncol = 800)
bt_curve_softImpute <- matrix(nrow = 10, ncol = 800)
bt_curve_ANCOM <- matrix(nrow = 10, ncol = 800)
bt_curve_Wilcoxon <- matrix(nrow = 10, ncol = 800)
bt_curve_zinb <- matrix(nrow = 10, ncol = 800)
set.seed(1234)
for(i in 1:800){
  print(i)
  sample_idx <- sample(1:dim(sim_tab_zi)[2], dim(sim_tab_zi)[2], replace = TRUE)
  
  bt_less_full <- sample_idx[which(sample_idx %in% less_identified_full)]
  bt_greater_full <- sample_idx[which(sample_idx %in% greater_identified_full)]
  
  ground_truth_less <- rep(1, 252)
  ground_truth_more <- rep(1, 252)
  ground_truth_less[which(sample_idx %in% less_identified_full)] = 0
  ground_truth_more[which(sample_idx %in% greater_identified_full)] = 0
  ground_truth <- c(ground_truth_less, ground_truth_more)
  
  
  # mbImpute
  bt_less_mbImpute <- sample_idx[which(sample_idx %in% less_identified_mbImpute)]
  bt_greater_mbImpute <- sample_idx[which(sample_idx %in% greater_identified_mbImpute)]
  bt_precision[1,i] <- (sum(bt_less_mbImpute %in% bt_less_full) +
        sum(bt_greater_mbImpute %in% bt_greater_full)) / (length(bt_less_mbImpute) + length(bt_greater_mbImpute))
  bt_Recall[1,i] <- (sum(bt_less_mbImpute %in% bt_less_full) +
        sum(bt_greater_mbImpute %in% bt_greater_full)) / (length(bt_less_full) + length(bt_greater_full))
  
  label_score <- c( less_pval_mbImpute[sample_idx], greater_pval_mbImpute[sample_idx])
  label_score[is.na(label_score)] = 0
  pROC_result <- roc( ground_truth, label_score , direction = "<")
  gp <- ggroc(pROC_result)
  idx_qt <-  quantile(1:length(pROC_result$sensitivities), c(0, 0.25, 0.5, 0.75, 1))
  bt_curve_mbImpute[,i] <- c(pROC_result$sensitivities[idx_qt], pROC_result$specificities[idx_qt]) 
  bt_auc[1,i] <-  auc(pROC_result)
  
  # softImpute
  bt_less_softImpute <- sample_idx[which(sample_idx %in% less_identified_softImpute)]
  bt_greater_softImpute <- sample_idx[which(sample_idx %in% greater_identified_softImpute)]
  bt_precision[2,i] <- (sum(bt_less_softImpute %in% bt_less_full) +
                            sum(bt_greater_softImpute %in% bt_greater_full)) / (length(bt_less_softImpute) + length(bt_greater_softImpute))
  bt_Recall[2,i] <- (sum(bt_less_softImpute %in% bt_less_full) +
                            sum(bt_greater_softImpute %in% bt_greater_full)) / (length(bt_less_full) + length(bt_greater_full))
  label_score <- c( less_pval_softImpute[sample_idx], greater_pval_softImpute[sample_idx])
  label_score[is.na(label_score)] = 0
  pROC_result <- roc( ground_truth, label_score , direction = "<")
  idx_qt <-  quantile(1:length(pROC_result$sensitivities), c(0, 0.25, 0.5, 0.75, 1))
  bt_curve_softImpute[,i] <- c(pROC_result$sensitivities[idx_qt], pROC_result$specificities[idx_qt]) 
  bt_auc[2,i] <-  auc(pROC_result)
  
  # ANCOM
  bt_ANCOM <- sample_idx[which(sample_idx %in% ancom_detected)]
  bt_precision[3,i]  <- (sum(ancom_detected %in% bt_less_full) +
                                 sum(ancom_detected %in% bt_greater_full)) / (length(bt_ANCOM))
  bt_Recall[3,i] <- (sum(ancom_detected %in% bt_less_full) +
                                 sum(ancom_detected %in% bt_greater_full)) / (length(bt_less_full) + length(bt_greater_full))
  label_score <- c(ANCOM_pval[sample_idx], ANCOM_pval[sample_idx])
  label_score[is.na(label_score)] = 0
  pROC_result <- roc( ground_truth, label_score , direction = "<")
  idx_qt <-  quantile(1:length(pROC_result$sensitivities), c(0, 0.25, 0.5, 0.75, 1))
  bt_curve_ANCOM[,i] <- c(pROC_result$sensitivities[idx_qt], pROC_result$specificities[idx_qt]) 
  bt_auc[3,i] <-  auc(pROC_result)
  
  
  # ZINB
  bt_less_ZINB <- sample_idx[which(sample_idx %in% ZINB_ident_less)]
  bt_greater_ZINB <- sample_idx[which(sample_idx %in% ZINB_ident_greater)]
  bt_precision[4,i] <- (sum(bt_less_ZINB %in% bt_less_full) +
                            sum(bt_greater_ZINB %in% bt_greater_full)) / (length(bt_less_ZINB) + length(bt_greater_ZINB))
  bt_Recall[4,i] <- (sum(bt_less_ZINB %in% bt_less_full) +
                            sum(bt_greater_ZINB %in% bt_greater_full)) / (length(bt_less_full) + length(bt_greater_full))
  label_score <- c( p_value_zinb_less[sample_idx], p_value_zinb_greater[sample_idx])
  label_score[is.na(label_score)] = 0
  pROC_result <- roc( ground_truth, label_score , direction = "<")
  idx_qt <-  quantile(1:length(pROC_result$sensitivities), c(0, 0.25, 0.5, 0.75, 1))
  bt_curve_zinb[,i] <- c(pROC_result$sensitivities[idx_qt], pROC_result$specificities[idx_qt]) 
  bt_auc[4,i] <-  auc(pROC_result)
  
  # Wilcoxon
  bt_less_Wilcoxon <- sample_idx[which(sample_idx %in% less_identified_zi)]
  bt_greater_Wilcoxon <- sample_idx[which(sample_idx %in% greater_identified_zi)]
  bt_precision[5,i] <- (sum(bt_less_Wilcoxon %in% bt_less_full) +
                          sum(bt_greater_Wilcoxon %in% bt_greater_full)) / (length(bt_less_Wilcoxon) + length(bt_greater_Wilcoxon))
  bt_Recall[5,i] <- (sum(bt_less_Wilcoxon %in% bt_less_full) +
                       sum(bt_greater_Wilcoxon %in% bt_greater_full)) / (length(bt_less_full) + length(bt_greater_full))
  label_score <- c( less_pval_zi[sample_idx], greater_pval_zi[sample_idx])
  label_score[is.na(label_score)] = 0
  pROC_result <- roc( ground_truth, label_score , direction = "<")
  idx_qt <-  quantile(1:length(pROC_result$sensitivities), c(0, 0.25, 0.5, 0.75, 1))
  bt_curve_Wilcoxon[,i] <- c(pROC_result$sensitivities[idx_qt], pROC_result$specificities[idx_qt]) 
  bt_auc[5,i] <-  auc(pROC_result)
}

bt_precision[is.nan(bt_precision)] = 0
apply(bt_precision, 1, sd)[-6] / precison_recall_comparison[c(2,3,4,5,1),1]
c(0.04324161, 0.03924473, 0.16956254, 0.03867931, 0.06655002)
apply(bt_Recall, 1, sd)[-6] / precison_recall_comparison[c(2,3,4,5,1),2]
c(0.041887898, 0.039921226, 0.008838368, 0.038971670, 0.032733448 )
bt_F1_score <- 2 * bt_precision * bt_Recall / (bt_precision + bt_Recall)
bt_F1_score[is.nan(bt_F1_score)] <- 0
rowMeans(bt_F1_score)
apply(bt_F1_score, 1, sd)[-6] / precison_recall_comparison[c(2,3,4,5,1),3]
c(0.03713716, 0.03453192, 0.01622628, 0.03377985, 0.04114185)
apply(bt_auc, 1, sd) / apply(bt_auc, 1, mean)
c(0.02414132, 0.02620456, NA, 0.02520467, 0.02992104)
apply(bt_auc, 1, mean)
c(0.8204205, 0.7722390, NA, 0.7366401, 0.6769746)

setwd("~/Dropbox/mbimpute/Simulation_reformulization/02012020")
saveRDS(bt_precision, "bt_precision.RData")
saveRDS(bt_Recall, "bt_recall.RData")
saveRDS(bt_F1_score, "bt_F1_score.RData")
saveRDS(bt_auc, "bt_auc.RData")
saveRDS(bt_roc_mbImpute, "bt_ROC_mbImpute.RData")
saveRDS(bt_roc_softImpute, "bt_ROC_softImpute.RData")
saveRDS(bt_roc_zinb, "bt_ROC_zinb.RData")
saveRDS(bt_roc_Wilcoxon, "bt_ROC_Wilcoxon.RData")


idx1 <- quantile(1:493, c(0, 0.25, 0.5, 0.75, 1))
idx2 <- quantile(1:505, c(0, 0.25, 0.5, 0.75, 1))
idx3 <- quantile(1:170, c(0, 0.25, 0.5, 0.75, 1))
idx4 <- quantile(1:225, c(0,0.25, 0.5, 0.75, 1))

df <- data.frame(cbind(c(pROC_result_mbImpute$sensitivities[idx1], pROC_result_softImpute$sensitivities[idx2], 
                       pROC_result_zinb$sensitivities[idx3], pROC_result_zi$sensitivities[idx2]), 
                       c(pROC_result_mbImpute$specificities[idx1], pROC_result_softImpute$specificities[idx2], 
                         pROC_result_zinb$specificities[idx3], pROC_result_zi$specificities[idx2]),
                       c(rep("MBImpute", length(idx1)), rep("softImpute", length(idx2)),
                         rep("ZINB", length(idx3)), rep("Wilcoxon", length(idx2)))))
df$X1 <- as.numeric(as.character(df$X1))
df$X2 <- as.numeric(as.character(df$X2))

colnames(df) <- c("sensitivity", "specificity", "Impute_method")

ggplot(df, aes(x = specificity, y = sensitivity, col = Impute_method)) + geom_line(size= 2) + scale_x_reverse() +
  scale_color_manual("Impute_method", values = c("softImpute" = "#6BAED6", "MBImpute" = "#33a02c", "ZINB" = "#FEE9AC", "Wilcoxon" = "#FDBB8A")) + 
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sd_mbImpute <- apply(bt_curve_mbImpute, 1, sd)
sd_softImpute <- apply(bt_curve_softImpute, 1, sd)
sd_zinb <- apply(bt_curve_zinb, 1, sd)
sd_Wilcoxon <- apply(bt_curve_Wilcoxon, 1, sd)
sd_ANCOM <- apply(bt_curve_ANCOM, 1, sd)

df <- data.frame(cbind(c(pROC_result_mbImpute$sensitivities[idx1], pROC_result_softImpute$sensitivities[idx2], 
                         pROC_result_zinb$sensitivities[idx3], pROC_result_zi$sensitivities[idx2], pROC_result_ANCOM$sensitivities[idx4]), 
                       c(pROC_result_mbImpute$specificities[idx1], pROC_result_softImpute$specificities[idx2], 
                         pROC_result_zinb$specificities[idx3], pROC_result_zi$specificities[idx2], pROC_result_ANCOM$specificities[idx4]),
                       c(rep("MBImpute", length(idx1)), rep("softImpute", length(idx2)),
                         rep("ZINB", length(idx3)), rep("Wilcoxon", length(idx2)), rep("ANCOM", length(idx4)))))
df$X1 <- as.numeric(as.character(df$X1))
df$X2 <- as.numeric(as.character(df$X2))

df1 <- data.frame("min" = pROC_result_mbImpute$sensitivities[idx1] - 2*max(sd_mbImpute[1:5]), "max" = pROC_result_mbImpute$sensitivities[idx1] + 2*max(sd_mbImpute[1:5]), "specificity" = pROC_result_mbImpute$specificities[idx1])
df2 <- data.frame("min" = pROC_result_softImpute$sensitivities[idx2] - 2*max(sd_softImpute[1:5]), "max" = pROC_result_softImpute$sensitivities[idx2] + 2*max(sd_softImpute[1:5]), "specificity" = pROC_result_softImpute$specificities[idx2])
df3 <- data.frame("min" = pROC_result_zinb$sensitivities[idx3] - 2*max(sd_zinb[1:5]), "max" = pROC_result_zinb$sensitivities[idx3] + 2*max(sd_zinb[1:5]), "specificity" = pROC_result_zinb$specificities[idx3])
df4 <- data.frame("min" = pROC_result_zi$sensitivities[idx2] - 2*max(sd_Wilcoxon[1:5]), "max" = pROC_result_zi$sensitivities[idx2] + 2*max(sd_Wilcoxon[1:5]), "specificity" = pROC_result_zi$specificities[idx2])
df5 <- data.frame("min" = pROC_result_ANCOM$sensitivities[idx4] - 2*max(sd_ANCOM[1:5]), "max" = pROC_result_ANCOM$sensitivities[idx4] + 2*max(sd_ANCOM[1:5]), "specificity" = pROC_result_ANCOM$specificities[idx4])


colnames(df) <- c("sensitivity", "specificity", "Impute_method")

setwd("~/Dropbox/mbimpute/Simulation_reformulization/02012020")
pdf("AUC_comparison.pdf", width = 10)
ggplot() +
  geom_ribbon(data = df1, aes(ymin=min,ymax=max,x=specificity), fill = alpha("#33a02c", 0.2)) +
  geom_ribbon(data = df2, aes(ymin=min,ymax=max,x=specificity), fill = alpha("#6BAED6", 0.2)) +
  geom_ribbon(data = df3, aes(ymin=min,ymax=max,x=specificity), fill = alpha("#FEE9AC", 0.4)) +
  geom_ribbon(data = df4, aes(ymin=min,ymax=max,x=specificity), fill = alpha("#FDBB8A", 0.2)) +
  geom_ribbon(data = df5, aes(ymin=min,ymax=max,x=specificity), fill = alpha("#FED193", 0.2)) +
  geom_line(data = df, aes(x = specificity, y = sensitivity, col = Impute_method), size = 1) + scale_x_reverse() + 
  scale_color_manual("Impute_method", values = c("softImpute" = "#6BAED6", "MBImpute" = "#33a02c", "ZINB" = "#FEE9AC", "Wilcoxon" = "#FDBB8A", "ANCOM" = "#FED193")) + 
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

auc_mbImpute
auc_softImpute
auc_zinb
auc_zi

apply(bt_auc, 1, sd)


df2 <- data.frame(cbind(c(auc_mbImpute, auc_softImpute, auc_zinb, auc_zi, auc_ANCOM), c("MBImpute", "softImpute", "ZINB", "Wilcoxon", "ANCOM")))
colnames(df2) <- c("auc", "method")
df2[,1] <- as.numeric(as.character(df2[,1]))
df2[,2] <- factor(df2[,2], levels = c("Wilcoxon", "ANCOM", "ZINB", "MBImpute", "softImpute")) 
pdf("auc summary.pdf")
ggplot(df2, aes(x = method, y = auc, fill = method)) + geom_bar(stat = "identity") + coord_cartesian(ylim=c(0,1))+
  scale_fill_manual("method", values = c("softImpute" = "#6BAED6", "MBImpute" = "#33a02c", "ZINB" = "#FEE9AC", "Wilcoxon" = "#FDBB8A", "ANCOM" = "#FED193")) + 
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

apply(bt_auc, 1, sd)[c(1,2,3,4,5)] / c(auc_mbImpute, auc_softImpute, auc_ANCOM, auc_zinb, auc_zi)



