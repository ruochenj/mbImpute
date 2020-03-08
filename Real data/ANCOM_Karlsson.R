setwd("~/Dropbox/mbimpute/code/Real data Karlsson")

otable <- read.csv("otu_real_data.csv", row.names = "X")
meta_data <- read.csv("meta_data.csv", row.names = "X")
D1 <- read.csv("D.csv", row.names = "X")

#process otu table and metagenomic data -- type 2 diabetes
library(curatedMetagenomicData)
setwd("~/Dropbox/mbimpute/code/Real data Karlsson")
# Karlsson_2015 <- Karlsson_2015.metaphlan_bugs_list.stool()
# suppressPackageStartupMessages(library(phyloseq))
# Karlsson.pseq = ExpressionSet2phyloseq( Karlsson_2015 )
#
# Karlsson.tree <- ExpressionSet2phyloseq( Karlsson.pseq, phylogenetictree = TRUE)
#
# wt = UniFrac(loman.tree, weighted=TRUE, normalized=FALSE,
#              parallel=FALSE, fast=TRUE)
# plot(hclust(wt), main="Weighted UniFrac distances")
#
# ZellerG_2014 <- ZellerG_2014.metaphlan_bugs_list.stool()
# ZellerG.pseq = ExpressionSet2phyloseq( ZellerG_2014 )
# sub_zeller <- subset_taxa(ZellerG.pseq, !is.na(Species))

Karlsson_cd <- curatedMetagenomicData("KarlssonFH_2013.metaphlan_bugs_list.stool", dryrun = FALSE, counts = TRUE, bugs.as.phyloseq = TRUE)
Karlsson_OTU <- Karlsson_cd$KarlssonFH_2013.metaphlan_bugs_list.stool@otu_table

pseq2 <- ExpressionSet2phyloseq(KarlssonFH_2013.metaphlan_bugs_list.stool(), phylogenetictree = TRUE)

rownames(pseq2@otu_table)
Karlsson_meta <- pseq2@sam_data
Karlsson_taxa_phy <- pseq2@phy_tree
Karlsson_OTU <- t(Karlsson_OTU)
Karlsson_Species_idx <- which(colnames(Karlsson_OTU) %in% rownames(pseq2@otu_table))
Karlsson_OTU_select <- Karlsson_OTU[,Karlsson_Species_idx]

#get a reasonable subset of meta features
head(Karlsson_meta)
########## CHANGE #############
apply(Karlsson_meta, 2, function(x){
  x <- as.numeric(x)
  var(x[!is.na(x)])/mean(x[!is.na(x)])
})
# which(! c("study_condition", "age", "number_reads", "triglycerides", "hba1c", "ldl", "c_ceptide", "cholesterol", "glucose", "adiponectin", "fgf_19", "hscrp", "leptin", "cd163") %in% colnames(Karlsson_meta))
# 
# which(colnames(Karlsson_meta) %in% c("study_condition", "age", "number_reads", "triglycerides", "hba1c", "ldl", "c_ceptide", "cholesterol", "glucose", "adiponectin", "fgf_19", "hscrp", "leptin", "cd163"))
chosen_meta <- Karlsson_meta[,c("study_condition", "age", "number_reads", "triglycerides", "hba1c", "ldl", "c_peptide", "cholesterol", "glucose", "adiponectin", "hscrp", "leptin")]
# chosen_meta$study_condition[chosen_meta$study_condition == "adenoma"]=1
# chosen_meta$study_condition[chosen_meta$study_condition == "control"]=0
# chosen_meta$study_condition[chosen_meta$study_condition == "CRC"]=2
chosen_meta$age <- (chosen_meta$age - min(chosen_meta$age)) / sqrt(var(chosen_meta$age))
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))
chosen_meta$triglycerides <- (chosen_meta$triglycerides - min(chosen_meta$triglycerides))/sqrt(var(chosen_meta$triglycerides))
chosen_meta$hba1c <- (chosen_meta$hba1c - min(chosen_meta$hba1c)) / sqrt(var(chosen_meta$hba1c))
chosen_meta$ldl <- (chosen_meta$ldl  - min(chosen_meta$ldl)) / sqrt(var(chosen_meta$ldl)) 
chosen_meta$c_peptide <- (chosen_meta$c_peptide - min(chosen_meta$c_peptide)) / sqrt(var(chosen_meta$c_peptide))
chosen_meta$cholesterol <- (chosen_meta$cholesterol - min(chosen_meta$cholesterol)) / sqrt(var(chosen_meta$cholesterol)) 
chosen_meta$glucose[is.na(chosen_meta$glucose)] = mean(chosen_meta$glucose[!is.na(chosen_meta$glucose)])
chosen_meta$glucose <- (chosen_meta$glucose - min(chosen_meta$glucose)) / sqrt(var(chosen_meta$glucose))
chosen_meta$adiponectin[is.na(chosen_meta$adiponectin)] <- mean(chosen_meta$adiponectin[!is.na(chosen_meta$adiponectin)])
chosen_meta$adiponectin <- (chosen_meta$adiponectin - min(chosen_meta$adiponectin)) / sqrt(var(chosen_meta$adiponectin))
chosen_meta$hscrp <- (chosen_meta$hscrp - min(chosen_meta$hscrp)) / sqrt( var(chosen_meta$hscrp))
chosen_meta$leptin[is.na(chosen_meta$leptin)] <- mean(chosen_meta$leptin[!is.na(chosen_meta$leptin)])
chosen_meta$leptin <- (chosen_meta$leptin - min(chosen_meta$leptin)) / sqrt( var(chosen_meta$leptin))

OTU_tab <- Karlsson_OTU_select

otable <- OTU_tab@.Data

study_con = meta_tab$study_condition
Vardat <- meta_tab[which(study_con == "control" | study_con == "T2D"),]
OTUdat <- otable[which(study_con == "control" | study_con == "T2D"),]
condition = study_con[which(study_con == "control" | study_con == "T2D")]

cov_mat <- Vardat[, -c(1,2)]

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

library(ggplot2)
library(exactRankTests)
library(nlme)



Vardat <- cbind(rownames(OTUdat), Vardat)
colnames(Vardat)[1] <- "Sample.ID"
OTUdat <- cbind(rownames(OTUdat), OTUdat)
colnames(OTUdat)[1] = "Sample.ID"

Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)

comparison_test=ANCOM.main(OTUdat,
                           Vardat,
                           adjusted=F,
                           repeated=F,
                           main.var="study_condition",
                           adj.formula=NULL,
                           repeat.var=NULL,
                           multcorr=2,
                           sig=0.1,
                           prev.cut=0.90,
                           longitudinal = F)

ANCOM_detected <- comparison_test$W.taxa$otu.names[comparison_test$W.taxa$detected_0.6]

write.csv(ANCOM_detected, file = "ANCOM_results.csv")

# ancom_detected <- unlist(lapply(strsplit(as.character(ANCOM_detected), split = ""), function(x){
#   paste0(x[5:length(x)], collapse = "")
# }))
# 
# # 0.05
# #ancom_detected <- c(20, 92, 40, 207, 121, 245, 116, 135, 286, 93, 290, 71, 118)
# sum(ancom_detected %in% less_identified_full) # 1
# sum(ancom_detected %in% less_identified_zi) # 1
# sum(ancom_detected %in% greater_identified_full) # 12
# sum(ancom_detected %in% greater_identified_zi)
# length(ancom_detected) # 13
# 
# sum(sim_tab_zi < log10(1.01) + 0.0001) / (dim(sim_tab)[1] * dim(sim_tab)[2])
# ancom_precison_recall <- c( (sum(ancom_detected %in% less_identified_full) +
#                                sum(ancom_detected %in% greater_identified_full)) / (length(ancom_detected)),
#                             (sum(ancom_detected %in% less_identified_full) +
#                                sum(ancom_detected %in% greater_identified_full)) / (length(less_identified_full) + length(greater_identified_full)))
# 

















