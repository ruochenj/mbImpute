# run softImpute on real data
library(softImpute)
si_impute <- function(sim_tab){
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
    if(se2 < se){
      se = se2
      j = i
    }
  }
  
  si1 <- softImpute(y_sim_cv_intro, rank.max = j, lambda = 0, trace.it = TRUE)
  impute_mat <- complete(y_sim_cv_intro, si1)
  
  y_imp[, filter_vec] = impute_mat
  
  return(y_imp)
}

# si_imp_crc = si_impute(sim_tab_zi_CRC)
# si_imp_control = si_impute(sim_tab_zi_control)

otu_real_data_T2D <- read.csv("otu_real_data_T2D.csv", row.names = "X")
otu_real_data_IGT <- read.csv("otu_real_data_IGT.csv", row.names = "X")
otu_real_data_control <- read.csv("otu_real_data_control.csv", row.names = "X")

T2D_mat_si <- si_impute(otu_real_data_T2D)
control_mat_si <- si_impute(otu_real_data_control)
IGT_mat_si <- si_impute(otu_real_data_IGT)


#regroup the conditions
mat_imp_si <- matrix(1, ncol = m, nrow = n)
mat_imp_si[T2D_idx,] = as.matrix(T2D_mat_si)
mat_imp_si[control_idx,] = as.matrix(control_mat_si)
mat_imp_si[IGT_idx,] = as.matrix(IGT_mat_si)

########imputation2#############
imp_mat2 <- mat_imp

set.seed(42) # Set a seed if you want reproducible results
library(Rtsne)
imp_rtsne_out <- Rtsne(imp_mat2) # Run TSNE

imp_mat2[imp_mat2 > 6] = 5

png("plots/si_rtsne_imp.png")
plot(imp_rtsne_out$Y, col = as.factor(study_con))
dev.off()
library(tsne)
imp_tsne_out <- tsne(imp_mat2)
png("plots/si_tsne_imp.png")
plot(imp_tsne_out, col = as.factor(study_con))
dev.off()
imp_pca_out <- prcomp(imp_mat2, center = TRUE)
par(mfrow = c(1,1))
png("plots/si_pca_1_2_imp.png")
plot(imp_pca_out$x[,1:2], col = as.factor(study_con))
dev.off()
png("plots/si_pca_1_3_imp2.png")
plot(imp_pca_out$x[,c(1,3)], col = as.factor(study_con))
dev.off()
png("plots/si_pca_2_3_imp.png")
plot(imp_pca_out$x[,2:3], col = as.factor(study_con))
dev.off()

png("plots/si_htmap_imp.png")
heatmap(as.matrix(imp_mat2), Colv = NA, Rowv = NA, scale="row")
dev.off()

png("plots/si_htmap_imp_ordered.png")
heatmap(as.matrix(imp_mat2), scale = "row")
dev.off()

png("plots/si_htmap_imp_dis.png")
distmat <- as.matrix(dist(as.matrix(imp_mat2)))
diag(distmat) = max(distmat)
heatmap(distmat)
dev.off()

selected_taxa <- c()
p_T2D_control <- rep(1, m)
p_T2D_IGT <- rep(1, m)
p_IGT_control <- rep(1, m)
for(j in 1:m){
  #test <- pairwise.wilcox.test(mat_imp[,j], g = study_con, p.adjust.method = "none", alternative = "greater")
  test <- pairwise.t.test(mat_imp_si[,j], g = study_con, p.adjust.method = "none", alternative = "greater")
  p_T2D_control[j] = test$p.value[2,1] #T2D > control
  p_T2D_IGT[j] <- test$p.value[2,2] #T2D > IGT
  p_IGT_control[j] <- test$p.value[1,1] #IGT > control
}

#install.packages("fdrtool")
grt_T2D_control <- which(p.adjust(p_T2D_control, method = "fdr") < 0.1)
grt_T2D_IGT <- which(p.adjust(p_T2D_IGT, method = "fdr") < 0.1)
grt_IGT_control <- which(p.adjust(p_IGT_control, method = "fdr") < 0.1)

grt_tc_df <- cbind(grt_T2D_control, p.adjust(p_T2D_control, method = "bonferroni")[grt_T2D_control])
grt_ti_df <- cbind(grt_T2D_IGT, p.adjust(p_T2D_IGT, method = "bonferroni")[grt_T2D_IGT])
grt_ic_df <- cbind(grt_IGT_control, p.adjust(p_IGT_control, method = "bonferroni")[grt_IGT_control])

p_T2D_control2 <- rep(1, m)
p_T2D_IGT2 <- rep(1, m)
p_IGT_control2 <- rep(1, m)
for(j in 1:m){
  #test <- pairwise.wilcox.test(mat_imp[,j], g = study_con, p.adjust.method = "none", alternative = "less")
  test <- pairwise.t.test(mat_imp_si[,j], g = study_con, p.adjust.method = "none", alternative = "less")
  p_T2D_control2[j] = test$p.value[2,1] #T2D < control
  p_T2D_IGT2[j] <- test$p.value[2,2] #T2D < IGT
  p_IGT_control2[j] <- test$p.value[1,1] #control < IGT
}

less_T2D_control <- which(p.adjust(p_T2D_control2, method = "fdr") < 0.1)
less_T2D_IGT <- which(p.adjust(p_T2D_IGT2, method = "fdr") < 0.1)
less_IGT_control <- which(p.adjust(p_IGT_control2, method = "fdr") < 0.1)

less_tc_df <- cbind(less_T2D_control, p.adjust(p_T2D_control2, method = "bonferroni")[less_T2D_control])
less_ti_df <- cbind(less_T2D_IGT, p.adjust(p_T2D_IGT2, method = "bonferroni")[less_T2D_IGT])
less_ic_df <- cbind(less_IGT_control, p.adjust(p_IGT_control2, method = "bonferroni")[less_IGT_control])

record_table_imputed <- rbind(cbind(length(less_tc_df)/2, length(less_ti_df)/2, length(less_ic_df)/2),
                              cbind(length(grt_tc_df)/2, length(grt_ti_df)/2, length(grt_ic_df)/2))
colnames(record_table_imputed) <- c("T2D VS Control", "T2D vs IGT", "IGT vs Control")
rownames(record_table_imputed) <- c("less", "greater")

intersect_table <- rbind(cbind(length(intersect(Karlsson_raw_counts_compare$less_tc_df[,1], less_tc_df[,1])), length(intersect(Karlsson_raw_counts_compare$less_ti_df[,1], less_ti_df[,1])), length(intersect(Karlsson_raw_counts_compare$less_ic_df[,1], less_ic_df[,1]))),
                         cbind(length(intersect(Karlsson_raw_counts_compare$greater_tc_df[,1], grt_tc_df[,1])), length(intersect(Karlsson_raw_counts_compare$greater_ti_df[,1], grt_ti_df[,1])), length(intersect(Karlsson_raw_counts_compare$greater_ic_df[,1], grt_ic_df[,1]))))

write.csv(intersect_table, "summary/si_intersect_table_wct.csv")
write.csv(record_table_imputed, "summary/si_record_table_imputed_wct.csv")

write.csv(colnames(otable)[ grt_tc_df[which(! grt_tc_df[,1] %in%  Karlsson_raw_counts_compare$greater_tc_df[,1]), 1] ]
          , "summary/si_newly_found_taxa_tc.csv")
write.csv(colnames(otable)[ grt_tc_df[which(grt_tc_df[,1] %in% Karlsson_raw_counts_compare$greater_tc_df[,1]),1] ]
          , "summary/si_overlapped_taxa_tc.csv")
write.csv(colnames(otable)[Karlsson_raw_counts_compare$greater_tc_df[,1]], "summary/si_raw_identified_taxa_tc.csv")

write.csv(colnames(otable)[ grt_ti_df[which(! grt_ti_df[,1] %in%  Karlsson_raw_counts_compare$greater_ti_df[,1]), 1] ]
          , "summary/si_newly_found_taxa_ti.csv")
write.csv(colnames(otable)[ grt_ti_df[which(grt_ti_df[,1] %in% Karlsson_raw_counts_compare$greater_ti_df[,1]),1] ]
          , "summary/si_overlapped_taxa_ti.csv")
write.csv(colnames(otable)[Karlsson_raw_counts_compare$greater_ti_df[,1]], "summary/si_raw_identified_taxa_ti.csv")

write.csv(colnames(otable)[ grt_ic_df[which(! grt_ic_df[,1] %in%  Karlsson_raw_counts_compare$greater_ic_df[,1]), 1] ]
          , "summary/si_newly_found_taxa_ic.csv")
write.csv(colnames(otable)[ grt_ic_df[which(grt_ic_df[,1] %in% Karlsson_raw_counts_compare$greater_ic_df[,1]),1] ]
          , "summary/si_overlapped_taxa_ic.csv")
write.csv(colnames(otable)[Karlsson_raw_counts_compare$greater_ic_df[,1]], "summary/si_raw_identified_taxa_ic.csv")

sum(length(grt_T2D_control), length(grt_IGT_control), length(grt_T2D_IGT))
# 203

