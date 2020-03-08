#Separate conditions and feed into mbImpute

study_con = x$study_condition
sum(study_con == "T2D")
sum(study_con == "control")
sum(study_con == "IGT")

n = length(study_con)

T2D_idx <- which(study_con == "T2D")
control_idx <- which(study_con == "control")
IGT_idx <- which(study_con == "IGT")

otable <- read.csv("otu_real_data.csv", row.names = "X")
met_data <- read.csv("meta_data.csv", row.names = "X")
D <- read.csv("D.csv", row.names = "X")

m = dim(otable)[2]
n = dim(otable)[1]

met_data <- met_data[,-2]

write.csv(otable[T2D_idx,], "otu_real_data_T2D.csv")
write.csv(met_data[T2D_idx,], "meta_data_T2D.csv")

write.csv(otable[control_idx,], "otu_real_data_control.csv")
write.csv(met_data[control_idx,], "meta_data_control.csv")

write.csv(otable[IGT_idx,], "otu_real_data_IGT.csv")
write.csv(met_data[IGT_idx,], "meta_data_IGT.csv")

#feed into the server

T2D_mat <- readRDS("RDatas/imputed2_T2D_mat.rds")
control_mat <- readRDS("RDatas/imputed1_control_mat.rds")
IGT_mat <- readRDS("RDatas/imputed1_IGT_mat.rds")


#regroup the conditions
mat_imp <- matrix(1, ncol = m, nrow = n)
mat_imp[T2D_idx,] = as.matrix(T2D_mat)
mat_imp[control_idx,] = as.matrix(control_mat)
mat_imp[IGT_idx,] = as.matrix(IGT_mat)

########imputation2#############
imp_mat2 <- mat_imp

set.seed(42) # Set a seed if you want reproducible results
library(Rtsne)
imp_rtsne_out <- Rtsne(imp_mat2) # Run TSNE

imp_mat2[imp_mat2 > 6] = 5.5

png("plots/rtsne_imp.png")
plot(imp_rtsne_out$Y, col = as.factor(study_con))
dev.off()
library(tsne)
imp_tsne_out <- tsne(imp_mat2)
png("plots/tsne_imp.png")
plot(imp_tsne_out, col = as.factor(study_con))
dev.off()
imp_pca_out <- prcomp(imp_mat2, center = TRUE)
par(mfrow = c(1,1))
png("plots/pca_1_2_imp.png")
plot(imp_pca_out$x[,1:2], col = as.factor(study_con))
dev.off()
png("plots/pca_1_3_imp2.png")
plot(imp_pca_out$x[,c(1,3)], col = as.factor(study_con))
dev.off()
png("plots/pca_2_3_imp.png")
plot(imp_pca_out$x[,2:3], col = as.factor(study_con))
dev.off()

png("plots/htmap_imp.png")
heatmap(as.matrix(imp_mat2), Colv = NA, Rowv = NA, scale="row")
dev.off()

png("plots/htmap_imp_ordered.png")
heatmap(as.matrix(imp_mat2), scale = "row")
dev.off()

png("plots/htmap_imp_dis.png")
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
  test <- pairwise.t.test(mat_imp[,j], g = study_con, p.adjust.method = "none", alternative = "greater")
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
  test <- pairwise.t.test(mat_imp[,j], g = study_con, p.adjust.method = "none", alternative = "less")
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

write.csv(intersect_table, "summary/intersect_table_wct.csv")
write.csv(record_table_imputed, "summary/record_table_imputed_wct.csv")

write.csv(colnames(otable)[ grt_tc_df[which(! grt_tc_df[,1] %in%  Karlsson_raw_counts_compare$greater_tc_df[,1]), 1] ]
          , "summary/newly_found_taxa_tc.csv")
write.csv(colnames(otable)[ grt_tc_df[which(grt_tc_df[,1] %in% Karlsson_raw_counts_compare$greater_tc_df[,1]),1] ]
          , "summary/overlapped_taxa_tc.csv")
write.csv(colnames(otable)[Karlsson_raw_counts_compare$greater_tc_df[,1]], "summary/raw_identified_taxa_tc.csv")

write.csv(colnames(otable)[ grt_ti_df[which(! grt_ti_df[,1] %in%  Karlsson_raw_counts_compare$greater_ti_df[,1]), 1] ]
          , "summary/newly_found_taxa_ti.csv")
write.csv(colnames(otable)[ grt_ti_df[which(grt_ti_df[,1] %in% Karlsson_raw_counts_compare$greater_ti_df[,1]),1] ]
          , "summary/overlapped_taxa_ti.csv")
write.csv(colnames(otable)[Karlsson_raw_counts_compare$greater_ti_df[,1]], "summary/raw_identified_taxa_ti.csv")

write.csv(colnames(otable)[ grt_ic_df[which(! grt_ic_df[,1] %in%  Karlsson_raw_counts_compare$greater_ic_df[,1]), 1] ]
          , "summary/newly_found_taxa_ic.csv")
write.csv(colnames(otable)[ grt_ic_df[which(grt_ic_df[,1] %in% Karlsson_raw_counts_compare$greater_ic_df[,1]),1] ]
          , "summary/overlapped_taxa_ic.csv")
write.csv(colnames(otable)[Karlsson_raw_counts_compare$greater_ic_df[,1]], "summary/raw_identified_taxa_ic.csv")

sum(length(grt_T2D_control), length(grt_IGT_control), length(grt_T2D_IGT))
# 117