## ------------------------------------------------------------------------
library(mbImpute)
# the OTU table
otu_tab[1:6, 1:6]
# the taxa distance matrix generated from phylogenetic tree 
D[1:6, 1:6]
# a numeric meta data corresponding to the otu table
meta_data[1:6, 1:6]
# get the condition from the meta data
condition = meta_data$study_condition
cond <- as.numeric(as.factor(condition))
meta_data[,1] <- as.numeric(as.factor(meta_data[,1]))
meta_data <- meta_data[,-1]

imputed_count_mat_lognorm <- mbImpute(condition = condition, otu_tab = otu_tab, meta_data = meta_data, D = D)$imp_count_mat_lognorm
# imputed_matrix <- mbImpute(condition = condition, otu_tab = otu_tab, meta_data = meta_data, D = D, parallel = TRUE, ncores = 4)

## ---- fig.align = "center", out.width='\\textwidth', fig.height = 8, fig.width = 8----
library(ggplot2)
# pca plot
raw_pca_out_full <- prcomp(otu_tab, center = TRUE)
df1 <- data.frame(raw_pca_out_full$x[,1:2], condition)
df1$condition <- as.factor(condition)
ggplot(data = df1, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 2) +
  scale_color_manual("condition", values = c("IGT" = "#DA5E03", "T2D" = "#389E78", "control"= "#2166ac"))

imp_pca_out_full <- prcomp(imputed_count_mat_lognorm, center = TRUE)
df2 <- data.frame(imp_pca_out_full$x[,1:2], condition)
df2$condition <- as.factor(condition)
ggplot(data = df2, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 2) +
  scale_color_manual("condition", values = c("IGT" = "#DA5E03", "T2D" = "#389E78", "control"= "#2166ac"))

## ---- fig.align = "center", out.width='\\textwidth', fig.height = 8, fig.width = 8----
# The histogram for each taxon after imputation. 
par(mfrow = c(3,3))
hist(imputed_count_mat_lognorm[condition == "T2D", 1], main = "T2D imputed 1")
hist(imputed_count_mat_lognorm[condition == "T2D", 2], main = "T2D imputed 2")
hist(imputed_count_mat_lognorm[condition == "T2D", 3], main = "T2D imputed 3")

hist(imputed_count_mat_lognorm[condition == "control", 1], main = "control imputed 1")
hist(imputed_count_mat_lognorm[condition == "control", 2], main = "control imputed 2")
hist(imputed_count_mat_lognorm[condition == "control", 3], main = "control imputed 3")


hist(imputed_count_mat_lognorm[condition == "IGT", 1], main = "IGT imputed 1")
hist(imputed_count_mat_lognorm[condition == "IGT", 2], main = "IGT imputed 2")
hist(imputed_count_mat_lognorm[condition == "IGT", 3], main = "IGT imputed 3")

