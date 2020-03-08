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

#get the distance matrix D base on phylogenetic tree
edges <- Karlsson_taxa_phy$edge

#generate phylogenetic tree
edges <- pseq2@phy_tree$edge

m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
#199, 449
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m

#tree_height = floor(log2(n_taxa))+1
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
  l <- rep(1,tree_height+1)
  l[1] = i
  for(j in 2:(tree_height+1)){
    if(sum(edges[,2] %in% l[j-1]) != 0){
      l[j] = edges[edges[,2] %in% l[j-1], 1]
    }
    else{
      l[j] = NA
    }
  }
  nd_mat[,i] = l[2:(tree_height+1)]
}


d1_mat <- matrix(0, nrow = m, ncol = m)
#records the position of 1:m in the edges set.
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0
#d1_mat denotes the distance for two taxa

#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

write.csv(otu_tab_processed, "otu_real_data.csv")
write.csv(x, "meta_data.csv")
write.csv(d1_mat, "D.csv")

otable <- read.csv("otu_real_data.csv", row.names = "X")
met_data <- read.csv("meta_data.csv", row.names = "X")
D1 <- read.csv("D.csv", row.names = "X")

