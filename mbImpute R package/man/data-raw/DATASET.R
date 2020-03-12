## code to prepare `DATASET` dataset goes here

#process otu table and metagenomic data -- type 2 diabetes
library(curatedMetagenomicData)
library(mbImpute)

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

D <- distance_mat_gen(edges, n_taxa = m)
otu_tab <- OTU_tab@.Data

meta_data <- data.frame(chosen_meta)

usethis::use_data(otu_tab)
usethis::use_data(D)
usethis::use_data(meta_data)
usethis::use_data(edges)


