# zinb
setwd("~/Dropbox/mbimpute/code/Real data Karlsson")

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
# otable <- read.csv("otu_real_data.csv", row.names = "X")
meta_data <- read.csv("meta_data.csv", row.names = "X")
D1 <- read.csv("D.csv", row.names = "X")

study_con = meta_tab$study_condition
Vardat <- meta_tab[which(study_con == "control" | study_con == "T2D"),]
OTUdat <- otable[which(study_con == "control" | study_con == "T2D"),]
condition = study_con[which(study_con == "control" | study_con == "T2D")]

cov_mat <- Vardat[, -c(1,2)]

library(pscl)

df_taxa <- c()
p_value_zinb <- c()
for(i in 1:dim(OTUdat)[2]){
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
    p_val <- sr$coefficients$count[2,4]
  }
  p_value_zinb <- c(p_value_zinb, p_val)
}

p_value_zinb_rec <- p_value_zinb
p_value_zinb[is.na(p_value_zinb_rec)] = 1
p.adjust(p_value_zinb, method = "fdr")
colnames(OTUdat)[which(p.adjust(p_value_zinb, method = "fdr") < 0.1)]
# colnames(OTUdat)[which(p_value_zinb < 0.1/300)]
write.csv( colnames(OTUdat)[which(p.adjust(p_value_zinb, method = "fdr") < 0.1)], "zinb_results.csv")

# ident <- which(p_value_zinb < 0.05/300)
# sum(ident %in% less_identified_full)
# sum(ident %in% greater_identified_full)
# 
# zinb_result <- c((29+24)/(84), (29+24)/(85+95))
