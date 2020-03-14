mbImpute: an accurate and robust imputation method for microbiome data
================
Ruochen Jiang, Wei Vivian Li, and Jingyi Jessica Li
2020-03-14

<!-- README.md is generated from README.Rmd. Please edit that file -->

# mbImpute

<!-- badges: start -->

<!-- badges: end -->

The goal of mbImpute is to impute false zero counts in microbiome
sequencing data, i.e., a sample-by-taxon count matrix, by jointly
borrowing information from similar samples, similar taxa and optional
metadata including sample covariates, and taxon phylogeny.

## Installation

You can use the following command in R to directly install the mbImpute
package from GitHub:

``` r
#Please install the following R packages
#install.pacakges("glmnet")
#install.packages("devtools")
#Then install the mbImpute package 
library(devtools)
install_github("ruochenj/mbImpute/mbImpute R package")
```

## Example

We use the microbiome dataset from Karlsson et al (2013) as an example
to demonstrate the use of mbImpute:

``` r
#Load the R packages
library(mbImpute)
library(glmnet)
#> Loading required package: Matrix
#> Loading required package: foreach
#> Loaded glmnet 2.0-18
# the OTU table
otu_tab[1:6, 1:6]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954419                         2602398
#> S118                       0                         1169731
#> S121                  550000                         3162050
#> S126                       0                          986563
#> S127                       0                         3940520
#> S131                       0                          502030
#>      s__Dialister_invisus s__Dorea_longicatena s__Ruminococcus_obeum
#> S112              1671620              1440718               1112728
#> S118              1412991               623343                190968
#> S121               827400               855614                274969
#> S126              2411487               163956                112493
#> S127                    0               554154                286616
#> S131                    0               159910                802850
#>      s__Coprococcus_comes
#> S112               947723
#> S118                58660
#> S121               281024
#> S126               148430
#> S127               210698
#> S131               450721
# the taxon phylogenetic distance matrix 
D[1:6, 1:6]
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    2    9   10   10    8
#> [2,]    2    0    9   10   10    8
#> [3,]    9    9    0    3    3    3
#> [4,]   10   10    3    0    2    4
#> [5,]   10   10    3    2    0    4
#> [6,]    8    8    3    4    4    0
# the (optional) meta data, i.e., sample covariate matrix with rows representing samples and corresponding to the rows in otu_tab
meta_data[1:6, 1:6]
#>      study_condition      age number_reads triglycerides     hba1c
#> S112             IGT 1.293993    0.6475183     0.9926486 1.2575721
#> S118         control 2.587987    1.4075527     0.2357540 0.8803004
#> S121         control 1.293993    0.9558486     0.1613054 1.2575721
#> S126         control 1.293993    1.2244933     0.7072621 1.3833293
#> S127         control 2.587987    0.9428587     0.9057918 1.2575721
#> S131             IGT 1.293993    1.0145444     1.8239918 1.5090865
#>            ldl
#> S112 1.0022890
#> S118 2.8883168
#> S121 0.9915117
#> S126 2.4895566
#> S127 3.2116358
#> S131 1.7351455
# obtain the sample conditions from the meta data (imputation will be performed within each condition)
condition <- meta_data$study_condition

# For all the categorical variables, make sure they are converted to numerical.
meta_data[,1] <- as.numeric(as.factor(meta_data[,1]))

# run mbImpute
imputed_count_mat_list <- mbImpute(condition = condition, otu_tab = otu_tab, meta_data = meta_data, D = D)
#> [1] "condition IGT is imputing"
#> [1] "Working on it!"
#> [1] "condition control is imputing"
#> [1] "Working on it!"
#> [1] "condition T2D is imputing"
#> [1] "Working on it!"
#> [1] "Finished."
# a glance at the imputed matrix
# First result is the count matrix at log scale, we recommend to perform downstream analysis on this data as the distribution for the values in each taxon follows approximately normal distribution (see our paper for more results)
imputed_count_mat_list$imp_count_mat_lognorm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                5.335660                        5.153952
#> S118                4.670297                        4.721291
#> S121                4.409327                        5.168918
# Second result is the imputed normalized count matrix with same library size set to 10^6 for each sample (subject/person).
imputed_count_mat_list$imp_count_mat_norm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                  216599                          142544
#> S118                   46804                           52635
#> S121                   25663                          147541
# Third result is the imputed count matrix with the same scale as the input count matrix.
imputed_count_mat_list$imp_count_mat_origlibsize[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954404                         2602397
#> S118                 1040127                         1169709
#> S121                  549997                         3162029
# If you have multiple cores and would like to do parallel computing, please use the following command
imputed_count_mat_list <- mbImpute(condition = condition, otu_tab = otu_tab, meta_data = meta_data, D = D, parallel = TRUE, ncores = 4)
#> [1] "condition IGT is imputing"
#> [1] "Working on it!"
#> [1] "condition control is imputing"
#> [1] "Working on it!"
#> [1] "condition T2D is imputing"
#> [1] "Working on it!"
#> [1] "Finished."
# a glance at the imputed matrix
# First result is the count matrix at log scale, we recommend to perform downstream analysis on this data as the distribution for the values in each taxon follows approximately normal distribution (see our paper for more results)
imputed_count_mat_list$imp_count_mat_lognorm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                5.335660                        5.153952
#> S118                4.670297                        4.721291
#> S121                4.409327                        5.168918
# Second result is the imputed normalized count matrix with same library size set to 10^6 for each sample (subject/person).
imputed_count_mat_list$imp_count_mat_norm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                  216599                          142544
#> S118                   46804                           52635
#> S121                   25663                          147541
# Third result is the imputed count matrix with the same scale as the input count matrix.
imputed_count_mat_list$imp_count_mat_origlibsize[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954404                         2602397
#> S118                 1040127                         1169709
#> S121                  549997                         3162029
# If you do not have meta data, or phylogenetic information, and the samples belong to one condition
otu_tab_T2D <- otu_tab[condition == "T2D",]
imputed_count_matrix_list <- mbImpute(otu_tab = otu_tab_T2D)
#> [1] "Meta data information unavailable"
#> [1] "Phylogenentic information unavailable"
#> [1] "condition 1 is imputing"
#> [1] "Working on it!"
#> [1] "Finished."
# a glance at the imputed matrix
# First result is the count matrix at log scale, we recommend to perform downstream analysis on this data as the distribution for the values in each taxon follows approximately normal distribution (see our paper for more results)
imputed_count_mat_list$imp_count_mat_lognorm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                5.335660                        5.153952
#> S118                4.670297                        4.721291
#> S121                4.409327                        5.168918
# Second result is the imputed normalized count matrix with same library size set to 10^6 for each sample (subject/person).
imputed_count_mat_list$imp_count_mat_norm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                  216599                          142544
#> S118                   46804                           52635
#> S121                   25663                          147541
# Third result is the imputed count matrix with the same scale as the input count matrix.
imputed_count_mat_list$imp_count_mat_origlibsize[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954404                         2602397
#> S118                 1040127                         1169709
#> S121                  549997                         3162029
```

Reference:

Karlsson, F. H., Tremaroli, V., Nookaew, I., Bergström, G., Behre, C.
J., Fagerberg, B., … & Bäckhed, F. (2013). Gut metagenome in European
women with normal, impaired and diabetic glucose control. Nature,
498(7452), 99-103.
