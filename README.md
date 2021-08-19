mbImpute: an accurate and robust imputation method for microbiome data
================
Ruochen Jiang, Wei Vivian Li, and Jingyi Jessica Li
2021-08-18

<!-- README.md is generated from README.Rmd. Please edit that file -->

# mbImpute

<!-- badges: start -->
<!-- badges: end -->

The goal of mbImpute is to impute false zero counts in microbiome
sequencing data, i.e., a sample-by-taxon count matrix, by jointly
borrowing information from similar samples, similar taxa and optional
metadata including sample covariates and taxon phylogeny.

## Installation

Please install the following R packages first.

``` r
install.packages("glmnet")
install.packages("devtools")
install.packages("Matrix")
```

Then you can use the following R command to directly install the
mbImpute package from GitHub:

``` r
library(devtools)
install_github("ruochenj/mbImpute/mbImpute R package")
```

## Example

We use the microbiome dataset from Karlsson et al (2013) as an example
to demonstrate the use of mbImpute:

``` r
# Load the R packages
library(mbImpute)
library(glmnet)
#> Loading required package: Matrix
#> Loaded glmnet 4.1-2
library(Matrix)

# Display part of the OTU table
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

# Display part of the taxon phylogenetic distance matrix, whose rows and columns correspond to the columns in otu_tab 
D[1:6, 1:6]
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    2    9   10   10    8
#> [2,]    2    0    9   10   10    8
#> [3,]    9    9    0    3    3    3
#> [4,]   10   10    3    0    2    4
#> [5,]   10   10    3    2    0    4
#> [6,]    8    8    3    4    4    0

# Display part of the (optional) meta data, i.e., the sample covariate matrix with rows representing samples and corresponding to the rows in otu_tab
meta_data[1:6, 1:6]
#>      study_condition      age number_reads triglycerides     hba1c       ldl
#> S112             IGT 1.293993    0.6475183     0.9926486 1.2575721 1.0022890
#> S118         control 2.587987    1.4075527     0.2357540 0.8803004 2.8883168
#> S121         control 1.293993    0.9558486     0.1613054 1.2575721 0.9915117
#> S126         control 1.293993    1.2244933     0.7072621 1.3833293 2.4895566
#> S127         control 2.587987    0.9428587     0.9057918 1.2575721 3.2116358
#> S131             IGT 1.293993    1.0145444     1.8239918 1.5090865 1.7351455
study_condition = meta_data[,1]
meta_data <- as.data.frame(unclass(meta_data))
meta_data <- meta_data[,-1]

# Demo 1: run mbImpute on a single core
imputed_count_mat_list <- mbImpute(condition = study_condition, otu_tab = otu_tab, metadata = meta_data, D = D)
#> [1] 1
#> [1] 0
#> [1] 11230 23280
#> [1] 0.7156728 5.0000000 5.0000000 5.0000000
#> [1] 11230 23280
#> [1] 0.7156728 0.5568270 5.0000000 5.0000000
#> [1] 11230 23280
#> [1] 0.7156728 0.5568270 0.5841656 5.0000000
#> [1] 11230 23280
#> [1] 0.7156728 0.5568270 0.5841656 0.6362045
#> [1] 10000
#> [1] "Finished."

# A glance at the imputed result, which includes three matrices
## The first is an imputed matrix on the log10 scale; we recommend users to perform downstream analysis based on normal distributions on this data, whose values in each taxon (column) follows an approximate normal distribution (see our paper for detail)
imputed_count_mat_list$imp_count_mat_lognorm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                5.335660                        5.153952
#> S118                4.000822                        4.721291
#> S121                4.409327                        5.168918
## The second is an imputed normalized count matrix, where each sample (row) is set to have the same total of a million reads
imputed_count_mat_list$imp_count_mat_norm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                  216599                          142544
#> S118                   10017                           52635
#> S121                   25663                          147541
## The third is an imputed count matrix on the original scale, with each sample (row) having the read count same as that in the original otu_tab
imputed_count_mat_list$imp_count_mat_origlibsize[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954404                         2602397
#> S118                  222608                         1169709
#> S121                  549997                         3162029

# Demo 2: if you have multiple (e.g., 4) cores and would like to do parallel computing
imputed_count_mat_list <- mbImpute(condition = meta_data$study_condition, otu_tab = otu_tab, metadata = meta_data, D = D, parallel = TRUE, ncores = 4)
#> [1] 1
#> [1] 0
#> [1] 11230 23280
#> [1] 0.7156728 5.0000000 5.0000000 5.0000000
#> [1] 11230 23280
#> [1] 0.7156728 0.5568270 5.0000000 5.0000000
#> [1] 11230 23280
#> [1] 0.7156728 0.5568270 0.5841656 5.0000000
#> [1] 11230 23280
#> [1] 0.7156728 0.5568270 0.5841656 0.6362045
#> [1] 10000
#> [1] "Finished."

# A glance at the imputed result, which includes three matrices
## The first is an imputed matrix on the log10 scale; we recommend users to perform downstream analysis based on normal distributions on this data, whose values in each taxon (column) follows an approximate normal distribution (see our paper for detail)
imputed_count_mat_list$imp_count_mat_lognorm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                5.335660                        5.153952
#> S118                4.000822                        4.721291
#> S121                4.409327                        5.168918
## The second is an imputed normalized count matrix, where each sample (row) is set to have the same total of a million reads
imputed_count_mat_list$imp_count_mat_norm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                  216599                          142544
#> S118                   10017                           52635
#> S121                   25663                          147541
## The third is an imputed count matrix on the original scale, with each sample (row) having the read count same as that in the original otu_tab
imputed_count_mat_list$imp_count_mat_origlibsize[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954404                         2602397
#> S118                  222608                         1169709
#> S121                  549997                         3162029

# Demo 3: if you do not have meta data or phylogenetic information, and the samples belong to one condition
otu_tab_T2D <- otu_tab[study_condition == "T2D",]
imputed_count_matrix_list <- mbImpute(otu_tab = otu_tab_T2D)
#> [1] "Meta data information unavailable"
#> [1] "Phylogenentic information unavailable"
#> [1] 4204 3707
#> [1] 1.055828 5.000000 5.000000 5.000000
#> [1] 4204 3707
#> [1] 1.055828 1.053811 5.000000 5.000000
#> [1] 4204 3707
#> [1] 1.055828 1.053811 0.962270 5.000000
#> [1] 4204 3707
#> [1] 1.0558281 1.0538114 0.9622700 0.7518937
#> [1] "Finished."

# A glance at the imputed result, which includes three matrices
## The first is an imputed matrix on the log10 scale; we recommend users to perform downstream analysis based on normal distributions on this data, whose values in each taxon (column) follows an approximate normal distribution (see our paper for detail)
imputed_count_mat_list$imp_count_mat_lognorm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                5.335660                        5.153952
#> S118                4.000822                        4.721291
#> S121                4.409327                        5.168918
## The second is an imputed normalized count matrix, where each sample (row) is set to have the same total of a million reads
imputed_count_mat_list$imp_count_mat_norm[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                  216599                          142544
#> S118                   10017                           52635
#> S121                   25663                          147541
## The third is an imputed count matrix on the original scale, with each sample (row) having the read count same as that in the original otu_tab
imputed_count_mat_list$imp_count_mat_origlibsize[1:3, 1:2]
#>      s__Clostridium_sp_L2_50 s__Faecalibacterium_prausnitzii
#> S112                 3954404                         2602397
#> S118                  222608                         1169709
#> S121                  549997                         3162029
```

Reference:

Karlsson, F. H., Tremaroli, V., Nookaew, I., Bergström, G., Behre, C.
J., Fagerberg, B., … & Bäckhed, F. (2013). Gut metagenome in European
women with normal, impaired and diabetic glucose control. Nature,
498(7452), 99-103.
