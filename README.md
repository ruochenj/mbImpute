
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mbImpute

<!-- badges: start -->

<!-- badges: end -->

The goal of mbImpute is to …

## Installation

You can install the released version of mbImpute using:

``` r
install.packages("mbImpute_0.1.0.tar.gz")
```

## Example

This is a basic example which shows you how to use mbImpute to perform
microbiome imputation:

``` r
library(mbImpute)
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 3.5.2

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
# the taxa distance matrix generated from phylogenetic tree 
D[1:6, 1:6]
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    2    9   10   10    8
#> [2,]    2    0    9   10   10    8
#> [3,]    9    9    0    3    3    3
#> [4,]   10   10    3    0    2    4
#> [5,]   10   10    3    2    0    4
#> [6,]    8    8    3    4    4    0
# a numeric meta data corresponding to the otu table
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
# get the condition from the meta data
condition = meta_data$study_condition
cond <- as.numeric(as.factor(condition))
meta_data[,1] <- as.numeric(as.factor(meta_data[,1]))
meta_data <- meta_data[,-1]

imputed_matrix <- mbImpute(condition = condition, otu_tab = otu_tab, meta_data = meta_data, D = D, k =5)
#> [1] "condition IGT is imputed"
#> [1]  49 344
#> [1] 0.4181722
#> [1] "design_mat generated"
#> [1] 3851
#> [1] 3934 3851
#> [1] 3934 3851
#> [1] 0.7068284
#> [1] 0.007268951
#> [1] "the preserved values in zero inflated matrix is: "
#> [1] 0.2185556
#> [1]  49 344
#> [1]  49 344
#> [1] "the mse for psi "
#> [1] 2
#> [1] "is"
#> [1] 0.7068284
#> [1] 0.4181548
#> [1]  49 344
#> [1] "impute_val generated"
#> [1]  49 344
#> [1] "condition control is imputed"
#> [1]  43 344
#> [1] 0.4879396
#> [1] "design_mat generated"
#> [1] 3079
#> [1] 3357 3079
#> [1] 3357 3079
#> [1] 0.6578296
#> [1] 0.02175325
#> [1] "the preserved values in zero inflated matrix is: "
#> [1] 0.1865
#> [1]  43 344
#> [1]  43 344
#> [1] "the mse for psi "
#> [1] 2
#> [1] "is"
#> [1] 0.6578296
#> [1] 0.487936
#> [1]  43 344
#> [1] "impute_val generated"
#> [1]  43 344
#> [1] "condition T2D is imputed"
#> [1]  53 344
#> [1] 0.4231258
#> [1] "design_mat generated"
#> [1] 4304
#> [1] 4328 4304
#> [1] 4328 4304
#> [1] 0.7480673
#> [1] 0.008594657
#> [1] "the preserved values in zero inflated matrix is: "
#> [1] 0.2404444
#> [1]  53 344
#> [1]  53 344
#> [1] "the mse for psi "
#> [1] 2
#> [1] "is"
#> [1] 0.7480673
#> [1] 0.4231108
#> [1]  53 344
#> [1] "impute_val generated"
#> [1]  53 344
# imputed_matrix <- mbImpute(condition = condition, otu_tab = otu_tab, meta_data = meta_data, D = D, k =5, parallel = TRUE, ncores = 4)
```

Visualize the imputation results

``` r
library(ggplot2)
# pca plot
raw_pca_out_full <- prcomp(otu_tab, center = TRUE)
df1 <- data.frame(raw_pca_out_full$x[,1:2], condition)
df1$condition <- as.factor(condition)
ggplot(data = df1, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 2) +
  scale_color_manual("condition", values = c("IGT" = "#DA5E03", "T2D" = "#389E78", "control"= "#2166ac"))
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r

imp_pca_out_full <- prcomp(imputed_matrix, center = TRUE)
df2 <- data.frame(imp_pca_out_full$x[,1:2], condition)
df2$condition <- as.factor(condition)
ggplot(data = df2, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 2) +
  scale_color_manual("condition", values = c("IGT" = "#DA5E03", "T2D" = "#389E78", "control"= "#2166ac"))
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

``` r
# The histogram for each taxon after imputation. 
hist(imputed_matrix[condition == "T2D", 1])
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
hist(imputed_matrix[condition == "T2D", 2])
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r
hist(imputed_matrix[condition == "T2D", 3])
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />

``` r

hist(imputed_matrix[condition == "control", 1])
```

<img src="man/figures/README-unnamed-chunk-3-4.png" width="100%" />

``` r
hist(imputed_matrix[condition == "control", 2])
```

<img src="man/figures/README-unnamed-chunk-3-5.png" width="100%" />

``` r
hist(imputed_matrix[condition == "control", 3])
```

<img src="man/figures/README-unnamed-chunk-3-6.png" width="100%" />

``` r

hist(imputed_matrix[condition == "IGT", 1])
```

<img src="man/figures/README-unnamed-chunk-3-7.png" width="100%" />

``` r
hist(imputed_matrix[condition == "IGT", 2])
```

<img src="man/figures/README-unnamed-chunk-3-8.png" width="100%" />

``` r
hist(imputed_matrix[condition == "IGT", 3])
```

<img src="man/figures/README-unnamed-chunk-3-9.png" width="100%" />
