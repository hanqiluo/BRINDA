
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRINDA

<!-- badges: start -->
<!-- badges: end -->

The BRINDA R package is a user-friendly all-in-one R package that uses a
series of functions to implement BRINDA adjustment
[method](https://brinda-nutrition.org/publications/) .

## Installation

You can install the development version of BRINDA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hanqiluo/BRINDA")

```

## Example

This is a basic example which shows you how to adjust retinol binding
protein, retinol, ferritin, soluble transferrin receptor, and zinc using
AGP and/or CRP.

``` r
library(BRINDA)
## basic example code
data(sample_data)

sample_data_adj <- BRINDA(dataset = sample_data,
       retinol_binding_protein_varname = rbp,
       retinol_varname = sr, ferritin_varname = sf,
       soluble_transferrin_receptor_varname = stfr,
       zinc_varname = zinc, crp_varname = crp,
       agp_varname = agp, population = Psc,
       crp_ref_value_manual = ,
       agp_ref_value_manual = ,
       output_format = )
#> [1] "*******************************************"
#> [1] "Initial data checks completed"
#> [1] "*******************************************"
#> [1] "                                           "
#> [1] "*******************************************"
#> [1] "** Overview of the dataset and BRINDA package inputs"
#> [1] "** Dataset Name: sample_data"
#> [1] "** Retinol Binding Protein Variable Name: rbp (n = 74)"
#> [1] "** Retinol Variable Name: sr (n = 76)"
#> [1] "** Ferritin Variable Name: sf (n = 74)"
#> [1] "** Soluble Transferrin Receptor Variable Name: stfr (n = 74)"
#> [1] "** Zinc Variable Name: zinc (n = 74)"
#> [1] "** AGP Variable Name: agp (n = 74)"
#> [1] "** CRP Variable Name: crp (n = 74)"
#> [1] "** Population Group: PSC"
#> [1] "** Output Format: SIMPLE"
#> [1] "*******************************************"
#> [1] "                                           "
#> [1] "*******************************************"
#> [1] "Proceed to the BRINDA adjustment"
#> [1] "*******************************************"
#> [1] "*******************************************"
#> [1] "Generated deciles of AGP/CRP based on inputs"
#> [1] "*******************************************"
#> [1] "-------------------------------------------"
#> [1] "Adjusting Retinol Binding Protein using both AGP and CRP"
#> [1] "-------------------------------------------"
#> [1] "-------------------------------------------"
#> [1] "Adjusting Retinol using both AGP and CRP"
#> [1] "-------------------------------------------"
#> [1] "-------------------------------------------"
#> [1] "Adjusting Ferritin using both AGP and CRP"
#> [1] "-------------------------------------------"
#> [1] "-------------------------------------------"
#> [1] "Adjusting Soluble Transferrin Receptor using AGP only"
#> [1] "--------------------------------------------"
#> [1] "-------------------------------------------"
#> [1] "Analyzing Serum Zinc"
#> [1] "No or weak correlation between Serum Zinc and AGP based on Spearman correlation measures"
#> [1] "No or weak correlation between Serum Zinc and CRP based on Spearman correlation measures"
#> [1] "BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and AGP/CRP"
#> [1] "Adjusted zinc values are equal to unadjusted zinc values"
#> [1] "--------------------------------------------"
#> [1] "*******************************************"
#> [1] "BRINDA adjustment completed"
#> [1] "Proceed to output dataset"
#> [1] "*******************************************"
#> [1] "-------------------------------------------"
#> [1] "variables rbp_adj, sr_adj, sf_adj, stfr_adj, zn_adj are generated by the BRINDA function"
#> [1] "-------------------------------------------"
#> [1] "*******************************************"
#> [1] "BRINDA adjustment function complete"
#> [1] "*******************************************"
```
