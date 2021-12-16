
[![DOI](https://sandbox.zenodo.org/badge/414102388.svg)](https://sandbox.zenodo.org/badge/latestdoi/414102388)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRINDA

<!-- badges: start -->
<!-- badges: end -->

The BRINDA R package is a user-friendly all-in-one R package that uses a
series of functions to implement BRINDA adjustment
[method](https://brinda-nutrition.org/publications/) .

## Installation

``` r
install.packages("BRINDA")
```

### Development version

Alternatively, you can install the development version of BRINDA from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("hanqiluo/BRINDA")
```

## Example

This is a basic example which shows you how to use the BRINDA package to
adjust retinol binding protein, retinol, ferritin, soluble transferrin
receptor, and zinc using AGP and/or CRP.

``` r
library(BRINDA)
## basic example code
data(sample_data)

sample_data_adj <- BRINDA(dataset = sample_data,
       retinol_binding_protein_varname = rbp,
       retinol_varname = sr, 
       ferritin_varname = sf,
       soluble_transferrin_receptor_varname = stfr,
       zinc_varname = zinc, 
       crp_varname = crp,
       agp_varname = agp, 
       population = Psc,
       crp_ref_value_manual = ,
       agp_ref_value_manual = ,
       output_format = )
#> -------------------------------------------
#> ** Initial data checks completed **
#> -------------------------------------------
#> ** Overview of the dataset and BRINDA package inputs **
#> **** Dataset Name: sample_data**
#> **** Retinol Binding Protein Variable Name: rbp (n = 74)
#> **** Retinol Variable Name: sr (n = 76)
#> **** Ferritin Variable Name: sf (n = 74)
#> **** Soluble Transferrin Receptor Variable Name: stfr (n = 74)
#> **** Zinc Variable Name: zinc (n = 74)
#> **** AGP Variable Name: agp (n = 74)
#> **** CRP Variable Name: crp (n = 74)
#> **** Population Group: PSC
#> **** Output Format: SIMPLE
#> -------------------------------------------
#> ** Generate deciles of AGP/CRP based on inputs **
#> **** log-AGP = -0.52
#> **** log-CRP = -2.26
#> -------------------------------------------
#> ** Proceed to the BRINDA adjustment **
#> **** Adjusting Retinol Binding Protein using both AGP and CRP
#> **** Adjusting Retinol using both AGP and CRP
#> **** Adjusting Ferritin using both AGP and CRP
#> **** Adjusting Soluble Transferrin Receptor using AGP only
#> **** Adjusted zinc values are equal to unadjusted zinc values
#> ****** No or weak correlation between Serum Zinc and AGP based on Spearman correlation measures
#> ****** No or weak correlation between Serum Zinc and CRP based on Spearman correlation measures
#> ****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and AGP/CRP
#> ** BRINDA adjustment completed **
#> -------------------------------------------
#> ** Proceed to output dataset **
#> variables rbp_adj, sr_adj, sf_adj, stfr_adj, zn_adj are generated by the BRINDA function
#> -------------------------------------------
#> ** BRINDA adjustment function complete **
#> -------------------------------------------
```

## Citation

Luo, H.; Addo, Y.; Jahan, A (2021) BRINDA: Computation of BRINDA
Adjusted Micronutrient Biomarkers for Inflammation. R package version
0.1.3, <https://CRAN.R-project.org/package=BRINDA>

## Contributing

If you would like to report bugs, suggest features, or leave comments,
please [create issues](https://github.com/hanqiluo/BRINDA/issues).

If you would like to contribute code, please fork the source code,
modify, and issue a [pull
request](https://github.com/hanqiluo/BRINDA/pulls).

## Acknowledgements

The Authors thank all members of the BRINDA working group who helped in
developing the BRINDA adjustment method. The authors also thank Charles
D. Arnold, Fanny Sandalinas, Kevin Tang, and Lucas Gosdin for the
extensive testing of the package, Joanne Arsenault and Christine
McDonald for their editing and comments, and Jae Yeon Kim for his
assistance with CRAN submission.
