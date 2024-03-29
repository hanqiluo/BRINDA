
[![DOI](https://zenodo.org/badge/414102388.svg)](https://zenodo.org/badge/latestdoi/414102388)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRINDA

<!-- badges: start -->

[![R-CMD-check](https://github.com/hanqiluo/BRINDA/workflows/R-CMD-check/badge.svg)](https://github.com/hanqiluo/BRINDA/actions)
<!-- badges: end -->

The BRINDA R package is a user-friendly all-in-one R package that uses a
series of functions to implement BRINDA adjustment
[method](https://www.brinda-nutrition.org/publications/) .

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

#
# Example 1: Preschool-age children (PSC) -----------------------------
# calculate BRINDA inflammation adjustment values for preschool-age children (Assuming the dataset contain information of preschool-age children)

#
# Example 1.1 PSC: Both AGP and CRP are available 
#
final_data_psc <- BRINDA(
    # Enter the name of the dataset 
    dataset = sample_data, 
    # Enter the variable name of retinol binding protein in your dataset (if available)
    retinol_binding_protein_varname = rbp, 
    # Enter the variable name of serum/plasma retinol in your dataset (if available)
    retinol_varname = sr,
    # Enter the variable name of serum/plasma ferritin in your dataset (if available)
    ferritin_varname = sf,
    # Enter the variable name of serum/plasma soluble transferrin receptor in your dataset (if available)
    soluble_transferrin_receptor_varname = stfr, 
    # Enter the variable name of serum/plasma zinc in your dataset (if available)
    zinc_varname = zinc,
    # Enter the variable name of CRP (unit must be mg/L) in your dataset (if available)
    crp_varname = crp, 
    # Enter the variable name of AGP (unit must be g/L) in your dataset (if available)
    agp_varname = agp,
    # Please write WRA, PSC, Other, or Manual. 
population = Psc, 
# leave crp_ref_value_manual empty, BRINDA R package will use an external crp  reference value for PSC
    crp_ref_value_manual = , 
    # leave agp_ref_value_manual empty, BRINDA R package will use an external agp reference value for PSC
    agp_ref_value_manual = , 
    # Please write FULL or SIMPLE (default output is SIMPLE, if users leave output_format empty)
output_format = full)
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
#> **** Output Format: FULL
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
#> variables log_agp_ref, log_crp_ref, log_agp, log_agp_diff, log_crp, log_crp_diff, zn_agp_cor, zn_agp_P_value, zn_crp_cor, zn_crp_P_value, rbp_adj, rbp_beta1, rbp_beta1_se, rbp_beta1_P_value, rbp_beta2, rbp_beta2_se, rbp_beta2_P_value, sr_adj, sr_beta1, sr_beta1_se, sr_beta1_P_value, sr_beta2, sr_beta2_se, sr_beta2_P_value, sf_adj, sf_beta1, sf_beta1_se, sf_beta1_P_value, sf_beta2, sf_beta2_se, sf_beta2_P_value, stfr_adj, stfr_beta1, stfr_beta1_se, stfr_beta1_P_value, zn_adj are generated by the BRINDA function
#> -------------------------------------------
#> ** BRINDA adjustment function complete **
#> -------------------------------------------

#
# Example 1.2 PSC: Only AGP is available 
#
final_data_psc <- BRINDA(
    # Enter the name of the dataset 
    dataset = sample_data, 
    # Enter the variable name of retinol binding protein in your dataset (if available)
    retinol_binding_protein_varname = rbp, 
    # Enter the variable name of serum/plasma retinol in your dataset (if available)
    retinol_varname = sr,
    # Enter the variable name of serum/plasma ferritin in your dataset (if available)
    ferritin_varname = sf,
    # Enter the variable name of serum/plasma soluble transferrin receptor in your dataset (if available)
    soluble_transferrin_receptor_varname = stfr, 
    # Enter the variable name of serum/plasma zinc in your dataset (if available)
    zinc_varname = zinc,
    # Enter the variable name of CRP (unit must be mg/L) in your dataset (if available)
    # crp_varname = crp, 
    # Enter the variable name of AGP (unit must be g/L) in your dataset (if available)
    agp_varname = agp,
    # Please write WRA, PSC, Other, or Manual. 
population = Psc, 
# leave crp_ref_value_manual empty, BRINDA R package will use an external crp  reference value for PSC
    crp_ref_value_manual = , 
    # leave agp_ref_value_manual empty, BRINDA R package will use an external agp reference value for PSC
    agp_ref_value_manual = , 
    # Please write FULL or SIMPLE (default output is SIMPLE, if users leave output_format empty)
output_format = full)
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
#> **** CRP: NA
#> **** Population Group: PSC
#> **** Output Format: FULL
#> -------------------------------------------
#> ** Generate deciles of AGP/CRP based on inputs **
#> **** log-AGP = -0.52
#> -------------------------------------------
#> ** Proceed to the BRINDA adjustment **
#> **** Adjusting Retinol Binding Protein using AGP only
#> **** Adjusting Retinol using AGP only
#> **** Adjusting Ferritin using AGP only
#> **** Adjusting Soluble Transferrin Receptor using AGP only
#> **** Adjusted zinc values are equal to unadjusted zinc values
#> ****** No or weak correlation between Serum Zinc and AGP based on Spearman correlation measures
#> ****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and AGP
#> ** BRINDA adjustment completed **
#> -------------------------------------------
#> ** Proceed to output dataset **
#> variables log_agp_ref, log_agp, log_agp_diff, zn_agp_cor, zn_agp_P_value, rbp_adj, rbp_beta1, rbp_beta1_se, rbp_beta1_P_value, sr_adj, sr_beta1, sr_beta1_se, sr_beta1_P_value, sf_adj, sf_beta1, sf_beta1_se, sf_beta1_P_value, stfr_adj, stfr_beta1, stfr_beta1_se, stfr_beta1_P_value, zn_adj are generated by the BRINDA function
#> -------------------------------------------
#> ** BRINDA adjustment function complete **
#> -------------------------------------------

#
# Example 1.3 PSC: Only CRP is available 
#
final_data_psc <- BRINDA(
    # Enter the name of the dataset 
    dataset = sample_data, 
    # Enter the variable name of retinol binding protein in your dataset (if available)
    retinol_binding_protein_varname = rbp, 
    # Enter the variable name of serum/plasma retinol in your dataset (if available)
    retinol_varname = sr,
    # Enter the variable name of serum/plasma ferritin in your dataset (if available)
    ferritin_varname = sf,
    # Enter the variable name of serum/plasma soluble transferrin receptor in your dataset (if available)
    soluble_transferrin_receptor_varname = stfr, 
    # Enter the variable name of serum/plasma zinc in your dataset (if available)
    zinc_varname = zinc,
    # Enter the variable name of CRP (unit must be mg/L) in your dataset (if available)
     crp_varname = crp, 
    # Enter the variable name of AGP (unit must be g/L) in your dataset (if available)
    # agp_varname = agp,
    # Please write WRA, PSC, Other, or Manual. 
population = Psc, 
# leave crp_ref_value_manual empty, BRINDA R package will use an external crp  reference value for PSC
    crp_ref_value_manual = , 
    # leave agp_ref_value_manual empty, BRINDA R package will use an external agp reference value for PSC
    agp_ref_value_manual = , 
    # Please write FULL or SIMPLE (default output is SIMPLE, if users leave output_format empty)
output_format = full)
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
#> **** AGP: NA
#> **** CRP Variable Name: crp (n = 74)
#> **** Population Group: PSC
#> **** Output Format: FULL
#> -------------------------------------------
#> ** Generate deciles of AGP/CRP based on inputs **
#> **** log-CRP = -2.26
#> -------------------------------------------
#> ** Proceed to the BRINDA adjustment **
#> **** Adjusting Retinol Binding Protein using CRP only
#> **** Adjusting Retinol using CRP only
#> **** Adjusting Ferritin using CRP only
#> **** Adjusted Soluble Transferrin Receptor values are equal to unadjusted Soluble Transferrin Receptor values
#> ****** BRINDA only uses AGP to adjust soluble transferrin receptor
#> ****** You did not provide information on AGP
#> **** Adjusted zinc values are equal to unadjusted zinc values
#> ****** No or weak correlation between Serum Zinc and CRP based on Spearman correlation measures
#> ****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and CRP
#> ** BRINDA adjustment completed **
#> -------------------------------------------
#> ** Proceed to output dataset **
#> variables log_crp_ref, log_crp, log_crp_diff, zn_crp_cor, zn_crp_P_value, rbp_adj, rbp_beta2, rbp_beta2_se, rbp_beta2_P_value, sr_adj, sr_beta2, sr_beta2_se, sr_beta2_P_value, sf_adj, sf_beta2, sf_beta2_se, sf_beta2_P_value, stfr_adj, zn_adj are generated by the BRINDA function
#> -------------------------------------------
#> ** BRINDA adjustment function complete **
#> -------------------------------------------

#
# Example 2: Women of reproductive age (WRA) -----------------------------
# calculate BRINDA inflammation adjustment values for non-pregnant women of reproductive age assuming the sample dataset contain information of women of reproductive age). 
final_data_wra <- BRINDA(dataset = sample_data,
    retinol_binding_protein_varname = rbp,
    retinol_varname = sr,
    ferritin_varname = sf,
    soluble_transferrin_receptor_varname = stfr, 
    zinc_varname = zinc,
    crp_varname = crp, 
    agp_varname = agp,
    # Please write WRA, PSC, Other, or Manual.
    population = WRA, 
# leave crp_ref_value_manual empty, BRINDA R package will use an external crp reference value for WRA
    crp_ref_value_manual = , 
   # leave agp_ref_value_manual empty, BRINDA R package will use an external agp reference value for WRA
    agp_ref_value_manual = , 
    output_format = simple)
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
#> **** Population Group: WRA
#> **** Output Format: SIMPLE
#> -------------------------------------------
#> ** Generate deciles of AGP/CRP based on inputs **
#> **** log-AGP = -0.63
#> **** log-CRP = -1.83
#> -------------------------------------------
#> ** Proceed to the BRINDA adjustment **
#> **** Adjusted Retinol Binding Protein values are equal to unadjusted Retinol Binding Protein values
#> ****** BRINDA does not adjust Retinol Binding Protein among women of reproductive age.
#> **** Adjusted retinol values are equal to unadjusted retinol values
#> ****** BRINDA does not adjust Serum Retinol among women of reproductive age
#> **** Adjusting Ferritin using both AGP and CRP
#> **** Adjusting Soluble Transferrin Receptor using AGP only
#> **** Adjusted zinc values are equal to unadjusted zinc values
#> ****** BRINDA does not adjust Serum Zinc among women of reproductive age
#> ** BRINDA adjustment completed **
#> -------------------------------------------
#> ** Proceed to output dataset **
#> variables rbp_adj, sr_adj, sf_adj, stfr_adj, zn_adj are generated by the BRINDA function
#> -------------------------------------------
#> ** BRINDA adjustment function complete **
#> -------------------------------------------

#
# Example 3: Other population groups ----------------------------------
# calculate BRINDA inflammation adjustment values for other population assuming the study population is neither women of reproductive age nor preschool-age children
#
final_data_other <- BRINDA(dataset = sample_data,
    retinol_binding_protein_varname = rbp,
    retinol_varname = sr,
    ferritin_varname = sf,
    soluble_transferrin_receptor_varname = stfr, 
    zinc_varname = zinc,
    crp_varname = crp, 
    agp_varname = agp,
    population = OTHER, 
# leave crp_ref_value_manual empty, BRINDA R package will calculate the lowest decile of CRP and use it as the crp reference value
    crp_ref_value_manual = , 
# leave agp_ref_value_manual empty, BRINDA R package will calculate the lowest decile of AGP and use it as the agp reference value
   agp_ref_value_manual = , 
   output_format = FULL)
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
#> **** Population Group: OTHER
#> **** Output Format: FULL
#> -------------------------------------------
#> ** Generate deciles of AGP/CRP based on inputs **
#> **** log-AGP = -0.7
#> **** log-CRP = -1.77
#> -------------------------------------------
#> ** Proceed to the BRINDA adjustment **
#> **** Adjusting Retinol Binding Protein using both AGP and CRP
#> **** Adjusting Retinol using both AGP and CRP
#> **** Adjusting Ferritin using both AGP and CRP
#> **** Adjusting Soluble Transferrin Receptor using both AGP and CRP
#> **** Adjusted zinc values are equal to unadjusted zinc values
#> ****** No or weak correlation between Serum Zinc and AGP based on Spearman correlation measures
#> ****** No or weak correlation between Serum Zinc and CRP based on Spearman correlation measures
#> ****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and AGP/CRP
#> ** BRINDA adjustment completed **
#> -------------------------------------------
#> ** Proceed to output dataset **
#> variables log_agp_ref, log_crp_ref, log_agp, log_agp_diff, log_crp, log_crp_diff, zn_agp_cor, zn_agp_P_value, zn_crp_cor, zn_crp_P_value, rbp_adj, rbp_beta1, rbp_beta1_se, rbp_beta1_P_value, rbp_beta2, rbp_beta2_se, rbp_beta2_P_value, sr_adj, sr_beta1, sr_beta1_se, sr_beta1_P_value, sr_beta2, sr_beta2_se, sr_beta2_P_value, sf_adj, sf_beta1, sf_beta1_se, sf_beta1_P_value, sf_beta2, sf_beta2_se, sf_beta2_P_value, stfr_adj, stfr_beta1, stfr_beta1_se, stfr_beta1_P_value, stfr_beta2, stfr_beta2_se, stfr_beta2_P_value, zn_adj are generated by the BRINDA function
#> -------------------------------------------
#> ** BRINDA adjustment function complete **
#> -------------------------------------------

#
# Example 4: User-defined CRP and AGP -------------------------------
# calculate BRINDA inflammation adjustment values for a population when users would like to apply user-defined CRP and AGP reference values 
final_data_manual <- BRINDA(dataset = sample_data,
    retinol_binding_protein_varname = rbp,
    retinol_varname = sr,
    ferritin_varname = sf,
    soluble_transferrin_receptor_varname = stfr, 
    zinc_varname = zinc,
    crp_varname = crp, 
    agp_varname = agp,
    # If users select Manual as the population group, 
    # users can select their own AGP and CRP reference values for the BRINDA adjustment. 
    population = MANUAL, 
    # Enter a user-specified CRP reference value
    crp_ref_value_manual = 0.2, 
    # Enter a user-specified AGP reference value 
    agp_ref_value_manual = 1.4, 
    output_format = SIMPLE)
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
#> **** Population Group: MANUAL
#> **** Output Format: SIMPLE
#> -------------------------------------------
#> ** Generate deciles of AGP/CRP based on inputs **
#> **** log-AGP = 0.34
#> **** log-CRP = -1.61
#> -------------------------------------------
#> ** Proceed to the BRINDA adjustment **
#> **** Adjusting Retinol Binding Protein using both AGP and CRP
#> **** Adjusting Retinol using both AGP and CRP
#> **** Adjusting Ferritin using both AGP and CRP
#> **** Adjusting Soluble Transferrin Receptor using both AGP and CRP
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

Luo, H.; Addo, Y.; Geng, J (2022) BRINDA: Computation of BRINDA Adjusted
Micronutrient Biomarkers for Inflammation. R package version 0.1.5,
<https://CRAN.R-project.org/package=BRINDA>

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
