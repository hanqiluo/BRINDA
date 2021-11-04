#' Micronutrient biomarker dataset
#'
#' A biomarker data set that was subset from a cross-sectional survey in
#' Malawi. It provides de-identified information for serum ferritin, soluble
#' transferrin receptor, retinol binding protein, retinol, zinc, C-reactive
#' protein (CRP), Alpha1-acid glycoprotein (AGP) to illustrate the use of the
#' package.
#'
#' @docType data
#'
#' @usage data(sample_data)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{id}{Unique identification numbers}
#'  \item{sf}{Serum ferritin, µg/l}
#'  \item{stfr}{Soluble transferrin receptor, mg/L}
#'  \item{rbp}{Retinol binding protein, µmol/L}
#'  \item{sr}{Serum retinol, µmol/L}
#'  \item{zinc}{Zinc, µg/L}
#'  \item{crp}{C-reactive Protein, mg/L}
#'  \item{agp}{Alpha1-acid glycoprotein, g/L}
#' }
#' @references National Statistical Office (NSO), Community Health Sciences Unit
#'  (CHSU) [Malawi], Centers for Disease Control and Prevention (CDC), and Emory
#'  University. 2017. Malawi Micronutrient Survey 2015-16. Atlanta, GA, USA:
#'  NSO, CHSU, CDC, and Emory University
#' @keywords datasets
"sample_data"
