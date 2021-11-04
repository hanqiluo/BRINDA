#' Micronutrient biomarker dataset
#'
#' Small,toy biomarker data set that was subset from a cross-sectional survey in Malawi
#' It provides artificial information for serum ferritin, soluble transferrin receptor,
#' retinol binding protein, retinol, zinc, C-reactive protein, Alpha1-acid glycoprotein to
#' illustrate the use of the package.
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
#' @references This data set was subset from the National Micronutrient Survey in Malawi (2016) for the BRINDA package.
#' @keywords datasets
#' @examples
#'
#' data(sample_data)
#' head(sample_data)
#'
"sample_data"
