#' @name BRINDA
#' @title Computation of BRINDA Adjusted Micronutrient Biomarkers for inflammation
#' @author Hanqi Luo, O.Yaw Addo, Afrin Jahan
#'
#' @description Inflammation can affect many micronutrient biomarkers and can thus lead to incorrect diagnosis of individuals and to over- or under-estimate the prevalence of deficiency in a population. Biomarkers Reflecting Inflammation and Nutritional Determinants of Anemia (BRINDA) is a multi-agency and multi-country partnership designed to improve the interpretation of nutrient biomarkers in settings of inflammation and to generate context-specific estimates of risk factors for anemia (Suchdev (2016) <doi:10.3945/an.115.010215>). In the past few years, BRINDA published a series of papers to provide guidance on how to adjust micronutrient biomarkers, retinol binding protein, serum retinol, serum ferritin by Namaste (2020), soluble transferrin receptor (sTfR), serum zinc, serum and Red Blood Cell (RBC) folate, and serum B-12, using inflammation markers, alpha-1-acid glycoprotein (AGP) and/or C-Reactive Protein (CRP) by Namaste (2020) <doi:10.1093/ajcn/nqaa141>, Rohner (2017) <doi:10.3945/ajcn.116.142232>, McDonald (2020) <doi:10.1093/ajcn/nqz304>, and Young (2020) <doi:10.1093/ajcn/nqz303>. The BRINDA inflammation adjustment method mainly focuses on Women of Reproductive Age (WRA) and Preschool-age Children (PSC); however, the general principle of the BRINDA method might apply to other population groups. The BRINDA R package is a user-friendly all-in-one R package that uses a series of functions to implement BRINDA adjustment method, as described above. The BRINDA R package will first carry out rigorous checks and provides users guidance to correct data or input errors (if they occur) prior to inflammation adjustments. After no errors are detected, the package implements the BRINDA inflammation adjustment for up to five micronutrient biomarkers, namely retinol-binding-protein, serum retinol, serum ferritin, sTfR, and serum zinc (when appropriate), using inflammation indicators of AGP and/or CRP for various population groups. Of note, adjustment for serum and RBC folate and serum B-12 is not included in the R package, since evidence shows that no adjustment is needed for these micronutrient biomarkers in either WRA or PSC groups (Young (2020) <doi:10.1093/ajcn/nqz303>).
#'
#' @param dataset Enter the name of the dataset (should be already loaded in the
#' R environment; micronutrient biomarkers should NOT be log-transformed).
#'
#' @param retinol_binding_protein_varname Enter the variable name of retinol
#' binding protein (if available) in your dataset. The variable can be in either
#' international or conventional units. The adjusted values in the output
#' dataset will be in the same unit as retinol binding protein variable in the
#' input dataset.
#'
#' @param retinol_varname Enter the variable name of serum/plasma retinol
#' (if available) in your dataset. The variable can be in either international
#' or conventional units. The adjusted values in the output dataset will be in
#' the same unit as retinol variable in the input dataset.
#'
#' @param ferritin_varname Enter the variable name of serum/plasma ferritin
#' (if available) in your dataset. The variable can be in either international
#' or conventional units. The adjusted values in the output dataset will be in
#' the same unit as serum ferritin in the input dataset.
#'
#' @param soluble_transferrin_receptor_varname Enter the variable name of
#' serum/plasma soluble transferrin receptor (if available) in your dataset.
#' The variable can be in either international or conventional units.
#' The adjusted values in the output dataset will be in the same unit as soluble
#' transferrin receptor in the input dataset.
#'
#' @param zinc_varname Enter the variable name of serum/plasma zinc
#' (if available) in your dataset. The variable can be in either international
#' or conventional units. The adjusted values in the output dataset will be in
#' the same unit as serum zinc in the input dataset.
#'
#' @param crp_varname Enter the variable name of CRP (if available) in your
#' dataset. Unit must be mg/L
#'
#' @param agp_varname Enter the variable name of AGP (if available)  in your
#' dataset (unit must be g/L)
#'
#' @param population_group Please write WRA, PSC, Other, or Manual. The BRINDA R
#' package can only analyze one population group at one time.
#'	If users select WRA or PSC, external CRP/AGP reference values will be used
#'	If users select Other, the lowest decile of the CRP and AGP will be
#'	calculated and used as CRP and AGP reference values
#'	If users select Manual as the population group, users can define their own
#'	AGP and CRP reference values for the BRINDA adjustment.
#'
#' @param crp_ref_value_manual Leave it empty if users select population_group
#' as WRA, PSC, or Other. If users select population_group as Manual,
#' and there is a CRP variable in the dataset, enter a user-specified CRP
#' reference value
#'
#' @param agp_ref_value_manual Leave it empty if users select population_group
#' as WRA, PSC, or Other If users select population_group as Manual, and there
#' is an AGP variable in the dataset, enter a user specified AGP reference value
#'
#' @param output_format Please write FULL or SIMPLE (SIMPLE by default if users
#' leave it empty). The SIMPLE output only provides users adjusted micronutrient
#' biomarker values. The FULL output provides users all the intermediate
#' parameters for the BRINDA adjustment, such as coefficients of log (AGP) and
#' log (CRP) and associated standard errors and P values, in addition to adjusted
#' micronutrient biomarker values.
#'
#' @return `brinda()` returns a data frame object that contains additional
#' variables of adjusted micronutrient biomarkers (by default). If users specify
#' output format = full, the output dataset will also include additional
#' variables such as coefficients of regressions of micronutrient biomarkers on
#' AGP and CRP, natural logs of AGP/CRP reference values.
#'
#' @examples
#' data(sample_data)
#'
#' # Example 1
#' # Calculate BRINDA inflammation adjustment values for preschool-age children
#' # (Assuming the data set contains information of preschool-age children)
#'
#' sample_data_adj <-
#'     BRINDA(dataset = sample_data,
#'        retinol_binding_protein_varname = rbp,
#'        retinol_varname = sr,
#'        ferritin_varname = sf,
#'        soluble_transferrin_receptor_varname = stfr,
#'        zinc_varname = zinc,
#'        crp_varname = crp,
#'        agp_varname = agp,
#'        population = Psc,
#'        crp_ref_value_manual = ,
#'        agp_ref_value_manual = ,
#'        output_format = )
#'
#' # Example 2
#' # Calculate BRINDA inflammation adjustment values for non-pregnant women of
#' # reproductive age assuming the sample data set contains information of women
#' # of reproductive age).
#'
#' sample_data_adj2 <-
#'     BRINDA(dataset = sample_data,
#'        retinol_binding_protein_varname = rbp,
#'        retinol_varname = sr,
#'        ferritin_varname = sf,
#'        soluble_transferrin_receptor_varname = stfr,
#'        zinc_varname = zinc,
#'        crp_varname = ,
#'        agp_varname = agp, population = WRA,
#'        crp_ref_value_manual = ,
#'        agp_ref_value_manual = ,
#'        output_format = )
#'
#' # Example 3
#' # Calculate BRINDA inflammation adjustment values for other population assuming
#' # the study population is neither women of reproductive age nor preschool-age
#' # children
#'
#' sample_data_adj3 <-
#'     BRINDA(dataset = sample_data,
#'        retinol_binding_protein_varname = rbp,
#'        retinol_varname = sr,
#'        ferritin_varname = sf,
#'        soluble_transferrin_receptor_varname = stfr,
#'        zinc_varname = zinc,
#'        crp_varname = crp,
#'        agp_varname = ,
#'        population = OTHER,
#'        crp_ref_value_manual = ,
#'        agp_ref_value_manual = ,
#'        output_format = FULL)
#'
#' # Example 4
#' # Calculate BRINDA inflammation adjustment values for a population when users
#' # would like to apply user-defined CRP and AGP reference values
#'
#' sample_data_adj4 <-
#'     BRINDA(dataset = sample_data,
#'        retinol_binding_protein_varname = rbp,
#'        retinol_varname = sr,
#'        ferritin_varname = sf,
#'        soluble_transferrin_receptor_varname = stfr,
#'        zinc_varname = zinc,
#'        crp_varname = crp,
#'        agp_varname = agp,
#'        population = MANUAL,
#'        crp_ref_value_manual = 0.2,
#'        agp_ref_value_manual = 1.4,
#'        output_format = FULL)
#'
#' @export BRINDA
#' @importFrom data.table "setnames"
#' @importFrom Hmisc "upData"
#' @importFrom berryFunctions "is.error"
#' @importFrom dplyr "%>%" "select" "mutate" "all_of"
#' @import stats
#' @importFrom rlang "quo_name" "enquo"
#' @importFrom rlang ".data"


#
# MAIN BRINDA function ---------------------------------------------------------
#
BRINDA <- function(
    dataset,
    retinol_binding_protein_varname,
    retinol_varname,
    ferritin_varname,
    soluble_transferrin_receptor_varname,
    zinc_varname,
    crp_varname,
    agp_varname,
    population_group,
    crp_ref_value_manual = NULL,
    agp_ref_value_manual = NULL,
    output_format
){

    # dataset: the name of the dataset
    # retinol_binding_protein: variable name of retinol_binding_protein
    # retinol: variable name of serum/plasma retinol
    # ferritin: variable name of ferritin
    # soluble_transferrin_receptor: variable name of soluble transferrin receptor
    # zinc: variable name of serum zinc
    # crp:  variable name of CRP (unit must be mg/L)
    # agp:  variable name of AGP (unit must be g/L)
    # population_group: Can be "WRA", "PSC", "Other", or "Manual"
    # crp_ref_value_manual: User specified CRP reference value if the
    # population_group is Manual
    # agp_ref_value_manual: User specified AGP reference value if the
    # population_group is Manual
    # output: can be "full" or "simple"

    #
    # Checking
    #

    #
    # check if the main dataset exists
    #

    dataset_name <- ifelse(quo_name(enquo(dataset)) == "", NA,
                           quo_name(enquo(dataset)))

    if(is.na(dataset_name)){
        stop("-------------------------------------------
                    You must specify an input dataset
             -------------------------------------------")
    }

    if(!(exists("dataset") && is.data.frame(get("dataset"))))
    {
        stop("-------------------------------------------
                    the input dataset does not exist
             -------------------------------------------")
    }


    #
    # unpack all the values to the quo values
    # for values of space, change into NA
    #
    rbp_quo <- ifelse(quo_name(enquo(retinol_binding_protein_varname)) == "", NA,
                     quo_name(enquo(retinol_binding_protein_varname)))
    sr_quo  <- ifelse(quo_name(enquo(retinol_varname)) == "", NA,
                    quo_name(enquo(retinol_varname)))
    sf_quo  <-  ifelse(quo_name(enquo(ferritin_varname)) == "", NA,
                    quo_name(enquo(ferritin_varname)))
    stfr_quo <- ifelse(quo_name(enquo(soluble_transferrin_receptor_varname)) == "", NA,
                      quo_name(enquo(soluble_transferrin_receptor_varname)))
    zn_quo <- ifelse(quo_name(enquo(zinc_varname)) == "", NA,
                    quo_name(enquo(zinc_varname)))
    agp_quo <- ifelse(quo_name(enquo(agp_varname)) == "", NA,
                     quo_name(enquo(agp_varname)))
    crp_quo <- ifelse(quo_name(enquo(crp_varname)) == "", NA,
                     quo_name(enquo(crp_varname)))
    population_quo <- toupper(quo_name(enquo(population_group)))
    output_format_quo <- toupper(quo_name(enquo(output_format)))
    # assign SIMPLE to output_format_quo
    output_format_quo <- ifelse(output_format_quo == "", "SIMPLE", output_format_quo)


    # check the dataset and compile a new dataset
    biodata <- biodata_gen(dataset = dataset,
                           rbp_quo = rbp_quo,
                           sr_quo = sr_quo,
                           sf_quo = sf_quo,
                           stfr_quo = stfr_quo,
                           zn_quo = zn_quo,
                           agp_quo = agp_quo,
                           crp_quo = crp_quo)

    # check the values of the character variable
    check_parameter(dataset = biodata,
                    population_quo = population_quo,
                    output_format_quo = output_format_quo,
                    crp_ref_value_manual = crp_ref_value_manual,
                    agp_ref_value_manual = agp_ref_value_manual
    )

    message("-------------------------------------------")
    message("** Initial data checks completed **")
    message("-------------------------------------------")

    display_summary(dataset_name = dataset_name,
                    dataset = dataset, rbp_quo = rbp_quo,
                    sr_quo = sr_quo, sf_quo = sf_quo,
                    stfr_quo = stfr_quo, zn_quo = zn_quo,
                    agp_quo = agp_quo, crp_quo = crp_quo,
                    population_quo = population_quo,
                    output_format_quo = output_format_quo)
    message("-------------------------------------------")

    #
    # BRINDA adjustment
    #

    # determine the reference values
    message("** Generate deciles of AGP/CRP based on inputs **")
    biodata <- brinda_decile(dataset = biodata,
                             crp_ref_value_manual = crp_ref_value_manual,
                             agp_ref_value_manual = agp_ref_value_manual,
                             population_quo = population_quo,
                             agp_quo = agp_quo,
                             crp_quo = crp_quo)
    message("-------------------------------------------")

    message("** Proceed to the BRINDA adjustment **")

    biodata <- brinda_adjustment (dataset = biodata,
                                  population_quo = population_quo,
                                  rbp_quo = rbp_quo,
                                  sr_quo = sr_quo,
                                  sf_quo = sf_quo,
                                  stfr_quo = stfr_quo,
                                  zn_quo = zn_quo,
                                  agp_quo = agp_quo,
                                  crp_quo = crp_quo)

    message("** BRINDA adjustment completed **")
    message("-------------------------------------------")

    message("** Proceed to output dataset **")

    #
    # Output dataset:
    #
    biodata <- output_dataset(output_format_quo = output_format_quo,
                              ori_dataset = dataset,
                              brinda_dataset = biodata,
                              rbp_quo = rbp_quo,
                              sr_quo = sr_quo,
                              sf_quo = sf_quo,
                              stfr_quo = stfr_quo,
                              zn_quo = zn_quo,
                              agp_quo = agp_quo,
                              crp_quo = crp_quo)
    message("-------------------------------------------")
    message("** BRINDA adjustment function complete **")
    message("-------------------------------------------")

    return(biodata)
}



#
# generate dataset--------------------------------------------------------------
#

# check all the variables are valid
# check if all variables are numeric
#
biodata_gen <- function(dataset, rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo,
                        agp_quo, crp_quo){
    biodata <- NULL

    mn_biomarker_list <- c(rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo)
    mn_biomarker_valid_list <- mn_biomarker_list[!is.na(mn_biomarker_list)]

    inflammation_marker_list <- c(agp_quo, crp_quo)
    inflammation_marker_valid_list <- inflammation_marker_list [!is.na(inflammation_marker_list)]

    quo_list <- c(rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo, agp_quo, crp_quo)
    quo_valid_list <- quo_list [!is.na(quo_list)]

    name_valid_list <- c("rbp", "sr", "sf", "stfr", "zn", "agp", "crp")[!is.na(quo_list)]


    # check if one biomarker and one inflammation marker exist
    if(length(mn_biomarker_valid_list) < 1){
        stop(paste0(
            "----------------------------------------------
            You must input at least one micronutrient biomarker
            ----------------------------------------------"))
    }

    if(length(inflammation_marker_valid_list) < 1){
        stop(paste0(
            "----------------------------------------------
            You must input at least one inflammation biomarker
            ----------------------------------------------"))
    }

    # check if variables exist in the dataset
    # check if variables are numeric
    # check if variables > 0
    for(biomarker_name in quo_valid_list) {
        if(!is.na(biomarker_name)){

            if(!exists(biomarker_name, dataset)) {
                stop(paste0(
                    "----------------------------------------------
                    variable ", biomarker_name, " does not  in the dataset
                    ----------------------------------------------"))
            }

            if(!is.numeric(as.vector(dataset[[biomarker_name]]))) {
                stop(paste0(
                    "----------------------------------------------
                    variable ", biomarker_name, " is not a numeric variable
                    ----------------------------------------------"))
            }

            if(sum((dataset[[biomarker_name]] < 0), na.rm = TRUE) > 0) {
                stop(paste0(
                    "----------------------------------------------
                    variable ", biomarker_name, " has negative values
                    ----------------------------------------------"))
            }
        }
    }

    biodata <- dataset %>%
        dplyr::select(quo_valid_list)

    colnames(biodata) <- name_valid_list

    return(biodata)
}

#
# Check if the input parameters are reasonable----------------------------------
#
check_parameter <- function(population_quo,
                            output_format_quo,
                            crp_ref_value_manual,
                            agp_ref_value_manual,
                            dataset){

    if(!(population_quo %in% c("WRA", "PSC", "OTHER", "MANUAL"))) {
        stop("----------------------------------------------
             population_group can only be WRA, PSC, OTHER, or MANUAL
             ----------------------------------------------")
    }

    if(!(output_format_quo %in% c("SIMPLE", "FULL"))) {
        stop("----------------------------------------------
             output_format can only be SIMPLE or FULL
             ----------------------------------------------")
    }

    if(population_quo == "MANUAL"){
        if(is.null(agp_ref_value_manual) & exists("agp", dataset)){
            stop("----------------------------------------------
                  When you identify the population_group as MANUAL and specify a AGP variable,
                  you must specify a reference value for parameter agp_ref_value_manual
                  -----------------------------------------------")
        }

        if(!is.null(agp_ref_value_manual) & !exists("agp", dataset)){
            stop("----------------------------------------------
                  You should not specify a agp reference value when there is no agp variable
                  -----------------------------------------------")
        }

        if(!is.null(agp_ref_value_manual)){
            if(!is.numeric(agp_ref_value_manual)){
                stop("----------------------------------------------
                     agp_ref_value_manual should be a numeric value
                     -----------------------------------------------")
            }

            if(agp_ref_value_manual <= 0){
                stop("----------------------------------------------
                     agp_ref_value_manual should be a positive value
                     -----------------------------------------------")
            }
        }

        if(is.null(crp_ref_value_manual) & exists("crp", dataset)){
            stop("----------------------------------------------
                  When you identify the population_group as MANUAL and specify a CRP variable,
                  you must specify a reference value for parameter crp_ref_value_manual
                  -----------------------------------------------")
        }

        if(!is.null(crp_ref_value_manual)) {
            if(!is.numeric(crp_ref_value_manual)){
                stop("----------------------------------------------
                     crp_ref_value_manual should be a numeric value
                     ----------------------------------------------")
            }

            if(crp_ref_value_manual <= 0){
                stop("----------------------------------------------
                     crp_ref_value_manual should be a postive value
                     ----------------------------------------------")
            }

            if(!is.null(crp_ref_value_manual) & !exists("crp", dataset)){
                stop("----------------------------------------------
                  You should not specify a CRP reference value when there is no CRP variable
                  -----------------------------------------------")
            }
        }
    }

    if(population_quo != "MANUAL"){
        if (!is.null(crp_ref_value_manual) | !is.null(agp_ref_value_manual)){
            stop("----------------------------------------------
                  You are only allowed to specify a value for crp_ref_value_manual or agp_ref_value_manual
                  if population_group is MANUAL
                     ----------------------------------------------")
        }
    }
}


#
# display summary --------------------------------------------------------------
#
display_summary <-
    function(dataset_name, dataset, rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo,
             agp_quo, crp_quo, population_quo,  output_format_quo){

        quo_list <- c(rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo, agp_quo, crp_quo)
        mn_biomarker_full_variable_name_list <- c("Retinol Binding Protein", "Retinol",
                                                  "Ferritin", "Soluble Transferrin Receptor", "Zinc",
                                                  "AGP", "CRP")

        message("** Overview of the dataset and BRINDA package inputs **")
        message(paste0("**** Dataset Name: ", dataset_name, "**"))

        # biomarkers
        for(number in 1:length(quo_list)){
            if(is.na(quo_list[number])) {
                message(paste0("**** ",  mn_biomarker_full_variable_name_list[number], ": NA"))
            } else {
                message(paste0("**** ",  mn_biomarker_full_variable_name_list[number], " Variable Name: ",
                             quo_list[number], " (n = ", sum(!is.na(dataset %>%
                                                                        dplyr::select(quo_list[number]))), ")"))
            }
        }
        # Population group and output format
        message(paste0("**** Population Group: ", population_quo))
        message(paste0("**** Output Format: ", output_format_quo))
    }

#
# BRINDA adjustment - finding decile--------------------------------------------
#

# WRA: lnAGP = -0.63 ln(g/L)
# WRA: lnCRP = -1.83 ln(mg/L)
# PSC: lnAGP = -0.52 ln(g/L)
# PSC: lnCRP = -2.26 ln(mg/L)

brinda_decile <- function(dataset, population_quo, crp_ref_value_manual,
                          agp_ref_value_manual, agp_quo, crp_quo) {
    if(population_quo == "WRA"){
        log_crp_ref <- -1.83
        log_agp_ref <- -0.63
    }

    if(population_quo == "PSC"){
        log_crp_ref <- -2.26
        log_agp_ref <- -0.52
    }

    if(population_quo == "OTHER"){
        log_crp_ref <- ifelse(exists("crp", dataset), log(quantile(dataset$crp, 0.1, na.rm = T)), NA)
        log_agp_ref <- ifelse(exists("agp", dataset), log(quantile(dataset$agp, 0.1, na.rm = T)), NA)
    }

    if(population_quo == "MANUAL"){
        log_crp_ref <- ifelse(is.null(crp_ref_value_manual), NA, log(crp_ref_value_manual))
        log_agp_ref <- ifelse(is.null(agp_ref_value_manual), NA, log(agp_ref_value_manual))
    }

    # save the variable in the new dataset
    # AGP
    if(!is.na(agp_quo)){
        dataset <- dataset %>%
            mutate(log_agp_ref = log_agp_ref)

        var.labels <- c(
            log_agp_ref   = "BRINDA output variable: natural log of the AGP reference value")

        dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)

        message(paste0("**** log-AGP = ", round(unique(dataset$log_agp_ref), 2)))
    }

    if(!is.na(crp_quo)) {
        dataset <- dataset %>%
            mutate(log_crp_ref = log_crp_ref)

        var.labels <- c(
            log_crp_ref   = "BRINDA output variable: natural log of the CRP reference value")

        dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)

        message(paste0("**** log-CRP = ", round(unique(dataset$log_crp_ref), 2)))
    }




    return(dataset)
}



#
# BRINDA adjustment - core -----------------------------------------------------
#
brinda_adjustment <- function(dataset, rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo,
                              agp_quo, crp_quo, population_quo) {


    mn_biomarker_list <- c(rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo)
    mn_biomarker_variable_name_list <- c("rbp", "sr", "sf", "stfr", "zn")[!is.na(mn_biomarker_list)]
    mn_biomarker_full_variable_name_list <- c("Retinol Binding Protein", "Retinol",
                                              "Ferritin", "Soluble Transferrin Receptor", "Zinc")[!is.na(mn_biomarker_list)]

    # CRP/AGP values
    if(exists("agp", dataset)){
        dataset <-
            dataset %>%
            mutate(
                log_agp      = log(ifelse(.data$agp == 0, .data$agp + 0.001, .data$agp)),
                log_agp_diff = pmax(.data$log_agp - .data$log_agp_ref, 0))

        # AGP/CRP labels
        var.labels <- c(log_agp = "BRINDA output variable: natural log of the AGP",
                        log_agp_diff = "BRINDA output variable: difference between natural log of the AGP and the AGP reference value")

        dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
    }

    if(exists("crp", dataset)){
        dataset <-
            dataset %>%
            mutate(
                log_crp      = log(ifelse(.data$crp == 0, .data$crp + 0.001, .data$crp)),
                log_crp_diff = pmax(.data$log_crp - .data$log_crp_ref, 0))

        # AGP/CRP labels
        var.labels <- c(log_crp = "BRINDA output variable: natural log of the CRP",
                        log_crp_diff = "BRINDA output variable: difference between natural log of the CRP and the CRP reference value")

        dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
    }

    # add zinc parameters if the CRP/AGP reference values
    dataset <- zinc_spearman_corr(dataset = dataset,
                                  zn_quo = zn_quo,
                                  crp_quo = crp_quo,
                                  agp_quo = agp_quo,
                                  population_quo = population_quo)

    for(biomarker in mn_biomarker_variable_name_list) {
        #print(biomarker)
        dataset$biomarker <- dataset[[biomarker]]

        dataset <-
            dataset %>%
            mutate(log_biomarker = log(ifelse(.data$biomarker == 0, .data$biomarker + 0.001, .data$biomarker)))

        mn_biomarker_full_name <- mn_biomarker_full_variable_name_list[which(mn_biomarker_variable_name_list  == biomarker)]

        # sTfR and AGP
        if(biomarker == "stfr" & exists("agp", dataset)){
            dataset <- brinda_adjustment_agp(
                dataset = dataset,
                mn_biomarker_full_name = mn_biomarker_full_name,
                biomarker = biomarker)
        }

        # sTfR and no AGP
        if(biomarker == "stfr" & !exists("agp", dataset)){
            dataset$biomarker_adj <- dataset$biomarker
            message("**** Adjusted Soluble Transferrin Receptor values are equal to unadjusted Soluble Transferrin Receptor values")
            message("****** BRINDA only uses AGP to adjust soluble transferrin receptor")
            message("****** You did not provide information on AGP")


            var.labels <- c(biomarker_adj = "BRINDA output variable: adjusted Soluble Transferrin Receptor")

            dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
            setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))
        }

        # Non stfr/zn data
        # Also exclude data with RBP, SR, ZN among WRA group
        if(!(biomarker %in% c("rbp", "sr", "zn") & population_quo == "WRA")){
            if(!(biomarker %in% c("stfr", "zn"))){
                # Both AGP and CRP
                if(exists("agp", dataset) & exists("crp", dataset)){
                    dataset <- brinda_adjustment_agp_crp(
                        dataset = dataset,
                        mn_biomarker_full_name = mn_biomarker_full_name,
                        biomarker = biomarker)
                }

                # Only AGP
                if(exists("agp", dataset) & !exists("crp", dataset)){

                    dataset <- brinda_adjustment_agp(
                        dataset = dataset,
                        mn_biomarker_full_name = mn_biomarker_full_name,
                        biomarker = biomarker)
                }

                # Only CRP
                if(!exists("agp", dataset) & exists("crp", dataset)){

                    dataset <- brinda_adjustment_crp(
                        dataset = dataset,
                        mn_biomarker_full_name = mn_biomarker_full_name,
                        biomarker = biomarker)
                }
            }
            # no stfr/zn done
        }

        # Fix for serum_retinol/RBP/zinc for WRA - No BRINDA adjustment
        if(population_quo == "WRA" & biomarker =="rbp"){
            message("**** Adjusted Retinol Binding Protein values are equal to unadjusted Retinol Binding Protein values")
            message("****** BRINDA does not adjust Retinol Binding Protein among women of reproductive age.")


            dataset$biomarker_adj <- dataset$biomarker

            var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name))
            dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
            setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))

        }

        if(population_quo == "WRA" & biomarker =="sr"){
            message("**** Adjusted retinol values are equal to unadjusted retinol values")
            message("****** BRINDA does not adjust Serum Retinol among women of reproductive age")


            dataset$biomarker_adj <- dataset$biomarker

            var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name))
            dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
            setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))
        }

        if(population_quo == "WRA" & biomarker =="zn"){
            message("**** Adjusted zinc values are equal to unadjusted zinc values")
            message("****** BRINDA does not adjust Serum Zinc among women of reproductive age")

            dataset$biomarker_adj <- dataset$biomarker

            var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name))
            dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
            setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))
        }

        if(biomarker == "zn" & population_quo != "WRA"){
            # when AGP and CRP are both available
            if(exists("agp", dataset) & exists("crp", dataset)){
                if((unique(dataset$zn_agp_cor) < -0.2 & unique(dataset$zn_agp_P_value) < 0.1) |
                   (unique(dataset$zn_crp_cor) < -0.2 & unique(dataset$zn_crp_P_value) < 0.1)){
                    dataset <- brinda_adjustment_agp_crp(
                        dataset = dataset,
                        mn_biomarker_full_name = mn_biomarker_full_name,
                        biomarker = biomarker)
                }else{
                    message("**** Adjusted zinc values are equal to unadjusted zinc values")
                    message("****** No or weak correlation between Serum Zinc and AGP based on Spearman correlation measures")
                    message("****** No or weak correlation between Serum Zinc and CRP based on Spearman correlation measures")
                    message("****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and AGP/CRP")

                    dataset$biomarker_adj <- dataset$biomarker

                    var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name))
                    dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
                    setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))
                }
            }

            # when only AGP available
            if(exists("agp", dataset) & !exists("crp", dataset)){
                if((unique(dataset$zn_agp_cor) < -0.2 & unique(dataset$zn_agp_P_value) < 0.1)){
                    dataset <- brinda_adjustment_agp_crp(
                        dataset = dataset,
                        mn_biomarker_full_name = mn_biomarker_full_name,
                        biomarker = biomarker)
                }else{
                    message("**** Adjusted zinc values are equal to unadjusted zinc values")
                    message("****** No or weak correlation between Serum Zinc and AGP based on Spearman correlation measures")
                    message("****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and AGP")
                    dataset$biomarker_adj <- dataset$biomarker

                    var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name))
                    dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
                    setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))
                }
            }

            # when only CRP
            if(!exists("agp", dataset) & exists("crp", dataset)){
                if((unique(dataset$zn_crp_cor) < -0.2 & unique(dataset$zn_crp_P_value) < 0.1)){
                    dataset <- brinda_adjustment_crp(
                        dataset = dataset,
                        mn_biomarker_full_name = mn_biomarker_full_name,
                        biomarker = biomarker)
                }else{
                    message("**** Adjusted zinc values are equal to unadjusted zinc values")
                    message("****** No or weak correlation between Serum Zinc and CRP based on Spearman correlation measures")
                    message("****** BRINDA does not adjust Serum Zinc because of no or weak correlation between Serum Zinc and CRP")

                    dataset$biomarker_adj <- dataset$biomarker

                    var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name))
                    dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
                    setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))
                }
            }
        } # finish zinc

        dataset$biomarker <- NULL
        dataset$log_biomarker <- NULL

    }
    return(dataset)
}



#
# Produce Zinc statistics ------------------------------------------------------
#
zinc_spearman_corr <- function(dataset, zn_quo, crp_quo, agp_quo, population_quo) {
    if(!is.na(zn_quo) & population_quo != "WRA"){
        if(!is.na(agp_quo)){
            spearman_agp_results <- cor.test(dataset$zn, dataset$agp, method=c("spearman"), exact=F)
            dataset <-
                dataset %>%
                mutate(zn_agp_cor = spearman_agp_results$estimate,
                       zn_agp_P_value = spearman_agp_results$p.value)

            var.labels <- c(zn_agp_cor = "BRINDA output variable: Spearman correlation coefficient between serum zinc and AGP",
                            zn_agp_P_value = "BRINDA output variable: P-value of Spearman correlation between serum zinc and AGP")

            dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
        }

        if(!is.na(crp_quo)){
            spearman_crp_results <- cor.test(dataset$zn, dataset$crp, method=c("spearman"), exact=F)
            dataset <-
                dataset %>%
                mutate(zn_crp_cor = spearman_crp_results$estimate,
                       zn_crp_P_value = spearman_crp_results$p.value)

            var.labels <- c(zn_crp_cor = "BRINDA output variable: Spearman correlation coefficient between serum zinc and CRP",
                            zn_crp_P_value = "BRINDA output variable: P-value of Spearman correlation between serum zinc and CRP")

            dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)
        }
    }
    return(dataset)
}


#
# Both AGP and CRP adjustment---------------------------------------------------
#
brinda_adjustment_agp_crp <- function(dataset, mn_biomarker_full_name, biomarker){
    message(paste0("**** Adjusting ",  mn_biomarker_full_name, " using both AGP and CRP"))

    if(!is.error(lm(log_biomarker ~ log_agp + log_crp, data = dataset, na.action=na.omit))){
        mysvyglm <- lm(log_biomarker ~ log_agp + log_crp, data = dataset, na.action=na.omit)
        beta1 <- coef(mysvyglm)[2]
        beta1_se <- summary(mysvyglm)$coefficients[2, 2]
        beta1_P_value <- summary(mysvyglm)$coefficients[2, 4]

        beta2 <- coef(mysvyglm)[3]
        beta2_se <- summary(mysvyglm)$coefficients[3, 2]
        beta2_P_value <- summary(mysvyglm)$coefficients[3, 4]

    } else {
        beta1 <- 0
        beta1_se <- 0
        beta1_P_value <- 0
        beta2 <- 0
        beta2_se <- 0
        beta2_P_value <- 0
    }
    dataset <- dataset %>%
        mutate(biomarker_adj = exp(.data$log_biomarker - beta1 * .data$log_agp_diff - beta2 * .data$log_crp_diff),
               beta1 = beta1,
               beta1_se = beta1_se,
               beta1_P_value = beta1_P_value,
               beta2 = beta2,
               beta2_se = beta2_se,
               beta2_P_value = beta2_P_value
        )

    var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name),
                   beta1         = paste0("BRINDA output variable: reg coeff of ln(AGP): ln(",  mn_biomarker_full_name, ") against ln(AGP) and ln(CRP)"),
                   beta1_se      = paste0("BRINDA output variable: standard error of reg coeff of ln(AGP): ln(",  mn_biomarker_full_name, ") against ln(AGP) and ln(CRP)"),
                   beta1_P_value = paste0("BRINDA output variable: P-value of reg coeff of ln(AGP): ln(",  mn_biomarker_full_name, ") against ln(AGP) and ln(CRP)"),
                   beta2         = paste0("BRINDA output variable: reg coeff of ln(CRP): ln(",  mn_biomarker_full_name, ") against ln(AGP) and ln(CRP)"),
                   beta2_se      = paste0("BRINDA output variable: standard error of reg coeff of  ln(CRP): ln(",  mn_biomarker_full_name, ") against ln(AGP) and ln(CRP)"),
                   beta2_P_value = paste0("BRINDA output variable: P-value of reg coeff of ln(CRP): ln(",  mn_biomarker_full_name, ") against ln(AGP) and ln(CRP)"))

    dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)

    setnames(dataset, "beta1", paste0(biomarker, "_beta1"))
    setnames(dataset, "beta1_se", paste0(biomarker, "_beta1_se"))
    setnames(dataset, "beta1_P_value", paste0(biomarker, "_beta1_P_value"))
    setnames(dataset, "beta2", paste0(biomarker, "_beta2"))
    setnames(dataset, "beta2_se", paste0(biomarker, "_beta2_se"))
    setnames(dataset, "beta2_P_value", paste0(biomarker, "_beta2_P_value"))
    setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))

    return(dataset)
}

#
# BRINDA adjustment - only AGP -------------------------------------------------
#
brinda_adjustment_agp <- function(dataset, mn_biomarker_full_name, biomarker){

    message(paste0("**** Adjusting ",  mn_biomarker_full_name, " using AGP only"))

    if(!is.error(lm(log_biomarker ~ log_agp, data = dataset, na.action=na.omit))){
        mysvyglm <- lm(log_biomarker ~ log_agp, data = dataset, na.action=na.omit)
        beta1 <- coef(mysvyglm)[2]
        beta1_se <- summary(mysvyglm)$coefficients[2, 2]
        beta1_P_value <- summary(mysvyglm)$coefficients[2, 4]
    } else {
        beta1 <- 0
        beta1_se <- 0
        beta1_P_value <- 0
    }
    dataset <- dataset %>%
        mutate(biomarker_adj = exp(.data$log_biomarker - beta1 * .data$log_agp_diff),
               beta1 = beta1,
               beta1_se = beta1_se,
               beta1_P_value = beta1_P_value)

    var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name),
                   beta1         = paste0("BRINDA output variable: reg coeff: ln(", mn_biomarker_full_name, ") against ln(AGP)"),
                   beta1_se      = paste0("BRINDA output variable: standard error of reg coeff of ln(AGP): ln(",  mn_biomarker_full_name, ") against ln(AGP)"),
                   beta1_P_value = paste0("BRINDA output variable: P-value of reg coeff of ln(AGP): ln(",  mn_biomarker_full_name, ") against ln(AGP)"))

    dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)

    setnames(dataset, "beta1", paste0(biomarker, "_beta1"))
    setnames(dataset, "beta1_se", paste0(biomarker, "_beta1_se"))
    setnames(dataset, "beta1_P_value", paste0(biomarker, "_beta1_P_value"))
    setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))

    return(dataset)
}

#
# BRINDA adjustment - CRP only -------------------------------------------------
#
brinda_adjustment_crp <- function(dataset, mn_biomarker_full_name, biomarker){

    message(paste0("**** Adjusting ",  mn_biomarker_full_name, " using CRP only"))

    if(!is.error(lm(log_biomarker ~ log_crp, data = dataset, na.action=na.omit))){
        mysvyglm <- lm(log_biomarker ~ log_crp, data = dataset, na.action=na.omit)
        beta2 <- coef(mysvyglm)[2]
        beta2_se <- summary(mysvyglm)$coefficients[2, 2]
        beta2_P_value <- summary(mysvyglm)$coefficients[2, 4]
    } else {
        beta2 <- 0
        beta2_se <- 0
        beta2_P_value <- 0
    }
    dataset <- dataset %>%
        mutate(biomarker_adj = exp(.data$log_biomarker  - beta2 * .data$log_crp_diff),
               beta2 = beta2,
               beta2_se = beta2_se,
               beta2_P_value = beta2_P_value)

    var.labels <- c(biomarker_adj = paste0("BRINDA output variable: adjusted ", mn_biomarker_full_name),
                   beta2         = paste0("BRINDA output variable: reg coeff of ln(CRP): ln(",  mn_biomarker_full_name, ") against ln(CRP)"),
                   beta2_se      = paste0("BRINDA output variable: standard error of reg coeff of ln(CRP): ln(",  mn_biomarker_full_name, ") against ln(CRP)"),
                   beta2_P_value = paste0("BRINDA output variable: P-value of reg coeff of ln(CRP): ln(",  mn_biomarker_full_name, ") against ln(CRP)"))

    dataset <- Hmisc::upData(dataset, labels = var.labels, print = F)

    setnames(dataset, "beta2", paste0(biomarker, "_beta2"))
    setnames(dataset, "beta2_se", paste0(biomarker, "_beta2_se"))
    setnames(dataset, "beta2_P_value", paste0(biomarker, "_beta2_P_value"))
    setnames(dataset, "biomarker_adj", paste0(biomarker, "_adj"))

    return(dataset)
}


#
# Output dataset ---------------------------------------------------------------
#
output_dataset <- function(output_format_quo, ori_dataset, brinda_dataset,
                           rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo, agp_quo, crp_quo){

    mn_biomarker_list <- c(rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo)
    mn_inflammation_biomarker_list <- c(rbp_quo, sr_quo, sf_quo, stfr_quo, zn_quo, agp_quo, crp_quo)

    mn_biomarker_variable_name_list <- c("rbp", "sr", "sf", "stfr", "zn")[!is.na(mn_biomarker_list)]
    mn_inflammation_biomarker_variable_name_list <- c("rbp", "sr", "sf", "stfr", "zn", "agp", "crp")[!is.na(mn_inflammation_biomarker_list)]

    mn_biomarker_adj_variable_name_list <- paste0(mn_biomarker_variable_name_list, "_adj")
    mn_biomarker_full_name_list <- names(brinda_dataset %>%
                                             dplyr::select(-all_of(mn_inflammation_biomarker_variable_name_list)))

    if(output_format_quo == "SIMPLE"){
        if(sum(names(ori_dataset) %in% mn_biomarker_adj_variable_name_list) > 0 ) {
            output_dataset <-
                cbind(
                    ori_dataset %>%
                        dplyr::select(-mn_biomarker_adj_variable_name_list[mn_biomarker_adj_variable_name_list %in% names(ori_dataset)]),
                    brinda_dataset %>%
                        dplyr::select(all_of(mn_biomarker_adj_variable_name_list)))
        } else {
            output_dataset <- cbind(ori_dataset,
                                    brinda_dataset %>%
                                        dplyr::select(mn_biomarker_adj_variable_name_list))
        }
        message(paste0("variables ", toString(mn_biomarker_adj_variable_name_list), " are generated by the BRINDA function"))
    }

    if(output_format_quo == "FULL"){
        if(sum(names(ori_dataset) %in% mn_biomarker_full_name_list) > 0 ) {
            output_dataset <-
                cbind(
                    ori_dataset %>%
                        dplyr::select(-mn_biomarker_full_name_list[mn_biomarker_full_name_list %in% names(ori_dataset)]),
                    brinda_dataset %>%
                        dplyr::select(all_of(mn_biomarker_full_name_list)))
        } else{
            output_dataset <- cbind(ori_dataset, brinda_dataset %>%
                                        dplyr::select(all_of(mn_biomarker_full_name_list)))
        }

        message(paste0("variables ", toString(mn_biomarker_full_name_list), " are generated by the BRINDA function"))
    }
    return(output_dataset)
}
# output function complete



