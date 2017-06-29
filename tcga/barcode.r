`%>%` = magrittr::`%>%`
.codes = import('./code_tables')

#' Regular expression to match TCGA barcodes
#'
#' eg. TCGA-02-0001-01
#'  TCGA - clinical centre - participant id - sample code
#' eg. TCGA-02-0001-01C-01D-0182-01
#'  (same as above)[D/R]NA,.. - 
regex_barcode = paste("^TCGA",                # TCGA identifer
                      "([a-zA-Z0-9]{2})",     # tissue source site (eg. GBM from MBA)
                      "([a-zA-Z0-9]{4})",     # participant id (4 digit alphanumeric)
                      "?([0-9]{2})?([A-Z])?", # tumor/normal id, number of vial
                      "?([0-9]{2})?([A-Z])?", # portion (numbered); analyte (eg. [D/R]NA)
                      "?([a-zA-Z0-9]{4})?",   # plate id (4 digit alphanumeric)
                      "?([0-9]+)?$",          # centre (eg. 01=BROAD GCC)
                      sep = "-")

#' Returns a logical whether the argument is a TCGA barcode
is_barcode = function(x, len=0) {
    grepl(regex_barcode, x) & nchar(x) >= len
}

#' Makes sure the vector argument is all barcode
must_barcode = function(x, len=0) {
    ibc = is_barcode(x, len)
    if (!all(ibc))
        stop("Is not a ", len, "+ char barcode: ", paste(x[!ibc], collapse=", "))
}

#' Regular expression to match TCGA unique IDs
#'
#' 32-bit random string
#' eg. ebf3e73f-41a0-4ca5-b608-fe1c629e16de
regex_id = paste("[a-zA-Z0-9]{8}",
                 "[a-zA-Z0-9]{4}",
                 "[a-zA-Z0-9]{4}",
                 "[a-zA-Z0-9]{4}",
                 "[a-zA-Z0-9]{12}",
                 sep = "-")

#' Take a barcode and convert to data.frame w/ actual info
barcode2index = function(ids) {
    m = stringr::str_match(ids, regex_barcode)

    dplyr::data_frame(Bio.ID=m[,1],
                      TSS.Code = m[,2],
                      Participant.ID = m[,3],
                      Code = m[,4],
                      Vial = m[,5],
                      Portion = m[,6],
                      Analyte = m[,7],
                      Plate.ID = m[,8],
                      Analysis.Center = m[,9]) %>%
        dplyr::left_join(.codes$tissue_source_site, by="TSS.Code") %>%
        dplyr::left_join(.codes$sample_type, by="Code") %>%
        dplyr::left_join(.codes$disease_study, by="Study.Name") %>%
        dplyr::rename(Sample.Code = Code,
                      Sample.Definition = Definition)
}

#' Take a barcode and extract the study from it
#'
#' @param ids  Character vector of TCGA barcodes
#' @return     Named (barcodes) character vector of TCGA study cohorts
barcode2study = function(ids) {
    barcode2index(ids) %>%
        dplyr::pull(Study.Abbreviation) %>%
        setNames(ids)
}

#' Take a barcode and extract the study from it
#'
#' @param ids         Character vector of TCGA barcodes
#' @param short       Use short type descriptors
#' @param factor_ref  Convert to factor with this reference level
barcode2type = function(ids, short=FALSE, factor_ref=NA) {
    if (short)
        ref = "Sample.Code"
    else
        ref = "Sample.Definition"

    re = barcode2index(ids) %>%
        dplyr::select_(ref) %>%
        unlist() %>%
        setNames(ids)

    if (!is.na(factor_ref))
        relevel(factor(re), factor_ref)
    else
        re
}

if (is.null(module_name())) {
    library(testthat)
    expect_true(is_barcode("TCGA-OR-A5J1-01A-11D-A29I-10"))
    expect_true(is_barcode("TCGA-OR-A5J1-01A-11D-A29I"))
    expect_true(is_barcode("TCGA-OR-A5J1-01A-11D"))
    expect_true(is_barcode("TCGA-OR-A5J1-01A"))
    expect_true(is_barcode("TCGA-OR-A5J1-01"))

    expect_false(is_barcode("TCGA-OR-A5J1-01A-11D-A29I-10A"))
    expect_false(is_barcode("TCGA-OR-A5J1-01A-11D-A29IA"))
    expect_false(is_barcode("TCGA-OR-A5J1-01", 16))
}
