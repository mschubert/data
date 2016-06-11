library(dplyr)
.io = import('ebits/io')
.p = import('../path')
.map_id = import('./map_id')$map_id

#' List all available tissues
tissues = function(id_type="specimen") {
    tt = list.files(.p$path("tcga"), pattern="_voom\\.RData")
    sub("_voom.RData", "", tt)
}

#' Get a data.frame listing all clinical data available
#'
#' @param tissue   Limit selection to a single/set of tissues
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the clinical data
clinical = function(tissue=NULL, id_type=NULL) {
    re = .p$load("tcga", "clinical.RData")
    if (!is.null(tissue))
        re = dplyr::filter(re, study %in% tissue)
    if (!is.null(id_type))
        rownames(re) = substr(re$Tumor_Sample_Barcode, 1,
                              setNames(.id_lengths, id_types)[id_type])
    re
}

#' Get a matrix for all RNA-seq measurements
#'
#' @param tissue   The tissue to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix with HGNC symbols x TCGA samples
rna_seq = function(tissue, id_type="specimen", ...) {
    .p$load("tcga", paste0(tissue, "_voom.RData")) %>%
        .map_id(id_type=id_type, ...)
}

#' Get a matrix for all RPPA measurements
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix with antibodies x TCGA samples
rppa = function(id_type="specimen", ...) {
    .p$load("tcga", "rppa.RData") %>%
        .map_id(id_type=id_type, ...)
}

#' Get a data.frame listing all mutations and types
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
mutations = function(id_type="specimen", ...) {
    # 19k/22k in READ and all 20k OV have no portion in barcode, assuming 'A'
    .p$load("tcga", "mutations.RData") %>%
        mutate(Tumor_Sample_Barcode =
                ifelse(nchar(Tumor_Sample_Barcode == 15),
                       paste0(Tumor_Sample_Barcode, "A"),
                       Tumor_Sample_Barcode)) %>%
        .map_id(id_type=id_type, along="Tumor_Sample_Barcode", ...)
}

#' Get a data.frame listing all GISTIC scores for CNAs
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
cna = function(id_type="specimen", ...) {
    .p$load("tcga", "cna.RData") %>%
        .map_id(id_type=id_type, ...)
}

#' Get a data.frame listing all GISTIC scores for CNAs
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
cna_thresholded = function(id_type="specimen", ...) {
    .p$load("tcga", "cna_thresholded.RData") %>%
        .map_id(id_type=id_type, ...)
}
