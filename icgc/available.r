b = import('ebits/base')
idx = import('./indexing')

#' Finds identifiers where data is available
#'
#' @param rna_seq   RNA sequencing data [T/F]
#' @param rppa      RPPA protein data [T/F]
#' @param clinical  Clinical data [T/F]
#' @param map_to    Which identifers to map to ('specimen' or 'donor')
#' @return          List of identifiers where all data types specified are available
available = function(clinical=NULL, rna_seq=NULL, rppa=NULL, mutations=NULL, map_to="icgc_specimen_id") {
    valid = list()
    if (grepl("specimen", map_to))
        to = "icgc_specimen_id"
    else if (grepl("donor", map_to))
        to = "icgc_donor_id"
    else
        stop("invalid map_to, need to be 'specimen' or 'donor'")

    if (!is.null(clinical))
        valid$clinical = clinical()[[to]]
    if (!is.null(rna_seq))
        valid$rna_seq = b$match(idx$getNames('expr_seq_norm')[[1]],
                                 from = "icgc_sample_id",
                                 to = to,
                                 data = clinical_sample())
    if (!is.null(rppa))
        valid$rppa = b$match(idx$getNames('protein')[[1]],
                              from = "icgc_sample_id",
                              to = to,
                              data = clinical_sample())

    if (!is.null(mutations))
        valid$mut = b$match(idx$getNames('mutations')[[1]],
                             from = "icgc_sample_id",
                             to = to,
                             data = clinical_sample())

    do.call(.b$intersect, valid)
}
