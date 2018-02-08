.io = import('ebits/io')

#' Get a data.frame listing all clinical data available
#'
#' @param tissue   Limit selection to a single/set of tissues
#' @return         A data.frame with data for all the clinical data
clinical = function(tissue=NULL, id_type="patient") {
    fpath = module_file("TCGAbiolinks-downloader/clinical")
    .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData")))
}
