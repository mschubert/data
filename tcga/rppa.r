`%>%` = magrittr::`%>%`
.io = import('ebits/io')
.map_id = import('./map_id')$map_id

#' Get a matrix for all RPPA measurements
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix with antibodies x TCGA samples
rppa = function(tissue, id_type="specimen") {
    fpath = module_file("TCGAbiolinks-downloader/rppa")
    rppa = .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData"))) %>%
        .map_id(id_type=id_type, along="Sample_ID") %>%
        dplyr::rename(Sample = Sample_ID)
}
