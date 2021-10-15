`%>%` = magrittr::`%>%`
.map_id = import('./map_id')$map_id

#' Get a data.frame listing all mutations and types
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
mutations = function(cohort, id_type="specimen", ...) {
    fdir = module_file("TCGAbiolinks-downloader", "snv_mutect2")
    fname = file.path(fdir, sprintf("TCGA-%s.rds", cohort))
    tibble::as_tibble(readRDS(fname)) %>%
        dplyr::mutate(Sample = .map_id(Tumor_Sample_Barcode, id_type=id_type)) %>%
        select(Sample, everything())
}
