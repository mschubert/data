`%>%` = magrittr::`%>%`
.io = import('io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id

#' Get copy number segments
#'
#' docs: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame or GRanges object (Segment_Mean is log2(copy-number)-1)
cna_segments = function(tissue, id_type="specimen", granges=FALSE) {
    fpath = module_file("TCGAbiolinks-downloader/cnv_segments")
    cna = .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData"))) %>%
        .map_id(id_type=id_type, along="Sample") %>%
        dplyr::mutate(ploidy = 1 + 2^Segment_Mean)

    if (granges)
        cna = GenomicRanges::makeGRangesFromDataFrame(cna,
            start.field="Start", end.field="End", keep.extra.columns=TRUE)

    cna
}

#' Get copy number for genes
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix of gene copy number variations
cna_genes = function(tissue, id_type="specimen") {
    stop("not implemented")
}
