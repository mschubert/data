#library(GenomicRanges)
.io = import('io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id
`%>%` = magrittr::`%>%`

#' Get a matrix for all RNA-seq measurements
#'
#' docs: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param trans    Value transformation: 'raw', 'log2cpm', 'vst' (default: raw read counts)
#' @param annot    Return SummarizedExperiment annotations or matrix w/ 'external_gene_name'
#' @return         A matrix with HGNC symbols x TCGA samples
rna_isoform = function(tissue, id_type="specimen", trans="raw", annot=FALSE) {
    fpath = module_file(paste0("TCGAbiolinks-downloader/rna_isoform_", trans))
    expr = .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData")))

    if (identical(annot, TRUE))
        stop("annotation not provided by TCGAbiolinks")

    if (trans == "raw") {
        expr = round(expr[,grepl("^raw_count_", colnames(expr))])
        colnames(expr) = sub("^raw_count_", "", colnames(expr))
    }

    .map_id(expr, id_type=id_type, along=2)
}
