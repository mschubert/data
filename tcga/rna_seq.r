.bc = import('./barcode')
.map_id = import('./map_id')$map_id
.coh = import('./cohorts')

#' Get a matrix for all RNA-seq measurements
#'
#' docs: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
#'
#' @param tissue   The tissue(s) to get expression for, or "pan" for all cohorts
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param trans    Value transformation: 'raw', 'log2cpm', 'vst' (default: raw read counts)
#' @param annot    Return SummarizedExperiment annotations or matrix w/ 'external_gene_name'
#' @param excl_ffpe  Exclude FFPE (paraffin embedded) samples
#' @param drop     Drop list if only one tissue requested
#' @return         A matrix with HGNC symbols x TCGA samples
rna_seq = function(tissue, id_type="specimen", trans="raw", annot=FALSE, excl_ffpe=TRUE, drop=TRUE) {
    proc_one = function(expr) {
        if (excl_ffpe)
            expr = expr[,!(is.na(expr$is_ffpe) | expr$is_ffpe)]
        if (identical(annot, TRUE)) {
            colnames(expr) = .map_id(colnames(expr), id_type=id_type)
            expr
        } else {
            re = .map_id(SummarizedExperiment::assay(expr), id_type=id_type)
            if (is.character(annot))
                rownames(re) = SummarizedExperiment::rowData(expr)[[annot]]
            re
        }
    }

    if (tissue == "pan")
        tissue = .coh$cohorts()

    fpath = module_file(paste0("TCGAbiolinks-downloader/rna_seq_", trans))
    expr = file.path(fpath, paste0("TCGA-", tissue, ".rds")) %>%
        lapply(readRDS) %>%
        lapply(proc_one)

    if (drop && length(expr) == 1)
        expr[[1]]
    else if (!annot)
        narray::stack(expr, along=2, allow_overwrite=TRUE)
    else
        expr
}
