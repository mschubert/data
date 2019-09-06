.io = import('io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id

#' Get a matrix for all CpG methylation beta values
#'
#' docs: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param annot    Return SummarizedExperiment annotations or matrix w/ 'external_gene_name'
#' @return         A matrix with CpG IDs symbols x TCGA samples
meth_cpg = function(tissue, id_type="specimen", annot=FALSE) {
    fpath = module_file("TCGAbiolinks-downloader/meth_beta")
    meth = .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData")))

    if (identical(annot, TRUE))
        meth
    else {
        re = .map_id(SummarizedExperiment::assay(meth), id_type=id_type)
        if (is.character(annot))
            rownames(re) = SummarizedExperiment::rowData(meth)[[annot]]
        re
    }
}
