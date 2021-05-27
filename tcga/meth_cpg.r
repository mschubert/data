.io = import('io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id
`%>%` = magrittr::`%>%`

#' Get gene annotations for CpGs for HM450
#'
#' @return  A data.frame with CpG annotations
meth_cpg2gene = function() {
    fname = "HM450.hg38.manifest.gencode.v37.tsv.gz"
    fpath = file.path(module_file("cache"), fname)
    if (!file.exists(fpath)) {
        url = paste0("http://zhouserver.research.chop.edu/InfiniumAnnotation/20210414/HM450/", fname)
        download.file(url, fpath)
    }
    readr::read_tsv(fpath)
}

#' Get annotations for CpG islands for HM450
#'
#' @return  A data.frame with CpG annotations
meth_cpg2island = function() {
    fname = "HM450.hg38.manifest.CpGIsland.tsv.gz"
    fpath = file.path(module_file("cache"), fname)
    if (!file.exists(fpath)) {
        url = paste0("http://zhouserver.research.chop.edu/InfiniumAnnotation/20210414/HM450/", fname)
        download.file(url, fpath)
    }
    readr::read_tsv(fpath)
}

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

#' Core promoter element of CpGs +/- 200 bp from TSS
#'
#' @inheritParams meth_cpg
meth_promoter_core = function(tissue, id_type="specimen", gene_type="ensembl_gene_id") {
    stop("not implemented")
}

#' Extended promoter element of CpGs +/- 1500 bp from TSS
#'
#' @inheritParams meth_cpg
meth_promoter_ext = function(tissue, id_type="specimen", gene_type="ensembl_gene_id") {
    stop("not implemented")
}

#' Gene body element of CpGs +1500/TTS
#'
#' @inheritParams meth_cpg
meth_gene_body = function(tissue, id_type="specimen", gene_type="ensembl_gene_id") {
    stop("not implemented")
}

#' GRanges object of islands mapped to genes
cpg_gene = function() {
    res = module_file("data", "hm450.hg38.manifest.gencode.v22.rds", mustWork=TRUE) %>%
        readRDS()
}

#' Tibble of CpGs mapped to TF binding sites
cpg_tf = function() {
    res = module_file("data", "hm450.hg19.ENCODE.TFBS.tsv.gz", mustWork=TRUE) %>%
        readr::read_tsv()
}
