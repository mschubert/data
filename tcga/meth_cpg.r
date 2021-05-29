import_package("plyranges", attach=TRUE)
.io = import('io')
.seq = import('seq')
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

#' Get CpGs for TF binding sites
meth_cpg2tf = function() {
    fname = "hm450.hg19.ENCODE.TFBS.tsv.gz"
    res = module_file("data", fname, mustWork=TRUE) %>%
        readr::read_tsv()
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
        re[rowSums(is.na(re)) == 0,] # only include cgs with full data
    }
}

#' Methylation level summary at core/extended promoter, gene body, or both
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Gene identifier: ensembl_{gene,transcript}_id, external_gene_name
#' @return         A matrix with summarized methylation: (core) promoter, body, etc.
meth_summary = function(tissue, id_type="ensembl_gene_id") {
    params = paste(c("meth_summary", tissue, id_type), collapse="-")
    fpath = file.path(module_file("cache"), paste0(params, ".rds"))
    if (file.exists(fpath))
        return(readRDS(fpath))
    cgmat = meth_cpg(tissue)

    c2g = meth_cpg2gene() %>%
        dplyr::select(seqnames=CpG_chrm, start=CpG_beg, end=CpG_end,
                      strand=probe_strand, probe_id=probeID) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    genes = .seq$gene_table() %>%
        dplyr::select(chromosome_name, strand, start=start_position,
                      end=end_position, ensembl_gene_id, external_gene_name,
                      ensembl_transcript_id) %>%
        dplyr::mutate(chromosome_name = paste0("chr", chromosome_name),
                      strand = setNames(c("+","-"),c(1,-1))[as.character(strand)]) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    .do = . %>%
        join_overlap_intersect(c2g) %>% as.data.frame() %>% dplyr::as_tibble() %>%
        select(probe_id, rlang::sym(id_type)) %>% dplyr::distinct()
    cgs = list(
        gene = genes %>% anchor_3p() %>% stretch(1500) %>% .do(),
        core = genes %>% anchor_5p() %>% mutate(width=400) %>% shift_upstream(200) %>% .do(),
        ext = genes %>% anchor_5p() %>% mutate(width=3000) %>% shift_upstream(1500) %>% .do(),
        body = genes %>% .do()
    )

    sums = function(cg) {
        cg2 = setNames(cg[[2]], cg[[1]])
        mat2 = cgmat
        narray::intersect(cg2, mat2, along=1)
        narray::map(mat2, along=1, mean, subsets=cg2)
    }
    res = lapply(cgs, sums)

    saveRDS(res, file=fpath)
    res
}

if (is.null(module_name())) {
    .cohorts = import('./cohorts')$cohorts
    for (tissue in .cohorts())
        for (id_type in c("ensembl_gene_id", "external_gene_name"))
            res = meth_summary(tissue, id_type)
}
