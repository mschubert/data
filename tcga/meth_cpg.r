import_package("plyranges", attach=TRUE)
.io = import('io')
.seq = import('seq')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id
`%>%` = magrittr::`%>%`

#' Get gene annotations for CpGs for HM450
#'
#' @return  A data.frame with CpG annotations
meth_cg2gene = function() {
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
meth_cg2island = function() {
    fname = "HM450.hg38.manifest.CpGIsland.tsv.gz"
    fpath = file.path(module_file("cache"), fname)
    if (!file.exists(fpath)) {
        url = paste0("http://zhouserver.research.chop.edu/InfiniumAnnotation/20210414/HM450/", fname)
        download.file(url, fpath)
    }
    readr::read_tsv(fpath)
}

#' Get CpGs for TF binding sites
meth_cg2tf = function() {
    fname = "HM450.hg19.ENCODE.TFBS.tsv.gz"
    fpath = file.path(module_file("cache"), fname)
    if (!file.exists(fpath)) {
        url = paste0("http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/", fname)
        download.file(url, fpath)
    }
    readr::read_tsv(fpath)
}

#' Define core, extended promoter; gene and gene body
meth_cg2custom = function() {
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
#' @param tissue  The tissue(s) to get expression for
#' @param gene    Gene identifier: ensembl_{gene,transcript}_id, external_gene_name
#' @return        A matrix with summarized methylation: (core) promoter, body, etc.
meth_summary = function(tissue, gene="ensembl_gene_id") {
    one_cohort = function(tissue) {
        if (gene %in% c("gene_name", "hgnc_symbol"))
            gene = "external_gene_name"
        sid = rlang::sym(gene)
        params = sprintf("%s-%s.rds", tissue, gene)
        fpath = file.path(module_file("cache", "meth_summary"), params)
        if (file.exists(fpath))
            return(readRDS(fpath))
        cgmat = meth_cpg(tissue)

        c2g = m2g %>%
            filter(probe_id %in% rownames(cgmat)) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
        .do = . %>%
            join_overlap_intersect(c2g) %>% as.data.frame() %>% dplyr::as_tibble() %>%
            select(probe_id, !!sid) %>% filter(!!sid != "") %>% dplyr::distinct() %>% unstack()
        cgs = list(
            gene = genes %>% anchor_3p() %>% stretch(1500) %>% .do(),
            core = genes %>% anchor_5p() %>% mutate(width=400) %>% shift_upstream(200) %>% .do(),
            ext = genes %>% anchor_5p() %>% mutate(width=3000) %>% shift_upstream(1500) %>% .do(),
            body = genes %>% .do()
        )

        sums = function(cg) {
            gmeans = function(cgi) colMeans(cgmat[rownames(cgmat) %in% cgi,,drop=FALSE])
    #        res = pbapply::pblapply(unstack(cg), gmeans)
            res = parallel::mclapply(cg, gmeans)
            do.call(rbind, res)
        }
        res = lapply(cgs, sums)

        saveRDS(res, file=fpath)
        res
    }

    m2g = meth_cg2gene() %>%
        dplyr::select(seqnames=CpG_chrm, start=CpG_beg, end=CpG_end,
                      strand=probe_strand, probe_id=probeID)

    genes = .seq$gene_table() %>%
        dplyr::select(chromosome_name, strand, start=start_position,
                      end=end_position, ensembl_gene_id, external_gene_name,
                      ensembl_transcript_id) %>%
        dplyr::mutate(chromosome_name = paste0("chr", chromosome_name),
                      strand = setNames(c("+","-"),c(1,-1))[as.character(strand)]) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    res = lapply(tissue, one_cohort)
    list(gene = lapply(res, function(r) r$gene) %>% narray::stack(along=2),
         core = lapply(res, function(r) r$core) %>% narray::stack(along=2),
         ext = lapply(res, function(r) r$ext) %>% narray::stack(along=2),
         body = lapply(res, function(r) r$body) %>% narray::stack(along=2))
}

if (is.null(module_name())) {
    .cohorts = import('./cohorts')$cohorts
    options(mc.cores = 10L)
    for (tissue in .cohorts())
        for (id_type in c("ensembl_gene_id", "external_gene_name"))
            res = meth_summary(tissue, id_type)
}
