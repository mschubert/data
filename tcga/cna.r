library(plyranges) # TODO: remove when sa-lee/plyranges#36 is fixed
`%>%` = magrittr::`%>%`
.io = import('io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id
.seq = import('ebits/seq')

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
#' @param gene     ID type of the gene (e.g. hgnc_symbol, ensembl_gene_id)
#' @param chr_excl  Chromosomes to exclude; default: Y,MT
#' @return         A matrix of gene copy number variations
cna_genes = function(tissue, id_type="specimen", gene="ensembl_gene_id",
                     chr_excl=c("Y","MT")) {
    if (id_type != "specimen")
        stop("only id_type='specimen' implemented with caching")

    fpath = .cache_file(tissue, gene)
    if (file.exists(fpath))
        obj = .io$load(fpath)
    else {
        warning("no cache file found, this may take a long time", immediate.=TRUE)
        obj = .make_cache(cohort=tissue, feat_id=gene, fpath=fpath,
                          feat_ranges=.seq$coords$gene(idtype=gene, granges=TRUE))
    }
    .filter(obj, gene, chr_excl=chr_excl)
}

#' Make copy number cache of features
#'
#' @param fpath        File path of where to store the cache
#' @param cohort       Character vector of TCGA cohort, e.g. 'BRCA'
#' @param feat_id      Character string to serve as identifier, e.g. 'ensembl_gene_id'
#' @param feat_ranges  GRanges object with `feat_id` as column
#' @return             Contents of the cached file
.make_cache = function(fpath, cohort, feat_id, feat_ranges) {
    fsym = rlang::sym(feat_id)

    cnas = cna_segments(cohort, id_type='specimen', granges=TRUE) %>%
        plyranges::select(Sample, ploidy)

    copy_ranges = feat_ranges %>%
        plyranges::select(!! fsym) %>%
        plyranges::join_overlap_intersect(cnas) %>%
        as.data.frame() %>% # plyranges takes >50GB of mem, why? (#36)
        dplyr::mutate(width = abs(start - end)) %>%
        dplyr::group_by(Sample, !! fsym) %>%
        dplyr::summarize(copies = weighted.mean(ploidy, width))

    fml = as.formula(sprintf("copies ~ %s + Sample", feat_id))
    copies = narray::construct(fml, data=copy_ranges)
    save(feat_ranges, copies, file=fpath)
    list(feat_ranges, copies)
}

#' Return cache file path
#'
#' @param cohort   Character vector of TCGA cohort, e.g. 'BRCA'
#' @param feat_id  Character string to serve as identifier, e.g. 'ensembl_gene_id'
#' @return         Character string of the cache file path
.cache_file = function(cohort, feat_id) {
    fdir = file.path(module_file(), "cache", sprintf("cna-%s", feat_id))
    if (!dir.exists(fdir))
        dir.create(fdir)
    file.path(fdir, sprintf("%s.RData", cohort))
}

#' Return (filtered) cache
#'
#' @param obj  List with columns: feat_ranges [GRanges], copies [matrix, feat x sample]
#' @param feat_id  Character string to serve as identifier, e.g. 'ensembl_gene_id'
#' @param chr_excl  Chromosomes to exclude; default: Y,MT
#' @return  A filtered copy number matrix [features x samples]
.filter = function(obj, feat_id, chr_excl=c('Y', 'MT')) {
    copies = obj$copies
    if (length(chr_excl) > 0) {
        keep = as.data.frame(obj$feat_ranges) %>%
            dplyr::filter(! seqnames %in% chr_excl)
        copies = copies[rownames(copies) %in% keep[[feat_id]],]
    }
    copies
}
