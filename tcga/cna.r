`%>%` = magrittr::`%>%`
.io = import('io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id
.seq = import('ebits/seq')

#' Get copy number segments
#'
#' docs: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
#'
#' @param cohort   The cohort(s) to get expression for
#' @param barcode  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame or GRanges object (Segment_Mean is log2(copy-number)-1)
cna_segments = function(cohort, barcode="specimen", granges=FALSE) {
    fpath = module_file("TCGAbiolinks-downloader/cnv_segments")
    cna = .io$load(file.path(fpath, paste0("TCGA-", cohort, ".RData"))) %>%
        .map_id(id_type=barcode, along="Sample") %>%
        dplyr::mutate(ploidy = 2^(Segment_Mean + 1))

    if (granges)
        cna = GenomicRanges::makeGRangesFromDataFrame(cna,
            start.field="Start", end.field="End", keep.extra.columns=TRUE)
    cna
}

#' Get copy number segments
#'
#' docs: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
#'
#' @param barcode  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame or GRanges object
cna_segments_ascat = function(barcode="specimen", granges=FALSE) {
    fpath = module_file("cache")
    cna = readRDS(file.path(fpath, "cnv_segments_ascat2.rds")) %>%
        .map_id(id_type=barcode, along="Sample")

    if (granges)
        cna = GenomicRanges::makeGRangesFromDataFrame(cna,
            start.field="Start", end.field="End", keep.extra.columns=TRUE)
    cna
}

#' Get copy number for chromosomes
#'
#' @param cohort   The cohort(s) to get expression for
#' @param barcode  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param chr_excl  Chromosomes to exclude; default: Y,MT
#' @return         A matrix of gene copy number variations
cna_chrs = function(cohort, barcode="specimen", chr_excl=c("Y","MT")) {
    feat_ranges =.seq$coords$chrs(granges=TRUE)
    obj = .cached(cohort, "chrs", feat_ranges)
    .filter(obj, "chrs", barcode, chr_excl=chr_excl)
}

#' Get copy number for chromosome arms
#'
#' @inheritParams  cna_chrs
cna_chr_arms = function(cohort, barcode="specimen", chr_excl=c("Y","MT")) {
    feat_ranges =.seq$coords$chr_arms(granges=TRUE)
    obj = .cached(cohort, "arm", feat_ranges)
    .filter(obj, "arm", barcode, chr_excl=chr_excl)
}

#' Get copy number for chromosome bands
#'
#' @inheritParams  cna_chrs
cna_bands = function(cohort, barcode="specimen", chr_excl=c("Y","MT")) {
    feat_ranges =.seq$coords$chr_bands(granges=TRUE)
    obj = .cached(cohort, "band", feat_ranges)
    .filter(obj, "band", barcode, chr_excl=chr_excl)
}

#' Get copy number for genes
#'
#' @param gene     ID type of the gene (e.g. hgnc_symbol, ensembl_gene_id)
#' @inheritParams  cna_chrs
cna_genes = function(cohort, barcode="specimen", gene="ensembl_gene_id",
                     chr_excl=c("Y","MT")) {
    feat_ranges =.seq$coords$gene(idtype=gene, granges=TRUE)
    obj = .cached(cohort, gene, feat_ranges)
    .filter(obj, gene, barcode, chr_excl=chr_excl)
}

#' Get copy numbers for a custom feature set
#'
#' @param cohort       Character vector of TCGA cohort, e.g. 'BRCA'
#' @param feat_id      Column name of GRanges to use as idenitifer
#' @param feat_ranges  GRanges object with `feat_id` as column
#' @param barcode      Where to cut the barcode, either "patient", "specimen", or "full"
#' @param as_array     logical, indicating on whether to convert result to matrix
#' @return             Copy number object
cna_custom = function(cohort, feat_id, feat_ranges, barcode="specimen", as_array=TRUE) {
    fsym = rlang::sym(feat_id)

    cnas = cna_segments(cohort, barcode=barcode, granges=TRUE) %>%
        plyranges::select(Sample, ploidy)

    copy_ranges = feat_ranges %>%
        plyranges::select(!! fsym) %>%
        plyranges::join_overlap_intersect(cnas) %>%
        as.data.frame() %>%
        dplyr::mutate(width = abs(start - end)) %>%
        dplyr::group_by(Sample, !! fsym) %>%
        dplyr::summarize(copies = weighted.mean(ploidy, width, na.rm=TRUE))

    if (as_array) {
        fml = as.formula(sprintf("copies ~ %s + Sample", feat_id))
        copy_ranges = narray::construct(fml, data=copy_ranges)
    }

    copy_ranges
}

#' Return a cached copy of request (and create if not available)
#'
#' @param fpath        File path of where to store the cache
#' @param cohort       Character vector of TCGA cohort, e.g. 'BRCA'
#' @param feat_id      Character string to serve as identifier, e.g. 'ensembl_gene_id'
#' @param feat_ranges  GRanges object with `feat_id` as column
#' @return             Contents of the cached file
.cached = function(cohort, feat_id, feat_ranges) {
    fdir = file.path(module_file(), "cache", sprintf("cna-%s", feat_id))
    if (!dir.exists(fdir))
        dir.create(fdir)
    fpath = file.path(fdir, sprintf("%s.rds", cohort))

    if (file.exists(fpath))
        return(readRDS(fpath))

    warning(fpath, " cache not found, this may take a long time", immediate.=TRUE)
    copies = cna_custom(cohort, feat_id, feat_ranges, as_array=TRUE)
    obj = list(feat_ranges=feat_ranges, copies=copies)
    saveRDS(obj, file=fpath)
    obj
}

#' Filter a copy object by TCGA barcode type and chromosomes
#'
#' @param obj  List with columns: feat_ranges [GRanges], copies [matrix, feat x sample]
#' @param feat_id  Character string to serve as identifier, e.g. 'ensembl_gene_id'
#' @param barcode  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param chr_excl  Chromosomes to exclude; default: Y,MT
#' @return  A filtered copy number matrix [features x samples]
.filter = function(obj, feat_id, barcode, chr_excl=c('Y', 'MT')) {
    keep = as.data.frame(obj$feat_ranges) %>%
        dplyr::filter(! seqnames %in% chr_excl)
    copies = .map_id(obj$copies, id_type=barcode, along=2)
    copies[rownames(copies) %in% keep[[feat_id]],]
}
