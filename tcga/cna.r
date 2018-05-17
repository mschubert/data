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
    fname = sprintf("cn_%s-%s.RData", tissue, gene)
    fpath = file.path(module_file(), "cache", fname)
    if (id_type != "specimen")
        stop("only id_type='specimen' implemented with caching")

    if (file.exists(fpath)) {
        dset = .io$load(fname)
        genes = dset$genes
        copies = dset$copies
    } else {
        warning("no cache file found, this may take a long time",
                immediate.=TRUE)

        genes = .seq$coords$gene(idtype=gene, granges=TRUE)
        cnas = cna_segments(tissue, id_type=id_type, granges=TRUE) %>%
            plyranges::select(Sample, ploidy)

        copy_ranges = genes %>%
            plyranges::select(!! rlang::sym(gene)) %>%
            plyranges::join_overlap_intersect(cnas) %>%
            as.data.frame() %>% # plyranges takes >50GB of mem, why? (#36)
            dplyr::mutate(width = abs(start - end)) %>%
            dplyr::group_by(Sample, ensembl_gene_id) %>%
            dplyr::summarize(copies = weighted.mean(ploidy, width))

        copies = narray::construct(copies ~ ensembl_gene_id + Sample, data=copy_ranges)
        save(genes, copies, file=fpath)
    }

    if (length(chr_excl) > 0) {
        keep = as.data.frame(GenomicRanges::mcols(genes)) %>%
            dplyr::filter(! seqnames %in% chr_excl)
        copies = copies[, colnames(copies) %in% keep[[gene]]]
    }
    copies
}
