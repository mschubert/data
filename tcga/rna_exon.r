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
rna_seq = function(tissue, id_type="specimen", trans="raw", annot=FALSE) {
    fpath = module_file(paste0("TCGAbiolinks-downloader/rna_exon_", trans))
    expr = .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData")))

    fannot = file.path(module_file("cache"), "exon_annot.rds")
    if (!file.exists(fannot)) {
        # load hg19 annotations from ensembl
        ensembl = biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", GRCh=37,
                                      host = "useast.ensembl.org")
        options(biomart.timeout = 3600)
        res = biomaRt::getBM(attributes=c("ensembl_gene_id", "ensembl_exon_id",
                                          "exon_chrom_start", "exon_chrom_end",
                                          "chromosome_name", "strand",
                                          "external_gene_name"),
                             mart=ensembl)
        # add rownames to exons, save in rds
        dannot = res %>%
            dplyr::mutate(chromosome_name = paste0("chr", chromosome_name),
                    strand = setNames(c("+", "-"), c(1,-1))[as.character(strand)]) #%>%
#            makeGRangesFromDataFrame(keep.extra.columns=TRUE)
#        ov = GenomicRanges::findOverlaps(dannot %>% filter(external_gene_name %in% c("RBM14")), expr@rowRanges)
#        names(dannot) = paste0("chr", seqnames(dannot), ":", start(dannot), "-", end(dannot), ":", strand(dannot))
        dannot$name = with(dannot,
                paste0(chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end, ":", strand))
        # save in rds
        saveRDS(dannot, file=fannot)
    } else
        dannot = readRDS(fannot)

    # add annotations to row data
    annot = dannot[match(rownames(expr), dannot$name),]
    GenomicRanges::mcols(expr@rowRanges) = annot[c("ensembl_exon_id", "ensembl_gene_id", "external_gene_name")]

    if (identical(annot, TRUE))
        expr
    else {
        re = .map_id(SummarizedExperiment::assay(expr), id_type=id_type)
        if (is.character(annot))
            rownames(re) = SummarizedExperiment::rowData(expr)[[annot]]
        re
    }
}
