# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.b = import('ebits/base')
.io = import('ebits/io')
.ar = import('ebits/array')
.p = import('../../path')
util = import('./util')

#' Regular expression for MAF files
archive_regex = "Mutation_Packager_Calls.Level_3.*\\.tar(\\.gz)?$"

#' Metadata fields to keep
keep_fields = c('Chromosome','UniProt_AApos','Tumor_Validation_Allele2',
    'Transcript_Position', 'Score', 'chromosome_name_WU',
    'ucsc_cons_WU', 'Tumor_Seq_Allele2', 'variant_WU',
    'ucsc_cons', 'chromosome_name', 'Validation_Method',
    'NVarCov', 'NVarRat', 'COSMIC_Gene_Freq', 'CGC_Mutation_Type',
    'Reference_Allele', 'Tumor_Seq_Allele1', 'Match_Norm_Seq_Allele1',
    'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
    'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
    'variant', 'Sequencing_Phase')

#' Read a single MAF file and return results
#'
#' @param fname  File name to load
#' @param        Print file name currently processing
#' @return       An expression matrix with genes x samples
file2mut = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = .io$read_table(fname, header=TRUE, quote="") %catch% NULL
    if (identical(re, NULL))
        warning("Failed: ", fname)
    else {
        re$study = .b$grep("gdac.broadinstitute.org_([A-Z]+)", fname)

        for (nn in keep_fields) {
            if (!is.null(re[[nn]]))
                re[[nn]] = as.character(re[[nn]])
        }
    }
    re
}

#' Process all MAF files with voom
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Mutation matrix if save is NULL
mutations = function(regex=archive_regex, dir=util$data_dir, save=NULL) {
    elist = util$list_files(regex, dir) %>%
        util$unpack() %>%
        util$select("\\.maf\\.txt")

    mut = elist %>%
        lapply(file2mut) %>%
        bind_rows()

    if (is.null(save))
        mut
    else
        save(mut, file=.p$file("tcga", save))
}
