# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.b = import('ebits/base')
.io = import('ebits/io')
.ar = import('ebits/array')
util = import('./util')

#' Regular expression for MAF files
archive_regex = "Mutation_Packager_Calls.Level_3.*\\.tar(\\.gz)?$"

#' Metadata fields to keep
keep_fields = c('Study', 'Hugo_Symbol', 'Entrez_Gene_Id', 'NCBI_Build',
        'Chromosome', 'Start_position', 'End_position', 'Strand',
        'Variant_Classification', 'Variant_Type',
        'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
        'dbSNP_RS', 'dbSNP_Val_Status',
        'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
        'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2',
        'Mutation_Status', 'Sequence_Source', 'COSMIC_Codon',
        'COSMIC_Gene', 'Transcript_Id', 'Exon', "ChromChange", 'AAChange')

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
        re$Study = .b$grep("gdac.broadinstitute.org_([A-Z]+)", fname)
        re = re[intersect(keep_fields, colnames(re))]
        for (nn in intersect(keep_fields, colnames(re)))
            re[[nn]] = as.character(re[[nn]])
    }
    re
}

#' Process all MAF files with voom
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Mutation matrix if save is NULL
mutations = function(regex=archive_regex, dir=util$data_dir) {
    elist = util$list_files(regex, dir) %>%
        util$unpack() %>%
        util$select("\\.maf\\.txt")

    mut = elist %>%
        lapply(file2mut) %>%
        bind_rows()

    mut
}
