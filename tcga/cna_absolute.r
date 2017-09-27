# this is from synapse
# https://www.synapse.org/#!Synapse:syn1710464

#' Get ABSOLUTE copy numbers from synapse
#'
#' @return  Data.frame with 'Sample' ID and genomic regions
cna_absolute = function() {
    fpath = module_file("data/pancan12_absolute.segtab.txt")
    stopifnot(file.exists(fpath))
    cna = readr::read_tsv(fpath)
}
