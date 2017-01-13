idx = import('./indexing')

#' Function to retrieve the RNA seq data from the processed ICGC object
#'
#' @param index       HDF5 index, either numerical or character;
#'                    If a character vector ICGC ids will be matched:
#'                    SA*: sample ID, SP*: specimen ID, SD*: donor ID
#' @param raw_counts  Get the raw counts (as opposed to normalized)
#' @param map_ids     character vector to map identifiers to: 'icgc_sample_id', 
#'                     'icgc_specimen_id', 'donor_id' [default: same as requested
#'                    identifiers]
#' @return            The requested sample matrix
rna_seq = function(index=available(rna_seq=TRUE), raw_counts=FALSE, voom=FALSE,
                   map_ids=TRUE) {
    if (voom)
        fname = "expr_seq_voom"
    else if (raw_counts)
        fname = "expr_seq_raw"
    else
        fname = "expr_seq_norm"

    idx$getHDF5(index=index, map_ids=map_ids, fname=fname)
}
