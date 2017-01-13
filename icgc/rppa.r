idx = import('./indexing')

#' Function to retrieve the RPPA protein data from the processed ICGC object
#'
#' @param index       HDF5 index, either numerical or character;
#'                    If a character vector ICGC ids will be matched:
#'                    SA*: sample ID, SP*: specimen ID, SD*: donor ID
#' @param map_ids     character vector to map identifiers to: 'icgc_sample_id', 
#' @return             The requested sample matrix
rppa = function(index=available(rppa=TRUE), map_ids=TRUE) {
    idx$getHDF5(index=index, map_ids=map_ids, fname="protein")
}
