idx = import('./indexing')

#' Function to retrieve the mutation data from the processed ICGC object
#'
#' Can do subsetting using either `index`, `samples`, `specimen`, or `donors`.
#'
#' @param index    HDF5 index, either numerical or character;
#'                 If a character vector ICGC ids will be matched:
#'                 SA*: sample ID, SP*: specimen ID, SD*: donor ID
#' @param map_ids  character vector to map identifiers to: 'icgc_sample_id', 
#' @param minN     Minimum number of mutated genes to be included
#' @return         The requested sample matrix
mutations = function(index, map_ids=TRUE, minN=0) {
    re = idx$getHDF5(index=index, map_ids=map_ids, fname="mutations")
    re[rowSums(re) > minN,]
}
