b = import('ebits/base')
config = import('./config')

#' Function to map the given variables to an HDF5 index and name mappings
#'
#' @param index      Numeric, or character vector of IDs
#' @param map_ids    Type of IDs to map to, default: same as requested
#' @param filter_to  Set of IDs to map to
#' @return           List with file index and mapped column names
varmap = function(valid, map_ids=TRUE, filter_to=NULL) {
    if (is.character(map_ids))
        to = map_ids
    else if (identical(map_ids, FALSE) || all(grepl("^SA", filter_to)))
        to = "icgc_sample_id"
    else if (all(grepl("^SP", filter_to)))
        to = "icgc_specimen_id"
    else if (all(grepl("^D", filter_to)))
        to = "icgc_donor_id"
    else
        stop("invalid ids - all must start with SA, D, or SP")

    .b$match(x = valid,
             from = "icgc_sample_id",
             to = to,
             filter_to = filter_to,
             data = file.path(config$cached_data, 'clinicalsample.RData'))
}
