b = import('base')
io = import('io')

#' Return a data.frame with all donors
#'
#' Fields include IDs, ages, gender, hardy death scale
#'
#' @param hardy  Which kinds of deaths to include (hardy scale)
#' @return       A data.frame
donors = function(hardy=0:4) {
	fname = "GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
    fpath = module_file("data", fname, mustWork=TRUE)
	io$read_table(fpath, header=TRUE) %>%
		dplyr::filter(DTHHRDY %in% hardy)
}
