.io = import('ebits/io')
.cohorts = import('./cohorts')$cohorts
`%>%` = magrittr::`%>%`

#' Get a data.frame listing all clinical data available
#'
#' @param tissue   Limit selection to a single/set of tissues
#' @return         A data.frame with data for all the clinical data
clinical = function(tissue=NULL, id_type="patient") {
    if (is.null(tissue))
        tissue = .cohorts()
    fpath = module_file("TCGAbiolinks-downloader/clinical")
    load_fun = function(t) .io$load(file.path(fpath, paste0("TCGA-", t, ".RData")))
    res = lapply(tissue, load_fun) %>%
        dplyr::bind_rows() %>%
        dplyr::as_tibble()
}
