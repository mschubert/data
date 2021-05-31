`%>%` = magrittr::`%>%`

#' Get Pan-cancer immune measures from Thorsson et al.
immune = function() {
    fpath = module_file("data", mustWork=TRUE)
    readxl::read_xlsx(file.path(fpath, "1-s2.0-S1074761318301213-mmc2.xlsx"), na="NA") %>%
        dplyr::rename(barcode = `TCGA Participant Barcode`,
                      cohort = `TCGA Study`)
}
