`%>%` = magrittr::`%>%`
.map_id = import('./map_id')$map_id

#' Get tumor purity estimates
#'
#' from: https://www.nature.com/articles/ncomms9971
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame or GRanges object (Segment_Mean is log2(copy-number)-1)
purity = function(tissue, id_type="specimen", granges=FALSE) {
    nan2na_num = function(x) as.numeric(ifelse(is.nan(x), NA, x))

    fpath = module_file("TCGAbiolinks-downloader")
    purity = readxl::read_excel(file.path(fpath, "ncomms9971-s2.xlsx"), skip=3) %>%
        dplyr::transmute(cohort = `Cancer type`,
                         Sample = `Sample ID`,
                         estimate = nan2na_num(ESTIMATE),
                         absolute = nan2na_num(ABSOLUTE),
                         lump = nan2na_num(LUMP),
                         IHC = nan2na_num(IHC),
                         CPE = nan2na_num(CPE)) %>%
        .map_id(id_type=id_type, along="sample")
}
