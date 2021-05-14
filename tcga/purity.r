`%>%` = magrittr::`%>%`
.bc = import('./barcode')
.map_id = import('./map_id')$map_id

purity = function(...) {
    message("[tcga/purity] using purity_aran2015")
    purity_aran2015(...)
}

#' Get tumor purity estimates from Aran & Butte (2015, Nat Comms)
#'
#' Suppl. Data 1 from https://www.nature.com/articles/ncomms9971
#'
#' Data sources are the following:
#' * estimate: ESTIMATE method doi.org/10.1038/ncomms3612 (with additional cohorts)
#' * absolute: ABSOLUTE method doi.org/10.1038/nbt.2203
#' * leuko_meth: avg_beta_cpg(over 30% meth in cancer, under 5% in blood) / 0.85
#' * he_stain: H&E staining as described in 9971
#' * consensus: computed consensus estimate as described in 9971
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with purity estimates
purity_aran2015 = function(tissue=NULL, id_type="specimen") {
    nan2na_num = function(x) as.numeric(ifelse(is.nan(x), NA, x))

    fpath = module_file("TCGAbiolinks-downloader")
    purity = readxl::read_excel(file.path(fpath, "ncomms9971-s2.xlsx"), skip=3) %>%
        dplyr::transmute(cohort = `Cancer type`,
                         Sample = `Sample ID`,
                         estimate = nan2na_num(ESTIMATE),
                         absolute = nan2na_num(ABSOLUTE),
                         leuko_meth = nan2na_num(LUMP),
                         he_stain = nan2na_num(IHC),
                         consensus = nan2na_num(CPE)) %>%
        .map_id(id_type=id_type, along="Sample")

    if (is.null(tissue))
        purity
    else
        dplyr::filter(purity, cohort %in% tissue)
}
