`%>%` = magrittr::`%>%`
io = import('io')
seq = import('seq')
cosmic = import('./cosmic')

#' Calculate and return aneuploidy scores for the GDSC panel
#'
#' @return  A data.frame with fields: cell_line, cosmic, aneuploidy, coverage
aneuploidy = function() {
    fname = "Summary_segmentation_data_994_lines_picnic.txt"
    fpath = module_file(file.path("data", fname))
    cnas = readr::read_tsv(fpath) %>%
        dplyr::mutate(width = abs(startpos - endpos)) %>%
        seq$aneuploidy(sample="cellLine", ploidy="totalCN", seqnames="chr") %>%
        dplyr::transmute(
            cell_line = cellLine,
            cosmic = suppressWarnings(cosmic$name2id(cell_line)), # all ok
            aneuploidy = aneuploidy)
}
