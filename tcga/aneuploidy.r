`%>%` = magrittr::`%>%`
sys = import('sys')
seq = import('seq')
.cohorts = import('./cohorts')$cohorts
.cna_segments = import('./cna')$cna_segments
.cna_absolute = import('./cna_old')$cna_absolute

#' Provide AneuFinder-like aneuploidy score
#'
#' @param ranges  GRanges object returned by to_granges function
#' @param per_chromosome  Calculate aneuploidy score per chromosome or genome
#' @param ...     Arguments passed to `seq$aneuploidy`
#' @return  ...
aneuploidy_log2seg = function(cohort=NULL, assembly="GRCh38", per_chromosome=FALSE, ...) {
    one_cohort = function(coh) {
        .cna_segments(coh) %>%
            dplyr::mutate(width = abs(End - Start)) %>%
            seq$aneuploidy(assembly=assembly, per_chromosome=per_chromosome,
                           sample="Sample", seqnames="Chromosome")
    }

    if (is.null(cohort))
        cohort = .cohorts()
    lapply(cohort, one_cohort) %>%
        dplyr::bind_rows()
}

aneuploidy_absolute = function(cohort=NULL) {
    if (is.null(cohort))
        cohort = .cohorts()

    # would be nice to get WGD here but no clear data release when mentioned
    # also ploidy is rather limited and taken from old pub
    .cna_absolute(cohort) %>%
        group_by(Sample) %>%
            summarize(ploidy_abs = weighted.mean(Modal_HSCN_1 + Modal_HSCN_2, Length),
                      aneup_abs = weighted.mean(abs(Modal_HSCN_1 + Modal_HSCN_2 - 2), Length)) %>%
        ungroup()
}

aneuploidy = function(cohort=NULL) {
    seg = aneuploidy_log2seg(cohort) %>% dplyr::select(Sample, aneup_log2seg=aneuploidy)
    abs = aneuploidy_absolute(cohort)
    dplyr::full_join(seg, abs, by="Sample")
}
