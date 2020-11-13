`%>%` = magrittr::`%>%`
sys = import('sys')
seq = import('seq')
.cna_segments = import('./cna')$cna_segments

#' Provide AneuFinder-like aneuploidy score
#'
#' @param ranges  GRanges object returned by to_granges function
#' @param per_chromosome  Calculate aneuploidy score per chromosome or genome
#' @param ...     Arguments passed to `seq$aneuploidy`
#' @return  ...
aneuploidy = function(cohort=NULL, assembly="GRCh38", per_chromosome=FALSE, ...) {
    cna = .cna_segments(cohort) %>%
        dplyr::mutate(width = abs(End - Start)) %>%
        seq$aneuploidy(assembly=assembly, per_chromosome=per_chromosome,
                       sample="Sample", seqnames="Chromosome")
}
