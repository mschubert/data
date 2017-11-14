library(dplyr)
.file = import('./file')
cosmic = import('./cosmic')
tissues = import('./tissues')$tissues
drug = import('./drug')
DRUG_PROPS = drug$DRUG_PROPS

#' Returns a drug response matrix using the filters specified
#'
#' @param metric         Either 'IC50s' or 'AUC'
#' @param filter_cosmic  A vector of COSMIC IDs to filter for (default: TRUE, include all)
#' @param drug_names     Boolean flag indicating whether to name drugs or use IDs
#' @param cell_names     Boolean flag indicating whether to use cell line names or COSMIC IDs
#' @param min_tissue_measured  Minimum number of measured responses per tissue, NA otherwise
#' @param drop           Remove columns that only contain NAs
#' @param median_top     Include only drug responses where the tissue median is in the top N tissues
#' @param stage          The minimum clinical stage; 0: experimental, 1: in dev; 2: approved
#' @return               A filtered and ID-mapped drug response matrix
drug_response = function(metric='IC50s', tissue=NULL, filter_cosmic=TRUE, drug_names=TRUE,
        cell_names=FALSE, min_tissue_measured=0, drop=FALSE, median_top=NA, stage=0) {
    if (grepl("IC50", metric))
        SCREENING = .file$get('DRUG_IC50')
    else if (grepl("AUC", metric))
        SCREENING = .file$get('DRUG_AUC')
    else
        stop("invalid metric")

    ar = import('ebits/array')
    io = import('ebits/io')
    tissues = tissues()
    stages = io$read_table(module_file('drugs_s1f.csv'), header=TRUE)

    if (is.numeric(median_top)) {
        ar$intersect(SCREENING, tissues, along=1)
        tissue_ranks = ar$map(SCREENING, along=1, subsets=tissues,
            function(x) median(x, na.rm=TRUE)) %>%
            ar$map(along=1, function(x) rank(x, ties.method="min"))

        for (tt in unique(tissues))
            for (did in colnames(SCREENING))
                if (tissue_ranks[tt,did] > median_top)
                    SCREENING[tt == tissues, did] = NA
    }

    if (!is.null(tissue)) {
        tissues = tissues(tissue)
        ar$intersect(SCREENING, tissues, along=1)
    }

    if (min_tissue_measured > 0) {
        if (grepl("AUC", metric))
            stop("concentration measurements only for IC50")

        for (tt in unique(tissues))
            for (did in colnames(SCREENING))
                if (sum(SCREENING[tt==tissues, did] <
                        drug$conc('max',ids=did), na.rm=TRUE) < min_tissue_measured)
                    SCREENING[tt==tissues, did] = NA
    }

    if (stage > 0) {
        if (stage == 1)
            keep = stages %>%
                filter(`Clinical Stage` != "experimental") %>%
                select(Identifier) %>%
                unlist() %>% unname()
        if (stage == 2)
            keep = stages %>%
                filter(`Clinical Stage` == "clinically approved") %>%
                select(Identifier) %>%
                unlist() %>% unname()

        SCREENING = SCREENING[,intersect(colnames(SCREENING), as.character(keep))]
    }

    if (drug_names)
        colnames(SCREENING) = drug$id2name(colnames(SCREENING))

    if (cell_names)
        rownames(SCREENING) = cosmic$id2name(rownames(SCREENING))

    if (drop)
        SCREENING = SCREENING[,apply(SCREENING, 2, function(x) !all(is.na(x)))]

    SCREENING
}
