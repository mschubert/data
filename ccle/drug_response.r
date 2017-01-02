b = import('ebits/base')
io = import('ebits/io')
ar = import('ebits/array')
df = import('ebits/data_frame')
util = import('./util')
index = import('./index')$index

#' Returns a drug response matrix (drugs x COSMIC IDs)
#'
#' Note that an IC50 value of 8 corresponds to the maximum assigned
#' resistance and that the real value might be higher.
#'
#' @param cosmic      Use COSMIC IDs as cell line identifiers (default: FALSE)
#' @param gdsc_names  Use GDSC drug names (default: FALSE)
#' @param drop_max    Remove capped response values (8 uM)
#' @return            The drug response matrix
drug_response = function(cosmic=FALSE, gdsc_names=FALSE, drop_max=FALSE) {
    if (cosmic)
        idx = select(index, `CCLE name`, index=COSMIC)
    else
        idx = select(index, `CCLE name`, index=`Expression arrays`)

    re = util$read("data", "CCLE_NP24.2009_Drug_data_2015.02.24.csv", header=TRUE) %>%
        select(`CCLE name`=`CCLE Cell Line Name`, Compound, IC50=`IC50 (uM)`) %>%
        df$update(idx, match_cols="CCLE name") %>%
        na.omit() %>%
        ar$construct(IC50 ~ index + Compound, data=.)

    if (drop_max)
        re[re >= 8 - .Machine$double.eps] = NA

    if (gdsc_names) {
        fname = io$file_path(module_file(), "drug_mapping.txt")
        lookup = io$read_table(fname, header=TRUE)
        colnames(re) = lookup$GDSC_NAME[match(colnames(re),
                lookup$CCLE_NAME)] %or% colnames(re)
    }
    re
}
