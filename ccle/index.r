.b = import('ebits/base')
.io = import('ebits/io')
.ar = import('ebits/array')
.df = import('ebits/data_frame')
.gdsc = import('../gdsc')

.exists = function(...) !is.na(module_file(...))
.file = function(...) file.path(module_file(), ...)
.read = function(...) io$read(module_file(..., mustWork=TRUE))

index = .read("data", "CCLE_sample_info_file_2012-10-18.txt", header=TRUE)
index$COSMIC = .gdsc$cosmic$name2id(index$Cell.line.primary.name, warn=FALSE)

#' Returns a gene expression matrix (genes x COSMIC IDs)
#'
#' @param index_type  The column of `index` to be used for column names
#' @return            The expression matrix
basal_expression = function(index_type="COSMIC") {
    fname = "CCLE_Expression.Arrays_2013-03-18"
    if (!.exists("cache", fname, ext=".RData")) {
        warning("data file does not exist, creating raw .CEL files")
        ma = import('ebits/process/microarray')
        fnames = list.files(fname, full.names=TRUE)
        expr = oligo::read.celfiles(fnames) %>%
            ma$normalize() %>%
            ma$annotate("hgnc_symbol")
        save(expr, file=.file("cache", fname, ext=".RData"))
    } else
        expr = .io$load(.file("cache", fname, ".RData"))

    expr = Biobase::exprs(expr)
    cn = .b$match(sub(".CEL$", "", colnames(expr)),
                  from = index$Expression.arrays,
                  to = index[[index_type]])
    
    colnames(expr) = unname(cn)
    nas = is.na(cn)
    if (any(nas)) {
        warning("dropping ", sum(nas), " identifiers for ", index_type, " mapping")
        expr = expr[,!nas]
        cn = cn[!nas]
    }

    expr
}

#' Returns a drug response matrix (drugs x COSMIC IDs)
#'
#' Note that an IC50 value of 8 corresponds to the maximum assigned
#' resistance and that the real value might be higher.
#'
#' @param index_type  The column of `index` to be used for column names
#' @return            The drug response matrix
drug_response = function(index_type="COSMIC", map_drugs=TRUE) {
    idx = select_(index, "CCLE.name", index=index_type)
    re = .read("data", "CCLE_NP24.2009_Drug_data_2015.02.24.csv", header=TRUE) %>%
        select(CCLE.name=CCLE.Cell.Line.Name, Compound, IC50=IC50..uM.) %>%
        .df$update(idx, match_cols="CCLE.name") %>%
        na.omit() %>%
        .ar$construct(IC50 ~ index + Compound, data=.)

    if (map_drugs) {
        fname = .io$file_path(module_file(), "drug_mapping.txt")
        lookup = .io$read_table(fname, header=TRUE)
        colnames(re) = lookup$GDSC_NAME[match(colnames(re),
                lookup$CCLE_NAME)] %or% colnames(re)
    }
    re
}
