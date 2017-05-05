.b = import('ebits/base')
.io = import('ebits/io')

.camel2underscore = function(string) {
    s = gsub("([a-z])([A-Z])", "\\1_\\L\\2", string, perl=TRUE)
    sub("^(.[a-z])", "\\L\\1", s, perl = TRUE)
}

#' List of code lookup tables
#'
#' wget https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables
#' unzip $_
.codes = list.files(module_file(), pattern="\\.tsv", full.names=TRUE) %>%
    setNames(., .camel2underscore(sub("\\.tsv", "", basename(.)))) %>%
    .b$lnapply(function(fname) {
       .io$read_table(fname, header=TRUE, quote=NULL,
           na.strings=NULL, colClasses='character', check.names=TRUE)
    })

for(.i in seq_along(.codes))
    assign(names(.codes)[.i], .codes[[.i]])
