.b = import('ebits/base')
.io = import('ebits/io')

#' List of code lookup tables
#'
#' see https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm
.codes = list(center_code = 'center_code.txt',
              disease_study = 'disease_study.txt',
              platform_code = 'platform_code.txt',
              sample_type = 'sample_type.txt',
              tissue_source_site = 'tissue_source_site.txt',
              portion_analyte = 'portion_analyte.txt') %>%
    .b$lnapply(function(x) {
       fname = file.path(module_file(), x)
       .io$read_table(fname, header=TRUE, quote=NULL,
           na.strings=NULL, colClasses='character', check.names=TRUE)
    })

for(.i in seq_along(.codes))
    assign(names(.codes)[.i], .codes[[.i]])
