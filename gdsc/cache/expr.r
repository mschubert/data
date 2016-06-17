library(dplyr)
ma = import('process/microarray')

expr = ArrayExpress::ArrayExpress("E-MTAB-3610", drop=TRUE) %>%
    ma$qc() %>%
    ma$normalize()
save(expr, file="expr_temp.RData")

#annot = data.frame(probe_id = rownames(expr)) %>%
#    mutate(hgnc = ma$annotate)
#
##TODO: add hgnc, ensembl, etc annotations in data.frame here for easy mapping
#annot = data.frame(
#    hgnc = ma$annotate(...)
#)

save(expr, file="expr_probes.RData")
