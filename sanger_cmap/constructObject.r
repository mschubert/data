library(lumi)
library(limma)
library(dplyr)
b = import('base')
io = import('io')
gdsc = import('data/gdsc')
# implicit: lumiHumanIDMapping

DATA = module_file("data/r_objects/1st_2nd_3rd_batch")

arrays = io$load(io$file_path(DATA, 'quantileNormLumi.ro'))
annotation = nuID2IlluminaID(arrays)
stopifnot(rownames(annotation) == rownames(arrays))
expr = exprs(arrays)
rownames(expr) = annotation[,'Symbol']
expr = avereps(expr)

targets = io$load(io$file_path(DATA, 'targets.ro'))
stopifnot(targets$name == colnames(expr))
colnames(expr) = targets$sampleID

index = targets %>%
    dplyr::select(name, treatment) %>%
    dplyr::mutate(cell_line = sub("_[^_]+$", "", treatment),
                  cosmic = gdsc$cosmic$name2id(cell_line),
                  drug_id = sub("^[^_]+_", "", treatment),
                  drug = gdsc$drug$id2name(drug_id))
