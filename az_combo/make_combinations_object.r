library(dplyr)
b = import('base')
io = import('io')
gdsc = import('data/gdsc')

# read all the .xls files and save them in a combined object
raw = list.files("raw/Combenefit_input", full.names=TRUE, recursive=TRUE) %>%
    lapply(function(fname) xlsx::read.xlsx(fname, 1, stringsAsFactors=FALSE, check.names=FALSE))

# generate index from fields, subset data, save both
data = lapply(raw, function(x) {
    re = data.matrix(x[1:6,2:7])
    rownames(re) = x[1:6,1]
    stopifnot(x[8:12,1,drop=TRUE] == c("Agent1", "Agent2", "Unit1", "Unit2", "Title"))
    names(dimnames(re)) = list(x[8,2], x[9,2])
    re
})

index = lapply(raw, function(x) {
    re = x[8:12,2] %>%
        sub("\\\\muM", "uM", .) %>%
        setNames(x[8:12,1]) %>%
        as.list() %>%
        as.data.frame()
    re$MaxConc1 = max(x[1:6,1])
    re$MaxConc2 = max(colnames(x)[2:7])
    re
}) %>% bind_rows()

# order data matrices the same as the index
filtered = index %>%
    mutate(index=1:nrow(.),
           Title = sub("[nN]ormal|[pP]arental", "", Title)) %>%
    filter(!grepl("FGF|LTED|F100|F\\.100|High", Title)) %>%
    cbind(., gdsc$cosmic$name2id(.$Title, table=TRUE)[,-1]) %>%
    na.omit() %>%
    transmute(COMPOUND_A = Agent1,
              COMPOUND_B = Agent2,
              UNIT_A = Unit1,
              UNIT_B = Unit2,
              MAX_CONC_A = MaxConc1,
              MAX_CONC_B = MaxConc2,
              CELL_LINE = from,
              COSMIC = to,
              index = index)

# combine with synergy score, etc table
syn = io$read_table("raw/MONO_HILL_COMBI_SYNERGY.csv", header=TRUE)
syn2 = syn
colnames(syn2) = sub("_A", "_X", colnames(syn2))
colnames(syn2) = sub("_B", "_A", colnames(syn2))
colnames(syn2) = sub("_X", "_B", colnames(syn2))
syn = bind_rows(syn, syn2)
combined = merge(filtered, syn, by=c('CELL_LINE','COMPOUND_A','COMPOUND_B','MAX_CONC_A','MAX_CONC_B'))

# subset data, index accordingly
data = data[combined$index]
names(data) = combined$index
index = combined
index$index = NULL
names(data) = 1:nrow(index)

# save r object
save(data, index, file="combinations.RData")
