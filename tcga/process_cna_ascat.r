library(dplyr)

meta = module_file("data", "bcr_biotab") %>%
    list.files("biospecimen_aliquot.*\\.txt$", recursive=TRUE, full.names=TRUE) %>%
    lapply(function(f) readr::read_tsv(f)[-1,] %>%
           transmute(Sample = bcr_aliquot_barcode,
                     GDC_Aliquot=bcr_aliquot_uuid)) %>%
    bind_rows()

segs = module_file("data", "ascat2_segments") %>%
    list.files("\\.seg\\.txt$", recursive=TRUE, full.names=TRUE) %>%
    lapply(readr::read_tsv) %>%
    bind_rows() %>%
    inner_join(meta) %>%
    select(-GDC_Aliquot) %>%
    select(Sample, everything()) %>%
    mutate(Chromosome = sub("^chr", "", Chromosome))

saveRDS(segs, file=file.path(module_file("cache"), "cnv_segments_ascat2.rds"))

genes = module_file("data", "ascat2_genes") %>%
    list.files("\\.tsv$", recursive=TRUE, full.names=TRUE) %>%
    tibble(cohort = sub("^([A-Z]+-[^.]+)\\..*", "\\1", basename(.)),
           checksum = basename(dirname(.)),
           fname = .)

for (i in seq_along(genes)) {
    cohort = genes[[i]] %>%
        lapply(readr::read_tsv) %>%
        bind_rows()
}

x = readr::read_tsv(genes[1])
