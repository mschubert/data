library(dplyr)
b = import('base')
io = import('io')
gdsc = import('data/gdsc')

dat = io$read_table('raw/All_Norm_Combo_data_PilotDream.csv', header=TRUE) %>%
    filter(!is.na(Synergy_Score)) %>%
    transmute(drug1 = CpdA,
              drug2 = CpdB,
              synergy_score = Synergy_Score,
              cell_line = Cell_Line,
              az_panel = Cell_Panel) %>%
    mutate(cell_line = sub("[nN]ormal|[pP]arental", "", cell_line)) %>%
    mutate(cell_line = sub("- Rerun", "", cell_line)) %>%
    filter(!grepl("FGF|LTED|F100|F\\.100|High", cell_line)) %>%
	mutate(cosmic = gdsc$cosmic$name2id(cell_line, warn=FALSE)) %>%
	mutate(cell_line = gdsc$cosmic$id2name(cosmic)) %>%
	filter(!is.na(cell_line))

#alldrugs = sort(unique(c(dat$drug1, dat$drug2)))

io$write_table(dat, file="combo_resp.txt")
