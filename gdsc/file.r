.io = import('ebits/io')

.files = list(
    MASTER_LIST = "20140320_MASTER_LIST.ro",
    ANALYSIS_SET = "",
    BASAL_EXPRESSION = "BASAL_EXPRESSION_12062013v2.ro",
    DRUG_PROPS = "DRUG_ANALYSIS_SET_20150123.rdata",
    DRUG_IC50 = "v17a_IC50s.Rdata",
    DRUG_AUC = "v17a_AUCs.Rdata",
    INTOGEN_DRIVERS = "cancer_drivers_5_2_2014.ro",
    CONC = "SCREEN_CONC.RData",
    NGS_BEM = "NGS_BEM_COSMIC_NURIAS_26022014.ro"
)

get = function(id) .io$load(module_file('cache', .files[[id]], mustWork=TRUE))
