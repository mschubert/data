.b = import('ebits/base')
.file = import('./file')

DRUG_PROPS = .file$get('DRUG_PROPS')

#' Converts a vector of drug names to IDs
#'
#' @param x            Vector of drug names
#' @param fuzzy_level  Fuzzy matching: 0 (exact), 1 (alphanumeric exact), 2 (closest match)
#' @param table        Return a data.frame instead of a named vector
#' @return             Vector with IDs as values and drug names as names
name2id = function(x, fuzzy_level=1, table=FALSE) {
    .b$match(x = x,
             from = DRUG_PROPS$DRUG_NAME,
             to = DRUG_PROPS$DRUG_ID,
             fuzzy_level = fuzzy_level, table = table)
}

#' Converts a vector of drug IDs to names
#'
#' @param x      Vector of drug IDs
#' @param table  Return a data.frame instead of a named vector
#' @return       Vector with names as values and IDs as names
id2name = function(id, table=FALSE) {
    .b$match(x = as.character(id),
             from = DRUG_PROPS$DRUG_ID,
             to = DRUG_PROPS$DRUG_NAME,
             fuzzy_level = 0, table = table)
}

#' Returns the maximum screening concentration given a drug name or ID
#'
#' @param type   'max' or 'min' for screening concentrations
#' @param drugs  Character vector of drug names
#' @param ids    Character vector of drug IDs
#' @param log    Return log uM (TRUE) or uM (FALSE)
#' @param ...    Arguments passed to drug/name2id
#' @return       Maximum concentration that was screened
conc = function(type="max", names=NULL, ids=NULL, log=TRUE, ...) {
    if (is.null(names) + is.null(ids) != 1)
        stop("Need either drug names or IDs")

    CONC = .file$get('CONC')

    if (!is.null(ids)) {
        mapping = as.character(CONC$DRUG_ID)
        query = as.character(ids)
    } else {
        mapping = id2name(as.character(CONC$DRUG_ID))
        query = names
    }

    if (log==TRUE)
        field = "_LOG_UM"
    else
        field = "_UMOLAR"

    conc = setNames(CONC[[paste0(toupper(type), field)]], mapping)
    setNames(conc[query], query)
}

#' Returns a data.frame with drug ID and name, and the target classes
#'
#' @return  The data.frame
targets = function() {
    lookup = setNames(c("targeted", "cytotoxic"),
                      c("targetted", "cytotoxic"))

    fac2int = function(x) as.integer(as.character(x))

    DRUG_PROPS %>%
        filter(QC_DRUG == "P") %>%
        transmute(drug = as.character(DRUG_NAME),
                  id = as.character(DRUG_ID),
                  type = unname(lookup[as.character(DRUG_TYPE)]),
                  synonyms = as.character(SYNONYMS),
                  brand_name = as.character(BRAND_NAME),
                  target = as.character(PUTATIVE_TARGET),
                  p53 = fac2int(p53.pathway),
                  egfr = fac2int(EGFR.signaling),
                  mapk_erk = fac2int(ERK.MAPK.signaling),
                  pi3k = fac2int(PI3K.signaling),
                  mtor = fac2int(TOR.signaling),
                  rtk = fac2int(RTK.signaling),
                  wnt = fac2int(WNT.signaling),
                  igfr = fac2int(IGFR.signaling),
                  abl = fac2int(ABL.signaling),
                  jnk_p38 = fac2int(JNK.and.p38.signaling),
                  metabolism = fac2int(metabolism),
                  dna_replication = fac2int(DNA.replication),
                  cell_cycle = fac2int(cell.cycle),
                  mitosis = fac2int(mitosis),
                  cytoskeleton = fac2int(cytoskeleton),
                  genome_integrity = fac2int(Genome.integrity),
                  chromatin_histone = fac2int(chromatin..histone.methylation),
                  chromatin_other = fac2int(chromatin..other),
                  other = fac2int(other))
}
