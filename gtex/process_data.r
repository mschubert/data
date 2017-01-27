io = import('io')
l1k = import('io/l1ktools_io')
idmap = import('process/idmap')
rnaseq = import('ebits/process/rna-seq')

#' Generate the RNA expression object
#'
#' Reads the Ensembl ID gene counts from the original GTEx .gct expression
#' file, maps the IDs to hgnc symbols, runs a variance stabilizing trans-
#' formation (DESeq2) on the resulting matrix
rna_seq = function() {
	fname = "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
    expr = l1k$read.gct(module_file("data", fname, mustWork=TRUE))@mat

	# sum all ensembl genes for each hgnc symbol to get counts
	expr = idmap$gene(expr, to="hgnc_symbol", summarize=sum)
	expr = rnaseq$vst(expr)

	fout = file.path(module_file("cache", mustWork=TRUE), "rna_seq_vst.gctx")
	io$save(expr, file=fout)
}

#' Read the text-based index and save a .RData object
samples = function() {
	fname = "GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
	fpath = module_file("data", fname, mustWork=TRUE)
	fout = file.path(module_file("cache", mustWork=TRUE), "samples.RData")
	
	samples = io$read_table(fpath, header=TRUE)
	io$save(samples, file=fout)
}
