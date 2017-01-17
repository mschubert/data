all:
	exit "need to specify target"

../cache/rna_seq2_vst.gctx: ./rna_seq2_vst.r
	Rscript $<

../cache/clinical.RData: ./clinical.r
	Rscript $<

../cache/mutations.RData: ./mutations.r
	Rscript $<

../cache/cna.RData: ./cna.r
	Rscript $<

../cache/cna_thresholded.RData: ./cna_thresholded.r
	Rscript $<
