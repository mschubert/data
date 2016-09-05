include ../Makefile

URLS = $(shell firehose_get -e -b data latest | grep "^Download" | cut -d' ' -f2,2)
COHORTS = $(shell echo $(URLS) | tr " " "\n" | cut -d'/' -f7,7)
#DIRS = 

# process all firehose data
firehose: data
	;

# download using the firehose tool
data:
	;

data/PrimTum_MoBEM_PanCan.tsv:
	wget -N -P data http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/BEMs/PrimaryTumours/PrimTum_MoBEM_PanCan.tsv

cache/BEM.RData: data/PrimTum_MoBEM_PanCan.tsv
	Rscript bem.r $^ $@	
