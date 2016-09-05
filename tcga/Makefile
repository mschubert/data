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
