source("directories.R")

args = c("COADREADRNASeq.pcl", 
	"COADREADarray.pcl", 
	initial_data, 
	normalized_data, 
	"COADREAD")

source("normalize_data.R")

args = c(args[4], 
		'COADREAD_ZEROONE.pcl', 
		'COADREADClin.tsv', 
		'COADREAD_TDM_ZEROONE.pcl', 
		'COADREADRNASeqClin.tsv', 
		'COADREAD_QN_ZEROONE.pcl',
		'COADREADRNASeqClin.tsv', 
		'COADREAD_LOG_ZEROONE.pcl', 
		'COADREADRNASeqClin.tsv', 
		paste0(args[4], "../output/"))

source("lasso_subtype_coadread.R")
