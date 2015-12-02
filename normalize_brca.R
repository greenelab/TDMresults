source("directories.R")

args = c("BRCARNASeq.pcl", 
	"BRCAarray.pcl", 
	initial_data, 
	normalized_data, 
	"BRCA")

source("normalize_data.R")

# args = c(args[4], 
# 		'BRCA_ZEROONE.pcl', 
# 		'BRCAClin.tsv', 
# 		'BRCA_TDM_ZEROONE.pcl', 
# 		'BRCARNASeqClin.tsv', 
# 		'BRCA_QN_ZEROONE.pcl',
# 		'BRCARNASeqClin.tsv', 
# 		'BRCA_LOG_ZEROONE.pcl', 
# 		'BRCARNASeqClin.tsv', 
# 		paste0(args[4], "../output/"))

# source("lasso_subtype_brca.r")
