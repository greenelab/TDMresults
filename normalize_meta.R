source("directories.R")

args = c("BRCARNASeq.pcl", 
	"METABRICFiltered.pcl", 
	initial_data, 
	normalized_data, 
	"META_BRCA",
	"BRCAarray.pcl")

source("normalize_data.R")

# args = c(args[4], 
# 		'METABRICFiltered_ZEROONE.pcl', 
# 		'combined_subtype.txt', 
# 		'META_BRCA_MA_ZEROONE.pcl', 
# 		'BRCAClin.tsv', 
# 		'META_BRCA_LOG_ZEROONE.pcl', 
# 		'BRCARNASeqClin.tsv', 
# 		'META_BRCA_QN_ZEROONE.pcl', 
# 		'BRCARNASeqClin.tsv', 
# 		'META_BRCA_TDM_ZEROONE.pcl', 
# 		'BRCARNASeqClin.tsv', 
# 		paste0(args[4],"../output/"))

# source("lasso_subtype_meta.r")
