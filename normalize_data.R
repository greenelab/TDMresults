#args = commandArgs(trailingOnly=TRUE)

target_file = args[1]
ref_file = args[2]
data_folder = args[3]
output_folder = args[4]
prefix = args[5]

if (length(args) > 5){
	tcga_file = args[6]
} else {
	tcga_file = NULL
}

source("package_loader.R")
source("tdm.R")

load_it(c("data.table", "gdata", "preprocessCore", "scales", "huge"))

message("Loading reference file...", appendLF=FALSE)

# Load ref_file.
# Read the first line of the reference file.
ref_head = readLines(paste0(data_folder, ref_file), n=1)
		
# Split all values in the header of the reference file by tabs.
split_head = unlist(strsplit(ref_head, "\t"))
	
# Read in the reference expression values.
ref_values = fread(paste0(data_folder, ref_file), header=F, skip=1, data.table=T)
		
# If there is no row name for the header, then handle it.
if(ncol(ref_values) == length(split_head)) {
	ref_head = split_head[2:length(split_head)]
} else {
	ref_head = split_head
}
rm(split_head)
		
# Add a rowname to the header, just to make things consistent.
setnames(ref_values, colnames(ref_values), c("gene", ref_head))
rm(ref_head)

message("loaded.")

message("Loading target file...", appendLF=FALSE)

# Load the target file.
# Get the header of the target file.
target_head = readLines(paste0(data_folder, target_file), n=1)
split_head = unlist(strsplit(target_head, "\t"))
		
# Read in the target expression file.
target_values = fread(paste0(data_folder, target_file), header=F, skip=1, data.table=T)
		
# Again, handle the case when there is no row name for the target
# file header.
if(ncol(target_values) == length(split_head)) {
	target_head = split_head[2:length(split_head)]
} else {
	target_head = split_head
}
rm(split_head)
setnames(target_values, colnames(target_values), c("gene", target_head))
rm(target_head)

message("loaded.")

message("Filtering genes in target to include only those in reference...", appendLF=FALSE)

target_values$gene = trim(target_values$gene)

# Figure out which genes are in the reference file, but not
# in the target file.
missing = ref_values$gene[!(ref_values$gene %in% target_values$gene)]

# Order the genes in the target file to be the same as those in
# the reference file.
#target_values = target_values[match(ref_values$gene, target_values$gene),1:ncol(target_values),with=FALSE]

if(nrow(target_values) < 1) {
	stop("No matching genes found between datasets.")
}

# Filter the genes in the target file to include only those in
# the reference file.
target_values = target_values[target_values$gene %in% ref_values$gene,1:ncol(target_values),with=FALSE]

# Add missing data entries of all 0's for the genes that were in the 
# reference file but not in the target file.
missing_matrix = matrix(rep(0.0, 
		(ncol(target_values)-1) * length(missing)), 
		nrow=length(missing))

if(nrow(missing_matrix) > 0) {
	missing_dt = data.table(cbind(gene = missing, data.frame(missing_matrix)))
	setnames(missing_dt, colnames(missing_dt), colnames(target_values))
	target_values = rbindlist(list(target_values, missing_dt))
}

message("complete.")

message("Checking reference values...", appendLF=FALSE)

# Determine if ref_file contains any negative values.
neg = any(as.vector(as.matrix(ref_values[,2:ncol(ref_values),with=FALSE])) < 0)

# If there are negative values in reference data, then inverse log transform and
# relog transform using x+1.
if(neg) {
	message("negative...transforming...", appendLF=FALSE)
	ref_values = inv_log_transform(ref_values)
	ref_values = log_transform_p1(ref_values)
}

# Order the genes in the target file to be the same as those in
# the reference file.
target_values = target_values[match(ref_values$gene, target_values$gene),1:ncol(target_values),with=FALSE]

message("completed.")

message("Converting NA's to 0's in reference data...", appendLF=FALSE)

# Convert NA to 0.
na_to_zero = function(dt, un = 0) gdata::NAToUnknown(dt, un)	
ref_values = na_to_zero(ref_values)

message("completed.")

message("Performing nonparanormal normalization...", appendLF=FALSE)

pnormal = data.frame(target_values, check.names=FALSE)
rownames(pnormal) = chartr('.', '-', target_values[[1]])
colnames(pnormal) = chartr('.', '-', colnames(target_values))

pnormal = pnormal[,-1]
pnormal = data.matrix(pnormal)

# Get the column names of the target file.
cols = colnames(target_values[,2:ncol(target_values), with=F])

# Create a nonparanormal normalized dataset.
pnormal = huge.npn(t(pnormal), npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)

write.table(t(pnormal), paste0(output_folder, prefix, "_NPN.pcl"), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

message("completed.")

message("Performing zero-to-one scaling of nonparanormal normalized data...", appendLF=FALSE)

pnormal = data.frame(t(pnormal))
rownames(pnormal) = chartr('.', '-', target_values[[1]])
pnormal = cbind(gene=target_values[[1]], pnormal)
pnormal = data.table(pnormal)

# Zero to one scale the paranormal normalized data.
pnormal_zeroone = zero_to_one_transform(pnormal)
pnormal_zeroone = data.frame(pnormal_zeroone, check.names=FALSE)
rownames(pnormal_zeroone) = chartr('.', '-', target_values[[1]])
colnames(pnormal_zeroone) = chartr('.', '-', colnames(target_values))

# Get the column names of the target file.
cols = colnames(target_values[,2:ncol(target_values), with=F])

# Round all entries.
for(j in cols) set(pnormal_zeroone, j=j, value=as.numeric(pnormal_zeroone[[j]]))

write.table(pnormal_zeroone, paste0(output_folder, prefix, "_NPN_ZEROONE.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Performing quantile normalization...", appendLF=FALSE)

# Create a quantile normalized dataset.
# Create a target object for the quantile normalization.
ref_df = data.frame(ref_values[,2:ncol(ref_values),with=FALSE])
rownames(ref_df) = ref_values[[1]]
target_df = data.frame(target_values[,2:ncol(target_values),with=FALSE])
rownames(target_df) = target_values[[1]]

targ = normalize.quantiles.determine.target(
    data.matrix(ref_df),
    target.length=nrow(ref_df))

# Quantile normalize the data, against the reference distribution,
# using replacement - not averaging.
qn = data.matrix(normalize.quantiles.use.target(
                   data.matrix(target_df),targ,copy=F))

qn = data.frame(qn, check.names=FALSE)
rownames(qn) = chartr('.', '-', target_values[[1]])
qn = cbind(gene=target_values[[1]],qn)
colnames(qn) = chartr('.', '-', colnames(target_values))

write.table(qn, paste0(output_folder, prefix, "_QN.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Performing zero-to-one scaling of quantile normalized data...", appendLF=FALSE)

qn = data.table(qn)

# Zero to one scale the quantile normalized data.
qn_zeroone = zero_to_one_transform(qn)
qn_zeroone = data.frame(qn_zeroone, check.names=FALSE)
rownames(qn_zeroone) = chartr('.', '-', target_values[[1]])
colnames(qn_zeroone) = chartr('.', '-', colnames(target_values))

# Get the column names of the target file.
cols = colnames(target_values[,2:ncol(target_values), with=F])

# Round all entries.
for(j in cols) set(qn_zeroone, j=j, value=as.numeric(qn_zeroone[[j]]))


write.table(qn_zeroone, paste0(output_folder, prefix, "_QN_ZEROONE.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Perform TDM normalization...", appendLF=FALSE)

# TDM normalized the data.
tdm_values = tdm_transform(
	target_data=target_values, 
	ref_data=ref_values, 
	negative=FALSE, 
	filter_p=FALSE, 
	inv_reference=TRUE, 
	log_target=TRUE)

tdm_values = data.frame(tdm_values, check.names=FALSE)
rownames(tdm_values) = chartr('.', '-', ref_values[[1]])
colnames(tdm_values) = chartr('.', '-', colnames(target_values))

write.table(tdm_values, paste0(output_folder, prefix, "_TDM.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Performing zero-to-one scaling of TDM normalized data...", appendLF=FALSE)

tdm_values = data.table(tdm_values)

# Zero to one scale the TDM normalized data.
tdm_zeroone = zero_to_one_transform(tdm_values)
tdm_zeroone = data.frame(tdm_zeroone, check.names=FALSE)
rownames(tdm_zeroone) = chartr('.', '-', ref_values[[1]])
colnames(tdm_zeroone) = chartr('.', '-', colnames(target_values))

# Get the column names of the target file.
cols = colnames(target_values[,2:ncol(target_values), with=F])

# Round all entries.
for(j in cols) set(tdm_zeroone, j=j, value=as.numeric(tdm_zeroone[[j]]))

write.table(tdm_zeroone, paste0(output_folder, prefix, "_TDM_ZEROONE.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Performing log transformation...", appendLF=FALSE)

# Log transform the data.
log = log_transform_p1(target_values)

log = data.frame(log, check.names=FALSE)
rownames(log) = chartr('.', '-', target_values[[1]])
colnames(log) = chartr('.', '-', colnames(target_values))
write.table(log, paste0(output_folder, prefix, "_LOG.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Performing zero-to-one scaling of log transformed data...", appendLF=FALSE)

log = data.table(log)

# Zero to one scale the LOG transformed data.
log_zeroone = zero_to_one_transform(log)
log_zeroone = data.frame(log_zeroone, check.names=FALSE)
rownames(log_zeroone) = chartr('.', '-', target_values[[1]])
colnames(log_zeroone) = chartr('.', '-', colnames(target_values))

# Get the column names of the target file.
cols = colnames(target_values[,2:ncol(target_values), with=F])

# Round all entries.
for(j in cols) set(log_zeroone, j=j, value=as.numeric(log_zeroone[[j]]))

write.table(log_zeroone, paste0(output_folder, prefix, "_LOG_ZEROONE.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

if(!is.null(tcga_file)) {
	message("Loading microarray data...", appendLF=FALSE)
	# Load the microarray file.
	# Get the header of the target file.
	tcga_head = readLines(paste0(data_folder, tcga_file), n=1)
	split_head = unlist(strsplit(tcga_head, "\t"))
	
	# Read in the microarray expression file.
	tcga_values = fread(paste0(data_folder, tcga_file), header=F, skip=1, data.table=T)
	
	# Again, handle the case when there is no row name for the 
	# file header.
	if(ncol(tcga_values) == length(split_head)) {
		tcga_head = split_head[2:length(split_head)]
	} else {
		tcga_head = split_head
	}
	rm(split_head)
	setnames(tcga_values, colnames(tcga_values), c("gene", tcga_head))
	rm(tcga_head)

	message("loaded.")

	message("Checking microarray data...", appendLF=FALSE)

	# Determine if ref_file contains any negative values.
	neg = any(as.vector(as.matrix(tcga_values[,2:ncol(tcga_values),with=FALSE])) < 0)

	# If there are negative values in reference data, then inverse log transform and
	# relog transform using x+1.
	if(neg) {
		message("negative...transforming...", appendLF=FALSE)
		tcga_values = inv_log_transform(tcga_values)
		tcga_values = log_transform_p1(tcga_values)
	}

	tcga_values = na_to_zero(tcga_values)

	# Filter the genes in the tcga file to include only those in
	# the reference file.
	tcga_values = tcga_values[tcga_values$gene %in% ref_values$gene,1:ncol(tcga_values),with=FALSE]
	# Order the genes in the tcga file to be the same as those in
	# the reference file.
	tcga_values = tcga_values[match(ref_values$gene, tcga_values$gene),1:ncol(tcga_values),with=FALSE]


	message("completed.")

	message("Performing zero-to-one scaling of microarray data...", appendLF=FALSE)

	# Zero to one scale the microarray data.
	tcga_zeroone = zero_to_one_transform(tcga_values)
	tcga_zeroone = data.frame(tcga_zeroone, check.names=FALSE)
	rownames(tcga_zeroone) = chartr('.', '-', tcga_values[[1]])
	colnames(tcga_zeroone) = chartr('.', '-', colnames(tcga_values))

	# Get the column names of the microarray file.
	cols = colnames(tcga_zeroone[,2:ncol(tcga_zeroone)])

	# Round all entries.
	for(j in cols) set(tcga_zeroone, j=j, value=as.numeric(tcga_zeroone[[j]]))

	write.table(tcga_zeroone, paste0(output_folder, prefix, "_MA_ZEROONE.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

	message("completed.")
}

message("Performing zero-to-one scaling of reference data...", appendLF=FALSE)

# Zero to one transform reference data.
ref_zeroone = zero_to_one_transform(ref_values)
ref_zeroone = data.frame(ref_zeroone, check.names=FALSE)
rownames(ref_zeroone) = chartr('.', '-', ref_values[[1]])
colnames(ref_zeroone) = chartr('.', '-', colnames(ref_values))

# Get the column names of the target file.
cols = colnames(ref_values[,2:ncol(ref_values), with=F])

# Round all entries.
for(j in cols) set(ref_zeroone, j=j, value=as.numeric(ref_zeroone[[j]]))


write.table(ref_zeroone, paste0(output_folder, prefix, "_ZEROONE.pcl"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

message("completed.")

message("Normalization complete.")
