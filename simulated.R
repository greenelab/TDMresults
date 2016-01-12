source('tdm.R')
source('package_loader.R')
source('directories.R')

run_clustering = function(
		num_trials = 20, 
		num_methods = 4, 
		noise_loop = TRUE, 
		noise = 0, 
		zero_data = FALSE, 
		noise_step = .5, 
		ref_file="ref_file.txt", 
		target_file="target_file.txt", 
		nsamples = 400,
		nconditions = 4,
		seed=7532) {
	
	set.seed(seed)

	# Set up a matrix to store the results in.
	#results = matrix(rep(0,num_trials * num_methods), ncol=num_trials)
	ref_cms = tdm_cms = log_cms = qn_cms = npn_cms = list()

	# Define lists to hold clustering stats in.
	ref_stats = tdm_stats = log_stats = qn_stats = npn_stats = list()

	# Define lists to hold correlations in.
	tdm_cor = log_cor = qn_cor = npn_cor = list()

	# Define lists for principal coordinates
	tdm_mds = log_mds = qn_mds = npn_mds = list()

	labels = list()

	# Read in the reference file.
	ref = fread(ref_file, header=T, sep="\t")

	# The reference file will be changed each time through the loop, so keep 
	# the original as well.
	orig_ref = ref

	ref_t = t(ref)
	
	# Rescale the data so that it goes from 4 to 14. The simulator
	# produces values from [0,1].
	orig_ref_scaled = data.frame(
		matrix(scales::rescale(
			as.vector(data.matrix(ref_t)), to=c(4,14)), ncol=nsamples))

	# Read in the target file to a data.table. This is the file we want
	# to normalize.
	targ = fread(target_file, header=F, sep="\t")
	
	orig_targ = targ	
	
	# Transpose the file so that genes are rows.
	targ_t = t(targ)
	
	# Rescale the data so that they go from [4,14].
	orig_targ_scaled = data.frame(matrix(scales::rescale(
					as.numeric(as.vector(
						data.matrix(
							targ_t[,2:ncol(targ_t)]))), to=c(4,14)), ncol=nsamples))

	# Now perform num_trials trials, possibly adding noise for each.
	for(i in 0:(num_trials-1)) {	
		# Get a permutation of the columns from the reference data.
		# They come from the simulator grouped together by condition,
		# so we want to break that up.
		orig_cols = sample(ncol(orig_ref_scaled))

		# Each time through this loop, reset the data to their
		# untransformed values.
		ref_scaled = orig_ref_scaled
		targ_scaled = orig_targ_scaled
		ref = orig_ref
		targ = orig_targ

		# The labels for the data are based on their original order,
		# i.e. their value in orig_cols. The first 100 are condition 1,
		# the next 100 are condition 2, etc.
		ref_labels = orig_cols

		# Divide the labels into nconditions groups.
		ref_labels = factor(as.numeric(cut(ref_labels, nconditions)))
		
		# Create a vector of labels for the data. The first label
		# to be encountered is 1, and so forth.
		ref_labels = as.numeric(factor(ref_labels,levels=unique(ref_labels)))
		
		# Now permute the reference data.
		cols = orig_cols
		ref_scaled = ref_scaled[,cols]
		
		# Name the samples in the reference data.
		colnames(ref_scaled) = paste0("JIHGFEDCBA", 1:ncol(ref_scaled))
		
		# The file with the reference data has samples represented as rows,
		# and genes as columns, so this creates a column with gene names.
		ref_values = data.table(cbind(colnames(ref), ref_scaled))
		
		# Set the name of the column with the genes.
		setnames(ref_values, colnames(ref_values)[1], "gene")

		# Permute the target file. We do this in the same order as 
		# the reference file so that we can keep track of which samples
		# are which. 
		targ_scaled = targ_scaled[,cols]
		targ_scaled = data.table(targ_scaled)
		
		# Add a column with the gene names.
		targ_scaled = cbind(as.character(targ[1,]), targ_scaled)
		
		# Inverse log transform the data
		targ_scaled = inv_log_transform(targ_scaled)
		
		# Rescale the data again, this time to an increased dynamic range
		# to simulate RNA-seq (the range, not so much the distribution).
		# Combined with additional noise, we'll get some of the kind of 
		# changes we see in real RNA-seq data.
		targ_scaled = round(data.frame(matrix(scales::rescale(
							as.numeric(as.vector(data.matrix(targ_scaled[,2:ncol(targ_scaled),with=F]))), 
							to=c(0,1000000)), ncol=nsamples)))
		targ_orig = data.table(targ_scaled)
		# Add noise to the target data. If NOISE_LOOP, then
		# the noise is proportional to the trial.
		if(noise_loop) {
			res = addNoise(targ_scaled, colnames(targ_scaled), i*noise_step,"additive")
			targ_scaled = res$xm
			targ_scaled = targ_scaled + abs(min(targ_scaled))
		} else if(noise > 0) {
			res = addNoise(targ_scaled, colnames(targ_scaled), noise,"additive")
			targ_scaled = res$xm
			targ_scaled = targ_scaled + abs(min(targ_scaled))
		}
		ref_values = data.table(cbind(colnames(ref), ref_scaled))
		
		# Set the name of the column with the genes.
		setnames(ref_values, colnames(ref_values)[1], "gene")
		
		# If ZERO_DATA, then zero out a number of genes proprotional
		# to the trial.
		if(zero_data) {
			targ_scaled[sample(nrow(targ_scaled), i/100*nrow(targ_scaled)),]=0
		}
		
		# Now add the gene names to the target data.
		targ_scaled = cbind(as.character(targ[1,]),targ_scaled)
		colnames(targ_scaled)[2:ncol(targ_scaled)] = paste0("ABCDEFGHIJ", 1:(ncol(targ_scaled)-1))
		colnames(targ_scaled)[1] = "gene"
		targ_scaled = data.table(targ_scaled)

		# Perform TDM transformation on the RNA-seq data.
		tdm_dt = tdm_transform(
			target_data=targ_scaled,
			ref_data=ref_values, 
			negative=TRUE, 
			filter_p=FALSE, 
			inv_reference=FALSE, 
			log_target=FALSE)
		
		# Filter the target data to include only those genes, and in
		# that order, found in the reference distribution.
		tdm_dt = tdm_dt[gene %in% ref_values$gene]    # gene filter
		tdm_dt = tdm_dt[match(ref_values$gene, gene)] # reorder
		
		# Convert reference values to be numeric instead of string.
		ref_num = apply(ref_values[,2:ncol(ref_values),with=F], 2, function(x) as.numeric(x))
		rownames(ref_num) = ref_values[[1]]
		
		targ_num = apply(targ_orig[,1:ncol(targ_orig),with=F], 2, function(x) as.numeric(x))
		rownames(targ_num) = targ_scaled[[1]]

		# Convert target data to be numeric instead of string.
		tdm_num = apply(tdm_dt[,2:ncol(tdm_dt),with=F], 2, function(x) as.numeric(x))
		rownames(tdm_num) = tdm_dt[[1]]
		
		# Create an additional target dataset that is log-transformed, but not TDM normalized.
		# We can use this as the 'naive' approach.
		log_dt = log_transform_p1(targ_scaled)
		setnames(log_dt, colnames(log_dt)[1], "gene")
		
		# Put genes in the order found in the reference distribution.
		log_dt = log_dt[match(ref_values$gene, gene)]
		
		# Make sure that the new target data and reference files are in the same order
		# and have the same genes.
		log_dt = log_dt[gene %in% ref_values$gene]
		log_dt = log_dt[match(ref_values$gene, gene)]
		
		# Convert the data to numeric.
		log_num = apply(log_dt[,2:ncol(log_dt),with=F], 2, function(x) as.numeric(x))
		rownames(log_num) = log_dt[[1]]
		
		# Create yet another target dataset that is log-transformed, but not TDM normalized.
		qn_dt = log_transform_p1(targ_scaled)
		
		# Create a target object for targeted quantile normalization.
		qn_targ = normalize.quantiles.determine.target(
				data.matrix(ref_values[,2:ncol(ref_values),with=F]), 
				target.length=nrow(ref_values))
		
		# Quantile normalize the data, against the reference distribution,
		# using replacement - not averaging. 
		qn_mt = normalize.quantiles.use.target(
				data.matrix(qn_dt[,2:ncol(qn_dt),with=F]),qn_targ,copy=T)
		
		# Add the gene names and turn back into a data.table.
		# It's unecessary to log transform or scale the data
		# because its distribution is now an exact match for
		# the reference data.
		qn_mt = cbind(qn_dt[[1]], qn_mt)
		qn_mt = data.table(qn_mt)
		setnames(qn_mt, colnames(qn_mt), colnames(qn_dt))
		
		# Put genes in the order found in the reference distribution.
		qn_mt = qn_mt[match(ref_values$gene, gene)]
		
		# Make sure that the RNA-seq and reference files are in the same order
		# and have the same genes.
		qn_mt = qn_mt[gene %in% ref_values$gene]
		qn_mt = qn_mt[match(ref_values$gene, gene)]
				
		# Convert the data to numeric.
		qn_num = apply(qn_mt[,2:ncol(qn_mt),with=F], 2, function(x) as.numeric(x))
		rownames(qn_num) = qn_mt[[1]]

		# Create an additional target dataset that is nonparanormal.
		npn_df = data.frame(targ_scaled)
		rownames(npn_df) = targ_scaled[[1]]
		npn_df = npn_df[,-1]
		npn_mt = t(data.matrix(npn_df))

		npn_mt = suppressMessages(huge.npn(npn_mt, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE))
		npn_mt = t(npn_mt)

		# Put genes in the order found in the reference distribution.
		npn_mt = npn_mt[match(ref_values$gene, rownames(npn_mt)),]
		
		# Make sure that the new target data and reference files are in the same order
		# and have the same genes.
		npn_mt = npn_mt[rownames(npn_mt) %in% ref_values$gene,]
		npn_mt = npn_mt[match(ref_values$gene, rownames(npn_mt)),]
		
		# Convert the data to numeric.
		npn_num = apply(npn_mt, 2, function(x) as.numeric(x))

		npn_ref_df = data.frame(ref_values)
		rownames(npn_ref_df) = ref_values[[1]]
		npn_ref_df = npn_ref_df[,-1]
		npn_ref_mt = t(data.matrix(npn_ref_df))

		npn_ref_mt = suppressMessages(huge.npn(npn_ref_mt, npn.func = "shrinkage", npn.thresh=NULL, verbose=TRUE))
		npn_ref_mt = t(npn_ref_mt)

		npn_ref_num = apply(npn_ref_mt, 2, function(x) as.numeric(x))

		# [0,1] transform everything
        ref_num = apply(ref_num, 1, function(x){scales::rescale(x, to=c(0,1))})
        targ_num = apply(targ_num, 1, function(x){scales::rescale(x, to=c(0,1))})
        tdm_num = apply(tdm_num, 1, function(x){scales::rescale(x, to=c(0,1))})
        log_num = apply(log_num, 1, function(x){scales::rescale(x, to=c(0,1))})
        qn_num = apply(qn_num, 1, function(x){scales::rescale(x, to=c(0,1))})
        npn_num = apply(npn_num, 1, function(x){scales::rescale(x, to=c(0,1))})
        npn_ref_num = apply(npn_ref_num, 1, function(x){scales::rescale(x, to=c(0,1))})

        tdm_cor[[i+1]] = sapply(1:nrow(tdm_num), function(x) cor.test(tdm_num[x,], ref_num[x,], method="kendall"))
        log_cor[[i+1]] = sapply(1:nrow(log_num), function(x) cor.test(log_num[x,], ref_num[x,], method="kendall"))
        qn_cor[[i+1]] = sapply(1:nrow(qn_num), function(x) cor.test(qn_num[x,], ref_num[x,], method="kendall"))
        npn_cor[[i+1]] = sapply(1:nrow(npn_num), function(x) cor.test(npn_num[x,], npn_ref_num[x,], method="kendall"))
        
        colnames(npn_num) = colnames(log_num)
        rownames(npn_num) = rownames(log_num)
        colnames(npn_ref_num) = colnames(ref_num)
        rownames(npn_ref_num) = rownames(ref_num)

		# PAM is like kmeans, but it seems to work a bit better.
		ref_p = pam(ref_num, k = nconditions)
		ref_dist = dist(ref_num, method="euclidean")
		
		npn_ref_p = pam(npn_ref_num, k=nconditions)
		
		# Get some stats on the clusters.
		ref_stats[[i+1]] = cluster.stats(ref_dist, as.numeric(ref_p$clustering))

		# Get the distance between all data points.
		tdm_dist = dist(tdm_num, method="euclidean")
		
		# Predict clusters for TDM normalized data, given reference clusters and
		# generate stats for them.
		tdm_pred = as.numeric(predict(as.kcca(ref_p, ref_num), newdata=tdm_num))
		tdm_stats[[i+1]] = get_cluster_stats(tdm_dist, tdm_pred)

		log_dist = dist(log_num, method="euclidean")
		log_pred = as.numeric(predict(as.kcca(ref_p, ref_num), newdata=log_num))
		log_stats[[i+1]] = get_cluster_stats(log_dist, log_pred)

		qn_dist = dist(qn_num, method="euclidean")
		qn_pred = as.numeric(predict(as.kcca(ref_p, ref_num), newdata=qn_num))
		qn_stats[[i+1]] = get_cluster_stats(qn_dist, qn_pred)

		npn_dist = dist(npn_num, method="euclidean")
		npn_pred = as.numeric(predict(as.kcca(npn_ref_p, npn_ref_num), newdata=npn_num))
		npn_stats[[i+1]] = get_cluster_stats(npn_dist, npn_pred)

		# Store the results for this trial.
		tdm_cms[[i+1]] = confusionMatrix(tdm_pred, ref_p$clustering)$overall
		log_cms[[i+1]] = confusionMatrix(log_pred, ref_p$clustering)$overall
		qn_cms[[i+1]] = confusionMatrix(qn_pred, ref_p$clustering)$overall
		npn_cms[[i+1]] = confusionMatrix(npn_pred, npn_ref_p$clustering)$overall
		
		tdm_mds[[i+1]] = plotMDS(t(tdm_num), labels=ref_labels, pch=8, dim.plot=c(1,2), col=c(ref_labels), gene.selection="pairwise", ndim=45)
		log_mds[[i+1]] = plotMDS(t(log_num), labels=ref_labels, pch=8, dim.plot=c(1,2), col=c(ref_labels), gene.selection="pairwise", ndim=45)
		qn_mds[[i+1]] = plotMDS(t(qn_num), labels=ref_labels, pch=8, dim.plot=c(1,2), col=c(ref_labels), gene.selection="pairwise", ndim=45)
		npn_mds[[i+1]] = plotMDS(t(npn_num), labels=ref_labels, pch=8, dim.plot=c(1,2), col=c(ref_labels), gene.selection="pairwise", ndim=45)

		labels[[i+1]] = ref_labels

		message("Step: ", i, appendLF=TRUE)
	}
	
	return(list(last_ref=ref_num, last_tdm=tdm_num, last_log=log_num, last_qn=qn_num, last_npn=npn_num, labels=labels, ref_stats=ref_stats, tdm_stats=tdm_stats, log_stats=log_stats, qn_stats=qn_stats, npn_stats=npn_stats, tdm_cor=tdm_cor, log_cor=log_cor, qn_cor=qn_cor, npn_cor=npn_cor, tdm_cms=tdm_cms, log_cms=log_cms, qn_cms=qn_cms, npn_cms=npn_cms, tdm_mds=tdm_mds, log_mds=log_mds, qn_mds=qn_mds, npn_mds=npn_mds))
}

get_cluster_stats = function(dist, pred) {
	if(length(unique(pred)) == 1) {
		return(NULL)
	} else {
		return(cluster.stats(dist, pred))
	}
} # end get_cluster_stats

# This function (from package_loader.R) will check all of the listed packages 
# (they don't all need to be done before the rest of the program) to see 
# if they are loaded. If not, it will load them if they are available, 
# else it will install them and then load them.
load_it(c("ggplot2",
		"reshape2",
		"Hmisc",
		"data.table",
		"scales",
		"sdcMicro",
		"flexclust",
		"fpc",
		"corrplot",
		"ape",
		"cluster",
		"plyr",
		"dplyr",
		"devtools",
		"quantro",
		"preprocessCore",
		"gridExtra",
		"huge",
		"caret",
		"limma"))

# This one is on github, so we will handle it separately.
if(!require("quantroSim", quietly=T)){
	install_github(repo = "quantroSim", username = "stephaniehicks")
	library(quantroSim, quiet=T)
}

ref_file = paste0(initial_data,
						"Rnn250_nbgr250_hop0.3_bionoise0.1_expnoise0.1",
						"_corrnoise0.1_clustAdd_dataset.txt")
target_file = ref_file

set.seed(2357)

NUM_TRIALS = 20 # How many runs do we want?
NUM_METHODS = 4 # How many competing methods are we trying?

NOISE_LOOP = TRUE # Should we loop while increasing noise?
NOISE = 0		  # How much noise?
ZERO_DATA = FALSE # Should we zero out some genes?

# Run the tests.
test_results = run_clustering(
		num_trials=NUM_TRIALS, 
		num_methods=NUM_METHODS, 
		noise_loop=NOISE_LOOP, 
		noise=NOISE, 
		zero_data=ZERO_DATA, 
		noise_step=.2,
		ref_file=ref_file,
		target_file=target_file)

# Get the correct labels for each run.
labels = test_results$labels

# Put the confusion matrices in data.frames.
tdm_results = do.call(rbind, test_results$tdm_cms)
log_results = do.call(rbind, test_results$log_cms)
qn_results = do.call(rbind, test_results$qn_cms)
npn_results = do.call(rbind, test_results$npn_cms)

# Make the data.frames numeric.
tdm_results = data.frame(apply(tdm_results, 2, as.numeric))
log_results = data.frame(apply(log_results, 2, as.numeric))
qn_results = data.frame(apply(qn_results, 2, as.numeric))
npn_results = data.frame(apply(npn_results, 2, as.numeric))

# Tag on the method name.
tdm_results = cbind(tdm_results, Method="TDM")
log_results = cbind(log_results, Method="LOG")
qn_results = cbind(qn_results, Method="QN")
npn_results = cbind(npn_results, Method="NPN")

# Tag on the noise levels.
tdm_results = cbind(tdm_results, Noise=1:20)
log_results = cbind(log_results, Noise=1:20)
qn_results = cbind(qn_results, Noise=1:20)
npn_results = cbind(npn_results, Noise=1:20)

# Put all the results in one data.frame.
results = rbind(tdm_results, log_results, qn_results, npn_results)

# Change the noise levels to the percent noise.
results$Noise = (results$Noise - 1) * .2

# tdm_silwidths = sapply(test_results$tdm_stats, function(x) x$avg.silwidth)
# log_silwidths = sapply(test_results$log_stats, function(x) x$avg.silwidth)
# qn_silwidths = sapply(test_results$qn_stats, function(x) x$avg.silwidth)
# npn_silwidths = sapply(test_results$npn_stats, function(x) x$avg.silwidth)

# tdm_sws_df = data.frame(AvgSilwidths=tdm_silwidths, Method="TDM")
# tdm_sws_df = cbind(tdm_sws_df, Noise=1:20)
# log_sws_df = data.frame(AvgSilwidths=log_silwidths, Method="LOG")
# log_sws_df = cbind(log_sws_df, Noise=1:20)
# qn_sws_df = data.frame(AvgSilwidths=qn_silwidths, Method="QN")
# qn_sws_df = cbind(qn_sws_df, Noise=1:20)
# npn_sws_df = data.frame(AvgSilwidths=npn_silwidths, Method="NPN")
# npn_sws_df = cbind(npn_sws_df, Noise=1:20)

# silwidths = rbind(tdm_sws_df, log_sws_df, qn_sws_df, npn_sws_df)

results$Method = factor(results$Method, levels=c("TDM", "QN", "LOG", "NPN"))
# Show the results.
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                    "#CC79A7", "#F0E442")
ggplot(data=results, aes(x=Noise, y=Accuracy, group=Method, color=Method, shape=Method)) +
		geom_line() +
		geom_point() + ylab("Proportion correctly classified") + 
		xlab("Percent noise added") +
		#geom_rect(aes(xmin=1.9, xmax=2.1, ymin=.74, ymax=.76),color="red") +
		#geom_rect(aes(xmin=13.9, xmax=14.1, ymin=.27, ymax=.29),color="red") +
		#geom_rect(aes(xmin=6.9, xmax=7.1, ymin=.32, ymax=.34),color="red") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1), 
				axis.title.y = element_text(size=10),
				axis.title.x = element_text(size=10)) +
		expand_limits(x = 0, y = 0) +
		theme(legend.position="bottom", legend.direction="horizontal") +
		ylim(.2,1) +
		scale_color_manual(values=cbbPalette)

ggsave(paste0(output, "simulated.pdf"), plot=last_plot(), width=6, height=3)

mds_log = data.frame(MDS1=test_results$log_mds[[10]]$x, MDS2=test_results$log_mds[[10]]$y, Method="LOG")
mds_tdm = data.frame(MDS1=test_results$tdm_mds[[10]]$x, MDS2=test_results$tdm_mds[[10]]$y, Method="TDM")
mds_qn = data.frame(MDS1=test_results$qn_mds[[10]]$x, MDS2=test_results$qn_mds[[10]]$y, Method="QN")
mds_npn = data.frame(MDS1=test_results$npn_mds[[10]]$x, MDS2=test_results$npn_mds[[10]]$y, Method="NPN")

mds_log = cbind(mds_log, labels = as.factor(test_results$labels[[10]]))
mds_tdm = cbind(mds_tdm, labels = as.factor(test_results$labels[[10]]))
mds_qn = cbind(mds_qn, labels = as.factor(test_results$labels[[10]]))
mds_npn = cbind(mds_npn, labels = as.factor(test_results$labels[[10]]))

mds_all = rbind(mds_tdm, mds_log, mds_qn, mds_npn)

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                    "#CC79A7", "#F0E442")

ggplot() + 
		geom_point(data=mds_all,aes(x=MDS1,y=MDS2,shape=labels,colour=labels),size=2) + 
		theme_bw() + 
		theme(axis.text.x = element_text(size=5),  
				axis.text.y = element_text(size=5), 
				axis.title.x = element_text(size=10), 
				axis.title.y = element_text(size=10), 
				panel.background = element_blank(), 
				panel.grid.major = element_blank(),  
				panel.grid.minor = element_blank(),  
				plot.background = element_blank(),
				plot.title = element_text(size=12)) +
		labs(x="Coordinate 1", y="Coordinate 2") +
		theme(legend.direction = "horizontal", 
				legend.position = "bottom", 
				legend.box="vertical") +
		scale_color_manual(values=cbbPalette) +
		facet_wrap(~Method, ncol=4)

ggsave(paste0(output, "mds_plot.pdf"), plot=last_plot(), height=3, width=6)

tdm_correlation = data.frame(Tau=sapply(test_results$tdm_cor, `[[`, 4), Dataset="TDM", Noise=1:20)
log_correlation = data.frame(Tau=sapply(test_results$log_cor, `[[`, 4), Dataset="LOG", Noise=1:20)
qn_correlation = data.frame(Tau=sapply(test_results$qn_cor, `[[`, 4), Dataset="QN", Noise=1:20)
npn_correlation = data.frame(Tau=sapply(test_results$npn_cor, `[[`, 4), Dataset="NPN", Noise=1:20)

correlation = rbind(tdm_correlation, log_correlation, qn_correlation, npn_correlation)

correlation$Noise = (correlation$Noise - 1) * .2

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                    "#CC79A7", "#F0E442")
ggplot(data=correlation, aes(x=Noise, y=Tau, color=Dataset, group=Dataset, shape=Dataset)) + 
	geom_line() +
	geom_point() +
	ylab("Kendall's Tau") +
	xlab("Noise") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1), 
			axis.title.y = element_text(size=10),
			axis.title.x = element_text(size=10)) +
	expand_limits(x = 0, y = 0) +
	theme(legend.position="bottom", legend.direction="horizontal") +
	scale_color_manual(values=cbbPalette)
ggsave(paste0(output, "correlations.pdf"), plot=last_plot(), width=6, height=3)

quantro(t(test_results$last_log), test_results$labels[[20]], B=100)

plot_mds = function(data) {
	the_plot = ggplot() + 
		geom_point(data=data,aes(x=MDS1,y=MDS2,shape=labels,colour=labels),size=2) + 
		theme_bw() + 
		theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),  
				axis.text.y = element_text(size=8, angle = 90, hjust = 1), 
				axis.title.x = element_text(size=10), 
				axis.title.y = element_text(size=10), 
				panel.background = element_blank(), 
				panel.grid.major = element_blank(),  
				panel.grid.minor = element_blank(),  
				plot.background = element_blank(),
				plot.title = element_text(size=12)) +
		labs(x="Coordinate 1", y="Coordinate 2") +
		theme(legend.direction = "horizontal", 
				legend.position = "bottom", 
				legend.box="vertical") +
		scale_color_manual(values=cbbPalette) 
		return(the_plot)
}

tdm_plot = plot_mds(mds_tdm)
log_plot = plot_mds(mds_log)
qn_plot = plot_mds(mds_qn)
npn_plot = plot_mds(mds_npn)

output = plot_grid(tdm_plot, log_plot, qn_plot, npn_plot, ncol = 4,
          align = "h", labels = c("A", "B", "C", "D"), label_size = 15)
save_plot("mds_plots.pdf", output, ncol = 4, base_width=2.5, base_height=3)