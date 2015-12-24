# Author: Jeffrey Thompson
# Purpose: Performs multiclass lasso logistic regression to classify cancer subtypes.

# Get command line arguments.
#args <- commandArgs(trailingOnly = T)

source("package_loader.R")
source("directories.R")

NFOLDS = 100
NSEEDS = 10
INITIAL_SEED = 2
ACCURACY_PLOT_FILENAME = "meta_accuracies.pdf"
KAPPA_PLOT_FILENAME = "meta_kappas.pdf"
ACCURACY_TABLE_FILENAME = "meta_accuracies.tsv"
BALANCED_PLOT_FILENAME = "meta_balanced_accuracies.pdf"

args = c(normalized_data, 
		'META_BRCA_ZEROONE.pcl', 
		'combined_subtype.txt', 
		'META_BRCA_MA_ZEROONE.pcl', 
		'BRCAClin.tsv', 
		'META_BRCA_LOG_ZEROONE.pcl', 
		'BRCARNASeqClin.tsv', 
		'META_BRCA_QN_ZEROONE.pcl', 
		'BRCARNASeqClin.tsv', 
		'META_BRCA_TDM_ZEROONE.pcl', 
		'BRCARNASeqClin.tsv', 
		'META_BRCA_NPN_ZEROONE.pcl',
    		'BRCARNASeqClin.tsv',
		'META_BRCA_REF_NPN_ZEROONE.pcl',
		'combined_subtype.txt',
		'META_BRCA_UN_ZEROONE.pcl',
		'BRCARNASeqClin.tsv',
		output)

input_dir = args[1]
ref_input = args[2]
ref_clin = args[3]
tcga_input = args[4]
tcga_clin = args[5]
log_input = args[6]
log_clin = args[7]
qn_input = args[8]
qn_clin = args[9]
tdm_input = args[10]
tdm_clin = args[11]
npn_input = args[12]
npn_clin = args[13]
ref_npn_input = args[14]
ref_npn_clin = args[15]
un_input = args[16]
un_clin = args[17]
output_dir = args[18]

load_it(c("glmnet", "caret", "e1071", "stringr", "plyr", "huge"))

# This function creates dataframes to hold the accuracies for each subtype across runs.
setup_acc_df = function(df) {
  df = data.frame(matrix(nrow=5, ncol=0))
  rownames(df) = c("Basal", "Her2", "LumA", "LumB", "Normal")
  return (df)
}

# Define data.frames to hold classification results.
metadf = setup_acc_df(metadf)
tcgadf = setup_acc_df(tcgadf)
logdf = setup_acc_df(logdf)
qndf = setup_acc_df(qndf)
tdmdf = setup_acc_df(tdmdf)
npndf = setup_acc_df(npndf)
undf = setup_acc_df(undf)
tdmnotshareddf = setup_acc_df(tdmnotshareddf)
qnnotshareddf = setup_acc_df(qnnotshareddf)
lognotshareddf = setup_acc_df(lognotshareddf)
npnnotshareddf = setup_acc_df(npnnotshareddf)
unnotshareddf = setup_acc_df(unnotshareddf)

# Define random seeds.
set.seed(INITIAL_SEED)
seeds = sample(1:10000, NSEEDS)

message(paste("Random seeds:", paste(seeds, collapse=", ")), appendLF=TRUE)

# Pre-processes test data to ready it for the LASSO.
# Returns a list with first entry the pre-processed dataset and second entry
# the classes for those data.
preprocessTCGA = function(dataFile, clinFile, RNAseq=TRUE) {
  data = read.table(paste0(input_dir, dataFile), 
    sep='\t', 
    row.names=1, 
    header=T)

  dataT = t(data)
  
  # Remove duplicate samples from the data.
  dataT = dataT[str_detect(substr(rownames(dataT),1,15), "(\\.01)|(\\.11)"),]  

  # Read clinical data.
  clin = read.table(paste0(input_dir, clinFile), sep='\t', row.names=1, header=T)

  if(!RNAseq) {
    clin[,1] = clin[,3]
  }

  # Remove samples without subtype labels.
  clin = subset(clin, clin[,1] != "")
  clin[,1] = droplevels(clin[,1])

  # Filter clinical data to include only those samples with expression data.
  clin.response = clin[,1][match(substr(chartr('-', '.', rownames(dataT)), 1, 15), substr(chartr('-', '.', rownames(clin)),1,15))]
  clin.response = na.omit(clin.response)

  # Filter expression data to include only those samples with ids also in clinical data.
  dataT = dataT[substr(chartr('-', '.', rownames(dataT)), 1, 15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]
  
  # Filter any data with missing clinical annotation for tumor subtype.
  dataT = dataT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(dataT)),1,15), 
    substr(chartr('-', '.', rownames(clin)),1,15))]),]
  
  return (list(data=dataT, response=clin.response))
} # end function preprocessTCGA

# Predict subtypes on TCGA data using previously trained model.
predictTCGA = function(title, data, model) { 
  # Predict classes based on the trained model.
  lasso.predict = predict(model, data$data, s = "lambda.1se", type="class")
  
  # Make sure all factors are included.
  lasso.predict = factor(lasso.predict, c("Basal", "Her2", "LumA", "LumB", "Normal"))
  
  # Build a contingency table of the results.
  con_table = table(pred = lasso.predict, true = data$response, exclude="none")
  
  # Generate statistics based on the contigency table.
  cm = confusionMatrix(con_table)
  
  return (cm)
}

# Load the training data.
train = read.table(paste0(input_dir, ref_input), sep='\t', row.names=1, header=T)

# Transpose rows and columns.
trainT = t(train)

# Load the clinical features.
clin = read.table(paste0(input_dir, ref_clin), sep='\t', row.names=1, header=T)

# Remove annotations that have no subtype label.
clin = subset(clin, clin[,1] != "NC")
clin[,1] = droplevels(clin[,1])

# Filter the expression data to remove those samples without subtype labels.
trainT = trainT[rownames(trainT) %in% chartr('-', '.', rownames(clin)),]

# Filter clinical data to include only those entries that also have expression data.
trainclin.subtypes = clin[,1][match(chartr('-', '.', rownames(trainT)), chartr('-', '.', rownames(clin)))]
trainclin.subtypes = na.omit(trainclin.subtypes) 
write.table(trainclin.subtypes, file=paste0(output_dir, "meta_ref_classes.txt"))

# Filter any expression data with missing clinical annotation for tumor subtype.
trainT = trainT[!is.na(clin[,1][match(chartr('-', '.', rownames(trainT)), chartr('-', '.', rownames(clin)))]),]

# Load the training data.
npn_train = read.table(paste0(input_dir, ref_npn_input), sep='\t', row.names=1, header=T)

# Transpose rows and columns.
npntrainT = t(npn_train)

# Filter the expression data to remove those samples without subtype labels.
npntrainT = npntrainT[rownames(npntrainT) %in% chartr('-', '.', rownames(clin)),]

# Filter any expression data with missing clinical annotation for tumor subtype.
npntrainT = npntrainT[!is.na(clin[,1][match(chartr('-', '.', rownames(npntrainT)), chartr('-', '.', rownames(clin)))]),]

# Pre-process the TCGA data.
tcga = preprocessTCGA(tcga_input, tcga_clin, FALSE)
log = preprocessTCGA(log_input, log_clin)
qn = preprocessTCGA(qn_input, qn_clin)
tdm = preprocessTCGA(tdm_input, tdm_clin)
npn = preprocessTCGA(npn_input, npn_clin)
un = preprocessTCGA(un_input, un_clin)

# Determine which samples are shared between TCGA RNA-seq and TCGA microarray data.
# Then filter the data to include only those that are shared.

tdm_shared = substr(chartr('-', '.', rownames(tdm$data)),1,15) %in% substr(chartr('-', '.', rownames(tcga$data)),1,15) 
tdm$data = tdm$data[tdm_shared,]
tdm$response = tdm$response[tdm_shared]

tcga_shared = substr(chartr('-', '.', rownames(tcga$data)),1,15) %in% substr(chartr('-', '.', rownames(tdm$data)),1,15) 
tcga$data = tcga$data[tcga_shared,]
tcga$response = tcga$response[tcga_shared]

qn_shared = substr(chartr('-', '.', rownames(qn$data)),1,15) %in% substr(chartr('-', '.', rownames(tcga$data)),1,15)
qn$data = qn$data[qn_shared,]
qn$response = qn$response[qn_shared]

log_shared = substr(chartr('-', '.', rownames(log$data)),1,15) %in% substr(chartr('-', '.', rownames(tcga$data)),1,15)
log$data = log$data[log_shared,]
log$response = log$response[log_shared]

npn_shared = substr(chartr('-','.', rownames(npn$data)),1,15) %in% substr(chartr('-','.',rownames(tcga$data)),1,15)
npn$data = npn$data[npn_shared,]
npn$response = npn$response[npn_shared]

un_shared = substr(chartr('-', '.', rownames(un$data)),1,15) %in% substr(chartr('-','.',rownames(tcga$data)),1,15)
un$data = un$data[un_shared,]
un$response = un$response[un_shared]

# Pre-allocate vectors to hold results.
tcga_accs = vector(mode="list", length=length(seeds))
tdm_accs = vector(mode="list", length=length(seeds))
qn_accs = vector(mode="list", length=length(seeds))
log_accs = vector(mode="list", length=length(seeds))
npn_accs = vector(mode="list", length=length(seeds))
un_accs = vector(mode="list", length=length(seeds))

# Perform a number of iterations equal to the number of random seeds.
# At each iteration, perform n-fold cross validation to build a model on the training data,
# then use that model to make predictions on the test data.
for(seednum in 1:length(seeds)) {
  set.seed(seeds[seednum])
  
  # Train a model for classifying cancer subtype on the METABRIC training data.
  set.seed(seeds[seednum])
  lasso.model=cv.glmnet(trainT, trainclin.subtypes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
  plot(lasso.model)

  set.seed(seeds[seednum])
  npn.lasso.model=cv.glmnet(npntrainT, trainclin.subtypes, family="multinomial", parallel=F, type.measure="class", nfolds=NFOLDS)

  acc = predictTCGA("TCGA Results", tcga, lasso.model)
  tcga_accs[[seednum]] = acc
  tcgadf = cbind(tcgadf, as.vector(acc$byClass[,8]))
      
  acc <- predictTCGA("TDM Results", tdm, lasso.model)
  tdm_accs[[seednum]] = acc
  tdmdf <- cbind(tdmdf,as.vector(acc$byClass[,8]))
  
  acc <- predictTCGA("QN Results", qn, lasso.model)
  qn_accs[[seednum]] = acc
  qndf <- cbind(qndf,as.vector(acc$byClass[,8]))

  acc <- predictTCGA("LOG Results", log, lasso.model)
  log_accs[[seednum]] = acc
  logdf <- cbind(logdf,as.vector(acc$byClass[,8]))

  acc = predictTCGA("NPN Results", npn, npn.lasso.model)
  npn_accs[[seednum]] = acc
  npndf = cbind(npndf,as.vector(acc$byClass[,8]))

  acc = predictTCGA("Untransformed Results", un, lasso.model)
  un_accs[[seednum]] = acc
  undf = cbind(undf,as.vector(acc$byClass[,8]))
}

# Build a table of accuracies across all datasets.
accuracies = data.frame(Basal=numeric(0), Her2=numeric(0), LumA=numeric(0), LumB=numeric(0), Normal=numeric(0))
accuracies = rbind(accuracies, apply(tcgadf, 1, mean))
accuracies = rbind(accuracies, apply(tdmdf,1,mean))
accuracies = rbind(accuracies, apply(qndf,1,mean))
accuracies = rbind(accuracies, apply(logdf,1,mean))
accuracies = rbind(accuracies, apply(npndf,1,mean))
accuracies = rbind(accuracies, apply(undf,1,mean))
rownames(accuracies) <- c("MA", "TDM", "QN", "LOG", "NPN", "UNTR")
colnames(accuracies) <- c("Basal", "Her2", "LumA", "LumB", "Normal")


# Aggregate accuracies:
tdm_tables = lapply(tdm_accs, function(x) x$table)
qn_tables = lapply(qn_accs, function(x) x$table)
log_tables = lapply(log_accs, function(x) x$table)
tcga_tables = lapply(tcga_accs, function(x) x$table)
npn_tables = lapply(npn_accs, function(x) x$table)
un_tables = lapply(un_accs, function(x) x$table)

tdm_reduced = Reduce("+", tdm_tables) / length(tdm_tables)
qn_reduced = Reduce("+", qn_tables) / length(qn_tables)
log_reduced = Reduce("+", log_tables) / length(log_tables)
tcga_reduced = Reduce("+", tcga_tables) / length(tcga_tables)
npn_reduced = Reduce("+", npn_tables) / length(npn_tables)
un_reduced = Reduce("+", un_tables) / length(un_tables)

tdm_reduced[1:5,1:5] = t(apply(tdm_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
qn_reduced[1:5,1:5] = t(apply(qn_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
log_reduced[1:5,1:5] = t(apply(log_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
tcga_reduced[1:5,1:5] = t(apply(tcga_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
npn_reduced[1:5,1:5] = t(apply(npn_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
un_reduced[1:5,1:5] = t(apply(un_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]

tdm_cm = confusionMatrix(tdm_reduced)
qn_cm = confusionMatrix(qn_reduced)
log_cm = confusionMatrix(log_reduced)
tcga_cm = confusionMatrix(tcga_reduced)
npn_cm = confusionMatrix(npn_reduced)
un_cm = confusionMatrix(un_reduced)

saveRDS(log_cm, file=paste0(output_dir, "meta_log_cm.RDS"))

all_accs = data.frame(rbind(tdm_cm$overall, qn_cm$overall, log_cm$overall, npn_cm$overall, un_cm$overall, tcga_cm$overall))
all_accs = cbind(all_accs, method=c("TDM", "QN", "LOG", "NPN", "UNTR", "MA"))

nir = all_accs[1,]$AccuracyNull

all_accs$method = factor(all_accs$method, levels=all_accs$method)

# Plot accuracies
ci = aes(ymin=AccuracyLower, ymax=AccuracyUpper)
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(all_accs, aes(x=factor(method), y=Accuracy, color=method)) +
geom_point(size=3) +
geom_errorbar(ci, width=.3) +
ylab("Total Accuracy") +
xlab("Normalization") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
geom_hline(aes(yintercept = nir), linetype="longdash") +
scale_color_manual(values=cbPalette) 

ggsave(paste0(output_dir, ACCURACY_PLOT_FILENAME), plot=last_plot(), width=3, height=2.5)

cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = all_accs) +
ylab("Kappa") +
xlab("Normalization") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
geom_point(aes(x=factor(method), y=Kappa, color=method), shape=18, size=4) +
scale_color_manual(values=cbPalette) 

ggsave(paste0(output_dir, KAPPA_PLOT_FILENAME), plot=last_plot(), width=3, height=2.5)

meta_all = rbind(tdm_cm$byClass[,8],
  qn_cm$byClass[,8],
  log_cm$byClass[,8],
  npn_cm$byClass[,8],
  un_cm$byClass[,8],
  tcga_cm$byClass[,8])

meta_melted = melt(meta_all)
colnames(meta_melted) = c("Dataset", "Subtype", "Accuracy")
meta_melted$Dataset = c("TDM", "QN", "LOG", "NPN", "UNTR", "MA")
meta_melted$Dataset = factor(meta_melted$Dataset, levels=c("TDM", "QN", "LOG", "NPN", "UNTR", "MA"))

cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data=meta_melted, aes(x=Dataset, y=Accuracy, ymax=1)) +
  ylim(0.45,1) +
  coord_flip() + 
  geom_point(size=3, aes(colour=Dataset)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title=element_text(hjust=0)) +
  ylab("Balanced Accuracy") + theme_bw() + facet_wrap(~Subtype, ncol=5) +
  theme(legend.position="none") +
  scale_color_manual(values=cbPalette)

ggsave(paste0(output_dir, BALANCED_PLOT_FILENAME), plot=last_plot(), height=2.5, width=6)

write.table(all_accs, file=paste0(output_dir, ACCURACY_TABLE_FILENAME), sep='\t', row.names=T, col.names=T)
