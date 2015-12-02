# Author: Jeffrey Thompson
# Purpose: Performs multiclass lasso logistic regression to classify cancer subtypes.

# Get command line arguments.

#args = commandArgs(trailingOnly = T)

source("package_loader.R")

NFOLDS = 100
NSEEDS = 10
INITIAL_SEED = 2
ACCURACY_PLOT_FILENAME = "coadread_accuracies.pdf"
KAPPA_PLOT_FILENAME = "coadread_kappas.pdf"
ACCURACY_TABLE_FILENAME = "coadread_accuracies.tsv"

args = c('/home/jeff/Repos/tdm2015test2/tdm_auto/normalized_data/', 
		'COADREAD_ZEROONE.pcl', 
		'COADREADClin.tsv', 
		'COADREAD_TDM_ZEROONE.pcl', 
		'COADREADRNASeqClin.tsv', 
		'COADREAD_QN_ZEROONE.pcl',
		'COADREADRNASeqClin.tsv', 
		'COADREAD_LOG_ZEROONE.pcl', 
		'COADREADRNASeqClin.tsv', 
    'COADREAD_NPN_ZEROONE.pcl',
    'COADREADRNASeqClin.tsv',
		'/home/jeff/Repos/tdm2015test2/tdm_auto/output/')

input_dir = args[1]
ref_input = args[2]
ref_clin = args[3]
tdm_input = args[4]
tdm_clin = args[5]
qn_input = args[6]
qn_clin = args[7]
log_input = args[8]
log_clin = args[9]
npn_input = args[10]
npn_clin = args[11]
output_dir = args[12]

load_it(c("glmnet", "caret", "e1071", "stringr", "plyr"))

# This function creates dataframes to hold the accuracies for each class across runs.
setup_acc_df = function(df) {
  df = data.frame(matrix(nrow=4, ncol=0))
  rownames(df) = c("CIMP", "CIMPL", "NCIMP", "Normal")
  return (df)
}

# Define data.frames to hold classification results.
refdf = setup_acc_df(refdf)
tdmdf = setup_acc_df(tdmdf)
qndf = setup_acc_df(qndf)
logdf = setup_acc_df(logdf)
npndf = setup_acc_df(npndf)
tdmnotshareddf = setup_acc_df(tdmnotshareddf)
qnnotshareddf = setup_acc_df(qnnotshareddf)
lognotshareddf = setup_acc_df(lognotshareddf)
npnnotshareddf = setup_acc_df(npnnotshareddf)

# Define random seeds.
set.seed(INITIAL_SEED)
seeds = sample(1:10000, NSEEDS)

message(paste("Random seeds:", paste(seeds, collapse=", ")), appendLF=TRUE)

# Pre-processes test data to ready it for the LASSO.
# Returns a list with first entry the pre-processed dataset and second entry
# the classes for those data.
preprocessTCGA = function(dataFile, clinFile) {
  data = read.table(paste0(input_dir, dataFile), 
    sep='\t', 
    row.names=1, 
    header=T)

  dataT = t(data)
  
  # Remove duplicate samples from the data.
  dataT = dataT[str_detect(substr(rownames(dataT),1,15), "(\\.01)|(\\.11)"),]  

  # Read clinical data.
  clin = read.table(paste0(input_dir, clinFile), sep='\t', row.names=1, header=T)

  # Make normal samples into their own class.
  clin[,2] = as.character(clin[,2])
  clin[clin[,1]=="Solid Tissue Normal",2] = "NORMAL"
  clin[,2] = as.factor(clin[,2])

  # Remove samples without cancer labels.
  clin = subset(clin, clin[,1] != "")
  clin[,1] = droplevels(clin[,1])

  # Remove samples without class labels.
  clin = subset(clin, clin[,2] != "")
  clin[,2] = droplevels(clin[,2])

  # Remove samples with unused class.
  clin = subset(clin, clin[,2] != "Low purity c2")
  clin[,2] = droplevels(clin[,2])

  # Map long class names to simple class names.
  old_names = c("COADREAD non-CIMP c11", "COADREAD CIMPL c10", "COADREAD CIMP c12", "NORMAL")
  new_names = c("NCIMP", "CIMPL", "CIMP", "NORMAL")
  map = setNames(new_names, old_names)
  clin[,2] = sapply(clin[,2], function(x) map[as.character(x)])

  # Filter clinical data to include only those samples with expression data.
  clin.response = clin[,2][match(substr(chartr('-', '.', rownames(dataT)), 1, 15), chartr('-', '.', rownames(clin)))]
  clin.response = na.omit(clin.response)

  # Filter expression data to include only those samples with ids also in clinical data.
  dataT = dataT[substr(chartr('-', '.', rownames(dataT)), 1, 15) %in% chartr('-', '.', rownames(clin)),]
  
  # Filter any data with missing clinical annotation for tumor class.
  dataT = dataT[!is.na(clin[,2][match(substr(chartr('-', '.', rownames(dataT)),1,15), 
    substr(chartr('-', '.', rownames(clin)),1,15))]),]
  
  return (list(data=dataT, response=clin.response))
} # end function preprocessTCGA

# Predict subtypes on TCGA data using previously trained model.
predictTCGA = function(title, data, model) { 
  # Predict classes based on the trained model.
  lasso.predict = predict(model, data$data, s = "lambda.1se", type = "class")
  
  # Make sure all factors are included.
  lasso.predict = factor(lasso.predict, c("CIMP", "CIMPL", "NCIMP", "NORMAL"))
  
  # Build a contingency table of the results.
  con_table = table(pred = lasso.predict, true = data$response, exclude="none")
  
  # Generate statistics based on the contigency table.
  cm = confusionMatrix(con_table)

  return(cm)
}

# Load the training data.
train = read.table(paste0(input_dir, ref_input), sep='\t', row.names=1, header=T)

# Transpose rows and columns.
trainT = t(train)

# Remove duplicate samples from the data.
trainT = trainT[str_detect(substr(rownames(trainT),1,15), "(\\.01)|(\\.11)"),]  

# Load the clinical features.
clin = read.table(paste0(input_dir, ref_clin), sep='\t', row.names=1, header=T)

# Make normal samples into their own class.
clin[,2] = as.character(clin[,2])
clin[clin[,1]=="Solid Tissue Normal",2] = "NORMAL"
clin[,2] = as.factor(clin[,2])

# Remove samples without cancer labels.
clin = subset(clin, clin[,1] != "")
clin[,1] = droplevels(clin[,1])

# Remove samples without class labels.
clin = subset(clin, clin[,2] != "")
clin[,2] = droplevels(clin[,2])

# Remove samples with unused class.
clin = subset(clin, clin[,2] != "Low purity c2")
clin[,2] = droplevels(clin[,2])

# Map long class names to simple class names.
old_names = c("COADREAD non-CIMP c11", "COADREAD CIMPL c10", "COADREAD CIMP c12", "NORMAL")
new_names = c("NCIMP", "CIMPL", "CIMP", "NORMAL")
map = setNames(new_names, old_names)
clin[,2] = sapply(clin[,2], function(x) map[as.character(x)])

# Filter the samples to remove those without classification.
trainT = trainT[substr(rownames(trainT),1,15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]

# Filter clinical data to include only those entries that also have expression data.
trainclin.classes = clin[,2][match(substr(chartr('-', '.', rownames(trainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]
trainclin.classes = na.omit(trainclin.classes)

# Filter any expression data with missing clinical annotation for tumor class.
trainT = trainT[!is.na(clin[,2][match(substr(chartr('-', '.', rownames(trainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]),]


# Pre-process the TCGA data.
tdm = preprocessTCGA(tdm_input, tdm_clin)
qn = preprocessTCGA(qn_input, qn_clin)
log = preprocessTCGA(log_input, log_clin)
npn = preprocessTCGA(npn_input, npn_clin)

# Determine which samples are shared between training and test data.
tdm_shared = substr(chartr('-', '.', rownames(tdm$data)),1,15) %in% substr(chartr('-', '.', rownames(trainT)),1,15)  
qn_shared = substr(chartr('-', '.', rownames(qn$data)),1,15) %in% substr(chartr('-', '.', rownames(trainT)),1,15)
log_shared = substr(chartr('-', '.', rownames(log$data)),1,15) %in% substr(chartr('-', '.', rownames(trainT)),1,15)
npn_shared = substr(chartr('-', '.', rownames(npn$data)),1,15) %in% substr(chartr('-', '.', rownames(trainT)),1,15)


# Filter data to include only those samples not shared between training and test.
tdm$data = tdm$data[!tdm_shared,]
tdm$response = tdm$response[!tdm_shared]
qn$data = qn$data[!qn_shared,]
qn$response = qn$response[!qn_shared]
log$data = log$data[!log_shared,]
log$response = log$response[!log_shared]
npn$data = npn$data[!npn_shared,]
npn$response = npn$response[!npn_shared]

# Pre-allocate vectors for classification results
tdm_accs = vector(mode="list", length=length(seeds))
qn_accs = vector(mode="list", length=length(seeds))
log_accs = vector(mode="list", length=length(seeds))
npn_accs = vector(mode="list", length=length(seeds))

# Perform a number of iterations equal to the number of random seeds.
# At each iteration, perform n-fold cross validation to build a model on the training data,
# then use that model to make predictions on the test data.
for(seednum in 1:length(seeds)) {
  # Train a model for classifying cancer subtype.
  set.seed(seeds[seednum])
  lasso.model=cv.glmnet(trainT, trainclin.classes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
  plot(lasso.model)
  
  acc = predictTCGA("TDM Results", tdm, lasso.model)
  tdm_accs[[seednum]] = acc
  tdmdf = cbind(tdmdf,as.vector(acc$byClass[,8]))
  
  acc = predictTCGA("QN Results", qn, lasso.model)
  qn_accs[[seednum]] = acc
  qndf = cbind(qndf,as.vector(acc$byClass[,8]))

  acc = predictTCGA("LOG Results", log, lasso.model)
  log_accs[[seednum]] = acc
  logdf = cbind(logdf,as.vector(acc$byClass[,8]))

  acc = predictTCGA("NPN Results", npn, lasso.model)
  npn_accs[[seednum]] = acc
  npndf = cbind(npndf, as.vector(acc$byClass[,8]))
}

# Build a table of accuracies across all datasets.
accuracies = data.frame(NCIMP=numeric(0), CIMPL=numeric(0), CIMP=numeric(0), NORMAL=numeric(0))
accuracies = rbind(accuracies, apply(tdmdf,1,function(x) mean(x)))
accuracies = rbind(accuracies, apply(qndf,1,function(x) mean(x)))
accuracies = rbind(accuracies, apply(logdf,1,function(x) mean(x)))
accuracies = rbind(accuracies, apply(npndf,1,function(x) mean(x)))
rownames(accuracies) = c("TDM", "QN", "LOG", "NPN")
colnames(accuracies) = c("CIMP", "CIMPL", "NCIMP", "NORMAL")

# write.table(accuracies, file=paste0(output_dir,"/COADlassoSubtypeAccuracies.txt"), sep='\t', row.names=T, col.names=T)

# Aggregate accuracies:
tdm_tables = lapply(tdm_accs, function(x) x$table)
qn_tables = lapply(qn_accs, function(x) x$table)
log_tables = lapply(log_accs, function(x) x$table)
npn_tables = lapply(npn_accs, function(x) x$table)

tdm_reduced = Reduce("+", tdm_tables) / length(tdm_tables)
qn_reduced = Reduce("+", qn_tables) / length(qn_tables)
log_reduced = Reduce("+", log_tables) / length(log_tables)
npn_reduced = Reduce("+", npn_tables) / length(npn_tables)

tdm_reduced[1:4,1:4] = t(apply(tdm_reduced, 1, function(x) as.integer(round(x))))[1:4,1:4]
qn_reduced[1:4,1:4] = t(apply(qn_reduced, 1, function(x) as.integer(round(x))))[1:4,1:4]
log_reduced[1:4,1:4] = t(apply(log_reduced, 1, function(x) as.integer(round(x))))[1:4,1:4]
npn_reduced[1:4,1:4] = t(apply(npn_reduced, 1, function(x) as.integer(round(x))))[1:4,1:4]

tdm_cm = confusionMatrix(tdm_reduced)
qn_cm = confusionMatrix(qn_reduced)
log_cm = confusionMatrix(log_reduced)
npn_cm = confusionMatrix(npn_reduced)

all_accs = data.frame(rbind(tdm_cm$overall, qn_cm$overall, log_cm$overall, npn_cm$overall))

all_accs = cbind(all_accs, method=c("TDM", "QN", "LOG", "NPN"))

nir = all_accs[1,]$AccuracyNull

# Plot accuracies
ci = aes(ymin=AccuracyLower, ymax=AccuracyUpper)
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cbPalette = c("#56B4E9", "#E69F00", "#000000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(all_accs, aes(x=factor(method), y=Accuracy, color=method)) +
geom_point(size=3) +
geom_errorbar(ci, width=.2) +
ylab("Total Accuracy") +
xlab("Normalization") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
geom_hline(aes(yintercept = nir), linetype="longdash") +
scale_color_manual(values=cbPalette)  

ggsave(paste0(output_dir, ACCURACY_PLOT_FILENAME), plot=last_plot(), width=3, height=3)

cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = all_accs) +
ylab("Kappa") +
xlab("Normalization") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
geom_point(aes(x=factor(method), y=Kappa, color=method), shape=18, size=4) +
scale_color_manual(values=cbPalette) 

ggsave(paste0(output_dir, KAPPA_PLOT_FILENAME), plot=last_plot(), width=3, height=3)

write.table(all_accs, file=paste0(output_dir, ACCURACY_TABLE_FILENAME), sep='\t', row.names=T, col.names=T)
