# Author: Jeffrey Thompson
# Purpose: Performs multiclass lasso logistic regression to classify cancer subtypes.

# Get command line arguments.
#args <- commandArgs(trailingOnly = T)

source("package_loader.R")

NFOLDS = 100
NSEEDS = 10
INITIAL_SEED = 2
ACCURACY_PLOT_FILENAME = "meta_accuracies.pdf"
KAPPA_PLOT_FILENAME = "meta_kappas.pdf"
ACCURACY_TABLE_FILENAME = "meta_accuracies.tsv"

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
output_dir = args[12]

load_it(c("glmnet", "caret", "e1071", "stringr", "plyr"))

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
tdmnotshareddf = setup_acc_df(tdmnotshareddf)
qnnotshareddf = setup_acc_df(qnnotshareddf)
lognotshareddf = setup_acc_df(lognotshareddf)

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

# Filter any expression data with missing clinical annotation for tumor subtype.
trainT = trainT[!is.na(clin[,1][match(chartr('-', '.', rownames(trainT)), chartr('-', '.', rownames(clin)))]),]

# Pre-process the TCGA data.
tcga = preprocessTCGA(tcga_input, tcga_clin, FALSE)
log = preprocessTCGA(log_input, log_clin)
qn = preprocessTCGA(qn_input, qn_clin)
tdm = preprocessTCGA(tdm_input, tdm_clin)

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

# Pre-allocate vectors to hold results.
tcga_accs = vector(mode="list", length=length(seeds))
tdm_accs = vector(mode="list", length=length(seeds))
qn_accs = vector(mode="list", length=length(seeds))
log_accs = vector(mode="list", length=length(seeds))

# Perform a number of iterations equal to the number of random seeds.
# At each iteration, perform n-fold cross validation to build a model on the training data,
# then use that model to make predictions on the test data.
for(seednum in 1:length(seeds)) {
  set.seed(seeds[seednum])
  
  # Train a model for classifying cancer subtype on the METABRIC training data.
  set.seed(seeds[seednum])
  lasso.model=cv.glmnet(trainT, trainclin.subtypes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
  plot(lasso.model)

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

}

# Aggregate accuracies:
tdm_tables = lapply(tdm_accs, function(x) x$table)
qn_tables = lapply(qn_accs, function(x) x$table)
log_tables = lapply(log_accs, function(x) x$table)
tcga_tables = lapply(tcga_accs, function(x) x$table)

tdm_reduced = Reduce("+", tdm_tables) / length(tdm_tables)
qn_reduced = Reduce("+", qn_tables) / length(qn_tables)
log_reduced = Reduce("+", log_tables) / length(log_tables)
tcga_reduced = Reduce("+", tcga_tables) / length(tcga_tables)

tdm_reduced[1:5,1:5] = t(apply(tdm_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
qn_reduced[1:5,1:5] = t(apply(qn_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
log_reduced[1:5,1:5] = t(apply(log_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
tcga_reduced[1:5,1:5] = t(apply(tcga_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]

tdm_cm = confusionMatrix(tdm_reduced)
qn_cm = confusionMatrix(qn_reduced)
log_cm = confusionMatrix(log_reduced)
tcga_cm = confusionMatrix(tcga_reduced)

all_accs = data.frame(rbind(log_cm$overall, qn_cm$overall, tdm_cm$overall, tcga_cm$overall))
all_accs = cbind(all_accs, method=c("LOG", "QN", "TDM", "MA"))

nir = all_accs[1,]$AccuracyNull

all_accs$method = factor(all_accs$method, levels=all_accs$method)

# Plot accuracies
ci = aes(ymin=AccuracyLower, ymax=AccuracyUpper)
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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
