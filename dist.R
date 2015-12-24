source("directories.R")
source("package_loader.R")

load_it(c("reshape2", "ggplot2", "plyr"))

coad_cm = readRDS(paste0(output, "coad_log_cm.RDS"))
brca_cm = readRDS(paste0(output, "brca_log_cm.RDS"))
meta_cm = readRDS(paste0(output, "meta_log_cm.RDS"))

coad_sample_counts = colSums(coad_cm$table)
coad_total = sum(coad_sample_counts)
brca_sample_counts = colSums(brca_cm$table)
brca_total = sum(brca_sample_counts)
meta_sample_counts = colSums(meta_cm$table)
meta_total = sum(meta_sample_counts)

coad_sample_counts = data.frame(Proportion=coad_sample_counts /  coad_total)
brca_sample_counts = data.frame(Proportion=brca_sample_counts / brca_total)
meta_sample_counts = data.frame(Proportion=meta_sample_counts / meta_total)

coad_sample_counts = cbind(coad_sample_counts, Class=rownames(coad_sample_counts), dataset="RNA-seq")
brca_sample_counts = cbind(brca_sample_counts, Subtype=rownames(brca_sample_counts), dataset="RNA-seq")
meta_sample_counts = cbind(meta_sample_counts, Subtype=rownames(meta_sample_counts), dataset="RNA-seq")

brca_ref_classes = read.table(paste0(output, "brca_ref_classes.txt"))
coad_ref_classes = read.table(paste0(output, "coad_ref_classes.txt"))
meta_ref_classes = read.table(paste0(output, "meta_ref_classes.txt"))

brca_ref_counts = data.frame(plyr::count(brca_ref_classes))
brca_ref_total = sum(brca_ref_counts$freq)
coad_ref_counts = data.frame(plyr::count(coad_ref_classes))
coad_ref_total = sum(coad_ref_counts$freq)
meta_ref_counts = data.frame(plyr::count(meta_ref_classes))
meta_ref_total = sum(meta_ref_counts$freq)

colnames(brca_ref_counts) = c("Subtype", "Proportion")
brca_ref_counts = brca_ref_counts[,c("Proportion", "Subtype")]
brca_ref_counts$Proportion = as.numeric(as.character(brca_ref_counts$Proportion)) / brca_ref_total
colnames(coad_ref_counts) = c("Class", "Proportion")
coad_ref_counts = coad_ref_counts[,c("Proportion", "Class")]
coad_ref_counts$Proportion = as.numeric(as.character(coad_ref_counts$Proportion)) / coad_ref_total
colnames(meta_ref_counts) = c("Subtype", "Proportion")
meta_ref_counts = meta_ref_counts[,c("Proportion", "Subtype")]
meta_ref_counts$Proportion = as.numeric(as.character(meta_ref_counts$Proportion)) / meta_ref_total

brca_ref_counts = cbind(brca_ref_counts, dataset="Microarray")
coad_ref_counts = cbind(coad_ref_counts, dataset="Microarray")
meta_ref_counts = cbind(meta_ref_counts, dataset="METABRIC")

brca = rbind(brca_ref_counts, brca_sample_counts)
coad = rbind(coad_ref_counts, coad_sample_counts)
meta = rbind(meta_ref_counts, meta_sample_counts)

cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7")
ggplot(brca) + theme_bw() + geom_bar(stat="identity", aes(x=Subtype, y=Proportion, fill=Subtype)) + facet_wrap(~dataset) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_fill_manual(values=cbPalette) 
ggsave(paste0(output, "brca_sample_distribution.pdf"), plot=last_plot(), height=2.5, width=6)

cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7")
ggplot(coad) + theme_bw() + geom_bar(stat="identity", aes(x=Class, y=Proportion, fill=Class)) + facet_wrap(~dataset) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_fill_manual(values=cbPalette)
ggsave(paste0(output, "coad_sample_distribution.pdf"), plot=last_plot(), height=2.5, width=6)

cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7")
ggplot(meta) + theme_bw() + geom_bar(stat="identity", aes(x=Subtype, y=Proportion, fill=Subtype)) + facet_wrap(~dataset) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_fill_manual(values=cbPalette) 
ggsave(paste0(output, "meta_sample_distribution.pdf"), plot=last_plot(), height=2.5, width=6)
