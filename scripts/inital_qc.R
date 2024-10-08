##Initial QC
library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(DropletUtils)
library(ggpubr)
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(scran) 
library(scuttle)
library(gridExtra)
library(ggridges)

###Perform Initial QC on each sample after CogentAP has been performed
initial_qc <- function(stats, sample_name){
  stats2 <- subset(stats, Sample != 'Non_sample')
  stats2$mito_percent <- (stats2$Mitochondrial_Reads / stats2$Gene_Reads) *100
  G1 <- ggplot(stats2, aes(y=Barcoded_Reads, x = Sample, col = Sample)) +  geom_boxplot()  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name) 
  G2 <- ggplot(stats2, aes(y=Mapped_Reads, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G3 <- ggplot(stats2, aes(y=Exon_Reads, x = Sample, col = Sample)) +  geom_boxplot()  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G4 <- ggplot(stats2, aes(y=Intron_Reads, x = Sample, col = Sample)) +  geom_boxplot()  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G5 <- ggplot(stats2, aes(y=Intergenic_Reads, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G6 <- ggplot(stats2, aes(y=Gene_Reads, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G7 <- ggplot(stats2, aes(y=Mitochondrial_Reads, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G8 <- ggplot(stats2, aes(y=Ribosomal_Reads, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G9 <- ggplot(stats2, aes(y=No_of_Genes, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  G10 <- ggplot(stats2, aes(y=mito_percent, x = Sample, col = Sample)) +  geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) + NoLegend() + ggtitle(sample_name)
  g <- arrangeGrob(G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, nrow = 2) 
  ggsave(g, filename = paste0('inital_QC_', sample_name, '.png'))
}

preprocess <- function(raw_counts, stats, well, gene_info, samples, sample_name){
  print(dim(raw_counts)) 
  raw_counts$GeneID <- gsub("_.*","",raw_counts$GeneID) ###Remove gene symbol to just show ENSG 
  sample <- subset(stats, Sample %in% c(samples)) 
  sample_barcodes <- sample$Barcode  
  sample_barcodes <- c(sample_barcodes, "GeneID") ##include GeneID as a column header for filtering 
  counts <- raw_counts[(names(raw_counts) %in%  sample_barcodes)] ##Filter the gene df so we use only sample barcodes  
  print(dim(counts)) 
  well$Barcode <- gsub("\\+", "", well$Barcode) ##remove the + from the barcode in the wellist so it will match the genedf 
  subset <- gene_info %>% dplyr::select(2,3) #subset gene info to just ensemble_Id and gene_name 
  colnames(subset) <- c("GeneID", 'symbol') ##change column names 
  counts <- left_join(counts, subset, by= c('GeneID')) ##add gene information to counts 
  countsb <- remove_duplicate_genes(counts, column_of_symbol = 'symbol', method = 'mean') ##remove duplicate genes that have the same symbol based on which one has the greatest mean 
  print(summary(duplicated(rownames(countsb))))  
  counts <- countsb[,-1]  
  print(dim(counts)) 
  keep <- colnames(counts)   #update stats file so same barcodes match gene matrix file 
  stats <- filter(stats, Barcode %in% keep)   #update stats file so same barcodes match gene matrix file  
  row <- rownames(counts)   ##Order gene info to match same order as count matrix 
  gene_info <- gene_info[match(row, gene_info$Gene_Name), ]     ###Convert gene matrices to dgcMatrix format ready to make single cell experiment 
  mat <- as.matrix(counts) 
  mat <- as(mat, "dgCMatrix")## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
  colnames(mat) <- paste(sample_name, colnames(mat), sep="_")   ###Paste in the experiment name into the column names to make barcodes/cell IDs unique for stats file  
  stats2 <- stats[,-1] 
  rownames(stats2) <- stats[,1]  
  rownames(stats2) <- paste(sample_name, rownames(stats2), sep = "_")  
  preprocessed <- list(mat = mat, gene_info = gene_info, stats2 = stats2)  
  return(preprocessed)
}

SCE_QC <- function(data, sample_name){
  sce <- SingleCellExperiment(list(counts=data$mat), colData = data$stats2, rowData = data$gene_info) ##create single cell experiment
  sce <- addPerCellQCMetrics(sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(sce), value = FALSE))) ##Calculates QC metrics per cell included the % of mito reads
  colData(sce)$Complexity <- log10(sce@colData$detected) / log10(sce@colData$sum) ##Add complexity measurement
  df <- as.data.frame(colData(sce))
  G1 <- ggplot(df, aes(x=subsets_Mito_percent)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  G2 <- ggplot(df, aes(x=sum)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  G3 <- ggplot(df, aes(x=detected)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  G4 <- ggplot(df, aes(x=Complexity)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  g <- arrangeGrob(G1, G2, G3, G4, nrow = 2) #generates g
  ggsave(g, filename = paste0('QC_', sample_name, '.png'))
  return(sce)
}


raw_counts <- read.table("sample_name_analysis_incl_introns_genematrix.csv", header = T, sep = ',') ##CogentAP output, obtain from GSE235529
stats <- read.csv("sample_name_analysis_incl_introns_stats.csv") ##CogentAP output, can be found in files directory
well <- read.table("sample_name_WellList.TXT",  header = T,  sep="\t") ##Output from CellSelect, obtain from GSE235529
gene_info <- read.csv("gene_info_incl_introns.csv") ##output from CogentAP, is the same for each sample, can be found in files directory

sample_name <- "sample_name" # i.e. Participant_1

###Initial QC compares metrics from CogentAP output between single nuclei, postive control and negative control
initial_qc(stats = stats, 
           sample_name = sample_name)

###Organizes data by removing positive and negative control from the counts matrix, converts Ensemble ID to gene Symbol
###Removes rows with the same gene symbol based on lowest average expression
###Organizes and orders the stats file to match the cell/nuclei barcodes that are retained in the counts matrix after removal of postive and negative controls
data <- preprocess(raw_counts = raw_counts, 
                   stats = stats,
                   well = well, 
                   gene_info = gene_info,
                   samples = c('sample'),
                   sample_name = sample_name)

###Plots QC data
sce <- SCE_QC(data = data, sample_name = sample_name)

