###Performing decontx on each sample after clustering
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(celda)
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(scran) 
library(scuttle)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

###Filter nuclei based on QC metrics, perform initial clustering and then decontx ambient RNA clean up 

remove_duplicate_genes<-function(eset, column_of_symbol, method = "mean"){
  eset<-as.data.frame(eset)
  rownames(eset)<-NULL
  dups <- dim(eset)[1] - length(unique(eset[,column_of_symbol]))
  if(dups==0){
    eset<-tibble:: column_to_rownames(eset,var = column_of_symbol)
    return(eset)
  }else{
    if(method=="mean"){
      order_index=apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) mean(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>%dplyr:: distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }else if(method == "sd"){
      order_index = apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) sd(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>% distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }
  }
}

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

filter_cells <- function(sce, mito_percent_no, sum_no, detected_low, detected_high, complexity_no, sample_name){
  sce <- sce[,!(sce@colData$subsets_Mito_percent == "NaN")] 
  sce <- sce[,!(sce@colData$subsets_Mito_percent > mito_percent_no)] ####Remove nanowells with a mitochondrial % above this number, 20 is typically default
  sce <- sce[,(sce@colData$sum >= sum_no)] ##Remove nanowells with number of reads below this amount, we have sometimes not used this filer, or used 20,000
  sce <- sce[,(sce@colData$detected >= detected_low & sce@colData$detected <= detected_high)]  ##Remove nanowells if they have a low number of a high number of genes or high number of genes
  sce <- sce[,(sce@colData$Complexity >= complexity_no)] ##Remove nanowells with a low complexity, this means they have a lot of reads but not that many genes detected, for full-length snRNA-seq 0.65 is typical
  qc.stats <- perCellQCFilters(sce@colData) 
  print(colSums(as.matrix(qc.stats)))
  sce<- quickPerCellQC(sce) #Further remove cells that are the log-total counts at 3 MADs  from the median
  ##Make QC plots
  df <- as.data.frame(colData(sce))
  G1 <- ggplot(df, aes(x=subsets_Mito_percent)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  G2 <- ggplot(df, aes(x=sum)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  G3 <- ggplot(df, aes(x=detected)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  G4 <- ggplot(df, aes(x=Complexity)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')
  g <- arrangeGrob(G1, G2, G3, G4, nrow = 2) #generates g
  ggsave(g, filename = paste0('AfterQC_', sample_name, '.png'))
  return(sce)
}

generate_seurat <- function(sce, counts_cut_off){
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 0.1  ##filter genes based on default >= 0.1. Log10(0.1) = -1
  sce1 <- sce
  sce1 <- sce1[keep,]
  seurat <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))  ####Convert to seurat object
  filtered <- assay(sce1) 
  filtered <- CreateAssayObject(counts = filtered) ##add unfiltered assay for data exploration, this is not used for clustering etc
  seurat[["filtered"]] <- filtered
  ##Add feature level meta data to both assay slots
  seurat[["filtered"]] <- AddMetaData(seurat[["filtered"]],
                                      as.data.frame(SingleCellExperiment::rowData(sce1))) 
  seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]],
                                 as.data.frame(SingleCellExperiment::rowData(sce)))
  
  ##Get  protein coding genes only for clustering
  genes <- seurat@assays$filtered@meta.features 
  PC <- genes %>% dplyr::filter(Gene_Biotype == 'protein_coding') ##Protein coding genes DF
  PCoding <- PC$Gene_Name #Protein coding genes list
  ##Make a new assay with protein coding genes only 
  filtered <- as.data.frame(seurat@assays$filtered@counts) 
  pcode <- filtered[rownames(filtered) %in% PCoding,]
  pcoding <- CreateAssayObject(counts = pcode)
  seurat[["Pcoding"]] <- pcoding
  print(Assays(seurat))
  return(seurat)
}

remove_genes <- function(seurat){
  seurat <- seurat[!grepl("^MT-", rownames(seurat)), ] #filter mitochondrial genes
  seurat <- seurat[!grepl("^HB[^(P)]", rownames(seurat)), ] 
  seurat <- seurat[!grepl("MALAT1", rownames(seurat)), ] 
  seurat <- seurat[!grepl("NEAT1", rownames(seurat)), ] 
  return(seurat)
}

Run_PCA <- function(seurat, assay, number_of_features, sample_name){
  DefaultAssay(object = seurat) <- assay
  seurat <- SCTransform(seurat, vst.flavor = "v2", assay = assay, verbose = FALSE, variable.features.n = number_of_features)
  seurat <- RunPCA(seurat, verbose = FALSE)
  G1 <- DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
  ggsave(G1, filename = paste0('PCA_Heatmap_', sample_name, '.png'))
  G2 <- ElbowPlot(seurat)
  ggsave(G2, filename = paste0('PCA_Elbow_', sample_name, '.png'))
  return(seurat)
}

cluster <- function(seurat, PCA_dim, res, sample_name){
  DefaultAssay(object = seurat) <- 'SCT'
  seurat <- RunUMAP(seurat, dims = PCA_dim)
  seurat <- RunTSNE(seurat, dims = PCA_dim)
  seurat <- FindNeighbors(seurat, dims = PCA_dim)
  seurat <- FindClusters(seurat, resolution = res)
  G1 <- DimPlot(seurat, reduction = "umap", label = TRUE)
  G2 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
  g <- arrangeGrob(G1, G2, ncol = 2) #generates g
  ggsave(g, filename = paste0('UMAP_TSNE_', sample_name, '.png'))
  return(seurat)
}

run_decontx <- function(seurat, data){
  sce <- SingleCellExperiment(list(counts = seurat@assays$RNA@counts), 
                              colData = seurat@meta.data, rowData = data$gene_info)
  sce <- decontX(sce, z = colData(sce)$CellType)
  print(colData(sce)$decontX_contamination)
  return(sce)
}

create_seurat_with_decont <- function(sce, counts_cut_off){
  seurat1 <- CreateSeuratObject(counts = decontXcounts(sce), meta.data = as.data.frame(colData(sce)))
  seurat1[["RNA"]] <- AddMetaData(seurat1[["RNA"]],
                                  as.data.frame(SingleCellExperiment::rowData(sce)))
  cols <- c('sum.1', 'detected.1', 'total.1', 'nCount_filtered', 'nFeature_filtered', 'nCount_Pcoding', 'nCount_SCT', 'nFeature_SCT', 
            'SCT_snn_res.0.2', 'SCT_snn_res.0.3', 'SCT_snn_res.0.4', 'SCT_snn_res.0.5', 'SCT_snn_res.0.6', 'SCT_snn_res.0.7', 
            'SCT_snn_res.0.8', 'SCT_snn_res.0.9', 'seurat_clusters', 'CellType', 'ident', 'sum.2', 'detected.2', 'subsets_Mito_sum.1', 
            'subsets_Mito_detected.1', 'subsets_Mito_percent.1', 'altexps_filtered_sum', 'altexps_filtered_detected', 'altexps_filtered_percent',
            'altexps_Pcoding_sum', 'altexps_Pcoding_detected','altexps_Pcoding_percent',  "altexps_SCT_sum", "altexps_SCT_detected", "altexps_SCT_percent",
            "total.2",  "decontX_clusters")
  seurat1@meta.data <- seurat1@meta.data[, ! colnames(seurat1@meta.data) %in% cols]
  ave.counts <- rowMeans(decontXcounts(sce))
  keep <- ave.counts >= counts_cut_off  ##filter genes based on default >= 0.1. Log10(0.1) = -1
  sce1 <- sce
  sce1 <- sce1[keep,]
  counts = decontXcounts(sce1)
  filtered <- CreateAssayObject(counts = counts) ##add unfiltered assay for data exploration, this is not used for clustering etc
  seurat1[["filtered"]] <- filtered
  seurat1@assays$filtered@counts[1:6, 1:6]
  ##Add feature level meta data to both assay slots
  seurat1[["filtered"]] <- AddMetaData(seurat1[["filtered"]],
                                       as.data.frame(SingleCellExperiment::rowData(sce1))) 
  seurat1@assays$filtered@counts[rownames(seurat1@assays$filtered@counts) %in% 'CD36', ][1:10] 
  ##Get  protein coding genes only for clustering
  genes <- seurat1@assays$filtered@meta.features 
  PC <- genes %>% dplyr::filter(Gene_Biotype == 'protein_coding') ##Protein coding genes DF
  PCoding <- PC$Gene_Name #Protein coding genes list
  ##Make a new assay with protein coding genes only 
  filtered <- as.data.frame(seurat1@assays$filtered@counts) 
  pcode <- filtered[rownames(filtered) %in% PCoding,]
  pcoding <- CreateAssayObject(counts = pcode)
  seurat1[["Pcoding"]] <- pcoding
  return(seurat1)
}
Add_Well_To_Meta <- function(well, seurat, sample_name){
  well$Barcode <- gsub("\\+", "", well$Barcode)
  well$Barcode  <- paste(sample_name, well$Barcode, sep = "_")
  well2 <- well %>% dplyr::select(Barcode, Sample, Size1, Circularity1, Integ.Signal1, SampleWell, Dispense.tip)
  barcodes <- rownames(seurat@meta.data)
  data <- filter(well2, Barcode %in% barcodes) 
  rownames(data) <- data$Barcode
  seurat <- AddMetaData(seurat, data)
  return(seurat)
}


##Load filtered data from cellranger for preliminary analysis
raw_counts <- read.table("sample_name_analysis_incl_introns_genematrix.csv", header = T, sep = ',') ##CogentAP output
stats <- read.csv("sample_name_analysis_incl_introns_stats.csv") ##CogentAP output
well <- read.table("sample_name_WellList.TXT",  header = T,  sep="\t") ##Output from CellSelect 
gene_info <- read.csv("sample_name_gene_info_incl_introns.csv")

sample_name <- "sample_name" # i.e. Participant_1

###Organizes data by removing positive and negative control from the counts matrix, converts Ensemble ID to gene Symbol
###Removes rows with the same gene symbol based on lowest average expression
###Organizes and orders the stats file to match the cell/nuclei barcodes that are retained in the counts matrix after removal of postive and negative controls
data <- preprocess(raw_counts = raw_counts, 
                    stats = stats,
                    well = well, 
                    gene_info = gene_info,
                    samples = c('sample'),
                    sample_name = sample_name)


###Create a single cell experiment and plots QC metrics
sce <- SCE_QC(data = data, sample_name = sample_name)

##Filter cells based on metric
sce <- filter_cells(sce = sce, 
                    mito_percent_no = 30, ####Remove nanowells with a mitochondrial % above this number
                    sum_no = 10000, ##Remove nanowells with number of reads below this amount
                    detected_low = 500, ##remove nanowells with that detect no of genes lower than this number, may indicate poor quality nuclei,
                    detected_high = 20000, ##remove nanowells with that detect no of genes higher than this number, could be an indicator of doublet or a lot of ambient RNA,
                    complexity_no = 0.65, ##remove nanowells with low complexity, i.e. a lot of reads but not many genes detected, for full-length snRNA-seq 0.65 is typical
                    sample_name = sample_name
)



###Plot how frequently the genes are expressed
ave.counts <- rowMeans(counts(sce))
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))


##Filter low expressed genes
seurat <- generate_seurat(sce = sce,
                          counts_cut_off = 0.1  ##filter genes based on default >= 0.1. Log10(0.1) = -1
)

##Run and plot PCA
seurat <- Run_PCA(seurat = seurat,
                  assay = "Pcoding", ##usually use used Protein coding assay
                  number_of_features = 5000, ##number of features used for clustering
                  sample_name = sample_name)

##Cluster seurat object 
seurat <- cluster(seurat = seurat,
                  PCA_dim = 1:8, ##number of PCA dimensions used
                  res = c(0.4, 0.6, 0.8, 1, 1.2, 1.4), ##Different resolutions to determine the different clusters
                  sample_name = sample_name
)

##Generic adipose tissue markers, list is not exhaustive 
AT <- c('PDGFRA',
        'PDGFRB',
        "DCN",
        'ZNF423',
        "CD38",
        'CTNNA2',
        "CEBPA",
        "STAT5A",
        'LIPE',
        'LEP',
        'ADIPOQ',
        "DGAT2",
        'PLIN1',
        'LPL',
        "PPARG",
        "TBX15",
        "IRS1",
        "HK2",
        'VWF',
        'CDH5',
        'PECAM1',
        'PTPRC',
        'HLA-DRB1',
        'HLA-DPA1',
        'ITGAX',
        'TREM2',
        'KIT',
        'TPSAB1'
)

##Normalize RNA assay slot
DefaultAssay(object = seurat) <- "RNA"
seurat <- NormalizeData(seurat)
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

###decide which resolution to use based on expression of makers
P1 <- DotPlot(seurat, features = AT, 
              scale.min = 0,
              group.by = "SCT_snn_res.0.8"
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(seurat, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.8" )

P1 + P2

###Preliminary classification of cell types. Please note this will change for sample to sample this is just an example.
seurat@meta.data <- seurat@meta.data %>% mutate(CellType = case_when(
  SCT_snn_res.0.8 == "0"  ~ "Pre_Ad", 
  SCT_snn_res.0.8 == "1"  ~ "EC",
  SCT_snn_res.0.8 == "2"  ~ "Adip_1",  
  SCT_snn_res.0.8 == "3"  ~ "Stem",
  SCT_snn_res.0.8 == "4"  ~ "Adip_2",
  SCT_snn_res.0.8 == "5"  ~ "Immune"
)) 

###Add nanowell information to meta data including nuclei size, well of origin etc
seurat <- Add_Well_To_Meta(well = well, seurat = seurat, sample_name = sample_name)

##Save RDS
saveRDS(seurat, file = paste0('seurat_', sample_name, '.RDS'))

# Create a SingleCellExperiment object and run decontX
sce <- run_decontx(seurat, data)

plotDecontXContamination(sce)

ave.counts <- rowMeans(decontXcounts(sce))
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))



##Create a seurat object with the decontx counts
seurat1 <- create_seurat_with_decont(sce, counts_cut_off = 0.1)

##Remove contaminant genes i.e. mitochondrial
seurat1 <- remove_genes(seurat1)


##Run and plot PCA
seurat1 <- Run_PCA(seurat = seurat1,
                   assay = "Pcoding", ##usually use used Protein coding assay 
                   number_of_features = 5000, ##number of features used for clustering
                   sample_name = sample_name)

##Cluster seurat object
seurat1 <- cluster(seurat = seurat1,
                   PCA_dim = 1:15, ##number of PCA dimensions used
                   res = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), ##Different resolutions to determine the different clusters
                   sample_name = sample_name
)


##Normalize RNA assay slot
DefaultAssay(object = seurat1) <- "RNA"
seurat1 <- NormalizeData(seurat1)
all.genes <- rownames(seurat1)
seurat1 <- ScaleData(seurat1, features = all.genes)

##Check data to see if decontx has improved clustering

P3 <- DimPlot(seurat1, reduction = "umap", group.by = 'SCT_snn_res.0.3')
AT <- c('PDGFRA',
        'PDGFRB',
        "DCN",
        'ATXN1',
        'ZNF423',
        "CD38",
        'CTNNA2',
        "CEBPA",
        "STAT5A",
        'LIPE',
        'LEP',
        'ADIPOQ',
        "DGAT2",
        'PLIN1',
        'LPL',
        'CD36',
        'FASN',
        "PPARG",
        'FABP4',
        'SAA1',
        "TBX15",
        "IRS1",
        "HK2",
        'VWF',
        'CDH5',
        'PECAM1',
        'PTPRC',
        'HLA-DRB1',
        'HLA-DPA1',
        'ITGAX',
        'TREM2',
        'KIT' #mast 
)

P4 <- DotPlot(seurat1, features = AT, scale.min = 0, group.by = 'SCT_snn_res.0.3')  +  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))
P3 + P4

###Add nanowell information to meta data including nuclei size, well of origin etc
seurat1 <- Add_Well_To_Meta(well = well, seurat = seurat1, sample_name = sample_name)
seurat1@meta.data

##Remove meta.data columns no longer needed
cols <- c('SCT_snn_res.0.2', "SCT_snn_res.0.4", 'SCT_snn_res.0.5', 'SCT_snn_res.0.6', 'SCT_snn_res.0.7', 'SCT_snn_res.0.8',
          'SCT_snn_res.0.9', 'seurat_clusters', 'CellType')
seurat1@meta.data <- seurat1@meta.data[, ! colnames(seurat1@meta.data) %in% cols]

##check for batch effects of sample well and dispense tip
DimPlot(seurat1, reduction = "umap", group.by = "SampleWell")
DimPlot(seurat1, reduction = "umap", group.by = "Dispense.tip")

##Save Seurat RDS file
saveRDS(seurat1, file = paste0('seurat_', sample_name, 'decont', '.RDS'))



