####Integtation comparisons
## Load libraries
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(Seurat) 
library(scran) 
library(scuttle)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library("patchwork")
library(STACAS)
library(SeuratWrappers)
library(harmony)
library(pdfCluster)
library(lisi)
library(kBET)


###Batch correction comparison
Participant <- "orig.ident"
AgeGroup <- 'AgeGroup'


###Load individual RDS seurat objects 
##Read RDS
B1 <- readRDS('B1_decont.RDS')
B6 <- readRDS('B6_decont.RDS')
B7 <- readRDS('B7_decont.RDS')
B8 <- readRDS('B8_decont.RDS')
B9 <- readRDS('B9_decont.RDS')
B10 <- readRDS('B10_decont.RDS')
B14 <- readRDS('B14_decont.RDS')
B15 <- readRDS('B15_decont.RDS')
B16 <- readRDS('B16_decont.RDS')
B17 <- readRDS('B17_decont.RDS')
B18 <- readRDS('B18_decont.RDS')
B19 <- readRDS('B19_decont.RDS')
B20 <- readRDS('B20_decont.RDS')
B22 <- readRDS('B22_decont.RDS')
B23 <- readRDS('B23_decont.RDS')
B24 <- readRDS('B24_decont.RDS')
B26 <- readRDS('B26_decont.RDS')
B27 <- readRDS('B27_decont.RDS')
B30 <- readRDS('B30_decont.RDS')
B31 <- readRDS('B31_decont.RDS')

##Change to protein coding assay
DefaultAssay(object = B1) <- "Pcoding"
DefaultAssay(object = B6) <- "Pcoding"
DefaultAssay(object = B7) <- "Pcoding"
DefaultAssay(object = B8) <- "Pcoding"
DefaultAssay(object = B9) <- "Pcoding"
DefaultAssay(object = B10) <- "Pcoding"
DefaultAssay(object = B14) <- "Pcoding"
DefaultAssay(object = B15) <- "Pcoding"
DefaultAssay(object = B16) <- "Pcoding"
DefaultAssay(object = B17) <- "Pcoding"
DefaultAssay(object = B18) <- "Pcoding"
DefaultAssay(object = B19) <- "Pcoding"
DefaultAssay(object = B20) <- "Pcoding"
DefaultAssay(object = B22) <- "Pcoding"
DefaultAssay(object = B23) <- "Pcoding"
DefaultAssay(object = B24) <- "Pcoding"
DefaultAssay(object = B26) <- "Pcoding"
DefaultAssay(object = B27) <- "Pcoding"
DefaultAssay(object = B30) <- "Pcoding"
DefaultAssay(object = B31) <- "Pcoding"


#Merge data sets
Merge <- merge(B1, y = c(B6, B7, B8, B9, B10, B14, B15, B16, B17, B18, B19, B20, B22, B23, B24, B26, B27, B30, B31), project = "BioAge") ##Variable number 8 pick project name
DefaultAssay(object = Merge) <- "Pcoding"

###Add Age group to meta data
Merge@meta.data <- Merge@meta.data %>% mutate(AgeGroup = case_when(
  orig.ident == "B1" ~ "Older",
  orig.ident == "B6" ~ "Younger", 
  orig.ident == "B7" ~ "Older",
  orig.ident == "B8" ~ "Younger",
  orig.ident == "B9" ~ "Older",
  orig.ident == "B10" ~ "Younger",
  orig.ident == "B14" ~ "Older",
  orig.ident == "B15" ~ "Younger",
  orig.ident == "B16" ~ "Younger",
  orig.ident == "B17" ~ "Older",
  orig.ident == "B18" ~ "Younger",
  orig.ident == "B19" ~ "Younger",
  orig.ident == "B20" ~ "Younger",
  orig.ident == "B22" ~ "Older",
  orig.ident == "B23" ~ "Younger",
  orig.ident == "B24" ~ "Older",
  orig.ident == "B26" ~ "Older",
  orig.ident == "B27" ~ "Older",
  orig.ident == "B30" ~ "Younger",
  orig.ident == "B31" ~ "Older",
))


unique(Merge@meta.data$AgeGroup)
unique(Merge@meta.data$orig.ident)

###save Merged seurat. 
saveRDS(Merge, 'merged.RDS')

rm(list = ls())

##Load merged seurat
Merge <- readRDS("merged.RDS")
Participant <- "orig.ident"
AgeGroup <- 'AgeGroup'

#######Looked at merged data without integration###########
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

Merge <- Run_PCA(Merge, assay = 'Pcoding', number_of_features = 5000, sample_name = sample_name)

##Cluster 
cluster <- function(seurat, PCA_dim, res, sample_name){
  DefaultAssay(object = seurat) <- 'SCT'
  seurat <- RunUMAP(seurat, dims = PCA_dim)
  seurat <- RunTSNE(seurat, dims = PCA_dim)
  seurat <- FindNeighbors(seurat, dims = PCA_dim)
  seurat <- FindClusters(seurat, resolution = res)
  return(seurat)
}

seurat <- cluster(seurat=Merge, PCA_dim = 1:10, res = c(0.2, 0.4, 0.6, 0.8, 1, 1.2), sample_name = sample_name)


DimPlot(seurat, split.by = 'orig.ident', group.by = 'orig.ident')

# plot merged data with no integration
P1 <- DimPlot(seurat, group.by = c(Participant)) 
P2 <- DimPlot(seurat, group.by = c(AgeGroup))
P3 <- DimPlot(seurat, group.by = c("SCT_snn_res.0.4")) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "NoBatchCorrection.SCT") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("uncorrected.SCT.png")

# Add umap coordinates, extract meta and save for later
umap <- data.frame(seurat@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("Uncorrected_SCT",  colnames(umap))
seurat <- AddMetaData(seurat, metadata = umap)
meta <- seurat@meta.data
write.csv(meta, 'meta.data.uncorrected.SCT.csv')

rm(list = ls())

################Run RPCA ######################
rm(list = ls())
Merge <- readRDS("merged.RDS")

start <- Sys.time()
###Split object by participant 
Merge.list <- SplitObject(Merge, split.by = "orig.ident") 


#normalize and identify variable features for each dataset independently
for (i in 1:length(Merge.list)) {
  suppressWarnings(Merge.list[[i]] <- SCTransform(Merge.list[[i]], vst.flavor = "v2", assay = 'Pcoding', method = "glmGamPoi",
                                                  variable.features.n = 5000, ##Can change this number if needed, default to 5000
                                                  verbose = FALSE, return.only.var.genes=FALSE))
}

features <- SelectIntegrationFeatures(object.list = Merge.list, nfeatures = 5000)
Merge.list <- PrepSCTIntegration(object.list = Merge.list, anchor.features = features)
Merge.list <- lapply(X = Merge.list, 
                     FUN = RunPCA,
                     features = features)

##find integration anchors
anchors <- FindIntegrationAnchors(object.list = Merge.list, normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:20, 
                                  reduction = "rpca", 
                                  k.anchor = 10, 
) 
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                              #k.weight = 100,
                              dims = 1:20)

DefaultAssay(combined.sct) <- "integrated"
combined.sct <- ScaleData(combined.sct, 
                          #assay = "SCT", 
                          verbose = FALSE)
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
ElbowPlot(combined.sct) 

####Keep the same numner of PCAs 1:10
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:10)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:10)
combined.sct <- FindClusters(combined.sct, resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2))

end <- Sys.time()
end - start # Time difference of 

Participant <- "orig.ident"
AgeGroup <- 'AgeGroup'


# plot data
P1 <- DimPlot(combined.sct, group.by = c(Participant)) 
P2 <- DimPlot(combined.sct, group.by = c(AgeGroup))
P3 <- DimPlot(combined.sct, group.by = c("integrated_snn_res.0.4")) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "RPCA_k.anchor.10") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("RPCA_k.anchor.10.png")

DimPlot(combined.sct, split.by = c(Participant))
DimPlot(combined.sct, group.by = c("integrated_snn_res.0.2"), split.by = Participant)
DimPlot(combined.sct, group.by = c("integrated_snn_res.0.2"), split.by = AgeGroup)

# Add umap coordinates, extract meta and save for later
umap <- data.frame(combined.sct@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("RPCA_k.anchor.10",  colnames(umap))
combined.sct <- AddMetaData(combined.sct, metadata = umap)
meta <- combined.sct@meta.data
write.csv(meta, 'meta.RPCA_k.anchor.SCT.10.csv')


######Run STACAS #######
rm(list = ls())
Merge <- readRDS("merged.RDS")

start <- Sys.time()
SCT <- replicate(20, "SCT") ##The numner of replicates needs to be how many samples you have. 
ref.list <- SplitObject(Merge, split.by = "orig.ident")
ref.list <- lapply(ref.list, FUN = SCTransform, variable.features.n = 5000)
features.sct <- SelectIntegrationFeatures(ref.list, nfeatures = 5000)

object_integrated <- Merge %>% SplitObject(split.by = "orig.ident") %>%
 Run.STACAS(dims = 1:20, anchor.features = features.sct, assay = SCT) %>%
  RunUMAP(dims = 1:10) 

object_integrated  <- FindNeighbors(object_integrated , reduction = "pca", dims = 1:10)
object_integrated  <- FindClusters(object_integrated , resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2))
object_integrated@meta.data

end <- Sys.time()
end - start # Time difference of 

Participant <- "orig.ident"
AgeGroup <- 'AgeGroup'


# plot data
P1 <- DimPlot(object_integrated, group.by = c(Participant)) 
P2 <- DimPlot(object_integrated, group.by = c(AgeGroup))
P3 <- DimPlot(object_integrated, group.by = c("integrated_snn_res.0.4")) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "STACAS.SCT.Anchor.integrate.dims.20") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("STACAS.SCT.Anchor.integrate.dims.20.png")

DimPlot(object_integrated, split.by = c(Participant), group.by = c("integrated_snn_res.0.4")) 


# Add umap coordinates, extract meta and save for later
umap <- data.frame(object_integrated@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("STACAS",  colnames(umap))
object_integrated <- AddMetaData(object_integrated, metadata = umap)
meta <- object_integrated@meta.data
write.csv(meta, 'meta.STACA.SCT.Anchor.integrate.dims.20.csv')

######Run fastMNN#########
rm(list = ls())
Merge <- readRDS("merged.RDS")

start <- Sys.time()
SCT <- replicate(20, "SCT")
ref.list <- SplitObject(Merge, split.by = "orig.ident")
ref.list <- lapply(ref.list, FUN = SCTransform, variable.features.n = 5000)
features.sct <- SelectIntegrationFeatures(ref.list, nfeatures = 5000)


merge.list <- SplitObject(Merge, split.by="orig.ident")
merge.list <- lapply(X = merge.list, 
                       FUN = SCTransform, 
                       method = "glmGamPoi", 
                       return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 5000)


Merge <- readRDS("seurat_Unadjusted_BioAgedecont.RDS")
ref.list <- SplitObject(Merge, split.by = "orig.ident")
ref.list <- lapply(ref.list, FUN = SCTransform, variable.features.n = 5000)

Merge <- RunFastMNN(object.list = ref.list, features = features.sct, assay ='SCT')
Merge <- RunUMAP(Merge, reduction = "mnn", dims = 1:10)
Merge <- FindNeighbors(Merge, reduction = "mnn", dims = 1:10)
Merge <- FindClusters(Merge)

end <- Sys.time()
end - start # Time difference of 

Participant <- "orig.ident"
AgeGroup <- 'AgeGroup'

# plot data
P1 <- DimPlot(Merge, group.by = c(Participant)) 
P2 <- DimPlot(Merge, group.by = c(AgeGroup))
P3 <- DimPlot(Merge) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "FastMNN.SCT") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("FastMNN.SCT.png")


# Add umap coordinates, extract meta and save for later
umap <- data.frame(Merge@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("FastMNN.SCT",  colnames(umap))
Merge <- AddMetaData(Merge, metadata = umap)
meta <- Merge@meta.data
write.csv(meta, 'meta.FastMNN.SCT.csv')

#####run harmony ######
rm(list = ls())
Merge <- readRDS("merged.RDS")

start <- Sys.time()

merge.list <- SplitObject(Merge, split.by="orig.ident")
merge.list <- lapply(X = merge.list, 
                       FUN = SCTransform, 
                       method = "glmGamPoi", 
                       return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 5000)

merge.sct <- merge(x = merge.list[[1]], y = merge.list[2:length(merge.list)], merge.data=TRUE)
VariableFeatures(merge.sct) <- var.features
merge.sct <- RunPCA(merge.sct, verbose = FALSE)
merge.sct <- RunHarmony(merge.sct, assay.use="SCT", group.by.vars = "orig.ident")
merge.sct <- RunUMAP(merge.sct, reduction = "harmony", dims = 1:10)
merge.sct <- FindNeighbors(merge.sct, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4)) 
end <- Sys.time()
end - start # Time difference of 

Participant <- "orig.ident"
AgeGroup <- 'AgeGroup'

# plot data
P1 <- DimPlot(merge.sct, group.by = c(Participant)) 
P2 <- DimPlot(merge.sct, group.by = c(AgeGroup))
P3 <- DimPlot(merge.sct) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "harmony.SCT") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("harmony.SCT.png")


# Add umap coordinates, extract meta and save for later
umap <- data.frame(merge.sct@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("harmony.SCT",  colnames(umap))
merge.sct <- AddMetaData(merge.sct, metadata = umap)
meta <- merge.sct@meta.data
write.csv(meta, 'meta.harmony.SCT.csv')

saveRDS(merge.sct, 'harmony.sct.RDS')

##Compute batch correction metrics using saved meta.data csv files. Examples below
######unadjusted
####Harmony
df <- read.csv('meta.harmony.SCT.csv')
head(df)

##compute Lisi
cols <- c('harmony.SCTUMAP_1', 'harmony.SCTUMAP_2')
df2 <- df[, colnames(df) %in% cols]

cols <- c('SCT_snn_res.0.3', 'orig.ident')
df3 <- df[, colnames(df) %in% cols]
table(df3$SCT_snn_res.0.3)
table(df3$orig.ident)

res <- compute_lisi(df2, df3, c('SCT_snn_res.0.3', 'orig.ident'))
head(res)
mean(res$orig.ident)

#compute ARI
adj.rand.index(df$SCT_snn_res.0.3, df$orig.ident)

###Hamrony UMAP optimization 
seurat <- readRDS('harmony.sct.RDS')

ElbowPlot(seurat)
DefaultAssay(object = seurat) <- 'SCT'
seurat <- RunUMAP(seurat, reduction = "harmony", min.dist = 0.1, spread = 2,  dims = 1:20)
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)) 


DimPlot(seurat, reduction = "umap", label = TRUE, group.by = 'orig.ident')
DimPlot(seurat, reduction = "umap", group.by = 'AgeGroup', split.by = 'AgeGroup')

##Check cell Types 
##Generic AT makers
AT <- c('PDGFRA',
        'PDGFRB',
        "DCN",
        'DLK1',
        'ZNF423',
        "CD38",
        'RBFOX1',
        'CTNNA2',
        'NRXN3',
        'CDH4',
        "CEBPA",
        "STAT5A",
        'LIPE',
        'LEP',
        'ADIPOQ',
        "DGAT2",
        'PLIN1',
        "PPARG",
        "TBX15",
        "IRS1",
        "HK2",
        'ATP5F1B',
        'GPX4',
        'VWF',
        'CDH5',
        'PECAM1',
        'PTPRC',
        'HLA-DRB1',
        'HLA-DPA1',
        'ITGAX',
        'TREM2',
        'KIT',
        'TPSAB1', 
        'GNLY',
        'NKG7',
        'CD8A'
)



P3 <- DimPlot(seurat, reduction = "umap", label = TRUE, group.by = 'SCT_snn_res.0.25')
P4 <- DotPlot(seurat, features = AT, scale.min = 0, group.by = 'SCT_snn_res.0.25')  +  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))
P3 + P4


###Normalize assays and remove protein coding because no longer needed

DefaultAssay(object = seurat) <- "RNA"
seurat <- NormalizeData(seurat)
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

DefaultAssay(object = seurat) <- "filtered"
seurat <- NormalizeData(seurat)
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat[['Pcoding']] <- NULL

##save final seurat
saveRDS(seurat, "harmony.sct.V2.RDS")





