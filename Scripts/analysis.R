###BioAge data exploration after decontx and harmony integration
library(Seurat)
library(ggplot2)
library(viridis)
library('org.Hs.eg.db')
library('Rgraphviz')
library(tidyverse)
library(ggpubr)
library(rstatix)
library(tidyr)
library(enrichplot)
library(ggnewscale)
library(psycho)
library(parameters)
library(effectsize)
library(RColorBrewer)
library(hypeR)
library("MetBrewer")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(plyr)
library(gridExtra)
library(dplyr)
library(UpSetR)
library("ComplexHeatmap")
library(ggbreak)
library(speckle)
library(limma)
library(plotrix)
library(patchwork)
library(SingleCellExperiment)
library(UpSetR)
library(escape)
library(dittoSeq)
library(data.table)
library(Hmisc)
library(ggradar)
library(ggradar2)
library(fmsb)
library(fgsea)
library(GSEAseuratse)
library(grid)

seurat <- readRDS('harmony.sct.V2.RDS')

Ren <- c("#355828", "#f6b3b0", "#ada43b", "#e69b00", "#bf3729", "#b0799a", "#2f357c" )

Project_Name_ <- "BioAge"

###remove meta.data columns that are no longer needed
cols <- c('SCT_snn_res.1.1', 'SCT_snn_res.0.6', 'SCT_snn_res.0.7', 'SCT_snn_res.0.8', 'SCT_snn_res.0.9', 'seurat_clusters',
          'seurat_clusters')
seurat@meta.data <- seurat@meta.data[, ! colnames(seurat@meta.data) %in% cols]
mean(seurat@meta.data$nFeature_filtered)

##Identify cell types seuratsed on known markers 

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
        'TPSAB1',  #mast 
        'GNLY',
        'NKG7',
        'CD8A'
)

P3 <- DimPlot(seurat, reduction = "umap", label = TRUE, group.by = 'SCT_snn_res.0.25')
P4 <- DotPlot(seurat, features = AT, scale.min = 0, group.by = 'SCT_snn_res.0.25')  +  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))
P3 + P4

seurat@meta.data <- seurat@meta.data %>% mutate(CellType = case_when(
  SCT_snn_res.0.25 == "0"  ~ "Adip_1", 
  SCT_snn_res.0.25 == "1"  ~ "Adip_2", 
  SCT_snn_res.0.25 == "2"  ~ "Pre_Ad",
  SCT_snn_res.0.25 == "3"  ~ "Vascular", 
  SCT_snn_res.0.25 == "4"  ~ "Stem", 
  SCT_snn_res.0.25 == "5"  ~ "Macrophages",
  SCT_snn_res.0.25 == "6"  ~ "Mast"
)) ##repeat for all clusters, the last arguement should not have a column.

seurat@meta.data <- seurat@meta.data %>% mutate(General_Cell = case_when(
  SCT_snn_res.0.25 == "0"  ~ "Adip", 
  SCT_snn_res.0.25 == "1"  ~ "Adip", 
  SCT_snn_res.0.25 == "2"  ~ "Pre_Ad",
  SCT_snn_res.0.25 == "3"  ~ "Vascular", 
  SCT_snn_res.0.25 == "4"  ~ "Stem", 
  SCT_snn_res.0.25 == "5"  ~ "Immune",
  SCT_snn_res.0.25 == "6"  ~ "Immune"
)) ##repeat for all clusters, the last arguement should not have a column.



Idents(seurat) <- seurat$CellType

my_levels2 <- c('Stem', 'Pre_Ad', 'Adip_1', 'Adip_2', 'Vascular', 'Macrophages', 'Mast')
levels(seurat) <- my_levels2

DefaultAssay(object = seurat) <- 'filtered' ##Make sure assay is on filtered

P3 <- DimPlot(seurat, reduction = "umap", label = TRUE)
P4 <- DotPlot(seurat, features = AT, scale.min = 0)  +  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))
P3 + P4

#####Find DEGS for each cluster
Adip_1 <- FindMarkers(seurat, ident.1 = "Adip_1", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(Adip_1, file = paste0(Project_Name_, "Adip_1", '.csv'))

Adip_2 <- FindMarkers(seurat, ident.1 = "Adip_2", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(Adip_2, file = paste0(Project_Name_, "Adip_2", '.csv'))

Pre_Ad <- FindMarkers(seurat, ident.1 = "Pre_Ad", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(Pre_Ad, file = paste0(Project_Name_, "Pre_Ad", '.csv'))

EC <- FindMarkers(seurat, ident.1 = "EC", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(EC, file = paste0(Project_Name_, "EC", '.csv'))

Stem <- FindMarkers(seurat, ident.1 = "Stem", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(Stem, file = paste0(Project_Name_, "Stem", '.csv'))

Macrophages <- FindMarkers(seurat, ident.1 = "Macrophages", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(Macrophages, file = paste0(Project_Name_, "Macrophages", '.csv'))

Mast <- FindMarkers(seurat, ident.1 = "Mast", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(Mast, file = paste0(Project_Name_, "Mast", '.csv'))

###find top upregualted markers for heatmap
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(Top5, n = 35)

top5 <- Top5$gene
str(top5)
top5 <- unique(top5)

P1 <- DimPlot(seurat, reduction = "umap", label = TRUE, cols = Ren, label.color = "white", repel = TRUE, label.box = TRUE) +
  NoLegend() + labs(tag = 'B')  + theme(text = element_text(size=10),
                                     axis.text.x = element_text(colour = "black", size = 10),
                                     axis.text.y = element_text(colour = "black", size = 10))

##Plot top 5 markers in dotplot

Top5 <- c('DCN', 'C3', 'LUM', 'COL1A2', 'C1S', ##Stem
          'RBFOX1', 'CNTNAP2', 'PTPRD', 'CTNNA2', 'NRXN3', ##Pre-AD
         'GPX4', 'PLIN4', "GPX1", "RARRES2", "SAA1", 'PDE3B', 'OGA', 'DDX5', 'WDPCP', 'EBF1', ##Adips
  'VWF', 'BTNL9', 'SPARCL1', 'RBP7', 'RGS5', ##EC
  'MS4A6A', 'HLA-DRA', 'F13A1', 'FRMD4B', 'MRC1', ##Macrophages
  'HDC', 'KIT', 'TPSB2', 'MS4A2', 'TPSAB1'#Mast cells
 )

P2 <- DotPlot(seurat, features = Top5, scale.min = 0)  +  coord_flip() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10)) +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG'))) + labs(tag = 'C')

###Add additional meta data
##Add sex
seurat@meta.data <- seurat@meta.data %>% mutate(Sex = case_when(
  ID == "B01" ~ "F",
  ID == "B06" ~ "M", 
  ID == "B07" ~ "F",
  ID == "B08" ~ "F",
  ID == "B09" ~ "F",
  ID == "B10" ~ "F",
  ID == "B14" ~ "M",
  ID == "B15" ~ "F",
  ID == "B16" ~ "M",
  ID == "B17" ~ "F",
  ID == "B18" ~ "F",
  ID == "B19" ~ "M",
  ID == "B20" ~ "F",
  ID == "B22" ~ "F",
  ID == "B23" ~ "M",
  ID == "B24" ~ "M",
  ID == "B26" ~ "M",
  ID == "B27" ~ "M",
  ID == "B30" ~ "M",
  ID == "B31" ~ "M",
))

seurat@meta.data$Age_Sex <- paste(seurat@meta.data$AgeGroup, sep = '_', seurat@meta.data$Sex)

##Add BMI
seurat@meta.data <- seurat@meta.data %>% mutate(BMI = case_when(
  ID == "B01" ~ "20.9",
  ID == "B06" ~ "34.5", 
  ID == "B07" ~ "29",
  ID == "B08" ~ "35.3",
  ID == "B09" ~ "35.6",
  ID == "B10" ~ "32.9",
  ID == "B14" ~ "25.4",
  ID == "B15" ~ "29.7",
  ID == "B16" ~ "24",
  ID == "B17" ~ "26.4",
  ID == "B18" ~ "27.2",
  ID == "B19" ~ "29.4",
  ID == "B20" ~ "22.4",
  ID == "B22" ~ "34.1",
  ID == "B23" ~ "24.9",
  ID == "B24" ~ "35.1",
  ID == "B26" ~ "32.2",
  ID == "B27" ~ "27.8",
  ID == "B30" ~ "22.9",
  ID == "B31" ~ "27.3",
))

##Add BMI category meta.data
seurat@meta.data <- seurat@meta.data %>% mutate(BMI_category = case_when(
  BMI < 25  ~ "lean",
  BMI < 30 & BMI > 25  ~ "overweight",
  BMI > 30  ~ "obese",
))

##Add Waist
seurat@meta.data <- seurat@meta.data %>% mutate(Waist = case_when(
  ID == "B01" ~ "84.55",
  ID == "B06" ~ "106.5", 
  ID == "B07" ~ "99.25",
  ID == "B08" ~ "93.25",
  ID == "B09" ~ "119.1",
  ID == "B10" ~ "98.5",
  ID == "B14" ~ "104.25",
  ID == "B15" ~ "86.25",
  ID == "B16" ~ "83.75",
  ID == "B17" ~ "85.25",
  ID == "B18" ~ "91",
  ID == "B19" ~ "95.1",
  ID == "B20" ~ "82.5",
  ID == "B22" ~ "89.85",
  ID == "B23" ~ "80",
  ID == "B24" ~ "119.85",
  ID == "B26" ~ "117.55",
  ID == "B27" ~ "106.25",
  ID == "B30" ~ "91.2",
  ID == "B31" ~ "101",
))

##Add WHR
seurat@meta.data <- seurat@meta.data %>% mutate(WHR = case_when(
  ID == "B01" ~ "0.91",
  ID == "B06" ~ "0.9", 
  ID == "B07" ~ "0.93",
  ID == "B08" ~ "0.82",
  ID == "B09" ~ "0.94",
  ID == "B10" ~ "0.88",
  ID == "B14" ~ "1",
  ID == "B15" ~ "0.77",
  ID == "B16" ~ "0.89",
  ID == "B17" ~ "0.82",
  ID == "B18" ~ "0.91",
  ID == "B19" ~ "0.98",
  ID == "B20" ~ "0.84",
  ID == "B22" ~ "0.77",
  ID == "B23" ~ "0.9",
  ID == "B24" ~ "1.05",
  ID == "B26" ~ "1.08",
  ID == "B27" ~ "1.04",
  ID == "B30" ~ "0.7",
  ID == "B31" ~ "1.1"
))


###Cell Type Composition
tab2 <- prop.table(table(seurat$CellType, seurat$orig.ident), margin = 2)
tab2 <- as.data.frame(tab2)
colnames(tab2) <- c('Cell_Type', 'Participant', 'Freq')
tab2

tab2 <- tab2 %>% mutate(AgeGroup = case_when(
  Participant == "B01" ~ "Older",
  Participant == "B06" ~ "Younger", 
  Participant == "B07" ~ "Older",
  Participant == "B08" ~ "Younger",
  Participant == "B09" ~ "Older",
  Participant == "B10" ~ "Younger",
  Participant == "B14" ~ "Older",
  Participant == "B15" ~ "Younger",
  Participant == "B16" ~ "Younger",
  Participant == "B17" ~ "Older",
  Participant == "B18" ~ "Younger",
  Participant == "B19" ~ "Younger",
  Participant == "B20" ~ "Younger",
  Participant == "B22" ~ "Older",
  Participant == "B23" ~ "Younger",
  Participant == "B24" ~ "Older",
  Participant == "B26" ~ "Older",
  Participant == "B27" ~ "Older",
  Participant == "B30" ~ "Younger",
  Participant == "B31" ~ "Older",
))

propeller(
  x = seurat,
  clusters = seurat@meta.data$CellType,
  sample = seurat@meta.data$ID,
  group = seurat@meta.data$AgeGroup,
  trend = FALSE,
  robust = FALSE,
  transform = "logit"
)

propeller(
  x = seurat,
  clusters = seurat@meta.data$CellType,
  sample = seurat@meta.data$ID,
  group = seurat@meta.data$Sex,
  trend = FALSE,
  robust = FALSE,
  transform = "logit"
)



###Make graph in ggplot2 with 3 different layers
my_dat <- dplyr::summarize(group_by(tab2, AgeGroup, Cell_Type),
                           my_mean = mean(Freq, na.rm = TRUE),
                           my_se = std.error(Freq, na.rm = TRUE))

tab2$Cell_Type <- factor(tab2$Cell_Type, levels = my_levels2)
my_dat$Cell_Type <- factor(my_dat$Cell_Type, levels = my_levels2)

levels <- c('Younger', 'Older')

tab2$AgeGroup <- factor(tab2$AgeGroup, levels = levels)
my_dat$AgeGroup <- factor(my_dat$AgeGroup, levels = levels)

P3 <- ggplot() + 
  geom_seuratr(data = my_dat,
           aes(y = my_mean, x = Cell_Type, fill = AgeGroup, 
               ymin = my_mean - my_se,
               ymax = my_mean + my_se), position = "dodge", stat="identity", width=0.75) + 
  scale_fill_manual(values=c("#87CEFA", "#EC7014"))+ 
  theme_classic() +
  geom_errorseuratr(data = my_dat,
                aes(y = my_mean,  x = Cell_Type, fill = AgeGroup, 
                    ymin = my_mean - my_se,
                    ymax = my_mean + my_se), stat="identity",  position = "dodge", width=0.75) + 
  geom_jitter(data = tab2, aes(Cell_Type, Freq, fill = AgeGroup), shape = 21, color = "black",
              position = position_jitterdodge(jitter.height = 0, jitter.width = .2)) +
labs(y = "% of cells ", x = "", fill = "Age Group", tag = "D") +
  ylim(0, 0.4) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10))  +
  annotate("segment", x = 5.75, xend = 6.25, y = 0.15, yend = 0.15,
           colour = "black", size = 1) +
  annotate("text", label = 'P = 0.059', x = 6, y = 0.175, size =3,
           colour = "black") +
  annotate("segment", x = 6.75, xend = 7.25, y = 0.05, yend = 0.05,
           colour = "black", size = 1) +
  annotate("text", label = '**', x = 7, y = 0.06, size =8,
           colour = "black") +
  theme(legend.position = "top") 


###Figure 1
P1 + P2 + P3 +
  plot_layout(design = "112
             332")

###plot proportions for sex
tab3 <- prop.table(table(seurat$CellType, seurat$orig.ident), margin = 2)
tab3 <- as.data.frame(tab3)
colnames(tab3) <- c('Cell_Type', 'Participant', 'Freq')
tab3

tab3 <- tab3 %>% mutate(Sex = case_when(
  Participant == "B01" ~ "F",
  Participant == "B06" ~ "M", 
  Participant == "B07" ~ "F",
  Participant == "B08" ~ "F",
  Participant == "B09"~ "F",
  Participant == "B10" ~ "F",
  Participant == "B14" ~ "M",
  Participant == "B15" ~ "F",
  Participant == "B16" ~ "M",
  Participant == "B17" ~ "F",
  Participant == "B18" ~ "F",
  Participant == "B19" ~ "M",
  Participant == "B20" ~ "F",
  Participant == "B22" ~ "F",
  Participant == "B23" ~ "M",
  Participant == "B24" ~ "M",
  Participant == "B26" ~ "M",
  Participant == "B27"~ "M",
  Participant == "B30" ~ "M",
  Participant == "B31" ~ "M",
))


my_dat <- dplyr::summarize(group_by(tab3, Sex, Cell_Type),
                           my_mean = mean(Freq, na.rm = TRUE),
                           my_se = std.error(Freq, na.rm = TRUE))

tab3$Cell_Type <- factor(tab3$Cell_Type, levels = my_levels2)
my_dat$Cell_Type <- factor(my_dat$Cell_Type, levels = my_levels2)

ggplot() + 
  geom_seuratr(data = my_dat,
           aes(y = my_mean, x = Cell_Type, fill = Sex, 
               ymin = my_mean - my_se,
               ymax = my_mean + my_se), position = "dodge", stat="identity", width=0.75) + 
  scale_fill_manual(values=c("#AA4499", "#44AA99"))+ 
  theme_classic() +
  geom_errorseuratr(data = my_dat,
                aes(y = my_mean,  x = Cell_Type, fill = Sex, 
                    ymin = my_mean - my_se,
                    ymax = my_mean + my_se), stat="identity",  position = "dodge", width=0.75) + 
  geom_jitter(data = tab3, aes(Cell_Type, Freq, fill = Sex), shape = 21, color = "black",
              position = position_jitterdodge(jitter.height = 0, jitter.width = .2)) +
  labs(y = "% of cells ", x = "", fill = "Sex") +
  ylim(0, 0.4) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10))


##plot individual umaps
seurat@meta.data <- seurat@meta.data %>% mutate(manu_ID = case_when(
  ID == "B01" ~ "Older Female 1",
  ID == "B06" ~ "Younger Male 1",
  ID == "B07" ~ "Older Female 2",
  ID == "B08" ~ "Younger Female 1",
  ID == "B09" ~ "Older Female 3",
  ID == "B10" ~ "Younger Female 2",
  ID == "B14" ~ "Older Male 1",
  ID == "B15" ~ "Younger Female 3",
  ID == "B16" ~ "Younger Male 2",
  ID == "B17" ~ "Older Female 4",
  ID == "B18" ~ "Younger Female 4",
  ID == "B19" ~ "Younger Male 3",
  ID == "B20" ~ "Younger Female 5",
  ID == "B22" ~ "Older Female 5",
  ID == "B23" ~ "Younger Male 4",
  ID == "B24" ~ "Older Male 2",
  ID == "B26" ~ "Older Male 3",
  ID == "B27" ~ "Older Male 4",
  ID == "B30" ~ "Younger Male 5",
  ID == "B31" ~ "Older Male 5"
))

Idents(seurat) <- seurat$manu_ID

levels <- c("Younger Female 1", "Younger Female 2", "Younger Female 3", "Younger Female 4", "Younger Female 5", "Younger Male 1", "Younger Male 2",
            "Younger Male 3", "Younger Male 4", "Younger Male 5", "Older Female 1", "Older Female 2", "Older Female 3", "Older Female 4", "Older Female 5",
            "Older Male 1", "Older Male 2", "Older Male 3", "Older Male 4", "Older Male 5")
levels(seurat) <- levels
Idents(seurat) <- seurat$CellType
levels(seurat) <- my_levels2
DimPlot(seurat, pt.size = 0.5, label = FALSE,  split.by = "manu_ID", reduction = "umap", cols = Ren, ncol = 5) +
  theme(text = element_text(size=10))

DimPlot(seurat, pt.size = 0.5, label = FALSE,  group.by = "manu_ID", reduction = "umap") +
  theme(text = element_text(size=10))

###DEGS between older and younger for main cell type

seurat@meta.data$Cell_Group <- paste0(seurat@meta.data$CellType,  sep = '_', seurat@meta.data$AgeGroup) ##use the column that your group data is stored under
seurat@meta.data
Idents(seurat) = seurat$Cell_Group

##Now we can run find markers but use ident.1 and ident.2 to define the direct comparisons.
####Example
Stem_group_differences <- FindMarkers(seurat, ##seurat object
                                    ident.1 = "Stem_Older", ##Cluster and first group #example
                                    ident.2 = "Stem_Younger", ##should match the first cluster but be the opposite group
                                    min.pct = 0.1, ##only test genes that are detected in a 10% of cells in either of the two populations. 
                                    logfc.threshold = 0.1, ##only include genes with a minimum logfc of 0.1
                                    only.pos = FALSE ##include both up and downregulated genes
)
write.csv(Stem_group_differences, file = paste0(Project_Name_, "Stem_group_differences", '.csv'))

Pre_Ad_group_differences <- FindMarkers(seurat, ident.1 = "Pre_Ad_Older",  ident.2 = "Pre_Ad_Younger",  min.pct = 0.1,  logfc.threshold = 0.1, only.pos = FALSE)
write.csv(Pre_Ad_group_differences, file = paste0(Project_Name_, "Pre_Ad_group_differences", '.csv'))

Adip_1_group_differences <- FindMarkers(seurat, ident.1 = "Adip_1_Older",  ident.2 = "Adip_1_Younger",  min.pct = 0.1,  logfc.threshold = 0.1, only.pos = FALSE)
write.csv(Adip_1_group_differences, file = paste0(Project_Name_, "Adip_1_group_differences", '.csv'))

Adip_2_group_differences <- FindMarkers(seurat, ident.1 = "Adip_2_Older",  ident.2 = "Adip_2_Younger",  min.pct = 0.1,  logfc.threshold = 0.1, only.pos = FALSE)
write.csv(Adip_2_group_differences, file = paste0(Project_Name_, "Adip_2_group_differences", '.csv'))

EC_group_differences <- FindMarkers(seurat, ident.1 = "EC_Older",  ident.2 = "EC_Younger",  min.pct = 0.1,  logfc.threshold = 0.1, only.pos = FALSE)
write.csv(EC_group_differences, file = paste0(Project_Name_, "EC_group_differences", '.csv'))

Macrophages_group_differences <- FindMarkers(seurat, ident.1 = "Macrophages_Older",  ident.2 = "Macrophages_Younger",  min.pct = 0.1,  logfc.threshold = 0.1, only.pos = FALSE)
write.csv(Macrophages_group_differences, file = paste0(Project_Name_, "Macrophages_group_differences", '.csv'))

Mast_group_differences <- FindMarkers(seurat, ident.1 = "Mast_Older",  ident.2 = "Mast_Younger",  min.pct = 0.1,  logfc.threshold = 0.1, only.pos = FALSE)
write.csv(Mast_group_differences, file = paste0(Project_Name_, "Mast_group_differences", '.csv'))


###Reload DEgs
Stem <- read.csv("Old_V_Young_DEG/Stem_group_differences.csv", row.names = 1)
Pre_Ad <- read.csv("Old_V_Young_DEG/Pre_Ad_group_differences.csv", row.names = 1)
Adip_1 <- read.csv("Old_V_Young_DEG/Adip_1_group_differences.csv", row.names = 1)
Adip_2 <- read.csv("Old_V_Young_DEG/Adip_2_group_differences.csv", row.names = 1)
EC <- read.csv("Old_V_Young_DEG/EC_group_differences.csv", row.names = 1)
Macrophages <- read.csv("Old_V_Young_DEG/Macrophages_group_differences.csv", row.names = 1)
Mast <- read.csv("Old_V_Young_DEG/Mast_group_differences.csv", row.names = 1)


##save background genes lists
backgroundStem <- rownames(Stem)
backgroundPre_Ad <- rownames(Pre_Ad)
backgroundAdip_1 <- rownames(Adip_1)
backgroundAdip_2 <- rownames(Adip_2)
backgroundEC <- rownames(EC)
backgroundMacrophages <- rownames(Macrophages)
backgroundMast <- rownames(Mast)


upregulated <- function(df, logfc){
  df <- df  %>% dplyr::select(2, 5)
  df <- subset(df, df$p_val_adj < 0.05 & df$avg_log2FC > logfc)
  df <- df[order(-df$avg_log2FC),] 
  return(df)
}

downregulated <- function(df, logfc){
  df <- df  %>% dplyr::select(2, 5)
  df <- subset(df, df$p_val_adj < 0.05 & df$avg_log2FC < logfc)
  df <- df[order(-df$avg_log2FC),] 
  return(df)
}

###seperate lists of upregulated and downregulated genes for each cluster
Stem_Up <- upregulated(df = Stem, logfc = 0.25)
Stem_Down <- downregulated(df = Stem, logfc = -0.25)

Pre_Ad_Up <- upregulated(df = Pre_Ad, logfc = 0.25)
Pre_Ad_Down <- downregulated(df = Pre_Ad, logfc = -0.25)

Adip_1_Up <- upregulated(df = Adip_1, logfc = 0.25)
Adip_1_Down <- downregulated(df = Adip_1, logfc = -0.25)

Adip_2_Up <- upregulated(df = Adip_2, logfc = 0.25)
Adip_2_Down <- downregulated(df = Adip_2, logfc = -0.25)

EC_Up <- upregulated(df = EC, logfc = 0.25)
EC_Down <- downregulated(df = EC, logfc = -0.25)

Macrophages_Up <- upregulated(df = Macrophages, logfc = 0.25)
Macrophages_Down <- downregulated(df = Macrophages, logfc = -0.25)

Mast_Up <- upregulated(df = Mast, logfc = 0.25)
Mast_Down <- downregulated(df = Mast, logfc = -0.25)


##Get list of upregulated genes for each cluster 
Stem_Up_genes <- rownames(Stem_Up)
Stem_Down_genes <- rownames(Stem_Down)

Pre_Ad_Up_genes <- rownames(Pre_Ad_Up)
Pre_Ad_Down_genes <- rownames(Pre_Ad_Down)

Adip_1_Up_genes <- rownames(Adip_1_Up)
Adip_1_Down_genes <- rownames(Adip_1_Down)

Adip_2_Up_genes <- rownames(Adip_2_Up)
Adip_2_Down_genes <- rownames(Adip_2_Down)

EC_Up_genes <- rownames(EC_Up)
EC_Down_genes <- rownames(EC_Down)

Macrophages_Up_genes <- rownames(Macrophages_Up)
Macrophages_Down_genes <- rownames(Macrophages_Down)

Mast_Up_genes <- rownames(Mast_Up)
Mast_Down_genes <- rownames(Mast_Down)


signatures_up <- list("Stem" = Stem_Up_genes,
                      "Pre_Ad" = Pre_Ad_Up_genes,
                      "Adip_1" = Adip_1_Up_genes,
                      "Adip_2" = Adip_2_Up_genes,
                      "Vascular" = EC_Up_genes,
                      "Macrophages" = Macrophages_Up_genes,
                      "Mast" = Mast_Up_genes
)

m <-  make_comb_mat(list_to_matrix(signatures_up))
UpSet(m)


UpSet(m[comb_size(m) >= 5], comb_col = "#EC7014")

##652 479
upset(fromList(signatures_up), nintersects = 15,
      order.by="freq", matrix.color="#EC7014", point.size=5, main.seuratr.color = "#EC7014",
      sets.seuratr.color=c( "#e69b00", "#f6b3b0",  "#ada43b", "#355828",  "#bf3729"), 
      sets.x.label = 'No. of DEGs',
      text.scale = 1.5) 


signatures_down <- list("Stem" = Stem_Down_genes,
                      "Pre_Ad" = Pre_Ad_Down_genes,
                      "Adip_1" = Adip_1_Down_genes,
                      "Adip_2" = Adip_2_Down_genes,
                      "Vascular" = EC_Down_genes,
                      "Macrophages" = Macrophages_Down_genes,
                      "Mast" = Mast_Down_genes
)

m2 <-  make_comb_mat(list_to_matrix(signatures_down))
UpSet(m2)

##652 479
upset(fromList(signatures_down), nintersects = 15,
      order.by="freq", matrix.color="#87CEFA", point.size=5, main.seuratr.color = "#87CEFA",
      sets.seuratr.color=c( "#e69b00",  "#bf3729", "#ada43b", "#355828",  "#f6b3b0"), 
      sets.x.label = 'No. of DEGs',
      text.scale = 1.5) 


###run Pathway analysis with HypeR
##Get genes sets
Hallmark <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
Reactome <- msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean=TRUE)

##run pathway individually for each cluster and save excel file for further probing
run_pathway <- function(topgenes, background, geneset, geneset2, fdr, clustername){
  hyp <- hypeR(topgenes, geneset, test="hypergeometric", background=background, absolute = FALSE, pval = 0.05, fdr = fdr)
  hyp_to_excel(hyp, file_path=paste0(clustername, geneset2, ".xlsx"))
  return(hyp)
}

###Run pathway for each cluster with upregulated (up in Older) and downregulated (up in Younger) for Hallmark and Reactome
#e.g.
path <- run_pathway(topgenes = Adip_1_Up_genes, ##Upregulated genes that have been filtered
                    background = backgroundAdip_1, ##background should be specific to the cluster
                    geneset = Hallmark, ##Hallmark, Reactome
                    geneset2 = 'Hallmark',
                    fdr = 0.05, 
                    clustername = 'Adip_1_Older'
)



##combine pathways and make a dotplot of top pathways in the Older group
combine_pathways <- function(reactome, hallmark, cluster_name){
  reactome$data$geneset <- 'Reactome'
  hallmark$data$geneset <- 'Hallmark'
  mega <- rbind(reactome$data,  hallmark$data)
  mega$CellType <- cluster_name
  mega <- mega[order(mega$fdr),]
  return(mega)
}

Adip1 <- combine_pathways(reactome = Adip1R, hallmark = Adip1H, cluster_name = 'Adip_1')
str(Adip1)
Adip1_top <- Adip1$label[1:8]

Adip2 <- combine_pathways(reactome = Adip2R, hallmark = Adip2H, cluster_name = 'Adip_2')
str(Adip2)
##remove irrelevant pathways 
labs <- c("Sensory Processing Of Sound By Outer Hair Cells Of The Cochlea", "Sensory Processing Of Sound")
Adip2<- Adip2[!Adip2$label %in% labs,]
Adip2_top <- Adip2$label[1:8]

###PreAd only had upregulated pathways in hallmark, organize seperately
Pre_AdH$data$geneset = 'Hallmark'
Pre_AdH$data$CellType = 'Pre_Ad'
Pre_Ad <- Pre_AdH$data
Pre_Ad <- Pre_Ad[order(Pre_Ad$fdr),]
Pre_Ad_top <- Pre_Ad$label[1:2]

###Stem only had upregulated pathways in Reactome, organize seperately
StemR$data$geneset = 'Reactome'
StemR$data$CellType = 'Stem'
Stem <- StemR$data
Stem <- Stem[order(Stem$fdr),]
Stem_top <- Stem$label[1:8]

Vascular <- combine_pathways2(reactome = VascularR, hallmark = VascularH, cluster_name = 'Vascular')
str(Vascular)
Vascular_top <- Vascular$label[1:4]

top_pathways <- unique(c(Stem_top, Pre_Ad_top, Adip1_top, Adip2_top,   Vascular_top))

#See if top pathway are significant in other clusters
Adip_1_subset <- Adip1[Adip1$label %in% top_pathways,]
Adip_2_subset <- Adip2[Adip2$label %in% top_pathways,]
Pre_Ad_subset <-  Pre_Ad[Pre_Ad$label %in% top_pathways,]
Stem_subset <-  Stem[Stem$label %in% top_pathways,]
Vascular_subset <-  Vascular[Vascular$label %in% top_pathways,]

#make dotplot of top pathways upregualted in older
mega <- rbind(Adip_1_subset, Adip_2_subset, Pre_Ad_subset, Stem_subset, Vascular_subset)
head(mega)

cols <- c("label", 'fdr', 'overlap', "CellType")

# extracting data frame cols
mod1 <- mega[, colnames(mega) %in% cols]
mod1$log10FDR <- -log10(mod1$fdr)
head(mod1)

colnames(mod1) <- c('Pathway', 'FDR', 'Count', 'CellType', 'log10FDR')

level3 <- c('Stem','Pre_Ad', 'Adip_1', 'Adip_2', 'Vascular')
mod1$CellType <- factor(mod1$CellType, levels= level3)
mod1$Pathway <- factor(mod1$Pathway, levels= top_pathways)

##647 542
ggplot(mod1, aes(x = CellType, Pathway, size = Count)) +
  geom_point(aes(colour = -log10FDR)) +
  scale_colour_gradient(high = "grey", low = "#EC7014", n.breaks = 5) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 40))



##combine pathways and make a dotplot of top pathways in the Younger group

Adip2d <- combine_pathways2(reactome = Adip2Rd, hallmark = Adip2Hd, cluster_name = 'Adip_2')
str(Adip2d)
Adip2d<- Adip2d[!Adip2d$label %in% labs,]
Adip2d_top <- Adip2d$label[1:8]

###No pathways upregulated in older for Pre_Ad 

Vasculard <- combine_pathways2(reactome = VascularRd, hallmark = VascularHd, cluster_name = 'Vascular')
str(Vasculard)
Vasculard_top <- Vasculard$label[1:8]

Stemd <- combine_pathways2(reactome = StemRd, hallmark = StemHd, cluster_name = 'Stem')
str(Stemd)
Stemd_top <- Stemd$label[1:8]


top_dpathways <- unique(c(Stemd_top, Adip2d_top, Vasculard_top))

#Pull out pathways from clusters and see if if they are significant.
Adip_2_subsetd <- Adip2d[Adip2d$label %in% top_dpathways,]
Stem_subsetd <-  Stemd[Stemd$label %in% top_dpathways,]
Vascular_subsetd <-  Vasculard[Vasculard$label %in% top_dpathways,]

#make big dotplot heatmap
mega <- rbind(Adip_2_subsetd, Stem_subsetd, Vascular_subsetd)
head(mega)

cols <- c("label", 'fdr', 'overlap', "CellType")

# extracting data frame cols
mod1 <- mega[, colnames(mega) %in% cols]
mod1$log10FDR <- -log10(mod1$fdr)
head(mod1)

colnames(mod1) <- c('Pathway', 'FDR', 'Count', 'CellType', 'log10FDR')

level3 <- c('Stem', 'Adip_2', 'Vascular')
mod1$CellType <- factor(mod1$CellType, levels= level3)
mod1$Pathway <- factor(mod1$Pathway, levels= top_dpathways)

ggplot(mod1, aes(x = CellType, Pathway, size = Count)) +
  geom_point(aes(colour = -log10FDR)) +
  scale_colour_gradient(high = "grey", low = "#87CEFA", n.breaks = 5) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 40))


###Make logFC heatmap of top contributing genes for each pathway upregulated in each cell type
Older <- c('RPS2', 'RPL28', 'RPL8', 'SLIT3', 'RPS19', 'RPS8', 'NCL', 'NNMT', 'RPS24', "RPL13A", ##stem
           "C1QA", "C3", "CTSD", "GNAI2", "LCK", "LGMN", "PIK3R5", "SH2B3", "CFD", "GSN", "MMP2", #Pre_Ad
            "C1R", "C1S", "CFH", "CTSS", "TIMP1", "TIMP2", "DCN", "EMP1", "IFITM3", "TGFBR3",  #Adip_1
           "CALD1", "COL6A3", "DAB2", "FBLN1", "FBN1", "FLNA", "FSTL1", "ITGB5", "LGALS1", "MFAP5", "MGP", "MYL9",   #Adip_1
           "PRRX1", "SPARC", "SPP1", "TAGLN",  #Adip_1
           "ADAM12", "COL12A1", "COL1A2", "COL3A1", "COL4A2", "COL5A1", "COL5A2", "COL6A2", "CXCL12", #Adip_2
          "ITGB5"," PCOLCE2", "TPM2", "VIM", "CLU", "GPX3", "IGFBP6", "LMNA", "SPTAN1", #Adip_2
          "DOCK4", 'DOCK9', "CBLB",  "CTSB", "LRP1", "LAMA2"  #vascular
)
           

Younger <- c('LUM', "IGFBP3", "COL15A1", "COL3A1", "FN1", "SPARCL1", 'SCD', "COL1A1", "COL1A2", "COL4A2", ##Stem
             "ACACA", "ACER3", "ACLY", "ACSL1", "DBI", "ELOVL5", "ELOVL6", "FASN", "G0S2", "GPAM", 'GLUL', "MGLL", ##Adip_2
             "PLA2G2A", "PNPLA3", "PPARG", 'PRDX6', 'SCD', "THRSP", "TSPO", ##Adip_2
             "RPL23", "RPL27", "RPL31", "RPL34", "RPL35A", "RPL37A", "RPL41", "RPL7", "RPL9", "RPS13", "RPS20", 
             "RPS25", "RPS26", "RPS27", "RPS29", "RPS3A", "Useurat52" #vascular 
             )


topDEGS <- unique(c(Younger, Older))
str(topDEGS)

df <- as.data.frame(topDEGS)
colnames(df) <- "gene"

Stem$gene <- rownames(Stem)
Stem <- subset(Stem, Stem$p_val_adj < 0.05)
Stem2 <- Stem[Stem$gene %in% topDEGS,]
Stem2 <- Stem2 %>% select(2, 6)
names(Stem2)[names(Stem2) == "avg_log2FC"] <- "Stem"
df <- left_join(df, Stem2, by = 'gene')

Pre_Ad$gene <- rownames(Pre_Ad)
Pre_Ad <- subset(Pre_Ad, Pre_Ad$p_val_adj < 0.05)
Pre_Ad2 <- Pre_Ad[Pre_Ad$gene %in% topDEGS,]
Pre_Ad2 <- Pre_Ad2 %>% select(2, 6)
names(Pre_Ad2)[names(Pre_Ad2) == "avg_log2FC"] <- "Pre_Ad"
df <- left_join(df, Pre_Ad2, by = 'gene')

Adip_1$gene <- rownames(Adip_1)
Adip_1 <- subset(Adip_1, Adip_1$p_val_adj < 0.05)
Adip_12 <- Adip_1[Adip_1$gene %in% topDEGS,]
Adip_12 <- Adip_12 %>% select(2, 6)
names(Adip_12)[names(Adip_12) == "avg_log2FC"] <- "Adip_1"
df <- left_join(df, Adip_12, by = 'gene')
head(df)

Adip_2$gene <- rownames(Adip_2)
Adip_2 <- subset(Adip_2, Adip_2$p_val_adj < 0.05)
Adip_22 <- Adip_2[Adip_2$gene %in% topDEGS,]
Adip_22 <- Adip_22 %>% select(2, 6)
names(Adip_22)[names(Adip_22) == "avg_log2FC"] <- "Adip_2"
df <- left_join(df, Adip_22, by = 'gene')
head(df)

EC$gene <- rownames(EC)
EC <- subset(EC, EC$p_val_adj < 0.05)
EC2 <- EC[EC$gene %in% topDEGS,]
EC2 <- EC2 %>% select(2, 6)
names(EC2)[names(EC2) == "avg_log2FC"] <- "Vascular"
df <- left_join(df, EC2, by = 'gene')
head(df)



#df[is.na(df)] <- 0

df.long <- pivot_longer(data = df, cols = -gene, names_to = "Cell_Type", values_to = "avg_log2FC")
head(df.long)

mylevs <- c('Stem', 'Pre_Ad', 'Adip_1', 'Adip_2', 'Vascular')

df.long$gene <- factor(df.long$gene, levels = topDEGS)
df.long$Cell_Type <- factor(df.long$Cell_Type, levels = mylevs)
min(df.long$avg_log2FC, na.rm = TRUE)
max(df.long$avg_log2FC, na.rm = TRUE)

df.long2 <- as.data.frame(df.long)

youn <- df.long2[df.long2$gene %in% Younger,]
old <- df.long2[df.long2$gene %in% Older,]

min(old$avg_log2FC, na.rm = TRUE)
max(old$avg_log2FC, na.rm = TRUE)

ggplot(data = youn, mapping = aes(x = gene, y = Cell_Type, fill = avg_log2FC)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colours = c("#87CEFA", "white", "#EC7014"),
                       values = scales::rescale(c(-0.9,  0,  0.4))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_fixed() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggplot(data = old, mapping = aes(x = gene, y = Cell_Type, fill = avg_log2FC)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colours = c("#87CEFA", "white", "#EC7014"),
                       values = scales::rescale(c(-0.9,  0,  1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_fixed() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


#####Run senescence analysis using senmayo gene list
senmayo <- read.csv('Senmayo.csv')


####run GSEA on DEG list between older and younger for each participant
##Cutomize enrichment plot
plotEnrichment <- function(pathway, stats,
                           gseaParam=1,
                           ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway) ###Sorts which genes belong to the pathway
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj) ##number of genes tested
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="darkgreen", size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color="darkgreen") + theme_bw() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize) +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    labs(x="rank", y="enrichment score") +
    theme(text = element_text(size=10),
          axis.text.x = element_text(colour = "black", size = 10),
          axis.text.y = element_text(colour = "black", size = 10),
          axis.title.y = element_text(colour = "black", size = 15),
          axis.title.x = element_text(colour = "black", size = 15))
  g
}


head(Adip_1)
df <- Adip_1

###Organize data
transform2 <- function(df){
  df$log10p <- -log10(df$p_val)
  df$rank <- df$log10p * df$avg_log2FC
  df <- df[order(-df$rank),] 
  df$Symbol <- rownames(df)
  df <- df  %>% dplyr::select(8, 7)
  names(df)[names(df) == "Symbol"] <- "SYMBOL"
  return(df)
}

set.seed(123)
res <- transform2(Adip_1)
ranks <- deframe(res)
head(ranks, 20)

Sen_gene <- senmayo$Gene

sengene <- list('SenMayo' = Sen_gene)
fgseaRes <- fgsea(pathways = sengene, 
                  stats    = ranks,
                  minSize  = 10,
                  maxSize  = 500,
                  gseaParam = 1
)

fgseaRes 

set.seed(123)
P <- plotEnrichment(pathway = sengene[["SenMayo"]],
               stats = ranks,
               gseaParam = 1) + labs(title="SenMayo Adip_1")
# Create a text
grob <- grobTree(textGrob("P = 0.042", x=0.8,  y=0.9,
                          gp=gpar(col="black", fontsize=13)))
grob2 <- grobTree(textGrob("NES = 1.53", x=0.8,  y=0.85,
                          gp=gpar(col="black", fontsize=13)))

# Plot
P <- P + annotation_custom(grob) + annotation_custom(grob2)

##Repeat for other cell types

#################################
###Run correlation analysis between CXCL14 and differentially expressed genes for each cell type
###Get average expression of filtered genes per participant and cell type
seurat@meta.data$Cell_ID <- paste(seurat@meta.data$CellType, sep = '_', seurat@meta.data$ID)

Older_DEGS <- c(Stem_Up_genes, Pre_Ad_Up_genes, Adip_1_Up_genes, Adip_2_Up_genes, EC_Up_genes, Macrophages_Up_genes)
Younger_DEGS <- c(Stem_Down_genes, Pre_Ad_Down_genes, Adip_1_Down_genes, Adip_2_Down_genes, EC_Down_genes, Macrophages_Down_genes)

genes <- unique(c(Younger_DEGS, Older_DEGS))

av <- AverageExpression(seurat, assays = 'filtered', group.by = 'Cell_ID', slot = 'data', features = genes)
av <- as.data.frame(av)

##Adip1 
Adip1 <- av %>% select(1:20)
colnames(Adip1) <- str_remove_all(colnames(Adip1), "filtered.Adip_1_")
Adip1 <- t(as.data.frame(Adip1))

cor <- rcorr(as.matrix(Adip1), type = 'pearson')

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

mat <- flattenCorrMatrix(cor$r, cor$P)

head(mat)
mat <- mat[order(mat$p),] 
head(mat)
min(mat$p)

mat$Pvalue_adjust <- p.adjust(mat$p, method="bonferroni")
significant <- subset(mat, p < 0.05)

CXCL14_Adip1 <- subset(significant, row == 'CXCL14')
str(CXCL14_Adip1)

##Adip2 
Adip2 <- av %>% select(21:40)
colnames(Adip2) <- str_remove_all(colnames(Adip2), "filtered.Adip_2_")
Adip2 <- t(as.data.frame(Adip2))

cor <- rcorr(as.matrix(Adip2), type = 'pearson')
mat <- flattenCorrMatrix(cor$r, cor$P)

mat$Pvalue_adjust <- p.adjust(mat$p, method="bonferroni")
significant <- subset(mat, p < 0.05)
str(significant)
CXCL14_Adip2 <- subset(significant, row == 'CXCL14')
str(CXCL14_Adip2)

##Pre_Ad
Pre_Ad <- av %>% select(81:100)
colnames(Pre_Ad) <- str_remove_all(colnames(Pre_Ad), "filtered.Pre_Ad_")
Pre_Ad <- t(as.data.frame(Pre_Ad))

cor <- rcorr(as.matrix(Pre_Ad), type = 'pearson')
mat <- flattenCorrMatrix(cor$r, cor$P)

mat$Pvalue_adjust <- p.adjust(mat$p, method="bonferroni")
significant <- subset(mat, p < 0.05)
str(significant)
CXCL14_Pre_Ad <- subset(significant, row == 'CXCL14')
str(CXCL14_Pre_Ad)

##Vascular
Vascular <- av %>% select(121:140)
colnames(Vascular) <- str_remove_all(colnames(Vascular), "filtered.Pre_Ad_")
Vascular <- t(as.data.frame(Vascular))

cor <- rcorr(as.matrix(Vascular), type = 'pearson')
mat <- flattenCorrMatrix(cor$r, cor$P)

mat$Pvalue_adjust <- p.adjust(mat$p, method="bonferroni")
significant <- subset(mat, p < 0.05)
str(significant)
CXCL14_Vascular <- subset(significant, row == 'CXCL14')
str(CXCL14_Vascular)

###Subset only posive correlations

CXCL14_Adip1_pos <- subset(CXCL14_Adip1, cor > 0)
CXCL14_Adip2_pos <- subset(CXCL14_Adip2, cor > 0)
CXCL14_Pre_Ad_pos <- subset(CXCL14_Pre_Ad, cor > 0)
CXCL14_Vascular_pos <- subset(CXCL14_Vascular, cor > 0)

colnames(CXCL14_Adip1_pos) = c('CXCL14', 'Gene', 'Adip1_cor', 'Adip_1_p', 'Adip_1_Padj')
colnames(CXCL14_Adip2_pos) = c('CXCL14', 'Gene', 'Adip2_cor', 'Adip_2_p', 'Adip_2_Padj')
colnames(CXCL14_Pre_Ad_pos) = c('CXCL14', 'Gene', 'Pre_Ad_cor', 'Pre_Ad_p', 'Pre_Ad_Padj')
colnames(CXCL14_Vascular_pos) = c('CXCL14', 'Gene', 'Vascular_cor', 'Vascular_p', 'Vascular_Padj')

Gene <- unique(c(CXCL14_Adip1_pos$Gene, CXCL14_Adip2_pos$Gene,  CXCL14_Pre_Ad_pos$Gene, CXCL14_Vascular_pos$Gene))

mega <- as.data.frame(Gene)

mega <- left_join(mega, CXCL14_Adip1_pos, by = 'Gene')
mega <- left_join(mega, CXCL14_Pre_Ad_pos, by = 'Gene')
mega <- left_join(mega, CXCL14_Adip2_pos, by = 'Gene')
mega <- left_join(mega, CXCL14_Vascular_pos, by = 'Gene')

head(mega)
mega2 <-  mega %>% select(1, 3, 6, 9, 12)
head(mega2)

colnames(mega2) <- str_remove_all(colnames(mega2), "_cor")


df.long <- pivot_longer(data = mega2, cols = -Gene, names_to = "CellType", values_to = "R_Value")
head(df.long)

ggplot(data = df.long, mapping = aes(x = Gene, y = CellType, fill = R_Value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option="plasma") +
  #coord_flip() +
  theme_classic() +
  theme(axis.title=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x =element_blank(),
        axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10, angle = 45, hjust = 1)) 


####Correlate genes expression with WHR for each cell type
av <- AverageExpression(seurat, assays = 'filtered', group.by = 'Cell_ID', slot = 'data')
av <- as.data.frame(av)

##Adip1 
Adip1 <- av %>% select(1:20)
colnames(Adip1) <- str_remove_all(colnames(Adip1), "filtered.Adip_1_")
Adip1 <- as.data.frame(t(Adip1))
Adip1[1:6, 1:6]
Adip1$WHR <- c(0.91, 0.9, 0.93, 0.82, 0.94, 0.88, 1, 0.77, 0.89, 0.82, 0.91, 0.98, 0.84, 0.77, 0.9, 1.05, 1.08, 1.04, 0.7, 1.1)

cor <- rcorr(as.matrix(Adip1), type = 'pearson')
tail(cor$P)


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

mat <- flattenCorrMatrix(cor$r, cor$P)

significant <- subset(mat, p < 0.05)
WHR_Adip1 <- subset(significant, column == 'WHR')
str(WHR_Adip1)

##Adip2 
Adip2 <- av %>% select(21:40)
colnames(Adip2) <- str_remove_all(colnames(Adip2), "filtered.Adip_2_")
Adip2 <- as.data.frame(t(Adip2))
Adip2[1:6, 1:6]
Adip2$WHR <- c(0.91, 0.9, 0.93, 0.82, 0.94, 0.88, 1, 0.77, 0.89, 0.82, 0.91, 0.98, 0.84, 0.77, 0.9, 1.05, 1.08, 1.04, 0.7, 1.1)


cor <- rcorr(as.matrix(Adip2), type = 'pearson')
mat <- flattenCorrMatrix(cor$r, cor$P)

significant <- subset(mat, p < 0.05)
WHR_Adip2 <- subset(significant, column == 'WHR')
str(WHR_Adip2)

##Pre_Ad
Pre_Ad <- av %>% select(81:100)
colnames(Pre_Ad) <- str_remove_all(colnames(Pre_Ad), "filtered.Pre_Ad_")
Pre_Ad <- as.data.frame(t(Pre_Ad))
Pre_Ad[1:6, 1:6]
Pre_Ad$WHR <- c(0.91, 0.9, 0.93, 0.82, 0.94, 0.88, 1, 0.77, 0.89, 0.82, 0.91, 0.98, 0.84, 0.77, 0.9, 1.05, 1.08, 1.04, 0.7, 1.1)


cor <- rcorr(as.matrix(Pre_Ad), type = 'pearson')
mat <- flattenCorrMatrix(cor$r, cor$P)

significant <- subset(mat, p < 0.05)
WHR_Pre_Ad <- subset(significant, column == 'WHR')
str(WHR_Pre_Ad)

##Vascular
Vascular <- av %>% select(121:140)
colnames(Vascular) <- str_remove_all(colnames(Vascular), "filtered.Vascular_")
Vascular <- as.data.frame(t(Vascular))
Vascular[1:6, 1:6]
Vascular$WHR <- c(0.91, 0.9, 0.93, 0.82, 0.94, 0.88, 1, 0.77, 0.89, 0.82, 0.91, 0.98, 0.84, 0.77, 0.9, 1.05, 1.08, 1.04, 0.7, 1.1)


cor <- rcorr(as.matrix(Vascular), type = 'pearson')
mat <- flattenCorrMatrix(cor$r, cor$P)

significant <- subset(mat, p < 0.05)
WHR_Vascular <- subset(significant, column == 'WHR')
str(WHR_Vascular)


###Subset only posive correlations

WHR_Adip1_pos <- subset(WHR_Adip1, cor > 0)
WHR_Adip2_pos <- subset(WHR_Adip2, cor > 0)
WHR_Pre_Ad_pos <- subset(WHR_Pre_Ad, cor > 0)
WHR_Vascular_pos <- subset(WHR_Vascular, cor > 0)


colnames(WHR_Adip1_pos) = c('Gene', 'WHR', 'Adip1_cor', 'Adip_1_p')
colnames(WHR_Adip2_pos) = c('Gene', 'WHR', 'Adip2_cor', 'Adip_2_p')
colnames(WHR_Pre_Ad_pos) = c('Gene', 'WHR', 'Pre_Ad_cor', 'Pre_Ad_p')
colnames(WHR_Vascular_pos) = c('Gene', 'WHR', 'Vascular_cor', 'Vascular_p')

Gene <- unique(c(WHR_Adip1_pos$Gene, WHR_Adip2_pos$Gene,  WHR_Pre_Ad_pos$Gene, WHR_Vascular_pos$Gene))

mega <- as.data.frame(Gene)

mega <- left_join(mega, WHR_Vascular_pos, by = 'Gene')
mega <- left_join(mega, WHR_Adip1_pos, by = 'Gene')
mega <- left_join(mega, WHR_Adip2_pos, by = 'Gene')
mega <- left_join(mega, WHR_Pre_Ad_pos, by = 'Gene')

head(mega)
mega2 <-  mega %>% select(1, 3, 6, 9, 12)
head(mega2)

write.csv(mega2, 'WHR_Gene_correlations.csv')

colnames(mega2) <- str_remove_all(colnames(mega2), "_cor")


df.long <- pivot_longer(data = mega2, cols = -Gene, names_to = "CellType", values_to = "R_Value")
head(df.long)



ggplot(data = df.long, mapping = aes(x = Gene, y = CellType, fill = R_Value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option="plasma") +
  #coord_flip() +
  theme_classic() +
  theme(axis.title=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x =element_blank(),
        axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10, angle = 45, hjust = 1)) 

####Compare WHR positive genes with genes upregulated in older cell populations of 
library(ggvenn)
Adip1 <- list(Older_Adip_1 = Adip_1_Up_genes,
              WHR_Adip1_pos = WHR_Adip1_pos$Gene)

P1 <- ggvenn(
  Adip1, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

Adip2 <- list(Older_Adip_2 = Adip_2_Up_genes,
              WHR_Adip2_pos = WHR_Adip2_pos$Gene)

P2 <- ggvenn(
  Adip2, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

Vascular <- list(Older_Vascular = EC_Up_genes,
                 WHR_Vascular_pos = WHR_Vascular_pos$Gene)

P3 <- ggvenn(
  Vascular, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

Pre_Ad <- list(Older_Pre_Ad = EC_Up_genes,
               WHR_Pre_Ad_pos = WHR_Pre_Ad_pos$Gene)

P4 <- ggvenn(
  Pre_Ad, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

P1 + P2 + P3 + P4



####Plot out fibrosis and ECM remodelling markers
##Markers
fib <- c('ITGB5', ##integrins
         'ITGA1',
         'ITGA7',
         'ITGA9',
         'ITGAE',
         'ITGAX',
         'ITGB2',
         'ITGB4',
         'ITGA10',
         'COL3A1', #Collagens
         'COL6A1', 
         'COL6A3',
         'COL1A2',
         'COL5A1',
         'COL5A2',
         'COL12A1',
         'LAMP1', #lyossome associated membrane protein 1
         'LAMA4', #Laminin subunit alpha 4
         'LAMB2', #Laminin subunit beta 2
         'LAMC3', #Laminin subunit gamma 3
         'LAMA2', #Laminin subunit alpha 2
         'PDGFRA', #known pro-fibrotic signal
         'LOXL2', #Lysyl Oxidase Like 2
         'LOXL3', #Lysyl Oxidase Like 3
         "F13A1", #Coagulation Factor XIII A Chain - metabolically unhealthy
         'FN1', #Fibronectin 1
         'SPP1', #ECM glycoprotein osteopontin (OPN), increased in obesity 
         'MMP2', ##Matrix metallopeptidase
         'CD44' ##OPN interagts with CD44 and integrins
         )


seurat@meta.data$Cell_Group <- paste0(seurat@meta.data$CellType,  sep = '_', seurat@meta.data$AgeGroup) ##use the column that your group data is stored under
av <- AverageExpression(seurat, assays = 'filtered', group.by = 'Cell_Group', slot = 'data', features = fib)
av2 <- AverageExpression(seurat, assays = 'filtered', group.by = 'Cell_Group', slot = 'scale.data', features = fib)


av <- as.data.frame(av2)

Stem <- av %>% select(11:12)
colnames(Stem) <- str_remove_all(colnames(Stem), "filtered.Stem_")
Stem <- Stem[order(-Stem$Younger),] 
Stem <- as.data.frame(t(Stem))
Stem <- Stem %>% rownames_to_column("group")
levels <- c('Younger', 'Older')
Stem$group <- factor(Stem$group, levels = levels)



P1 <- ggradar(
  Stem, 
  values.radar = c("-0.5", "0", "2.2"),
  grid.min = -0.5, grid.mid = 0, grid.max = 2.2,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c( "#87CEFA", "#EC7014"),
  fill = TRUE,
  fill.alpha = 0.5,
  # seuratckground and grid lines
  seuratckground.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  axis.label.size = 3,
  grid.label.size = 4,
  legend.text.size = 10
)
 P1 

 Pre_Ad <- av %>% select(9:10)
 colnames(Pre_Ad) <- str_remove_all(colnames(Pre_Ad), "filtered.Pre_Ad_")
 Pre_Ad <- Pre_Ad[order(-Pre_Ad$Older),] 
 Pre_Ad <- as.data.frame(t(Pre_Ad))
 min(Pre_Ad)
 max(Pre_Ad)
 Pre_Ad <- Pre_Ad %>% rownames_to_column("group")
 levels <- c('Younger', 'Older')
 Pre_Ad$group <- factor(Pre_Ad$group, levels = levels)

  P2 <- ggradar(
   Pre_Ad, 
   values.radar = c("-0.7", "0", "0.7"),
   grid.min = -0.7, grid.mid = 0, grid.max = 0.7,
   # Polygons
   group.line.width = 1, 
   group.point.size = 3,
   group.colours = c( "#87CEFA", "#EC7014"),
   fill = TRUE,
   fill.alpha = 0.5,
   # seuratckground and grid lines
   seuratckground.circle.colour = "white",
   gridline.mid.colour = "grey",
   legend.position = "bottom",
   axis.label.size = 3,
   grid.label.size = 4,
   legend.text.size = 10
 )
 P2


Adip1 <- av %>% select(1:2)
colnames(Adip1) <- str_remove_all(colnames(Adip1), "filtered.Adip_1_")
Adip1 <- Adip1[order(-Adip1$Older),] 
Adip1 <- as.data.frame(t(Adip1))
Adip1 <- Adip1 %>% rownames_to_column("group")
levels <- c('Younger', 'Older')
Adip1$group <- factor(Adip1$group, levels = levels)



P3 <- ggradar(
  Adip1, 
  values.radar = c("-0.4", "0", "0.2"),
  grid.min = -0.4, grid.mid = 0, grid.max = 0.2,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c( "#87CEFA", "#EC7014"),
  fill = TRUE,
  fill.alpha = 0.5,
  # seuratckground and grid lines
  seuratckground.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  axis.label.size = 3,
  grid.label.size = 4,
  legend.text.size = 10
)
P3

Adip_2 <- av %>% select(3:4)
colnames(Adip_2) <- str_remove_all(colnames(Adip_2), "filtered.Adip_2_")
Adip_2 <- Adip_2[order(-Adip_2$Older),] 
Adip_2 <- as.data.frame(t(Adip_2))
min(Adip_2)
max(Adip_2)
Adip_2 <- Adip_2 %>% rownames_to_column("group")
levels <- c('Younger', 'Older')
Adip_2$group <- factor(Adip_2$group, levels = levels)

P4 <- ggradar(
  Adip_2, 
  values.radar = c("-0.6", "0", "1.1"),
  grid.min = -0.6, grid.mid = 0, grid.max = 1.1,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c( "#87CEFA", "#EC7014"),
  fill = TRUE,
  fill.alpha = 0.5,
  # seuratckground and grid lines
  seuratckground.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  axis.label.size = 3,
  grid.label.size = 4,
  legend.text.size = 10
)
P4


Vascular <- av %>% select(13:14)
colnames(Vascular) <- str_remove_all(colnames(Vascular), "filtered.Vascular_")
Vascular <- Vascular[order(-Vascular$Older),] 
Vascular <- as.data.frame(t(Vascular))
min(Vascular)
max(Vascular)
Vascular <- Vascular %>% rownames_to_column("group")
levels <- c('Younger', 'Older')
Vascular$group <- factor(Vascular$group, levels = levels)

P5 <- ggradar(
  Vascular, 
  values.radar = c("-0.6", "0", "0.4"),
  grid.min = -0.6, grid.mid = 0, grid.max = 0.4,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c( "#87CEFA", "#EC7014"),
  fill = TRUE,
  fill.alpha = 0.5,
  # seuratckground and grid lines
  seuratckground.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  axis.label.size = 3,
  grid.label.size = 4,
  legend.text.size = 10,
)
P5

P1 + P2 + P3 + P4 + P5

ggarrange(
 P1, P2, P3, P4, P5, labels = c("Stem", "Pre_Ad", 'Adip_1', 'Adip2', 'Vascular'),
  common.legend = TRUE, legend = "bottom"
)


saveRDS(seurat, 'harmony.sct.V3.RDS')


