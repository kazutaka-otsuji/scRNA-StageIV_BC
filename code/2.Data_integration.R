#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 2.Data integration

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)

set.seed(1)

### Merge 4 samples (w/o integration)
ump_Pre_CCR@meta.data$batch <- "1st"
ump_Post1_CCR@meta.data$batch <- "2nd"
ump_Post2_CCR@meta.data$batch <- "1st"
ump_Meta_CCR@meta.data$batch <- "3rd"

obj_all <- merge(ump_Pre_CCR, 
                 y = c(ump_Post1_CCR, ump_Post2_CCR, ump_Meta_CCR), 
                 add.cell.ids = c("Pre", "Post1", "Post2", "Meta"))

### Normalize, Find variable features, scale data with cell cycle regression, Run PCA, and Run UMAP
obj_all <- NormalizeData(obj_all, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:50) %>%
  FindClusters(., resolution = 0.8) %>%
  RunUMAP(., dims=1:50)

obj_all$Clusters <- factor(paste0("C", as.numeric(obj_all@meta.data$seurat_clusters)))
obj_all@meta.data$Clusters <- factor(obj_all@meta.data$Clusters, 
                                         levels = c(paste0("C", 1:18)))

# Plots
Idents(obj_all)<-"Clusters"
cols <- c("#386CB0", "#F0027F", "#B2DF8A", "#FFFF99", "#7570B3", "#7FC97F", "#6A3D9A", "#B3CDE3", "#A6761D", "#D95F02",
          "#E6AB02", "#FB9A99", "#CAB2D6", "#1B9E77", "#FF7F00", "#666666", "#FDC086", "#B15928")
p <- DimPlot(obj_all, reduction = "umap",label = T, label.size = 6, pt.size = 0.2, repel=T, cols=cols)+ 
  NoAxes() +
  ggtitle("Clusters (unintegrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0, size=20),
        legend.key.size = unit(0, 'lines'))
ggsave(file="Figures/Fig.S3A_All_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)

Idents(obj_all)<-"orig.ident"
cols <- c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB")
p <- DimPlot(obj_all, reduction = "umap",label = F,  pt.size = 0.1, cols=cols) +
  NoAxes() +
  ggtitle("Sample identification (unintegrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S3A_All_Sample.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)

Idents(obj_all)<-"batch"
p <- DimPlot(obj_all, reduction = "umap", group.by = "batch", label = F, cols = "Set1")+ 
  NoAxes() +
  ggtitle("Batch (unintegrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S3A_All_Batch.png", plot=p,  dpi=300, bg="white",units="in", width=5, height=4)


Idents(obj_all)<-"seurat_clusters"
Celltype <- c("Cancer-cell", "Cancer-cell", "Cancer-cell", "T/NK-cell", "Myeloid", "T/NK-cell", 
              "Cancer-cell", "Cancer-cell", "T/NK-cell", "T/NK-cell", "Cancer-cell",
              "Cancer-cell", "FB/ET", "B-cell", "T/NK-cell", "T/NK-cell", "Cancer-cell", "B-cell")
current.cluster.ids <- c(0:17)
obj_all@meta.data[,'Celltype'] <- plyr::mapvalues(x = obj_all@meta.data$seurat_clusters,
                                                      from = current.cluster.ids,
                                                      to = Celltype)
obj_all@meta.data$Celltype <- as.factor(obj_all@meta.data$Celltype)

Idents(obj_all)<-"Celltype"
p <- DimPlot(obj_all, reduction = "umap",label = T, label.size = 7, pt.size = 0.1, repel=T, cols="Spectral")+ 
  NoAxes() +
  ggtitle("Cell types (unintegrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S3A_All_Celltype.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


### Data integration with STACAS
library(STACAS)
nfeatures <- 1000
ndim <- 15

object <- obj_all
object = JoinLayers(object)
obj.list <- SplitObject(object, split.by = "batch")
stacas_anchors <- FindAnchors.STACAS(obj.list, 
                                     anchor.features = nfeatures,
                                     dims = 1:ndim)

st1 <- SampleTree.STACAS(
  anchorset = stacas_anchors,
  obj.names = names(obj.list)
)

object_integrated <- IntegrateData.STACAS(stacas_anchors,
                                          sample.tree = st1,
                                          dims=1:ndim)

object_integrated <- object_integrated %>%
  ScaleData() %>%
  RunPCA(npcs = ndim)%>%
  RunUMAP(dims = 1:ndim) %>%
  FindNeighbors(dims = 1:ndim) %>%
  FindClusters(resolution = 0.7) 

# Plots  
p <- DimPlot(object_integrated, reduction = "umap", group.by = "batch", label = F, cols = "Set1")+ 
  NoAxes() +
  ggtitle("Batch (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S3B_STACAS_All_Batch.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)

object_integrated@meta.data$Celltype <- as.factor(object_integrated@meta.data$Celltype) %>%
  factor(object_integrated@meta.data$Celltype,
         levels = c("Cancer-cell","T/NK-cell","Myeloid","B-cell", "FB/ET"))

p <- DimPlot(object_integrated, reduction = "umap", group.by = "Celltype",label = T, label.size = 7, pt.size = 0.1, repel=T, cols="Spectral")+ 
  NoAxes() +
  ggtitle("Celltype (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.1B_STACAS_All_Celltype.png", plot=p,  dpi=300, bg="white", units="in", width=5.5, height=4)

object_integrated@meta.data$orig.ident <- as.factor(object_integrated@meta.data$orig.ident) %>%
  factor(object_integrated@meta.data$orig.ident,
         levels = c("Pre","Post1","Post2","Meta"))
cols <- c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB")

p <- DimPlot(object_integrated, reduction = "umap", group.by = "orig.ident", label = F, cols=cols)+ 
  NoAxes() +
  ggtitle("Sample identification (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.1C_STACAS_All_Sample.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)

object_integrated@meta.data$Clusters <- factor(1+as.numeric(object_integrated@meta.data$seurat_clusters)) %>%
  factor(paste0("C", as.numeric(object_integrated@meta.data$Clusters)), levels = c(paste0("C", 1:18)))

cols <- c("#386CB0",  "#B2DF8A","#F0027F", "#FFFF99",  "#7570B3", 
          "#7FC97F", "#6A3D9A", "#A6761D", "#B3CDE3", "#FB9A99", 
          "#E6AB02", "#D95F02", "#CAB2D6", "#1B9E77", "#FF7F00",
          "#666666")

p <- DimPlot(object_integrated, reduction = "umap", group.by = "Clusters", label = T, label.size = 6, pt.size = 0.1, repel=T, cols=cols)+ 
  NoAxes() +
  ggtitle("Clusters (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20),
        legend.key.size = unit(0, 'lines'))
ggsave(file="Figures/Fig.S3B_STACAS_All_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


DefaultAssay(object_integrated) <- "RNA"
p <- FeaturePlot(object_integrated,
                 features = c("ESR1", "FOXA1","ERBB2", "MKI67",
                              "EPCAM", "KRT8", "KRT14", "VIM"), ncol=4, order=T, min.cutoff = 0)
ggsave(file="Figures/Fig.S3D.png", plot=p, dpi=300, width=12,height=5)


## Extract cancer cells
## Pre-treatment vs Post-treatment
Idents(object_integrated)<-"Clusters"
obj_cancer <- subset(object_integrated, idents = c("C1", "C2", "C3", "C5", "C6", "C9", "C10", "C14"))

Idents(obj_cancer)<-"orig.ident"
sample.id2 <- c("Pretreatment", "Posttreatment", "Posttreatment", "Posttreatment")
current.cluster.ids <- c("Pre", "Post1", "Post2", "Meta")
obj_cancer@meta.data[,"sample.id2"] <- plyr::mapvalues(x = obj_cancer@meta.data$orig.ident,
                                                       from = current.cluster.ids, to = sample.id2)


# Find differentially expressed features between "Pre-treatment" and "Post-treatment"
Idents(obj_cancer)<-"sample.id2"

PrevsPost.de.markers <- FindMarkers(obj_cancer, 
                                    ident.1 = "Posttreatment", ident.2 = "Pretreatment", 
                                    verbose = FALSE) #Wilcoxon Rank Sum Test
PrevsPost.de.markers <- PrevsPost.de.markers %>% 
  tibble::rownames_to_column(var = "gene")
head(PrevsPost.de.markers)

# The basic scatter plot: x is "avg_log2FC", y is "p_val_adj"
# add a column of NAs
PrevsPost.de.markers$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
PrevsPost.de.markers$diffexpressed[PrevsPost.de.markers$avg_log2FC >1  & PrevsPost.de.markers$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
PrevsPost.de.markers$diffexpressed[PrevsPost.de.markers$avg_log2FC < -1 & PrevsPost.de.markers$p_val_adj < 0.05] <- "DOWN"

PrevsPost.de.markers$delabel <- NA
PrevsPost.de.markers$delabel[PrevsPost.de.markers$diffexpressed != "NO"] <- PrevsPost.de.markers$gene[PrevsPost.de.markers$diffexpressed != "NO"]

library(ggrepel)
genes.to.label = c("CALCA", "BMP5", "MARCOL", "NRXN2", "FABP4", "PPDPFL", "H19", 
                   "ACAN",  "SNAI2",  "PDE3A",
                   "KLK6", "KLK7", "LTF", "ELF5", "KLK5", "PROM1", "TMEM213", "IRX1", "KLK8", "SLC38A3", 
                   "OBP2B", "FGF13", "CLDN3", "CST6", "TACSTD2")


min_positive_pval <- min(PrevsPost.de.markers$p_val_adj[PrevsPost.de.markers$p_val_adj > 0], na.rm = TRUE)
replacement_value <- min_positive_pval / 2
PrevsPost.de.markers$p_val_adj[PrevsPost.de.markers$p_val_adj == 0] <- replacement_value


PrevsPost.de.markers$label <- ifelse(PrevsPost.de.markers$gene %in% genes.to.label, PrevsPost.de.markers$gene, NA)

# Plot
p <- ggplot(data=PrevsPost.de.markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 10, na.rm = TRUE) +
  scale_color_manual(values=c("#004488", "#BBBBBB", "darkred"))

ggsave(file="Figures/Fig.1D.png", plot=pC4, bg="white", dpi=300, width=8,height=5)




### Merge Post1 and Meta ###
obj <- merge(ump_Post1_CCR, 
             y = c(ump_Meta_CCR), 
             add.cell.ids = c("Post1", "Meta"))

### Normalize, Find variable features, scale data with cell cycle regression, Run PCA, and Run UMAP
obj_Post1_Meta <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:40) %>%
  FindClusters(., resolution = 0.5) %>%
  RunUMAP(., dims=1:40)


# Annotation
Idents(obj_Post1_Meta)<-"seurat_clusters"
Celltype <- c("T/NK-cell", "Cancer-cell", "Cancer-cell", "Cancer-cell", "B-cell/Myeloid", "Cancer-cell", 
              "Cancer-cell", "T/NK-cell", "T/NK-cell", "B-cell/Myeloid", "FB/ET",
              "T/NK-cell")
current.cluster.ids <- c(0:11)
obj_Post1_Meta@meta.data[,'Celltype'] <- plyr::mapvalues(x = obj_Post1_Meta@meta.data$seurat_clusters,
                                                              from = current.cluster.ids,
                                                              to = Celltype)

obj_Post1_Meta@meta.data$Celltype <- as.factor(obj_Post1_Meta@meta.data$Celltype)



#### STACAS ####
library(STACAS)

nfeatures <- 1000
ndim <- 15

object <- obj_Post1_Meta
object = JoinLayers(object)

obj.list <- SplitObject(object, split.by = "batch")

stacas_anchors <- FindAnchors.STACAS(obj.list, 
                                     anchor.features = nfeatures,
                                     dims = 1:ndim)

st1 <- SampleTree.STACAS(
  anchorset = stacas_anchors,
  obj.names = names(obj.list)
)

object_integrated <- IntegrateData.STACAS(stacas_anchors,
                                          sample.tree = st1,
                                          dims=1:ndim)

object_integrated <- object_integrated %>%
  ScaleData() %>%
  RunPCA(npcs = ndim) %>%
  FindNeighbors(dims = 1:ndim) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:ndim)


object_integrated@meta.data$Celltype <- factor(object_integrated@meta.data$Celltype,
                                               levels = c("Cancer-cell","T/NK-cell","B-cell/Myeloid","FB/ET"))


p <- DimPlot(object_integrated, reduction = "umap", group.by = "Celltype",
             label = T, label.size = 7, pt.size = 0.1,
             repel=T, cols="Spectral")+ 
  NoAxes() +
  ggtitle("Celltype (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S9A_All_Celltype.png", plot=p,  dpi=300, bg="white", 
       units="in", width=5.5, height=4)


object_integrated@meta.data$orig.ident <- factor(object_integrated@meta.data$orig.ident,
                                                 levels = c("Post1","Meta"))

cols <- c("#66C2A5", "#8DA0CB")
p <- DimPlot(object_integrated, reduction = "umap", group.by = "orig.ident",
             label = F, cols=cols)+ 
  NoAxes() +
  ggtitle("Sample identification (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S9A_All_Sample.png", plot=p,  dpi=300, bg="white", 
       units="in", width=5, height=4)


object_integrated@meta.data$Clusters <- factor(1+as.numeric(object_integrated@meta.data$seurat_clusters))
object_integrated@meta.data$Clusters <- factor(paste0("C", as.numeric(object_integrated@meta.data$Clusters)), 
                                               levels = c(paste0("C", 1:14)))


cols <- c("#386CB0", "#F0027F", "#B2DF8A","#7570B3", "#FFFF99",  
          "#7FC97F", "#A6761D", "#6A3D9A", "#B3CDE3", 
          "#FB9A99", "#1B9E77", "#FF7F00",
          "#666666", "#FDC086")　

Idents(object_integrated)<-"Clusters"
p <- DimPlot(object_integrated, reduction = "umap", group.by = "Clusters",
             label = T, label.size = 6, pt.size = 0.1,
             repel=T, cols=cols)+ 
  NoAxes() +
  ggtitle("Clusters (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20),
        legend.key.size = unit(0, 'lines'))
ggsave(file="Figures/Fig.S9A_All_Cluster.png", plot=p,  dpi=300, bg="white", 
       units="in", width=5, height=4)




## Extract cancer cells
Idents(object_integrated)<-"Clusters"
obj_cancer <- subset(object_integrated, idents = c("C2", "C3", "C4", "C6", "C9", "C10", "C12"))
obj_cancer@meta.data$Clusters <- droplevels(obj_cancer@meta.data$Clusters)

Idents(obj_cancer)<-"Clusters"
cols <- c("#F0027F","#B2DF8A",  "#7570B3", "#B3CDE3",  "#FB9A99")

p <- DimPlot(obj_cancer, reduction = "umap", group.by = "Clusters",
             label = T, label.size = 6, pt.size = 0.1,
             repel=T, cols=cols)+ 
  NoAxes() +
  ggtitle("Clusters of cancer (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20),
        legend.key.size = unit(0, 'lines'))
ggsave(file="Figures/Fig.4A_Cancer_Cluster.png", plot=p,  dpi=300, bg="white", 
       units="in", width=5, height=4)


cols <- c("#66C2A5", "#8DA0CB")
p <- DimPlot(obj_cancer, reduction = "umap", group.by = "orig.ident",
             label = F, pt.size = 0.1,cols=cols)+ 
  NoAxes() +
  ggtitle("Samples of cancer (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.4B_Cancer_Sample.png", plot=p,  dpi=300, bg="white", 
       units="in", width=5, height=4)


obj_cancer@meta.data$Phase <- factor(obj_cancer@meta.data$Phase,
                                     levels = c("G1", "S", "G2M"))

p <- DimPlot(obj_cancer, reduction = "umap", group.by = "Phase",
             label = F, pt.size = 0.1,cols="Dark2")+ 
  NoAxes() +
  ggtitle("Cell cycle phase of cancer (integrated)")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=18))
ggsave(file="Figures/Fig.S9B_STACAS_Cancer_Phase.png", plot=p,  dpi=300, bg="white", 
       units="in", width=5, height=4)
