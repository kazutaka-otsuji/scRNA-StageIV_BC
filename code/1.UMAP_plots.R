#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 1. UMAP plots

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)

set.seed(1)

### Read 10x data for each sample
data_Pre <- Read10X(data.dir = "path/Pre")
data_Pre <- CreateSeuratObject(counts = data_Pre, min.cells = 3, min.features = 200, project = "Pre")

data_Post1 <- Read10X(data.dir = "path/Post1")
data_Post1 <- CreateSeuratObject(counts = data_Post1, min.cells = 3, min.features = 200, project = "Post1")

data_Post2 <- Read10X(data.dir = "path/Post2")
data_Post2 <- CreateSeuratObject(counts = data_Post2, min.cells = 3, min.features = 200, project = "Post2")

data_Meta <- Read10X(data.dir = "path/Meta")
data_Meta <- CreateSeuratObject(counts = data_Meta, min.cells = 3, min.features = 200, project = "Meta")

### Normalize, Find variable features, scale data, Run PCA, and Run UMAP
ump_Pre <- NormalizeData(data_Pre, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
  ScaleData(., features = rownames(.)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:35) %>%
  FindClusters(., resolution = 0.5) %>%
  RunUMAP(., dims=1:35)

ump_Post1 <- NormalizeData(data_Post1, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
  ScaleData(., features = rownames(.)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:25) %>%
  FindClusters(., resolution = 0.5) %>%
  RunUMAP(., dims=1:25)

ump_Post2 <- NormalizeData(data_Post2, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
  ScaleData(., features = rownames(.)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:35) %>%
  FindClusters(., resolution = 0.5) %>%
  RunUMAP(., dims=1:35)

ump_Meta <- NormalizeData(data_Meta, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
  ScaleData(., features = rownames(.)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:10) %>%
  FindClusters(., resolution = 1) %>%
  RunUMAP(., dims=1:10)

### Cell cycle scoring
exp.mat <- read.table(file = "path/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

ump_Pre <- CellCycleScoring(ump_Pre, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ump_Post1 <- CellCycleScoring(ump_Post1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ump_Post2 <- CellCycleScoring(ump_Post2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ump_Meta<- CellCycleScoring(ump_Meta, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

### cell cycle regress, Run PCA, and Run UMAP
ump_Pre_CCR <- ScaleData(ump_Pre, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ump_Pre)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:25) %>%
  FindClusters(., resolution = 0.8) %>%
  RunUMAP(., dims=1:25)

ump_Post1_CCR <- ScaleData(ump_Post1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ump_Post1)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:40) %>%
  FindClusters(., resolution = 0.5) %>%
  RunUMAP(., dims=1:40)

ump_Post2_CCR <- ScaleData(ump_Post2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ump_Post2)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:35) %>%
  FindClusters(., resolution = 0.5) %>%
  RunUMAP(., dims=1:35)

ump_Meta_CCR <- ScaleData(ump_Meta, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ump_Meta)) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:11) %>%
  FindClusters(., resolution = 0.6) %>%
  RunUMAP(., dims=1:11)

### UMAP plot
ump_Pre_CCR$Clusters <- factor(paste0("C", as.numeric(ump_Pre_CCR@meta.data$seurat_clusters)))
ump_Pre_CCR@meta.data$Clusters <- factor(ump_Pre_CCR@meta.data$Clusters, 
                                             levels = c(paste0("C", 1:15)))

p <- DimPlot(ump_Pre_CCR, reduction = "umap",label = T, label.size = 6, pt.size = 0.2, repel=T)+ 
  NoAxes() +
  ggtitle("Clusters")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Pre_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


ump_Post1_CCR$Clusters <- factor(paste0("C", as.numeric(ump_Post1_CCR@meta.data$seurat_clusters)))
ump_Post1_CCR@meta.data$Clusters <- factor(ump_Post1_CCR@meta.data$Clusters, 
                                         levels = c(paste0("C", 1:10)))

p <- DimPlot(ump_Post1_CCR, reduction = "umap",label = T, label.size = 6, pt.size = 0.2, repel=T)+ 
  NoAxes() +
  ggtitle("Clusters")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Post1_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


ump_Post2_CCR$Clusters <- factor(paste0("C", as.numeric(ump_Post2_CCR@meta.data$seurat_clusters)))
ump_Post2_CCR@meta.data$Clusters <- factor(ump_Post2_CCR@meta.data$Clusters, 
                                           levels = c(paste0("C", 1:9)))

p <- DimPlot(ump_Post2_CCR, reduction = "umap",label = T, label.size = 6, pt.size = 0.2, repel=T)+ 
  NoAxes() +
  ggtitle("Clusters")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Post2_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


ump_Meta_CCR$Clusters <- factor(paste0("C", as.numeric(ump_Meta_CCR@meta.data$seurat_clusters)))
ump_Meta_CCR@meta.data$Clusters <- factor(ump_Meta_CCR@meta.data$Clusters, 
                                           levels = c(paste0("C", 1:4)))

p <- DimPlot(ump_Meta_CCR, reduction = "umap",label = T, label.size = 6, pt.size = 0.2, repel=T)+ 
  NoAxes() +
  ggtitle("Clusters")+
  theme(text=element_text(family="Arial", size=12),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Meta_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


### Cell-type annotation
obj <- ump_Pre_CCR # example

# 1. FindAllMarkers
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1) %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(obj, features = top10$gene)


# 2. Marker gene expression
DotPlot(obj,
        features = c("THBD", "PECAM1", # endothelial
                     "MYL9", "ACTA2", "FAP", # fibroblast
                     "FCER1G", "CD1C", "CLEC10A", # myeloid
                     "KLRD1", "KLRB1", # NK cell
                     "CD3D", # T cell
                     "CD79A", # B cell
                     "KRT14", "KRT8", "EPCAM" # breast epithelial)
        ))+coord_flip()&
          theme(text=element_text(family="Arial"))+ 
          RotatedAxis()
        
# 3. SingleR
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
        
library(SingleR)
data_for_singleR <- GetAssayData(obj, assay="RNA", layer="counts")
pred <- SingleR(test = data_for_singleR, ref = hpca.se, labels = hpca.se$label.main)
obj$cell_type <- pred$labels
table(obj@meta.data$cell_type, obj@meta.data$Clusters)
        
# 4. Add "Celltype" into meta.data and plot
Idents(ump_Pre_CCR)<-"seurat_clusters"
Celltype <- c("Cancer-cell", "T/NK-cell", "T/NK-cell", "Cancer-cell", "T/NK-cell", "Cancer-cell",
              "Myeloid", "FB/ET", "Cancer-cell", "T/NK-cell", "B-cell",
              "Myeloid", "B-cell", "Cancer-cell", "CMP")
current.cluster.ids <- c(0:14)
ump_Pre_CCR@meta.data[,'Celltype'] <- plyr::mapvalues(x = ump_Pre_CCR@meta.data$seurat_clusters,
                                                      from = current.cluster.ids,
                                                      to = Celltype)
ump_Pre_CCR@meta.data$Celltype <- as.factor(ump_Pre_CCR@meta.data$Celltype)
        
Idents(ump_Pre_CCR)<-"Celltype"
p <- DimPlot(ump_Pre_CCR, reduction = "umap",label = T, label.size = 7, pt.size = 0.2, repel=T, cols="Spectral")+
     NoAxes() +
     ggtitle("Cell type identification")+
     theme(text=element_text(family="Arial", size=12),
     legend.text = element_text(size = 15),
     plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Pre_Celltype.png", plot=p,  dpi=300, bg="white",units="in", width=5, height=4)
        
        
Idents(ump_Post1_CCR)<-"seurat_clusters"
Celltype <- c("T/B/NK-cell", "Cancer-cell", "Cancer-cell", "Cancer-cell","Myeloid", "Cancer-cell",
              "T/B/NK-cell", "T/B/NK-cell", "Myeloid","FB/ET")
current.cluster.ids <- c(0:9)
ump_Post1_CCR@meta.data[,'Celltype'] <- plyr::mapvalues(x = ump_Post1_CCR@meta.data$seurat_clusters,
                                                        from = current.cluster.ids,
                                                        to = Celltype)
        
ump_Post1_CCR@meta.data$Celltype <- as.factor(ump_Post1_CCR@meta.data$Celltype)
ump_Post1_CCR@meta.data$Celltype <- factor(ump_Post1_CCR@meta.data$Celltype,
                                           levels = c("Cancer-cell","T/B/NK-cell","Myeloid","FB/ET"))
        
Idents(ump_Post1_CCR)<-"Celltype"
p <- DimPlot(ump_Post1_CCR, reduction = "umap",label = T, label.size = 7, pt.size = 0.2, repel=T, cols="Spectral")+ 
     NoAxes() +
     ggtitle("Cell type identification")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Post1_Celltype.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
Idents(ump_Post2_CCR)<-"seurat_clusters"
Celltype <- c("Cancer-cell","Cancer-cell", "Cancer-cell", "Cancer-cell","Cancer-cell","Myeloid",
              "Cancer-cell","FB/ET", "T/NK-cell")
current.cluster.ids <- c(0:8)
ump_Post2_CCR@meta.data[,'Celltype'] <- plyr::mapvalues(x = ump_Post2_CCR@meta.data$seurat_clusters,
                                                        from = current.cluster.ids,
                                                        to = Celltype)
        
Idents(ump_Post2_CCR)<-"Celltype"
p <- DimPlot(ump_Post2_CCR, reduction = "umap",label = T, label.size = 7, pt.size = 0.2, repel=T, cols="Spectral")+ 
     NoAxes() +
     ggtitle("Cell type identification")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Post2_Celltype.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
Idents(ump_Meta_CCR)<-"seurat_clusters"
Celltype <- c("Cancer-cell","Cancer-cell","T/NK-cell", "B-cell/Myeloid")
current.cluster.ids <- c(0:3)
ump_Meta_CCR@meta.data[,'Celltype'] <- plyr::mapvalues(x = ump_Meta_CCR@meta.data$seurat_clusters,
                                                           from = current.cluster.ids,
                                                           to = Celltype)
        
Idents(ump_Meta_CCR)<-"Celltype"
cols =c("#d7191c", "#fdae61", "#abdda4")
p <- DimPlot(ump_Meta_CCR, reduction = "umap",label = T, label.size = 7, pt.size = 1, repel=T, cols=cols)+ 
     NoAxes() +
     ggtitle("Cell type identification")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2A_Meta_Celltype.png", plot=p,  dpi=300, bg="white",  units="in", width=6, height=4)


# Save files
save(ump_Pre_CCR, file = "Files/ump_Pre_CCR.rda")
save(ump_Post1_CCR, file = "Files/ump_Post1_CCR.rda")
save(ump_Post2_CCR, file = "Files/ump_Post2_CCR.rda")
save(ump_Meta_CCR, file = "Files/ump_Meta_CCR.rda")


        
### Merge 4 samples (w/o integration)
ump_Pre_CCR@meta.data$batch <- "1st"
ump_Post1_CCR@meta.data$batch <- "2nd"
ump_Post2_CCR@meta.data$batch <- "1st"
ump_Meta_CCR@meta.data$batch <- "3rd"
        
obj_all <- merge(ump_Pre_CCR, 
                 y = c(ump_Post1_CCR, ump_Post2_CCR, ump_Meta_CCR), 
                 add.cell.ids = c("Pre", "Post1", "Post2", "Meta"))
        
### Normalize, Find variable features, scale data with cell cycle regression, Run PCA, and Run UMAP
ump_all_CCR <- NormalizeData(obj_all, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
               FindVariableFeatures(., selection.method = "vst",nfeatures = 1000) %>% 
               ScaleData(., vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
               RunPCA(.) %>%
               FindNeighbors(., dims = 1:50) %>%
               FindClusters(., resolution = 0.8) %>%
               RunUMAP(., dims=1:50)

ump_all_CCR$Clusters <- factor(paste0("C", as.numeric(ump_all_CCR@meta.data$seurat_clusters)))
ump_all_CCR@meta.data$Clusters <- factor(ump_all_CCR@meta.data$Clusters, 
                                          levels = c(paste0("C", 1:18)))

# Plots
Idents(ump_all_CCR)<-"Clusters"
cols <- c("#386CB0", "#F0027F", "#B2DF8A", "#FFFF99", "#7570B3", "#7FC97F", "#6A3D9A", "#B3CDE3", "#A6761D", "#D95F02",
          "#E6AB02", "#FB9A99", "#CAB2D6", "#1B9E77", "#FF7F00", "#666666", "#FDC086", "#B15928")
p <- DimPlot(ump_all_CCR, reduction = "umap",label = T, label.size = 6, pt.size = 0.2, repel=T, cols=cols)+ 
     NoAxes() +
     ggtitle("Clusters (unintegrated)")+
      theme(text=element_text(family="Arial", size=12),
            legend.text = element_text(size = 14),
            plot.title = element_text(hjust = 0, size=20),
            legend.key.size = unit(0, 'lines'))
ggsave(file="Figures/Fig.S2B_All_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
Idents(ump_all_CCR)<-"orig.ident"
cols <- c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB")
p <- DimPlot(ump_all_CCR, reduction = "umap",label = F,  pt.size = 0.1, cols=cols) +
     NoAxes() +
     ggtitle("Sample identification (unintegrated)")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2B_All_Sample.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
Idents(ump_all_CCR)<-"batch"
p <- DimPlot(ump_all_CCR, reduction = "umap", group.by = "batch", label = F, cols = "Set1")+ 
     NoAxes() +
     ggtitle("Batch (unintegrated)")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2_All_Batch.png", plot=p,  dpi=300, bg="white",units="in", width=5, height=4)
        
ump_all_CCR
Idents(ump_all_CCR)<-"seurat_clusters"
Celltype <- c("Cancer-cell", "Cancer-cell", "Cancer-cell", "T/NK-cell", "Myeloid", "T/NK-cell", 
              "Cancer-cell", "Cancer-cell", "T/NK-cell", "T/NK-cell", "Cancer-cell",
              "Cancer-cell", "FB/ET", "B-cell", "T/NK-cell", "T/NK-cell", "Cancer-cell", "B-cell")
 current.cluster.ids <- c(0:17)
 ump_all_CCR@meta.data[,'Celltype'] <- plyr::mapvalues(x = ump_all_CCR@meta.data$seurat_clusters,
                                                       from = current.cluster.ids,
                                                       to = Celltype)
ump_all_CCR@meta.data$Celltype <- as.factor(ump_all_CCR@meta.data$Celltype)
        
Idents(ump_all_CCR)<-"Celltype"
p <- DimPlot(ump_all_CCR, reduction = "umap",label = T, label.size = 7, pt.size = 0.1, repel=T, cols="Spectral")+ 
     NoAxes() +
     ggtitle("Cell types (unintegrated)")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20))
ggsave(file="Figures/Fig.S2B_All_Celltype.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
        
### Data integration with STACAS
library(STACAS)
nfeatures <- 1000
ndim <- 15
        
object <- ump_all_CCR
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
        
object_integrated <- object_integrated %>% ScaleData() %>%
                     RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim)

# Plots  
p <- DimPlot(object_integrated, reduction = "umap", group.by = "batch", label = F, cols = "Set1")+ 
     NoAxes() +
     ggtitle("Batch (integrated)")+
     theme(text=element_text(family="Arial", size=12),
                legend.text = element_text(size = 15),
                plot.title = element_text(hjust = 0, size=20))
     ggsave(file="Figures/Fig.S2C_STACAS_All_Batch.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
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
ggsave(file="Figures/Fig.1B_STACAS_All_Sample.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)
        
object_integrated@meta.data$Clusters <- factor(1+as.numeric(object_integrated@meta.data$seurat_clusters)) %>%
                                        factor(paste0("C", as.numeric(object_integrated@meta.data$Clusters)), levels = c(paste0("C", 1:18)))
        
cols <- c("#386CB0", "#F0027F", "#B2DF8A", "#FFFF99", "#7570B3", "#7FC97F", "#6A3D9A", "#B3CDE3", "#A6761D", "#D95F02", 
          "#E6AB02", "#FB9A99", "#CAB2D6", "#1B9E77", "#FF7F00", "#666666", "#FDC086", "#B15928")
p <- DimPlot(object_integrated, reduction = "umap", group.by = "Clusters", label = T, label.size = 6, pt.size = 0.1, repel=T, cols=cols)+ 
     NoAxes() +
     ggtitle("Clusters (integrated)")+
     theme(text=element_text(family="Arial", size=12),
           legend.text = element_text(size = 15),
           plot.title = element_text(hjust = 0, size=20),
           legend.key.size = unit(0, 'lines'))
ggsave(file="Figures/Fig.S2C_STACAS_All_Cluster.png", plot=p,  dpi=300, bg="white", units="in", width=5, height=4)


# Save files
save(object_integrated, file = "Files/object_integrated.rda")
