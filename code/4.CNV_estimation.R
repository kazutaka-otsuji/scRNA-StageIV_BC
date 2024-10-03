#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 4. CNV estimation

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

set.seed(1)

##### InferCNV #####
library("infercnv")
library("tibble")
library("tsvio")
library("NGCHM")
library("infercnvNGCHM")
library("rjags")
library("ape")


obj <- ump_Pre_CCR # example

### CreateInfercnvObject
mtx <- as.data.frame(obj@assays$RNA@counts)
cells <- Cells(obj)
V1 <- as.character(obj@meta.data$Clusters)

annot <- data.frame(row.names = cells, V1 = V1)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mtx,
                                        annotations_file=annot,
                                        delim="\t",
                                        gene_order_file=gene_order_file,
                                        ref_group_names=c("C7", "C8", "C12"),
                                        min_max_counts_per_cell=c(1e3,1e8),
                                        chr_exclude=c("chrY", "chrM"))


### perform infercnv operations to reveal cnv signal
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              min_cells_per_gene = 10,
                              out_dir="output_dir", 
                              analysis_mode = "subcluster",
                              cluster_by_groups=F, 
                              tumor_subcluster_partition_method='random_trees',
                              tumor_subcluster_pval=0.05,
                              denoise=TRUE,
                              HMM=F,
                              no_plot=F) 