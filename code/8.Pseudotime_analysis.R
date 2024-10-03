#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 8.Pseudotime_analysis
#### Monocle3 ####
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(magrittr)

set.seed(1234)

Idents(obj_cancer) <- "Phase" # obj_cancer: integrated cancer clusters of Post1 and Meta (Fig. 4A). 
obj_G1 <- subset(obj_cancer, idents="G1") # Only cancer cells in G1 cell-cycle phase are contained in this analysis.

cds <- obj_G1 %>%
  as.cell_data_set() %>%
  cluster_cells(resolution = 1e-3) %>%
  learn_graph(use_partition = TRUE, verbose = FALSE) %>%
  order_cells()

cds <- estimate_size_factors() %>%
  preprocess_cds()


rowData(cds)$gene_short_name <- rownames(cds)


cols <- c("Post1-CNV1"="#C2ADC0", "Post1-CNV2"="#FFED6F","Post1-CNV3"= "#C7D98C", "Post1-CNV4"="#DED7DA",
          "Meta" = "#C77CFF")

p5C <- plot_cells(cds,
                   color_cells_by="Post1_CNV1_Meta",
                   label_groups_by_cluster=F,
                   label_leaves=F,
                   label_cell_groups = F,
                   graph_label_size = 4,
                   cell_size = 0.5)+NoAxes()+
  theme(text=element_text(family="Arial"),
        legend.text=element_text(size=12),
        legend.key.size = unit(0, 'lines'),
        legend.position=c(0.2,0.3),
        legend.title=element_blank()) +
  scale_color_manual(values = cols)+ 
  guides(color = guide_legend(override.aes = list(size = 2))) 

ggsave(file="Figures/Fig.5C.png", plot=p5C, dpi=300, width=4,height=3.2)



p5D <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 group_cells_by = "cluster",
                 label_cell_groups = FALSE,
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 label_roots = FALSE,
                 trajectory_graph_segment_size = 1,
                 cell_size = 0.5)+NoAxes()+
  theme(text=element_text(family="Arial"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size = unit(0.2, 'inch'),
        legend.position=c(0.2,0.3)) 

ggsave(file="Figures/Fig.5D.png", plot=p5D, dpi=300, width=4,height=3.2)



p5E <- plot_cells(cds,
                 genes=c("CALCA", "CALCB","ADM", "CALCRL",  
                         "RAMP1", "RAMP2"),
                 show_trajectory_graph = FALSE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size = 0.4)+ NoAxes()+
  theme(text=element_text(family="Arial"),
        legend.key.size = unit(0.2, 'inch'),
        legend.position="bottom")

ggsave(file="Figures/Fig.5E.png", plot=p5E, dpi=300, width=4,height=3.2)


CALCA_genes <- c("CALCA", "CALCB","ADM", "CALCRL",  
                 "RAMP1", "RAMP2")
CALCA_lineage_cds <- cds[rowData(cds)$gene_short_name %in% CALCA_genes,
                         colData(cds)$Phase %in% c("G1")]

p5F <- plot_genes_in_pseudotime(CALCA_lineage_cds,
                                color_cells_by="Post1_CNV1_Meta",
                                min_expr=0.1,
                                cell_size = 1)+
  scale_color_manual(values = cols)+
  theme(text=element_text(family="Arial", size = 16),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.1, 'inch'),
        legend.position="bottom",
        legend.title=element_blank())+
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file="Figures/Fig.5F.png", plot=p13, dpi=300, width=4,height=8)