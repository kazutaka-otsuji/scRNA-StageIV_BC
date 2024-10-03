#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 6. TF activity estimation

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(BITFAM)

set.seed(1)

# obj_cancer: integrated cancer clusters of Post1 and Meta (Fig. 4A).
integrated_norm_data <- as.matrix(obj_cancer[["integrated"]]@data)
BITFAM_res_Post <- BITFAM(data = integrated_norm_data, species = "human", scATAC_obj = NA, ncores = parallel::detectCores())
Z_Post <- BITFAM_activities(BITFAM_res_Post)

Idents(obj_cancer) <- "Clusters"
Annotation <- as.data.frame(obj_cancer@active.ident)
Annotation$Clusters <- Annotation$Annotation


cluster_ids <- c("C2", "C3", "C4", "C6", "C9", "C10", "C12")
unique_clusters_df <- as.data.frame(sapply(cluster_ids, function(cluster) ifelse(Annotation$Clusters == cluster, 1, 0)))
colnames(unique_clusters_df) <- paste0("Integrated_Cluster", sub("C", "", cluster_ids))

# random forest model
rf_results <- lapply(1:ncol(unique_clusters_df), function(i) {
  cluster_name <- colnames(unique_clusters_df)[i]
  Z_Post_temp <- cbind(Z_Post, unique_clusters_df[, i])
  colnames(Z_Post_temp)[28] <- cluster_name
  Z_Post_temp[[cluster_name]] <- factor(Z_Post_temp[[cluster_name]])
  fit_rf <- randomForest(as.formula(paste(cluster_name, "~ .")), data = Z_Post_temp)
  importance(fit_rf)[order(importance(fit_rf)[, 1], decreasing = TRUE), ][1:10]
})


selected_tf_names <- unique(unlist(sapply(rf_results, rownames)))
selected_TF_Z <- Z_Post[, selected_tf_names, drop = FALSE]
selected_TF_Z <- selected_TF_Z[, !grepl("\\.", colnames(selected_TF_Z))]
selected_TF_Z <- apply(selected_TF_Z, 2, function(x) (x - min(x)) / (max(x) - min(x)))
selected_TF_Z2 <- t(selected_TF_Z)
colnames(selected_TF_Z2) <- Annotation$Clusters


# Heatmap
library(ComplexHeatmap)
library(tidyr)
library(circlize)
library(tibble)

colnames(selected_TF_Z2) <- make.names(colnames(selected_TF_Z2), unique = TRUE)
selected_TF_Z2_df <- selected_TF_Z2 %>%
  as.data.frame() %>%
  rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "Integrated Cluster", values_to = "Activity")

top_TFs <- selected_TF_Z2_df %>%
  group_by(`Integrated Cluster`, TF) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE), .groups = 'drop') %>%
  group_by(`Integrated Cluster`) %>%
  top_n(10, mean_activity) %>%
  ungroup()

heatmap_data <- selected_TF_Z2_df %>%
  filter(TF %in% top_TFs$TF) %>%
  pivot_wider(names_from = `Integrated Cluster`, values_from = Activity) %>%
  column_to_rownames("TF")

cluster_info <- gsub("\\..*", "", colnames(heatmap_data))
sorted_columns <- order(cluster_info)
heatmap_data <- heatmap_data[, sorted_columns]
cluster_info <- cluster_info[sorted_columns]

desired_cluster_order <- c("C3", "C10", "C12", "C6", "C2", "C9", "C4")
desired_meta_order <- c("Post1-CNV1", "Post1-CNV2", "Post1-CNV3", "Post1-CNV4", "Meta")

combined_df$Integrated_Clusters <- factor(combined_df$Integrated_Clusters, levels = desired_cluster_order)
combined_df$Post1_CNV1_Meta <- factor(combined_df$Post1_CNV1_Meta, levels = desired_meta_order)

cluster_annotation <- HeatmapAnnotation(
  `Integrated Clusters` = combined_df$Integrated_Clusters,
  `Post1-CNV&Meta` = combined_df$Post1_CNV1_Meta,
  col = list(
    `Integrated Clusters` = c("C2" = "#F0027F", "C3" = "#B2DF8A", "C4" = "#7570B3", "C6" = "#7FC97F",
                              "C9" = "#B3CDE3", "C10" = "#FB9A99", "C12" = "#FF7F00"),
    `Post1-CNV&Meta` = c("Post1-CNV1" = "#C2ADC0", "Post1-CNV2" = "#FFED6F", "Post1-CNV3" = "#C7D98C", 
                         "Post1-CNV4" = "#DED7DA", "Meta" = "#C77CFF")
  ),
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10)
)

col_fun <- colorRamp2(c(0, 1), c("white", "red"))


ht <- Heatmap(
  selected_TF_Z2, 
  name = "Activity",
  col = col_fun,
  top_annotation = cluster_annotation,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = combined_df$Integrated_Clusters,
  border = TRUE,
  show_heatmap_legend = FALSE
)


color_legend <- Legend(
  title = "Activity",
  at = seq(0, 1, by = 0.25),
  labels = seq(0, 1, by = 0.25),
  col_fun = colorRamp2(c(0, 1), c("white", "red")),
  direction = "horizontal",
  title_position = "topcenter",
  border = TRUE
)

lgd1 <- Legend(labels = names(colors)[1:7],
               title = "Post1/Meta-Integrated Clusters",
               legend_gp = gpar(fill = unname(colors[1:7]), col = "black"),
               title_position = "topcenter",
               nrow = 2,
               border = TRUE,
               direction = "horizontal")

lgd2 <- Legend(labels = names(colors)[8:12],
               title = "Post1-CNV & Meta",
               legend_gp = gpar(fill = unname(colors[8:12]), col = "black"),
               title_position = "topcenter",
               nrow = 2,
               border = TRUE,
               direction = "horizontal")


pdf("Figures/Fig.4D.pdf", width = 7.5, height = 7.5, fonts = "ArialMT", colormodel = "srgb")

pushViewport(viewport(layout = grid.layout(2, 1, heights = c(9, 1))))


pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
pushViewport(viewport(layout = grid.layout(1, 3, widths = c(2, 3, 4))))

pushViewport(viewport(layout.pos.col = 1))
draw(color_legend, just = "center")
popViewport()

pushViewport(viewport(layout.pos.col = 2))
draw(lgd1, just = "center")
popViewport()

pushViewport(viewport(layout.pos.col = 3))
draw(lgd2, just = "center")
popViewport()


popViewport()
popViewport()
dev.off()




# TF activities on the UMAP plot.
selected_TF_Z <- as.data.frame(selected_TF_Z)
selected_TF_Z_MPC <- selected_TF_Z[, c("ATF3", "SNAI2", "KLF4", "EGR1", "ID1", "SPI1", "SPIB", "MEF2C")]
colnames(selected_TF_Z_MPC) <- paste0("TF_", c("ATF3", "SNAI2", "KLF4", "EGR1", "ID1", "SPI1", "SPIB", "MEF2C"))

obj_cancer <- do.call(AddMetaData, c(list(obj_cancer), as.list(selected_TF_Z_MPC)))

feature_names <- colnames(selected_TF_Z_MPC)
plots <- lapply(feature_names, function(feature) {
  FeaturePlot(obj_cancer, features = feature,
              cols = c("lightgray", "red"), order = TRUE, pt.size = 0.25) +
    NoAxes() +
    theme(text = element_text(family = "Arial", size = 12),
          legend.position = c(0.05, 0.25),
          legend.key.size = unit(0.2, 'inch'))
})

p <- gridExtra::grid.arrange(grobs = plots, nrow = 3)
ggsave(file = "Figures/Fig.4E.png", plot = p, dpi = 300, width = 8, height = 6.4)