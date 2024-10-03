#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 7.Expression correlations

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(Rmagic)

set.seed(1)


# Prepare magic data
prepare_magic_data <- function(obj) as.matrix(Rmagic::magic(t(as.matrix(obj[["integrated"]]@data)))$result)

# obj_cancer: integrated cancer clusters of Post1 and Meta (Fig. 4A).
Idents(obj_cancer) <- "orig.ident"
obj_list <- list(
  Post1 = subset(obj_cancer, idents = "Post1"),
  Meta = subset(obj_cancer, idents = "Meta")
)

Idents(obj_cancer) <- "Post1_CNV1_Meta"
obj_list <- c(obj_list, list(
  Post1_CNV1 = subset(obj_cancer, idents = "Post1-CNV1"),
  Post1_CNV2 = subset(obj_cancer, idents = "Post1-CNV2"),
  Post1_CNV3 = subset(obj_cancer, idents = "Post1-CNV3"),
  Post1_CNV4 = subset(obj_cancer, idents = "Post1-CNV4")
))


colors <- c("#66C2A5", "#8DA0CB", "#C2ADC0", "#FFED6F", "#C7D98C", "#DED7DA")
names(colors) <- names(obj_list)

gene_pairs <- list(
  c("ADM", "CALCA"),
  c("ADM", "CALCRL"),
  c("CALCA", "CALCRL"),
  c("CALCA", "RAMP1"),
  c("ADM", "RAMP2")
)

plot_list <- unlist(
  lapply(names(obj_list), function(name) {
    data <- prepare_magic_data(obj_list[[name]])
    lapply(gene_pairs, function(pair) {
      correlation <- cor(data[, pair[1]], data[, pair[2]], method = "pearson")
      print(correlation)
      ggplot(data = as.data.frame(data), aes_string(x = pair[1], y = pair[2])) +
        geom_point(color = colors[name], size = 0.5) +
        labs(x = pair[1], y = pair[2]) +
        theme_classic() +
        theme(
          text = element_text(family = "Arial", size = 16),
          axis.text = element_text(size = 14)
        ) +
        annotate("text", x = mean(range(data[, pair[1]])), y = max(data[, pair[2]]) * 1.05,
                 label = paste0("Correlation: ", round(correlation, 2)), size = 5, color = "black", vjust = 0.5)
    })
  }),
  recursive = FALSE
)


plot_grid_list <- list(
  plot_grid(plotlist = plot_list[1:5], ncol = 1),
  plot_grid(plotlist = plot_list[6:10], ncol = 1),
  plot_grid(plotlist = plot_list[11:15], ncol = 1),
  plot_grid(plotlist = plot_list[16:20], ncol = 1),
  plot_grid(plotlist = plot_list[21:25], ncol = 1),
  plot_grid(plotlist = plot_list[26:30], ncol = 1)
)


P10A <- plot_grid(plot_grid_list[[1]], plot_grid_list[[2]], align = "hv", ncol = 1)
ggsave(file = "Figures/Fig.S10A_cor.png", plot = P10A, dpi = 300, width = 14.4, height = 4.8)


P10B <- plot_grid(plot_grid_list[[3]], plot_grid_list[[4]], plot_grid_list[[5]], plot_grid_list[[6]], align = "hv", ncol = 1)
ggsave(file = "Figures/Fig.S10B_cor.png", plot = P8, bg = "white", dpi = 300, width = 9.6, height = 9.6)