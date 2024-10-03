#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 3. Pathway enrichment analysis

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)

set.seed(1)

### Seurat enrichment ###
library(enrichR)
# gene-set library
# https://maayanlab.cloud/Enrichr/#libraries

Idents(obj_cancer)<-"sample.id2" # integrated cancer clusters of all 4 samples (Fig. 1B, C, Supplementary Fig. 3B).

# modify the source code
DEenrichRPlot <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    balanced = TRUE,
    logfc.threshold = 0.25,
    assay = NULL,
    max.genes,
    test.use = 'wilcox',
    p.val.cutoff = 0.05,
    cols = NULL,
    enrich.database = NULL,
    num.pathway = 10,
    return.gene.list = FALSE,
    ...
) {
  enrichr.installed <- PackageCheck("enrichR", error = FALSE)
  if (!enrichr.installed[1]) {
    stop(
      "Please install the enrichR package to use DEenrichRPlot",
      "\nThis can be accomplished with the following command: ",
      "\n----------------------------------------",
      "\ninstall.packages('enrichR')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  if (is.null(x = enrich.database)) {
    stop("Please specify the name of enrichR database to use")
  }
  if (!is.numeric(x = max.genes)) {
    stop("please set max.genes")
  }
  assay <- assay %||% DefaultAssay(object = object)
  
  DefaultAssay(object = object) <- assay
  
  all.markers <- FindMarkers(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    only.pos = FALSE,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    assay = assay
  )
  
  pos.markers <- all.markers[all.markers[, 2] > logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
  
  if(nrow(pos.markers) == 0){
    message("No positive markers pass the logfc.thershold")
    pos.er <- c()
  }
  
  else{
    pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
    pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
    pos.er <- do.call(what = cbind, args = pos.er)
    pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
    pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
    pos.er <- pos.er[1:num.pathway, ]
    pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
    gene.list <- list(pos = pos.er)
  }
  
  if (isTRUE(x = balanced)) {
    neg.markers <- all.markers[all.markers[, 2] < -logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    Sys.sleep(1)
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er <- neg.er[1:num.pathway, ]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])
    
    if(isTRUE(length(neg.er$term) == 0) & isTRUE(length(pos.er == 0))){
      stop("No positive or negative marker genes identified")
    }
    
    else{
      if(isTRUE(length(neg.er$term) == 0)){
        
        gene.list <- list(pos = pos.er)
        
      }
      else{
        gene.list <- list(pos = pos.er, neg = neg.er)
      }
    }
    
  }
  if (return.gene.list) {
    return(gene.list)
  }
  
  if(nrow(pos.markers) == 0){
    message("No positive markers to plot")
    
    if (isTRUE(x = balanced)) {
      
      p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
        geom_bar(stat = "identity", fill = "dodgerblue") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                  size = 5,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      p <- p2
      
    }
    else{
      stop("Nothing to plot")
    }
  }
  
  else {
    p <- ggplot(data = pos.er, aes_string(x = "term", y = "log10pval")) +
      geom_bar(stat = "identity", fill = "orange") +
      coord_flip() + xlab("Pathway") +
      scale_fill_manual(values = cols, drop = FALSE) +
      ylab("-log10(pval)") +
      ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
      theme_classic() +
      geom_text(aes_string(label = "term", y = 0),
                size = 5,
                color = "black",
                position = position_dodge(1),
                hjust = 0)+
      theme(axis.title.y= element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    if (isTRUE(x = balanced)) {
      
      p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
        geom_bar(stat = "identity", fill = "dodgerblue") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                  size = 5,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      p <- p+p2
      
    }
  }
  
  return(p)
}


p_post <- DEenrichRPlot(
  obj_cancer,
  ident.1 = "Posttreatment",
  balanced = TRUE,
  logfc.threshold = 1,
  max.genes=50,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = c("orange", "dodgerblue"),
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,
  return.gene.list = TRUE,
  repel=TRUE
)&theme(text=element_text(family="Arial"),
        plot.title = element_text(size=14))

ggsave(file="Figures/Fig.1E.png", plot=p_post, bg="white", dpi=300, width=11,height=4)


## Clustergram ##
library(pheatmap)


markers <- FindMarkers(obj_cancer, ident.1 = "Posttreatment", logfc.threshold = 1, test.use = "wilcox", min.pct = 0.25)
filtered_markers <- markers[markers$p_val_adj < 0.05, ]

top_pos_genes <- rownames(filtered_markers[filtered_markers$avg_log2FC > 1, ])[1:min(50, sum(filtered_markers$avg_log2FC > 1))]
top_neg_genes <- rownames(filtered_markers[filtered_markers$avg_log2FC < -1, ])[1:min(50, sum(filtered_markers$avg_log2FC < -1))]


dbs <- c("MSigDB_Hallmark_2020")
pos_msigdb_results <- enrichr(top_pos_genes, dbs)$MSigDB_Hallmark_2020
neg_msigdb_results <- enrichr(top_neg_genes, dbs)$MSigDB_Hallmark_2020


create_gene_term_matrix <- function(terms, genes, enrich_results) {
  gene_term_matrix <- matrix(0, nrow = length(genes), ncol = length(terms), dimnames = list(genes, terms))
  enrich_data <- enrich_results[enrich_results$Term %in% terms, ]
  genes_list <- strsplit(enrich_data$Genes, ";")
  term_indices <- match(enrich_data$Term, terms)
  sapply(seq_along(genes_list), function(i) {
    gene_indices <- match(genes_list[[i]], genes)
    gene_term_matrix[gene_indices, term_indices[i]] <<- -log10(enrich_data$Adjusted.P.value[i])
  })
  return(gene_term_matrix[rowSums(gene_term_matrix, na.rm = TRUE) > 0, , drop = FALSE])
}

pos_terms <- factor(p_post$pos$Term, levels = unique(p_post$pos$Term))
neg_terms <- factor(p_post$neg$Term, levels = unique(p_post$neg$Term))

pos_gene_term_matrix <- create_gene_term_matrix(levels(pos_terms), top_pos_genes, pos_msigdb_results)
neg_gene_term_matrix <- create_gene_term_matrix(levels(neg_terms), top_neg_genes, neg_msigdb_results)


create_heatmap <- function(gene_term_matrix, terms, color) {
  Heatmap(
    as.matrix(gene_term_matrix),
    name = "Enrichment -log10(p-value)",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 45,
    col = colorRamp2(c(0, max(gene_term_matrix, na.rm = TRUE)), c("white", color)),
    rect_gp = gpar(col = "#999999"),
    column_names_gp = grid::gpar(fontsize = 14),
    heatmap_legend_param = list(direction = "horizontal")
  )
}

htmp_pos <- create_heatmap(pos_gene_term_matrix, levels(pos_terms), "orange")
htmp_neg <- create_heatmap(neg_gene_term_matrix, levels(neg_terms), "dodgerblue")


ae1_data <- as.matrix(AverageExpression(obj_cancer, features = rownames(filtered_markers), group.by = "orig.ident", assays = "RNA")$RNA)

annotation_col <- data.frame(Sample = factor(c("Pre", "Post1", "Post2", "Meta"), levels = c("Pre", "Post1", "Post2", "Meta")))
rownames(annotation_col) <- colnames(ae1_data)

ha <- HeatmapAnnotation(df = annotation_col, col = list(Sample = c("Pre" = "#E78AC3", "Post1" = "#66C2A5", "Post2" = "#FC8D62", "Meta" = "#8DA0CB")), show_legend = FALSE)
color_mapping <- viridis(10)

create_expression_heatmap <- function(ae_data, gene_term_matrix, ha) {
  scaled_data <- t(apply(ae_data[rownames(gene_term_matrix), , drop = FALSE], 1, function(x) (x - min(x)) / (max(x) - min(x))))
  Heatmap(
    scaled_data,
    name = "Scaled Expression",
    top_annotation = ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = color_mapping,
    column_names_side = "top",
    column_names_rot = 45,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 14),
    heatmap_legend_param = list(direction = "horizontal")
  )
}

ht_pos <- create_expression_heatmap(ae1_data, pos_gene_term_matrix, ha)
ht_neg <- create_expression_heatmap(ae1_data, neg_gene_term_matrix, ha)


pdf("Figures/Fig.1F.pdf", width = 6.6, height = 7.5, fonts = "ArialMT", colormodel = "srgb", useDingbats = FALSE)
draw(htmp_pos + ht_pos, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()

pdf("Figures/Fig.S3F.pdf", width = 6.6, height = 7.5, fonts = "ArialMT", colormodel = "srgb", useDingbats = FALSE)
draw(htmp_neg + ht_neg, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()