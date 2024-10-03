#---------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#---------------------------------------------------------

# 5. Cell-to-cell communication

library(SeuratData)
library(Connectome)
library(dplyr)

# obj: SeuratObject of Post1 cancer cells

# convert a v5 assay to a v3 assay
obj[["RNA3"]] <- as(object = obj[["RNA"]], Class = "Assay")
DefaultAssay(obj) <- "RNA3"

conn <- obj %>%
  NormalizeData() %>%
  { ScaleData(., features = union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol, Connectome::ncomms8866_human$Receptor.ApprovedSymbol) %>% intersect(rownames(.))) } %>%
  CreateConnectome(assay = "RNA3", species = 'human', min.cells.per.ident = 75, p.values = FALSE, calculate.DOR = FALSE) %>%
  FilterConnectome(min.pct = 0.1, min.z = 0.25, remove.na = TRUE)



## Centrality analysis
cols = c("#C2ADC0", "#FFED6F", "#C7D98C", "#DED7DA")

# modify the source code
Centrality <- function(connectome,
                       cols.use = NULL,
                       weight.attribute = 'weight_sc',
                       min.z = NULL,
                       normalize = T,
                       group.by = 'mode',...){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)
  
  if(weight.attribute == 'weight_sc' & is.null(min.z)){
    message("\nWeight attribute is 'weight_sc', recommend also setting min.z = 0 to avoid negative ligand and receptor scores")
  }
  
  # Store
  master <- connectome
  
  # Filter as demanded (remove NAs at minimum)
  master_sub <- FilterConnectome(master,min.z = min.z,remove.na = T,...)
  
  # Set up to plot ModalDotPlot
  
  modes <- as.character(unique(master_sub$mode))
  ligands <- as.character(unique(master_sub$ligand))
  recepts <- as.character(unique(master_sub$receptor))
  genes <- as.character(union(unique(master_sub$ligand),unique(master_sub$receptor)))
  pairs <- as.character(unique(master_sub$pair))
  
  cells <- as.character(unique(union(master_sub$source,master_sub$target)))
  
  # Determine which type of network plot:
  if (group.by == 'mode'){
    groups <- modes
  }
  if (group.by == 'ligand'){
    groups <- ligands
  }
  if (group.by == 'receptor'){
    groups <- recepts
  }
  if (group.by == 'gene'){
    groups <- genes
  }
  if (group.by == 'mechanism'){
    groups <- pairs
  }
  
  # Process data
  df <- data.frame()
  
  for (i in 1:length(groups)){
    if (group.by == 'mode'){
      temp <- subset(master_sub,mode == groups[[i]])
    }
    if (group.by == 'ligand'){
      temp <- subset(master_sub,ligand == groups[[i]])
    }
    if (group.by == 'receptor'){
      temp <- subset(master_sub,receptor == groups[[i]])
    }
    if (group.by == 'gene'){
      temp <- subset(master_sub,ligand == groups[[i]] | receptor == groups[[i]])
    }
    if (group.by == 'mechanism'){
      temp <- subset(master_sub,pair == groups[[i]])
    }
    
    # Make network for analysis
    net <- igraph::graph_from_data_frame(temp, directed = T)
    hub <- hub_score(net,weights = temp[,weight.attribute], scale = T)$vector
    auth <- authority_score(net,weights = temp[,weight.attribute], scale = T)$vector
    total.edgeweight <- sum(temp[,weight.attribute])
    # Get cumulative source/sink edgeweights per cell type
    for (j in 1:length(cells)){
      temp2 <- subset(temp,source == cells[[j]])
      wt.source <- sum(temp2[,weight.attribute])
      temp2 <- subset(temp,target == cells[[j]])
      wt.sink <- sum(temp2[,weight.attribute])
      if (normalize == T){
        wt.source <- wt.source/total.edgeweight
        wt.sink <- wt.sink/total.edgeweight
      }
      # Compile info for single network
      row <- data.frame(
        group = groups[[i]],
        cells = cells[[j]],
        hub.score = hub[cells[[j]]],
        auth.score = auth[cells[[j]]],
        wt.source = wt.source,
        wt.sink = wt.sink,
        row.names = NULL)
      
      df <- rbind(df,row)
    }
  }
  
  #Alphabetize
  df$group <- factor(df$group,levels = sort(as.character(unique(df$group)),decreasing = T))
  
  # Plots
  p1 <- ggplot(df,aes(x=group,y=wt.source,color = as.factor(cells)))+
    geom_point(size = df$hub.score*10,alpha = 0.6)+
    coord_flip()+
    theme_dark()+
    theme(legend.position="none") + ggtitle('Outgoing Centrality')+
    geom_text(data=df %>% group_by(group) %>% top_n(1,hub.score),aes(group,wt.source,label=cells))+
    ylab('Outgoing Edgeweight by Cell Type')
  if (normalize == T){
    p1 <- p1+
      ylab('Outgoing Edgeweight Fraction by Cell Type')+
      ylim(0,1)
  }
  
  p2 <- ggplot(df,aes(group,wt.sink,color = as.factor(cells)))+
    geom_point(size = df$auth.score*10,alpha = 0.6)+
    coord_flip()+
    theme_dark()+
    theme(legend.position="none") + ggtitle('Incoming Centrality')+
    geom_text(data=df %>% group_by(group) %>% top_n(1,auth.score),aes(group,wt.sink,label=cells))+
    ylab('Incoming Edgeweight by Cell Type')
  if (normalize == T){
    p2 <- p2+
      ylab('Incoming Edgeweight Fraction by Cell Type')+
      ylim(0,1)
  }
  
  # Define vertical axis label
  if (group.by == 'mode'){
    p1 <- p1 + xlab('Network (by Family)')
    p2 <- p2 + xlab('Network (by Family)')
  }
  if (group.by == 'ligand'){
    p1 <- p1 + xlab('Network (by Ligand)')
    p2 <- p2 + xlab('Network (by Ligand)')
  }
  if (group.by == 'receptor'){
    p1 <- p1 + xlab('Network (by Receptor)')
    p2 <- p2 + xlab('Network (by Receptor)')
  }
  if (group.by == 'gene'){
    p1 <- p1 + xlab('Network (by Gene)')
    p2 <- p2 + xlab('Network (by Gene)')
  }
  if (group.by == 'mechanism'){
    p1 <- p1 + xlab('Network (by Mechanism)')
    p2 <- p2 + xlab('Network (by Mechanism)')
  }
  
  # Modify colors if desired
  if (!is.null(cols.use)){
    p1 <- p1 + scale_colour_manual(values = cols.use)
    p2 <- p2 + scale_colour_manual(values = cols.use)
  }
  # Put legend on bottom
  legend <- get_legend(
    p1 +
      guides(color = guide_legend(nrow = 2,byrow=TRUE,override.aes = list(size=5))) +
      theme(legend.position = "bottom")
  )
  # Assemble plot
  plot.top <- plot_grid(p1, p2,nrow = 1)
  return(plot_grid(plot.top,legend,ncol = 1,rel_heights = c(1, .1)))
}


p <- Centrality(conn,
                cols.use = cols,
                modes.include = NULL,
                min.z = NULL,
                weight.attribute = 'weight_sc',
                group.by = 'mode')


ggsave(file="Figures/Fig.S8A_Connectome_Centrality.png", plot=p, dpi=300, bg="white",
       units="in", width=9, height=6)




## Exploring connectomic data using CircosPlot
test <- conn
test <- data.frame(test %>% group_by(vector) %>% top_n(5,weight_sc))

cells.of.interest <- c("Post1-CNV1", "Post1-CNV2", "Post1-CNV3", "Post1-CNV4")
grid.col <- as.vector(c("#C2ADC0", "#FFED6F", "#C7D98C", "#DED7DA"))
names(grid.col) <- c("Post1-CNV1", "Post1-CNV2", "Post1-CNV3", "Post1-CNV4")


# CircosPlot
# modify the source code
CircosPlot <- function(connectome,
                       weight.attribute = 'weight_sc',
                       cols.use = NULL,
                       min.z = NULL,
                       lab.cex = 1,
                       balanced.edges = T,
                       edge.color.by.source = T,
                       small.gap = 1,
                       big.gap = 10,
                       title = NULL,...){
  library(tidyverse)
  library(circlize)
  library(dplyr)
  library(scales)
  library(ComplexHeatmap)
  
  # If (weight.attribute != 'weight_norm'){
  if (weight.attribute == 'weight_sc' & is.null(min.z)){
    connectome <- FilterConnectome(connectome, remove.na = T,min.z = 0,...)
  }else{
    connectome <- FilterConnectome(connectome,remove.na = T,min.z = min.z,...)
  }
  #}
  # Pull the dataframe of interest for plotting and format with weight as third column
  connectome$lig.stash <- as.character(connectome$ligand)
  connectome$rec.stash <- as.character(connectome$receptor)
  df <- data.frame(connectome %>% select(ligand,receptor))
  df$ligand <- make.unique(as.character(df$ligand))
  df$receptor <- make.unique(as.character(df$receptor))
  #df$weight <- connectome[,weight.attribute]
  temp <- connectome[,!colnames(connectome) %in% colnames(df)]
  df <- cbind(df,temp)
  
  # Squash ligands back together to single name if they are duplicates (different edges on same cell type)
  for (i in 1:length(unique(df$lig.stash))){
    temp <- subset(df,lig.stash == unique(df$lig.stash)[i])
    for (j in 1:length(unique(temp$source))){
      temp2 <- subset(temp,source == unique(temp$source)[j])
      dummy <- paste(rep(' ',j-1),collapse = '') # Add number of spaces corresponding to number of unique sources
      df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
    }
    #if(length(unique(temp$source)) == 1){
    #  df[rownames(temp),]$ligand <- as.character(temp$lig.stash)
    #}
  }
  
  # Squash receptors back together to single name if they are duplicates (different edges on same cell type)
  for (i in 1:length(unique(df$rec.stash))){
    temp <- subset(df,rec.stash == unique(df$rec.stash)[i])
    for (j in 1:length(unique(temp$target))){
      temp2 <- subset(temp,target == unique(temp$target)[j])
      dummy <- paste(rep(' ',j-1),collapse = '') # Add number of spaces corresponding to number of unique targets
      df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
    }
    #if(length(unique(temp$target)) == 1){
    #  df[rownames(temp),]$receptor <- as.character(temp$rec.stash)
    #}
  }
  
  # Squash ligands back together, by cell type, if they are expressed on multiple cell types
  #temp <- subset(df,ligand != lig.stash) # this is a problem
  #if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$source))){
  #    temp2 <- subset(temp,source == unique(temp$source)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique sources
  #    df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
  #  }
  #}
  # Squash receptors back together, by cell type, if they are expressed on multiple cell types
  #temp <- subset(df,receptor != rec.stash) # this is a problem
  #if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$target))){
  #    temp2 <- subset(temp,target == unique(temp$target)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique targets
  #    df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
  #  }
  #}
  
  #Establish ordering (order) so that genes are grouped nicely by celltype
  source.order <- df[order(df$source), ]
  target.order <- df[order(df$target), ]
  source.order.un <- unique(source.order[,c('ligand','source')])
  target.order.un <- unique(target.order[,c('receptor','target')])
  
  source.order$id <- 1:nrow(source.order)
  target.order$id <- 1:nrow(target.order)
  source.order.un$id <- 1:nrow(source.order.un)
  target.order.un$id <- 1:nrow(target.order.un)
  
  sector.order.un <- c(as.character(source.order.un$ligand),
                       as.character(target.order.un$receptor))
  
  # Coloring setup
  if (is.null(cols.use)){
    nodes <- as.character(unique(union(df$source,df$target)))
    cols.use <- hue_pal()(length(nodes))
    names(cols.use) <- nodes
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }else{
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }
  
  
  # Map to get ligand colorings (edges)
  map <- base::merge(source.order, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.edges <- as.character(map$cols.use)
  names(lig.cols.edges) <- map$ligand
  
  # Map to get receptor colorings (edges) # this does not work
  map <- base::merge(target.order, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.edges <- as.character(map$cols.use)
  names(rec.cols.edges) <- map$receptor
  
  # Map to get ligand colorings (sectors)
  map <- base::merge(source.order.un, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.sect <- as.character(map$cols.use)
  names(lig.cols.sect) <- map$ligand
  
  # Map to get receptor colorings (sectors)
  map <- base::merge(target.order.un, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.sect <- as.character(map$cols.use)
  names(rec.cols.sect) <- map$receptor
  
  # Make sector colors (grid.cols)
  sectors <- c(source.order.un$ligand,target.order.un$receptor)
  sector.cols <- c(as.character(lig.cols.sect),as.character(rec.cols.sect))
  
  # Plotting
  # Decide edge order and edge color order
  if (edge.color.by.source == T){
    edge.color <- lig.cols.edges
    df.plot <- source.order
  }else{
    edge.color <- rec.cols.edges
    df.plot <- target.order
  }
  # Decide weight attributes and balanced vs. not
  if (weight.attribute == 'weight_norm'){
    if (balanced.edges == T){
      df.plot <- df.plot[,c('ligand','receptor','weight_norm')]
    }else{
      df.plot <- df.plot[,c('ligand','receptor','ligand.expression','recept.expression')]
    }
  }
  if (weight.attribute == 'weight_sc'){
    if (balanced.edges == T){
      df.plot <- df.plot[,c('ligand','receptor','weight_sc')]
    }else{
      df.plot <- df.plot[,c('ligand','receptor','ligand.scale','recept.scale')]
    }
  }
  #if (weight.attribute == 'score'){
  #  if (balanced.edges == T){
  #    df.plot <- df.plot[,c('ligand','receptor','score')]
  #  }else{
  #    df.plot <- df.plot[,c('ligand','receptor','ligand.norm.lfc','recept.norm.lfc')]
  #  }
  #}
  
  circos.clear()
  #circos.par(gap.degree = gap.degree)
  chordDiagram(df.plot,
               order = sector.order.un,
               col = edge.color,
               grid.col = sector.cols,
               directional = 1,
               direction.type = "arrows",
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = 1,
               small.gap = small.gap,
               big.gap = big.gap)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
    #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  # Make and add legend
  legend <- Legend(at = as.character(unique(union(df$source,df$target))),
                   type = "grid",
                   legend_gp = gpar(fill = as.character(cols.use[as.character(unique(union(df$source,df$target))),]$cols.use)),
                   title_position = "topleft",
                   title = "Cell Type", 
                   labels_gp = gpar( fontsize = 15))
  draw(legend, x = unit(10, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
  if(!is.null(title)){title(title)}
  
  p1.base <- recordPlot()
  return(p1.base)
}


P3G <- CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 1,title = 'Ligand vs. receptor expression (from scaled slot)', cols.use = grid.col)

pdf("Figures/Fig.3G.pdf")
P3G
dev.off()


# Network plot
features <- c('VEGFA', 'FGF2', 'COL1A2', 'LAMC2')

plots <- lapply(features, function(feature) {
  NetworkPlot(test, features = feature, min.pct = 0.1, weight.attribute = 'weight_sc', include.all.nodes = TRUE,
              title = paste("Network plot involving", feature), cols.use = cols)
})

pdf("Figures/Fig.S8B.pdf", width = 3, height = 3)
invisible(lapply(plots, print)) 
dev.off()

