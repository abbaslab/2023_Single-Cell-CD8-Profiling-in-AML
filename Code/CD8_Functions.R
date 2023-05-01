############################################
## Project: AML scRNAseq
## Script Purpose: Define functions and marker lists that are going to be used in multiple scripts
## Date: 1/25/2022
## Author: Poonam Desai
## Notes: 
############################################

# Load required packages
library(Seurat); library(readxl); library(dplyr); library(ggplot2); library(clustree)
library(data.table); library(ggpubr); library(ArchR); library(harmony)
library(RColorBrewer); library(pheatmap); library(purrr)
library(ComplexHeatmap)
library(grid); library(gridExtra)

group.colors <- c("#FFBB78", "#C5B0D5", "#FF9896", "#AEC7E8", "#C7C7C7", "#DBDB8D")
names(group.colors) <- c("N", "Mem", "EM", "CTL", "MAIT", "Ex") 

class13.colors <- c("#F28E2B", "#B07AA1", "#D4A6C8", "#E15759", "#FF9D9A", "#D37295", "#FABFD2", 
                             "#4E79A7", "#A0CBE8", "#86BCB6", "#5FA2CE", "#499894",  "#59A14F", "#8CD17D", "#BAB0AC", "#F1CE63")

names(class13.colors) <- c("N", "Mem1", "Mem2", "EM1", "EM2", "INF1", "INF2", 
                           "CTL1", "CTL2", "CTL INF1", "CTL INF2", "CTL INF3", "NK-like1", "NK-like2", "MAIT", "Ex")

anno_colors = list(type=c("Healthy"= "seagreen", "NewlyDx"="dodgerblue", "RelRef"="firebrick3"), group=group.colors,
                   class13 = class13.colors)              


# Analyze the Seurat Object
analyze_seurat <- function(seurat_object, number_features=2500, res=0.8, md=0.01, k=20){
  TCR.genes <- grep("^TR[ABGD][VJ]",rownames(seurat_object),value = T)
  seurat_object <- FindVariableFeatures(seurat_object,nfeatures = number_features)
  var.genes <- VariableFeatures(seurat_object)
  var.genes <- setdiff(var.genes,TCR.genes)
  VariableFeatures(seurat_object) <- var.genes
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(object = seurat_object, verbose = FALSE,features = VariableFeatures(seurat_object),npcs=50)
  seurat_object <- RunHarmony(seurat_object, group.by.vars="orig.ident",assay.use ="RNA", max.iter.cluster=500)
  seurat_object <- FindNeighbors(seurat_object, dims=1:30,reduction = "harmony",k.param = k)
  seurat_object <- FindClusters(seurat_object, resolution = seq(0,1,0.1), random.seed=123)
  seurat_object <- FindClusters(object = seurat_object,resolution=res,random.seed=123,graph.name = 'RNA_snn')
  seurat_object <- RunUMAP(seurat_object,reduction = "harmony",seed.use = 123,dims=1:50,
                           a=0.5,b=1.5, verbose=FALSE)
  return(seurat_object)
}

# Calculate group means for heatmap
groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}


#Dot Plots
dot_plot <- function(seurat_object, features, lim=NULL, feats=NULL, group, min=-1.5, max=1.5){
  dot_plot <- DotPlot(object = seurat_object, features=features, group.by=group, col="RdYlBu", 
                      col.min = min, col.max = max, scale.max = 50) + 
    theme(legend.position="none",axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y=element_text(size=10),axis.text.y=element_text(size=10),
          axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    labs(x=feats, y=NULL) +  scale_y_discrete(limits=lim) + coord_flip()
  return(dot_plot)
}

dot_plot_leg <- function(seurat_object, features, lim=NULL, group="class13", feats=NULL, min=-1.5, max=1.5){
  dot_plot <- DotPlot(object = seurat_object, features=features, group.by=group,col="RdYlBu", col.min = min, col.max=max, scale.max = 50) + 
    theme(legend.position="right",axis.ticks.x = element_blank(), axis.text.x = element_text(size=10, angle=60, hjust=1), 
          axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    labs(x=feats, y=NULL) +  scale_y_discrete(limits=lim) + coord_flip()
  return(dot_plot)
}

#function for a grouped dot plot (no legend)
group_dot_plot<- function(seurat_object, features=list(feat_naive,feat_ctl,feat_ex),lim=NULL, group="class9", 
                          feats=NULL, min=-1.5, max=1.5){
  feat_len <- length(features)
  plot_list <- vector("list", feat_len)
  n <- vector("list", feat_len)
  g_build <- vector("list", feat_len)
  g_table <- vector("list", feat_len)
  for(i in 1:(feat_len-1)){
    plot_list[[i]]<- dot_plot(seurat_object, features[[i]], lim=lim, group=group, feats=feats[[i]], min=min, max=max)
    n[i] <- length(features[[i]])
    g_build[[i]] <- ggplot_build(plot_list[[i]])
    g_table[[i]] <- ggplot_gtable(g_build[[i]])
    if (i==1){ g <- g_table[[1]]} else {g <- gtable_rbind(g,g_table[[i]])}
  }
  plot_list[[feat_len]] <- dot_plot_x(seurat_object, features[[feat_len]], lim=lim, group=group, feats=feats[[feat_len]], min=min, max=max)
  n[feat_len] <- length(features[[feat_len]])
  g_build[[feat_len]] <- ggplot_build(plot_list[[feat_len]])
  g_table[[feat_len]] <- ggplot_gtable(g_build[[feat_len]])
  
  g <- gtable_rbind(g, g_table[[feat_len]])
  panels <- g$layout$t[grepl("panel", g$layout$name)]
  a <- unlist(n)
  g$heights[panels] <- unit(a,"null")
  return(g)
}


bk <- seq(-2,2,by=0.05)
colour_bk <- c(colorRampPalette(c("#4575b4","#e0f3f8"))(35),
               colorRampPalette(c("#e0f3f8","#ffffbf"))(10),
               colorRampPalette(c("#ffffbf","#d73027"))(36))


theme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.ticks = element_line(color = "white", size  =  0.2),
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,   
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "none",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_rect(fill = NA, color = "white"),
      ##panel.grid.major = element_line(color = "grey35"),
      panel.grid.major = element_blank(),
      ##panel.grid.minor = element_line(color = "grey20"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white", face="bold"),
      strip.text.y = element_text(size = base_size, color = "white",angle = -90),
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")
    )
}

########## Cluster Membership ##########

cluster_membership <- function(x,y){
  a <- table(x,y)
  b <- t(a/rowSums(a))
  p <- pheatmap::pheatmap(mat = b, color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100), 
                          cluster_cols = FALSE, cellheight=15, cellwidth = 15)
  return(p)
}


########## Marker Lists ##########

canonical_marker <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CCR7", "LEF1", "SELL", "IL7R", "STAT1", 
                      "STAT3", "FOXP3", "MX1", "IFI44L", "IFIT3", "GZMK", "CD27", "LTB", "GNLY", "PRF1", "GZMB", 
                      "NKG7", "NCAM1","SLC4A10", "TRAV1-2")

