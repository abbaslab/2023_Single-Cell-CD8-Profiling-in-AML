############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Figure 3
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: 
############################################

file_prefix <- ""
source(paste0(file_prefix,"CD8_Functions.R"))
cd8 <- readRDS(paste0(file_prefix,"cd8.rds"))

library(ComplexHeatmap); library(SCENIC)
library(patchwork)

########## A: Monocle3 Trajectory ##########
library(monocle3)
library(SeuratWrappers)

cd8.cds <- as.cell_data_set(cd8)
cd8.cds <- preprocess_cds(cd8.cds, num_dim = 100)
cd8.cds <- align_cds(cd8.cds, alignment_group = "orig.ident")
cd8.cds <- reduce_dimension(cd8.cds)
plot_cells(cd8.cds, color_cells_by="class13")

cd8.cds <- cluster_cells(cds=cd8.cds)
cd8.cds <- learn_graph(cd8.cds, close_loop = FALSE)
plot_cells(cd8.cds, color_cells_by = "class13", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)

cd8.cds <- order_cells(cd8.cds)
monocle_umap <- plot_cells(cd8.cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves=TRUE)

pdf(paste0(file_prefix,"cd8_monocle_trajectory.pdf"), height=6, width=8)
monocle_umap & NoAxes()
dev.off()

########## B/C: Monocle 3 Group/Type #########

pdf(paste0(file_prefix,"cd8_monocle_type.pdf"), height=6, width=8)
plot_cells(cd8.cds, color_cells_by = "type", label_branch_points = FALSE, label_leaves=FALSE, label_cell_groups = FALSE) + 
  scale_color_manual(values=anno_colors$type, name="Patient Type") & NoAxes()
dev.off()

pdf(paste0(file_prefix,"cd8_monocle_group.pdf"), height=6, width=8)
plot_cells(cd8.cds, color_cells_by = "group", label_branch_points = FALSE, label_leaves=FALSE, label_cell_groups = FALSE) + 
  scale_color_manual(values=anno_colors$group, name="Cell Type") & NoAxes()
dev.off()

########## D: Monocle Reduction Feature Plots ##########
pdf(paste0(file_prefix,"CX3CR1_expression.pdf"))
FeaturePlot(cd8, "CX3CR1", reduction="MON", max.cutoff = 2.5) + scale_color_viridis(option="H") & NoAxes() & NoLegend()
dev.off()

pdf(paste0(file_prefix,"KLRG1_expression.pdf"))
FeaturePlot(cd8, "KLRG1", reduction="MON", max.cutoff = 3) + scale_color_viridis(option="H") & NoAxes() & NoLegend()
dev.off()

pdf(paste0(file_prefix,"CXCR4_expression.pdf"))
FeaturePlot(cd8, "CXCR4", reduction="MON") + scale_color_viridis(option="H") & NoAxes() & NoLegend()
dev.off()

pdf(paste0(file_prefix,"FOSB_expression.pdf"))
FeaturePlot(cd8, "FOSB", reduction="MON", max.cutoff = 4) + scale_color_viridis(option="H") & NoAxes() & NoLegend()
dev.off()

########## E: Correlation Plots ##########

# obtain expression values of the genes to test from the Seurat object
expression_matrix <- FetchData(cd8, vars = c("CX3CR1", "CXCR4"))
expression_matrix <- round(expression_matrix , 1)
expression_matrix$class13 <- cd8$class13
expression_matrix$type <- cd8$type
expression_matrix$sample <- cd8$Patient

#Gene Correlation
ColAssign <- function(Var){
  require(ggthemes);require(RColorBrewer)
  pal <- tableau_color_pal(palette = "Classic 20",direction = 1,type="regular")
  if (length(Var) > 20) {
    palOut <- colorRampPalette(pal(20))(length(Var))
    names(palOut) <- Var
  } else if (length(Var) == 20) {
    palOut <- pal(20)
    names(palOut) <- Var
  } else if (length(Var) < 20) {
    palOut <- pal(20)
    palOut <- setdiff(palOut,c("#7f7f7f","#c7c7c7"))# remove grey colors
    #palOut <- sample(palOut)
    palOut <- c(palOut,c("#7f7f7f","#c7c7c7"))
    palOut <- palOut[1:length(Var)]
    names(palOut) <- Var
  }
  return(palOut)
}
df <- expression_matrix %>% group_by(class13,type,.drop = FALSE) %>% summarise(CXCR4=mean(CXCR4),CX3CR1=mean(CX3CR1))
table(df$type)
palOut <- ColAssign(unique(df$class13))

pdf(paste0(file_prefix,"CX3CR1_CXCR4_correlation.pdf"))
ggplot(df, aes(y=CX3CR1, x=CXCR4,shape=type)) + geom_point(aes(size=2,color=class13))+
  geom_smooth(se=FALSE,method = lm, linetype="dashed", size=0.5) + scale_color_manual(values=anno_colors$class13, name="Clusters")+ 
  theme_classic() + stat_cor(label.x=1.5) + facet_wrap(.~type)
dev.off()

cor.test(df$CX3CR1, df$CXCR4)



########## F: Type Histogram over Pseudotime ##########

pdf(paste0(file_prefix,"pseudotime_histogram.pdf"))
ggplot(cd8@meta.data, aes(x=pseudotime, fill=type)) + geom_histogram(alpha=0.4, bins=20) + theme_classic() +
  scale_fill_manual(values=anno_colors$type, name="Patient Type") +
  labs(x="Pseudotime", y="Frequency")
dev.off()


########## G: SCENIC ##########
setwd(paste0(file_prefix,"/scenic/"))
scenicOptions <- readRDS(paste0(file_prefix,"scenicOptions.Rds"))
load(paste0(file_prefix,"cellInfo.rda"))
cellInfo$class13 =ifelse(cellInfo$class12 %in% c("NKT TIGIT+", "GZMK+EM EOMES+", "CTL CX3CR1+"),"Group1","Group2")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
AUC.Mat <- getAUC(regulonAUC)

rownames(AUC.Mat)=unlist(lapply(rownames(AUC.Mat), function(x) {unlist(strsplit(x," "))[[1]]}))
rownames(AUC.Mat)=unlist(lapply(rownames(AUC.Mat), function(x) {unlist(strsplit(x,"_"))[[1]]}))
#heatmap for Field
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$class12),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c("NKT TIGIT+","GZMK+EM TGFB1+","CTL TGFB1+","NKT IFNG+","CTL CX3CR1+","GZMK+EM EOMES+")]
rownames(regulonActivity_byCellType_Scaled)=unlist(lapply(rownames(regulonActivity_byCellType_Scaled), function(x) {unlist(strsplit(x," "))[[1]]}))
rownames(regulonActivity_byCellType_Scaled)=unlist(lapply(rownames(regulonActivity_byCellType_Scaled), function(x) {unlist(strsplit(x,"_"))[[1]]}))

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=class13_cells[colnames(regulonAUC)])

pdf(paste0(file_prefix,"SCENIC_RSS.pdf", height=14, width=10))
plotRSS_oneSet(rss, setName = "NK-like1") + plotRSS_oneSet(rss, setName = "NK-like2") +
plotRSS_oneSet(rss, setName = "CTL1") + plotRSS_oneSet(rss, setName = "CTL2") +
plotRSS_oneSet(rss, setName = "EM1") + plotRSS_oneSet(rss, setName = "EM2") + plot_layout(ncol=2)
dev.off()

expr=as.data.frame(t(AUC.Mat))
expr$class13=cellInfo$class13

p=ggboxplot(data=expr,x="class13", y="TBX21", color="black", xlab ="",ylab="Regulon activity",
            fill = "class13",bxp.errorbar = FALSE) + 
  scale_fill_manual(values=c("gold2","deepskyblue"))+
  theme(axis.text.x = element_text(size=10, angle=45,hjust=1),axis.text.y = element_text(size=10),
        axis.title.y=element_text(size=12),legend.position = 'none')
library(rstatix)
wilcox_stat <- expr %>% group_by("class13") %>% pairwise_wilcox_test(TBX21~class13, p.adjust.method = "BH") %>% add_significance("p.adj") %>% add_xy_position(x="class13", fun="max", step.increase = 0.1)
pdf(paste0(file_prefix,"TBX21.pdf", height=5, width=3))
p + stat_pvalue_manual(wilcox_stat, label="p.adj.signif", tip.length = 0.0,hide.ns = T) 
dev.off()

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

pdf(paste0(file_prefix,"EOMES_AUC.pdf"))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, AUC.Mat, regulonAUC["EOMES (923g)",], plots=c("AUC"))
dev.off()

pdf(paste0(file_prefix,"TBX21_AUC.pdf"))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, AUC.Mat, regulonAUC["TBX21 (188g)",], plots=c("AUC"))
dev.off()

pdf(paste0(file_prefix,"JUN_AUC.pdf"))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, AUC.Mat, regulonAUC["JUN (50g)",], plots=c("AUC"))
dev.off()

pdf(paste0(file_prefix,"FOSB_AUC.pdf"))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, AUC.Mat, regulonAUC["FOSB (158g)",], plots=c("AUC"))
dev.off()
