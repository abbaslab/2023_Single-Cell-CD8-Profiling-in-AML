############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for SCENIC Analysis
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#SCENIC
setwd(paste0(file_prefix,"/scenic/"))
scenicOptions <- readRDS(paste0(file_prefix,"scenicOptions.Rds"))
load(paste0(file_prefix,"cellInfo.rda"))
cellInfo$class13 =ifelse(cellInfo$class12 %in% c("NKT TIGIT+", "GZMK+EM EOMES+", "CTL CX3CR1+"),"Group1","Group2")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
AUC.Mat <- getAUC(regulonAUC)

rownames(AUC.Mat)=unlist(lapply(rownames(AUC.Mat), function(x) {unlist(strsplit(x," "))[[1]]}))
rownames(AUC.Mat)=unlist(lapply(rownames(AUC.Mat), function(x) {unlist(strsplit(x,"_"))[[1]]}))

#RSS
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

#TBX21 Plot
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

#TSNE for individual markers
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

#Top Regulons
tf=c() 
for (i in 1:ncol(regulonActivity_byCellType_Scaled)){
  sel.TFs <- data.frame(TF=rownames(AUC.Mat),p=NA,row.names = rownames(AUC.Mat),cluster=regulonActivity_byCellType_Scaled[,i],other=rowMeans(regulonActivity_byCellType_Scaled[,-c(i)]))
  for (j in sel.TFs$TF) {
    test.x <- wilcox.test(AUC.Mat[j,cellInfo$class12==colnames(regulonActivity_byCellType_Scaled)[i]],AUC.Mat[j,cellInfo$class12!=colnames(regulonActivity_byCellType_Scaled)[i]])
    sel.TFs[j,"p"] <- test.x$p.value
  }
  sel.TFs <- sel.TFs[order(sel.TFs$p),]
  sel.TFs$tumor_hi <- NA
  sel.TFs$tumor_hi[sel.TFs$TF %in% rownames(regulonActivity_byCellType_Scaled)[regulonActivity_byCellType_Scaled[,i]>0]] <- "Primitive"
  sel.TFs2 <- sel.TFs[!is.na(sel.TFs$tumor_hi),]
  sel.TFs2=sel.TFs2[sel.TFs2$cluster>0,]
  sel.TFs2=sel.TFs2[sel.TFs2$p<0.01,]
  sel.TFs2$difference=sel.TFs2$cluster - sel.TFs2$other
  sel.TFs2=sel.TFs2[order(sel.TFs2$difference,decreasing = T),]
  tf=c(tf,rownames(sel.TFs2)[1:10])
}


png("/Users/bofeiwang/Downloads/poonam_paper/scenic/class12_top15.png",res=200,width=8, height=10,units="in")
p=ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[tf,],name="Regulon activity",cluster_rows = F,
                          column_title=NULL,show_heatmap_legend = T,column_names_rot = 45,cluster_columns = F,
                          row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 16),column_names_centered = T)
draw(p, annotation_legend_side = 'top',heatmap_legend_side="right")
dev.off()  
