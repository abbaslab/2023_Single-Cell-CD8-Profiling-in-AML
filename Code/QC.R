############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: QC of data
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#Doublet score
pdf(paste0(file_prefix,"Doublet_Score.pdf", height=5, width=5))
FeaturePlot(cd8, "Doublet_Score") + 
  scale_color_gradient(low="gray", high="purple",limits=c(0,1)) + 
  ggtitle("Doublet Score") & NoAxes()
dev.off()

# Filtering information/violin plots
pdf(paste0(file_prefix,"features_counts.pdf",  height=4, width=9))
VlnPlot(cd8, c("nFeature_RNA", "nCount_RNA"), pt.size=0, cols=anno_colors$class13) & 
  xlab(NULL)
dev.off()

#Mitochondrial and Ribosomal percentage
pdf(paste0(file_prefix,"mito_ribo_ratios.pdf",  height=4, width=9))
VlnPlot(cd8, c("mitoRatio", "riboRatio"), pt.size=0, cols=anno_colors$class13) & 
  ylim(0,1)
dev.off()

#Cell Cycle Score
pdf(paste0(file_prefix,"cell_cycle.pdf", height=4, width=9))
VlnPlot(cd8, c("S.Score", "G2M.Score"), pt.size=0, cols=anno_colors$class13) &
  xlab(NULL)
dev.off()

#CD3,CD4,CD8 Vln plots
pdf(paste0(file_prefix,"CD3_CD8.pdf",  height=4, width=10))
VlnPlot(cd8, c("CD3D","CD3E","CD3G","CD8A", "CD8B"), cols=anno_colors$class13, pt.size=0) & 
  NoLegend() & 
  theme(axis.title=element_blank(), axis.text=element_text(size=10))
dev.off()

#Cells per Individual
a <- data.frame(table(cd8$Patient))
my.colors <- paletteDiscrete(sort(unique(a$Var1)),set = "stallion2")
pdf(paste0(file_prefix,"cd8_patient_bar.pdf", width=9, height=6))
ggbarplot(a, x="Var1", y="Freq",fill="Var1", palette = my.colors) +
  theme(axis.text.x=element_text(angle=90, size=8), legend.position = "none") + 
  xlab(NULL) + 
  ylab("# of Cells")
dev.off()

#Cells by Patient Type
a <- data.frame(table(cd8$type))
pdf(paste0(file_prefix,"cd8_type_bar.pdf", width=4, height=4))
ggbarplot(a, x="Var1", y="Freq",fill="Var1", palette = anno_colors$type) +
  theme(axis.text.x=element_text(angle=90, size=8)) + 
  xlab(NULL) + 
  ylab("# of Cells")
dev.off()