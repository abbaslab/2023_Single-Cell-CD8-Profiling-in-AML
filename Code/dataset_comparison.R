############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for comparison with other datasets
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#Load Datasets
vg_cd8 <- readRDS("~/OneDrive - Inside MD Anderson/LabMembers/Resources/datasets/RNA_sc/2019_vanGalen_Cell/Objects/vg_cd8.rds")
petti_cd8 <- readRDS("/Users/poonamdesai/Library/CloudStorage/OneDrive-InsideMDAnderson/LabMembers/Resources/datasets/RNA_sc/2021_Petti_BCD/Objects/petti_cd8.rds")

#Define Markers
zhang_ex <- c("PDCD1", "TOX", "CXCL13", "TIGIT", "CTLA4", "TNFRSF9", "HAVCR2", "LAG3")

inhibitors <- c("ADORA2A", "ARG1", "BTLA", "BTN2A2", "CD160", "CD200", "CD200R1", "CD274", "CD276",
                "CD44", "CD47", "CD96", "CEACAM1", "CTLA4", "CYBB", "HAVCR2", "HLA-DPA1", "HLA-G", "ICOSLG",
                "IDO1", "IL2RB", "KIR2DL1", "KIR2DL3", "KIR2DL4",  "KIR3DL1", "KIR3DL2", "KIR3DL3", 
                "KLRC1", "LAG3", "LAIR1", "LGALS3", "LGALS9", "LILRB1", "LILRB2", "LYN", "NECTIN2", "NRP1", "PDCD1",
                "PDCD1LG2", "PTPN1", "PVRIG", "REL", "SIGLEC7", "SIGLEC9", "SIRPA", "TDO2", "TIGIT", "TNFRSF14",
                "TNFRSF8", "VSIR")

cd8_signature <- c("RPL41", "MTRNR2L12", "TXNIP", "CD2", "PRF1", "PSMB9", "LIMD2", "LTB", "GIMAP1", "GIMAP4", 
                           "CX3CR1", "TRAF3IP3", "GIMAP7", "GIMAP5", "DENND2D", "CLIC1", "MYL12A", "S100A11", "RAC2", 
                           "PSMB10", "ARPC1B", "CORO1A", "PCBP1", "UCP2", "PLEK")


#Van Galen data with Zheng markers
mat <- vg@assays$RNA@data[zhang_ex,]
mat <- groupMeans(mat,unlist(vg$RNA_snn_res.1))
pdf("./Figures/Additional/zhang_ex_vg_t_exp.pdf", height=5, width=4)
Heatmap(t(mat), col = colour_bk, cluster_rows = FALSE, cluster_columns = TRUE,
        column_title = NULL,  rect_gp = gpar(col = "white", lwd = 1),
        width = nrow(mat)*unit(5, "mm"), height = ncol(mat)*unit(5, "mm"), row_names_gp =gpar(fontsize=8),
        column_names_gp =gpar(fontsize=8),heatmap_legend_param = list(title="Avg. Exp."))
dev.off()

#Van Galen immune inhibitory
mat <- vg_cd8@assays$RNA@data[inhibitors,]
mat <- t(scale(t(mat)))
mat <- groupMeans(mat,unlist(vg_cd8$RNA_snn_res.1))
mat <- mat[complete.cases(mat),]
pdf("./Figures/Additional/inh_vg_cd8_z.pdf", height=12, width=7)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, heatmap_legend_param = list(title="Z-score"), na_col = "black",
        column_title = NULL,  width = ncol(mat)*unit(5, "mm"), height = nrow(mat)*unit(5, "mm"))
dev.off()

#Petti data with Zheng markers
mat <- petti_cd8@assays$RNA@data[zhang_ex,]
mat <- groupMeans(mat,unlist(petti_cd8$RNA_snn_res.1))
pdf("./Figures/Additional/zhang_ex_petti_cd8_exp.pdf", height=4, width=4)
Heatmap(t(mat), col = colour_bk, cluster_rows = FALSE, cluster_columns = TRUE,
        column_title = NULL,  rect_gp = gpar(col = "white", lwd = 1),
        width = nrow(mat)*unit(5, "mm"), height = ncol(mat)*unit(5, "mm"), row_names_gp =gpar(fontsize=8),
        column_names_gp =gpar(fontsize=8),heatmap_legend_param = list(title="Avg. Exp."))
dev.off()

#Petti immune inhibitory
mat <- petti@assays$RNA@data[inhibitors,]
mat <- t(scale(t(mat)))
mat <- groupMeans(mat,unlist(petti$RNA_snn_res.1))
mat <- mat[complete.cases(mat),]
pdf("./Figures/Additional/inh_petti_cd8_z.pdf", height=12, width=7)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, heatmap_legend_param = list(title="Z-score"), na_col = "black",
        column_title = NULL,  width = ncol(mat)*unit(5, "mm"), height = nrow(mat)*unit(5, "mm"))
dev.off()

#CD8-derived score
petti_eff <- readRDS()

eff_df <- data.frame(petti_eff$patient, petti_eff$type, petti_eff$dysfunction)
colnames(eff_df) <- c("patient", "type", "score")
library(rstatix)
wilcox_stat_eff<-eff_df %>% group_by(patient) %>% rstatix::pairwise_wilcox_test(score ~ type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="patient", fun="max", step.increase = 0.1)
wilcox_sig_eff <- subset(wilcox_stat_eff, p.adj < 0.05)
eff_box <- ggboxplot(data = eff_df, x="patient", y="dysfunction",fill="type", outlier.shape = NA, remove=c("boxplot"),
                     size=0.3) + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1), legend.position = "top", text = element_text(size=6)) + 
  guides(fill=guide_legend(title="Type")) 

pdf("./Plots/patient_BH_boxplot.pdf", height=4, width=6)
eff_box +stat_pvalue_manual(wilcox_sig_eff, label="p.adj.signif", tip.length = 0.0)
dev.off()

eff_df_mean <- eff_df %>% group_by(type,patient) %>% summarise_at(vars(cd8_signature), list(name = mean))

dys <- eff_df_mean$name
names(dys) <- paste0(eff_dys_mean$patient,"_",eff_dys_mean$type)
dys <- dys[order(names(dys))]
names(dys) <- rep(c("NewDx","RelRef"),6)

jci <- read_xlsx("~/Downloads/159579-JCI-CMED-RV-3_sd_617420.xlsx", sheet = "Supplementary Table 2", skip=2)
jci <- jci[-1,]

table(dysfunction_signature %in% jci$`IED172 (RNA-seq)`)
table(dysfunction_signature %in% jci$`IED68 (NanoString)`)
table(dysfunction_signature %in% jci$`LASSO penalized regression, PI20 (NanoString)`)
table(dysfunction_signature %in% jci$`LASSO penalized regression, PI24 (RNA-seq)`)

IED172 <- jci$`IED172 (RNA-seq)`
IED172 <- IED172[1:172]
IED68 <- jci$`IED68 (NanoString)`
IED68 <- IED68[complete.cases(IED68)]
PI20 <- jci$`LASSO penalized regression, PI20 (NanoString)`
PI20 <- PI20[complete.cases(PI20)]
PI24  <- jci$`LASSO penalized regression, PI24 (RNA-seq)`
PI24 <- PI24[complete.cases(PI24)]

library(ggVennDiagram); library(patchwork)
x <- list("Dysfunction" = dysfunction_signature,"IED172"=IED172)
x1 <- ggVennDiagram(x, label_alpha = 0, edge_size = 0, label_size = 4) + 
  scale_fill_distiller(palette = "Reds", direction = 1) + scale_color_manual(values=c("black", "black")) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

x <- list("Dysfunction" = dysfunction_signature,"IED68"=IED68)
x2 <- ggVennDiagram(x, label_alpha = 0, edge_size = 0, label_size = 4) + 
  scale_fill_distiller(palette = "Reds", direction = 1) + scale_color_manual(values=c("black", "black")) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 

x <- list("Dysfunction" = dysfunction_signature, "PI20"=PI20)
x3 <- ggVennDiagram(x, label_alpha = 0, edge_size = 0, label_size = 4) + 
  scale_fill_distiller(palette = "Reds", direction = 1) + scale_color_manual(values=c("black", "black")) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

x <- list("Dysfunction" = dysfunction_signature,"PI24"=PI24)
x4 <- ggVennDiagram(x, label_alpha = 0, edge_size = 0, label_size = 4) + 
  scale_fill_distiller(palette = "Reds", direction = 1) + scale_color_manual(values=c("black", "black")) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 

x <- list("IED172"=IED172,"IED68"=IED68, "PI20"=PI20, "PI24"=PI24)
x5 <- ggVennDiagram(x, label_alpha = 0, edge_size = 0, label_size = 4) + 
  scale_fill_distiller(palette = "Purples", direction = 1) + scale_color_manual(values=c("black", "black","black","black")) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pdf("./Figures/Additional/jci_sigs_venn.pdf", height=7,width=7)
(x1 + x2) /(x3 + x4) &  plot_layout(guides = 'collect')
dev.off()

jci_sigs <- list("IED172"=IED172, "IED68"=IED68, "PI20"=PI20, "PI24"=PI24) 

library(AUCell)
exprMatrix=GetAssayData(cd8, slot="counts")
cells_rankings <- AUCell_buildRankings(as.matrix(exprMatrix), plotStats=TRUE, splitByBlocks = TRUE) #build the rank for genes
cells_AUC <- AUCell_calcAUC(jci_sigs, cells_rankings)
score_mat <- t(SummarizedExperiment::assay(cells_AUC, 'AUC'))
colnames(score_mat)<-make.names(colnames(score_mat))
cd8=AddMetaData(cd8, as.data.frame(score_mat))

rm(exprMatrix, cells_rankings,cells_AUC,score_mat)

cd8_eff <- subset(cd8, subset = class13 %in% c("EM1", "EM2", "CTL1", "CTL2", "NK-like1", "NK-like2"))
cd8_eff$groupEff <- NA
cd8_eff$groupEff[colnames(cd8_eff)[which(cd8_eff$class13 %in% c("EM1", "CTL1", "NK-like1"))]] <-"Group1"
cd8_eff$groupEff[colnames(cd8_eff)[which(cd8_eff$class13 %in% c("EM2", "CTL2", "NK-like2"))]] <-"Group2"

c <- data.frame(cd8_eff$dysfunction, cd8_eff$IED172, cd8_eff$Patient, cd8_eff$class13)
colnames(c) <- c("Dysfunction", "IED172", "Patient", "class13")
c2 <- c %>% group_by(Patient,class13) %>% summarise_all(mean)
ggscatter(c, x="Dysfunction", y="IED172", add="loess",conf.int = TRUE,cor.coef.size = 6, 
          cor.coef = TRUE, cor.method = "pearson", size=0.1, xlab = "CD8-derived signature", 
          ylab = "IED172", add.params = list(color="steelblue", fill="skyblue"))
pdf("./Figures/Additional/IED172_correlation.pdf", width=4, height=4)
ggscatter(c2, x="Dysfunction", y="IED172", add="reg.line",conf.int = TRUE,cor.coef.size = 6, 
          cor.coef = TRUE, cor.method = "pearson", size=1, xlab = "CD8-derived signature", 
          ylab = "IED172", add.params = list(color="steelblue", fill="skyblue"))
dev.off()

c <- data.frame(cd8_eff$dysfunction, cd8_eff$IED68, cd8_eff$Patient, cd8_eff$class13)
colnames(c) <- c("Dysfunction", "IED68", "Patient", "class13")
c2 <- c %>% group_by(Patient,class13) %>% summarise_all(mean)
ggscatter(c, x="Dysfunction", y="IED68", add="loess",conf.int = TRUE,cor.coef.size = 6, 
          cor.coef = TRUE, cor.method = "pearson", size=0.1, xlab = "CD8-derived signature", 
          ylab = "IED68", add.params = list(color="steelblue", fill="skyblue"))
pdf("./Figures/Additional/IED68_correlation.pdf", width=4, height=4)
ggscatter(c2, x="Dysfunction", y="IED68", add="reg.line",conf.int = TRUE,cor.coef.size = 6, 
          cor.coef = TRUE, cor.method = "pearson", size=1, xlab = "CD8-derived signature", 
          ylab = "IED68", add.params = list(color="steelblue", fill="skyblue"))
dev.off()

pdf("./Plots/heatmap_patient_score.pdf", height=1.5, width=7)
Heatmap(t(dys), cluster_columns = F, width = 17*unit(7, "mm"), height = 1*unit(7, "mm"),
        rect_gp = gpar(col = "white", lwd = 1), col=colour_bk, column_gap = unit(7,"mm"),
        column_split = c(rep("115225",2), rep("220882",2),rep("452198",2), rep("508084",2), rep("823477",2), rep("869586",2)))
dev.off()

