############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for MiloR analysis
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#Exhaustion Expression
zhang_ex <- c("PDCD1", "TOX", "CXCL13", "TIGIT", "CTLA4", "TNFRSF9", "HAVCR2", "LAG3")

mat <- FetchData(cd8, zhang_ex)
mat <- groupMeans(t(mat),unlist(cd8$group))

pdf(paste0(file_prefix,"zhang_ex_group.pdf"), height=3.5, width=3)
ht <- Heatmap(t(mat), col = colour_bk, show_row_dend = F,rect_gp = gpar(col = "white", lwd = 1),
              width = nrow(mat)*unit(5, "mm"), height = ncol(mat)*unit(5, "mm"), row_names_gp =gpar(fontsize=8),
              column_names_gp =gpar(fontsize=8),heatmap_legend_param = list(title="Avg. Exp.", direction="horizontal"))
draw(ht, heatmap_legend_side ="bottom")
dev.off()

#Exhaustion Marker Expression (all clusters)
mat <- FetchData(cd8, zhang_ex)
mat <- groupMeans(t(mat),unlist(cd8$class13))

pdf("./Figures/Additional/zhang_ex_class13.pdf", height=5, width=4)
Heatmap(t(mat), col = colour_bk, show_row_dend = F,rect_gp = gpar(col = "white", lwd = 1),
        width = nrow(mat)*unit(5, "mm"), height = ncol(mat)*unit(5, "mm"), row_names_gp =gpar(fontsize=8),
        column_names_gp =gpar(fontsize=8),heatmap_legend_param = list(title="Avg. Exp."))
dev.off()

#Exhaustion Score vs. Pan-Cancer
panc_cd8 <- readRDS("~/OneDrive - Inside MD Anderson/Codes/Papers_and_Github/ZeminZhang_PanTCellScience2021/data/metaInfo/panC.freq.CD8.ht.tb.rds")

panc_cd8$exhausted <- "No"
panc_cd8$exhausted[grepl("Tex", panc_cd8$group.var)] <- "Yes"

#add aml to this
x <- table(cd8$Patient, cd8$class13)
y <- as.data.frame(x/rowSums(x))
colnames(y) <- c("Patient", "Cluster", "Freq")
y_ex <- subset(y, subset= Cluster %in% "Ex")
y_ex <- y_ex[!(y_ex$Patient %in% c("NLBM4", "NLBM5", "NLBM6")),]
colnames(y_ex) <- c("donor.var", "exhausted", "freq")
y_ex$exhausted <- NULL
y_ex$cancer <- "AML"
#order from lowest to highest

panc_ex <- subset(panc_cd8, subset = exhausted %in% "Yes")
panc_cd8_ex <- panc_ex %>% group_by(donor.var) %>% dplyr::summarise(freq=sum(freq), cancer=dplyr::first(cmp.var))
panc_df <- data.frame(panc_cd8_ex)
panc_aml_ex <- bind_rows(panc_df,y_ex)

mean_ex <- panc_aml_ex %>% group_by(cancer) %>% dplyr::summarise(mean = mean(freq))
med_ex <- panc_aml_ex %>% group_by(cancer) %>% dplyr::summarise(median = median(freq))
med_ex <- data.frame(med_ex)
med_ex <- med_ex[order(med_ex$median),]
med_ex_aml <- med_ex[c(2,1,3:21),]


pdf(paste0(file_prefix,"pancancer_exh.pdf"), height=3, width=4)
ggboxplot(data = panc_aml_ex, x="cancer", y="freq",fill="skyblue", outlier.shape = NA, remove=c("boxplot"),
          order=med_ex_aml$cancer, size=0.3) + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1), legend.position = "top", text = element_text(size=6)) + 
  stat_compare_means(label.x = 2, size=2, label.y=0.7) + guides(fill=guide_legend(title="Location"))
dev.off()

#Pan-Cancer Significance
kruskal.test(freq~cancer, panc_aml_ex)

#Pan Cancer exhaustion
pdf("./Figures/Figure1/pancancer_exh.pdf", height=3, width=4)
ggboxplot(data = panc_aml_ex, x="cancer", y="freq",fill="skyblue", outlier.shape = NA, remove=c("boxplot"),
          order=med_ex_aml$cancer, size=0.3) + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1), legend.position = "top", text = element_text(size=6)) + 
  stat_compare_means(label.x = 2, size=2, label.y=0.7) + guides(fill=guide_legend(title="Location"))
dev.off()

#Barplot for Exhaustion Relative Abundance
sample_levels <- c("NLBM4", "NLBM5", "NLBM6", "PT1A",   "PT2A",   "PT3A",   "PT4A",
                   "PT5A",   "PT6A",   "PT7A",   "PT8A",   "PT9A",  "PT10A",  "PT11A",  "PT12A",  "PT13A",  "PT14A",
                   "PT15A",  "PT16A",  "PT17A",  "PT18A",  "PT19A",  "PT20A",  "PT21A",  "PT22A",  "PT23A",  "PT24A",
                   "PT25A",  "PT26A",  "PT27A",  "PT28A",  "PT29A",  "PT30A", "PT32A")
sample_type <- c("Healthy", "Healthy", "Healthy", "RelRef",   "RelRef",   "RelRef",   "RelRef",
                 "RelRef",  "RelRef",   "RelRef",  "RelRef",   "NewlyDx",  "NewlyDx", "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",
                 "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",
                 "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx", "NewlyDx")
names(sample_type) <- sample_levels

x <- table(cd8$Patient, cd8$class13)
y <- as.data.frame(x/rowSums(x))
colnames(y) <- c("Patient", "Cluster", "Freq")
#Add type information
y$Type <- sample_type[match(y$Patient, patient_levels)]
y$Type <- factor(y$Type, levels=c("Healthy","NewlyDx", "RelRef"))
y$group[y$Cluster %in% "Ex"] <-"Ex"

y_group <- y[y$group =="Ex",]
y_group$group <- droplevels(y_group$group)
y_group$Cluster <- droplevels(y_group$Cluster)
wilcox_stat_cd8 <- y_group %>% group_by(Cluster) %>% pairwise_wilcox_test(Freq~Type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="Cluster", fun="max", step.increase = 0.1)
wilcox_sig_cd8 <- subset(wilcox_stat_cd8, p.adj < 0.05)
ex_box <- ggboxplot(data=y_group,x="Cluster", y="Freq", color="Type", palette = anno_colors$type,bxp.errorbar = TRUE)  
ex_box <- ex_box +stat_pvalue_manual(wilcox_stat_cd8, label="p.adj", tip.length = 0.0)

pdf(paste0(file_prefix,"Exhaustion_box.pdf"))
ex_box
dev.off()

# Corrplot exhaustion markers
ex <- FetchData(cd8, c(zhang_ex,"class13"))
ex <- subset(ex, ex$class13 == "Ex")
ex$class13 <- NULL
C <- cor(ex)
testRes = cor.mtest(ex, conf.level = 0.95)

pdf(paste0(file_prefix,"corr_ex_circles.pdf", height=4, width=4))
corrplot(C, tl.cex = 0.5, tl.col="black", col=COL2('BrBG'), tl.pos = 'd', 
         type='lower', method = 'circle')
dev.off()

no_ex <- FetchData(cd8, c(zhang_ex,"class13"))
no_ex <- subset(no_ex, no_ex$class13 != "Ex")
no_ex$class13 <- NULL
c <- cor(no_ex)
testRes1 = cor.mtest(ex, conf.level = 0.95)

pdf(paste0(file_prefix,"corr_no_ex_circles.pdf", height=4, width=4))
corrplot(c, tl.cex = 0.5, tl.col="black", col= COL2('BrBG'), tl.pos = 'd',
         type='upper')
dev.off()

#Immune Inhibitory Marker Expression
inhibitors <- c("ADORA2A", "ARG1", "BTLA", "BTN2A2", "CD160", "CD200", "CD200R1", "CD274", "CD276",
                "CD44", "CD47", "CD96", "CEACAM1", "CTLA4", "CYBB", "HAVCR2", "HLA-DPA1", "HLA-G", "ICOSLG",
                "IDO1", "IL2RB", "KIR2DL1", "KIR2DL3", "KIR2DL4",  "KIR3DL1", "KIR3DL2", "KIR3DL3", 
                "KLRC1", "LAG3", "LAIR1", "LGALS3", "LGALS9", "LILRB1", "LILRB2", "LYN", "NECTIN2", "NRP1", "PDCD1",
                "PDCD1LG2", "PTPN1", "PVRIG", "REL", "SIGLEC7", "SIGLEC9", "SIRPA", "TDO2", "TIGIT", "TNFRSF14",
                "TNFRSF8", "VSIR")

mat <- cd8@assays$RNA@data[inhibitors,]
mat <- groupMeans(mat,unlist(cd8[["group_type"]]))
mat <- mat[,c("N Healthy", "N NewlyDx", "N RelRef", "Mem Healthy", "Mem NewlyDx", "Mem RelRef",
              "EM Healthy", "EM NewlyDx", "EM RelRef", "CTL Healthy", "CTL NewlyDx", "CTL RelRef", "MAIT Healthy", "MAIT NewlyDx", 
              "MAIT RelRef", "Ex Healthy", "Ex NewlyDx", "Ex RelRef")]
ht1 <- Heatmap(mat, cluster_rows = TRUE, cluster_columns = FALSE, heatmap_legend_param = list(title="Expression"),
               column_split = c(rep("A",3), rep("B",3),rep("C",3), rep("D",3), rep("E",3), rep("F",3)), cluster_column_slices = FALSE,
               column_title = NULL,  width = ncol(mat)*unit(5, "mm"), height = nrow(mat)*unit(5, "mm"))

pdf(paste0(file_prefix,"inh_type_ht_exp.pdf", width=7, height=12))
draw(ht1)
dev.off()


