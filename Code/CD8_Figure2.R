############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Figure 2
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: 
############################################

file_prefix <- ""
source(paste0(file_prefix,"CD8_Functions.R"))
cd8 <- readRDS("cd8.rds")
library(miloR); library(patchwork)

########## A: Dot Plot All Clusters ##########
feat_naive <- c("LTB","IL7R", "SELL", "LEF1", "TCF7", "CCR7")
feat_ctl <- c("NKG7", "PRF1", "GZMH", "GZMB", "GNLY")
feat_mait <- c("SLC4A10", "TRAV1-2")
feat_cd8 <- c("CX3CR1", "KLRG1", "FOSB", "CXCR4", "CD44", "EOMES", "GZMK")
feat_cyt <- c("STAT1", "IFNG", "IFIT3", "TNF")
feat_ex <- c("PDCD1", "TIGIT", "TOX","LAG3")
cd8_dot <- group_dot_plot(cd8, features=list(feat_naive, feat_ctl, feat_cd8,feat_cyt, feat_mait, feat_ex), 
                          group = "class13",feats=c("Naive/Memory", "Cytotoxic", "Differential", "Inflammatory", "MAIT", "Exhaustion"),
                          lim = c("N", "Mem1", "Mem2", "EM1", "EM2", "CTL1", "CTL2", "NK-like1", "NK-like2", "INF1", "INF2", "CTL INF1",
                                  "CTL INF2", "CTL INF3", "MAIT", "Ex"))
leg <- get_legend(dot_plot_leg(cd8, feat_naive))
pdf(paste0(file_prefix,"cd8_dot_plot.pdf"), height=8.25, width=4.25)
ggarrange(cd8_dot)
dev.off()

pdf(paste0(file_prefix,"cd8_dot_legend.pdf"), height=4, width=2)
ggarrange(leg)
dev.off()


########## B: Density Plot #########
gdat<-Embeddings(cd8, reduction = 'umap')
gdat<-as.data.frame(gdat)
gdat$type <- cd8$type

pdf(paste0(file_prefix,"density_plot_vertical.pdf"), width=4, height=10)
ggplot(data = gdat,aes(x = UMAP_1, y = UMAP_2)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F) + labs(x="UMAP1", y="UMAP2") +
  geom_point(color = 'white',size = .005) + scale_fill_viridis(option="magma") + theme_black() + 
  facet_wrap(.~type,ncol = 1) 
dev.off()


########## C: Box Plot ##########

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

y$group[y$Cluster %in% "N"] <-"N"
y$group[y$Cluster %in% c("Mem1","Mem2")] <- "Mem"
y$group[y$Cluster %in% "MAIT"] <-"MAIT"
y$group[y$Cluster %in% "Ex"] <-"Ex"
y$group[y$Cluster %in% c("EM1", "EM2", "INF1", "INF2")] <-"EM"
y$group[y$Cluster %in% c("CTL1", "CTL2", "NK-like1", "NK-like2")] <-"CTL"
y$group[y$Cluster %in% c("CTL INF1", "CTL INF2", "CTL INF3")] <- "CTL INF"
y$group <- factor(y$group, levels=c("N", "Mem", "MAIT", "Ex", "EM", "CTL", "CTL INF"))

y_group <- y[y$group =="N",]
y_group$group <- droplevels(y_group$group)
y_group$Cluster <- droplevels(y_group$Cluster)
wilcox_stat_cd8 <- y_group %>% group_by(Cluster) %>% pairwise_wilcox_test(Freq~Type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="Cluster", fun="max", step.increase = 0.1)
wilcox_sig_cd8 <- subset(wilcox_stat_cd8, p.adj < 0.05)
naive_box <- ggboxplot(data=y_group,x="Cluster", y="Freq", color="Type", palette = anno_colors$type,bxp.errorbar = TRUE)  
naive_box <- naive_box +stat_pvalue_manual(wilcox_sig_cd8, label="p.adj.signif", tip.length = 0.0)

y_group <- y[y$group =="Mem",]
y_group$group <- droplevels(y_group$group)
y_group$Cluster <- droplevels(y_group$Cluster)
wilcox_stat_cd8 <- y_group %>% group_by(Cluster) %>% pairwise_wilcox_test(Freq~Type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="Cluster", fun="max", step.increase = 0.3)
wilcox_sig_cd8 <- subset(wilcox_stat_cd8, p.adj < 0.05)
memory_box <- ggboxplot(data=y_group,x="Cluster", y="Freq", color="Type", palette = anno_colors$type,bxp.errorbar = TRUE)  
memory_box <- memory_box +stat_pvalue_manual(wilcox_sig_cd8, label="p.adj.signif", tip.length = 0.0)

y_group <- y[y$group =="EM",]
y_group$group <- droplevels(y_group$group)
y_group$Cluster <- droplevels(y_group$Cluster)
wilcox_stat_cd8 <- y_group %>% group_by(Cluster) %>% pairwise_wilcox_test(Freq~Type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="Cluster", fun="median_iqr", step.increase = 0.2)
wilcox_sig_cd8 <- subset(wilcox_stat_cd8, p.adj < 0.05)
em_box <- ggboxplot(data=y_group,x="Cluster", y="Freq", color="Type", palette = anno_colors$type,bxp.errorbar = TRUE)  
em_box <- em_box +stat_pvalue_manual(wilcox_sig_cd8, label="p.adj.signif", tip.length = 0.0)

y_group <- y[y$group =="CTL",]
y_group$group <- droplevels(y_group$group)
y_group$Cluster <- droplevels(y_group$Cluster)
wilcox_stat_cd8 <- y_group %>% group_by(Cluster) %>% pairwise_wilcox_test(Freq~Type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="Cluster", fun="median_iqr", step.increase = 0.15)
wilcox_sig_cd8 <- subset(wilcox_stat_cd8, p.adj < 0.05)
ctl_box <- ggboxplot(data=y_group,x="Cluster", y="Freq", color="Type", palette = anno_colors$type,bxp.errorbar = TRUE)  
ctl_box <- ctl_box +stat_pvalue_manual(wilcox_sig_cd8, label="p.adj.signif", tip.length = 0.0)

y_group <- y[y$group =="CTL INF",]
y_group$group <- droplevels(y_group$group)
y_group$Cluster <- droplevels(y_group$Cluster)
wilcox_stat_cd8 <- y_group %>% group_by(Cluster) %>% pairwise_wilcox_test(Freq~Type, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% add_xy_position(x="Cluster", fun="median_iqr", step.increase = 1)
wilcox_sig_cd8 <- subset(wilcox_stat_cd8, p.adj < 0.05)
ctl_inf_box <- ggboxplot(data=y_group,x="Cluster", y="Freq", color="Type", palette = anno_colors$type,bxp.errorbar = TRUE)  
ctl_inf_box <- ctl_inf_box +stat_pvalue_manual(wilcox_sig_cd8, label="p.adj.signif", tip.length = 0.0)


cd8_box <- (naive_box + memory_box + plot_layout(ncol=2, widths=c(1,2))) / em_box / ctl_box / ctl_inf_box
cd8_box <- cd8_box & xlab(NULL) & theme(axis.text.x = element_text(size=8, angle=60, hjust=1),
                                                     axis.text.y = element_text(size=8))
cd8_box <- cd8_box + plot_layout(guides="collect") & theme(legend.position = "top")
pdf(paste0(file_prefix,"cd8_box_bh_scaled.pdf"), width=4, height=10)
cd8_box & ylim(0,0.7)
dev.off()

########## D: Milo UMAP ##########

cd8_milo <- as.SingleCellExperiment(cd8)
cd8_milo <- Milo(cd8_milo)
cd8_milo <- buildGraph(cd8_milo, k = 20, d = 50)
cd8_milo <- makeNhoods(cd8_milo, prop = 0.1, k = 20, d=50, refined = TRUE)

#As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples.
#plotNhoodSizeHist(cd8_milo)
cd8_milo <- countCells(cd8_milo, meta.data = data.frame(colData(cd8_milo)), sample="orig.ident")
cd8_milo <- calcNhoodDistance(cd8_milo, d=50)

cd8_design <- data.frame(colData(cd8_milo))[,c("orig.ident", "type")]
cd8_design$type <- as.factor(cd8_design$type)
cd8_design <- distinct(cd8_design)
rownames(cd8_design) <- c(cd8_design$orig.ident)

da_results <- testNhoods(cd8_milo, design = ~ type, design.df = cd8_design)
da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "class13")
da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "type")
#ggplot(da_results, aes(type_fraction)) + geom_histogram(bins=50)
da_results$type <- ifelse(da_results$type_fraction < 0.7, "Mixed", da_results$type)
da_results$type <- factor(da_results$type, rev(c("Healthy", "NewlyDx", "Mixed", "RelRef")))
da_results$class13 <- factor(da_results$class13, rev(c("N", "Mem1", "Mem2", "EM1", "EM2","INF1", "INF2",
                                                       "CTL1", "CTL2", "CTL INF1", "CTL INF2", "CTL INF3", 
                                                       "NK-like1", "NK-like2", "MAIT", "Ex")))

cd8_milo <- buildNhoodGraph(cd8_milo)
milo <- plotNhoodGraphDA(cd8_milo, da_results, alpha=1,size_range=c(1,5)) + guides(edge_width="none")
milo <- milo + guides(size=guide_legend(title="Neighborhood Size"))

pdf(paste0(file_prefix,"milo_umap.pdf"), height=5, width=7)
milo
dev.off()


########## E: Milo BeeSwarm ##########
pdf(paste0(file_prefix,"milo_type_DAbeeswarm.pdf"), width = 4, height=1.5)
plotDAbeeswarm(da_results, group.by = "type") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf(paste0(file_prefix,"milo_class13_beeswarm.pdf"), width = 4, height=4.5)
plotDAbeeswarm(da_results, "class13") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

