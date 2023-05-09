############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for cluster figure generation
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#Group UMAP
pdf(paste0(file_prefix,"group_umap.pdf"), width=4,height=4)
DimPlot(cd8, group.by="group", cols=anno_colors$group, label=TRUE, label.box = TRUE, label.size = 3) + ggtitle(NULL) & NoLegend() & NoAxes()
dev.off()

#All Clusters UMAP
pdf(paste0(file_prefix,"class13_umap.pdf"), width=3, height=3)
DimPlot(cd8, group.by="class13", cols=anno_colors$class13, label=FALSE) + ggtitle(NULL) & NoLegend() & NoAxes()
dev.off()

#Group Dot Plot
can_markers <- c("CCR7", "LEF1", "SELL", "IL7R", "LTB", "GZMB", "GNLY", "PRF1", "NKG7", "GZMK", 
                 "PDCD1", "LAG3", "TIGIT", "SLC4A10", "TRAV1-2")

pdf(paste0(file_prefix,"canonical_group_dot_plot.pdf"), height=4, width=2.5)
dot_plot_x(cd8, rev(can_markers), group = "group")
dev.off()

#Group Distribution by Patient Type
group_type_df <- as.data.frame(table(cd8$group, cd8$type))
kruskal.test(Freq~Var1, data=group_type_df)

pdf(paste0(file_prefix,"group_stacked_bar.pdf"), width=3, height=4)
ggplot(cd8@meta.data, aes(x=type, fill=group, color=group)) + geom_bar(position = position_fill(reverse = TRUE), alpha=0.8, width=0.5) + 
  theme_classic() + scale_fill_manual(values=group.colors, name= "Cell Types") + scale_color_manual(values=anno_colors$group, name= "Cell Types") +
  theme( axis.text.x = element_text(size=10, angle=90), plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  labs(x=NULL, y="Percent")
dev.off()


#Pie Chart
group <- data.frame(table(cd8$group))
colnames(group) <- c("cluster", "n")
group$cluster <- factor(group$cluster)
group$prop<-round(group$n/sum(group$n)*100,2)
group <- group %>% arrange(cluster,dplyr::desc(n)) %>% mutate(lab.ypos=100-(cumsum(prop)-0.5*prop))

pdf(paste0(file_prefix,"group_pie_chart.pdf"), height=3, width=3)
ggplot(group, aes(x = "", y = prop, fill = cluster)) + geom_bar(width = 1, stat = "identity", color = "white") + 
  scale_fill_manual(values = anno_colors$group)  + geom_label_repel(aes(label=prop, y=lab.ypos), size=3, show.legend = F) + 
  coord_polar("y")  + theme_void() + theme(legend.position = "none")
dev.off()

#Groups by Patient Pie Chart
a <- data.frame(table(cd8$group, cd8$Patient))
pdf(paste0(file_prefix,"cd8_patient_group_pie.pdf", width=9, height=9))
ggplot(data=a, aes(x=1, y=Freq,fill=Var1, color=Var1)) + 
  scale_fill_manual(values=anno_colors$group) + 
  scale_color_manual(values = anno_colors$group) +
  geom_col(position='fill', color="white") + coord_polar(theta="y", start = 0) +
  facet_wrap(.~ Var2) + 
  theme_void()
dev.off()

#Cluster Membership by Patient Type/Group
pdf(paste0(file_prefix,"group_cluster_membership.pdf", width=3, height=3))
cluster_membership(cd8$type, cd8$group)
dev.off()

#Cluster Membership by Patient Type/Cluster
pdf(paste0(file_prefix,"type_cluster_membership.pdf", width=3, height=3))
cluster_membership(cd8$type, cd8$class13)
dev.off()

#Dot Plot All Clusters
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


##Density Plot
gdat<-Embeddings(cd8, reduction = 'umap')
gdat<-as.data.frame(gdat)
gdat$type <- cd8$type

pdf(paste0(file_prefix,"density_plot_vertical.pdf"), width=4, height=10)
ggplot(data = gdat,aes(x = UMAP_1, y = UMAP_2)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F) + labs(x="UMAP1", y="UMAP2") +
  geom_point(color = 'white',size = .005) + scale_fill_viridis(option="magma") + theme_black() + 
  facet_wrap(.~type,ncol = 1) 
dev.off()


#Box Plot

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

#KIR gene expression
pdf(paste0(file_prefix,"kir_expression.pdf"), width=4, height=10)
dot_plot_x(cd8, c("KIR2DL1", "KIR2DL3", "KIR2DL4", "KIR3DL1", "KIR3DL2", "KIR3DL3",), group = "class13")
dev.off()

#Regulatory Gene expression
pdf(paste0(file_prefix,"reg_expression.pdf"), width=4, height=10)
dot_plot_x(cd8, c("FOXP3", "ITGAE", "CD28", "IL2RB", "ITGA4", "IKZF2", "KIR3DL1", "KIR2DL3"), group = "class13")
dev.off()

#Downsampling
cd8_newdx <- subset(cd8, subset = type %in% "NewlyDx")
cd8_hr <- subset(cd8, subset = type %in% c("Healthy", "RelRef"))
newdx_sample <- sample(colnames(cd8_newdx), size=4397)

cd8_sub <- subset(cd8_newdx, cells = newdx_sample)
cd8_sub <- merge(cd8_sub, cd8_hr)
cd8_sub <- analyze_seurat(cd8_sub)

pdf("./Figures/Additional/cd8_sample_umap_type.pdf", height=3, width=7)
DimPlot(cd8_sub, group.by = "group", cols = anno_colors$group, split.by = "type") 
dev.off()

DimPlot(cd8_sub, group.by = "RNA_snn_res.1", label=T)
dot_plot_x(cd8_sub, canonical_marker, group="RNA_snn_res.1")
dot_plot_x(cd8_sub, feat_inf, group="RNA_snn_res.1")

feat_naive <- c("LTB","IL7R", "SELL", "LEF1", "TCF7", "CCR7")
feat_ctl <- c("NKG7", "PRF1", "GZMH", "GZMB", "GNLY")
feat_mait <- c("SLC4A10", "TRAV1-2")
feat_cd8 <- c("CX3CR1", "KLRG1", "FOSB", "CXCR4", "CD44", "EOMES", "GZMK")
feat_cyt <- c("STAT1", "IFNG", "IFIT3", "TNF")
feat_ex <- c("PDCD1", "TIGIT", "TOX","LAG3")
cd8_dot <- group_dot_plot(cd8_sub, features=list(feat_naive, feat_ctl, feat_cd8,feat_cyt, feat_mait, feat_ex), 
                          group = "RNA_snn_res.1")
pdf("./Figures/Additional/cd8_subsample_dot_plot.pdf", height=8.25, width=4.25)
ggarrange(cd8_dot)
dev.off()


#Correlation Plots
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

