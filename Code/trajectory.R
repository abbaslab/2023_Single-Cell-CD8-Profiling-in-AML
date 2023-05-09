############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Trajectory Analysis
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#Monocle3 Trajectory Analysis
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


#Monocle 3 Group/Type
pdf(paste0(file_prefix,"cd8_monocle_type.pdf"), height=6, width=8)
plot_cells(cd8.cds, color_cells_by = "type", label_branch_points = FALSE, label_leaves=FALSE, label_cell_groups = FALSE) + 
  scale_color_manual(values=anno_colors$type, name="Patient Type") & NoAxes()
dev.off()

pdf(paste0(file_prefix,"cd8_monocle_group.pdf"), height=6, width=8)
plot_cells(cd8.cds, color_cells_by = "group", label_branch_points = FALSE, label_leaves=FALSE, label_cell_groups = FALSE) + 
  scale_color_manual(values=anno_colors$group, name="Cell Type") & NoAxes()
dev.off()


#Monocle Reduction Feature Plots
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


#Type Histogram over Pseudotime
pdf(paste0(file_prefix,"pseudotime_histogram.pdf"))
ggplot(cd8@meta.data, aes(x=pseudotime, fill=type)) + geom_histogram(alpha=0.4, bins=20) + theme_classic() +
  scale_fill_manual(values=anno_colors$type, name="Patient Type") +
  labs(x="Pseudotime", y="Frequency")
dev.off()

#Faceted Type Histogram over Pseudotime
pdf("./Figures/Additional/pseudotime_histogram_faceted.pdf", height=3,width=7)
ggplot(cd8@meta.data, aes(x=pseudotime, fill=type)) + geom_histogram(alpha=0.6, bins=20,color="white") + theme_classic() +
  scale_fill_manual(values=anno_colors$type, name="Patient Type") +
  labs(x="Pseudotime", y="Frequency") + facet_wrap(.~type)
dev.off()

#Genes over Pseudotime
cd8_cds_pr_test_res <- graph_test(cd8.cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(cd8_cds_pr_test_res, q_value == 0 & morans_I > 0.25))
genes <- pr_deg_ids[!startsWith(pr_deg_ids, "RP")]
genes <- c(genes, "KLRG1")

pt.matrix <- normalized_counts(cd8.cds, norm_method = "log")[match(genes,rownames(rowData(cd8.cds))),order(pseudotime(cd8.cds))]
png("./Figures/Figure3/pseudotime_heatmap.png", res=200, bg="white", units="in", width=4, height=6)
Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 4),
  clustering_method_rows = "average",
  clustering_distance_rows = "pearson",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE, column_title = "Pseudotime", column_title_side = "bottom")
dev.off()

#Cluster Vln Plot across Pseudotime
VlnPlot(cd8, "pseudotime", pt.size=0, group.by="class13", cols=anno_colors$class13) + coord_flip() & NoLegend()
ggsave("./Figures/Figure3/class13_pseudotime.png", bg="white")


#correlation between slingshot and monocle3
library(slingshot)
sce <- slingshot(as.SingleCellExperiment(cd8),clusterLabels = "cluster_group",reducedDim = 'PCA')
ggplot(as.data.frame(colData(as.SingleCellExperiment(cd8))), aes(x = sce$slingPseudotime_1, y = cluster_group, colour = cluster_group)) +
  geom_quasirandom(groupOnX = FALSE) + scale_color_tableau() + theme_classic() +
  xlab("Pseudotime") + ylab("Clusters") + ggtitle("")

pseudotime=data.frame(slingshot_pseudotime=sce$slingPseudotime_1,slingshot2=sce$slingPseudotime_2,slingshot3=sce$slingPseudotime_3,type=sce$type,cluster_group=sce$cluster_group,patient=sce$Patient)
pseudotime$slingshot_pseodutime=ifelse(is.na(pseudotime$slingshot_pseodutime),pseudotime$slingshot2, pseudotime$slingshot_pseodutime)
pseudotime$id=colnames(sce)

temp=data.frame(monocle=cd8$pseudotime,id=colnames(cd8))
pseudotime=merge(pseudotime,temp,by='id')
ggscatter(pseudotime, x = "slingshot_pseodutime", y = "monocle", add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",xlab = "slingshot", ylab = "monocle",size=0.4) +
  theme(legend.position = "none", axis.title = element_text(size = 26),
        axis.text.x = element_text(size = 23),axis.text.y = element_text(size = 23),
        axis.ticks = element_blank()) 
rm(temp,pseudotime)