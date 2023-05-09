############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for GSEA
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: file_prefix should be the input/output directory
############################################

library(tidyverse); library(ggrepel)
library(ggVennDiagram); library(fgsea)
library(enrichR); library(stringr); library(msigdbr)
library(patchwork); library(ComplexHeatmap)

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(file_prefix,"cd8.rds")

#Preliminary Analysis
Idents(cd8) <- cd8$class13
#Differences in CTL, NKT, EM
ctl_markers <- FindMarkers(cd8, ident.1 = "CTL1", ident.2 = "CTL2", logfc.threshold = 0)
ctl_frame <- ctl_markers %>% mutate(Expression = case_when(avg_log2FC >= log(2) & p_val_adj <= 0.05 ~ "Upregulated",
                                                           avg_log2FC <= -log(2) & p_val_adj <= 0.05 ~ "Downregulated",
                                                           TRUE ~ "Unchanged"))
pval_0 <- subset(ctl_frame, ctl_frame$p_val_adj == 0)
pval_jitter <- round(rnorm(nrow(pval_0), mean =295, sd = 5),0) 
pval_jitter <- 1*10^-pval_jitter
pval_0$pvalue <- pval_jitter
ctl_frame[rownames(pval_0),]$p_val_adj <- pval_0$pvalue


nkt_markers <- FindMarkers(cd8, ident.1 = "NK-like1", ident.2 = "NK-like2", logfc.threshold = 0)
nkt_frame <- nkt_markers %>% mutate(Expression = case_when(avg_log2FC >= log(2) & p_val_adj <= 0.05 ~ "Upregulated",
                                                           avg_log2FC <= -log(2) & p_val_adj <= 0.05 ~ "Downregulated",
                                                           TRUE ~ "Unchanged"))
pval_0 <- subset(nkt_frame, nkt_frame$p_val_adj == 0)
pval_jitter <- round(rnorm(nrow(pval_0), mean =295, sd = 5),0) 
pval_jitter <- 1*10^-pval_jitter
pval_0$pvalue <- pval_jitter
nkt_frame[rownames(pval_0),]$p_val_adj <- pval_0$pvalue


em_markers <- FindMarkers(cd8, ident.1 = "EM1", ident.2 = "EM2", logfc.threshold = 0)
em_frame <- em_markers %>% mutate(Expression = case_when(avg_log2FC >= log(2) & p_val_adj <= 0.05 ~ "Upregulated",
                                                         avg_log2FC <= -log(2) & p_val_adj <= 0.05 ~ "Downregulated",
                                                         TRUE ~ "Unchanged"))
pval_0 <- subset(em_frame, em_frame$p_val_adj == 0)
pval_jitter <- round(rnorm(nrow(pval_0), mean =295, sd = 5),0) 
pval_jitter <- 1*10^-pval_jitter
pval_0$pvalue <- pval_jitter
em_frame[rownames(pval_0),]$p_val_adj <- pval_0$pvalue


#Pathway Plots
markers_all <- FindMarkers(cd8, ident.1 = c("EM1", "CTL1", "NK-like1"), ident.2 = c("EM2", "CTL2", "NK-like2"),
                           logfc.threshold = 0)

hallmark_gene_sets <- msigdbr(species="human", category="H")
pathwaysH = split(x = hallmark_gene_sets$human_gene_symbol, f = hallmark_gene_sets$gs_name)

go_bp <- msigdbr(species="human", category="C5", subcategory = c("GO:BP"))
pathwaysGO_BP = split(x = go_bp$human_gene_symbol, f = go_bp$gs_name)

ranks_for_fgsea <- function(find_markers_output){
  all.markers <- find_markers_output[order(-find_markers_output$avg_log2FC),]
  all.markers.list <- c(all.markers$avg_log2FC)
  names(all.markers.list) <- rownames(all.markers)
  return(all.markers.list)
}

all_ranks <- ranks_for_fgsea(markers_all)
all_fgsea <- fgsea(pathways=pathwaysH, stats=all_ranks,minSize=15)
all_fgsea <- data.frame(all_fgsea)

go_bp_fgsea <- fgsea(pathways=pathwaysGO_BP, stats=all_ranks,minSize=15)
go_bp_fgsea <- data.frame(go_bp_fgsea)
go_bp_fgsea_sig <- go_bp_fgsea[go_bp_fgsea$padj<0.05,]
go_bp_fgsea_sig <- go_bp_fgsea_sig[abs(go_bp_fgsea_sig$NES)>1.5,]

go_bp_pathways <- c("GOBP_CELLULAR_RESPONSE_TO_STRESS","GOBP_REGULATION_OF_MRNA_METABOLIC_PROCESS","GOBP_CELL_CYCLE", "GOBP_REGULATION_OF_CELL_DIFFERENTIATION","GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS")

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space="Lab") ##color scheme

hallmark_paths <- bind_rows(all_fgsea %>% filter(NES >0) %>% arrange(padj) %>% head(4),
                            all_fgsea %>% filter(NES<0) %>% arrange(padj) %>% head(5), all_fgsea[17,]
) %>% mutate(pathway = str_remove(pathway, "HALLMARK_")) %>% mutate(pathway=str_replace_all(pathway, "_", " ")) %>%
  arrange(NES) %>% mutate(pathway=factor(pathway,levels=pathway)) %>% 
  ggplot(aes(x=pathway, y=-log10(padj))) + labs(x="Hallmark", y=NULL) +
  geom_point(aes(y = -log10(padj),x = pathway, color = 'black',fill =NES,size = abs(NES)),shape = 21,
             color = 'black',stroke = .2) + ylim(0,28) + scale_size(name="NES",breaks=c(1.5,1.75,2,2.25,2.5), limits=c(1,3))+
  scale_fill_gradientn(colours =myPalette(100), breaks=c(-2,-1,0,1,2), limits=c(-3,2.5)) + coord_flip() + theme_bw()

go_bp_paths <- bind_rows(go_bp_fgsea %>% filter(NES >0) %>% arrange(padj) %>% head(5),
                         go_bp_fgsea[go_bp_fgsea$pathway %in% go_bp_pathways,] %>% arrange(padj),
) %>% mutate(pathway = str_remove(pathway, "GOBP_")) %>% mutate(pathway=str_replace_all(pathway, "_", " ")) %>%
  arrange(NES) %>% mutate(pathway=factor(pathway,levels=pathway)) %>% 
  ggplot(aes(x=pathway, y=-log10(padj))) + labs(x="GO Biological Process", y="-log10(FDR)") +
  geom_point(aes(y = -log10(padj),x = pathway, color = 'black',fill =NES, size = abs(NES)),shape = 21,
             color = 'black',stroke = .2) + ylim(0,28) + scale_size(name="NES", breaks=c(1.5,1.75,2,2.25,2.5), limits=c(1,3)) + 
  scale_fill_gradientn(colours =myPalette(100), breaks=c(-2,-1,0,1,2), limits=c(-3,2.5)) + coord_flip() + theme_bw()


pdf(paste0(file_prefix,"pathway_plots.pdf", height=4.5, width=7))
hallmark_paths / go_bp_paths + plot_layout(guides = "collect")
dev.off()

#Bubble Plots

all_top_pathways <- bind_rows(all_fgsea %>% filter(NES >0) %>% arrange(padj, NES) %>% head(4),all_fgsea[17,],
                              all_fgsea %>% filter(NES<0) %>% arrange(padj, NES) %>% head(top)
) %>% mutate(pathway_name = str_remove(pathway, "HALLMARK_")) %>% mutate(pathway_name=str_replace_all(pathway_name, "_", " "))



pdf("./Figures/Supplementary/Supp4_hallmark_fgsea.pdf")
ggplot(all_fgsea, aes(y = -log10(padj),x = NES,label=pathway)) +
  geom_point(aes(y = -log10(padj),x = NES, color = 'black',fill =NES,size = padj),shape = 21,
             color = 'black',stroke = .2) + geom_label_repel(data = all_top_pathways, mapping = aes(NES, -log(padj,10), label = pathway_name),
                                                             size = 2, max.overlaps = 10, min.segment.length=0, force=5, xlim = c(NA, Inf),
                                                             force_pull=0.25) + scale_fill_gradientn(colours =myPalette(100)) +
  scale_size(range = c(5,1),limits = c(0, 1), breaks = c(0,0.001,0.01,0.05,1)) + theme_classic() + coord_cartesian(clip = "off") +
  xlab(expression("log"[2]*"NES")) + ylab(expression("-log"[10]*"FDR")) + ggtitle("FGSEA of CD8 Clusters for Hallmark Pathways") +
  theme(legend.justification = "top")
dev.off()

gobp_top_pathways <- bind_rows(go_bp_fgsea_sig %>% filter(NES >0) %>% arrange(padj, NES) %>% head(top),
                               go_bp_fgsea_sig %>% filter(NES<0) %>% arrange(padj, NES) %>% head(top)
) %>% mutate(pathway_name = str_remove(pathway, "HALLMARK_")) %>% mutate(pathway_name=str_replace_all(pathway_name, "_", " "))


pdf("./Figures/Supplementary/Supp4_gobp_fgsea.pdf")
ggplot(go_bp_fgsea, aes(y = -log10(padj),x = NES,label=pathway)) +
  geom_point(aes(y = -log10(padj),x = NES, color = 'black',fill =NES,size = padj),shape = 21,
             color = 'black',stroke = .2) + geom_label_repel(data = gobp_top_pathways, mapping = aes(NES, -log(padj,10), label = pathway_name),
                                                             size = 2, max.overlaps = 10, min.segment.length=0, force=5, xlim = c(NA, Inf),
                                                             force_pull=0.25) + scale_fill_gradientn(colours =myPalette(100)) +
  scale_size(range = c(5,1),limits = c(0, 1), breaks = c(0,0.001,0.01,0.05,1)) + theme_classic() + coord_cartesian(clip = "off") +
  xlab(expression("log"[2]*"NES")) + ylab(expression("-log"[10]*"FDR")) + ggtitle("FGSEA of CD8 Top Clusters for GO BP Pathways") +
  theme(legend.justification = "top")
dev.off()


#Cytosig
cytosig_class9 <- data.frame(fread("cytosig_class9.csv"))
rownames(cytosig_class9) <- c("CTL INF2", "CTL INF3", "NK-like1", "EM1", 
                              "CTL1", "N", "Mem1", "INF1", "NK-like2", "CTL2", 
                              "EM2", "Mem2", "MAIT", "CTL INF1", "INF2", "Ex")
cytosig_class9$ID <- NULL

cytosig <- cytosig_class9[rev(c("EM1", "CTL1", "NK-like1", "EM2", "CTL2", "NK-like2")),]

pdf(paste0(file_prefix,"cytosig_heatmap.pdf", height=2, width=9.5))
Heatmap(cytosig, show_column_dend = FALSE, show_row_dend = FALSE, rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(title="Activity"))
dev.off()

#Correlation Plots
expression_matrix <- FetchData(cd8, vars = c("GOBP_AEROBIC_RESPIRATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                                             "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "class13", "type"))
df <- expression_matrix %>% group_by(class13,type,.drop = FALSE) %>% summarise(AEROBIC=mean(GOBP_AEROBIC_RESPIRATION),OXPHOS=mean(HALLMARK_OXIDATIVE_PHOSPHORYLATION),
                                                                               HYPOXIA = mean(HALLMARK_HYPOXIA), GLYCOLYSIS=mean(HALLMARK_GLYCOLYSIS))


hyp_gly <- ggplot(df, aes(y=GLYCOLYSIS, x=HYPOXIA,shape=type)) + geom_point(aes(size=2,color=class13))+
  geom_smooth(se=FALSE,method = lm, linetype="dashed", size=0.5) + scale_color_manual(values=class13.colors, name="Clusters")+ 
  theme_classic() + stat_cor() + theme(axis.title = element_text(size=12))

aero_op <- ggplot(df, aes(y=OXPHOS, x=AEROBIC,shape=type)) + geom_point(aes(size=2,color=class13))+
  geom_smooth(se=FALSE,method = lm, linetype="dashed", size=0.5) + scale_color_manual(values=class13.colors, name="Clusters")+ 
  theme_classic() + stat_cor() + xlim(0.15,0.23) + ylim(0.1,0.2) + theme(axis.title = element_text(size=12))


pdf(paste0(file_prefix,"pathway_correlation_plots.pdf", height=7, width=13))
hyp_gly + aero_op + plot_layout(guides="collect")
dev.off()

#Load differential expression markers

#Volcano Plots
p1 <- ggplot(ctl_frame, aes(avg_log2FC, -log(p_val_adj,10))) +  geom_point(size = 1) + xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"FDR"))
p2 <- p1 + geom_point(aes(color = Expression), size = 1) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
  geom_vline(xintercept=c(-log(2),log(2)),lty=2,col="black",lwd=0.25) + geom_hline(yintercept = -log(0.05,10),lty=2,col="black",lwd=0.25)
top <- 10
top_genes <- bind_rows(ctl_frame %>% filter(Expression == 'Upregulated') %>%  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
                         head(top),
                       ctl_frame %>% filter(Expression == 'Downregulated') %>% arrange(p_val, desc(abs(avg_log2FC))) %>% head(top)
)
p3 <-  p2 + geom_label_repel(data = top_genes, mapping = aes(avg_log2FC, -log(p_val_adj,10), label = rownames(top_genes)),
                             size = 2, max.overlaps = 50, min.segment.length=0, force = 2) +theme_classic() + ggtitle("Differentially Expressed Genes between CTL clusters")

pdf("./Figures/Supplementary/Supp4_CTL_volcano.pdf")
p3 + theme(legend.position = "none")
dev.off()

p1 <- ggplot(nkt_frame, aes(avg_log2FC, -log(p_val_adj,10))) +  geom_point(size = 1) + xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"FDR"))
p2 <- p1 + geom_point(aes(color = Expression), size = 1) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
  geom_vline(xintercept=c(-log(2),log(2)),lty=2,col="black",lwd=0.25) + geom_hline(yintercept = -log(0.05,10),lty=2,col="black",lwd=0.25)
top_genes <- bind_rows(nkt_frame %>% filter(Expression == 'Upregulated') %>%  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
                         head(top),
                       nkt_frame %>% filter(Expression == 'Downregulated') %>% arrange(p_val, desc(abs(avg_log2FC))) %>% head(top)
)
p3 <-  p2 + geom_label_repel(data = top_genes, mapping = aes(avg_log2FC, -log(p_val_adj,10), label = rownames(top_genes)),
                             size = 2, max.overlaps = 50, min.segment.length=0, force = 2) +theme_classic() + ggtitle("Differentially Expressed Genes between NK-like clusters")
pdf("./Figures/Supplementary/Supp4_NK_like_volcano.pdf")
p3 + theme(legend.position = "none")
dev.off()

p1 <- ggplot(em_frame, aes(avg_log2FC, -log(p_val_adj,10))) +  geom_point(size = 1) + xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"FDR"))
p2 <- p1 + geom_point(aes(color = Expression), size = 1) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
  geom_vline(xintercept=c(-log(2),log(2)),lty=2,col="black",lwd=0.25) + geom_hline(yintercept = -log(0.05,10),lty=2,col="black",lwd=0.25)
top_genes <- bind_rows(em_frame %>% filter(Expression == 'Upregulated') %>%  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
                         head(top),
                       em_frame %>% filter(Expression == 'Downregulated') %>% arrange(p_val, desc(abs(avg_log2FC))) %>% head(top)
)
p3 <-  p2 + geom_label_repel(data = top_genes, mapping = aes(avg_log2FC, -log(p_val_adj,10), label = rownames(top_genes)),
                             size = 2, max.overlaps = 50, min.segment.length=0, force = 2) +theme_classic() + ggtitle("Differentially Expressed Genes between EM clusters")

pdf("./Figures/Supplementary/Supp4_EM_volcano.pdf")
p3 + theme(legend.position = "none")
dev.off()


#Top Clusters Venn Diagram/Heatmap
ctl_top50 <- ctl_markers %>% top_n(50, avg_log2FC) %>% rownames()
nkt_top50 <- nkt_markers %>% top_n(50, avg_log2FC) %>% rownames()
em_top50 <- em_markers %>% top_n(50, avg_log2FC) %>% rownames()

top50 <- list("EM1"=em_top50, "NK-like1"=nkt_top50,"CTL1"=ctl_top50)

pdf(paste0(file_prefix,"top50_venn.pdf", height=5, width = 5))
ggVennDiagram(top50, label_alpha = 0, edge_size = 0, label_size = 4, label="count") + 
  scale_fill_distiller(palette = "Reds", direction = 1) + scale_color_manual(values=c("black", "black", "black")) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
  labs(title="Common DEGs between CD8 subsets", subtitle = "Higher Expression in Group1 Clusters")
dev.off()

top_overlap <- intersect(intersect(em_top50, nkt_top50),ctl_top50)
top_overlap <- top_overlap[order(top_overlap)]
top_overlap_genes <- GetAssayData(cd8)[top_overlap,]
mat <- cd8@assays$RNA@data[intersect(top_overlap,rownames(cd8)),]
mat <- groupMeans(mat,unlist(cd8$class13))
mat2 <- mat[,c("N", "Mem1", "Mem2", "MAIT","INF1", "INF2", "CTL INF1", "CTL INF2", "CTL INF3", "Ex", "EM1",
               "CTL1", "NK-like1", "EM2", "CTL2", "NK-like2")]

pdf(paste0(file_prefix,"top_overlap_class13.pdf", width=6, height=3))
pheatmap::pheatmap(t(mat2), cluster_rows=T, show_rownames=T,cluster_cols=F,show_colnames = T,
                   scale="column",col=colour_bk,clustering_method = "complete", cellheight=10, cellwidth = 10, fontsize = 8)
dev.off()

#Density Plot
pdf(paste0(file_prefix,"type_score.pdf", width=6, height=4))
ggplot(cd8@meta.data, aes(x=score, fill=type)) + geom_density(alpha=0.4) + theme_classic() +
  scale_fill_manual(values=anno_colors$type)
dev.off()

#Correlation with CD8-derived score
cd8_hallmark <- cd8@meta.data[,startsWith(colnames(cd8@meta.data), "HALLMARK")]
corr_list <- list()
corr_sig <- list()
for (i in colnames(cd8_hallmark)){
  a <-cor.test(cd8$dysfunction, cd8_hallmark[,i])
  corr_list[[i]] <- a
  if (a$p.value < 0.05){
    corr_sig[[i]] <- a$estimate
  }
}

corr_list <- unlist(corr_sig)
corr_list <- corr_list[order(corr_list)]
corr_frame <- as.data.frame(corr_list)
colnames(corr_frame) <- "Cor"
corr_frame$pathway <- rownames(corr_frame)
corr_frame$pathway <- str_remove(corr_frame$pathway, ".cor")
corr_frame$pathway <- str_remove(corr_frame$pathway, "HALLMARK_")
corr_frame$pathway <- str_replace_all(corr_frame$pathway, "_", " ")


cd8_go_bp <- cd8@meta.data[,startsWith(colnames(cd8@meta.data), "GOBP")]
corr_list <- list()
corr_sig <- list()
for (i in colnames(cd8_go_bp)){
  a <-cor.test(cd8$dysfunction, cd8_go_bp[,i])
  corr_list[[i]] <- a
  if (a$p.value < 0.01){
    corr_sig[[i]] <- a$estimate
  }
}

corr_list <- unlist(corr_sig)
corr_list <- corr_list[order(corr_list)]
corr_frame <- as.data.frame(corr_list)
colnames(corr_frame) <- "Cor"
corr_frame$pathway <- rownames(corr_frame)
corr_frame$pathway <- str_remove(corr_frame$pathway, ".cor")
corr_frame$pathway <- str_remove(corr_frame$pathway, "GOBP_")
corr_frame$pathway <- str_replace_all(corr_frame$pathway, "_", " ")

corr_sub <- corr_frame[abs(corr_frame$Cor) > 0.5,]


h_cor <- hallmark_cor %>% arrange(cor) %>% mutate(pathway=factor(pathway,levels=pathway)) %>% 
  ggplot(aes(x=pathway, y=cor, fill=cor)) + labs(x="Hallmark", y=NULL) +
  geom_bar(stat="identity") + ylim(-0.62,0.6) +
  scale_fill_gradientn(colours =myPalette(100), breaks=c(-1,0,1), limits=c(-1,1)) + coord_flip() + theme_bw()

go_cor <- gobp_cor %>% arrange(cor) %>% mutate(pathway=factor(pathway,levels=pathway)) %>% 
  ggplot(aes(x=pathway, y=cor, fill=cor)) + labs(x="GO Biological Processes", y=NULL) +
  geom_bar(stat="identity") + ylim(-0.62,0.6) + 
  scale_fill_gradientn(colours =myPalette(100), breaks=c(-1,0,1), limits=c(-1,1)) + coord_flip() + theme_bw()


pdf(paste0(file_prefix,"cd8_score_correlations.pdf", width=12, height=5.5))
h_cor / go_cor + plot_layout(guides="collect")
dev.off()


#Correlation with Pseudotime

expression_matrix <- FetchData(cd8, vars = c("score", "pseudotime", "class13", "type", "Patient"))
expression_matrix <- expression_matrix[!is.infinite(expression_matrix$pseudotime),]

#Gene Correlation
df <- expression_matrix %>% group_by(Patient,.drop = FALSE) %>% dplyr::summarise(cd8_score=mean(score),pseudotime=mean(pseudotime))
table(df$type)
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
palOut <- ColAssign(unique(df$class13))

pdf(paste0(file_prefix,"cd8_score_pseudotime.pdf", height=7, width=5.5))
ggplot(df, aes(y=cd8_score, x=pseudotime,shape=type)) + geom_point(aes(size=2,color=class13))+
  geom_smooth(se=FALSE,method = lm, linetype="dashed", size=0.5) + scale_color_manual(values=class13.colors, name="Clusters")+ 
  theme_classic() + stat_cor() + xlim(0,40)+ facet_grid(type~.) + theme(axis.title=element_text(size=18))
dev.off()
#cor.test(df$dysfunction,df$pseudotime)

df_sub <- df[df$type %in% c("NewlyDx", "RxNaive"),]
w <- wilcox.test(cd8_score ~ type, data=df_sub)
p.adjust(w$p.value,n=length(df_sub$cd8_score), method="BH")

#GSEA plots
pdf("./Figures/Supplementary/Supp4_hallmark_infa.pdf")
plotEnrichment(pathwaysH[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], all_ranks) + labs(title="Hallmark IFNa Response")
dev.off()

pdf("./Figures/Supplementary/Supp4_hallmark_hypoxia.pdf")
plotEnrichment(pathwaysH[["HALLMARK_HYPOXIA"]], all_ranks) + labs(title="Hallmark Hypoxia Response")
dev.off()

pdf("./Figures/Supplementary/Supp4_gobp_oxphos_fgsea.pdf")
plotEnrichment(pathwaysGO_BP[["GOBP_OXIDATIVE_PHOSPHORYLATION"]], all_ranks) + labs(title="GO Biological Processes Oxidative Phosphorylation")
dev.off()

pdf("./Figures/Supplementary/Supp4_gobp_aerobic_fgsea.pdf")
plotEnrichment(pathwaysGO_BP[["GOBP_AEROBIC_RESPIRATION"]], all_ranks) + labs(title="GO Biological Processes Aerobic Respiration")
dev.off()
