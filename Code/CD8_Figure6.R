############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Figure 6
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: All code should be in here except for print to pdf commands
## Need to change the TCR folder to point to TCR on the OneDrive
############################################

library(immunarch); library(Platypus); library(alakazam)
library(data.table); library(stringr)
library(circlize); library(scales)
library(Seurat)

setwd("~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/")
source("CD8_Functions.R")
cd8 <- readRDS("./20220912_cd8.rds")

########## Preliminary Analysis ##########
setwd("~/Documents/T_cells/")
tcr_list <- list()
for (i in patient_levels){
  tcr_list[[i]] <- read.csv(paste0("./TCR/",i,".csv"))
  }
setwd("~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/")

for(patient in names(tcr_list)){
  if (startsWith(patient, "NLBM")){tcr_list[[patient]]$barcode <- paste0(patient, "_",tcr_list[[patient]]$barcode)}
  else {tcr_list[[patient]]$barcode <- paste0(patient, "A_",tcr_list[[patient]]$barcode)}
}

tcr_df <- do.call(rbind, tcr_list)

data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "__")
}

tcr_dt <- as.data.table(tcr_df)
tcr_dt_collapsed <- tcr_dt[, {lapply(.SD, data_concater)} , by=barcode]
rownames(tcr_dt_collapsed) <- tcr_dt_collapsed$barcode

barcode_split<-unlist(str_split(tcr_dt_collapsed$barcode, "_"))
barcode_split <- barcode_split[seq(1,2*nrow(tcr_dt_collapsed),by=2)]
tcr_dt_collapsed$Patient <- barcode_split
# 
tcr_dt_collapsed$raw_clonotype_id <- str_remove(tcr_dt_collapsed$raw_clonotype_id, "__")
#tcr_dt_collapsed$barcode <- sub(".*?_","",tcr_dt_collapsed$barcode)

sample_levels <- c("NLBM4", "NLBM5", "NLBM6", "PT1A",   "PT2A",   "PT3A",   "PT4A",
                   "PT5A",   "PT6A",   "PT7A",   "PT8A",   "PT9A",  "PT10A",  "PT11A",  "PT12A",  "PT13A",  "PT14A",
                   "PT15A",  "PT16A",  "PT17A",  "PT18A",  "PT19A",  "PT20A",  "PT21A",  "PT22A",  "PT23A",  "PT24A",
                   "PT25A",  "PT26A",  "PT27A",  "PT28A",  "PT29A",  "PT30A", "PT32A")
sample_type <- c("Healthy", "Healthy", "Healthy", "RelRef",   "RelRef",   "RelRef",   "RelRef",
                 "RelRef",  "RelRef",   "RelRef",  "RelRef",   "NewlyDx",  "NewlyDx", "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",
                 "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",
                 "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx",  "NewlyDx", "NewlyDx")
names(sample_type) <- sample_levels
tcr_dt_collapsed$type <- sample_type[match(tcr_dt_collapsed$Patient, names(sample_type))]

#subset barcodes of CD8+ TCRs from the CD8 (w/ TCR) object
cd8_barcodes <- cd8$barcode[complete.cases(cd8$barcode)]
#subset data table used for analysis
tcr_df_cd8 <- subset(tcr_df, subset = barcode %in% cd8_barcodes)
tcr_dt_cd8 <- subset(tcr_dt_collapsed, subset = barcode %in% cd8_barcodes)

tcr_dt_cd8$type <- sample_type[match(tcr_dt_cd8$Patient,names(sample_type))]
tcr_dt_cd8$class13 <- cd8$class13[match(tcr_dt_cd8$barcode, cd8$barcode)]
tcr_dt_cd8$pt_clon <- paste0(tcr_dt_cd8$Patient, "_", tcr_dt_cd8$raw_clonotype_id)

samples <- sample_levels
patient_levels <- c("NLBM4", "NLBM5", "NLBM6", "PT1",   "PT2",   "PT3",   "PT4",
             "PT5",   "PT6",   "PT7",   "PT8",   "PT9",  "PT10",  "PT11",  "PT12",  "PT13",  "PT14",
             "PT15",  "PT16",  "PT17",  "PT18",  "PT19",  "PT20",  "PT21",  "PT22",  "PT23",  "PT24",
             "PT25",  "PT26",  "PT27",  "PT28",  "PT29",  "PT30", "PT32")
names(samples) <- patient_levels

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


########## A: Alpha Diversity Curve ##########
type_curve <- alphaDiversity(data.frame(tcr_dt_cd8[,c("barcode", "raw_clonotype_id", "type")]), 
                             group="type", clone="raw_clonotype_id")
pdf("./Figures/Figure6/alpha_diversity.pdf", height=3,width=4)
plot(type_curve, colors=anno_colors$type)
dev.off()

########## B: Rank Abundance #########
abundance <- estimateAbundance(data.frame(tcr_dt_cd8[,c("barcode", "raw_clonotype_id", "type")]), 
                               clone="raw_clonotype_id", group="type")
pdf("./Figures/Figure6/rank_aubndance.pdf", height=3, width=4)
plot(abundance, legend_title="Type",xlim=c(1,150), colors=anno_colors$type)
dev.off()

########## C: Expansion Stacked Bar Plot ##########
cd8_tcr_df <- cd8@meta.data[c("type","cloneType")]
cd8_tcr_df <- cd8_tcr_df[complete.cases(cd8_tcr_df$cloneType),]

CT_type <- table(cd8_tcr_df$type, cd8_tcr_df$cloneType)
CT_type <- as.data.frame(CT_type/rowSums(CT_type))
#CT_type_melt <- melt(CT_type)

pdf("./Figures/Figure6/expansion_stacked_bar.pdf", width=2, height=4)
ggbarplot(CT_type, x="Var1", y="Freq", fill="Var2", color="Var2", palette=colorblind_vector(5)) + 
  theme(legend.position = "none", axis.text.x=element_text(angle=90)) + xlab(NULL) + ylab("Percent")
dev.off()

#use similar strategy to the boxplot for cd8 cluster percentage
#need the dataframe with ALL cells not just the summary table!!!!
wilcox_stat <- CT_type %>% group_by(Var2) %>% pairwise_wilcox_test(Freq~Var1, p.adjust.method = "BH") %>% 
  add_significance("p.adj")

########## D: Clonal Expansion UMAP ##########
pdf("./Figures/Figure6/expansion_umap.pdf", width=8, height=2.5)
DimPlot(cd8, group.by = "cloneType", shuffle=T, split.by="type") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + theme_classic() + 
  theme(plot.title = element_blank(), legend.position = "right")
dev.off()

df <- data.frame(FetchData(cd8, c("type", "cloneType")))
k <- kruskal.test(cloneType ~ type, data=df)
p.adjust(k$p.value)

########## E: CD8-Score/Expansion Density Plot ##########

pdf("./Figures/Figure6/expansion_density.pdf")
ggplot(cd8@meta.data, aes(x=score, fill=cloneType)) + geom_density(alpha=0.4) + theme_classic() +
  scale_fill_manual(values=colorblind_vector(5)) + theme(legend.position = "none") 
dev.off()


df <- data.frame(FetchData(cd8, c("score", "cloneType")))
k <- kruskal.test(cloneType ~ score, data=df)
p.adjust(k$p.value)

########## F: Circos Plot ##########
library(circlize); library(scRepertoire)
cd8$class13 <- factor(cd8$class13, levels=c(c("N", "Mem1", "Mem2","EM1", "EM2", "INF1", "INF2", 
                                              "CTL1", "CTL2", "CTL INF1", "CTL INF2", "CTL INF3", 
                                              "NK-like1", "NK-like2", "MAIT", "Ex")))

cd8_rr <- subset(cd8, subset=type %in% "RelRef")
circles1 <- getCirclize(cd8_rr, group.by = "class13")

cd8_newdx <- subset(cd8, subset=type %in% "NewlyDx")
circles2 <- getCirclize(cd8_newdx, group.by = "class13")

cd8_healthy <- subset(cd8, subset=type %in% "Healthy")
circles3 <- getCirclize(cd8_healthy, group.by = "class13")

col_fun = colorRamp2(c(1,300), c("gray", "red"), transparency = 0.5)

pdf("./Figures/Figure6/circos_type.pdf", width=18, height=6)
l <- layout(matrix(c(1,2,3),ncol=3))
circos.par(gap.after = c("N"=10, "Mem1"=1, "Mem2"=10, "EM1"=1,"EM2"=1,"INF1"=1,"INF2"=10,
                         "CTL1"=1,"CTL2"=1,"CTL INF1"=1, "CTL INF2"=1,"CTL INF3"=10, "NK-like1"=1,
                         "NK-like2"= 10, "MAIT"= 10, "Ex"= 10))
circlize::chordDiagram(circles3,self.link = 1, grid.col=anno_colors$class13, col=col_fun, directional = 1, 
                       direction.type = "arrows",  link.arr.type = "big.arrow", link.sort = TRUE, link.decreasing = TRUE)
circlize::chordDiagram(circles2,self.link = 1, grid.col=anno_colors$class13, col=col_fun, directional = 1, 
                       direction.type = "arrows",  link.arr.type = "big.arrow", link.sort = TRUE, link.decreasing = TRUE)
circlize::chordDiagram(circles1,self.link = 1, grid.col=anno_colors$class13, col=col_fun, directional = 1, 
                       direction.type = "arrows",  link.arr.type = "big.arrow", link.sort = TRUE, link.decreasing = TRUE)
dev.off()

########## G: Cluster Overlap by Type ##########

#clonal overlap
cH <- expression2List(cd8_healthy, split.by = "class13")
cH_df <- clonalOverlap(cH, cloneCall="aa", method="morisita", exportTable = TRUE)
cH_df <- cH_df %>% pivot_longer(cols=c(1:16))
cH_df <- cH_df[complete.cases(cH_df),]

cR <- expression2List(cd8_rr, split.by = "class13")

cR_df <- clonalOverlap(cR, cloneCall="aa", method="morisita", exportTable = TRUE)
cR_df <- cR_df %>% pivot_longer(cols=c(1:14))
cR_df <- cR_df[complete.cases(cR_df),]

cN <- expression2List(cd8_newdx, split.by = "class13")
cN_df <- clonalOverlap(cN, cloneCall="aa", method="morisita", exportTable = TRUE)
cN_df <- cN_df %>% pivot_longer(cols=c(1:16))
cN_df <- cN_df[complete.cases(cN_df),]

pdf("./Figures/Additional/clonal_overlap_morisita.pdf", width=24, height=8)
clonalOverlap(cH, cloneCall="aa", method="morisita") + clonalOverlap(cN, cloneCall="aa", method="morisita") + 
  clonalOverlap(cR, cloneCall="aa", method="morisita") + plot_layout(guides="collect") & 
  scale_fill_gradientn(colors = c("blue","white","red"), breaks=c(0,0.25,0.5,0.75,1),
                       limits=c(0,1), na.value = "white") & theme(axis.text.x=element_text(angle=90))
dev.off()


########## I: STARTRAC analysis ##########
Idents(cd8) <- cd8$class13
startrac <- StartracDiversity(cd8, 
                              type = "type", 
                              sample = "orig.ident", 
                              by = "overall", exportTable = TRUE)


startrac$type <- "RxNaive"
startrac[startrac$aid %in% c("NLBM4", "NLBM5", "NLBM6"),]$type <- "Healthy"
startrac[startrac$aid %in% c("PT1A",   "PT2A",   "PT3A",   "PT4A",
                             "PT5A",   "PT6A",   "PT7A",   "PT8A"),]$type <- "RelRef"

melted <- melt(startrac)
melted <- melted[melted$variable %in% c("tran"),]
melted$type <- factor(melted$type, c("Healthy", "RxNaive", "RelRef"))

melted_sub <- melted[complete.cases(melted$value),]
melted_sub$majorCluster <- factor(melted_sub$majorCluster, c("N", "Mem1", "Mem2", "EM1", "EM2","INF1", "INF2",
                                     "CTL1", "CTL2", "CTL INF1", "CTL INF2", "CTL INF3", 
                                     "NK-like1", "NK-like2", "MAIT", "Ex"))


pdf("./Figures/Figure6/startrac.pdf", height=2, width=6)
ggplot(melted_sub, aes(x = majorCluster, y = value)) + 
  geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0) + theme_classic() +
  facet_grid(.~ type) + ylab("Index Score") + scale_fill_manual(values=class13.colors) + 
  guides(fill = "none") + theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=60, hjust=1, size=6), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                strip.background = element_blank(),
                                panel.border = element_rect(colour = "black", fill = NA))
dev.off()


k <- kruskal.test(value ~ type, data=melted)
p.adjust(k$p.value)
