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

library(miloR)

#Milo UMAP

cd8_milo <- as.SingleCellExperiment(cd8)
cd8_milo <- Milo(cd8_milo)
cd8_milo <- buildGraph(cd8_milo, k = 20, d = 50)
cd8_milo <- makeNhoods(cd8_milo, prop = 0.1, k = 20, d=50, refined = TRUE)

#Nhood size
pdf(paste0(file_prefix,"milo_nhood_size.pdf"), height=5, width=7)
plotNhoodSizeHist(cd8_milo)
dev.off()

cd8_milo <- countCells(cd8_milo, meta.data = data.frame(colData(cd8_milo)), sample="orig.ident")
cd8_milo <- calcNhoodDistance(cd8_milo, d=50)

cd8_design <- data.frame(colData(cd8_milo))[,c("orig.ident", "type")]
cd8_design$type <- as.factor(cd8_design$type)
cd8_design <- distinct(cd8_design)
rownames(cd8_design) <- c(cd8_design$orig.ident)

da_results <- testNhoods(cd8_milo, design = ~ type, design.df = cd8_design)
da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "class13")
da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "type")

pdf(paste0(file_prefix,"milo_type_fraction.pdf"), height=5, width=7)
ggplot(da_results, aes(type_fraction)) + geom_histogram(bins=50)
dev.off()

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


#Milo BeeSwarm Plots
pdf(paste0(file_prefix,"milo_type_DAbeeswarm.pdf"), width = 4, height=1.5)
plotDAbeeswarm(da_results, group.by = "type") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf(paste0(file_prefix,"milo_class13_beeswarm.pdf"), width = 4, height=4.5)
plotDAbeeswarm(da_results, "class13") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

