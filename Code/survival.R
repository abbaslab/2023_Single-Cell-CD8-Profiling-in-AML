############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Survival Analysis
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas, Bofei Wang
## Notes: 
############################################

library(GSVA); library(survminer); library(survival)
library(Seurat); library(stringr)
library(ComplexHeatmap)

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"
source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS(paste0(file_prefix,"cd8.rds"))

cd8_derived_sig <- c("RPL41", "MTRNR2L12", "TXNIP", "CD2", "PRF1", "PSMB9", "LIMD2", "LTB", "GIMAP1", "GIMAP4", 
                           "CX3CR1", "TRAF3IP3", "GIMAP7", "GIMAP5", "DENND2D", "CLIC1", "MYL12A", "S100A11", "RAC2", 
                           "PSMB10", "ARPC1B", "CORO1A", "PCBP1", "UCP2", "PLEK")

ref_meta <- read.csv("~/Documents/T_cells/ssGSEA/BeatAML/MetaData_All_v4.csv", row.names = 1)

#Preliminary Analysis
tcga_meta<- ref_meta %>% filter(Cohort=="TCGA") %>% filter(Adult == TRUE) %>% 
  filter(Healthy_Disease == "D" & At_Diagnosis != FALSE)  %>%
  filter(FAB_Simple != "M3")

tcga_meta$short_name <- str_sub(tcga_meta$Sample,1,15)
tcga_meta$OS.time<- cbio_tcga$`Overall Survival (Months)`[match(tcga_meta$short_name, cbio_tcga$`Sample ID`)]
tcga_meta[tcga_meta$OS.time == 0 & !is.na(tcga_meta$OS.time),]$OS.time <- NA
tcga_meta$OS.time2 <- tcga_meta$OS.time * 30.4
tcga_meta$OS <- ifelse(tcga_meta$Dead_Alive=="Dead",1,0)

cd8_sig2 <- c(cd8_derived_sig, "CD8A")

cancer_gene_exp <- tcga_zeng[,tcga_meta$Sample]
keep <- apply(cancer_gene_exp, 1, function(x) length(unique(x)) != 1) #removes anything that had constant expression
cancer_gene_exp <- cancer_gene_exp[keep, ]
cancer_gene_exp <- t(cancer_gene_exp)
cancer_gene_exp <- t(scale(cancer_gene_exp))
gsva_es <- gsva(cancer_gene_exp, list("Score" = cd8_derived_sig, "CD8A-norm" =  cd8_sig2), method="ssgsea")


#Correlation & Heatmap

score_frame <- data.frame(gsva_es)
score_frame <- score_frame["Score",]
colnames(dys) <- gsub("\\.","-",colnames(score_frame))
score_frame <- score_frame[order(score_frame[1,], decreasing = TRUE)]


ha_score <- HeatmapAnnotation(Score = anno_barplot(as.numeric(score_frame[1,]), bar_width = 0.75, 
                                                 border = FALSE, 
                                                 axis = TRUE,
                                                 axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                 gp = gpar(col = NA, fill = "grey20"), show_annotation_name = F, height = unit(4, "cm")))

cancer_score <- cancer_gene_exp[cd8_derived_sig[cd8_derived_sig %in% rownames(tcga_zeng)],tcga_meta$Sample]
cancer_score <- cancer_score[,names(score_frame)]

ht_score <- Heatmap(cancer_score, show_column_names = FALSE, cluster_columns = FALSE, top_annotation = ha_score)

pdf(paste0(file_prefix,"score_heatmap_tcga.pdf", height=10, width=10))
ht_score                        
dev.off()

C <- cor(t(cancer_score))
library(corrplot)

pdf(paste0(file_prefix,"score_gene_correlation.pdf", height=5, width=5))
corrplot(C, tl.cex = 0.5, tl.col="black", col=colorRampPalette(c("blue","white","red"))(100), order='FPC')
dev.off()


#TCGA Survival

tcga_meta$score <- gsva_es["Score",]
tcga_meta$score_norm <- gsva_es["CD8A-norm",]
split.exp <- as.numeric(quantile(tcga_meta$dysfunction,probs=0.5))
tcga_meta$dys <- ifelse(tcga_meta$dysfunction> as.numeric(quantile(tcga_meta$dysfunction,probs=0.5)), 1,0)
norm_fit<- survfit(Surv(OS.time2,Dead_Alive=="Dead") ~ dys, tcga_meta)
coxph(Surv(OS.time2, OS==1) ~ dys, data = tcga_meta)

pdf("./Figures/Survival/cd8_score_tcga.pdf", width=5, height = 6)
ggsurvplot(norm_fit, data=tcga_meta, pval=T, palette = c("blue", "red"),
           legend.labs = c("Low Score", "High Score"), xlab= "Survival Time (Days)",
           pval.coord = c(50,0.05), legend=c(0.75,0.9), legend.title="", title="TCGA LAML Survival Probability", conf.int = FALSE,
           risk.table = TRUE) 
dev.off()

split_norm <- as.numeric(quantile(tcga_meta$dys_norm,probs=0.5))
tcga_meta$dys2 <- ifelse(tcga_meta$dys_norm> as.numeric(quantile(tcga_meta$dys_norm,probs=0.5)), 1,0)
norm_fit2<- survfit(Surv(OS.time2,Dead_Alive=="Dead") ~ dys2, tcga_meta)

pdf("./Figures/Survival/cd8_score_tcga_norm.pdf", width=5, height = 6, onefile = FALSE)
ggsurvplot(norm_fit2, data=tcga_meta, pval=T, palette = c("blue", "red"),
           legend.labs = c("Low Score", "High Score"), xlab= "Survival Time (Days)",
           pval.coord = c(50,0.05), legend=c(0.75,0.9), legend.title="", title="TCGA LAML Survival Probability", conf.int = FALSE,
           risk.table = TRUE)
dev.off()

coxph(Surv(OS.time2, OS==1) ~ dys2, data = tcga_meta)


# Multivariate Analysis

tcga_meta$age_group <- ifelse(tcga_meta$Age > 60, "Greater", "Less")
tcga_meta$age_group <- factor(tcga_meta$age_group, levels=c("Less", "Greater"))
#cancer_survival$cyto.risk <- tcga_meta[match(cancer_survival$bcr_patient_barcode, tcga_meta$Sample)]
tcga_meta$Cytogenetic_risk_incomplete=factor(tcga_meta$Cytogenetic_risk_incomplete,levels = c("Favorable","Intermediate","Unfavorable" ))

covariate_names <- c(`age_group:Greater` ="Age > 60", Cytogenetic_risk_incomplete="Risk Level", dys_norm="CD8-Derived Score", 
                     Percent_Blasts_in_BM = "Blast %", `Cytogenetic_risk_incomplete:Intermediate` = "ELN 2017 Intermediate Cytogenetics",
                     `Cytogenetic_risk_incomplete:Unfavorable` = "ELN2017 Adverse Cytogenetics")

tcga_meta %>%
  analyse_multivariate(vars(OS.time2, OS==1),
                       covariates = vars(age_group, Cytogenetic_risk_incomplete, dys_norm, Percent_Blasts_in_BM),
                       covariate_name_dict = covariate_names) -> result
forest_plot(result,
            factor_labeller = covariate_names,
            endpoint_labeller = c(time="OS"),
            orderer = ~order(HR),
            labels_displayed = c("factor"),
            ggtheme = ggplot2::theme_bw(base_size = 8),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2))


#Hazard Ratio Plot

library(phenoTest)

dys_sig <- c("RPL41", "TXNIP", "CD2", "PRF1", "PSMB9", "LIMD2", "LTB", "GIMAP1", "GIMAP4", 
             "CX3CR1", "TRAF3IP3", "GIMAP7", "GIMAP5", "DENND2D", "CLIC1", "MYL12A", "S100A11", "RAC2", 
             "PSMB10", "ARPC1B", "CORO1A", "PCBP1", "UCP2", "PLEK", "CD8A")

cancer_gene_exp <- tcga_zeng[dys_sig,tcga_meta$Sample]
rownames(tcga_meta) <- tcga_meta$Sample

pheno_surv <- matrix(c("OS", "OS.time2"), ncol=2, byrow=TRUE)
colnames(pheno_surv) <- c('event', 'time')
vars2test <- list(survival = pheno_surv, categorical='Sex', continuous='Age')
exp_set <- ExpressionSet(as.matrix(cancer_gene_exp), phenoData = AnnotatedDataFrame(tcga_meta))
epheno <- phenoTest::ExpressionPhenoTest(exp_set, vars2test, p.adjust.method = "none")
epheno <- pAdjust(epheno,method='BH')

barplotSignifSignatures(epheno[,'OS'],list(dys_sig),alpha=0.99)
barplotSignatures(epheno[,'OS'],list(dys_sig), ylim=c(0,1))

pvals <- heatmapPhenoTest(exp_set,dys_sig,vars2test=vars2test[1],heat.kaplan='heat')

pdf("./Figures/Survival/tcga_smooth_cox.pdf", width=6, height = 6, onefile = FALSE)
smoothCoxph(pData(exp_set)$OS.time2, pData(exp_set)$OS, pData(exp_set)$dysfunction, xlab="CD8-derived Score")
dev.off()
cp <- coxph(Surv(OS.time2, OS==1) ~ dysfunction, data = pData(exp_set))
summary(cp)

#HLA gene correlation

hla=list()
hla$hla1=c("B2M","HLA-A","HLA-B","HLA-C");hla$hla2=c("HLA-DMB","HLA-DMA","HLA-DRB1","HLA-DPA1","HLA-DRA","HLA-DPB1")
gsva_es <- gsva(cancer_gene_exp, hla, method="ssgsea", abs.ranking=F)
tcga_meta$hla1_ssgsea_group <- ifelse(t(gsva_es)[,1]> as.numeric(quantile(t(gsva_es)[,1],probs=0.5)), 1,0)
tcga_meta$hla1_ssgsea_score <- t(gsva_es)[,1]
tcga_meta$hla2_ssgsea_group <- ifelse(t(gsva_es)[,2]> as.numeric(quantile(t(gsva_es)[,2],probs=0.5)), 1,0)
tcga_meta$hla2_ssgsea_score <- t(gsva_es)[,2]
plot.cor(y="hla2_ssgsea_score",ylab="HLA2_ssgsea_score",filename = "HLA2_ssGSEA_Score")


pdf("./Figures/Survival/hla1_score_scatter.pdf", width=7, height=7)
ggscatter(tcga_meta, x="dys_norm", y="hla1_ssgsea_score",add="reg.line",conf.int = TRUE,cor.coef.size = 6, 
          cor.coef = TRUE, cor.method = "pearson",xlab = "CD8-derived Score", ylab = "HLA1 ssGSEA Score",size=1.5,
          add.params = list(color="darkblue", fill="darkblue")) +
  theme(legend.position = "none", axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.ticks = element_blank()) 
dev.off()

#Abbas Survival Analysis
abbas_meta<- ref_meta %>% filter(Cohort=="Abbas") %>% filter(Adult == TRUE) %>% 
  filter(Healthy_Disease == "D" & At_Diagnosis != FALSE)  %>%
  filter(FAB_Simple != "M3" & FAB_Simple != "M6" & FAB_Simple != "M7")

ac <- read.csv("~/OneDrive - Inside MD Anderson/LabMembers/Resources/BulkRNA/AML_logTPM/Abbas_logTPM.csv")
colnames(ac) <- colnames(ac) %>% str_replace_all("\\.","-") %>% str_replace_all("X","")
rownames(ac) <- ac$Genes
abbas_TPM <- ac


ac <- t(ac)
ac <- ac/ac[,"CD8A"]
abbas_counts_scaled <- t(scale(ac))
ac <- abbas_TPM[,abbas_meta$Sample]
abbas_counts_scaled <- t(scale(t(ac)))
keep <- apply(abbas_counts_scaled, 1, function(x) length(unique(x)) != 1) #removes anything that had constant expression
abbas_counts_scaled <- abbas_counts_scaled[keep, ]
gsva_es <- gsva(abbas_counts_scaled, list(dysfunction_signature),method="ssgsea",abs.ranking=F)
abbas_meta$dysfunction <- t(gsva_es)
abbas_meta$dys <- ifelse(t(gsva_es)> as.numeric(quantile(t(gsva_es),probs=0.5)), 1,0)
split.exp=as.numeric(quantile(abbas_meta$dysfunction,probs=0.5, na.rm=TRUE))
#gene.cox <- coxph(Surv(OS.time, OS == 1) ~ escore < split.exp, abbas_meta,x=TRUE)
#hr1 <- smoothHR(data=phei, coxfit=gene.cox)
#plot(hr1, predictor="gene.expi", prob=0, conf.level=0.95)
fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ dysfunction < split.exp, abbas_meta)

pdf("~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/Figures/Additional/abbas_survival.pdf")
ggsurvplot(fit,data=abbas_meta,pval = T,legend.labs = c("High Score", "Low Score"), xlab= "Time (days)",
           pval.coord = c(4000,0.15), legend=c(0.7,0.9), legend.title="", title="Abbas AML Data Survival Probability",
           risk.table = TRUE, conf.int = TRUE)
dev.off()