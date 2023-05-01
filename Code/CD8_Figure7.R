############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Figure 7
## Date: 3/2/2023
## Author: Poonam Desai, Andre Faustino
## Notes: 
############################################

library(data.table); library(stringr)
library(circlize); library(scales)
library(Seurat); library(readr)
library(ggplot2); library(dplyr)
library(patchwork); library(UpSetR)
library(scRepertoire)

setwd("~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/")
source("CD8_Functions.R")
cd8 <- readRDS("./CD8_paper/20220912_cd8.rds")

########## Preliminary Analysis ##########

gliph <- read.table(file="~/Downloads/antigen_specificity/gliph_poonam_members_annotated.tsv", sep = "\t", header=T)
rownames(gliph) <- gliph$sample_id

data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "__")
}

gliph_dt <- as.data.table(gliph)
gliph_dt_collapsed <- gliph_dt[, {lapply(.SD, data_concater)} , by=sample_id]
rownames(gliph_dt_collapsed) <- gliph_dt_collapsed$barcode

gliph_dt_collapsed$group <- str_replace_all(gliph_dt_collapsed$group, "Pathogen-related__UnknownSpecificity", "Pathogen-related")
gliph_dt_collapsed$group <- str_replace_all(gliph_dt_collapsed$group, "Cancer-related/Pathogen-related__Pathogen-related", "Cancer-related/Pathogen-related")
gliph_dt_collapsed$group <- str_replace_all(gliph_dt_collapsed$group, "Pathogen-related__Pathogen-related/UnknownSpecificity", "Pathogen-related")
gliph_dt_collapsed$group <- str_replace_all(gliph_dt_collapsed$group, "Pathogen-related/UnknownSpecificity__UnknownSpecificity", "Pathogen-related/UnknownSpecificity")
gliph_dt_collapsed$tcr_group <- gliph_dt_collapsed$group
gliph_dt_collapsed$group <- NULL

gliph_dt_collapsed$group_umap <- str_replace_all(gliph_dt_collapsed$group_umap, "Ambiguous__Pathogen-related", "Ambiguous")
gliph_dt_collapsed$group_umap <- str_replace_all(gliph_dt_collapsed$group_umap, "Pathogen-related__UnknownSpecificity", "Pathogen-related")

gliph_dt_collapsed$condition <- NULL
gliph_dt_collapsed$patient <- NULL

rownames(gliph_dt_collapsed) <- gliph_dt_collapsed$sample_id
gliph_dt_collapsed$sample_id <- NULL

cd8 <- AddMetaData(cd8, gliph_dt_collapsed)

gliph_dt$class13 <- cd8$class13[match(gliph_dt$sample_id, colnames(cd8))]
gliph_dt_sum <- data.frame(gliph_dt) %>%
  group_by(group_umap, class13) %>%
  dplyr::summarise(median_clonotypes = median(Freq))

########## A: Clonotype Specificity ##########

gliph_cluster_ann <- read_tsv(file = '~/Downloads/antigen_specificity/gliph_cluster_ann.tsv')
gliph_members_ann_summarize<- read_tsv(file = "~/Downloads/antigen_specificity/gliph_members_annotated.tsv")

#

gliph_clusters_ann_summarize <- gliph_cluster_ann %>%
  group_by(group_umap) %>%
  dplyr::summarise(n_cluster = n())

# Annotation barplot

p1 <- ggplot(gliph_clusters_ann_summarize, aes(x = group_umap, y = n_cluster, fill=factor(n_cluster))) +
  geom_col() + geom_text(aes(label=n_cluster)) + scale_fill_manual(values = c("gold", "firebrick4","steelblue", gray(0.3))) +
  labs(x = NULL, y = "# of clusters") +
  scale_y_continuous(trans = 'sqrt') +
  theme_classic() + theme(axis.text.x = element_blank(), legend.position = "none")

#

gliph_members_ann_matrix <- gliph_members_ann_summarize%>%
  group_by(group_umap, condition) %>%
  dplyr::summarise(
    median_clonotypes = median(Freq)
  )

gliph_members_ann_matrix <- gliph_members_ann_matrix %>%
  group_by(group_umap) %>%
  mutate(
    zscore = scale(median_clonotypes)[,1]
  )

# Heatmap abundance
p2 <- ggplot(gliph_members_ann_matrix, aes(x = group_umap, y = condition, fill = zscore)) +
  geom_tile(colour = "grey", size = 0.45) +
  labs(x = NULL, y = NULL, fill = 'Z-score') +
  scale_y_discrete(
    limits = c('RelRef', 'RxNaive', 'Healthy'),
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(low = muted("blue"),
                       mid = "white",
                       high = muted("red")) +
  scale_x_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = 'black', size = 12),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 0.9,
      size = 12
    ),
    axis.ticks = element_line(size = 0.4),
    legend.position = 'bottom'
  )

pdf("./Figures/Additional/gliph_spec.pdf", height=6, width=3)
p1/p2 +
  plot_layout(heights = c(3,1))
dev.off()


########## B: Clonotypes by Cell Type##########

t <- table(cd8$group, cd8$group_umap)
t <- as.data.frame(t)
u <- data.frame(pivot_wider(t, values_from=Freq, names_from = c(Var1)))
rownames(u) <- u$Var2
u$Var2 <- NULL
u <- as.matrix(u)
fisher_test(t(u), simulate.p.value = TRUE)

pdf("./Figures/Additional/spec_type_group.pdf", width=4, height=4)
ggbarplot(t, x="Var1", y="Freq", fill="Var2", color="Var2") + scale_color_manual(values = c("gold", "firebrick4","steelblue", gray(0.3))) +
  scale_fill_manual(values = c("gold", "firebrick4","steelblue", gray(0.3)))+
  theme(axis.text.x=element_text(angle=90)) + xlab(NULL) + ylab("Percent") + scale_x_discrete(limits=c("N", "Mem", "EM", "CTL", "MAIT", "Ex")) +
  stat_compare_means()
dev.off()

########## C: Clonotypes by Expansion ##########

t <- table(cd8$cloneType, cd8$group_umap)
t <- as.data.frame(t)

xtab <- table(cd8$cloneType, cd8$group_umap)
fisher_test(xtab)

pdf("./Figures/Additional/spec_clonType.pdf", width=5, height=6)
ggbarplot(t, x="Var1", y="Freq", fill="Var2", color="Var2") + scale_color_manual(values = c("gold", "firebrick4","steelblue", gray(0.3))) +
  scale_fill_manual(values = c("gold", "firebrick4","steelblue", gray(0.3)))+
  theme(axis.text.x=element_text(angle=90), legend.position = "right") + xlab(NULL) + ylab("Percent")
dev.off()
########## D: Clonotype Sharing Upset ##########

gliph_cdr3_ann_summarise <- gliph_members_ann_summarise %>%
  ungroup() %>%
  distinct(TcRb, condition, group_umap) %>%
  mutate(
    Healthy = ifelse(condition == "Healthy", 1, 0),
    RelRef = ifelse(condition == "RelRef", 1, 0),
    RxNaive = ifelse(condition == "RxNaive", 1, 0),
  )

gliph_cdr3_ann_summarise <- gliph_cdr3_ann_summarise %>%
  group_by(TcRb, group_umap) %>%
  dplyr::summarise(
    Healthy = sum(Healthy),
    RelRef = sum(RelRef),
    RxNaive = sum(RxNaive),
  )

gliph_cdr3_ann_summarise <- gliph_cdr3_ann_summarise %>%
  as.data.frame()

categories <- unique(gliph_cdr3_ann_summarise$group_umap)

#

pdf("./Figures/Additional/tcr_upset.pdf", width=7, height=9)
UpSetR::upset(
  gliph_cdr3_ann_summarise,
  query.legend = "top",
  sets.bar.color = 'grey26',
  text.scale = 2,
  nintersects = 10,
  sets = c('RelRef', 'RxNaive', 'Healthy'),
  mainbar.y.label = '# of Clonotypes',
  keep.order = TRUE,
  queries = list(
    list(query = elements,
         params = list(
           "group_umap", categories), color = grey(0.3), query.name = "Unknown Specificity", active = T),
    list(query = elements,
         params = list(
           "group_umap", categories[-1]), color = "steelblue", query.name = "Pathogen-related", active = T),
    list(query = elements,
         params = list(
           "group_umap", categories[-c(1:2)]), color = "gold", query.name = "Ambiguous", active = T),
    list(query = elements,
         params = list(
           "group_umap", "Cancer-related"), color = "firebrick4", query.name = "Cancer-related", active = T)
  )
)
dev.off()

# Barplot shareability

gliph_cdr3_shareability <- gliph_cdr3_ann_summarise %>%
  mutate(
    is_shared = ifelse(
      Healthy + RelRef + RxNaive > 1, 'Shared', 'Exclusive'
    )
  )

gliph_cdr3_shareability <- gliph_cdr3_shareability %>%
  group_by(is_shared) %>%
  dplyr::summarise(
    n_cases = n()
  )

#

pdf("./Figures/Additional/Shareability.pdf", width=3, height=4)
ggplot(gliph_cdr3_shareability, aes(x = is_shared, y = n_cases)) +
  geom_col() + geom_text(aes(label=n_cases)) +
  labs(x = NULL, y = "# of CDR3β clonotypes") +
  scale_y_continuous(trans = 'sqrt') +
  theme_classic() +
  theme(
    axis.text = element_text(colour = 'black', size = 12),
    axis.title.y = element_text(size = 14)
  )
dev.off()


########## E: Clonotype Sharing RelRef ##########

gliph_rr <- subset(gliph_poonam_members_ann_summarise, condition == "RelRef")

gliph_rr_summarise <- gliph_rr %>%
  ungroup() %>%
  distinct(TcRb, patient, group_umap) %>%
  mutate(
    PT1 = ifelse(patient == "PT1A", 1, 0),
    PT2 = ifelse(patient == "PT2A", 1, 0),
    PT3 = ifelse(patient == "PT3A", 1, 0),
    PT4 = ifelse(patient == "PT4A", 1, 0),
    PT5 = ifelse(patient == "PT5A", 1, 0),
    PT6 = ifelse(patient == "PT6A", 1, 0),
    PT7 = ifelse(patient == "PT7A", 1, 0),
    PT8 = ifelse(patient == "PT8A", 1, 0)
  )

gliph_rr_summarise <- gliph_rr_summarise %>%
  group_by(TcRb, group_umap) %>%
  dplyr::summarise(
    PT1 = sum(PT1),
    PT2 = sum(PT2),
    PT3 = sum(PT3),
    PT4 = sum(PT4),
    PT5 = sum(PT5),
    PT6 = sum(PT6),
    PT7 = sum(PT7),
    PT8 = sum(PT8)
  )

gliph_rr_summarise <- gliph_rr_summarise %>%
  as.data.frame()

categories <- unique(gliph_rr_summarise$group_umap)
#PT2-PT6 = Pathogen
#PT1-PT2 = Unknown

pdf("./Figures/Additional/tcr_rr_upset.pdf", width=9, height=4.5)
UpSetR::upset(
  gliph_rr_summarise,
  query.legend = "top",
  sets.bar.color = 'grey26',
  text.scale = 1,
  nintersects = 10,
  sets = c("PT1", "PT2", "PT3", "PT4", "PT5", "PT6", "PT7", "PT8"),
  mainbar.y.label = '# of Clonotypes',
  keep.order = FALSE,
  queries = list(
    list(query = elements,
         params = list(
           "group_umap", categories), color = grey(0.3), query.name = "Unknown Specificity", active = T),
    list(query = elements,
         params = list(
           "group_umap", categories[-1]), color = "steelblue", query.name = "Pathogen-related", active = T),
    list(query = elements,
         params = list(
           "group_umap", categories[-c(1:2)]), color = "gold", query.name = "Ambiguous", active = T),
    list(query = elements,
         params = list(
           "group_umap", "Cancer-related"), color = "firebrick4", query.name = "Cancer-related", active = T)
  )
)
dev.off()

########## F: Clonotype Sharing NewDx ##########

gliph_new <- subset(gliph_poonam_members_ann_summarise, condition == "RxNaive")

gliph_new_summarise <- gliph_new %>%
  ungroup() %>%
  distinct(TcRb, patient, group_umap) %>%
  mutate(
    PT9 = ifelse(patient == "PT9A", 1, 0),
    PT10 = ifelse(patient == "PT10A", 1, 0),
    PT11 = ifelse(patient == "PT11A", 1, 0),
    PT12 = ifelse(patient == "PT12A", 1, 0),
    PT13 = ifelse(patient == "PT13A", 1, 0),
    PT14 = ifelse(patient == "PT14A", 1, 0),
    PT15 = ifelse(patient == "PT15A", 1, 0),
    PT16 = ifelse(patient == "PT16A", 1, 0),
    PT17 = ifelse(patient == "PT17A", 1, 0),
    PT18 = ifelse(patient == "PT18A", 1, 0),
    PT19 = ifelse(patient == "PT19A", 1, 0),
    PT20 = ifelse(patient == "PT20A", 1, 0),
    PT21 = ifelse(patient == "PT21A", 1, 0),
    PT22 = ifelse(patient == "PT22A", 1, 0),
    PT23 = ifelse(patient == "PT23A", 1, 0),
    PT24 = ifelse(patient == "PT24A", 1, 0),
    PT25 = ifelse(patient == "PT25A", 1, 0),
    PT26 = ifelse(patient == "PT26A", 1, 0),
    PT27 = ifelse(patient == "PT27A", 1, 0),
    PT28 = ifelse(patient == "PT28A", 1, 0),
    PT29 = ifelse(patient == "PT29A", 1, 0),
    PT30 = ifelse(patient == "PT30A", 1, 0),
    PT32 = ifelse(patient == "PT32A", 1, 0)
  )

gliph_new_summarise <- gliph_new_summarise %>%
  group_by(TcRb, group_umap) %>%
  dplyr::summarise(
    PT9 = sum(PT9),
    PT10 = sum(PT10),
    PT11 = sum(PT11),
    PT12 = sum(PT12),
    PT13 = sum(PT13),
    PT14 = sum(PT14),
    PT15 = sum(PT15),
    PT16 = sum(PT16),
    PT17 = sum(PT17),
    PT18 = sum(PT18),
    PT19 = sum(PT19),
    PT20 = sum(PT20),
    PT21 = sum(PT21),
    PT22 = sum(PT22),
    PT23 = sum(PT23),
    PT24 = sum(PT24),
    PT25 = sum(PT25),
    PT26 = sum(PT26),
    PT27 = sum(PT27),
    PT28 = sum(PT28),
    PT29 = sum(PT29),
    PT30 = sum(PT30),
    PT32 = sum(PT32)
  )

gliph_new_summarise <- gliph_new_summarise %>%
  as.data.frame()

categories <- unique(gliph_new_summarise$group_umap)
#PT9-PT25 = Pathogen
#PT9-PT21 = Pathogen
#PT26-PT30 = Pathogen
#PT10-PT12 = Pathogen
#PT10-PT30 = Pathogen
#PT29-PT14 = Pathogen
#PT29-PT26 = Ambiguous
#PT26-PT9-PT25 = Pathogen

pdf("./Figures/Additional/tcr_newdx_upset.pdf", width=14, height=8)
UpSetR::upset(
  gliph_new_summarise,
  query.legend = "top",
  sets.bar.color = 'grey26',
  text.scale = 1,
  nintersects = 35,
  sets = c("PT9",  "PT10",  "PT11",  "PT12",  "PT13",  "PT14",
           "PT15",  "PT16",  "PT17",  "PT18",  "PT19",  "PT20",  "PT21",  "PT22",  "PT23",  "PT24",
           "PT25",  "PT26",  "PT27",  "PT28",  "PT29",  "PT30", "PT32"),
  mainbar.y.label = '# of Clonotypes',
  keep.order = FALSE,
  queries = list(
    list(query = elements,
         params = list(
           "group_umap", categories), color = grey(0.3), query.name = "Unknown Specificity", active = T),
    list(query = elements,
         params = list(
           "group_umap", categories[-1]), color = "steelblue", query.name = "Pathogen-related", active = T),
    list(query = elements,
         params = list(
           "group_umap", categories[-c(1:2)]), color = "gold", query.name = "Ambiguous", active = T),
    list(query = elements,
         params = list(
           "group_umap", "Cancer-related"), color = "firebrick4", query.name = "Cancer-related", active = T)
  )
)
dev.off()
