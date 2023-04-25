############################################
## Project: CD8+ T-cell Analysis
## Script Purpose: Code for Figure 1
## Date: 8/26/2022
## Author: Poonam Desai, Hussein Abbas
## Notes: 
############################################

file_prefix <- "~/OneDrive - Inside MD Anderson/LabMembers/Poonam/Grants_Manuscripts/Manuscripts/CD8_paper/"

source(paste0(file_prefix,"Code/CD8_Functions.R"))
cd8 <- readRDS("cd8.rds")

########## B: Oncoprint #########
cytogenetics <- read_xlsx(paste0(file_prefix,"Metadata_20220708.xlsx"), sheet="Poonam_Metadata")
cytogenetics <- data.frame(cytogenetics)
rownames(cytogenetics) <- cytogenetics$SampleID

cyto_ht <- cytogenetics[,c(3,4,6,7,8,9,11,12,13,14,15,16,17,18,22,23,24,39,40)]
cyto_ht <- subset(cyto_ht, subset= Type %in% c("NewlyDx", "RelRef"))
cyto_ht$Type <- NULL
cyto_ht <- t(cyto_ht)
cyto_ht <- cyto_ht[,c(1,7,2,4,8,5,6,3,12,19,30,18,20,23,25,31,9,13,26,17,22,11,28,14,15,21,24,27,16,29,10)]
rownames(cyto_ht) <- c("Karyotype", "Del(5/5q)", "Del(7/7q)", "Trisomy 8", "Del(17/17p)", "EZH2",
                       "TP53", "FLT3", "NPM1", "RAS", "SF3B1", "IDH1", "IDH2", "DNMT3A", "TET2", "ASXL1", "ELN 2017", "Status")

lgd1 = Legend(labels = c("Mutation Present", "No Mutation", "Not Tested"), legend_gp = gpar(fill = c("plum3", "lightgray", "snow4")), 
              title = "Mutation")
lgd2 = Legend(labels=c("Diploid", "Complex", "Not Complex"), legend_gp = gpar(fill=c("rosybrown","palevioletred", "lightpink")), title="Karyotype")
lgd3 = Legend(labels=c("Chromosome 5q Deletion", "Chromosome 7 Deletion", "Chromosome 7q Deletion", "Chromosome 17 Deletion", "Chromosome 8 Trisomy"), 
              legend_gp=gpar(fill=c("peachpuff1", "steelblue3", "lightskyblue", "indianred2", "tan2")), title= "Chromosome Alterations")
lgd4= Legend(labels=c("Favorable", "Intermediate", "Adverse"), legend_gp = gpar(fill=c("darkseagreen3","khaki", "paleturquoise")), title="ELN 2017 Criteria")
lgd5= Legend(labels=c("Secondary", "Denovo"), legend_gp = gpar(fill=c("wheat3", "seashell")), title="AML Status")
pd = packLegend(lgd2, lgd3, lgd1, lgd4, lgd5, direction = "horizontal", gap = unit("2", "cm"))

pdf(paste0(file_prefix,"patient_cytogenetics.pdf"), width=12, height=7.5)
Heatmap(cyto_ht, col = c("1" = "plum3", "0"="lightgray", "NP" = "snow4", "Diploid"="rosybrown",
                         "Not Complex"="lightpink","Complex"="palevioletred", "del7"="steelblue","del7q"="steelblue", 
                         "del5q"="steelblue", "Trisomy8"="steelblue", "del17_17p"="steelblue", "No"="lightgray", "Favorable"="darkseagreen3", 
                         "Intermediate"="khaki", "Adverse"="paleturquoise", "Secondary"= "wheat3", "Denovo"="seashell1"), column_title = "Patients", row_title=" ", 
        row_split=c(rep("A", 5), rep("B",11), rep("C",2)),column_split = c(rep("A",8), rep("B",23)),column_gap = unit(3, "mm"),
        border=FALSE, width=unit(7, "in"),height=unit(4,"in"), rect_gp = gpar(col = "white", lwd = 2), show_heatmap_legend = FALSE)
draw(pd, just=c("center","bottom"), y=unit(0.5, "cm"))
dev.off()

########## C: Group UMAP ##########

pdf(paste0(file_prefix,"group_umap.pdf"), width=4,height=4)
DimPlot(cd8, group.by="group", cols=anno_colors$group, label=TRUE, label.box = TRUE, label.size = 3) + ggtitle(NULL) & NoLegend() & NoAxes()
dev.off()

########## D: All Clusters UMAP ##########

pdf(paste0(file_prefix,"class13_umap.pdf"), width=3, height=3)
DimPlot(cd8, group.by="class13", cols=anno_colors$class13, label=FALSE) + ggtitle(NULL) & NoLegend() & NoAxes()
dev.off()

########## E: Group Dot Plot ##########
can_markers <- c("CCR7", "LEF1", "SELL", "IL7R", "LTB", "GZMB", "GNLY", "PRF1", "NKG7", "GZMK", 
                 "PDCD1", "LAG3", "TIGIT", "SLC4A10", "TRAV1-2")

pdf(paste0(file_prefix,"canonical_group_dot_plot.pdf"), height=4, width=2.5)
dot_plot_x(cd8, rev(can_markers), group = "group")
dev.off()

########## F: Group Distribution by Patient Type ##########

group_type_df <- as.data.frame(table(cd8$group, cd8$type))
kruskal.test(Freq~Var1, data=group_type_df)


pdf(paste0(file_prefix,"group_stacked_bar.pdf"), width=3, height=4)
ggplot(cd8@meta.data, aes(x=type, fill=group, color=group)) + geom_bar(position = position_fill(reverse = TRUE), alpha=0.8, width=0.5) + 
  theme_classic() + scale_fill_manual(values=group.colors, name= "Cell Types") + scale_color_manual(values=anno_colors$group, name= "Cell Types") +
  theme( axis.text.x = element_text(size=10, angle=90), plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  labs(x=NULL, y="Percent")
dev.off()


########## G: Pie Chart ##########
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

########## H: Exhaustion Expression ##########
zhang_ex <- c("PDCD1", "TOX", "CXCL13", "TIGIT", "CTLA4", "TNFRSF9", "HAVCR2", "LAG3")

mat <- FetchData(cd8, zhang_ex)
mat <- groupMeans(t(mat),unlist(cd8$group))

pdf(paste0(file_prefix,"zhang_ex_group.pdf"), height=3.5, width=3)
ht <- Heatmap(t(mat), col = colour_bk, show_row_dend = F,rect_gp = gpar(col = "white", lwd = 1),
              width = nrow(mat)*unit(5, "mm"), height = ncol(mat)*unit(5, "mm"), row_names_gp =gpar(fontsize=8),
              column_names_gp =gpar(fontsize=8),heatmap_legend_param = list(title="Avg. Exp.", direction="horizontal"))
draw(ht, heatmap_legend_side ="bottom")
dev.off()

########## I: Exhaustion Score vs. Pan-Cancer ##########
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

