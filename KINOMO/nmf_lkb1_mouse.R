#!/usr/bin/env Rscript

## NMF


## NMF

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)
#install.packages("NMF")
library(NMF)
library(RColorBrewer)
#pat <- 'PA060'

library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

#NMF global clustering all factors # top100



#########################lkb1_mouse-rt###################################


#LAC3
LAC3_rank_4 <- readRDS(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC3/LAC3_KINOMO_nmf_rank_4.rds")
LAC3_top100 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC3/LAC3_KINOMO_nmf_rank_4_top100_W.csv")
dim(LAC3_top100)

LAC3_top100
rownames(LAC3_top100) <- LAC3_top100[,1]
LAC3_top100 <- LAC3_top100[,-1]
rownames(LAC3_top100)

LAC3_top100

LAC3_genelist_top100 <- c(LAC3_top100[,1],LAC3_top100[,2],LAC3_top100[,3],LAC3_top100[,4])
LAC3_genelist_top100 <- LAC3_genelist_top100[!duplicated(LAC3_genelist_top100)]
LAC3_genelist_top100
saveRDS(LAC3_genelist_top100,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC3/LAC3_genelist_top100_rank4_unique.rds")

LAC3_top100_rank4 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC3/LAC3_KINOMO_rank_4_W.csv")
rownames(LAC3_top100_rank4) <- LAC3_top100_rank4[,1]
LAC3_top100_rank4 <- LAC3_top100_rank4[,-1]
LAC3_top100_rank4_unique <- LAC3_top100_rank4[LAC3_genelist_top100,]
LAC3_top100_rank4_unique
colnames(LAC3_top100_rank4_unique)<-c('LAC3_R4_F1','LAC3_R4_F2','LAC3_R4_F3','LAC3_R4_F4')
head(LAC3_top100_rank4_unique)

saveRDS(LAC3_top100_rank4_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC3/LAC3_top100_rank4_unique.rds")
write.csv(LAC3_top100_rank4_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC3/LAC3_top100_rank4_unique.csv")


#LAC11
LAC11_rank_3 <- readRDS(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC11/LAC11_KINOMO_nmf_rank_3.rds")
LAC11_top100 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC11/LAC11_KINOMO_nmf_rank_3_top100_W.csv")
dim(LAC11_top100)

LAC11_top100
rownames(LAC11_top100) <- LAC11_top100[,1]
LAC11_top100 <- LAC11_top100[,-1]
rownames(LAC11_top100)

LAC11_top100

LAC11_genelist_top100 <- c(LAC11_top100[,1],LAC11_top100[,2],LAC11_top100[,3])#,LAC11_top100[,4])
LAC11_genelist_top100 <- LAC11_genelist_top100[!duplicated(LAC11_genelist_top100)]
LAC11_genelist_top100
saveRDS(LAC11_genelist_top100,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC11/LAC11_genelist_top100_rank3_unique.rds")

LAC11_top100_rank3 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC11/LAC11_KINOMO_rank_3_W.csv")
rownames(LAC11_top100_rank3) <- LAC11_top100_rank3[,1]
LAC11_top100_rank3 <- LAC11_top100_rank3[,-1]
LAC11_top100_rank3_unique <- LAC11_top100_rank3[LAC11_genelist_top100,]
LAC11_top100_rank3_unique
colnames(LAC11_top100_rank3_unique)<-c('LAC11_R3_F1','LAC11_R3_F2','LAC11_R3_F3')#,'LAC11_R4_F4')
head(LAC11_top100_rank3_unique)

saveRDS(LAC11_top100_rank3_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC11/LAC11_top100_rank3_unique.rds")
write.csv(LAC11_top100_rank3_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC11/LAC11_top100_rank3_unique.csv")


#LAC2
LAC2_rank_5 <- readRDS(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC2/LAC2_KINOMO_nmf_rank_5.rds")
LAC2_top100 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC2/LAC2_KINOMO_nmf_rank_5_top100_W.csv")
dim(LAC2_top100)

LAC2_top100
rownames(LAC2_top100) <- LAC2_top100[,1]
LAC2_top100 <- LAC2_top100[,-1]
rownames(LAC2_top100)

LAC2_top100

LAC2_genelist_top100 <- c(LAC2_top100[,1],LAC2_top100[,2],LAC2_top100[,3],LAC2_top100[,4],LAC2_top100[,5])
LAC2_genelist_top100 <- LAC2_genelist_top100[!duplicated(LAC2_genelist_top100)]
LAC2_genelist_top100
saveRDS(LAC2_genelist_top100,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC2/LAC2_genelist_top100_rank5_unique.rds")

LAC2_top100_rank5 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC2/LAC2_KINOMO_rank_5_W.csv")
rownames(LAC2_top100_rank5) <- LAC2_top100_rank5[,1]
LAC2_top100_rank5 <- LAC2_top100_rank5[,-1]
LAC2_top100_rank5_unique <- LAC2_top100_rank5[LAC2_genelist_top100,]
LAC2_top100_rank5_unique
colnames(LAC2_top100_rank5_unique)<-c('LAC2_R5_F1','LAC2_R5_F2','LAC2_R5_F3','LAC2_R5_F4','LAC2_R5_F5')
head(LAC2_top100_rank5_unique)

saveRDS(LAC2_top100_rank5_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC2/LAC2_top100_rank5_unique.rds")
write.csv(LAC2_top100_rank5_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC2/LAC2_top100_rank5_unique.csv")


#LAC1
LAC1_rank_6 <- readRDS(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC1/LAC1_KINOMO_nmf_rank_6.rds")
LAC1_top100 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC1/LAC1_KINOMO_nmf_rank_6_top100_W.csv")
dim(LAC1_top100)

LAC1_top100
rownames(LAC1_top100) <- LAC1_top100[,1]
LAC1_top100 <- LAC1_top100[,-1]
rownames(LAC1_top100)

LAC1_top100

LAC1_genelist_top100 <- c(LAC1_top100[,1],LAC1_top100[,2],LAC1_top100[,3],LAC1_top100[,4],LAC1_top100[,5],LAC1_top100[,6])
LAC1_genelist_top100 <- LAC1_genelist_top100[!duplicated(LAC1_genelist_top100)]
LAC1_genelist_top100
saveRDS(LAC1_genelist_top100,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC1/LAC1_genelist_top100_rank6_unique.rds")

LAC1_top100_rank6 <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC1/LAC1_KINOMO_rank_6_W.csv")
rownames(LAC1_top100_rank6) <- LAC1_top100_rank6[,1]
LAC1_top100_rank6 <- LAC1_top100_rank6[,-1]
LAC1_top100_rank6_unique <- LAC1_top100_rank6[LAC1_genelist_top100,]
LAC1_top100_rank6_unique
colnames(LAC1_top100_rank6_unique)<-c('LAC1_R6_F1','LAC1_R6_F2','LAC1_R6_F3','LAC1_R6_F4','LAC1_R6_F5','LAC1_R6_F6')
head(LAC1_top100_rank6_unique)

saveRDS(LAC1_top100_rank6_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC1/LAC1_top100_rank6_unique.rds")
write.csv(LAC1_top100_rank6_unique,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/LAC1/LAC1_top100_rank6_unique.csv")





#Integrate factors

library("dplyr")
library("plyr")
library("readr")
library("purrr")
library(ComplexHeatmap)
library(circlize)
data_all <- list.files(path = ".", pattern = "*.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "Genes")
data_all[is.na(data_all)] <- 0

write.csv(data_all,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_mouse_tumor_gene_top100_correlation.csv")


set.seed(123)

#NSCLC_top100_genes_unique_updated <- read.csv(file="NSCLC_best_rank_top30_genes_unique_most_upd_2.csv")
#lkb1_mouse_top100_genes_unique_updated <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_mouse_tumor_gene_top100_correlation.csv")
lkb1_mouse_top100_genes_unique_updated <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_mouse_tumor_gene_top100_correlation_1.csv")
row.names(lkb1_mouse_top100_genes_unique_updated) <- lkb1_mouse_top100_genes_unique_updated[,1]
lkb1_mouse_top100_genes_unique_updated <-lkb1_mouse_top100_genes_unique_updated[,-1]
#rownames(NSCLC_top200_genes_unique)
lkb1_mouse_top100_genes_unique_updated[1:5,1:5]

colnames(lkb1_mouse_top100_genes_unique_updated)
mydata.cor = cor(lkb1_mouse_top100_genes_unique_updated, method = c("pearson"))
#mydata.cor = cor(lkb1_mouse_top100_genes_unique_updated, method = c("kendal"))

#lkb1_mouse_metadata <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/LKB1_metadata.csv")
lkb1_mouse_metadata <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/LKB1_metadata_1.csv")
#row.names(lkb1_mouse_metadata) <- lkb1_mouse_metadata[,1]
#lkb1_mouse_metadata <-lkb1_mouse_metadata[,-1]
#rownames(NSCLC_top200_genes_unique)
lkb1_mouse_metadata[1:5,]
lkb1_mouse_metadata <- as.data.frame(lkb1_mouse_metadata)
ann <- data.frame(lkb1_mouse_metadata$ID)
colnames(ann) <- c('ID')
table(ann$ID)
colours <- list('ID' = c('cGASKO-iso' = 'limegreen', 'cGASKO-PD1' = 'gold','iso' = 'red2', 'MCAK-iso' = 'royalblue',
                         'MCAK-PD1' = 'purple', 'PD1' = 'pink','UNK' = 'grey'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

#pdf("lkb1_mouse_kinomo_heatmap_v1.pdf",height = 20,width = 20)
#pdf("lkb1_mouse_kinomo_heatmap_v2.pdf",height = 20,width = 20)
pdf("lkb1_mouse_kinomo_heatmap_v3.pdf",height = 20,width = 20)
#pdf("testv3.2.pdf",height = 30,width = 30)
#pa = cluster::pam(mydata.cor, k = 6)
#pa = cluster::pam(mydata.cor, k = 7)
#pa = cluster::pam(mydata.cor, k =9)
pa = cluster::pam(mydata.cor, k =11)
#split <- paste0("MP", pa$clustering)
#split <- factor(paste0("MP", pa$clustering), levels=c("MP1","MP2","MP3","MP4","MP5","MP6"))
#split <- factor(paste0("MP", pa$clustering), levels=c("MP1","MP2","MP3","MP4","MP5","MP6","MP7"))
#split <- factor(paste0("MP", pa$clustering), levels=c("MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9"))
# split <- factor(paste0("MP", pa$clustering), levels=c("MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9","MP10","MP11","MP12",
#                                                       "MP13","MP14","MP15"))
split <- factor(paste0("MP", pa$clustering), levels=c("MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9","MP10","MP11"))
#split <- pa$clustering
#Heatmap(mydata.cor, name = "Corr", row_split = paste0("MP", pa$clustering),column_split = paste0("MP", pa$clustering))
#heatmap(mydata.cor, cluster_rows = T, cluster_cols = T, show_rownames = TRUE,   main = 'Heatmap',row_km = 2,column_km = 2)
#Heatmap(mydata.cor, name = "mat", column_km = 3, row_km = 3)
Heatmap(
  mydata.cor,
 # col = colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9, 1), c("blue", "white","white","white","white", "red", "red", "red", "red","red","red")),
  col=custom_magma,
  name = "Corr", #row_split = factor(paste0("MP", pa$clustering)),factor(column_split = paste0("MP", pa$clustering)),
  row_split=split, column_split = split,cluster_row_slices = FALSE,cluster_column_slices = FALSE,
  top_annotation=colAnn)
# decorate_heatmap_body("Corr", {
#   grid.text("outlier", 1.5/10, 2.5/4, default.units = "npc")
#   grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 2))
# }, slice = 2)
dev.off()

write.csv(pa$clustering,file="clustering.csv")
#write.csv(pa$clustering,file="clustering_1.csv")
write.csv(pa$clustering,file="clustering_2.csv")
#stacked bar plots

library(ggplot2)

#lkb1_mouse_metaprograms <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_metaprograms.csv")
#lkb1_mouse_metaprograms <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_metaprograms_1.csv")
lkb1_mouse_metaprograms <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_metaprograms_2.csv")
#row.names(lkb1_mouse_metadata) <- lkb1_mouse_metadata[,1]
#lkb1_mouse_metadata <-lkb1_mouse_metadata[,-1]
#rownames(NSCLC_top200_genes_unique)
#lkb1_mouse_metadata[1:5,]
lkb1_mouse_metaprograms <- as.data.frame(lkb1_mouse_metaprograms)

# # create a dataset
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# data <- data.frame(specie,condition,value)

# Grouped
#pdf("stacked.plots.kinomo.pdf")
#pdf("~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/stacked.plots.kinomo.1.pdf",height=10,width=10)
#pdf("~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/stacked.plots.kinomo.2.pdf",height=7,width=7)
pdf("~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/stacked.plots.kinomo.3.pdf",height=5,width=6)
split <- factor(paste0("MP", pa$clustering), levels=c("MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9","MP10","MP11"))
ggplot(lkb1_mouse_metaprograms, aes(fill=Sample_Code,  x=split)) + 
  geom_bar(position="stack", stat="count")+theme_bw()
ggplot(lkb1_mouse_metaprograms, aes(fill=ID,  x=split)) + 
  geom_bar(position="stack", stat="count")+theme_bw()
# ggplot(lkb1_mouse_metaprograms, aes(fill=EGFR_Status, y=EGFR_Status, x=MP)) + 
#   geom_bar(position="stack", stat="identity")+theme_bw()
# ggplot(lkb1_mouse_metaprograms, aes(fill=Percent_Tumor_PDL1_Expression, y=Percent_Tumor_PDL1_Expression, x=MP)) + 
#   geom_bar(position="stack", stat="identity")+theme_bw()
# ggplot(lkb1_mouse_metaprograms, aes(fill=Path_Response_Classification, y=Path_Response_Classification, x=MP)) + 
#   geom_bar(position="stack", stat="identity")+theme_bw()
# ggplot(lkb1_mouse_metaprograms, aes(fill=Most_Recent_Disease_Status, y=Most_Recent_Disease_Status, x=MP)) + 
#   geom_bar(position="stack", stat="identity")+theme_bw()
# ggplot(lkb1_mouse_metaprograms, aes(fill=Recurrence, y=Recurrence, x=MP)) + 
#   geom_bar(position="stack", stat="identity")+theme_bw()

dev.off()

#write.csv(lkb1_mouse_top100_genes_unique_updated,file='lkb1_mouse_tumor_gene_top100_correlation_transform_3.csv')
#lkb1_mouse_tumor_gene_top100_correlation <- read.csv(file="lkb1_mouse_tumor_gene_top100_correlation.csv")
#lkb1_mouse_tumor_gene_top100_correlation <- read.csv(file="lkb1_mouse_tumor_gene_top100_correlation_1.csv")
#lkb1_mouse_tumor_gene_top100_correlation <- read.csv(file="lkb1_mouse_tumor_gene_top100_correlation_transform_1.csv")
#lkb1_mouse_tumor_gene_top100_correlation <- read.csv(file="lkb1_mouse_tumor_gene_top100_correlation_transform_3.csv")
lkb1_mouse_tumor_gene_top100_correlation <- read.csv(file="lkb1_mouse_tumor_gene_top100_correlation_transform_4.csv")
lkb1_mouse_tumor_gene_top100_correlation <- as.data.frame(lkb1_mouse_tumor_gene_top100_correlation)
lkb1_mouse_tumor_gene_top100_correlation <- t(lkb1_mouse_tumor_gene_top100_correlation)
#write.csv(lkb1_mouse_tumor_gene_top100_correlation,file="lkb1_mouse_tumor_gene_top100_correlation_transform.csv")
#write.csv(lkb1_mouse_tumor_gene_top100_correlation,file="lkb1_mouse_tumor_gene_top100_correlation_transform_1.csv")
#write.csv(lkb1_mouse_tumor_gene_top100_correlation,file="lkb1_mouse_tumor_gene_top100_correlation_transform_2.csv")
#write.csv(lkb1_mouse_tumor_gene_top100_correlation,file="lkb1_mouse_tumor_gene_top100_correlation_transform_4.csv")
write.csv(lkb1_mouse_tumor_gene_top100_correlation,file="lkb1_mouse_tumor_gene_top100_correlation_transform_5.csv")

# meta-program hierarchies

lkb1_mouse_metaprograms <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_mouse_kinomo_enrichment_v1.csv")
row.names(lkb1_mouse_metaprograms) <- lkb1_mouse_metaprograms[,1]
lkb1_mouse_metaprograms <-lkb1_mouse_metaprograms[,-1]
#rownames(NSCLC_top200_genes_unique)
#lkb1_mouse_metadata[1:5,]
lkb1_mouse_metaprograms <- as.data.frame(lkb1_mouse_metaprograms)
lkb1_mouse_metaprograms <- t(lkb1_mouse_metaprograms)
dist_mat <- dist(lkb1_mouse_metaprograms, method = 'euclidean')

hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

cut_avg <- cutree(hclust_avg, k = 4)

rect.hclust(hclust_avg , k = 4, border = 2:6)
abline(h = 4, col = 'red')
plot(cophenetic(hclust_avg) ~ dist_mat)
cor(cophenetic(hclust_avg), dist_mat)

#NSCLC metaprograms hierarchy
lkb1_mouse_metaprograms <- read.csv(file="~/Documents/Izar_Group/NSCLC/nmf_genes_selected_new_new_color_gradient_v2_nsclc.csv")
row.names(lkb1_mouse_metaprograms) <- lkb1_mouse_metaprograms[,1]
lkb1_mouse_metaprograms <-lkb1_mouse_metaprograms[,-1]
#rownames(NSCLC_top200_genes_unique)
#lkb1_mouse_metadata[1:5,]
lkb1_mouse_metaprograms <- as.data.frame(lkb1_mouse_metaprograms)
lkb1_mouse_metaprograms <- t(lkb1_mouse_metaprograms)
dist_mat <- dist(lkb1_mouse_metaprograms, method = 'euclidean')

hclust_avg <- hclust(dist_mat, method = 'average')
pdf("NSCLC.KINOMO.MPs.Hierarchy.pdf")
plot(hclust_avg)

cut_avg <- cutree(hclust_avg, k = 5)

rect.hclust(hclust_avg , k = 5, border = 2:6)
abline(h = 5, col = 'red')
plot(cophenetic(hclust_avg) ~ dist_mat)
cor(cophenetic(hclust_avg), dist_mat)
dev.off()

# install.packages('dendextend')
# library(dendextend)
# avg_dend_obj <- as.dendrogram(hclust_avg)
# avg_col_dend <- color_branches(avg_dend_obj, h = 4)
# plot(avg_col_dend)

# library(dplyr)
# seeds_df_cl <- mutate(seeds_df, cluster = cut_avg)
# count(seeds_df_cl,cluster)

# mydata.cor.1 <- mydata.cor
# mydata.cor.1[mydata.cor.1 < 0.6] <- 0
# 
# 
# mydata.cor.2 <- mydata.cor
# 
# cor_mat <- mydata.cor.2
# cor_mat[!lower.tri(cor_mat)] <- NA # remove diagonal and redundant values
# data.frame( cor_mat) %>%
#               rownames_to_column() %>%
#               gather(key="variable", value="correlation", -rowname) %>%
#               filter(abs(correlation) > 0.6)
# 
# any_over_90 <- function(x) any(x > .8, na.rm = TRUE)
# cor_mat_1 <- focus_if(cor_mat,any_over_90, mirror = TRUE)
# 
# #library(corrplot)
# write.csv(mydata.cor , file= "NSCLC_best_rank_top30_genes_unique_most_upd_2.csv")
# #mydata.cor<-read.csv(file= "NSCLC_best_rank_top30_genes_unique_most_upd_2.csv")
# #mydata.cor[is.na(mydata.cor)] <- 0
# # mydata.cor<-read.csv(file= "NSCLC_best_rank_top30_genes_unique_most_upd_2.csv")
# row.names(mydata.cor) <- mydata.cor[,1]
# mydata.cor <-mydata.cor[,-1]
palette = colorRampPalette(c("blue", "white", "red")) (50)
#palette = colorRampPalette(c("white", "red")) (10)

test <- read.csv(file="~/Documents/Ben_Izar_Project/NMF/test.csv")
test <- as.data.frame(test)
test[,2]

test.cin70 <- read.csv(file="~/Documents/Ben_Izar_Project/NMF/CIN70.csv")
test.cin70 <- as.data.frame(test.cin70)
#test[,2]
colnames(test.cin70)
library(dplyr)
library(ComplexHeatmap)

test.cin70.df <- test.cin70%>%
  group_by(Barcode) %>%
  summarise_at(vars(CIN70_Signature), list(name = mean))
write.csv(test.cin70.df,file="test.cin70.df.csv")
#install.packages("ggplots")
#library(ggplots)
#pdf("NSCLC_top100_genes_unique_best_rank_corelation_1.pdf",height=30,width=30)
#pdf("NSCLC_top100_genes_unique_best_rank_corelation_2.pdf",height=30,width=30)
#pdf("NSCLC_top100_genes_unique_best_rank_corelation_2_stk_vs_nonstk.pdf",height=30,width=30)
#pdf("NSCLC_top100_genes_unique_best_rank_corelation_2_prim_stk_vs_bm_nonstk.pdf",height=30,width=30)
#pdf("NSCLC_top100_genes_unique_best_rank_corelation_2_prim_vs_bm_CIN70.pdf",height=30,width=30)
#pdf("test.2.pdf",height=30,width=30)
pdf("test.v3.pdf",height=30,width=30)
#pdf("NSCLC_top30_genes_unique_best_rank_corelation_gt_0.6_1.pdf",height=30,width=30)
#pdf("NSCLC_top30_genes_unique_best_rank_corelation_gt_0.6_2.pdf",height=30,width=30)
#pdf("NSCLC_top30_genes_unique_best_rank_corelation_gt_0.4.pdf",height=30,width=30)
#pdf("NSCLC_top30_genes_unique_best_rank_corelation_gt_0.9.pdf",height=30,width=30)

#pdf("NSCLC_top100_genes_unique_2nd_best_rank_corelation_1.pdf",height=30,width=30)
par(mar=c(0,0,0,0)+0.1)
#heatmap(x = mydata.cor.1, col = palette, symm = TRUE, labRow=rownames(mydata.cor.1), labCol=colnames(mydata.cor.1),margins=c(10,10))

# heatmap(x = mydata.cor, col = palette, symm = TRUE, labRow=rownames(mydata.cor), labCol=colnames(mydata.cor),margins=c(10,10),
#         RowSideColors = test[,2])
heatmap(x = mydata.cor, col = palette, symm = TRUE, labRow = FALSE, labCol = FALSE, margins=c(10,10))#,
       # RowSideColors = test[,2])
# heatmap.2(x = mydata.cor, col = palette, dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none',margins=c(10,10),
#           RowSideColors = test[,5])
#Heatmap(mydata.cor,  left_annotation = test$CIN70)
# heatmap.2(x = mydata.cor, col = palette, dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none',margins=c(10,10),
#Heatmap(mat, name = "mat", left_annotation = test$CIN70)

# legend("topright",      # location of the legend on the heatmap plot
#        #legend = c('Non-STK11-mut (BRAIN_METS)','Non-STK11-mut (PRIMARY)','STK11-mut (BRAIN_METS)','STK11-mut (PRIMARY)'), # category labels
#        #col = c('#00CED1','#696969','#FF4500','#4B0082'),  # color key
#        legend = c('BRAIN_METS','PRIMARY'), # category labels
#        col = c('#FF0000','#0000FF'),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
#)

#heatmap(x = mydata.cor, col = palette, symm = TRUE, labRow=rownames(mydata.cor), labCol=colnames(mydata.cor),margins=c(10,10))
dev.off()

colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')

colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colPrimvsBMmainTweakBM <- c('BRAIN_METS'='#FFF0F0','PRIMARY'='#0000FF')
colPrimvsBMmainTweakPRIMARY<- c('BRAIN_METS'='#FF0000','PRIMARY'='#EBF7FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colSTKvsNonSTKTweakSTK <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#F1F1F1')
colSTKvsNonSTKTweakNonSTK <- c('Non-STK11-mut'='#EBFFE3', 'STK11-mut'='#000000')

heatmap.3(x = mydata.cor, col = palette, dendrogram="both", #margins=c(6,12),
          symm = TRUE, labRow=rownames(mydata.cor), labCol=colnames(mydata.cor),margins=c(10,10),
          RowSideColors = rsc)

legend("topright",      # location of the legend on the heatmap plot
       legend = c("MBM", "MPM"), # category labels
       col = c("#A80D11", "#008DB8"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)




install.packages('RCurl')
library('RCurl') 
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
cell_cycle_genes <- read.csv(file="cell_cycle_genes.csv")

GeneListConvert <- function(gene.list, species, from, to) {
  attr.list <- list('ensg' = 'ensembl_gene_id', 'gn' = 'hgnc_symbol', 'entrez' = 'entrezgene_id')
  # get mart for appropriate species
  if (species == 'human') {
    ensembl.mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if (species == 'mouse') {
    ensembl.mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    attr.list['gn'] <- 'mgi_symbol'
  } else {
    print(paste('Species', species, 'not supported. Exiting'))
    return()
  }
  # get name.map based on from and to
  attr.set <- unlist(attr.list[c(from, to)])
  name.map <- biomaRt::getBM(attributes = attr.set, filters = attr.set[1],
                             values = gene.list, mart = ensembl.mart)
  # filter and return
  ret.list <- name.map[,2]
  ret.list <- ret.list[ which(ret.list != '') ]
  ret.list <- ret.list[ which(!is.na(ret.list)) ]
  return(ret.list)
}

#BiocManager::install("biomaRt")
#BiocManager::install("AnnotationDbi")

mat<-NSCLC_top100_genes_unique_updated[intersect(rownames(NSCLC_top100_genes_unique_updated),cell_cycle_genes[,4]),]
any(is.na(mat))
dim(mat)
palette = colorRampPalette(c("blue", "white", "red")) (50)

bk <- c(-100,seq(0,100,by=20))
#colors (one less than breaks
mycols <- c("white",colorRampPalette(colors = c("white","red"))(length(bk)-2))

#library('grid')
library(ComplexHeatmap)
#cell_cycle_genes_gn<-GeneListConvert(cell_cycle_genes[,2],'human','ensg','gn')
#library('RColorBrewer')
#pdf("NSCLC_top300_genes_unique_2nd_best_rank_heatmap_cellcycle_1.pdf",height=30,width=30)
#pdf("NSCLC_top300_genes_unique_best_rank_heatmap_cellcycle_1.pdf",height=30,width=30)
pdf("NSCLC_top100_genes_unique_best_rank_heatmap_cellcycle_02Dec2021.pdf",height=30,width=30)
par(mar=c(0,0,0,0)+0.1)
#heatmap(x = as.matrix(mat), col = palette, labRow=rownames(mat), labCol=colnames(mat),margins=c(10,10))
#heatmap(x = as.matrix(mat), col = brewer.pal(11,"RdBu"), labRow=rownames(mat), labCol=colnames(mat),margins=c(10,10))

#heatmap(x = as.matrix(mat), col = mycols, breaks=bk, scale="none", cexRow = 1.5, labRow=rownames(mat), 
        
heatmap(x = as.matrix(mydata.cor), col = mycols, breaks=bk, scale="none", cexRow = 1.5, labRow=rownames(mydata.cor), 
labCol=colnames(mydata.cor),margins=c(10,10))

dev.off()


library(circlize)
library(pheatmap)

#cell_cycle_genes <- read.csv(file="cell_cycle_genes_G2M.csv")
cell_cycle_genes <- read.csv(file="cell_cycle_genes_S.csv")
mat<-NSCLC_top100_genes_unique_updated[intersect(rownames(NSCLC_top100_genes_unique_updated),cell_cycle_genes[,4]),]
any(is.na(mat))
dim(mat)

AT2_signatures <- c('SFTPB','SFTPC','SFTPD','MUC1','ETV5')
Club_signatures <- c('CYP2F2','SCGB3A2','CCKAR')
CIN_signature <- c('PELI2','BMP2','SHH','TNS4','RAB3B','ROBO1','ARHGAP28','CHN2','CST1','F13A1','CPVL','SEMA6D','NHSL2','GTF2IP7','DPYSL3','PCDH7','KHDRBS3','TRAC','TMEM156','CST4','CD24','FGF5','NTN4')
NFKB_signature <- c('PPARG','DDIT3','NUPR1','RAB3B','IGFBP4','LRRC8C','TCP11L2','MAFK','NRG1','F2R','KRT19','CTGF','ZFC3H1')
ISG_signature <- c('CCL5','CXCL10','IFNA1','IFNB1','ISG15','IFI27L2','SAMD9L')
Yamanaka_signature <- c('POU5F1', 'SOX2', 'KLF4', 'MYC', 'c-Myc','OCT3', 'OCT4', 'OTF-3', 'OTF3', 'OTF4', 'Oct-3', 'Oct-4')
FB1_0.4_signature_v2 <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/NMF/signatures/FB1_0.4_signature_v2.csv")
FB2_0.4_signature_v2 <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/NMF/signatures/FB2_0.4_signature_v2.csv")
FB3_0.4_signature_v2 <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/NMF/signatures/FB3_0.4_signature_v2.csv")
FB4_0.4_signature_v2 <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/NMF/signatures/FB4_0.4_signature_v2.csv")
FB5_0.4_signature_v2 <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/NMF/signatures/FB5_0.4_signature_v2.csv")
FB10_0.4_signature_v2 <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/NMF/signatures/FB10_0.4_signature_v2.csv")

NSCLC_top100_genes_unique_updated <- read.csv(file="NSCLC_best_rank_top300_genes_unique_most_updated.csv")
row.names(NSCLC_top100_genes_unique_updated) <- NSCLC_top100_genes_unique_updated[,1]
NSCLC_top100_genes_unique_updated <-NSCLC_top100_genes_unique_updated[,-1]
#rownames(NSCLC_top200_genes_unique)
NSCLC_top100_genes_unique_updated[1:5,1:5]

mat<-NSCLC_top100_genes_unique_updated[intersect(rownames(NSCLC_top100_genes_unique_updated), FB10_0.4_signature_v2[,1]),]
any(is.na(mat))
dim(mat)

sample_categories <- read.csv(file="sample_categories.csv")

#col_fun = colorRamp2(c(-100, 0, 100), c("white", "red"))
#col_fun(seq(-20, 10))


bk <- c(-100,seq(0,100,by=20))
#colors (one less than breaks
mycols <- c("white",colorRampPalette(colors = c("white","red"))(length(bk)-2))

ann <- data.frame(sample_categories$PRIM_vs_BM,sample_categories$STK11mut_vs_NonSTK11mut)
colnames(ann) <- c('PRIM_vs_BM','STK11mut_vs_NonSTK11mut')
colours <- list(bar = c('PRIMARY' = 'red', 'BRAIN_METS' = 'blue','STK11_mut' = 'green', 'Non_STK11_mut' = 'yellow', 'STK11_WT' = 'black'))
colAnn <- HeatmapAnnotation(df = ann,
  which = 'col',
  col = list('PRIM_vs_BM'= c('PRIMARY' = 'red', 'BRAIN_METS' = 'blue'),'STK11mut_vs_NonSTK11mut'= c('STK11_mut' = 'green', 'Non_STK11_mut' = 'yellow', 'STK11_WT' = 'black')))#,
  #annotation_width = unit(c(1, 4), 'cm'),
  #gap = unit(1, 'mm'))

#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_cellcycle_G2M_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_cellcycle_S_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_CIN_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_NFKB_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_AT2_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_Club_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_ISG_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_FB1_0.4_signature_v2_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_FB2_0.4_signature_v2_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_FB3_0.4_signature_v2_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_FB4_0.4_signature_v2_02Dec2021.pdf",height=20,width=20)
#pdf("NSCLC_top100_genes_unique_best_rank_heatmap_FB5_0.4_signature_v2_02Dec2021.pdf",height=20,width=20)
pdf("NSCLC_top100_genes_unique_best_rank_heatmap_FB10_0.4_signature_v2_09Dec2021.pdf",height=20,width=20)
par(mar=c(0,0,0,0)+0.1)
Heatmap(
as.matrix(mat),
name = "Metagene Expr",
#column_title = "Cell cycle gene signature: G2/M",
#column_title = "Cell cycle gene signature: S",
#column_title = "CIN signature",
#column_title = "NFKB_signature",
#column_title = "AT2_signature",
#column_title = "Club_signature",
#column_title = "ISG_signature",
#column_title = "FB1_0.4_signature",
#column_title = "FB2_0.4_signature",
#column_title = "FB3_0.4_signature",
#column_title = "FB4_0.4_signature",
#column_title = "FB5_0.4_signature",
column_title = "FB10_0.4_signature",
col = mycols, 
show_row_names = TRUE,
show_column_names = TRUE,
cluster_rows = TRUE,
cluster_columns = TRUE,
show_column_dend = TRUE,
show_row_dend = TRUE,
row_dend_reorder = TRUE,
column_dend_reorder = TRUE,
column_names_gp = grid::gpar(fontsize = 7),
row_names_gp = grid::gpar(fontsize = 16),
clustering_method_rows = "ward.D2",
clustering_method_columns = "ward.D2", 
top_annotation=colAnn)
dev.off()

Heatmap(mat, name = "mat", col = mycols)

pheatmap(mat)

cols = colorRampPalette(c("white", "red"))(100)

pheatmap(mat, cluster_rows = T, cluster_cols = T, show_rownames = TRUE, 
color = cols, main = 'Heatmap')

pheatmap(mat,
              annotation_colors = palette)

Heatmap(x = as.matrix(mat), col = mycols, breaks=bk, scale="none", cexRow = 1.5)# labRow=rownames(mat), 

#heatmap.2(data_matrix,col=brewer.pal(11,"RdBu"),scale="row", trace="none")


library(pheatmap)
pdf("NSCLC_top100_genes_unique_2nd_best_rank_heatmap_cellcycle_2.pdf",height=30,width=30)
par(mar=c(0,0,0,0)+0.1)
#heatmap(x = as.matrix(mat), col = palette, labRow=rownames(mat), labCol=colnames(mat),margins=c(10,10))
#heatmap(x = as.matrix(mat), col = brewer.pal(11,"RdBu"), labRow=rownames(mat), labCol=colnames(mat),margins=c(10,10))
pheatmap(as.matrix(mat), col = mycols, breaks=bk, cluster_rows = TRUE,
  cluster_cols = TRUE, show_rownames = T, show_colnames = T,fontsize_row = 20)
dev.off()

#NSCLC_top300_genes_unique_2nd_best_rank_corelation <- #read.csv(file="NSCLC_top300_genes_unique_2nd_best_rank_corelation.csv")
NSCLC_top50_genes_unique_best_rank_corelation <- read.csv(file="NSCLC_best_rank_top50_genes_unique_most_upd_1.csv")
#row.names(NSCLC_top300_genes_unique_2nd_best_rank_corelation) <- NSCLC_top300_genes_unique_2nd_best_rank_corelation[,1]
#NSCLC_top300_genes_unique_2nd_best_rank_corelation <-NSCLC_top300_genes_unique_2nd_best_rank_corelation[,-1]
#rownames(NSCLC_top200_genes_unique)
#NSCLC_top300_genes_unique_2nd_best_rank_corelation[1:5,1:5]

#colnames(mydata.cor) <- colnames(NSCLC_top300_genes_unique_2nd_best_rank_corelation)
#rownames(mydata.cor) <- rownames(NSCLC_top300_genes_unique_2nd_best_rank_corelation)

#pdf("NSCLC_top300_genes_unique_rank_corelation_1.pdf",height=30,width=30)
#par(mar=c(0,0,0,0)+0.1)
#heatmap(x = mydata.cor, col = palette, symm = TRUE, labRow=rownames(mydata.cor), labCol=colnames(mydata.cor),margins=c(10,10))
#dev.off()

pa = cluster::pam(mydata.cor, k = 2:5)

#pdf("test.pdf")
#pdf("testv2.pdf",height = 30,width = 30)

install.packages("maptree")
library(maptree)
library(fpc)
fit2 <- hclust(dist(1-mydata.cor), method="ward.D2")
pdf("test.1.pdf")
plot(fit2,cex=0.3)
#rect.hclust(fit2,h = 50, which = c(1,2), border = 3:4)
#rect.hclust(fit2,k = 7, border = "red")
dev.off()

# a <- agnes (mydata.cor, method="ward")
# b <- kgs (fit2, fit2$height, maxclust=20)
# pdf("test.2.pdf")
# plot (names (b), b, xlab="# clusters", ylab="penalty")
# dev.off()

ges.clust <- pamk(data = dist(1-mydata.cor), krange = 2:15)$pamobject
ges.clust$silinfo

pdf("test.3.pdf")
#pa = cluster::pam(mydata.cor, k = 7)
Heatmap(mydata.cor, name = "mat", row_split = paste0("gr", ges.clust$silinfo),column_split = paste0("gr", ges.clust$silinfo))
#heatmap(mydata.cor, cluster_rows = T, cluster_cols = T, show_rownames = TRUE,   main = 'Heatmap',row_km = 2,column_km = 2)
#Heatmap(mydata.cor, name = "mat", column_km = 3, row_km = 3)
dev.off()

library(mclust)
guess <- Mclust(mydata.cor)
summary(guess)

clusters <- 9
k_fit <- kmeans(mydata.cor, 
                clusters,     # number of clusters
                iter.max=100, # max steps in the fitting process
                nstart=100)   # number of random starts

library(cluster)
pdf("test.4.pdf")
Heatmap(mydata.cor, 
         k_fit$cluster, 
         color=TRUE, 
         shade=TRUE, 
         labels=2, 
         lines=0)
dev.off()

##### Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours
# Custom color palette
library(RColorBrewer)
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

ggplot(data = mydata, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

