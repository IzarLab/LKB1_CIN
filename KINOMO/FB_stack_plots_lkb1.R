
library(ggplot2)

# create a dataset
tumor_nontumor <- read.csv(file="~/Documents/Izar_Group/durva/kinomo/downstream/durva_kinomo_enrichment_v1.csv")
#tumor_nontumor <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/melanoma/NMF_melanoma/melanoma_stack_plots_1.csv")
row.names(tumor_nontumor) <- tumor_nontumor[,1]
tumor_nontumor<-tumor_nontumor[,-1]
rownames(tumor_nontumor)
tumor_nontumor<-as.data.frame(tumor_nontumor)

# FB_CIN70 <- read.csv(file="~/Documents/Ben_Izar_Project/NMF/FB_CIN70.csv")
# FB_CIN70 <- as.data.frame(FB_CIN70)
# # Color palette v2
# colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
# colPrimvsBMmainTweakBM <- c('BRAIN_METS'='#FFF0F0','PRIMARY'='#0000FF')
# colPrimvsBMmainTweakPRIMARY<- c('BRAIN_METS'='#FF0000','PRIMARY'='#EBF7FF')
# colSTKvsNonSTKmain <- c('STK11-WT'='#3CB371', 'STK11-MUT'='#000000')
# colSTKvsNonSTKTweakSTK <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#F1F1F1')
# colSTKvsNonSTKTweakNonSTK <- c('Non-STK11-mut'='#EBFFE3', 'STK11-mut'='#000000')
# 
# colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')
# 
# #colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut_BRAIN_METS'='#3CB371','Non-STK11-mut_PRIMARY'='#2F4F4F','STK11-mut_BRAIN_METS'='#FF4500','STK11-mut_PRIMARY'='#4B0082')
# 
# colPrimvsBMvsSTKvsNonSTKTweakNonSTKPrim.STKPrim.STKBM <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#F2F2F2','STK11-mut (BRAIN_METS)'='#FCF2EE','STK11-mut (PRIMARY)'='#F1E9F7')
# 
# colPrimvsBMvsSTKvsNonSTKTweakNonSTKBM.STKPrim.STKBM <- c('Non-STK11-mut (BRAIN_METS)'='#E3FDEF','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FCF2EE','STK11-mut (PRIMARY)'='#F1E9F7')
# 
# colPrimvsBMvsSTKvsNonSTKTweakSTKPrim.NonSTKPrim.NonSTKBM <- c('Non-STK11-mut (BRAIN_METS)'='#E3FDEF','Non-STK11-mut (PRIMARY)'='#F2F2F2','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#F1E9F7')
# 
# colPrimvsBMvsSTKvsNonSTKTweakSTKBM.NonSTKPrim.NonSTKBM <- c('Non-STK11-mut (BRAIN_METS)'='#E3FDEF','Non-STK11-mut (PRIMARY)'='#F2F2F2','STK11-mut (BRAIN_METS)'='#FCF2EE','STK11-mut (PRIMARY)'='#4B0082')
# 
# c
# cols = c(
#   'Primary'='#A80D11',
#   'Brain_Mets'='#008DB8')
# 
# cols = c(
#   'STK11-MUT'='#A80D11',
#   'STK11-WT'='#008DB8')
# 
# cols = c(
#   'PA001' = '#808080', 'PA004' = '#d3d3d3', 'PA005' = '#2f4f4f',
#    'PA019' = '#556b2f', 'PA025' = '#8b4513','PA034' = '#2e8b57', 'PA042' = '#228b22',
#  'PA043' = '#7f0000','PA048' = '#191970', 'PA054' = '#808000', 'PA056' = '#b8860b',
#    'PA060' = '#008b8b',
#     'PA067' = '#4682b4', 'PA068' = '#d2691e','PA070' = '#9acd32','PA072' = '#cd5c5c',
#     'PA076' = '#00008b', 'PA080' = '#32cd32','PA104' = '#8fbc8f','PA125' = '#8b008b',
#    'PA141' = '#b03060','N254' = '#ff4500',#'N561' = '#00ced1',
#  'N586' = '#ffa500',
#     'KRAS_10' = '#ffd700', 'KRAS_11' = '#6a5acd', 'KRAS_12' = '#deb887',
#   'KRAS_13' = '#00ff00', 'KRAS_17' = '#00fa9a', 'KRAS_4' = '#dc143c',
#     'KRAS_6' = '#0000ff', 'KRAS_7' = '#a020f0', 'KRAS_8' = '#adff2f',
#    'STK_1' = '#da70d6', 'STK_14' = '#ff00ff', 'STK_15' = '#1e90ff',
#   'STK_18' = '#f0e68c', 'STK_2' = '#dda0dd','STK_20' = '#90ee90',
#    'STK_21' = '#ffa07a', 'STK_22dot2' = '#87cefa', 'STK_3' = '#7fffd4',
#   'STK_5dot1' = '#ff69b4','STK_5dot1' = '#ffb6c1'
# )

library(dplyr)

# ggplot(data= iris, aes(x= Sepal.Width, fill= Species))+
#   geom_bar() +
#   geom_text(stat = "count", aes(label =..count..), position=position_stack(vjust=0.5))+
#   scale_y_continuous(limits = c(0,20))
# iris
cols <- colSTKvsNonSTKmain
cols <- colPrimvsBMmain
tumor_nontumor <- tumor_nontumor %>%
  group_by(Sample) %>%
  #  mutate(label_y = cumsum(Number) - 0.5 * Weight)
  mutate(ypos = cumsum(Number) - 0.5 * Number,
         label = paste0(Type, ":", Number, "-factor")) 

#tumor_nontumor$Sample <- factor(tumor_nontumor$Sample, levels=tumor_nontumor$Sample[order(tumor_nontumor$Sample)])
tumor_nontumor$Sample = factor(tumor_nontumor$Sample, levels = unique(tumor_nontumor$Sample))

# a <- tumor_nontumor
# b <- FB_CIN70
# 
# df <- merge(a,b,by=tumor_nontumor$Sample)

ggplot(tumor_nontumor, aes(x = factor(Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP11','MP12','MP13','MP14','MP15')), y = Number, fill = Type, group=1,label=Number)) + 
  geom_col(aes(y = Number, fill = Type)) + ggtitle("NSCLC: KINOMO Metaprograms (Sample-wise)")+
  geom_bar( stat='identity')+ geom_text(
  size = 3, position = position_stack(vjust = 0.5)) 
  #geom_line(aes(y=FB_CIN70$CIN70)) + 
  #scale_y_continuous(labels = scales::format_format(big.mark = ".", decimal.mark = ",", scientific = FALSE),
   #                  sec.axis = sec_axis(~ max(FB_CIN70$CIN70)))

ggplot(data = tumor_nontumor, aes(x = factor(Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP11','MP12','MP13','MP14','MP15')), y = Number, fill = Type, group=1,label=Number)) + 
  geom_bar( stat='identity')+ geom_text(
    size = 3, position = position_stack(vjust = 0.5))   +
  
  # geom_line( aes(y=tumor_nontumor$CIN70_11), size=2, color="Green") +
  # # Custom the Y scales:
  # scale_y_continuous(
  #   
  #   # Features of the first axis
  #   name = "Count",
  #   
  #   # Add a second axis and specify its features
  #   sec.axis = sec_axis(~./5, name="CIN70")
  # )+ #+
theme_bw()+
  scale_fill_manual(values=cols) +
  # theme(
  #   #axis.title.y = element_text(color = cols, size=12),
  #   axis.title.y.right = element_text(color = "Green", size=12)
  # ) +
  
  #ggtitle("NSCLC: KINOMO Metaprograms (Primary vs BM) and CIN70")
 # ggtitle("NSCLC: KINOMO Metaprograms (STK11-MUT vs STK11-WT) and CIN70")
  ggtitle("NSCLC: KINOMO Metaprograms (STK11-MUT vs STK11-WT) ")
#ggtitle("NSCLC: KINOMO Metaprograms (Sample-wise)")
#theme_ipsum()

# yplot + stat_compare_means()+theme(text = element_text(size = 16))+   # Add counts by group to boxplot
#   annotate("text",
#            #x = 1:length(table(ichorcna_data$Sample.Type)),
#            x = 1:length(table(tumor_nontumor$Sample)),
#            y = aggregate(Number ~ Sample, tumor_nontumor, median)[ , 2],
#            label = table(tumor_nontumor$Sample),
#            col = "White",
#            vjust = - 1)

#df$derma <- factor(df$derma, levels = df$derma[order(df$prevalence)])
#df1$col1 = factor(df1$col1)
#p <- ggplot()
ggplot(data = tumor_nontumor, aes(x = factor(Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP11','MP12','MP13','MP14','MP15')),y = Number, fill = Type)) + 
  geom_bar( stat='identity')+
  #  geom_bar(position = "fill", stat='identity')+
#  geom_line(aes(y = FB_CIN70$CIN70), color = "red") 
  theme_bw()+
  scale_fill_manual(values=cols)+
  ggtitle("NSCLC KINOMO FB-sample contributions") +
  #ggtitle("NSCLC KINOMO FB: Primary vs Brain_Mets") +
  #ggtitle("NSCLC KINOMO FB: STK11-mut vs Non-STK11-mut") +
  xlab("Samples") + 
  # ylab("Percentage")+
  ylab("Count") + theme(text = element_text(size = 20)) # +#theme(axis.text.x=element_text(size=rel(1.1))

#aes(x = factor(Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP11','MP12','MP13','MP14','MP15')),y = Number, fill = Type)) + 
level_order <- factor(tumor_nontumor$Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9',
                                                       'MP10','MP11','MP12','MP13','MP14','MP15'))

ggplot(tumor_nontumor, aes(x=level_order,  y=Number,fill=Type )) + 
#ggplot(tumor_nontumor, aes(fill=Type, y=Number, x=Sample)) + 
  geom_bar(stat="identity")+
  #geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=cols)

# df <- merge(a,b,by="x")
#   
# ggplot(df, aes(x=x, y=y.x, fill=z)) +
#     geom_bar(position="stack",stat="identity") + 
#     geom_line(aes(x=x, y=y.y)) + 
#     ylab("") + xlab("x")
# geom_text(position="stack",aes(Sample,Number,label=Type),size=5)
# geom_text(aes(label=Type), position=position_dodge(cumsum(Number) - 0.5 * Number))
# geom_text(stat = "Number", aes(label =..Number..), position=position_stack(vjust=0.5))
#ggsave("Melanoma_KINOMO_FB_sample_contributions_percent.pdf",height = 10,width = 10)
#ggsave("Melanoma_KINOMO_FB_sample_contributions_count.pdf",height = 10,width = 10)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_all_samples.pdf",height = 10,width = 10)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_all_samples_v3.pdf",height = 10,width = 15)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_stk11_vs_nonstk11.pdf",height = 10,width = 10)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_stk11_vs_nonstk11_v2.pdf",height = 10,width = 15)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_prim_vs_bm.pdf",height = 10,width = 10)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_prim_vs_bm_v2.pdf",height = 10,width = 15)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_prim_vs_bm_v3.pdf",height = 10,width = 15)
#ggsave("NSCLC_KINOMO_FB_sample_contributions_count_stk11_vs_nonstk11_v3.pdf",height = 10,width = 15)
ggsave("NSCLC_KINOMO_FB_sample_contributions_count_stk11_vs_nonstk11_v5.pdf",height = 10,width = 15)


############ selected genes
library("scales")

# nmf_selected_genes
#nmf_selected_genes<- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_v1.csv")
#nmf_selected_genes<- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_v2.csv")
#nmf_selected_genes<- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_v5.csv")
#nmf_selected_genes<- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_v6.csv")
nmf_selected_genes<- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_v2.1.csv")
#nmf_selected_genes<- read.csv(file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/nmf_genes_selected_new_new_color_gradient_v2_nsclc.csv")
row.names(nmf_selected_genes) <- nmf_selected_genes[,1]
nmf_selected_genes <-nmf_selected_genes[,-1]
#rownames(Melanoma_top200_genes_unique)
#Melanoma_FB_stouff[1:5,1:5]

vec_range3 <- rescale(as.matrix(nmf_selected_genes)) 

#vec_range3 <- nmf_selected_genes
library(pheatmap)

# save_pheatmap_pdf_2 <- function(x, filename, width=4, height=24) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }

save_pheatmap_pdf_2 <- function(x, filename, width=3, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#set.seed(1)
#m1<-matrix(c(rnorm(1000)), ncol=100)
#cols = colorRampPalette(c("white", "black"))(2)
library(RColorBrewer)
# p1<-pheatmap(nmf_selected_genes,
#          cluster_rows = F,
#          cluster_cols = F,
#          show_rownames = TRUE, 
#          color = cols,
#          breaks = c(0, 0.5,1),  # distances 0 to 3 are red, 3 to 9 black
#          main = 'KINOMO Factor Blocks \n Selected metagenes (Top100)')
#col.pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(1:100))
# redmono = c("#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FFF5F0")
# mat_breaks <- quantile_breaks(vec_range3,n = 100)
# quantile_breaks <- function(xs, n = 10){
#   breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
#   unique(breaks)
# }
# 
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
#col.pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(mat_breaks))

#pdf(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/kinomo_selected_genes_lkb1_v5.pdf",width=10, height=30)
#pdf(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/kinomo_selected_genes_lkb1_v6.pdf",width=5, height=7)
p1<-pheatmap(vec_range3,#clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan",
             #            p1<-pheatmap(nmf_selected_genes,#clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan",
             #           cutree_cols = 5,cutree_rows = 6,s
             show_rownames = T,show_colnames = T, #legend_breaks = c(0,1),
             # legend_labels = c("No","Yes"),
             color = colorRampPalette(c("grey100","grey0"))(50),#scale = "none",treeheight_row = 100,
             #color = rev(redmono),# color = colorRampPalette(c('white','black'))(2),#scale = "none",treeheight_row = 100,
            # color = col.pal,
             # color = colorRampPalette(c('white','black'))(2),#scale = "none",treeheight_row = 100,
             #ccolor = colorRampPalette(c("white", "red"))(20),# color = colorRampPalette(c('white','black'))(2),#scale = "none",treeheight_row = 100,
             cluster_cols = F,cluster_rows = F,main = 'KINOMO MPs \n Metagenes\n (LKB1-mouse)')
#cluster_cols = F,cluster_rows = F,main = 'KINOMO Factor Blocks \n Selected metagenes (Top100)')
#save_pheatmap_pdf_2(p1,file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/melanoma/NMF_melanoma/nmf_selected_genes_v1.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/melanoma/NMF_melanoma/nmf_selected_genes_v2.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/melanoma/NMF_melanoma/nmf_selected_genes_v3.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/melanoma/NMF_melanoma/nmf_selected_genes_v4.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/kinomo_selected_genes_nsclc_v1.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/kinomo_selected_genes_nsclc_v2.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/kinomo_selected_genes_nsclc_v3.1.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/kinomo_selected_genes_nsclc_v3.2.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/kinomo_selected_genes_nsclc_v3.3.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Ben_Izar_Project/NMF/Ranked_metagenes_FBs/kinomo_selected_genes_nsclc_v3.3.refined.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/kinomo_selected_genes_lkb1_v1.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/kinomo_selected_genes_lkb1_v2.pdf")
save_pheatmap_pdf_2(p1,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/kinomo_selected_genes_lkb1_v2.1.pdf")
#save_pheatmap_pdf_2(p1,file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/kinomo_selected_genes_lkb1_v5.pdf")
dev.off()

agg.CIN70 <- read.csv(file="Agg_CIN70.csv")
agg.CIN70.df <- agg.CIN70%>%
  group_by(Sample) %>%
  summarise_at(vars(Agg_CIN70), list(name = mean))
write.csv(agg.CIN70.df,file="agg.CIN70.df.csv")
