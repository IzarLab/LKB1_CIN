
library(ggplot2)

# create a dataset
tumor_nontumor <- read.csv(file="kinomo_enrichment_v1.csv")
row.names(tumor_nontumor) <- tumor_nontumor[,1]
tumor_nontumor<-tumor_nontumor[,-1]
rownames(tumor_nontumor)
tumor_nontumor<-as.data.frame(tumor_nontumor)


library(dplyr)

cols <- colSTKvsNonSTKmain
cols <- colPrimvsBMmain
tumor_nontumor <- tumor_nontumor %>%
  group_by(Sample) %>%
  #  mutate(label_y = cumsum(Number) - 0.5 * Weight)
  mutate(ypos = cumsum(Number) - 0.5 * Number,
         label = paste0(Type, ":", Number, "-factor")) 

tumor_nontumor$Sample = factor(tumor_nontumor$Sample, levels = unique(tumor_nontumor$Sample))

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
