## Karan's version

library(tidyverse)
library(reshape2)
library(matrixStats)
library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(purrr)
library(DropletUtils)
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(celldex)
library(gridExtra)
library(mixtools)
library(enrichR)
library(stringr)
#library(hash)
library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)
# Load necessary libraries
library(ggplot2)
library(viridis)

#df = read.csv("ranked_interactions-cGASKO-PD1_GSEA.txt", sep = "\t")
df = read.csv("ranked_interactions-MSA-KO_GSEA.txt", sep = "\t")

df = subset(df,fdr_ligand <= 0.05)
df = subset(df,numSigI1_fdr05_receptor >= 5)
df = subset(df, log2FC_ligand >=0.01 | log2FC_ligand <= -0.01)
df$interaction <- apply(df[, c("cell_type_ligand", "cell_type_receptor")], 1, function(x) {
  paste(sort(x), collapse = ";")
})


#pdf("cGASKO-PD1_GSEA.pdf", width = 12, height = 6)
#pdf("cGASKO-PD1_GSEA_v2.pdf", width = 10, height = 10)
df_total <- NULL
#pdf("MSA-KO_GSEA_v2.pdf", width = 10, height = 10)
pdf("MSA-KO_GSEA_v2.1.pdf", width = 10, height = 10)

for (x in unique(df$interaction)){
  print(x)
  
  interaction_df = subset(df,interaction == x)
  pos_df <- subset(interaction_df, log2FC_ligand > 0 )
  pos_genes <- unique(c(pos_df$ligand,pos_df$receptor))
  
  
  # Load Enrichr
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (websiteLive) head(dbs)
  # read signature
  dbs <- c("MSigDB_Hallmark_2020","GO_Biological_Process_2021")
  term_list = c()
  if (websiteLive & length(pos_genes) > 2) {
    enriched <- enrichr(pos_genes, dbs)
    
    mut_enr_ch<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))
    mut_enr_go<-mutate(enriched[[2]], qscore = -log(Adjusted.P.value, base=10))
    dim(mut_enr_ch)
    dim(mut_enr_go)
    mut_enr <- rbind.data.frame(mut_enr_ch,mut_enr_go)
    dim(mut_enr)
    write.csv(mut_enr,file="mut_enr_metagenes.csv")
    mp = "Meta-program n"
    cohort = "Sarcoma"
    enrichdbs = "Cancer_Hallmarks_and_Gene_Ontology"
    topn = "Top100"
    pos_mut_enr <- mut_enr
    mut_enr_p <- subset(mut_enr, Adjusted.P.value < 0.05)
    
    h_mut_enr1 <- mut_enr_p[order(mut_enr_p$qscore, decreasing = TRUE), ][1:10,]#[1:500,]
    
    pos_gene_sets <- h_mut_enr1
    term_list <- c(term_list,pos_gene_sets$Term)
  }
  
  
  neg_df <- subset(interaction_df, log2FC_ligand < 0 )
  neg_genes <- unique(c(neg_df$ligand,neg_df$receptor))
  print(neg_genes)
  
  
  # Load Enrichr
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (websiteLive) head(dbs)
  # read signature
  dbs <- c("MSigDB_Hallmark_2020","GO_Biological_Process_2021")
  
  if (websiteLive & length(neg_genes) > 2) {      
    enriched <- enrichr(neg_genes, dbs)
    
    mut_enr_ch<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))
    mut_enr_go<-mutate(enriched[[2]], qscore = -log(Adjusted.P.value, base=10))
    dim(mut_enr_ch)
    dim(mut_enr_go)
    mut_enr <- rbind.data.frame(mut_enr_ch,mut_enr_go)
    dim(mut_enr)
    write.csv(mut_enr,file="mut_enr_metagenes.csv")
    mp = "Meta-program n"
    cohort = "Sarcoma"
    enrichdbs = "Cancer_Hallmarks_and_Gene_Ontology"
    topn = "Top100"
    neg_mut_enr <- mut_enr
    mut_enr_p <- subset(mut_enr, Adjusted.P.value < 0.05)
    
    h_mut_enr1 <- mut_enr_p[order(mut_enr_p$qscore, decreasing = TRUE), ][1:10,]#[1:500,]
    neg_genes_sets <- h_mut_enr1
    term_list <- c(term_list,rev(neg_genes_sets$Term))
  }
  print(term_list)
  term_list <- term_list[!duplicated(term_list)]
  neg_subset = subset(neg_mut_enr, Term %in% term_list)
  pos_subset = subset(pos_mut_enr, Term %in% term_list)
  
  # Ensure the term_list has no duplicates
  term_list <- term_list[!duplicated(term_list)]
  
  # Create the pos_subset and neg_subset dataframes
  neg_subset <- subset(neg_mut_enr, Term %in% term_list)
  pos_subset <- subset(pos_mut_enr, Term %in% term_list)
  
  # Merge the two dataframes on the Term column
  combined_df <- merge(pos_subset[, c("Term", "qscore", "Adjusted.P.value")], 
                       neg_subset[, c("Term", "qscore", "Adjusted.P.value")], 
                       by = "Term", 
                       all = TRUE, 
                       suffixes = c("_pos", "_neg"))
  # Replace NA values with 0
  combined_df[is.na(combined_df)] <- 0
  
  # Set Term as a factor with levels in the order of term_list
  combined_df$Term <- factor(combined_df$Term, levels = term_list)
  
  # Order the dataframe by the Term factor
  combined_df <- combined_df[order(combined_df$Term), ]
  combined_df <- combined_df[rev(rownames(combined_df)), ]
  
  
  # Print the combined dataframe
  print(combined_df)
  df_total <- rbind(df_total,combined_df)
  
  data_long <- data.frame(
    Term = rep(combined_df$Term, 2),
    qscore = c(combined_df$qscore_pos,-combined_df$qscore_neg),
    pvalue = c(combined_df$Adjusted.P.value_pos,combined_df$Adjusted.P.value_neg),
    group = rep(c("Positive", "Negative"), each = nrow(combined_df))
  )
  
  # Adjust the p-values: take -log(pvalue) and make negative values for the "Negative" group
  data_long$pvalue <- -log10(data_long$pvalue)
  data_long$pvalue[is.na(data_long$pvalue) | is.infinite(data_long$pvalue)] <- 0
  
  data_long$pvalue[data_long$group == "Negative"] <- -data_long$pvalue[data_long$group == "Negative"]
  
  # Plot the data
  data_long$color_group <- with(data_long, ifelse(group == "Positive", "red", "blue"))
  max_pvalue <- max(abs(data_long$pvalue))
  
  #write.csv(data_long,'data_long.csv')
  # Plot the data
  p<- ggplot(data_long, aes(fill = pvalue, y = qscore, x = Term)) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_gradientn(colors = c("#000080", "#4682B4", "#B6D0E2", "#FAA0A0","#A95C68", "#800020"),
                         values = scales::rescale(c(-max_pvalue, 0, max_pvalue)),
                         name = "-log10(P-Value)", limits = c(-max_pvalue, max_pvalue)) +
    geom_hline(yintercept = 0, color = "black", size = 1) + 
    
    ggtitle(paste(x, " interactions, GSEA MSA-KO"))+
    ylab("-log10(Adj.Pvalue)") + 
    xlab("Term") +
    #coord_flip() +  # Flip coordinates if you want horizontal bars
    theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  print(p)
  
}
dev.off()
#cGASKO-PD1_GSEA
#write.csv(df_total,'df_total.csv')
#write.csv(df_total,'df_total_cGASKO-PD1_GSEA.csv')
write.csv(df_total,'df_total_MSA_KO_GSEA.csv')

#df_total <- read.csv('MSA-KO_GSEA_v2_likelihood.csv')
df_total <- read.csv('df_total_cGASKO-PD1_GSEA.csv')
df_total <- as.data.frame(df_total)
likelihoods_pos <- 10^(df_total$qscore_pos)

likelihoods_neg <- 10^(df_total$qscore_neg)
df_total$likelihoods_pos <- likelihoods_pos
df_total$likelihoods_neg <- likelihoods_neg
#write.csv(df_total,'MSA-KO_GSEA_v2_likelihood_v2.csv')

new_min <- 0  # Desired minimum value
new_max <- 1  # Desired maximum value
old_min <- min(df_total$likelihoods_pos)
old_max <- max(df_total$likelihoods_pos)

df_total$likelihoods_pos_scaled <- (df_total$likelihoods_pos - old_min) / (old_max - old_min) * (new_max - new_min) + new_min

old_min <- min(df_total$likelihoods_neg)
old_max <- max(df_total$likelihoods_neg)

df_total$likelihoods_neg_scaled <- (df_total$likelihoods_neg - old_min) / (old_max - old_min) * (new_max - new_min) + new_min

#write.csv(df_total,'MSA-KO_GSEA_v2_likelihood_v2.csv')
write.csv(df_total,'df_total_cGASKO-PD1_GSEA_v2.csv')

# Load necessary libraries
library(ggplot2)
# Example data
#data <- read.csv('MSA-KO_GSEA_v2_likelihood_v3.csv')
#data <- read.csv('df_total_cGASKO-PD1_GSEA_v3.csv')
#data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_cGASKO-PD1_GSEA_v5_Malignant_Tcells.csv')
#data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_cGASKO-PD1_GSEA_v6_Malignant_Tcells.csv')
#data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_cGASKO-PD1_GSEA_v7_Malignant_Myeloid.csv')
data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_cGASKO-PD1_GSEA_v7_Malignant_Tcells.csv')
#data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_cGASKO-PD1_GSEA_v6_Malignant_Myeloid.csv')
#data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_KLisovsKLcGASKOiso_GSEA_v6_Malignant_Tcells.csv')
#data <- read.csv('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/contact.tracing/df_total_KLisovsKLcGASKOiso_GSEA_v6_Malignant_Myeloid.csv')
data <- as.data.frame(data)

head(data)
# data <- data.frame(
#   Pathway = c("Pathway X", "Pathway Y", "Pathway A", "Pathway B", "Pathway C", "Pathway D", "Pathway E","Pathway A", "Pathway B"),
#   p_value = c(0.001, 0.05, 0.001, 0.05, 0.01, 0.0001, 0.02,0.001, 0.05),
#   oddsRatio = c(-2.5, -1.8, 2.5, 1.8, 3.0, 2.2, 1.5,-0.5, -0.8)
# )
# Calculate Combined Score
##data$CombinedScore <- -log(data$p_value) * data$oddsRatio
# Plot
# Define colors based on Combined Score
#data$color <- ifelse(data$CombinedScore > 0, "#87CEEB", "#FF00FF")

# Plot
#pdf("MSA-KO_GSEA_v3.pdf", width = 10, height = 10)
#pdf("MSA-KO_GSEA_Malignant_Tcellsv5.pdf", width = 10, height = 10)
#pdf("MSA-KO_GSEA_Malignant_Tcellsv6.pdf", width = 10, height = 10)
#pdf("cGASKO-PD1_GSEA_v7_Malignant_Myeloid.pdf", width = 10, height = 10)
pdf("cGASKO-PD1_GSEA_v7_Malignant_Tcells.pdf", width = 10, height = 10)
#pdf("MSA-KO_GSEA_Malignant_Myeloidv6.pdf", width = 10, height = 10)
#pdf("KLisovsKLcGASKOiso_GSEA_Malignant_Tcells6.pdf", width = 10, height = 10)
#pdf("KLisovsKLcGASKOiso_GSEA_Malignant_Myeloid6.pdf", width = 10, height = 10)
#pdf("cGASKO-PD1_GSEA_v3.pdf", width = 10, height = 10)
ggplot(data, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Combined.Score)) +
#  ggplot(data, aes(x = reorder(Term, Likelihood), y = Likelihood, fill = Likelihood)) +
  geom_bar(stat = "identity") +
  #coord_flip() +
  #scale_fill_gradientn(colors = c("#000080", "#4682B4", "#A95C68", "#800020"))+
 # scale_fill_gradientn(colors = c("#000080",  "#A95C68", "#800020"))+
  scale_fill_gradient2(low = "#000080", mid = "#FFFFFF", high = "#800020", midpoint = 0) +
  
  labs(
    #title = "Likelihood",
    title = "cGASKO_PD1 vs PD1: Malignant-Tcells",
    #title = "cGASKO_PD1 vs PD1: Malignant-Myeloid",
    #title = "KL-iso vs KL-cGAS-KO-iso: Malignant-Tcells",
   # title = "KL-iso vs KL-cGAS-KO-iso: Malignant-Myeloid",
    #title = "cGASKO_PD1 vs PD1: Malignant-Myeloid",
    y = "-log(p)*oddsRatio",
    x = "Pathway"
  ) +
  theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# #cGAS-STING
# 
# # Load necessary libraries
# # Load necessary libraries
# library(ggplot2)
# library(ggsignif)
# library(ggbreak)
# 
# # Sample data (replace with your actual data)
# data <- data.frame(
#   Group = c("PT", "BM", "Group C", "Group D"),  # Replace with actual group names
#   Value = c(0.05, 0.1, 1.5, 2.5),  # Replace with actual values
#   SE = c(0.01, 0.015, 0.1, 0.15)  # Replace with actual standard errors
# )
# 
# # Create the plot
# p <- ggplot(data, aes(x = Group, y = Value)) +
#   geom_bar(stat = "identity", fill = "skyblue", color = "black") +
#   geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE), width = 0.2) +
#   theme_minimal() +
#   ylab("Your Y-axis Label") +
#   xlab("Your X-axis Label")
# 
# # Add significance comparisons
# p <- p + geom_signif(comparisons = list(c("PT", "BM")), 
#                      map_signif_level = TRUE, 
#                      y_position = c(2.7),  # Adjust y position for p-value annotation
#                      textsize = 3.5)
# 
# # Apply the axis break
# p_broken <- p + scale_y_break(c(0.2, 1.0))  # Adjust the break points as needed
# 
# # Display the plot
# print(p_broken)
# 
# 
# library(ggplot2)
# 
# # Example data frame (replace with your actual data)
# data <- data.frame(
#   Pathway = c("Pathway1", "Pathway2", "Pathway3", "Pathway4", "Pathway5",
#               "Pathway6", "Pathway7", "Pathway8", "Pathway9", "Pathway10"),
#   value = c(-3200, -1500, -1000, -1200, -1100, 500, 700, 800, 1500, 2700)
# )
# 
# # Generate the plot with a color gradient and legend
# ggplot(data, aes(x = reorder(Pathway, value), y = value, fill = value)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # Flip the axes
#   scale_fill_gradient2(low = "navy", mid = "white", high = "darkred", midpoint = 0,
#                        name = "-log(p) * oddsRatio") + # Adding the legend
#   labs(x = "", y = "-log(p) * oddsRatio", 
#        title = "GSEA on Cancer cell : T cells interactions") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         legend.position = "right") # Adjust legend position if necessary
