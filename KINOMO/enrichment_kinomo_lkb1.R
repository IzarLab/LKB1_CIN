#library(rmarkdown)
# Author: Somnath Tagore, Ph.D. 
# Title: Enrichment Analysis (GSEA, Pathway)
# Script Name: enrichment.R
# Last Updated: 09/01/2022

#Instructions
#install.packages("enrichR")
library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# combine both Cancer hallmarks (ch) and GO

#tumor
#mylist <- read.csv(file="~/Documents/Izar_Group/durva/kinomo/downstream/durva_kinomo_enrichment_pathway.csv")
#mylist <- read.csv(file="~/Documents/Izar_Group/sotr_notr/KINOMO/downstream/withfilter/sotr_withfilter_kinomo_enrichment_pathway.csv")
#mylist <- read.csv(file="~/Documents/Izar_Group/sotr_notr/KINOMO/downstream/withoutfilter/sotr_withoutfilter_kinomo_enrichment_pathway.csv")
#mylist <- read.csv(file="~/Documents/Izar_Group/sotr_notr/KINOMO/downstream/withoutfilter/sotr_withoutfilter_kinomo_enrichment_pathway.csv")
#mylist <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_pathway.csv")
mylist <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_pathway_1.csv")
#mylist <- read.csv(file="~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/lkb1_kinomo_enrichment_pathway_2.csv")
mylist <- as.data.frame(mylist)
#mylist <- toupper(mylist)
# mylist <- c(
#   
# )
#tumor.mylist <- mylist

dbs <- c("MSigDB_Hallmark_2020","GO_Biological_Process_2021")
#dbs <- c("MSigDB_Hallmark_2020")
if (websiteLive) {
  #enriched <- enrichr(unique(mylist[1:100,1]), dbs)
  #enriched <- enrichr(mylist[,1], dbs)
  #enriched <- enrichr(mylist, dbs)
  enriched <- enrichr(mylist$MP7, dbs)
}

#mut_enr<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))

mut_enr_ch<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))
mut_enr_go<-mutate(enriched[[2]], qscore = -log(Adjusted.P.value, base=10))

dim(mut_enr_ch)
dim(mut_enr_go)
mut_enr <- rbind.data.frame(mut_enr_ch,mut_enr_go)
dim(mut_enr)

celltype = "Tumor"

cohort = "LKB1-Mouse"
enrichdbs = "Cancer_Hallmarks_and_Gene_Ontology"
# topn = "Top100"
# mp = "MP15"
topn = "All"
mp = "MP7"
#h_mut_enr1 <- mut_enr[1:100,]#[1:500,] 
h_mut_enr1 <- mut_enr
h_mut_enr <- h_mut_enr1
#df[df$col1 == 1, ]

write.csv(h_mut_enr,paste0(cohort,"_",mp,"_",celltype,"_Enrichment_",enrichdbs,"_",topn,"_v5.csv"))
          
ggp<- h_mut_enr %>%
  #extract(Term, into = c("Term", "id"), "(.*)\\s\\((.*)\\)") %>%
  #mutate(Term=stri_trans_totitle(Term)) %>%
  ggplot(aes(qscore, reorder(Term, qscore), fill = P.value)) +
  scale_fill_gradient(low = "red", high = "blue") +
  geom_bar(stat = "identity") +
  
  #facet_wrap( ~ Class, nrow = 3, , scales = "free") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank()
  ) +
  geom_text(
    aes(label = paste("P.val=", round(P.value,3))),
    color = "black",
    size = 4,
    hjust = 1, nudge_x = 2
    #+ 
    # xlab("qscore") + ylab("Description")
    #position = position_dodge(width = 1.5)
  ) + theme_bw()
myggp<-ggp + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))+ggtitle(paste(cohort,"-",mp,": ",enrichdbs,"\ntop = ",topn,"_genes")) +xlab("qscore") + ylab("Description")
#pdf(paste0(cohort,"_",mp,"_",celltype,"_Enrichment_",enrichdbs,"_",topn,"_v1.pdf"), width = 10, height = 10)
#pdf(paste0(cohort,"_",mp,"_",celltype,"_Enrichment_",enrichdbs,"_",topn,"_v2.pdf"), width = 10, height = 10)
#pdf(paste0(cohort,"_",mp,"_",celltype,"_Enrichment_",enrichdbs,"_",topn,"_v3.pdf"), width = 10, height = 10)
pdf(paste0(cohort,"_",mp,"_",celltype,"_Enrichment_",enrichdbs,"_",topn,"_v5.pdf"), width = 10, height = 10)
print(myggp)
#print(mygg1p)
#print(mygg2p)
dev.off()


#LKB1-mouse data
# mouse to human gene conversion
`%nin%` = Negate(`%in%`)
#mylist.copy <- read.csv(file="LKB1.all.samples.integrated.v1.2.dge.iso.vs.mcak.iso.all.pathways.csv")
#mylist.copy <- read.csv(file="LKB1.all.samples.integrated.v1.2.dge.iso.vs.cGASKO-iso.all.pathways.csv")
mylist.copy <- read.csv(file="LKB1.all.samples.integrated.v1.2.dge.PD1.vs.MCAK-PD1.all.pathways.csv")
#mylist.copy <- read.csv(file="LKB1.all.samples.integrated.v1.2.dge.PD1.vs.cGASKO-PD1.all.pathways.csv")

mylist.copy <- as.data.frame(mylist.copy)
musGenes <- mylist.copy$gene_symbol[1:100]
musGenes <- toupper(musGenes)
#musGenes.iso.vs.mcak.iso <- musGenes
#musGenes.iso.vs.cGASKO.iso <- musGenes
#musGenes.PD1.vs.MCAK.PD1 <- musGenes
#musGenes.PD1.vs.cGASKO.PD1 <- musGenes

#musGenes <- musGenes.iso.vs.cGASKO.iso[!(musGenes.iso.vs.cGASKO.iso %in% musGenes.iso.vs.mcak.iso)]

musGenes <- musGenes.PD1.vs.cGASKO.PD1[!(musGenes.PD1.vs.cGASKO.PD1 %in% musGenes.PD1.vs.MCAK.PD1)]

# library("biomaRt")
# human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 
# genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
#                   values = musGenes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
# 
# humanx <- unique(genesV2[, 2])
# 
# convertMouseGeneList(musGenes)
# 

library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

gene_list <- musGenes
convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
        
      }
    }
  }
  
 # return (output)
}
write.csv(output,file="output.csv")
dbs <- c("MSigDB_Hallmark_2020")
#dbs <- c("MSigDB_Hallmark_2020")
if (websiteLive) {
  #enriched <- enrichr(unique(mylist[1:100,1]), dbs)
  #enriched <- enrichr(mylist[,1], dbs)
  #enriched <- enrichr(mylist, dbs)
  enriched <- enrichr(musGenes, dbs)
}

#mut_enr<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))

mut_enr_ch<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))
#mut_enr_go<-mutate(enriched[[2]], qscore = -log(Adjusted.P.value, base=10))

dim(mut_enr_ch)
dim(mut_enr_go)
mut_enr <- rbind.data.frame(mut_enr_ch)
dim(mut_enr)
celltype = "Malignant"

cohort = "LKB1-mouse"
enrichdbs = "Cancer_Hallmarks"
# topn = "Top100"
# mp = "MP15"
topn = "Top100 genes - Enriched in PD1"
#mp = "MCAK-ISO_vs_ISO"
#mp = "cGASKO-ISO_vs_ISO"
mp = "MCAK-PD1_vs_PD1"
#mp = "cGASKO-PD1_vs_PD1"
#PD1.vs.cGASKO.PD1
h_mut_enr1 <- mut_enr[1:30,]#[1:500,] 
h_mut_enr <- h_mut_enr1
#df[df$col1 == 1, ]

ggp<- h_mut_enr %>%
  #extract(Term, into = c("Term", "id"), "(.*)\\s\\((.*)\\)") %>%
  #mutate(Term=stri_trans_totitle(Term)) %>%
  ggplot(aes(qscore, reorder(Term, qscore), fill = P.value)) +
  scale_fill_gradient(low = "red", high = "blue") +
  geom_bar(stat = "identity") +
  
  #facet_wrap( ~ Class, nrow = 3, , scales = "free") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank()
  ) +
  geom_text(
    aes(label = paste("P.val=", round(P.value,3))),
    color = "black",
    size = 4,
    hjust = 1, nudge_x = 2
    #+ 
    # xlab("qscore") + ylab("Description")
    #position = position_dodge(width = 1.5)
  ) + theme_bw()
myggp<-ggp + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))+ggtitle(paste("Cohort: ",cohort,"\n","Database: ",enrichdbs,"\n","Comparison: ",mp,"\n","Enrichment: ",topn)) +xlab("qscore") + ylab("Description")
#pdf(paste0(cohort,"_",mp,"_",celltype,"_Enrichment_",enrichdbs,"_",topn,"_v1.pdf"), width = 10, height = 10)
pdf(paste0(cohort,"_",mp,"_","Enrichment_",enrichdbs,"_",topn,"_v2.pdf"), width = 10, height = 10)
print(myggp)
#print(mygg1p)
#print(mygg2p)
dev.off()


# dotplot

library(ggplot2)
#data <- read.csv("~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/S1.csv", sep =";", header = TRUE, stringsAsFactors = FALSE)
#data <- read.csv("~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/s1_lkb1_dotplot.csv", sep =",", header = TRUE, stringsAsFactors = FALSE)
data <- read.csv("~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/s1_lkb1_dotplot.v2.csv", sep =",", header = TRUE, stringsAsFactors = FALSE)

category_order <- unique(data$Pathways)

# Convert category to a factor with defined order
data$category <- factor(data$Pathways, levels = category_order)

# ggplot(data, aes(x = value, y = category)) +
#   geom_point() +
#   labs(x = "Value", y = "Category") +
#   ggtitle("Dotplot with Ordered Categories")

S1<- ggplot(data, aes(x= MPs, y=category, size=Enrichment, color=Adj.Pval, group=MPs)) + geom_point(alpha = 0.8) + 
  theme_classic()
S1 

# S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(0.000000000000000007, 
#                                                                                           0.002))
#pdf('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/s1_lkb1_dotplot.pdf',height=7,width=10)
pdf('~/Documents/Izar_Group/lkb1-cin-cgas-mouse/malignant/kinomo/downstream/s1_lkb1_dotplot.v2.pdf',height=5,width=7)
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
S1+scale_size(breaks=c(0,1,2,3),range = c(0, 3))
#S1
dev.off()


