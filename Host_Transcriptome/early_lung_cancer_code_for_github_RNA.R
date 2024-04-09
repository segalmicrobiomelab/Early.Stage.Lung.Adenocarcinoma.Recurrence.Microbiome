# 07 APR 2024 ; Fares Darawshy ; Segal Lab 
#  Signatures of Early Lung Cancer. RNA sequencing

#Load/Save wd ----- 
setwd() # set your wd 
#load libraries 
library(devtools)
library(fgsea)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(ggalt)
library(forcats)
library(car)
library(magrittr)
library(lubridate)
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library("vegan")
library(ade4)
library("reshape2")
library(dplyr)	
library(Matrix)
library(scales)
library(cowplot)
library(randomForest)
library(caret)
library("mlbench")
library(RCurl)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(AnnotationHub)
library(ensembldb)
library(VennDiagram)
library(ranger)
library(xgboost)
library(rmarkdown)
library(glue)
library(ggtext)

#Load RNA Seq file 
mycounts <-read.delim2("NovaSeq_Raw_Counts_Gene_Symbols.txt", sep="\t", row.names=1)

#load meadata 
early_lung_cancer_metadata <- read.csv(file="LNSmetafile_91sub.csv")
rownames(early_lung_cancer_metadata) <- early_lung_cancer_metadata$X
early_lung_cancer_metadata$sample_name <- early_lung_cancer_metadata$X
early_lung_cancer_metadata <- early_lung_cancer_metadata %>% select(-X)


#explore 
head(mycounts)
rownames(mycounts)
colnames(mycounts)

#Keep Only MetaData of RNASeq Samples
#Keep only samples with RNA ID that matches the RNA count data 
RNA.data = early_lung_cancer_metadata[early_lung_cancer_metadata$RNA_Seq_done_Final==1,]
RNA.data <- RNA.data[RNA.data$RNA_ID!= "n.a",]

#Pick up the column that have RNA ID in metadata 
#We want to match ID names in each element 
RNA.data$RNA_ID

#Order Meta Data by SampleId
RNA.data <- RNA.data[order(RNA.data$RNA_ID),]

#Order Count Data by SampleID and replace the X before each name 
colnames(mycounts) <- gsub("X", "", colnames(mycounts))
mycounts <- mycounts [, order(colnames(mycounts))]

mycounts_mod <- mycounts[,RNA.data$RNA_ID]

#Confirm Sample IDs match for Count and Meta Data
table(colnames(mycounts_mod)==as.character(RNA.data$RNA_ID))
rownames(RNA.data) <- RNA.data$RNA_ID

#confirm columns names of data and rownamnes of metadata are similar 
table(colnames(mycounts_mod)==rownames(RNA.data))


#############################Creating DESeq Object ############################

#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts_mod, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts_mod))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

#Convert Model Variable into Factor
RNA.data$Sample_Type_Involved <- as.factor(RNA.data$Sample_Type_Involved)


dds <- DESeqDataSetFromMatrix(countData = d1, colData = RNA.data, design = ~ Sample_Type_Involved)

#Normalization Step 
dds <- estimateSizeFactors(dds)

#Retrive normalized counts matrix 
normalized_counts <- counts(dds, normalized=TRUE)
#save it 
write.table(normalized_counts, file="new_results/normalized_counts_mod.txt", sep="\t", quote=F, col.names=NA)

#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]

#Transform Data
vsd <- varianceStabilizingTransformation(dds)

#Drop Levels
dds$Sample_Type_Involved   <- droplevels(dds$Sample_Type_Involved)
vsd$Sample_Type_Involved   <- droplevels(vsd$Sample_Type_Involved)


###############################################################################
################Running Differential Analysis##################
################Involved Vs Uninvolved (Supp figure 7 ########################

######################Plot PCoA################################

#Create Distance Matrix
vegdist   = vegdist(t(assay(vsd)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(vsd), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Sample_Type_Involved, data = newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type_Involved",suffixes=c("",".centroid"))

pdf("new_figures/NSCLC.Lung.Tissue.Involved.vs.UnIn.Bray.RNA.PCoA_mod.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= Sample_Type_Involved)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  scale_color_manual(values=c("blue", "red")) + 
  #plot ellipse
  #stat_ellipse(type = "t") + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Sample_Type_Involved), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample_Type_Involved))+ 
  #labels centroids 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Lung", "Tumor")), size=10) +
  #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
  scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")+
  ggtitle("PCoA plot,Lung vs. Tumor, p=0.001")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Stats 
adonis2(vegdist ~ vsd$Sample_Type_Involved)
write.table(adonis2(vegdist ~ vsd$Sample_Type_Involved), file="new_results/Bray.Lung.Cohort.Stat.RNA_mod.txt", sep = "\t", row.names = TRUE)

################################################################################
######### PCOA of tumor only: recurrence vs no recurrence (Figure 3A) #####


dds <- DESeqDataSetFromMatrix(countData = d1, colData = RNA.data, design = ~ Progression)
vsd <- varianceStabilizingTransformation(dds)


dds.involved <- dds [, dds$Sample_Type_Involved %in% "Lung.Tissue.In"]
vsd.involved <- vsd[, vsd$Sample_Type_Involved  %in% "Lung.Tissue.In"]

dds.involved$Progression <- factor(dds.involved$Progression)
vsd.involved$Progression <- factor(vsd.involved$Progression)


######################Plot PCoA################################

#Create Distance Matrix
vegdist   = vegdist(t(assay(vsd.involved)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(vsd.involved), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression, data = newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression",suffixes=c("",".centroid"))

pdf("new_figures/Lung.Tissue.Involved.Recu.vs.noRecu.Bray.RNA.PCoA.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= Progression)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  scale_color_manual(values=c("orange", "#E0115F")) + 
  #plot ellipse
  #stat_ellipse(type = "t") + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression))+ 
  #labels centroids 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence")), size=10) +
  #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
  scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")+
  ggtitle("PCoA plot, Tumor: Recurrence vs. No recurrence, p=0.283")+
  theme(plot.title = element_text(hjust = 0.2))
dev.off()

#Stats 
adonis2(vegdist ~ vsd.involved$Progression)
write.table(adonis2(vegdist ~ vsd.involved$Progression), file="new_results/Bray.Lung.Cohort.Inv.Recu.vs.norecur.Stat.RNA.txt", sep = "\t", row.names = TRUE)




########################## recurrent vs non recurrent UnInvolved  PCOA (Figure 3C###############


dds.Uninvolved <- dds [, dds$Sample_Type_Involved %in% "Lung.Tissue.UnIn"]
vsd.Uninvolved <- vsd[, vsd$Sample_Type_Involved  %in% "Lung.Tissue.UnIn"]

dds.Uninvolved$Progression

dds.Uninvolved$Progression <- factor(dds.Uninvolved$Progression)
vsd.Uninvolved$Progression <- factor(vsd.Uninvolved$Progression)

######################Plot PCoA################################

#Create Distance Matrix
vegdist   = vegdist(t(assay(vsd.Uninvolved)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(vsd.Uninvolved), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression, data = newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression",suffixes=c("",".centroid"))

pdf("new_figures/NSCLC.Lung.Tissue.UnInvolved.Recu.vs.noRecu.Bray.RNA.PCoA.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= Progression)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  scale_color_manual(values=c("#92A1CF", "blue")) + 
  #plot ellipse
  #stat_ellipse(type = "t") + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression))+ 
  #labels centroids 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence")), size=10) +
  #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
  scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")+
  ggtitle("PCoA plot,Lung Tissue, Recurrence vs. No Recurrence, p=0.005")+
  theme(plot.title = element_text(hjust = 0.216))
dev.off()

#Stats 
adonis2(vegdist ~ vsd.Uninvolved$Progression)
write.table(adonis2(vegdist ~ vsd.Uninvolved$Progression), file="new_results/Bray.Lung.Cohort.UnInv.Recu.vs.norecur.Stat.RNA.txt", sep = "\t", row.names = TRUE)




################################################################################
########################### edgeR Analysis (Figures 3B, 3D) #####################################
################################################################################

########################### edgeR Tumor vs Lung ################################


#get genes table and add 1
x= as(mycounts_mod, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA.data, "Sample_Type_Involved")
group <- factor(group, levels = c("Lung.Tissue.UnIn", "Lung.Tissue.In"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here#############

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = "new_results/edgeR.results.Tumor.vs.Lung.csv")

#save for GSEA a ranked file

#Histogram for quality assesment 
ggplot(res, aes(x=logFC)) +
  geom_histogram(color="blue", fill="white")

# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file="new_results/Tumor_vs_Lung_for_GSEA.rnk", col_names=FALSE)



###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < 0.2 ] <- "red"
cols[res$logFC < 0 & res$FDR < 0.2 ] <- "blue"


#plot with labels 
pdf(file = "new_figures/edgeR_Volcano_Tumor_vs_Lung_FDR_0.2.pdf", height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>1 & res$FDR < 0.2 , as.character(res$Gene.symbol),
                               ifelse(res$logFC<-1 & res$FDR < 0.2, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




########################### edgeR Tumor vs Recurrence vs No recurrence  #################

#first select tumor samples only from metadata and counts table 
mycounts.Tumor <- mycounts_mod %>% 
  dplyr::select(., starts_with("NYU"))

#get genes table and add 1
x= as(mycounts.Tumor, "matrix")
x= x+1 

#get tumor only table 
RNA.data.Tumor <- RNA.data %>% filter(Sample_Type_Involved =="Lung.Tissue.In")

#overall you have 146 Tumor samples with RNA data 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA.data.Tumor, "Progression_Lab_Inv")

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here#############

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = "new_results/edgeR.results.Tumor_Rec_vs_noRec_.csv")

#save for GSEA a ranked file

#Histogram for quality assesment 
ggplot(res, aes(x=logFC)) +
  geom_histogram(color="blue", fill="white")

# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write.csv(y, file = "new_results/tumor_test.csv", col.names = FALSE)
write_tsv(y, file="new_results/Tumor_Rec_vs_no.Rec_for_GSEA.rnk", col_names=FALSE)



###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "#E0115F"
cols[res$logFC < 0 & res$FDR < alpha ] <- "orange"


#plot with labels 
pdf(file = paste0("new_figures/edgeR_Tumor_Rec_vs_no.Rec_FDR_", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>1 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<-1 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







########################### edgeR Lung vs Recurrence vs No recurrence  #################

#first select Lung samples only from metadata and counts table 
mycounts.Lung <- mycounts_mod %>% 
  dplyr::select(., -starts_with("NYU"))

#get genes table and add 1
x= as(mycounts.Lung, "matrix")
x= x+1 

#get Lung only table 
RNA.data.Lung <- RNA.data %>% filter(Sample_Type_Involved =="Lung.Tissue.UnIn")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA.data.Lung, "Progression_Lab_Inv")

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here#############

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = "new_results/edgeR.results.Lung_Rec_vs_noRec_.csv")

#save for GSEA a ranked file

#Histogram for quality assesment 
ggplot(res, aes(x=logFC)) +
  geom_histogram(color="blue", fill="white")

# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file="new_results/Lung_Rec_vs_no.Rec_for_GSEA.rnk", col_names=FALSE)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "blue"
cols[res$logFC < 0 & res$FDR < alpha ] <- "#92A1CF"


#plot with labels 
pdf(file = paste0("new_figures/edgeR_Lung_Rec_vs_no.Rec_FDR_", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>1 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<-1 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()



###################################################
######### IPA hearmap tumor and lung , rec vs no rec (Figure 3E) ######


IPA_res <- read.csv("new_figures/IPA_res_tumor_lung_rec_vs_norec.csv")

#get max and min values 
max(IPA_res$Tumor, na.rm = TRUE)
min(IPA_res$Tumor, na.rm = TRUE)
max(IPA_res$Lung, na.rm = TRUE)
min(IPA_res$Lung, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-1.213, 0,  7.635), c("#FED976", "white", "#6A51A3"))
col_fun(seq(-1.213, 0,  7.635))

#convert data to matrix 
IPA_res_mat <- IPA_res
IPA_res_mat <-IPA_res_mat %>% arrange(desc(Lung))
#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "new_figures/ipa_for_figure_3E.pdf", height = 30, width = 16)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        column_names_gp = gpar(fontace="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(1.5, "cm"), height = unit(28, "cm"))
dev.off()







########### GSEA results Hallmark heatmap #####



gsea_res <- read.csv("new_results/gsea_hallmark_for_heatmap.csv")

#get max and min values 
max(gsea_res$Tumor, na.rm = TRUE)
min(gsea_res$Tumor, na.rm = TRUE)
max(gsea_res$Lung, na.rm = TRUE)
min(gsea_res$Lung, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-2.237938, 0,  1.765604), c("#FED976", "white", "#6A51A3"))
col_fun(seq(-2.237938, 0,  1.765604))

#convert data to matrix 
gsea_res_mat <- gsea_res
gsea_res_mat <-gsea_res_mat %>% arrange(desc(Lung))
#set pathways as rownames
rownames(gsea_res_mat) <- gsea_res_mat$NAME

#get rid of extra columns
gsea_res_mat <- gsea_res_mat[,-1]

#convert to matrix 
gsea_res_mat <- as.matrix(gsea_res_mat)
library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "new_figures/gsea_ehatmap_hallmark.pdf", height = 16, width = 16)
ComplexHeatmap::Heatmap(gsea_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        column_names_gp = gpar(fontace="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =gsea_res$NAME, 
                        width = unit(1.5, "cm"))
dev.off()

