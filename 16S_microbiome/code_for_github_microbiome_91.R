# 07 APR 2024 ; Fares Darawshy ; Segal Lab 
# Microbial Signatures of Early Lung Cancer 

#Load/Save wd ----- 
setwd(#set your working directory )

#load libraries 
library(qiime2R)
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
library(knitr)
library(decontam)
library(microbiomeMarker)
library(caret)
library(randomForest)
library(matrixStats)
library(PMCMRplus)
library(ggtext)
library(rmarkdown)
library(glue)
library(cowplot)
library(patchwork)
library(table1)

#set seed 
set.seed(1234)


#import data from QIIME2 using q2R

#reading sequence variants 
SVs <- read_qza("table.filtered.qza")

#reading Metadata of 91 patients only 
early_lung_cancer_metadata <- read.csv(file="LNSmetafile_91sub.csv")
rownames(early_lung_cancer_metadata) <- early_lung_cancer_metadata$X
early_lung_cancer_metadata$sample_name <- early_lung_cancer_metadata$X
early_lung_cancer_metadata <- early_lung_cancer_metadata %>% select(-X)


#read metadata of all cohort 
map = "Surgical.Cohort.Map.txt"
mapping.table = sample_data(read.table(map,header = T, sep = "\t", row.names = 1))

#reading Taxonomy 
taxonomy <- read_qza("taxonomy.filtered.qza")
head(taxonomy$data)

#you can see in the output that taxonomy improted is a single string, along with confidence level
#you need to break this string according to level. 
taxonomy <- parse_taxonomy(taxonomy$data)
head(taxonomy)

#creating a phyloseq object from all these files 
physeq <- qza_to_phyloseq(features = "table.filtered.qza",
                          tree = "rooted-tree_quality.filtered.qza",
                          taxonomy =  "taxonomy.filtered.qza",
                          metadata = "Surgical.Cohort.Map.txt")
#inspect 
physeq
otu_table(physeq)
sample_data(physeq)
colnames(tax_table(physeq))

#create sub-tables to use later 
#remove all taxa with 0 abundance.
#none removed in this case 
physeq = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)  

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) x/sum(x)
physeq.rel.table= transform_sample_counts(physeq, normalizeSample)

# Create phyllum and order tables (do it after normalization and out of the relative table)
Phylum.rel.table = tax_glom(physeq.rel.table, taxrank = "Phylum")
Class.rel.table = tax_glom(physeq.rel.table, taxrank = "Class")
Order.rel.table = tax_glom(physeq.rel.table, taxrank = "Order")
Family.rel.table = tax_glom(physeq.rel.table, taxrank = "Family")
Genus.rel.table = tax_glom(physeq.rel.table, taxrank = "Genus")
Species.rel.table = tax_glom(physeq.rel.table, taxrank = "Species")

#table without DFW and mock 
physeq.3 = subset_samples(physeq, Sample_Type_Involved %in% c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"))

physeq.3@sam_data$sample_name <- rownames(sample_data(physeq.3))

#rearrange metadata to include those 91 patients only from metadta 
temp <- data.frame(sample_data(physeq.3))

BKG_data <- temp %>% filter(Sample_Type=="BKG")
lung_tissue_data <- temp %>% filter(Sample_Type=="Lung.Tissue")

temp_2 <- inner_join(lung_tissue_data, early_lung_cancer_metadata, by="sample_name")

temp_2 <- temp_2 %>% select(-ends_with(".y"))

colnames(temp_2) <- gsub(".x", "", colnames(temp_2))

rownames(temp_2) <- temp_2$sample_name

temp_3 <- bind_rows(temp_2, BKG_data)

#create new physeq object 
physeq_mod <- subset_samples(physeq.3, sample_name %in% temp_3$sample_name)

#get relative table 
physeq_mod_rel <- transform_sample_counts(physeq_mod, normalizeSample)

###############################################################################  
###################################Topographical Analysis #####################
###############################################################################  

#########Calculating sequence depth and plot###########
sample_sums(physeq_mod)
sort(sample_sums(physeq_mod))

#you can omit taxa seen in 0 of the samples 
#by comparing number of rows to physeq object, you can see there is no such taxa 
taxa.sums <- taxa_sums(physeq_mod)>0
taxa.sums <- as.data.frame(taxa.sums)
nrow(taxa.sums)

#convert total counts into data frame
library(data.table)
sequence.depth = data.table::data.table(as(sample_data(physeq_mod), "data.frame"), TotalReads=sample_sums(physeq_mod), keep.rownames = TRUE)
sequence.depth

#change 1st column form rn to sample_ID
setnames(sequence.depth, "rn", "Sample_ID")

#convert total counts into data frame
sequence.depth.3 = data.table::data.table (as(sample_data(physeq_mod), "data.frame"), TotalReads=sample_sums(physeq_mod), keep.rownames = TRUE)
sequence.depth.3
#change 1st column form rn to sample_ID
setnames(sequence.depth.3, "rn", "Sample_ID")

#Stats. Pairwise comparison using Wilcoxon 
compare_means(TotalReads ~ Sample_Type_Involved, data = sequence.depth.3)
write.table(compare_means(TotalReads ~ Sample_Type_Involved, data = sequence.depth.3), file ="sequence.depth.stats.new_modified.txt", sep = "\t", row.names = TRUE)
my_comparisons <- list( c("Lung.Tissue.In", "BKG"), c("Lung.Tissue.In", "Lung.Tissue.UnIn"), c("BKG", "Lung.Tissue.UnIn") )

#plot Boxplot ( Figure Supp 2A)
pdf(file = "Surg.Cohort.Sequence.Depth.BKG.Tumor.Lung.boxplot_modified.pdf", height = 7, width = 6)
ggplot(data=sequence.depth.3, mapping = aes(x=Sample_Type_Involved, y=TotalReads))+
  geom_boxplot(fill=c("grey", "red", "blue"))+
  geom_jitter()+
  labs(x= '', 
       y='Reads Count')+
  scale_x_discrete(labels = c("BKG", "Tumor", "Lung"))+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(mapping=aes(label=..p.value..), comparisons=my_comparisons, method = "wilcox.test")+
  theme_classic()+
  theme(axis.text.x = element_text(size=20, face = "bold", colour = "black"), 
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 20, face = "bold", colour = "black"), 
        legend.position = "none")
dev.off()

#median of all samples sequence depth 
sequence.depth.3 %>% 
  summarise(mean=mean(TotalReads), 
            median=median(TotalReads), 
            q1 = quantile(TotalReads, 0.25), 
            q3 = quantile(TotalReads, 0.75), 
            iqr = IQR(TotalReads))

sequence.depth.3.measures <- sequence.depth.3 %>% 
  group_by(Sample_Type_Involved) %>% 
  summarise(mean=mean(TotalReads), 
            median=median(TotalReads), 
            q1 = quantile(TotalReads, 0.25), 
            q3 = quantile(TotalReads, 0.75), 
            iqr = IQR(TotalReads))

write.csv(sequence.depth.3.measures, file = "sequence.depth.statistical.measures_modified.csv")

################################################################################
################################Rarefaction (Suppfigure 2D)################################### 
#calculation and plot 
set.seed(1234)

psdata.three.groups <- physeq_mod
sample_sums(psdata.three.groups)

calculate_rarefaction_curves <- function(psdata.three.groups, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata.three.groups, measures, depth) {
    if(max(sample_sums(psdata.three.groups)) < depth) return()
    psdata.three.groups <- prune_samples(sample_sums(psdata.three.groups) >= depth, psdata.three.groups)
    
    rarified_psdata <- rarefy_even_depth(psdata.three.groups, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata.three.groups, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata.three.groups, c('Shannon'), rep(c(1, 10, 100, 200, 400, 600, 800, 1000, 1:100 * 10000), each = 100))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata.three.groups)), by.x = 'Sample', by.y = 'row.names')

library(scales)
library(ggforce)
pdf (file = "rarefaction.curve.by.sample.group.three.groups.Shannon.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line() + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange() + 
  ylab("Shannon Index")+
  xlab("Sequence Depth")+
  facet_wrap(facets = ~ Measure, scales = 'free_y')+ 
  #facet_zoom(xlim=c(0,10000))+
  #scale_x_continuous(limits = c(0,10000))+
  theme_bw()+
  theme(legend.title = element_blank())
dev.off()

#zoom on low depth: 1000
pdf (file = "rarefaction.curve.by.sample.group.Zoomed.1000.Shannon.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line(alpha=0.7) + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange(shape=1, alpha=0.5) + 
  ylab("Shannon Index")+
  xlab("Sequence Depth")+
  facet_wrap(facets = ~ Measure, scales = 'free_y')+ 
  facet_zoom(xlim=c(0,1000))+
  scale_x_continuous(limits = c(0,50000))+
  theme_bw()+
  theme(legend.title = element_blank())
dev.off()

#zoom on low depth: 30K
pdf (file = "rarefaction.curve.by.sample.group.Zoomed.30K.Shannon.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line() + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange() + 
  ylab("Shannon Index")+
  xlab("Sequence Depth")+
  facet_wrap(facets = ~ Measure, scales = 'free_y')+ 
  facet_zoom(xlim=c(0,30000))+
  #scale_x_continuous(limits = c(0,2000))+
  theme_bw()+
  theme(legend.title = element_blank())
dev.off()

#zoom on low depth: 50K
pdf (file = "rarefaction.curve.by.sample.group.Zoomed.50K.Shannon.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line() + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange() + 
  ylab("Shannon Index")+
  xlab("Sequence Depth")+
  facet_wrap(facets = ~ Measure, scales = 'free_y')+ 
  facet_zoom(xlim=c(0,50000))+
  #scale_x_continuous(limits = c(0,2000))+
  theme_bw()+
  theme(legend.title = element_blank())
dev.off()


#repeat for Observed 

#calculation and plot 
set.seed(1234)

psdata.three.groups <- physeq_mod
sample_sums(psdata.three.groups)

calculate_rarefaction_curves <- function(psdata.three.groups, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata.three.groups, measures, depth) {
    if(max(sample_sums(psdata.three.groups)) < depth) return()
    psdata.three.groups <- prune_samples(sample_sums(psdata.three.groups) >= depth, psdata.three.groups)
    
    rarified_psdata <- rarefy_even_depth(psdata.three.groups, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata.three.groups, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata.three.groups, c('Observed'), rep(c(1, 10, 100, 200, 400, 600, 800, 1000, 1:100 * 10000), each = 100))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary$Sample <- gsub("X", "",rarefaction_curve_data_summary$Sample)

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata.three.groups)), by.x = 'Sample', by.y = 'row.names')

library(scales)
library(ggforce)
pdf (file = "rarefaction.curve.by.sample.group.three.groups.Observed.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line() + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange() + 
  ylab("ASVs")+
  xlab("Reads Count")+
  #facet_zoom(xlim=c(0,50000))+
  #scale_x_continuous(limits = c(0,2000))+
  scale_y_continuous(breaks = seq(0,350,25))+
  theme_classic()+
  theme(legend.title = element_blank(), 
        legend.position = "top",
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
dev.off()

pdf (file = "rarefaction.curve.by.sample.group.three.groups.Zoomed.2000.Observed.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line(alpha =0.5) + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange(alpha =0.5) + 
  ylab("ASVs")+
  xlab("Reads Count")+
  facet_zoom(xlim=c(0,2000))+
  scale_x_continuous(limits = c(0,50000))+
  scale_y_continuous(breaks = seq(0,350,25))+
  theme_bw()+
  theme(legend.title = element_blank(), 
        legend.position = "top",
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"), 
        strip.background = element_rect(color = "lightgrey", fill="white"), 
        panel.background = element_rect(fill = "white"))
dev.off()

pdf (file = "rarefaction.curve.by.sample.group.three.groups.Zoomed.10000.Observed.new_mod.pdf", width = 8, height = 8)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    color=Sample_Type_Involved,
    group = Sample)) + 
  geom_line(alpha =0.5) + 
  scale_color_manual(labels = c("BKG","Tumor", "Lung")  ,values = c("grey", "red", "blue"))+
  geom_pointrange(alpha =0.5) + 
  ylab("ASVs")+
  xlab("Sequence Depth")+
  facet_zoom(xlim=c(0,10000))+
  scale_x_continuous(limits = c(0,50000))+
  scale_y_continuous(breaks = seq(0,350,25))+
  theme_classic()+
  theme(legend.title = element_blank(), 
        legend.position = "top",
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
dev.off()


###############################################################################
###############Alpha Diversity of all cohort including BKG (Supp Figure 2B 2C) data used to plot supp figure 6 in Prism ###################


#calculate alpha diversity measures 
alpha.measures <- data.frame("Shannon" = estimate_richness(physeq_mod, measures = "Shannon"), 
                             "Observed"= estimate_richness(physeq_mod, measures = "Observed"), 
                             "Chao1"= estimate_richness(physeq_mod, measures = "Chao1"),
                             "Sample_Type" = sample_data(physeq_mod)$Sample_Type_Involved, 
                             "Progression_Lab_Inv" = sample_data(physeq_mod)$Progression_Lab_Inv)
#explore your alpha.measures
alpha.measures
alpha.measures$Sample_Type <- as.factor(alpha.measures$Sample_Type)

###import alpha measures to plot in prism 
write.csv(alpha.measures, file = "new_results/alpha_measures_for_prism.csv")

#stats and plot for Shannon
compare_means(Shannon ~ Sample_Type, data = alpha.measures)
write.table(compare_means(Shannon ~ Sample_Type, data = alpha.measures), file ="Shannon.stats.new_mod.txt", sep = "\t", row.names = TRUE)
my_comparisons <- list( c("Lung.Tissue.In", "BKG"), c("Lung.Tissue.In", "Lung.Tissue.UnIn"), c("BKG", "Lung.Tissue.UnIn") )

#Plot Shannon
pdf(file = "Shannon.Index.By.Sample.Type.new_mod.pdf", height = 7, width = 6)
ggplot(alpha.measures, aes(x=Sample_Type, y=Shannon))+
  geom_boxplot(fill=c("grey", "red", "blue"), outlier.colour = NA)+
  geom_jitter()+
  xlab("")+ylab("Shannon Index")+
  scale_x_discrete(labels = c("BKG", "Tumor", "Lung"))+
  stat_compare_means(mapping=aes(label=p.value), comparisons = my_comparisons, method = "wilcox.test")+
  theme_classic()+
  theme(axis.text.x = element_text(size=20, face = "bold", colour = "black"), 
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 20, face = "bold", colour = "black"), 
        legend.position = "none")
dev.off()

#stats and plot for observed 
compare_means(Observed ~ Sample_Type, data = alpha.measures)
write.table(compare_means(Observed ~ Sample_Type, data = alpha.measures), file ="Observed.stats.new_mod.txt", sep = "\t", row.names = TRUE)
my_comparisons <- list( c("Lung.Tissue.In", "BKG"), c("Lung.Tissue.In", "Lung.Tissue.UnIn"), c("BKG", "Lung.Tissue.UnIn") )

#Plot Observed
pdf(file = "Alpha.Observed.By.Sample.Type.new_mod.pdf", height = 7, width = 6)
ggplot(alpha.measures, aes(x=Sample_Type, y=Observed))+
  geom_boxplot(fill=c("grey", "red", "blue"), outlier.colour = NA)+
  geom_jitter(width = 0.2)+
  xlab("")+ylab("Richness")+
  scale_x_discrete(labels = c("BKG", "Tumor", "Lung"))+
  stat_compare_means(mapping=aes(label=..p.adj..), comparisons = my_comparisons, method = "wilcox.test")+
  theme_classic()+
  theme(axis.text.x = element_text(size=20, face = "bold", colour = "black"), 
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 20, face = "bold", colour = "black"), 
        legend.position = "none")
dev.off()

##################################Beta Diversity PCOA (Supp Figure 3) ################################


library(vegan)
####Creating PCOA plot for all samples#######
#Create Distance Matrix
vegdist   = vegdist(t(otu_table(physeq_mod_rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(physeq_mod_rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Sample_Type_Involved,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type_Involved",suffixes=c("",".centroid"))

pdf("Lung.Cohort.Site.three.groups.Bray.PCoA.new_mod.pdf", height = 10, width = 10)
ggplot(newResults, aes(PC1, PC2, color= Sample_Type_Involved)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("grey", "red", "blue")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Sample_Type_Involved), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample_Type_Involved))+ 
  #If you want to identify specific samples use the code bellow
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BKG", "Tumor", "Lung")), size=12) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
  #Use the code bellow if you want to switch the X axis around
  #scale_x_reverse() +
  #ggtitle("PCoA plot of Bray Distance. p<0.001")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#stats of bray 
#overall p value
adonis2(vegdist~ sample_data(physeq_mod_rel)$Sample_Type_Involved)

#BKG vs lung 
ps.re.1 <- subset_samples(physeq_mod_rel, Sample_Type_Involved %in% c("BKG", "Lung.Tissue.UnIn"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$Sample_Type_Involved), file = "Bray.Stats.all.cohort.BKG.vs.Lung_mod.csv")

#BKG vs Tumor 
ps.re.1 <- subset_samples(physeq_mod_rel, Sample_Type_Involved %in% c("BKG", "Lung.Tissue.In"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$Sample_Type_Involved), file = "Bray.Stats.all.cohort.BKG.vs.Tumor_mod.csv")

#Tumor vs Lung
ps.re.1 <- subset_samples(physeq_mod_rel, Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$Sample_Type_Involved), file = "Bray.Stats.all.cohort.Tumor.vs.Lung_mod.csv")



##### write 16s  plot in prism ### 
data_16s <- sample_data(physeq_mod) %>% data.frame() %>% select(c(sample_name, Sample_Type_Involved, Progression_Lab_Inv, Final.Conce_per_uL))
write.csv(data_16s, file = "new_results/16S_data_for_prism.csv")

################################################################################
#############################Decontam using KW function (Supp Figure 4) #########################
#################################################################################


decontaminant_KW(input_phyloseq = physeq_mod, 
                 SampleID.unique ="sample_name",
                 sample_type_var_name = "Sample_Type_Involved",
                 sample_types=c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"), 
                 negative_sample_type=c("BKG"), 
                 sample_type_color=c("grey" ,"red", "blue"),
                 stat_option="mean", 
                 compare_type = c("Lung.Tissue.In_and_Lung.Tissue.UnIn", "Lung.Tissue.In_or_Lung.Tissue.UnIn"),
                 display_contam_method = "none",
                 graph_option ="boxplot",
                 test_threshold=0.1, 
                 log_scale = "yes", 
                 taxa_genus_output = "no") 



decontaminant_KW(input_phyloseq = physeq_mod, 
                 SampleID.unique ="sample_name",
                 sample_type_var_name = "Sample_Type_Involved",
                 sample_types=c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"), 
                 negative_sample_type=c("BKG"), 
                 sample_type_color=c("grey" ,"red", "blue"),
                 stat_option="mean", 
                 compare_type = c("Lung.Tissue.In_and_Lung.Tissue.UnIn", "Lung.Tissue.In_or_Lung.Tissue.UnIn"),
                 display_contam_method = "none",
                 graph_option ="boxplot",
                 test_threshold=0.25, 
                 log_scale = "yes", 
                 taxa_genus_output = "no") 

decontaminant_KW(input_phyloseq = physeq_mod, 
                 SampleID.unique ="sample_name",
                 sample_type_var_name = "Sample_Type_Involved",
                 sample_types=c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"), 
                 negative_sample_type=c("BKG"), 
                 sample_type_color=c("grey" ,"red", "blue"),
                 stat_option="mean", 
                 compare_type = c("Lung.Tissue.In_and_Lung.Tissue.UnIn", "Lung.Tissue.In_or_Lung.Tissue.UnIn"),
                 display_contam_method = "none",
                 graph_option ="boxplot",
                 test_threshold=0.5, 
                 log_scale = "yes", 
                 taxa_genus_output = "no") 


decontaminant_KW(input_phyloseq = physeq_mod, 
                 SampleID.unique ="sample_name",
                 sample_type_var_name = "Sample_Type_Involved",
                 sample_types=c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"), 
                 negative_sample_type=c("BKG"), 
                 sample_type_color=c("grey" ,"red", "blue"),
                 stat_option="mean", 
                 compare_type = c("Lung.Tissue.In_and_Lung.Tissue.UnIn", "Lung.Tissue.In_or_Lung.Tissue.UnIn"),
                 display_contam_method = "none",
                 graph_option ="boxplot",
                 test_threshold=0.75, 
                 log_scale = "yes", 
                 taxa_genus_output = "no") 

decontaminant_KW(input_phyloseq = physeq_mod, 
                 SampleID.unique ="sample_name",
                 sample_type_var_name = "Sample_Type_Involved",
                 sample_types=c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"), 
                 negative_sample_type=c("BKG"), 
                 sample_type_color=c("grey" ,"red", "blue"),
                 stat_option="mean", 
                 compare_type = c("Lung.Tissue.In_and_Lung.Tissue.UnIn", "Lung.Tissue.In_or_Lung.Tissue.UnIn"),
                 display_contam_method = "none",
                 graph_option ="boxplot",
                 test_threshold=0.9, 
                 log_scale = "yes", 
                 taxa_genus_output = "no") 



####### final figures of decontam and contamlist #####


#use subplot function 
decontaminant_subplot_KW(input_phyloseq = physeq_mod, 
                         SampleID.unique ="sample_name",
                         sample_types= c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"),
                         sample_type_var_name="Sample_Type_Involved",                                       
                         sample_type_color=c("grey", "red", "blue"), 
                         negative_sample_type="BKG", ###the sample type that you want to be as negative control
                         compare_type = "Lung.Tissue.In_and_Lung.Tissue.UnIn",
                         stat_option = "mean", 
                         method_type = "preval",
                         test_threshold = 0.5,
                         graph_option = "boxplot", 
                         taxa_genus_output="no", 
                         log_scale = "yes")

phyloseq::tax_table(physeq_mod)


#get final contamlist (as Supp Table 1)
contamlist <- read.csv(file = "contamlist_new.csv")
contamlist <- data.frame(contamlist)
contamlist$asv <- make.unique(contamlist$asv)
rownames(contamlist) <- contamlist$asv
contamlist <- contamlist %>% select(asv, contaminant)

##### differntial edgeR Tumor vs lung analysis (remove BKG) Supp Figure 5################

ps <- subset_samples(physeq_mod, Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"))
ps.rel <- transformSampleCounts(ps, normalizeSample)

ps@sam_data$Sample_Type_Involved <- factor(ps@sam_data$Sample_Type_Involved, levels = c("Lung.Tissue.UnIn", "Lung.Tissue.In"))
ps.rel@sam_data$Sample_Type_Involved <- factor(ps.rel@sam_data$Sample_Type_Involved, levels = c("Lung.Tissue.UnIn", "Lung.Tissue.In"))


#####################edgeR ##########################################

#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "Sample_Type_Involved", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)


# make taxa name column and name asv column 

res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, Sample_Type_Involved=="Lung.Tissue.In")
ps.rel.2 <- subset_samples(ps.rel, Sample_Type_Involved=="Lung.Tissue.UnIn")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.2 <- meanRA.save

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

#get abundance of each category alone to save it in the table 
contamlist
#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "tumor", "lung"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.1, abundance.2,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "new_results/edgeR.results_Sample_Type_Involved_tumor_vs_lung.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="tumor", res$abundance.1, 
                        res$abundance.2)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 1 & res$FDR < 0.2 ] <- "red"
cols[res$logFC < -1 & res$FDR < 0.2 ] <- "blue"

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(category=="tumor") %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(category=="lung") %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "blue","darkgrey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "blue", "red"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="new_figures/edgeR_Bubble_Sample_Type_Involved_tumor_vs_lung.pdf", height = 11, width = 12)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="new_figures/edgeR_Bubble_Sample_Type_Involved_tumor_vs_lung_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 



####### Table 1 ########

Table.1.data <- sample_data(physeq_mod) %>% data.frame() %>% filter(Sample_Type_Involved =="Lung.Tissue.In")

#add statistical analysis column 
#define p value function 

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

#create table 1 using table 1 package 
table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         factor(Male_1) + factor(RACE) + as.numeric(Age) + factor(Smoking_status) + as.numeric(Pack_Years) 
       + factor(Stage) +factor(Outcome) + factor(Progression_Lab) + factor(ProgType_Lab)+ as.numeric(Size_cm)
       | Progression, 
       #set data 
       data=Table.1.data, 
       #set stat display options for continous variables (options from stat.default)
       render.continuous = c((.="Median [Q1, Q3]")))
#now add stats (same as above except where notes added ) for these two groups only 
table1(~ 
         factor(Male_1) + factor(RACE) + as.numeric(Age) + factor(Smoking_status) + as.numeric(Pack_Years) 
       + factor(Stage) +factor(Outcome) + factor(Progression_Lab) + factor(ProgType_Lab)+ as.numeric(Size_cm)
       | Progression, data=Table.1.data, 
       # don't display overall column 
       overall = F,
       render.continuous = c((.="Median [Q1, Q3]")),
       #add p value as extra column (defined in function above)
       extra.col = list("P-value" = pvalue))

#get p value of numerical variables 
t.test(as.numeric(Table.1.data$Age) ~ Table.1.data$Progression_Lab) # p = 0.09795
t.test(as.numeric(Table.1.data$Pack_Years) ~ Table.1.data$Progression_Lab) #p-value = 0.4096
t.test(as.numeric(Table.1.data$Size_cm) ~ Table.1.data$Progression_Lab) #p-value = 0.4096

#copy and paste the above results and edit in excel / word 
#5 year mortality data 
mort_dat <- read_xlsx("5_year_survival_updated.xlsx")
early_lung_cancer_metadata$Sample_ID<- rownames(early_lung_cancer_metadata)

mort_dat <- inner_join(mort_dat, early_lung_cancer_metadata, by="Sample_ID")
mort_dat$outcome_in_5_years
mort_dat %>% 
  group_by(group, outcome_in_5_years) %>% 
  summarise(count=n())

contingency_table <- table(mort_dat$group, mort_dat$outcome_in_5_years)

# Perform chi-square test
chi_square_test <- chisq.test(contingency_table)

# Print the results
print(chi_square_test)

########## KM recurrence curve (Figure 1A) #######
########################Kaplan Meir Curve for Recurrence########################

#Load relevant packages 
library(dplyr)
library(survival)
library(survminer)

#use Genus level table
#Use subset of samples with tumor only 
survival_data <- sample_data(physeq_mod) %>% data.frame() %>% filter(Sample_Type_Involved =="Lung.Tissue.In")

#search for columns to keep 
survival_data$New_TTP
survival_data$Progression
survival_data$New_time_followup.or.death

#Column to be used as time (in days) - New_TTP 
# column to be used as status is Progression Lab
survival.data <- survival_data[, c("New_TTP","Progression", "Progression_Lab_Inv" )]

#rename the new df 
colnames(survival.data) <- c("time", "status", "Recurrence")

# turn the data into data frame so it can be used in survival function. 
#turn the class of columns so its either factor or numeric to be used in survival functions 
survival.data<- data.frame(survival.data)
survival.data$time <- as.numeric(survival.data$time)
survival.data$status <- as.numeric(survival.data$status)
survival.data$Recurrence <- as.factor(survival.data$Recurrence)
#get median survival in each group 
survival.data %>% 
  group_by(Recurrence) %>% 
  summarise(median=median(time), 
            Q1=quantile(time, 0.25), 
            Q3=quantile(time, 0.75))
#create survival object

surv <- Surv(survival.data$time, survival.data$status)

#create survival df
sfit <- survfit(Surv(time, status)~Recurrence, data=survival.data)
sfit

survival.plot <- ggsurvplot(sfit, data = survival.data,
                            fun = "pct", ggtheme = theme_pubr(), 
                            conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
                            legend.labs=c("No", "Yes"), legend.title="Recurrence",  
                            palette=c("dodgerblue2", "red"), size=1,
                            title="Kaplan-Meier Curve for Recurrence Following Surgery",
                            ylab="Percent Without Recurrence",
                            xlab="Time(Days)",
                            # tables.col = "strata", #color of numbers inside risk table
                            tables.height = 0.15,
                            fontsize = 4,
                            risk.table.y.text.col = T,
                            cumcensor = TRUE, 
                            tables.theme = theme_cleantable(),
                            # risk.table.y.text = FALSE,
                            risk.table.title = "Number at risk (cumulative number of recurrence)",
                            cumcensor.title = "Cumulative number of censored subjects"
)


#save the plot 
pdf(file = "Early.Lung.Cancer.Kaplan.Meier.Plot.pdf", width = 8, height = 8)
survival.plot
dev.off()








######## Prune the data to continue analysis ################
physeq_mod.noBKG <- subset_samples(physeq_mod, Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"))
physeq_mod.noBKG.rel <- transform_sample_counts(physeq_mod.noBKG, normalizeSample)
#Pruning taxa strategy 
pruned.taxa <- genefilter_sample(physeq_mod.noBKG.rel, filterfun_sample(function(x) x> 0.001), A = 0.01 * nsamples(physeq_mod.noBKG.rel))

#prune relative table 
physeq_mod.noBKG.rel_pruned <- prune_taxa(pruned.taxa, physeq_mod.noBKG.rel)

physeq_mod.noBKG.rel_pruned

phyloseq::plot_bar(physeq_mod.noBKG.rel_pruned)+geom_hline(yintercept = 0.7)


#prune count table 
physeq_mod.noBKG_pruned <- prune_taxa(pruned.taxa, physeq_mod.noBKG)


#prepare Tumor only tables 
In.Lung.Tissue.pruned <- subset_samples(physeq_mod.noBKG_pruned, sample_data(physeq_mod.noBKG_pruned)$Sample_Type_Involved == "Lung.Tissue.In")
In.Lung.Tissue.pruned.rel <- transform_sample_counts(In.Lung.Tissue.pruned, normalizeSample)

#prepare Lung only tables 
UnIn.Lung.Tissue.pruned <- subset_samples(physeq_mod.noBKG_pruned, sample_data(physeq_mod.noBKG_pruned)$Sample_Type_Involved == "Lung.Tissue.UnIn")
UnIn.Lung.Tissue.pruned.rel <- transform_sample_counts(UnIn.Lung.Tissue.pruned, normalizeSample)





############ Involved subset analysis by recurrence (Figure 1) ##################

ps <- In.Lung.Tissue.pruned
ps.rel <- In.Lung.Tissue.pruned.rel

#ddPCR 
ddpcr_dat <- data.frame(ps@sam_data)

ddpcr_dat <- ddpcr_dat %>% dplyr::mutate(Final.Conce_per_uL=as.numeric(Final.Conce_per_uL))

#calculate stats for ddpcr
my.comparisons <- compare_means(Final.Conce_per_uL~Progression_Lab_Inv, data = ddpcr_dat)
#save it
write.csv(my.comparisons, file = "new_results/ddPCR.stats.tumor_recurrence.csv")


# plot it
pdf(file = "new_figures/ddPCR_tumor_recurrence.pdf", width=6, height=9)
ggplot(ddpcr_dat, aes(x=Progression_Lab_Inv, y=Final.Conce_per_uL, fill=Progression_Lab_Inv))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange", "#E0115F"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("16S rRNA copies/uL")+
  scale_x_discrete(labels = c("No Recurrence", "Recurrence"))+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(comparisons = list(c("In.Recurrence", "In.No.Recurrence")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



#calculate alpha diversity measures 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed"= estimate_richness(ps, measures = "Observed"), 
                             "Chao1"= estimate_richness(ps, measures = "Chao1"),
                             "Progression_Lab_Inv" = sample_data(ps)$Progression_Lab_Inv)
#explore your alpha.measures
alpha.measures
alpha.measures$Progression_Lab_Inv <- as.factor(alpha.measures$Progression_Lab_Inv)

#stats and plot for Shannon
compare_means(Shannon ~ Progression_Lab_Inv, data = alpha.measures)
write.table(compare_means(Shannon ~ Progression_Lab_Inv, data = alpha.measures), file ="new_results/Shannon.stats.recurrence.tumor.new_mod.txt", sep = "\t", row.names = TRUE)

#Plot Shannon
pdf(file = "new_figures/Shannon.Index.By.recurrence.tumor.new_mod.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=Progression_Lab_Inv, y=Shannon, fill=Progression_Lab_Inv))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange", "#E0115F"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon Index")+
  scale_x_discrete(labels = c("No Recurrence", "Recurrence"))+
  stat_compare_means(comparisons = list(c("In.Recurrence", "In.No.Recurrence")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#stats and plot for observed 
compare_means(Observed ~ Progression_Lab_Inv, data = alpha.measures)
write.table(compare_means(Observed ~ Progression_Lab_Inv, data = alpha.measures), file ="new_results/Observed.stats.recurrence.tumor.new_mod.txt", sep = "\t", row.names = TRUE)

#Plot Observed
pdf(file = "new_figures/Observed.Index.By.recurrence.tumor.new_mod.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=Progression_Lab_Inv, y=Observed, fill=Progression_Lab_Inv))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange", "#E0115F"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("No Recurrence", "Recurrence"))+
  stat_compare_means(comparisons = list(c("In.Recurrence", "In.No.Recurrence")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

##################################Beta Diversity PCOA################################


library(vegan)
####Creating PCOA plot for all samples#######
#Create Distance Matrix
vegdist   = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression_Lab_Inv,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression_Lab_Inv",suffixes=c("",".centroid"))

pdf("new_figures/beta_diversity_PCOA_tumor_recurrence.pdf", height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Progression_Lab_Inv)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange", "#E0115F")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression_Lab_Inv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression_Lab_Inv)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  ggtitle("Beta Diversity, Bray, Progression_Lab_Inv, p=0.12")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
write.table(adonis2(vegdist ~ newResults$Progression_Lab_Inv), 
            file = "new_results/Beta.Diversity.Bray.tumor_Progression_Lab_Inv.txt", sep = "\t", row.names = T)




###edgeR involved subset ###
dge <- phyloseq_to_edgeR(physeq = ps, group = "Progression_Lab_Inv", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)


# make taxa name column and name asv column 

res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, Progression_Lab_Inv=="In.Recurrence")
ps.rel.2 <- subset_samples(ps.rel, Progression_Lab_Inv=="In.No.Recurrence")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.2 <- meanRA.save

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

#get abundance of each category alone to save it in the table 
contamlist
#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "2"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.1, abundance.2,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "new_results/edgeR.results_tumor_Progression_Lab_In.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 1 & res$FDR < 0.2 ] <- "#E0115F"
cols[res$logFC < -1 & res$FDR < 0.2 ] <- "orange"

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(category=="1") %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(category=="2") %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"#E0115F", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange","darkgrey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "#E0115F", "orange"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="new_figures/edgeR_Bubble_tumor_Progression_Lab_In.pdf", height = 11, width = 12)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="new_figures/edgeR_Bubble_tumor_Progression_Lab_In_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 



















############ UnInvolved subset analysis by recurrence Figure 1 ##################

ps <- UnIn.Lung.Tissue.pruned
ps.rel <- UnIn.Lung.Tissue.pruned.rel

#ddPCR 
ddpcr_dat <- data.frame(ps@sam_data)

ddpcr_dat <- ddpcr_dat %>% dplyr::mutate(Final.Conce_per_uL=as.numeric(Final.Conce_per_uL))

#calculate stats for ddpcr
my.comparisons <- compare_means(Final.Conce_per_uL~Progression_Lab_Inv, data = ddpcr_dat)
#save it
write.csv(my.comparisons, file = "new_results/ddPCR.stats.lung_recurrence.csv")


# plot it
pdf(file = "new_figures/ddPCR_lung_recurrence.pdf", width=6, height=9)
ggplot(ddpcr_dat, aes(x=Progression_Lab_Inv, y=Final.Conce_per_uL, fill=Progression_Lab_Inv))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("#92A1CF", "blue"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("16S rRNA copies/uL")+
  scale_x_discrete(labels = c("No Recurrence", "Recurrence"))+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(comparisons = list(c("UnIn.Recurrence", "UnIn.No.Recurrence")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



#calculate alpha diversity measures 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed"= estimate_richness(ps, measures = "Observed"), 
                             "Chao1"= estimate_richness(ps, measures = "Chao1"),
                             "Progression_Lab_Inv" = sample_data(ps)$Progression_Lab_Inv)
#explore your alpha.measures
alpha.measures
alpha.measures$Progression_Lab_Inv <- as.factor(alpha.measures$Progression_Lab_Inv)

#stats and plot for Shannon
compare_means(Shannon ~ Progression_Lab_Inv, data = alpha.measures)
write.table(compare_means(Shannon ~ Progression_Lab_Inv, data = alpha.measures), file ="new_results/Shannon.stats.recurrence.lung.new_mod.txt", sep = "\t", row.names = TRUE)

#Plot Shannon
pdf(file = "new_figures/Shannon.Index.By.recurrence.lung.new_mod.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=Progression_Lab_Inv, y=Shannon, fill=Progression_Lab_Inv))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("#92A1CF", "blue"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon Index")+
  scale_x_discrete(labels = c("No Recurrence", "Recurrence"))+
  stat_compare_means(comparisons = list(c("UnIn.Recurrence", "UnIn.No.Recurrence")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#stats and plot for observed 
compare_means(Observed ~ Progression_Lab_Inv, data = alpha.measures)
write.table(compare_means(Observed ~ Progression_Lab_Inv, data = alpha.measures), file ="new_results/Observed.stats.recurrence.lung.new_mod.txt", sep = "\t", row.names = TRUE)

#Plot Observed
pdf(file = "new_figures/Observed.Index.By.recurrence.lung.new_mod.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=Progression_Lab_Inv, y=Observed, fill=Progression_Lab_Inv))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("#92A1CF", "blue"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("No Recurrence", "Recurrence"))+
  stat_compare_means(comparisons = list(c("UnIn.Recurrence", "UnIn.No.Recurrence")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

##################################Beta Diversity PCOA################################


library(vegan)
####Creating PCOA plot for all samples#######
#Create Distance Matrix
vegdist   = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression_Lab_Inv,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression_Lab_Inv",suffixes=c("",".centroid"))

pdf("new_figures/beta_diversity_PCOA_lung_recurrence.pdf", height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Progression_Lab_Inv)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("#92A1CF", "blue")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression_Lab_Inv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression_Lab_Inv)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  ggtitle("Beta Diversity, Bray, Progression_Lab_Inv, p=0.001")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
write.table(adonis2(vegdist ~ newResults$Progression_Lab_Inv), 
            file = "new_results/Beta.Diversity.Bray.lung_Progression_Lab_Inv.txt", sep = "\t", row.names = T)




###edgeR involved subset ###
dge <- phyloseq_to_edgeR(physeq = ps, group = "Progression_Lab_Inv", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)


# make taxa name column and name asv column 

res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, Progression_Lab_Inv=="UnIn.Recurrence")
ps.rel.2 <- subset_samples(ps.rel, Progression_Lab_Inv=="UnIn.No.Recurrence")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.2 <- meanRA.save

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

#get abundance of each category alone to save it in the table 
contamlist
#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "2"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.1, abundance.2,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "new_results/edgeR.results_lung_Progression_Lab_In.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 1 & res$FDR < 0.2 ] <- "blue"
cols[res$logFC < -1 & res$FDR < 0.2 ] <- "#92A1CF"

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(category=="1") %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(category=="2") %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"blue", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "#92A1CF","darkgrey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "#92A1CF", "blue"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="new_figures/edgeR_Bubble_lung_Progression_Lab_In.pdf", height = 11, width = 12)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="new_figures/edgeR_Bubble_lung_Progression_Lab_In_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 













######## Forest Plot of Cox Model Results (Figure 2C, E) #############

library(readxl)
library(ggtext)


####Tumor results##### 

#only significant results with padj < 0.2 were choosen 

res <- readxl::read_excel("COX_PH_tumor.xlsx")

#match asv column from taxatable 
taxtable <- data.frame(tax_table(physeq))
#now create a column called names and paste name at Genus or species level 
taxtable$names <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",taxtable$Species,sep=""),
                               ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,sep=""),
                                      ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",sep=""),
                                             ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",sep=""),
                                                    ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",sep=""),
                                                           ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",sep=""),
                                                                  ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",sep=""), paste(rownames(taxtable))))))))))

#create a column "name ASV" that include the ASV (or OTU)
taxtable$nameASV <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",rownames(taxtable),sep=""),
                                 ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,"_",rownames(taxtable),sep=""),
                                        ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",rownames(taxtable),sep=""),
                                               ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",rownames(taxtable),sep=""),
                                                      ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",rownames(taxtable),sep=""),
                                                             ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",rownames(taxtable),sep=""),
                                                                    ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",rownames(taxtable),sep=""), paste(rownames(taxtable))))))))))


#create a column of asv which is the row name 
taxtable$asv <- rownames(taxtable)

#now merge the table with your results according to asv column 
Tumor.forest.df <- inner_join(res, taxtable, by="asv")
Tumor.forest.df$nameASV.trim <- paste0(Tumor.forest.df$names,paste0(".."), paste0(str_sub(Tumor.forest.df$asv, - 3, - 1)))
Tumor.forest.df$FDR <- as.numeric(Tumor.forest.df$FDR)

#set color vector accroding to contam
Tumor.forest.df <- Tumor.forest.df %>% 
  mutate(contaminant= factor(contaminant)) %>% 
  mutate(color=ifelse(Tumor.forest.df$contaminant =="TRUE", "red", "black")) %>% 
  mutate(upper_CI= as.numeric(upper_CI)) %>% 
  mutate(lower_CI= as.numeric(lower_CI)) %>% 
  mutate(HR= as.numeric(HR)) %>% 
  mutate(nameASV.trim.colored = Tumor.forest.df$nameASV.trim) %>% 
  filter(FDR<=0.2) %>% 
  arrange(., desc(HR)) %>% 
  dplyr::mutate(ord = order(as.numeric(rownames(.))))

#add color of potential contaminants 
Tumor.forest.df$nameASV.trim.colored <- paste0("<span style=\"color: ", Tumor.forest.df$color, "\">", Tumor.forest.df$nameASV.trim.colored, "</span>")

#plot 
Tumor.forest.plot <- ggplot(data=Tumor.forest.df, aes(x=HR, y=fct_reorder(nameASV.trim.colored, -ord), xmin=lower_CI, xmax=upper_CI))+
  geom_point(size=5)+
  geom_errorbarh(height=0.2)+
  geom_vline(xintercept = 1, color="black", linetype="dashed", alpha=0.5)+
  ylab("")+xlab("Hazard Ratio")+
  xlim(0,3.5)+
  theme_classic()+
  theme(axis.text.y= element_markdown(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20), 
        axis.title.x=element_text(face = "bold", size = 24))
#save 
pdf(file = "Forest.Plot.Tumor.COX.PH.pdf", height = 14, width = 10)
Tumor.forest.plot
dev.off()



####Lung samples#### 
#only significant results with padj < 0.2 were choosen 

res <- readxl::read_excel("COX_PH_lung.xlsx")

#match asv column from taxatable 
taxtable <- data.frame(tax_table(physeq))
#now create a column called names and paste name at Genus or species level 
taxtable$names <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",taxtable$Species,sep=""),
                               ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,sep=""),
                                      ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",sep=""),
                                             ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",sep=""),
                                                    ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",sep=""),
                                                           ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",sep=""),
                                                                  ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",sep=""), paste(rownames(taxtable))))))))))

#create a column "name ASV" that include the ASV (or OTU)
taxtable$nameASV <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",rownames(taxtable),sep=""),
                                 ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,"_",rownames(taxtable),sep=""),
                                        ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",rownames(taxtable),sep=""),
                                               ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",rownames(taxtable),sep=""),
                                                      ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",rownames(taxtable),sep=""),
                                                             ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",rownames(taxtable),sep=""),
                                                                    ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",rownames(taxtable),sep=""), paste(rownames(taxtable))))))))))


#create a column of asv which is the row name 
taxtable$asv <- rownames(taxtable)

#now merge the table with your results according to asv column 
lung.forest.df <- inner_join(res, taxtable, by="asv")
lung.forest.df$nameASV.trim <- paste0(lung.forest.df$names,paste0(".."), paste0(str_sub(lung.forest.df$asv, - 3, - 1)))
lung.forest.df$FDR <- as.numeric(lung.forest.df$FDR)

#set color vector accroding to contam
lung.forest.df <- lung.forest.df %>% 
  mutate(contaminant= factor(contaminant)) %>% 
  mutate(color=ifelse(lung.forest.df$contaminant =="TRUE", "red", "black")) %>% 
  mutate(upper_CI= as.numeric(upper_CI)) %>% 
  mutate(lower_CI= as.numeric(lower_CI)) %>% 
  mutate(HR= as.numeric(HR)) %>% 
  mutate(nameASV.trim.colored = lung.forest.df$nameASV.trim) %>% 
  filter(FDR<=0.2) %>% 
  arrange(., desc(HR)) %>% 
  dplyr::mutate(ord = order(as.numeric(rownames(.))))

#add color of potential contaminants 
lung.forest.df$nameASV.trim.colored <- paste0("<span style=\"color: ", lung.forest.df$color, "\">", lung.forest.df$nameASV.trim.colored, "</span>")

#plot only top 40 so figure can be similar to tumor figure
lung.forest.df <- lung.forest.df %>% dplyr::slice(1:40)
lung.forest.df <- lung.forest.df[-1,]

lung.forest.plot <- ggplot(data=lung.forest.df, aes(x=HR, y=fct_reorder(nameASV.trim.colored, -ord), xmin=lower_CI, xmax=upper_CI))+
  geom_point(size=5)+
  geom_errorbarh(height=0.2)+
  geom_vline(xintercept = 1, color="black", linetype="dashed", alpha=0.5)+
  ylab("")+xlab("Hazard Ratio")+
  xlim(0,3.5)+
  theme_classic()+
  theme(axis.text.y= element_markdown(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20), 
        axis.title.x=element_text(face = "bold", size = 24))
#save 
pdf(file = "Forest.Plot.lung.COX.PH.pdf", height = 14, width = 10)
lung.forest.plot
dev.off()


###### comparison of AUC generated by elastic net for every dataset (Figure 4E, 4f) ######

#tumor samples 
AUC_values <- c(0.57, 0.53, 0.6)
lower_ci <- c(0.44, 0.39, 0.47)
upper_ci <- c(0.71, 0.67, 0.73)
data_type <- c("Host", "Microbiome", "Host + Microbiome")
plot_data <- data.frame(data_type, AUC_values, lower_ci, upper_ci)
plot_data <- plot_data %>% mutate(data_type=factor(data_type, levels=c("Host + Microbiome", "Host", "Microbiome")))

pdf(file = "new_figures/forest_plot_tumor_model_mod.pdf", width = 8, height = 5)
ggplot(plot_data, aes(x=AUC_values, y=data_type, xmin=lower_ci, xmax=upper_ci))+
  geom_point(size=5)+
  geom_errorbarh(height=0.2)+
  geom_vline(xintercept = 0.5, color="red", linetype="dashed", alpha=0.5)+
  ylab("")+xlab("AUC (95% CI)")+
  theme_classic()+
  theme(axis.text.y= element_text(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20), 
        axis.title.x=element_text(face = "bold", size = 20))
dev.off()


#lung samples 
AUC_values <- c(0.8, 0.65, 0.83)
lower_ci <- c(0.71, 0.52, 0.74)
upper_ci <- c(0.89, 0.78, 0.92)
data_type <- c("Host", "Microbiome", "Host + Microbiome")
plot_data <- data.frame(data_type, AUC_values, lower_ci, upper_ci)
plot_data <- plot_data %>% mutate(data_type=factor(data_type, levels=c("Host + Microbiome", "Host", "Microbiome")))

pdf(file = "new_figures/forest_plot_lung_model_mod.pdf", width = 8, height = 5)
ggplot(plot_data, aes(x=AUC_values, y=data_type, xmin=lower_ci, xmax=upper_ci))+
  geom_point(size=5)+
  geom_errorbarh(height=0.2)+
  geom_vline(xintercept = 0.5, color="red", linetype="dashed", alpha=0.5)+
  ylab("")+xlab("AUC (95% CI)")+
  theme_classic()+
  theme(axis.text.y= element_text(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20), 
        axis.title.x=element_text(face = "bold", size = 20))
dev.off()







######## Random forest model and generation of figures 2A, B, D##########

### use data of upregualted taxa 
#read results from your edgeR analysis 
taxa_upreg <- read.csv(file = "new_results/edgeR.results_tumor_Progression_Lab_In.csv")

#select only significant and upregualted taxa 
asv_to_save <- taxa_upreg %>% filter(FDR<=0.2 & logFC>0) %>% select(asv)
#####tumor samples only #####
#use relative table 
ps.rel <- In.Lung.Tissue.pruned.rel
ps.rel <- subset_taxa(ps.rel, rownames(tax_table(ps.rel)) %in% asv_to_save$asv)




#####tumor samples only #####
#use relative table 
ps.rel <- ps.rel
#prepare taxatable 
#now get taxa names at genus or species level from taxa table 
tax_table(ps.rel)

taxtable <- data.frame(tax_table(ps.rel))

#now create a column called names and paste name at Genus level 
taxtable$names <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",taxtable$Species,sep=""),
                               ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,sep=""),
                                      ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",sep=""),
                                             ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",sep=""),
                                                    ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",sep=""),
                                                           ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",sep=""),
                                                                  ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",sep=""), paste(rownames(taxtable))))))))))

#create a column of asv / OTU that is from the row name 
taxtable$asv <- rownames(taxtable)
taxtable$nameASV.trim <- paste0(taxtable$names, paste0("..", paste0(str_sub(taxtable$asv, -3, -1))))

#prepeare contamlist 
contamlist

#run the model. change top_n every time you want to examine top taxa
set.seed(1234)

#get metadata 
metadata.for.RF <- sample_data(ps.rel) %>% 
  data.frame()

#get features (taxa) 
features <- otu_table(ps.rel)
#transpose
features.transposed<- data.frame(t(features)) 

#define outcome 
outcome <- as.factor(sample_data(ps.rel)$Progression_Lab)

#build a data frame 
rf.data_complete <- data.frame(features.transposed, outcome)

#prepare data 
#remove na 
rf.data_complete <- rf.data_complete %>% 
  drop_na()

#be sure data is numerical 
lapply(rf.data_complete, as.integer)

#replace infinite data 
rf.data_complete[mapply(is.infinite, rf.data_complete)] <- NA

rf.data_complete <- rf.data_complete %>% 
  mutate_if(is.character, as.factor)


# Split the data into training and testing sets
set.seed(1234)
train_index <- sample(nrow(rf.data_complete), 0.8 * nrow(rf.data_complete))
train_data <- rf.data_complete[train_index, ]
test_data <- rf.data_complete[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)


#get variable importance 
# Extract the mean decrease impurity (MDI) values for each feature
mdi <- randomForest::importance(fit)
mdi <- mdi %>% 
  data.frame() %>% 
  dplyr::select(MeanDecreaseGini)

# Normalize the MDI values
Gini <- mdi %>% 
  dplyr::mutate(Norm.Gini= MeanDecreaseGini/sum(MeanDecreaseGini)) %>% 
  dplyr::mutate(Gini = Norm.Gini/max(Norm.Gini)) %>% 
  dplyr::arrange(desc(Gini)) %>% 
  dplyr::select(Gini)


#add enrich group (or outcome)
features.means <- by(t(features), outcome, colMeans)
features.means <- do.call(cbind, features.means)
idx_enrich <- apply(features.means, 1, which.max)
group_enrich <- colnames(features.means)[idx_enrich]

Gini$enrich_group <- group_enrich

#leave out genes with 0 value of Gini that don't affct the model (need to change every time you run it)
Gini.100.per<- Gini %>% 
  dplyr::filter(Gini>0)

#from here subset genes of increasing percentage so you will use later to fit the model 

#top 1%
features.1.per <- Gini %>% 
  dplyr::slice(1:round(0.01*nrow(Gini.100.per)))

#top 5%
features.5.per <- Gini %>% 
  dplyr::slice(1:round(0.05*nrow(Gini.100.per)))

#top 10%
features.10.per <- Gini %>% 
  dplyr::slice(1:round(0.1*nrow(Gini.100.per)))

#top 20%
features.20.per <- Gini %>% 
  dplyr::slice(1:round(0.2*nrow(Gini.100.per)))

#top 50%
features.50.per <- Gini %>% 
  dplyr::slice(1:round(0.5*nrow(Gini.100.per)))

#top 75%
features.75.per <- Gini %>% 
  dplyr::slice(1:round(0.75*nrow(Gini.100.per)))

#Finally get one df to export 
Ginidf.to.save <- Gini %>% 
  dplyr::mutate(one_per = Gini) %>% 
  dplyr::mutate(five_per = Gini) %>% 
  dplyr::mutate(ten_per = Gini) %>% 
  dplyr::mutate(twent_per = Gini) %>% 
  dplyr::mutate(fifty_per = Gini) %>% 
  dplyr::mutate(sevent_fiv_per = Gini) %>% 
  dplyr::mutate(one_hund = Gini) %>% 
  dplyr::select(-Gini)
#replace with NA in each column 
nrow(Ginidf.to.save)
Ginidf.to.save$one_per [nrow(features.1.per) :nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$five_per [nrow(features.5.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$ten_per [nrow(features.10.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$twent_per [nrow(features.20.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$fifty_per [nrow(features.50.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$sevent_fiv_per [nrow(features.75.per):nrow(Ginidf.to.save)] <- NA

#adding taxa names and contamlist 
Ginidf.to.save$asv <- rownames(Ginidf.to.save)
Ginidf.to.save$asv <- gsub("X", "", Ginidf.to.save$asv)
Ginidf.to.save <- inner_join(Ginidf.to.save, contamlist, by="asv")
Ginidf.to.save <- inner_join(Ginidf.to.save, taxtable, by="asv")
Ginidf.to.save <- Ginidf.to.save %>% dplyr::select(names, asv, one_per, five_per, ten_per, twent_per, fifty_per, sevent_fiv_per, one_hund,
                                                   enrich_group, contaminant)
Ginidf.to.save$SN <- rownames(Ginidf.to.save)
write.csv(Ginidf.to.save, file = "new_results/RF_results_Tumor_samples_taxa_modified_upreg_edgeR.csv")

######get ROC plor and AUC
#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence
library(pROC)

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = "new_figures/roc.plot.Tumor_samples_100_per_taxa_modified_upreg_edgeR.pdf", height = 8, width = 8)
roc.plot
dev.off()



####repeat model for 1%####
n <- 1
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.1.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+  
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.Tumor_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()


##### repeat for 5 % 
n <- 5
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.5.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.Tumor_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()


#10% 
n <- 10
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.10.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.Tumor_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()



### 20% 
n <- 20
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.20.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.Tumor_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()



### 50% 
n <- 50
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.50.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.Tumor_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()



#### 75% 
n <- 75
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.75.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.Tumor_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()





#########unaffected lung subset 

#read results from your edgeR analysis 
taxa_upreg <- read.csv(file = "new_results/edgeR.results_lung_Progression_Lab_In.csv")
#filter for upregulated and significant taxa 
asv_to_save <- taxa_upreg %>% filter(FDR<=0.2 & logFC>0) %>% select(asv)
#####lung samples only #####
#use relative table 
ps.rel <- UnIn.Lung.Tissue.pruned.rel
ps.rel <- subset_taxa(ps.rel, rownames(tax_table(ps.rel)) %in% asv_to_save$asv)


#prepare taxatable 
#now get taxa names at genus or species level from taxa table 
tax_table(ps.rel)

taxtable <- data.frame(tax_table(ps.rel))

#now create a column called names and paste name at Genus level 
taxtable$names <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",taxtable$Species,sep=""),
                               ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,sep=""),
                                      ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",sep=""),
                                             ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",sep=""),
                                                    ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",sep=""),
                                                           ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",sep=""),
                                                                  ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",sep=""), paste(rownames(taxtable))))))))))

#create a column of asv / OTU that is from the row name 
taxtable$asv <- rownames(taxtable)
taxtable$nameASV.trim <- paste0(taxtable$names, paste0("..", paste0(str_sub(taxtable$asv, -3, -1))))

#prepeare contamlist 
contamlist

#run the model. change top_n every time you want to examine top taxa
set.seed(1234)

#get metadata 
metadata.for.RF <- sample_data(ps.rel) %>% 
  data.frame()

#get features (taxa) 
features <- otu_table(ps.rel)
#transpose
features.transposed<- data.frame(t(features)) 

#define outcome 
outcome <- as.factor(sample_data(ps.rel)$Progression_Lab)

#build a data frame 
rf.data_complete <- data.frame(features.transposed, outcome)

#prepare data 
#remove na 
rf.data_complete <- rf.data_complete %>% 
  drop_na()

#be sure data is numerical 
lapply(rf.data_complete, as.integer)

#replace infinite data 
rf.data_complete[mapply(is.infinite, rf.data_complete)] <- NA

rf.data_complete <- rf.data_complete %>% 
  mutate_if(is.character, as.factor)


# Split the data into training and testing sets
set.seed(1234)
train_index <- sample(nrow(rf.data_complete), 0.8 * nrow(rf.data_complete))
train_data <- rf.data_complete[train_index, ]
test_data <- rf.data_complete[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get variable importance 
# Extract the mean decrease impurity (MDI) values for each feature
mdi <- randomForest::importance(fit)
mdi <- mdi %>% 
  data.frame() %>% 
  dplyr::select(MeanDecreaseGini)

# Normalize the MDI values
Gini <- mdi %>% 
  dplyr::mutate(Norm.Gini= MeanDecreaseGini/sum(MeanDecreaseGini)) %>% 
  dplyr::mutate(Gini = Norm.Gini/max(Norm.Gini)) %>% 
  dplyr::arrange(desc(Gini)) %>% 
  dplyr::select(Gini)


#add enrich group (or outcome)
features.means <- by(t(features), outcome, colMeans)
features.means <- do.call(cbind, features.means)
idx_enrich <- apply(features.means, 1, which.max)
group_enrich <- colnames(features.means)[idx_enrich]

Gini$enrich_group <- group_enrich

#leave out genes with 0 value of Gini that don't affct the model (need to change every time you run it)
Gini.100.per<- Gini %>% 
  dplyr::filter(Gini>0)

#from here subset genes of increasing percentage so you will use later to fit the model 

#top 1%
features.1.per <- Gini %>% 
  dplyr::slice(1:round(0.01*nrow(Gini.100.per)))

#top 5%
features.5.per <- Gini %>% 
  dplyr::slice(1:round(0.05*nrow(Gini.100.per)))

#top 10%
features.10.per <- Gini %>% 
  dplyr::slice(1:round(0.1*nrow(Gini.100.per)))

#top 20%
features.20.per <- Gini %>% 
  dplyr::slice(1:round(0.2*nrow(Gini.100.per)))

#top 50%
features.50.per <- Gini %>% 
  dplyr::slice(1:round(0.5*nrow(Gini.100.per)))

#top 75%
features.75.per <- Gini %>% 
  dplyr::slice(1:round(0.75*nrow(Gini.100.per)))

#Finally get one df to export 
Ginidf.to.save <- Gini %>% 
  dplyr::mutate(one_per = Gini) %>% 
  dplyr::mutate(five_per = Gini) %>% 
  dplyr::mutate(ten_per = Gini) %>% 
  dplyr::mutate(twent_per = Gini) %>% 
  dplyr::mutate(fifty_per = Gini) %>% 
  dplyr::mutate(sevent_fiv_per = Gini) %>% 
  dplyr::mutate(one_hund = Gini) %>% 
  dplyr::select(-Gini)
#replace with NA in each column 
nrow(Ginidf.to.save)
Ginidf.to.save$one_per [nrow(features.1.per) :nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$five_per [nrow(features.5.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$ten_per [nrow(features.10.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$twent_per [nrow(features.20.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$fifty_per [nrow(features.50.per):nrow(Ginidf.to.save)] <- NA
Ginidf.to.save$sevent_fiv_per [nrow(features.75.per):nrow(Ginidf.to.save)] <- NA

#adding taxa names and contamlist 
Ginidf.to.save$asv <- rownames(Ginidf.to.save)
Ginidf.to.save$asv <- gsub("X", "", Ginidf.to.save$asv)
Ginidf.to.save <- inner_join(Ginidf.to.save, contamlist, by="asv")
Ginidf.to.save <- inner_join(Ginidf.to.save, taxtable, by="asv")
Ginidf.to.save <- Ginidf.to.save %>% dplyr::select(names, asv, one_per, five_per, ten_per, twent_per, fifty_per, sevent_fiv_per, one_hund,
                                                   enrich_group, contaminant)
Ginidf.to.save$SN <- rownames(Ginidf.to.save)
write.csv(Ginidf.to.save, file = "new_results/RF_results_lung_samples_taxa_modified_upreg_edgeR.csv")

######get ROC plor and AUC
#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = "new_figures/roc.plot.lung_samples_100_per_taxa_modified_upreg_edgeR.pdf", height = 8, width = 8)
roc.plot
dev.off()



####repeat model for 1%####
n <- 1
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.1.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.lung_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()


##### repeat for 5 % 
n <- 5
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.5.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.lung_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()


#10% 
n <- 10
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.10.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.lung_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()



### 20% 
n <- 20
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.20.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.lung_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()



### 50% 
n <- 50
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.50.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.lung_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()



#### 75% 
n <- 75
# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(features.75.per), outcome)

set.seed(1234)
train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = train_data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions
predictions <- predict(fit, newdata = test_data[,-ncol(test_data)], type = "prob")
predictions <- data.frame(predictions)

#extract prediciton of. the outcome of interset (recurrence in this case)
predictions <- predictions$Recurrence

#build ROC object 
roc <- roc(test_data[, ncol(test_data)], predictions)

library(pROC)
#using ggplot 
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(face="bold",size=18))+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

pdf(roc.plot, file = paste0("new_figures/roc.plot.lung_samples_", n, "_per_taxa_mod_upreg_edgeR.pdf"), height = 8, width = 8)
roc.plot
dev.off()




##### prepare final figures 


# tumor samples ASV level
AUC <- c(0.65, 0.67, 0.67, 0.64, 0.62, 0.56, 0.55)
perc.taxa <- c(1, 5, 10, 20, 50, 75, 100)

AUC.all.taxa.tumor <- data.frame(AUC, perc.taxa)

AUC.all.taxa.tumor.plot <- ggplot(AUC.all.taxa.tumor, aes(x=perc.taxa, y=AUC))+
  geom_line()+
  geom_point()+
  ylim(c(0.5,0.85))+
  xlab("% Top Taxa")+
  theme_classic()

pdf(file = "new_figures/AUC.tumor.samples.ASV.plot_modified.pdf", width = 6, height = 6)
AUC.all.taxa.tumor.plot
dev.off()

#all taxa lung samples 
AUC <- c(0.42, 0.42, 0.63, 0.81, 0.83, 0.82, 0.87)
AUC.all.taxa.lung <- data.frame(AUC, perc.taxa)

AUC.all.taxa.lung.plot <- ggplot(AUC.all.taxa.lung, aes(x=perc.taxa, y=AUC))+
  geom_line()+
  geom_point()+
  ylim(c(0.5,0.85))+
  xlab("% Top Taxa")+
  theme_classic()

pdf(file = "new_figures/AUC.lung.samples.ASV.plot_modified.pdf", width = 6, height = 6)
AUC.all.taxa.lung.plot
dev.off()

#Combined plot
AUC.tumor <- c(0.65, 0.67, 0.67, 0.64, 0.62, 0.56, 0.55)
AUC.lung <- c(0.42, 0.42, 0.63, 0.81, 0.83, 0.82, 0.87)

AUC.df <- data.frame(perc.taxa, AUC.tumor, AUC.lung)
AUC.df
colors <- c( "AUC.tumor" = "red", "AUC.lung" = "blue")

AUC.combined.plot <- ggplot(AUC.df, aes(x=perc.taxa))+
  geom_line(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_point(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_line(aes(y=AUC.lung, color="AUC.lung"))+
  geom_point(aes(y=AUC.lung, color="AUC.lung"))+
  ylim(c(0.4, 0.9))+
  scale_x_continuous(breaks = seq(0,100,10))+
  xlab("% Top Taxa")+
  ylab("AUC")+
  labs(color="Legend")+
  scale_color_manual(values = colors, labels = c("Lung", "Tumor"))+
  ggtitle("AUC according to sample type at ASV level")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 26, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "top")
#save 
pdf(file = "new_figures/AUC.combined.ASV.level_modified_upreg_edgeR.pdf", width = 10, height = 12)
AUC.combined.plot
dev.off()

#limit figure to start from 10% 
AUC.combined.plot <- ggplot(AUC.df, aes(x=perc.taxa))+
  geom_line(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_point(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_line(aes(y=AUC.lung, color="AUC.lung"))+
  geom_point(aes(y=AUC.lung, color="AUC.lung"))+
  ylim(c(0.4, 0.9))+
  scale_x_continuous(breaks = seq(10,100,10), limits = c(10,100))+
  xlab("% Top Taxa")+
  ylab("AUC")+
  labs(color="Legend")+
  scale_color_manual(values = colors, labels = c("Lung", "Tumor"))+
  ggtitle("AUC according to sample type at ASV level")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 26, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "top")
#save 
pdf(file = "new_figures/AUC.combined.ASV.level_modified_upreg_edgeR_limits_10_100.pdf", width = 10, height = 6)
AUC.combined.plot
dev.off()



#### best AUC results plots 

#tumor 

plot_df <- read.csv(file = "new_results/RF_results_Tumor_samples_taxa_modified_upreg_edgeR.csv")

plot_df$nameASV.trim <- paste0(plot_df$names, paste0("..", paste0(str_sub(plot_df$asv, -3, -1))))

#change color of contaminant to purple 
plot_df <- plot_df %>% dplyr::mutate(color= ifelse(contaminant=="TRUE", "red", "black"))
plot_df$nameASV.trim.colored <- paste0("<span style=\"color: ", plot_df$color, "\">", plot_df$nameASV.trim, "</span>")
plot_df <- plot_df %>% filter(ten_per!="NA")

#plot_df <- plot_df %>% dplyr::slice(1:30)

#plot
p <- ggplot(plot_df, aes(x=ten_per, y=fct_reorder(nameASV.trim.colored, ten_per, max)))+
  geom_col(position = position_dodge(0.9), fill="red", color="black", alpha = 0.7, width = 0.5)+
  ylab("")+
  xlab("Gini Index")+
  ggplot2::ggtitle("Top Taxa Achieving best AUC in Tumor Samples at ASV level (10%)")+
  theme_bw()+
  theme(axis.text.y = element_markdown(size =20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 26, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))

#save the plot 
pdf(file="new_figures/10_per_best_auc_Tumor.samples.ASV.level_modified_edgeR.pdf", height = 11, width = 10)
p
dev.off()



#lung
plot_df <- read.csv(file = "new_results/RF_results_lung_samples_taxa_modified_upreg_edgeR.csv")

plot_df$nameASV.trim <- paste0(plot_df$names, paste0("..", paste0(str_sub(plot_df$asv, -3, -1))))

#change color of contaminant to purple 
plot_df <- plot_df %>% dplyr::mutate(color= ifelse(contaminant=="TRUE", "red", "black"))
plot_df$nameASV.trim.colored <- paste0("<span style=\"color: ", plot_df$color, "\">", plot_df$nameASV.trim, "</span>")
plot_df <- plot_df %>% filter(one_hund!="NA")

#plot top 30 
plot_df <- plot_df %>% dplyr::slice(1:30)

#plot
p <- ggplot(plot_df, aes(x=one_hund, y=fct_reorder(nameASV.trim.colored, one_hund, max)))+
  geom_col(position = position_dodge(0.9), fill="blue", color="black", alpha = 0.7, width = 0.5)+
  ylab("")+
  xlab("Gini Index")+
  ggplot2::ggtitle("Top Taxa Achieving best AUC in Tumor Samples at ASV level (100%)")+
  theme_bw()+
  theme(axis.text.y = element_markdown(size =20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 26, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))

#save the plot 
pdf(file="new_figures/100_per_best_auc_lung.samples.ASV.level_modified_edgeR.pdf", height = 11, width = 12)
p
dev.off()

