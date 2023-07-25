# 05/01/2023 ; Fares Darawshy ; Segal Lab 
# Microbial and Host Transcriptomic Signatures of Early Lung Cancer 
# for github publication 

#Load/Save wd ----- 
setwd() #set your own wd
#save.image(file="16S.Early.Lung.Cancer.RData")
load(file = "16S.Early.Lung.Cancer.RData" )

###install needed packages  
#install_packages: this function take a list of packages and install them
install_packages <- function(package_list){
  list_installed <- installed.packages()
  new_pkgs <- subset(package_list, !(package_list %in% list_installed[, "Package"]))
  already_install_pkgs <- subset(package_list, (package_list %in% list_installed[, "Package"]))
  print(c(" packages already installed:",already_install_pkgs))
  print(strrep("_", 70))
  if(length(new_pkgs)!=0){
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(new_pkgs, update=F)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}


#packages needed for the script
pkg_list <- c("devtools", "tidyverse", "readr", "readtext", "vegan", "ade4", "biomformat", "cachem", "utf8", "backports", "colorspace", "rhdf5", "DelayedArray","Biobase", "Biostrings", "magick",
              "phyloseq", "ape", "phangorn", "ggpubr", "decontam", "ggplot2", "reshape2", "caret", "randomForest", "survival", "survminer", "ggforce", "scales", "pROC",
              "tidyr", "matrixStats", "DESeq2", "edgeR", "limma", "Glimma", "RColorBrewer", "pheatmap",
              "pheatmap", "ggrepel", "scales", "data.table", "fBasics", "forcats", "maptools", 
              "lubridate", "boot", "table1", "stringr", "papaja", "flextable",
              "Gmisc", "glue", "htmlTable", "grid", "magrittr", "rmarkdown", "plotly",
              "microbiome", "DT", "webshot", "lubridate", "png", "RCurl", "cowplot", "janitor",
              "optmatch", "MatchIt", "decontam", "qdap", "stringr","openxlsx", "chisq.posthoc.test", 
              "FSA", "cobalt", "ggplotify", "grid", "gridExtra", "combinat", "mixOmics", "gplots", "plyr", "readxl", "DESeq2", "mia", "microbiomeMarker", "ggpubr", "jpeg", "openxlsx")

#installing all packages needed  (run this if you need to install packages)
#install_packages(pkg_list)

#loading packages (run this to load packages)
#for (i in pkg_list){
#  eval(bquote(library(.(i))))
#}

#load libraries 
library(pacman)
pacman::p_load(qiime2R,ggpubr, dplyr, ggalt, forcats, car, magrittr, lubridate, DESeq2, edgeR, limma, Glimma, 
               RColorBrewer, pheatmap, ggplot2, gplots, ggrepel, pathfindR, scales, data.table, fBasics, forcats, 
               omu, ade4, "vegan", phyloseq, maptools, "reshape2", survival, survminer, cowplot, ggforce, scales,pROC, 
               dplyr, knitr, decontam, tidyverse, microbiomeMarker, caret, randomForest, install = FALSE, update = FALSE)

## To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}

#reading Metadata 
map = "Surgical.Cohort.Map.txt"
mapping.table = sample_data(read.table(map,header = T, sep = "\t", row.names = 1))


################################16S Microbial Analysis##################################

#creating a phyloseq object from metadata and qiime output
physeq <- qza_to_phyloseq(features = "table.filtered.qza",
                          tree = "rooted-tree_quality.filtered.qza",
                          "taxonomy.filtered.qza",
                          metadata = "Surgical.Cohort.Map.txt")

#remove all taxa with 0 abundance.
physeq = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)  

#create relative table 
physeq.rel.table= transform_sample_counts(physeq, normalizeSample)

#create a table with background, tumor and lung samples 
physeq.3 = subset_samples(physeq, Sample_Type_Involved %in% c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"))
physeq.3.rel <- transform_sample_counts(physeq.3, normalizeSample)

#create subset tables 
#a table of lower airway tissue samples 
physeq.noBKG <- subset_samples(physeq, Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"))
physeq.noBKG.rel <- transform_sample_counts(physeq.noBKG, normalizeSample)


# filter and remove taxa not seen more than *0.002* times in at least *1%* of the samples. This protects against an ASV with small mean & trivially large C.V. ----
pruned.taxa <- genefilter_sample(physeq.noBKG.rel, filterfun_sample(function(x) x> 0.002), 
                                 A = 0.01 * nsamples(physeq.noBKG.rel))
#prune relative table 
physeq.noBKG.pruned.rel <- prune_taxa(pruned.taxa, physeq.noBKG.rel)
#prune count table 
physeq.noBKG.pruned <- prune_taxa(pruned.taxa, physeq.noBKG)

#prepare Tumor only tables 
In.Lung.Tissue.pruned <- subset_samples(physeq.noBKG.pruned, sample_data(physeq.noBKG.pruned)$Sample_Type_Involved == "Lung.Tissue.In")
In.Lung.Tissue.pruned.rel <- transform_sample_counts(In.Lung.Tissue.pruned, normalizeSample)

#prepare Lung only tables 
UnIn.Lung.Tissue.pruned <- subset_samples(physeq.noBKG.pruned, sample_data(physeq.noBKG.pruned)$Sample_Type_Involved == "Lung.Tissue.UnIn")
UnIn.Lung.Tissue.pruned.rel <- transform_sample_counts(UnIn.Lung.Tissue.pruned, normalizeSample)

################################################################################
#Figure 1A
######Kaplan Meir Curve for Recurrence##########

#Load relevane packages 

#Use subset of samples with tumor only 

#We need data from our metadata that include time to recurrence, status of recurrence
metadata <- data.frame(sample_data(In.Lung.Tissue.pruned))

#decide which columns to keep: Column to be used as time (in days) - New_TTP, column to be used as status is Progression_Lab_Inv
survival.data <- metadata[, c("New_TTP","Progression", "Progression_Lab_Inv" )]

#rename the new df 
colnames(survival.data) <- c("time", "status", "Recurrence")

# turn the data into data frame so it can be used in survival function. and turn the class of columns so its either factor or numeric to be used in survival functions 
survival.data<- data.frame(survival.data)
survival.data$time <- as.numeric(survival.data$time)
survival.data$status <- as.numeric(survival.data$status)
survival.data$Recurrence <- as.factor(survival.data$Recurrence)

#create survival object
surv <- Surv(survival.data$time, survival.data$status)

#create survival df
sfit <- survfit(Surv(time, status)~Recurrence, data=survival.data)

#plot
survival.plot <- ggsurvplot(sfit, data = survival.data,
                            fun = "pct", ggtheme = theme_pubr(), 
                            conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
                            legend.labs=c("No", "Yes"), legend.title="Recurrence",  
                            palette=c("dodgerblue2", "red"), size=1,
                            title="Kaplan-Meier Curve for Recurrence Following Surgery",
                            ylab="Percent Without Recurrence",
                            xlab="Time(Days)",
                            censor = TRUE,
                            censor.shape = 124,
                            # cumevents = TRUE,
                            cumcensor= TRUE,
                            # tables.col = "strata", #color of numbers inside risk table
                            tables.height = 0.15,
                            fontsize = 4,
                            risk.table.y.text.col = T,
                            tables.theme = theme_cleantable(),
                            # risk.table.y.text = FALSE,
                            risk.table.title = "Number at risk (cumulative number of recurrence)",
                            cumcensor.title = "Cumulative number of censored subjects"
)

#save the plot 
pdf(file = "Early.Lung.Cancer/Kaplan.Meier.Plot.pdf", width = 8, height = 8)
survival.plot
dev.off()

#stats
#overall p value 
surv_pvalue(sfit)

#Cox regression analysis 
cox <- coxph(Surv(time, status)~Recurrence, data = survival.data)
summary(cox)

######################Tumor vs Lung vs BKG #####################################
# alpha diversity (Supp Figure 2A )

#calculate alpha diversity measures 
alpha.measures <- data.frame("Shannon" = estimate_richness(physeq.3, measures = "Shannon"), 
                             "Observed"= estimate_richness(physeq.3, measures = "Observed"), 
                             "Chao1"= estimate_richness(physeq.3, measures = "Chao1"),
                             "Sample_Type" = sample_data(physeq.3)$Sample_Type_Involved)
#explore your alpha.measures
alpha.measures
alpha.measures$Sample_Type <- as.factor(alpha.measures$Sample_Type)

#stats and plot for Shannon
compare_means(Shannon ~ Sample_Type, data = alpha.measures)
my_comparisons <- list( c("Lung.Tissue.In", "BKG"), c("Lung.Tissue.In", "Lung.Tissue.UnIn"), c("BKG", "Lung.Tissue.UnIn") )

#Plot Shannon
ggplot(alpha.measures, aes(x=Sample_Type, y=Shannon))+
  geom_boxplot(fill=c("grey", "red", "blue"), outlier.colour = NA)+
  geom_jitter(width = 0.2)+
  xlab("")+ylab("")+
  scale_x_discrete(labels = c("BKG", "Tumor", "Lung"))+
  stat_compare_means(mapping=aes(label=..p.adj..), comparisons = my_comparisons, method = "wilcox.test")+
  theme_bw()

#repeat for Richness (observed) to generate Supp Figure 2B

#################################################################################
# Tumor, Recu vs no Recu (similar figues generated using Prism) (Supp Figures 1A, 1B)

#calculate alpha diversity measures 
alpha.measures <- data.frame("Shannon" = estimate_richness(In.Lung.Tissue.pruned, measures = "Shannon"), 
                             "Observed"= estimate_richness(In.Lung.Tissue.pruned, measures = "Observed"), 
                             "Chao1"= estimate_richness(In.Lung.Tissue.pruned, measures = "Chao1"),
                             "Recurrence" = sample_data(In.Lung.Tissue.pruned)$Progression_Lab_Inv)
#explore your alpha.measures
alpha.measures
alpha.measures$Recurrence <- as.factor(alpha.measures$Progression_Lab_Inv)

#stats and plot for Shannon
compare_means(Shannon ~ Recurrence, data = alpha.measures)
my_comparisons <- list( c("In.Recurrence", "In.No.Recurrence"))

#Plot Shannon
ggplot(alpha.measures, aes(x=Recurrence, y=Shannon))+
  geom_boxplot(fill=c("orange", "#E0115F"), outlier.colour = NA)+
  geom_jitter(width = 0.2)+
  xlab("")+ylab("")+
  scale_x_discrete(labels = c("Recurrence", "No Recurrence"))+
  stat_compare_means(mapping=aes(label=..p.adj..), comparisons = my_comparisons, method = "wilcox.test")+
  theme_bw()

#repeat for Richness (Observed)

#### Lung, Recurrence vs no Recu (similar figues generated using Prism)

#calculate alpha diversity measures 
alpha.measures <- data.frame("Shannon" = estimate_richness(UnIn.Lung.Tissue.pruned, measures = "Shannon"), 
                             "Observed"= estimate_richness(UnIn.Lung.Tissue.pruned, measures = "Observed"), 
                             "Chao1"= estimate_richness(UnIn.Lung.Tissue.pruned, measures = "Chao1"),
                             "Recurrence" = sample_data(UnIn.Lung.Tissue.pruned)$Progression_Lab_Inv)
#explore your alpha.measures
alpha.measures
alpha.measures$Recurrence <- as.factor(alpha.measures$Progression_Lab_Inv)

#stats and plot for Shannon
compare_means(Shannon ~ Recurrence, data = alpha.measures)
my_comparisons <- list( c("UnIn.Recurrence", "UnIn.No.Recurrence"))

#Plot Shannon
ggplot(alpha.measures, aes(x=Recurrence, y=Shannon))+
  geom_boxplot(fill=c("#92A1CF", "blue"), outlier.colour = NA)+
  geom_jitter(width = 0.2)+
  xlab("")+ylab("")+
  scale_x_discrete(labels = c("Recurrence", "No Recurrence"))+
  stat_compare_means(mapping=aes(label=..p.adj..), comparisons = my_comparisons, method = "wilcox.test")+
  theme_bw()

################################Alpha diversity Rarefaction################################### 

#Supp Figure 1C

#calculation and plot 
set.seed(1234)

psdata.three.groups <- physeq.3
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

rarefaction_curve_data <- calculate_rarefaction_curves(psdata.three.groups, c('Observed'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 100))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata.three.groups)), by.x = 'Sample', by.y = 'row.names')

#plot 
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
  facet_wrap(facets = ~ Measure, scales = 'free_y')+ 
  facet_zoom(xlim=c(0,10000))+
  #scale_x_continuous(limits = c(0,10000))+
  theme_bw()


############Beta Diversity 

#Tumor vs Lung (Supp Fig 3)

#Calculate Bray distance
bray.distance = distance(physeq.3.rel, "bray")
bray.pco = dudi.pco(d = cailliez(bray.distance), scannf = FALSE, nf = 2)

#plot 
s.class(bray.pco.3$li, fac = as.factor(sample_data(physeq.3.rel)$Sample_Type_Involved), col=c("grey", "red", "blue"))

#stats
adonis2(bray.distance.3~ sample_data(physeq.3.rel)$Sample_Type_Involved)


#Tumor, Recurrence vs no Recurrence (Figure 1C)

#Calculate Bray distance
bray.distance = distance(In.Lung.Tissue.pruned.rel, "bray")
bray.pco = dudi.pco(d = cailliez(bray.distance), scannf = FALSE, nf = 2)

#plot 
s.class(bray.pco.3$li, fac = as.factor(sample_data(In.Lung.Tissue.pruned.rel)$Progression_Lab_Inv), col=c("orange", "#E0115F"))

#stats
adonis2(bray.distance.3~ sample_data(In.Lung.Tissue.pruned.rel)$Progression_Lab_Inv)


#Lung Recurrence vs no Recurrence (Figure 1D)

#Calculate Bray distance
bray.distance = distance(UnIn.Lung.Tissue.pruned.rel, "bray")
bray.pco = dudi.pco(d = cailliez(bray.distance), scannf = FALSE, nf = 2)

#plot 
s.class(bray.pco.3$li, fac = as.factor(sample_data(In.Lung.Tissue.pruned.rel)$Progression_Lab_Inv), col=c("#92A1CF", "blue"))

#stats
adonis2(bray.distance.3~ sample_data(In.Lung.Tissue.pruned.rel)$Progression_Lab_Inv)


##################################################################################
##################################################################################
################################Differential Enrichment analyses: edgeR ####################################

########Compare tumor vs Lung  (Supp Figure 5)

#define function to convert a phyloseq object into edgeR 
phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, 
              ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method = method)
  # Check for division by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data,\n         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = physeq.noBKG.pruned, group = "Sample_Type_Involved", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res <- TT@.Data[[1]]

#Reverse Directionality if you need to
res$logFC <- res$logFC*(-1)

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

#create a column 
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
RA.df.1 <- data.frame(otu_table(physeq.noBKG.pruned.rel))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Tumor <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(physeq.noBKG.pruned.rel))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Lung <- meanRA.save

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
res<- res %>% 
  mutate(category= ifelse(logFC>0, "Tumor", "Lung"))

# save results as a table 
write.csv(res, file = "edgeR.results.tumor.vs.lung_Pruned.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="Tumor", res$abundance.Tumor, res$abundance.Lung)

#plot as bubble plot 

res.bubble <- res
        
#get top 20 taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:20)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:20)
        
#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)
        
#create a color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "blue","grey"))
        
#Now add ASV name so you have last three digits of ASV string  
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))
  
#arrange data frame while ordering according to logFC and creating numbered order column       
res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))
        
        
#size based on relative abundance. Color based on category 
        
bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("blue", "red"))+
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
  geom_vline(xintercept=0, color="red",linetype="dashed")+
#save the plot         
pdf(file="edgeR.Bubble.Plot.Tumor.vs.Lung.pdf", height = 11, width = 12)
bubble.plot

##################################################################################
##################################################################################

####Tumor samples comparing recurrence vs no recurrence groups (Figure 1E)

#### use In.Lung.Tissue.pruned to compare Recurrence and no Recurrence groups in tumor samples (group = Progression_Lab_Inv) (Figure 1E)

#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = In.Lung.Tissue.pruned, group = "Progression_Lab_Inv", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res <- TT@.Data[[1]]

#Reverse Directionality if you need to
res$logFC <- res$logFC*(-1)

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

#create a column 
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
RA.df.1 <- data.frame(otu_table(In.Lung.Tissue.pruned.rel))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Tumor.rec <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(In.Lung.Tissue.pruned.rel))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Tumor.no.rec <- meanRA.save

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
res<- res %>% 
  mutate(category= ifelse(logFC>0, "Rec", "No_Rec"))

# save results as a table 
write.csv(res, file = "edgeR.results.tumor_Rec_vs_no_Rec_Pruned.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="Rec", res$abundance.Tumor.rec, res$abundance.Tumor.no.rec)

#plot as bubble plot 
res.bubble <- res

#get top 20 taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:20)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:20)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#create a color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"#E0115F", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange","grey"))

#Now add ASV name so you have last three digits of ASV string  
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

#arrange data frame while ordering according to logFC and creating numbered order column       
res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("orange", "#E0115F"))+
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
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  #save the plot         
  pdf(file="edgeR.Bubble.Plot.Tumor.Rec_vs_no_Rec.pdf", height = 11, width = 12)
bubble.plot

##################################################################################
##################################################################################

#Lung samples comparing recurrence vs no recurrence groups  (Figure 1F)

#### use UnIn.Lung.Tissue.pruned to compare Recurrence and no Recurrence groups in lung samples (group = Progression_Lab_Inv) (Figure 1F)

#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = In.Lung.Tissue.pruned, group = "Progression_Lab_Inv", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res <- TT@.Data[[1]]

#Reverse Directionality if you need to
res$logFC <- res$logFC*(-1)

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

#create a column 
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
RA.df.1 <- data.frame(otu_table(In.Lung.Tissue.pruned.rel))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Lung.rec <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(In.Lung.Tissue.pruned.rel))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Lung.no.rec <- meanRA.save

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
res<- res %>% 
  mutate(category= ifelse(logFC>0, "Rec", "No_Rec"))

# save results as a table 
write.csv(res, file = "edgeR.results.Lung_Rec_vs_no_Rec_Pruned.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="Rec", res$abundance.Lung.rec, res$abundance.Lung.no.rec)

#plot as bubble plot 
res.bubble <- res

#get top 20 taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:20)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:20)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#create a color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"blue", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "#92A1CF","grey"))

#Now add ASV name so you have last three digits of ASV string  
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

#arrange data frame while ordering according to logFC and creating numbered order column       
res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("#92A1CF", "blue"))+
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
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  #save the plot         
  pdf(file="edgeR.Bubble.Plot.Lung.Rec_vs_no_Rec.pdf", height = 11, width = 12)
bubble.plot

##################################################################################


##################################################################################
##################################################################################
################################Random Forest (Figure 2) ####################################

############################################## Tumor samples 

set.seed(1234)

#prepare taxonomy table  
tax_table(In.Lung.Tissue.pruned)
taxtable <- data.frame(tax_table(In.Lung.Tissue.pruned))

#now create a column called names and paste name at Genus level and species level if available 
taxtable$names <- paste(ifelse(!is.na(taxtable$Species), paste("g_",taxtable$Genus,"_s_",taxtable$Species,sep=""),
                               ifelse(!is.na(taxtable$Genus), paste("g_",taxtable$Genus,sep=""),
                                      ifelse(!is.na(taxtable$Family), paste("f_",taxtable$Family,"_g__",sep=""),
                                             ifelse(!is.na(taxtable$Order), paste("o_",taxtable$Order, "_f__g__",sep=""),
                                                    ifelse(!is.na(taxtable$Class), paste("c_",taxtable$Class, "_o__f__g__",sep=""),
                                                           ifelse(!is.na(taxtable$Phylum), paste("p_",taxtable$Phylum, "_c__o__f__g__",sep=""),
                                                                  ifelse(!is.na(taxtable$Kingdom), paste("k_",taxtable$Kingdom, "_p__c__o__f__g__",sep=""), paste(rownames(taxtable))))))))))


#create a column of asv / OTU that is from the row name 
taxtable$asv <- rownames(taxtable)
#create a column that has taxa name and last three digits of ASV to be used as labels 
taxtable$nameASV.trim <- paste0(taxtable$names, paste0("..", paste0(str_sub(taxtable$asv, -3, -1))))

#prepare contamlist from decontamination analysis 
contamlist 


#run the model. change top_n every time you want to examine top taxa

#Top_n = 8
set.seed(1234)

microbiomeMarker_RF_full <- run_sl(
  In.Lung.Tissue.pruned,
  group = "Progression_Lab",
  nfolds = 5,
  nrepeats = 3,
  taxa_rank = "none",
  top_n = 8,
  norm = "none",
  method = "RF", importance = "permutation", 
  transform = "identity")


#build results dataframe from combining all above dataframes
microbiomeMarker_RF_full_results <-  marker_table(microbiomeMarker_RF_full) %>% 
  data.frame() %>% 
  dplyr::rename(., asv = feature) %>% 
  inner_join(., taxtable, by= "asv") %>% 
  inner_join(., contamlist, by="asv") %>% 
  dplyr::select(asv, names, nameASV.trim, ef_imp, enrich_group, color, 
                Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%   
  mutate(potential.contaminant = ifelse(color=="black", "FALSE", "TRUE"))

#export as a table 
write.csv(microbiomeMarker_RF_full_results, file = "RF.Tumor.Samples_Pruned_1_Per_8.csv")

# plot ROC 
p <- plot_sl_roc(microbiomeMarker_RF_full, group = "Progression_Lab") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")+  # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")

#save it 
pdf(file = "RF.Tumor.Samples_Pruned_1_Per_8.pdf", width = 8, height = 8)
p
dev.off()


####repeat the above using different top_n
#top_n to be used: 1%, 5%, 10%, 20%, 50%, 75%, and 100% of taxa N 
# taxa N = 776, so topn will be: 8, 39, 78, 155, 388, 582, 776


################################################################################


#repeat the above for Lung samples using UnIn.Lung.Tissue.pruned table: 

#Top_n = 8
set.seed(1234)

microbiomeMarker_RF_full <- run_sl(
  UnIn.Lung.Tissue.pruned,
  group = "Progression_Lab",
  nfolds = 5,
  nrepeats = 3,
  taxa_rank = "none",
  top_n = 8,
  norm = "none",
  method = "RF", importance = "permutation", 
  transform = "identity")


#build results dataframe from combining all above dataframes
microbiomeMarker_RF_full_results <-  marker_table(microbiomeMarker_RF_full) %>% 
  data.frame() %>% 
  dplyr::rename(., asv = feature) %>% 
  inner_join(., taxtable, by= "asv") %>% 
  inner_join(., contamlist, by="asv") %>% 
  dplyr::select(asv, names, nameASV.trim, ef_imp, enrich_group, color, 
                Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%   
  mutate(potential.contaminant = ifelse(color=="black", "FALSE", "TRUE"))

#export as a table 
write.csv(microbiomeMarker_RF_full_results, file = "RF.Lung.Samples_Pruned_1_Per_8.csv")

# plot ROC 
p <- plot_sl_roc(microbiomeMarker_RF_full, group = "Progression_Lab") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")+  # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")

#save it 
pdf(file = "RF.Lung.Samples_Pruned_1_Per_8.pdf", width = 8, height = 8)
p
dev.off()

####repeat the above using different top_n
#top_n to be used: 1%, 5%, 10%, 20%, 50%, 75%, and 100% of taxa N 
# taxa N = 776, so topn will be: 8, 39, 78, 155, 388, 582, 776


################################################################################


####Figure 2A####
# From the plots generated above, get the AUC value for each top_n (or use the following code)

# tumor samples ASV level
AUC <- c(0.75, 0.73, 0.73, 0.68, 0.62, 0.58, 0.53)
AUC.all.taxa.tumor <- data.frame(AUC, perc.taxa)

# Lung samples ASV level
AUC <- c(0.66, 0.77, 0.78, 0.73, 0.7, 0.67, 0.63)
AUC.all.taxa.lung <- data.frame(AUC, perc.taxa)

#Combined plot
AUC.tumor <- c(0.75, 0.73, 0.73, 0.68, 0.62, 0.58, 0.53)
AUC.lung <- c(0.66, 0.77, 0.78, 0.73, 0.7, 0.67, 0.63)

AUC.df <- data.frame(perc.taxa, AUC.tumor, AUC.lung)
AUC.df
colors <- c( "AUC.tumor" = "red", "AUC.lung" = "blue")

#plot AUC according to Gini index (Figure 2A)
ggplot(AUC.df, aes(x=perc.taxa))+
  geom_line(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_point(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_line(aes(y=AUC.lung, color="AUC.lung"))+
  geom_point(aes(y=AUC.lung, color="AUC.lung"))+
  ylim(c(0.5, 0.8))+
  scale_x_continuous(breaks = seq(0,100,10))+
  xlab("% Top Taxa")+
  ylab("AUC")+
  labs(color="Legend")+
  scale_color_manual(values = colors, labels = c("Lung", "Tumor"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 26, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "top")



#####plot best AUC among each category of samples. #Figure 2B

#build a dataframe out of Random forest CSV output (or just pull the results above if saved in specific vector)
RF_Tumor.samples_ASV_best_auc <- read.csv("RF.Tumor.Samples_Pruned_1_Per.csv")

#create Taxa name wih ASV
RF_Tumor.samples_ASV_best_auc$nameASV.trim.colored <- paste0("<span style=\"color: ", RF_Tumor.samples_ASV_best_auc$color, "\">", RF_Tumor.samples_ASV_best_auc$nameASV.trim, "</span>")

#create Gini index column that will be in 0-1 range 
RF_Tumor.samples_ASV_best_auc<- RF_Tumor.samples_ASV_best_auc %>% 
  dplyr::mutate(gini.index = ef_imp/100)

#plot barplot 
ggplot(RF_Tumor.samples_ASV_best_auc, aes(x=gini.index, y=fct_reorder(nameASV.trim.colored, gini.index, max)))+
  geom_col(position = position_dodge(0.9), fill="red", color="black", alpha = 0.7)+
  ylab("")+
  xlab("Gini Index")+
  ggplot2::ggtitle("Top Taxa Achieving best AUC in Tumor Samples at ASV level")+
  theme_bw()+
  theme(axis.text.y = element_markdown(size =20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 26, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))

################################################################################

####repeat for lung samples (Figure 2C)

#build a df out of csv output (or just pull the results above if saved in specific vector)
RF_Lung.samples_ASV_best_auc <- read.csv("RF.Lung.Samples_Pruned_10_Per_78.csv")

#build the df using Taxa column, Gini index, and % of the taxa that have best AUC. 
# create colored nameASV column 
RF_Lung.samples_ASV_best_auc$nameASV.trim.colored <- paste0("<span style=\"color: ", RF_Lung.samples_ASV_best_auc$color, "\">", RF_Lung.samples_ASV_best_auc$nameASV.trim, "</span>")

#create Gini index column that will be in 0-1 range 
RF_Lung.samples_ASV_best_auc<- RF_Lung.samples_ASV_best_auc %>% 
  dplyr::mutate(gini.index = ef_imp/100)

##plot barplot 
ggplot(RF_Lung.samples_ASV_best_auc, aes(x=gini.index, y=fct_reorder(nameASV.trim.colored, gini.index, max)))+
  geom_col(position = position_dodge(0.9), fill="dodgerblue", color="black", alpha = 0.7, width = 0.5)+
  ylab("")+
  xlab("Gini Index")+
  ggplot2::ggtitle("Top Taxa Achieving best AUC in Lung Samples at ASV level")+
  theme_bw()+
  theme(axis.text.y = element_markdown(size =20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 26, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))


################################################################################
################################################################################
################################################################################

#############Background Identification Analysis ( Supp Fig 4)###############

#define needed functions used in this analysis 
####wilcox_contam_test: Wilcoxon comparsion used to run decontam objects 
wilcox_contam_test <- function(relab,neg,threshold){
  p.value <- apply(relab,1,function(x){return(wilcox.test(as.numeric(as.character(x[!neg])),
                                                          as.numeric(as.character(x[neg])),
                                                          alternative="less",exact=FALSE)$p.value)})
  fdr <- unlist(lapply(p.value,function(x){return(p.adjust(x, method="BH", n = length(p.value)))}))
  contam <- as.data.frame(cbind(p.value, fdr, ifelse(p.value<threshold, 'TRUE', 'FALSE')))
  rownames(contam) <- rownames(relab)
  colnames(contam) <- c("p.value","fdr","contaminant")
  return(contam)
}

####decontaminant_KW: identify decontaminant based on different method
decontaminant_KW<-function(input_phyloseq, 
                           SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                           sample_type_var_name,                              
                           sample_types=list(), 
                           sample_type_color=list(), 
                           sample_type_color_2nd=list(), 
                           negative_sample_type, ###the sample type that you want to be as negative control
                           compare_type=list(), 
                           stat_option=c("mean", "median"), ### the statistics to determine rank for the boxplot
                           display_contam_method=c("MannWhit","preval","freq","combin","none"), #if none is selected then there will be no red label
                           graph_option=c("boxplot","mean_SE", "mean_SD"), 
                           test_threshold,
                           output_suffix) {
  ################################################################################################################################################################################
  ##############################################################setting up all the function argument and data frames############################################################## 
  ################################################################################################################################################################################
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  display_contam_method<- match.arg(display_contam_method)
  print(stat_option)
  
  if (missing(output_suffix)) { output_suffix<-"OP"} #if output_suffix is missing, then will just give a suffix _OP at the end of the file
  if (missing(graph_option)) { graph_option<-"boxplot"} #if graph_option is missing, then will just have boxplot as default
  
  ##making sure sample_types and sample_type_color and sample_type_color_2nd has same number of elements
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #if sample_type_color_2nd is not missing then make sure it has same number of element as sample_type_color
  if (!missing(sample_type_color_2nd)) {
    stopifnot("sample_types need to have same number of elements as sample_type_color_2nd"= length(sample_types)==length(sample_type_color_2nd))
  } #if graph_option is missing, then will just have boxplot as default
  if (missing(sample_type_color_2nd)) {
    sample_type_color_2nd<-rep("black",length(sample_types)) #make the secondary color to be black by default, unless specify
  } 
  
  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
    assign(paste0("sample_type_color_2nd",i), sample_type_color_2nd[[i]])
    assign(paste0("sampletype",i), sample_types[[i]])
  }  
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <<- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))
  
  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #this is to make taxa name with rank order in them
  for (i in 1:length(sample_types)) { 
    j <- i+length(sample_types)
    k <- i+2*length(sample_types)  
    column2<-top.relab.stats[,j]
    top.relab.stats[k] <- paste0(rownames(top.relab.stats)," ","(",column2,")")
    colnames(top.relab.stats)[k] <- c(paste0(paste0("rank.order.",sample_types[i]),".name"))
  }  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  
  #top.relab.stats now has the summary statistics for each sample type and the corresponding rank order. each row is a specific taxa 
  for (rank_sample_type in sample_types) { #this determines the rank order 
    print(paste0("ranking:", rank_sample_type))
    l<-match(rank_sample_type,sample_types) ###look at the element in rank_sample_type to see what their index at sample_types
    r<-l+length(sample_types) 
    # Re-rank based upon the rank group decreasing from highest to lowest
    column3<-top.relab.stats[,r] ###r would be the column number which consider the rank order for the ith sample_types
    top.relab.stats <- top.relab.stats[order(column3, decreasing = F), , drop = FALSE ]
    top.2 <- top.relab.stats
    top.3 <- head(top.2, 100) #top.3 would contain the top 100 taxa by the rank order of ith sample_types 
    
    # Filtering ONLY the top 100 taxa 
    relab.subset.top <- relab.subset %>% 
      dplyr::filter(rownames(relab.subset) %in% rownames(top.3))
    
    # Re ordered and saved into a new group (previously filtered)
    #relab.subset.top.match has same order as top.3; otherwise relab.subset.top.match is same as relab.subset.top
    relab.subset.top.match <- relab.subset.top[match(rownames(top.3),rownames(relab.subset.top)),]  
    
    for (h in 1:length(compare_type)) {
      name_compare_type<-get(paste0("compare_type",h))
      print(paste0(negative_sample_type," compare:"))
      print(name_compare_type)
      AND<-grepl("_and_",name_compare_type, fixed=T)
      OR<-grepl("_or_",name_compare_type, fixed=T)
      COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
      
      contam<-list()
      
      if ((AND==TRUE) | (OR==TRUE)) {
        if (AND==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
        }else if (OR==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
        }
        #loop through the two sample types 
        for (i in 1:length(compare_sublist)) { 
          assign(paste0("compare_sublist",i), compare_sublist[i])
        } 
        
        for (b in 1:length(compare_sublist)) {
          # only keeping the compare list  
          selected_sample_type<-get(paste0("compare_sublist",b))
          p<-match(selected_sample_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[p]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          
          # Remove unused from conc, in this case libsizes 
          libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
          
          counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
          relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
          #remove the new column "Row.names" created when merging by 0
          rownames(counts.subset_comb) <- counts.subset_comb$Row.names
          rownames(relab.subset_comb) <- relab.subset_comb$Row.names
          counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
          relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
          
          negative_comb<- append(comparison.list_select,comparison.list_negative)
          libsizes_comb <- append(libsizes_select,libsizes_negative )
          
          contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="prevalence",
                                       neg=negative_comb,
                                       threshold=test_threshold)
          contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="frequency",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
          contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                           method="combined",
                                           neg=negative_comb,
                                           conc=libsizes_comb,
                                           threshold=test_threshold)
          contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                              neg=negative_comb,
                                              threshold=test_threshold)
          contam_dataframe<- contam %>% data.frame #compress the lists from the different results
          contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
          contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
          assign(paste0("contam_result",b),contam_result)
        }
        #combine the results for the two sample types
        contam_result_final<-merge(contam_result1,contam_result2, by=0)
        rownames(contam_result_final) <- contam_result_final$Row.names
        contam_result_final <- subset(contam_result_final, select=-c(Row.names))
        if (AND==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }else if (OR==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }
      }else{
        if (COMBINE==TRUE) { ###if compare_type is a combination set 
          ###this works if A_combine_B... it basically combine subset A and subset B
          combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
          f<-match(combine_list,sample_types) 
          comparison.list_select<- unlist(comparison.list[f]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
          #Setting up the dataset so that it include both sample types that are combine =
          libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
        } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
          if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
            f<-match(name_compare_type,sample_types) 
            comparison.list_select<- unlist(comparison.list[f]) 
            comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
            #keeping only the sample type of interest
            libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
            counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
            relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
          } else { stop(paste0(paste0(name_compare_type," must be in sample type list: "),sample_types)) }
          
        }
        #combining with the negative control group
        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
        #remove the new column "Row.names" created when merging by 0
        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
        
        negative_comb<- append(comparison.list_select,comparison.list_negative)
        libsizes_comb <- append(libsizes_select,libsizes_negative )
        
        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="prevalence",
                                     neg=negative_comb,
                                     threshold=test_threshold)
        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="frequency",
                                     neg=negative_comb,
                                     conc=libsizes_comb,
                                     threshold=test_threshold)
        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="combined",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                            neg=negative_comb,
                                            threshold=test_threshold)
        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
        
        contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
      lowerbound<-length(sample_types)+1
      upperbound<-2*length(sample_types)
      rankorder <- subset(top.relab.stats[c(lowerbound:upperbound)]) ##keeping the rank order only
      csv_output_contam <- merge(contam_result_final,rankorder,by="row.names")
      
      contam_result_final_t<-contam_result_final
      #adding the contaminant result back to the relative abundance table (top.3 which contains the top 100 taxa)
      if (display_contam_method=="MannWhit"){
        contam_result_final_t$contaminant <- contam_result_final_t$wilcox.contaminant  
      }
      else if (display_contam_method=="preval") {
        contam_result_final_t$contaminant <- contam_result_final_t$prev.contaminant
      }
      else if (display_contam_method=="freq") {
        contam_result_final_t$contaminant <- contam_result_final_t$freq.contaminant
      }
      else if (display_contam_method=="combin") {
        contam_result_final_t$contaminant <- contam_result_final_t$combined.contaminant
      }
      else if (display_contam_method=="none") {
        contam_result_final_t$contaminant <- "FALSE" #if display_contam_method is none then turn all the contaminant value to FALSE so all the label would be black
      }
      contam_result_final_t<-subset(contam_result_final_t, select=c(contaminant))
      top.4<-merge(top.3,contam_result_final_t, by=0, keep.x=T)
      rownames(top.4) <- top.4$Row.names
      top.4 <- subset(top.4, select=-c(Row.names))
      top.4$color <- ifelse(top.4$contaminant, "red", "black")
      
      for (a in sample_types) {
        assign(paste0("rank_order_",a),subset(top.4, select=c(paste0("rank.order.",a), "color", paste0(paste0("rank.order.",a),".name") )))
      }      
      
      # Transposition of the figure, maybe at this point replace the column names 
      relab.subset.top.match.transpose <- as.data.frame(t(relab.subset.top.match))
      
      # name a new column with categories 
      relab.subset.top.match.transpose$category <- rownames(relab.subset.top.match.transpose)
      
      for (a in sample_types) {
        relab.subset.top.match.transpose$category[grepl(a,relab.subset.top.match.transpose$category)] <- a
      }  
      
      # Reverse column order 
      relab.subset.top.match.transpose.reverse <- relab.subset.top.match.transpose[,order(ncol(relab.subset.top.match.transpose):1)]
      for (a in sample_types) {
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(relab.subset.top.match.transpose.reverse, category == c(a) ))
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(get(paste0("relab.subset.top.match.transpose.reverse.",a)), select=-c(category)))
      }  
      
      ###results plot
      plot.df <- as.data.frame(rbind(cbind("Mann-Whitney U test",rownames(contam_result_final)[contam_result_final$wilcox.contaminant == T]),
                                     cbind("decontam prevalence",rownames(contam_result_final)[contam_result_final$prev.contaminant == T]),
                                     cbind("decontam frequency",rownames(contam_result_final)[contam_result_final$freq.contaminant == T]),
                                     cbind("decontam combined",rownames(contam_result_final)[contam_result_final$combined.contaminant == T])))
      
      ### Plotting the results ### 
      # plot_contam_plot 
      top.4<-top.4[order(top.4[[paste0("rank.order.",rank_sample_type)]]),]
      plot.df <- plot.df[plot.df[,2] %in% rownames(top.3),]
      
      plot.df[,1] <- factor(plot.df[,1],levels=c("Mann-Whitney U test","decontam prevalence","decontam frequency","decontam combined"))
      plot.df[,2] <- factor(plot.df[,2],levels=rev(rownames(top.3)))
      plot1 <- ggplot()+
        geom_point(mapping=aes(x=!!plot.df[,1],y=!!plot.df[,2],color=!!plot.df[,1]),size=2)+
        guides(color = "none")+
        theme_bw()+  
        scale_color_discrete(drop=FALSE) +
        scale_y_discrete(drop = FALSE, labels = function(y) {
          y_org<-y
          is_long <- nchar(y) > 35
          y[is_long] <- paste0(substr(sapply(strsplit(y[is_long],".",fixed=T),"[[",4),1,10),
                               paste0("..",paste0(substr(sapply(strsplit(y[is_long],".",fixed=T),"[[",5),1,10),".."),
                                      paste0(substr(sapply(strsplit(y[is_long],".",fixed=T),"[[",6),1,10),
                                             paste0("..",str_sub(y[is_long],-4,-1)))))
          is_short <- nchar(y) < 17
          y[is_short] <- paste0(paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",4),1,5),
                                       paste0("..",paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                              paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",6),1,6),
                                                     paste0("..",str_sub(y_org[is_short],-4,-1))))))
          y
        }
        ) +
        scale_x_discrete(drop = FALSE)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      myplots <- list()  # new empty list
      myplots[[1]]<- plot1
      
      ###looping thro each sample type 
      number<-1
      number2<-2
      testing1<<-relab.subset.top.match.transpose.reverse.BKG
      testing2<<-relab.subset.top.match.transpose.reverse.Upper
      testing3<<-relab.subset.top.match.transpose.reverse.Lower
      
      for (a in sample_types) {
        #this gets relab.subset.top.match.transpose.reverse.____ from wide to long 
        temp <- get(paste0("relab.subset.top.match.transpose.reverse.",a)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
        temp <- merge(temp, top.4, by.x="id", by.y="row.names")
        temp$id2<- temp[[paste0(paste0("rank.order.",a),".name")]]
        
        p <- temp %>%   mutate(id2 = fct_reorder(id2, -get(paste0("rank.order.",rank_sample_type)))) %>%
          ggplot(., aes(x=id2, y=value)) + 
          ylab("Relative Abundance")+
          scale_x_discrete(
            labels = function(x) {
              x_org<-x
              is_long <- nchar(x) > 35
              x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                   paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                          paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                                 paste0("..",str_sub(x[is_long],-4,-1)))))
              is_short <- nchar(x) < 17
              x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                    paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                           paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                         paste0("..",str_sub(x_org[is_short],-4,-1))))))
              x
            }
          ) +
          coord_flip() +
          theme_bw() + 
          theme(axis.text.y=element_text(color=rev(top.4$color)),axis.title.y=element_blank()) +
          theme(axis.text.x=element_text(angle=90,hjust=1))+
          ggtitle(a)
        if (graph_option=="boxplot") {
          p<-p+geom_boxplot(color=sample_type_color[[number]], alpha=0.2) 
        } else if(graph_option=="mean_SD") {
          p<-p+geom_jitter(color=sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
            stat_summary(fun= mean, 
                         geom="pointrange", 
                         fun.max = function (x) mean(x)+sd(x),
                         fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                         color=sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
        } else if(graph_option=="mean_SE") {
          p<-p+geom_jitter(color=sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
            stat_summary(fun= mean, 
                         geom="pointrange", 
                         fun.max = function (x) mean(x)+sd(x),
                         fun.min = function (x) ifelse(mean(x)-sd(x)/sqrt(length(x)) < 0, 0, mean(x)-sd(x)/sqrt(length(x))),
                         color=sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
        }
        myplots[[number2]] <- p
        number<-number+1
        number2<-number2+1
      }
      
      ### GGARRANGE ALL PLOTS ### 
      plot.all <- ggarrange(plotlist=myplots, ncol=length(myplots), nrow=1, align="h")
      
      pdf_output=paste0(paste0(paste0(paste0(paste0("contaminant",paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("_rank_",rank_sample_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".pdf"))
      pdf(file=pdf_output, width=20+(length(sample_types)-3)*5, height=10)
      show(plot.all)
      dev.off()        
      
      ####export csv
      csv_output=paste0(paste0(paste0(paste0("contaminant",paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".csv"))
      xlsx_output=paste0(paste0(paste0(paste0("contaminant",paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".xlsx"))
      write.csv(csv_output_contam,csv_output, row.names = FALSE)
      write.xlsx(csv_output_contam, xlsx_output, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
    }
  }
}

####decontaminant_testing_KW: identify decontaminant based on a particular method at a range of test threshold level 
decontaminant_testing_KW <- function (input_phyloseq, 
                                      SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                                      sample_type_var_name, 
                                      sample_types=list(), 
                                      negative_sample_type, ###the sample type that you want to be as negative control
                                      compare_type=list(), 
                                      method_type=c("MannWhit","preval","freq","combin"),
                                      stat_option=c("mean", "median"), ### the statistics to determine rank for the boxplot
                                      test_threshold_all_level=list(), ### test_threshold_all_level=c(0.1,0.3,0.5,0.7,0.9) would test for the level specified and create an excel spread sheet with the parameters
                                      output_excel)  { 
  
  ################################################################################################################################################################################
  ##############################################################setting up all the function argument and data frames############################################################## 
  ################################################################################################################################################################################
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  method_type <- match.arg(method_type)
  print(stat_option)
  print(method_type)
  
  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sampletype",i), sample_types[[i]])
  }
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <<- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))
  
  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  ###keeping only the statistics for each sample type
  upperlimitchart<-length(sample_types)*2
  top.relab.stats_table <- top.relab.stats[c(1:upperlimitchart)]
  relab.subset_final<-relab.subset
  
  for (h in 1:length(compare_type)) {
    name_compare_type<-get(paste0("compare_type",h))
    print(paste0(negative_sample_type," compare:"))
    print(name_compare_type)
    AND<-grepl("_and_",name_compare_type, fixed=T)
    OR<-grepl("_or_",name_compare_type, fixed=T)
    COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
    contam<-list()
    top.relab.stats_table_new<-top.relab.stats_table
    for (test_threshold in test_threshold_all_level) {
      if ((AND==TRUE) | (OR==TRUE)) {
        if (AND==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
        }else if (OR==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
        }
        #loop through the two sample types 
        for (i in 1:length(compare_sublist)) { 
          assign(paste0("compare_sublist",i), compare_sublist[i])
        } 
        
        for (b in 1:length(compare_sublist)) {
          # only keeping the compare list  
          selected_sample_type<-get(paste0("compare_sublist",b))
          p<-match(selected_sample_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[p]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          
          # Remove unused from conc, in this case libsizes 
          libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
          
          counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
          relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
          #remove the new column "Row.names" created when merging by 0
          rownames(counts.subset_comb) <- counts.subset_comb$Row.names
          rownames(relab.subset_comb) <- relab.subset_comb$Row.names
          counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
          relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
          
          negative_comb <- append(comparison.list_select,comparison.list_negative)
          libsizes_comb <- append(libsizes_select,libsizes_negative )
          
          contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="prevalence",
                                       neg=negative_comb,
                                       threshold=test_threshold)
          contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="frequency",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
          contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                           method="combined",
                                           neg=negative_comb,
                                           conc=libsizes_comb,
                                           threshold=test_threshold)
          contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                              neg=negative_comb,
                                              threshold=test_threshold)
          contam_dataframe<- contam %>% data.frame #compress the lists from the different results
          contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
          contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
          assign(paste0("contam_result",b),contam_result)
        }
        #combine the results for the two sample types
        contam_result_final<-merge(contam_result1,contam_result2, by=0)
        rownames(contam_result_final) <- contam_result_final$Row.names
        contam_result_final <- subset(contam_result_final, select=-c(Row.names))
        if (AND==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }else if (OR==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }
      }else{
        if (COMBINE==TRUE) { ###if compare_type is a combination set 
          ###this works if A_combine_B... it basically combine subset A and subset B
          combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
          f<-match(combine_list,sample_types) 
          comparison.list_select<- unlist(comparison.list[f]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
          #Setting up the dataset so that it include both sample types that are combine =
          libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
        } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
          
          if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
            f<-match(name_compare_type,sample_types) 
            comparison.list_select<- unlist(comparison.list[f]) 
            comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
            #keeping only the sample type of interest
            libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
            counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
            relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
          } else { stop(paste0(paste0(name_compare_type," must be in sample type list:"),sample_types)) }
          
        }
        #combining with the negative control group
        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
        #remove the new column "Row.names" created when merging by 0
        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
        
        negative_comb<- append(comparison.list_select,comparison.list_negative)
        libsizes_comb <- append(libsizes_select,libsizes_negative )
        
        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="prevalence",
                                     neg=negative_comb,
                                     threshold=test_threshold)
        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="frequency",
                                     neg=negative_comb,
                                     conc=libsizes_comb,
                                     threshold=test_threshold)
        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="combined",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                            neg=negative_comb,
                                            threshold=test_threshold)
        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
        
        contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
      ### Come back from the DECONTAM section to add this in ###    
      # Adding contam data (saved to the end to add)
      ###figure out which taxa to highlight in red based on which test you choice 
      if (method_type == "MannWhit") {
        contam_result_final$contamiant <- contam_result_final$wilcox.contaminant 
        print("MannWhit")
      }
      else if (method_type == "preval") {
        contam_result_final$contamiant <- contam_result_final$prev.contaminant 
        print("preval")
      }
      else if (method_type == "freq") {
        contam_result_final$contamiant <- contam_result_final$freq.contaminant 
        print("freq")
      }
      else if (method_type == "combin") {
        contam_result_final$contamiant <- contam_result_final$combined.contaminant
        print("combin")
      }
      contam_result_final <- subset(contam_result_final, select=c("contamiant"))
      print(test_threshold)
      colnames(contam_result_final)<- paste0("level_",test_threshold)
      top.relab.stats_table_new<-merge(top.relab.stats_table_new,contam_result_final, by="row.names")
      row.names(top.relab.stats_table_new)<-top.relab.stats_table_new$Row.names #previous merge step generate a new column call Row.names
      top.relab.stats_table_new<- subset(top.relab.stats_table_new,select=-c(Row.names))
    }
    
    ####export csv
    csv_output_name<-paste0(output_excel,paste0(paste0("contam",paste0(paste0('_negC_',negative_sample_type),paste0("_compare_",name_compare_type))),paste0(paste0("_method_",method_type),"_all_level.csv")))
    write.csv(top.relab.stats_table_new,csv_output_name)
    
  }
}

####decontaminant_sensitivity_KW: 
decontaminant_sensitivity_KW <- function (input_phyloseq, 
                                          SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                                          sample_type_var_name, 
                                          sample_types=list(), 
                                          negative_sample_type, ###the sample type that you want to be as negative control
                                          compare_type=list(), 
                                          stat_option=c("mean", "median"), ### the statistics to determine rank for the boxplot
                                          test_threshold_all_level=list(), ### test_threshold_all_level=c(0.1,0.3,0.5,0.7,0.9) would test for the level specified and create an excel spread sheet with the parameters
                                          expected_contam_taxa=list(), ### list of taxa you expect to be contaminant
                                          expected_NOT_contam_taxa=list(), ### list of taxa you do NOT expect to be contaminant
                                          output_excel)  { 
  
  ################################################################################################################################################################################
  ##############################################################setting up all the function argument and data frames############################################################## 
  ################################################################################################################################################################################
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  print(stat_option)
  
  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sampletype",i), sample_types[[i]])
  }
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <<- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))
  
  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  ###keeping only the statistics for each sample type
  upperlimitchart<-length(sample_types)*2
  top.relab.stats_table <- top.relab.stats[c(1:upperlimitchart)]
  relab.subset_final<-relab.subset
  
  for (h in 1:length(compare_type)) {
    name_compare_type<-get(paste0("compare_type",h))
    print(paste0(negative_sample_type," compare:"))
    print(name_compare_type)
    AND<-grepl("_and_",name_compare_type, fixed=T)
    OR<-grepl("_or_",name_compare_type, fixed=T)
    COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
    contam<-list()
    top.relab.stats_table_new<-top.relab.stats_table
    test_num<-1
    for (test_threshold in test_threshold_all_level) {
      if ((AND==TRUE) | (OR==TRUE)) {
        if (AND==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
        }else if (OR==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
        }
        #loop through the two sample types 
        for (i in 1:length(compare_sublist)) { 
          assign(paste0("compare_sublist",i), compare_sublist[i])
        } 
        
        for (b in 1:length(compare_sublist)) {
          # only keeping the compare list  
          selected_sample_type<-get(paste0("compare_sublist",b))
          p<-match(selected_sample_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[p]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          
          # Remove unused from conc, in this case libsizes 
          libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
          
          counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
          relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
          #remove the new column "Row.names" created when merging by 0
          rownames(counts.subset_comb) <- counts.subset_comb$Row.names
          rownames(relab.subset_comb) <- relab.subset_comb$Row.names
          counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
          relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
          
          negative_comb <- append(comparison.list_select,comparison.list_negative)
          libsizes_comb <- append(libsizes_select,libsizes_negative )
          
          #generate output in global environment to check
          assign(paste0("counts.checking",b),counts.subset_comb, envir=.GlobalEnv)
          assign(paste0("relab.checking",b),relab.subset_comb,)
          assign(paste0("negative_comb.checking",b),negative_comb,envir=.GlobalEnv)
          assign(paste0("libsizes_comb.checking",b),libsizes_comb,envir=.GlobalEnv)
          
          contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="prevalence",
                                       neg=negative_comb,
                                       threshold=test_threshold)
          contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="frequency",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
          contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                           method="combined",
                                           neg=negative_comb,
                                           conc=libsizes_comb,
                                           threshold=test_threshold)
          contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                              neg=negative_comb,
                                              threshold=test_threshold)
          contam_dataframe<- contam %>% data.frame #compress the lists from the different results
          contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
          contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
          assign(paste0("contam_result",b),contam_result)
        }
        #combine the results for the two sample types
        contam_result_final<-merge(contam_result1,contam_result2, by=0)
        rownames(contam_result_final) <- contam_result_final$Row.names
        contam_result_final <- subset(contam_result_final, select=-c(Row.names))
        if (AND==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }else if (OR==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }
      }else{
        if (COMBINE==TRUE) { ###if compare_type is a combination set 
          ###this works if A_combine_B... it basically combine subset A and subset B
          combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
          f<-match(combine_list,sample_types) 
          comparison.list_select<- unlist(comparison.list[f]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
          #Setting up the dataset so that it include both sample types that are combine =
          libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
        } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
          
          if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
            f<-match(name_compare_type,sample_types) 
            comparison.list_select<- unlist(comparison.list[f]) 
            comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
            #keeping only the sample type of interest
            libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
            counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
            relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
          } else { stop(paste0(paste0(name_compare_type," must be in sample type list:"),sample_types)) }
          
        }
        #combining with the negative control group
        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
        #remove the new column "Row.names" created when merging by 0
        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
        
        negative_comb<- append(comparison.list_select,comparison.list_negative)
        libsizes_comb <- append(libsizes_select,libsizes_negative)
        
        #generate output in global environment to check
        counts.checking<<-counts.subset_comb
        relab.checking<<-relab.subset_comb
        negative_comb.checking<<-negative_comb
        libsizes_comb.checking<<-libsizes_comb
        
        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="prevalence",
                                     neg=negative_comb,
                                     threshold=test_threshold)
        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="frequency",
                                     neg=negative_comb,
                                     conc=libsizes_comb,
                                     threshold=test_threshold)
        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="combined",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                            neg=negative_comb,
                                            threshold=test_threshold)
        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
        
        contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
      
      #### add the test level as suffix to the contaminant variables
      colnames(contam_result_final)<-paste(colnames(contam_result_final),test_threshold,sep="_")
      assign(paste0("contam_result_final",test_num),contam_result_final)
      test_num<-test_num+1
    }
    
    #merge all of the contaminant results with different testing level
    contam_result_final_combine<-contam_result_final1
    for (u in 2:length(test_threshold_all_level)) {
      contam_result_final_combine<-merge(contam_result_final_combine,get(paste0("contam_result_final",u)),by="row.names")
      row.names(contam_result_final_combine)<-contam_result_final_combine$Row.names #previous merge step generate a new column call Row.names
      contam_result_final_combine<- subset(contam_result_final_combine,select=-c(Row.names))
    }
    top.relab.stats_table_combine<-merge(top.relab.stats_table_new,contam_result_final_combine, by="row.names")
    row.names(top.relab.stats_table_combine)<-top.relab.stats_table_combine$Row.names #previous merge step generate a new column call Row.names
    top.relab.stats_table_combine<- subset(top.relab.stats_table_combine,select=-c(Row.names))
    top.relab.stats_table_combine<<-top.relab.stats_table_combine  
    ####export csv
    csv_output_name<-paste0(output_excel,paste0(paste0("contam",paste0(paste0('_negC_',negative_sample_type),paste0("_compare_",name_compare_type))),"_all_level_SENSITIVITY.csv"))
    write.csv(top.relab.stats_table_combine,csv_output_name)
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    ###########################################generate plot###########################################
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    top.relab.stats_table_combine$taxa_full_name<-row.names(top.relab.stats_table_combine)
    total_column<-length(test_threshold_all_level)*4
    
    ################################
    ######EXPECTED TO BE CONTAMINANT
    ################################
    #this step keep the taxa which are listed in expected_contam_taxa. And then it check to see if the sample is labeled as TRUE for each method (identify as contaminant)
    #It then goes through each taxa and calcuates a percentage of correctly identified contaminant
    top.relab.stats_table_combine$expected<-NULL
    for (taxa in expected_contam_taxa) {
      top.relab.stats_table_combine[,taxa]<-grepl(taxa, top.relab.stats_table_combine$taxa_full_name , fixed = TRUE)
      top.relab.stats_table_combine$expected[top.relab.stats_table_combine[,taxa]==TRUE]<-TRUE
    }
    top.relab.stats_table_combine$expected_NO<-NULL
    for (taxa2 in expected_NOT_contam_taxa) {
      top.relab.stats_table_combine[,taxa2]<-grepl(taxa2, top.relab.stats_table_combine$taxa_full_name , fixed = TRUE)
      top.relab.stats_table_combine$expected_NO[top.relab.stats_table_combine[,taxa2]==TRUE]<-TRUE
    }
    top.relab.stats_table_combine_t<-top.relab.stats_table_combine
    top.relab.stats_table_combine_expected<-top.relab.stats_table_combine_t[(!is.na(top.relab.stats_table_combine_t$expected)), ]
    top.relab.stats_table_combine_expected_org<-top.relab.stats_table_combine_expected
    for (level in test_threshold_all_level){ 
      for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
        top.relab.stats_table_combine_expected[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_combine_expected[,paste0(method,level)]))
      }
    }
    top.relab.stats_table_combine_expected <- top.relab.stats_table_combine_expected[, -c(1:6)]# remove the relab abundance and rank columns
    top.relab.stats_table_combine_expected <- top.relab.stats_table_combine_expected[, c(1:total_column)] #keeping only the testings results 
    top.relab.stats_table_combine_total_expected<-top.relab.stats_table_combine_expected %>% dplyr::summarise_all(funs(mean))
    top.relab.stats_table_combine_total_expected$taxa_name<-"Contaminant"
    #for individual taxa of interest that is in the expected contaminant group
    for (taxa in expected_contam_taxa){
      top.relab.stats_table_combine_taxa<-top.relab.stats_table_combine_expected_org[(top.relab.stats_table_combine_expected_org[,taxa]), ]
      for (level in test_threshold_all_level){ 
        for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
          top.relab.stats_table_combine_taxa[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_combine_taxa[,paste0(method,level)]))
        }
      }
      top.relab.stats_table_combine_taxa <- top.relab.stats_table_combine_taxa[, -c(1:6)]# remove the relab abundance and rank columns
      top.relab.stats_table_combine_taxa <- top.relab.stats_table_combine_taxa[, c(1:total_column)] #keeping only the testings results 
      top.relab.stats_table_combine_taxa <- top.relab.stats_table_combine_taxa %>% dplyr::summarise_all(funs(mean))
      top.relab.stats_table_combine_taxa$taxa_name<-taxa
      assign(paste0("top.relab.stats_table_combine_taxa",taxa),top.relab.stats_table_combine_taxa)
    }
    ####################################
    ######EXPECTED NOT TO BE CONTAMINANT
    ####################################
    #this step keep the taxa which are listed in expected_NOT_contam_taxa And then it check to see if the sample is labeled as FALSE for each method (identify as NOT being contaminant)
    #It then goes through each taxa and calculates a percentage of correctly identified NON-contaminant
    top.relab.stats_table_combine_expected_NO<-top.relab.stats_table_combine_t[(!is.na(top.relab.stats_table_combine_t$expected_NO) ), ]
    top.relab.stats_table_combine_expected_NO_org<-top.relab.stats_table_combine_expected_NO
    for (level in test_threshold_all_level){ 
      for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
        top.relab.stats_table_combine_expected_NO[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_combine_expected_NO[,paste0(method,level)])) #see how many are not label as contaminant for the ones which SHOULD NOT be contaminant
      }
    }
    top.relab.stats_table_combine_expected_NO <- top.relab.stats_table_combine_expected_NO[, -c(1:6)]# remove the relab abundance and rank columns
    top.relab.stats_table_combine_expected_NO <- top.relab.stats_table_combine_expected_NO[, c(1:total_column)] #keeping only the testings results 
    top.relab.stats_table_combine_total_expected_NO<-top.relab.stats_table_combine_expected_NO %>% dplyr::summarise_all(funs(mean))
    for (level in test_threshold_all_level){ 
      for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
        top.relab.stats_table_combine_total_expected_NO[,paste0(method,level)]<-1-top.relab.stats_table_combine_total_expected_NO[,paste0(method,level)] #see how many are not label as contaminant for the ones which SHOULD NOT be contaminant
      }
    }
    top.relab.stats_table_combine_total_expected_NO$taxa_name<-"Non-contaminant"
    #for individual taxa of interest that is in the NOT expected contaminant group
    for (taxa in expected_NOT_contam_taxa){
      top.relab.stats_table_NO_taxa<-top.relab.stats_table_combine_expected_NO_org[(top.relab.stats_table_combine_expected_NO_org[,taxa]), ]
      dsfsdf<<-top.relab.stats_table_NO_taxa
      for (level in test_threshold_all_level){ 
        for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
          top.relab.stats_table_NO_taxa[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_NO_taxa[,paste0(method,level)]))
        }
      }
      top.relab.stats_table_NO_taxa <- top.relab.stats_table_NO_taxa[, -c(1:6)]# remove the relab abundance and rank columns
      top.relab.stats_table_NO_taxa <- top.relab.stats_table_NO_taxa[, c(1:total_column)] #keeping only the testings results 
      top.relab.stats_table_NO_taxa <- top.relab.stats_table_NO_taxa %>% dplyr::summarise_all(funs(mean))
      for (level in test_threshold_all_level){ 
        for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
          top.relab.stats_table_NO_taxa[,paste0(method,level)]<-1-top.relab.stats_table_NO_taxa[,paste0(method,level)]
        }
      }
      top.relab.stats_table_NO_taxa$taxa_name<-taxa
      assign(paste0("top.relab.stats_table_NO_taxa",taxa),top.relab.stats_table_NO_taxa)
    }
    #####combine all the files
    TOTAL_dataframe<-top.relab.stats_table_combine_total_expected
    for (taxa in expected_contam_taxa){
      TOTAL_dataframe<-rbind(TOTAL_dataframe, get(paste0("top.relab.stats_table_combine_taxa",taxa)))
    }
    TOTAL_dataframe<-rbind(TOTAL_dataframe, top.relab.stats_table_combine_total_expected_NO)
    for (taxa in expected_NOT_contam_taxa){
      TOTAL_dataframe<-rbind(TOTAL_dataframe, get(paste0("top.relab.stats_table_NO_taxa",taxa)))
    }
    #change wide to long for each method and stack them to graph 
    TOTAL_prev<-TOTAL_dataframe %>% dplyr::select(starts_with(c("prev","taxa")))
    TOTAL_prev_t<-reshape2::melt(TOTAL_prev, id.vars=c("taxa_name"))
    TOTAL_prev_t$variable<-as.numeric(str_sub(TOTAL_prev_t$variable, -3,-1))
    TOTAL_prev_t$method<-"prev"
    
    TOTAL_freq<-TOTAL_dataframe %>% dplyr::select(starts_with(c("freq","taxa")))
    TOTAL_freq_t<-reshape2::melt(TOTAL_freq, id.vars=c("taxa_name"))
    TOTAL_freq_t$variable<-as.numeric(str_sub(TOTAL_freq_t$variable, -3,-1))
    TOTAL_freq_t$method<-"freq"
    
    TOTAL_combined<-TOTAL_dataframe %>% dplyr::select(starts_with(c("combined","taxa")))
    TOTAL_combined_t<-reshape2::melt(TOTAL_combined, id.vars=c("taxa_name"))
    TOTAL_combined_t$variable<-as.numeric(str_sub(TOTAL_combined_t$variable, -3,-1))
    TOTAL_combined_t$method<-"combined"
    
    TOTAL_wilcox<-TOTAL_dataframe %>% dplyr::select(starts_with(c("wilcox","taxa")))
    TOTAL_wilcox_t<-reshape2::melt(TOTAL_wilcox, id.vars=c("taxa_name"))
    TOTAL_wilcox_t$variable<-as.numeric(str_sub(TOTAL_wilcox_t$variable, -3,-1))
    TOTAL_wilcox_t$method<-"wilcox"
    
    #########################################################################################################################################################
    #########################################################################################################################################################
    #####determine percent of contaminant identify (this step go back to the initial dataset with ALL taxa and then see how many taxa is identified as contaminant to get a percentage)
    top.relab.stats_table_combine_TOTAL <- top.relab.stats_table_combine[, -c(1:6)]# remove the relab abundance and rank columns
    top.relab.stats_table_combine_TOTAL <- top.relab.stats_table_combine_TOTAL[, c(1:total_column)] #keeping only the testings results 
    top.relab.stats_table_combine_TOTAL<-top.relab.stats_table_combine_TOTAL %>% dplyr::summarise_all(funs(mean))
    top.relab.stats_table_combine_TOTAL$taxa_name2<-"TOTAL"
    
    TOTAL2_prev<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("prev","taxa")))
    TOTAL2_prev_t<-reshape2::melt(TOTAL2_prev, id.vars=c("taxa_name2"))
    TOTAL2_prev_t$variable<-as.numeric(str_sub(TOTAL2_prev_t$variable, -3,-1))
    TOTAL2_prev_t$method<-"prev"
    
    TOTAL2_freq<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("freq","taxa")))
    TOTAL2_freq_t<-reshape2::melt(TOTAL2_freq, id.vars=c("taxa_name2"))
    TOTAL2_freq_t$variable<-as.numeric(str_sub(TOTAL2_freq_t$variable, -3,-1))
    TOTAL2_freq_t$method<-"freq"
    
    TOTAL2_combined<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("combined","taxa")))
    TOTAL2_combined_t<-reshape2::melt(TOTAL2_combined, id.vars=c("taxa_name2"))
    TOTAL2_combined_t$variable<-as.numeric(str_sub(TOTAL2_combined_t$variable, -3,-1))
    TOTAL2_combined_t$method<-"combined"
    
    TOTAL2_wilcox<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("wilcox","taxa")))
    TOTAL2_wilcox_t<-reshape2::melt(TOTAL2_wilcox, id.vars=c("taxa_name2"))
    TOTAL2_wilcox_t$variable<-as.numeric(str_sub(TOTAL2_wilcox_t$variable, -3,-1))
    TOTAL2_wilcox_t$method<-"wilcox"
    #########################################################################################################################################################                  
    #########################################################################################################################################################                 
    
    for (method in c("prev", "freq","combined","wilcox")){
      g<-ggplot()+
        geom_line(data=get(paste0(paste0("TOTAL_",method),"_t")) %>% filter(taxa_name!="Contaminant" & taxa_name!="Non-contaminant" &  taxa_name!="TOTAL"), 
                  aes(x=variable, y=value, linetype =taxa_name), linewidth=0.3)+ 
        geom_line(data=get(paste0(paste0("TOTAL_",method),"_t")) %>% filter(taxa_name=="Contaminant" | taxa_name=="Non-contaminant" ), 
                  aes(x=variable, y=value, color=taxa_name), linewidth=2)+
        geom_line(data=get(paste0(paste0("TOTAL2_",method),"_t")), 
                  aes(x=variable, y=value), colour= "darkorchid1", linewidth=1, alpha=0.5)+
        labs(title=method,x="test threshold",y="percent identified")+
        scale_linetype_discrete(name = "GENUS")+ 
        scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
        scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
        scale_color_discrete(name = "Correctly Identified As Expected")+ 
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
              panel.grid.minor = element_blank(),strip.background=element_blank(),
              plot.title=element_text(face="bold",hjust=0.5, size = 28), 
              plot.subtitle = element_text(hjust=0.5),
              axis.title=element_text(face="bold", size = 18),
              axis.text.x=element_text(face="bold",size = 16),
              axis.text.y=element_text(face="bold",size = 16),
              axis.ticks=element_blank(),
              plot.margin=unit(c(1,1,1,1),"line"))
      assign(paste0("graph_",method),g)
    }
    
    remove_y <- theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank())
    
    p<-list(graph_prev,
            graph_freq+remove_y,
            graph_combined+remove_y,
            graph_wilcox+remove_y)
    #put the graphs next to each other. add a title and collec the legend- place it at bottom and horizontal
    combine_p<-wrap_plots(p, nrow = 1) + plot_layout(guides = "collect")  &  theme(legend.position='bottom',legend.direction = "horizontal") &
      plot_annotation(title = paste0(paste0(negative_sample_type," compare with "),name_compare_type),
                      theme= theme(title = element_text(color = "red", face = "bold", size = 24))) 
    #output graphs in pdf format
    graph_output_title<-paste0(output_excel,paste0(paste0(paste0(negative_sample_type,"_compare_with_"),name_compare_type),".pdf"))
    pdf(graph_output_title, width=22, height=7)
    show(combine_p)
    dev.off()
    
  }
}


####decontaminant_subplot_KW: identify decontaminant based on a particular method and create boxplot for each sample type  
decontaminant_subplot_KW <- function (input_phyloseq, 
                                      SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                                      sample_type_var_name,                                       
                                      sample_types=list(), 
                                      sample_type_color=list(), 
                                      sample_type_color_2nd=list(), 
                                      negative_sample_type=list(), ###the sample type that you want to be as negative control
                                      compare_type=list(), 
                                      method_type=c("MannWhit","preval","freq","combin"),
                                      stat_option=c("mean", "median"), ### the statistics to determine rank for the plot
                                      graph_option=c("boxplot","mean_SE", "mean_SD"), 
                                      test_threshold, ###prevalence level for the contaminant test
                                      output_suffix)  {  #
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  method_type <- match.arg(method_type)
  print(method_type)
  print(stat_option)  
  
  if (missing(output_suffix)) { output_suffix<-"OP"} #if output_suffix is missing, then will just give a suffix _OP at the end of the file
  if (missing(graph_option)) { graph_option<-"boxplot"} #if graph_option is missing, then will just have boxplot as default
  
  ##making sure sample_types and sample_type_color has same number of elements
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #if sample_type_color_2nd is not missing then make sure it has same number of element as sample_type_color
  if (!missing(sample_type_color_2nd)) {
    stopifnot("sample_types need to have same number of elements as sample_type_color_2nd"= length(sample_types)==length(sample_type_color_2nd))
  } #if graph_option is missing, then will just have boxplot as default
  if (missing(sample_type_color_2nd)) {
    sample_type_color_2nd<-rep("black",length(sample_types)) #make the secondary color to be black by default, unless specify
  } 
  
  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
    assign(paste0("sample_type_color_2nd",i), sample_type_color_2nd[[i]])
    assign(paste0("sampletype",i), sample_types[[i]])
  }  
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <<- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))
  
  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #this is to make taxa name with rank order in them
  for (i in 1:length(sample_types)) { 
    j <- i+length(sample_types)
    k <- i+2*length(sample_types)  
    column2<-top.relab.stats[,j]
    top.relab.stats[k] <- paste0(rownames(top.relab.stats)," ","(",column2,")")
    colnames(top.relab.stats)[k] <- c(paste0(paste0("rank.order.",sample_types[i]),".name"))
  }  
  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  
  ###keeping only the statistics for each sample type 
  top.relab.stats_table <- top.relab.stats[c(1:length(sample_types))]
  relab.subset_final<<-relab.subset
  
  for (h in 1:length(compare_type)) {
    name_compare_type<-get(paste0("compare_type",h))
    print(paste0(negative_sample_type," compare:"))
    print(name_compare_type)
    AND<-grepl("_and_",name_compare_type, fixed=T)
    OR<-grepl("_or_",name_compare_type, fixed=T)
    COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
    contam<-list()
    
    if ((AND==TRUE) | (OR==TRUE)) {
      if (AND==TRUE) {
        compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
      }else if (OR==TRUE) {
        compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
      }
      #loop through the two sample types 
      for (i in 1:length(compare_sublist)) { 
        assign(paste0("compare_sublist",i), compare_sublist[i])
      } 
      
      for (b in 1:length(compare_sublist)) {
        # only keeping the compare list  
        selected_sample_type<-get(paste0("compare_sublist",b))
        p<-match(selected_sample_type,sample_types) 
        comparison.list_select<- unlist(comparison.list[p]) 
        comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
        
        # Remove unused from conc, in this case libsizes 
        libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
        counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
        relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
        
        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
        #remove the new column "Row.names" created when merging by 0
        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
        
        negative_comb<- append(comparison.list_select,comparison.list_negative)
        libsizes_comb <- append(libsizes_select,libsizes_negative )
        
        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="prevalence",
                                     neg=negative_comb,
                                     threshold=test_threshold)
        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="frequency",
                                     neg=negative_comb,
                                     conc=libsizes_comb,
                                     threshold=test_threshold)
        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="combined",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                            neg=negative_comb,
                                            threshold=test_threshold)
        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
        
        contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        assign(paste0("contam_result",b),contam_result)
      }
      #combine the results for the two sample types
      
      contam_result1a<<-contam_result1
      contam_result2a<<-contam_result2
      
      contam_result_final<-merge(contam_result1,contam_result2, by=0)
      rownames(contam_result_final) <- contam_result_final$Row.names
      contam_result_final <- subset(contam_result_final, select=-c(Row.names))
      if (AND==TRUE) {
        contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }else if (OR==TRUE) {
        contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
        contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
    }else{
      if (COMBINE==TRUE) { ###if compare_type is a combination set 
        ###this works if A_combine_B... it basically combine subset A and subset B
        combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
        f<-match(combine_list,sample_types) 
        comparison.list_select<- unlist(comparison.list[f]) 
        comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
        combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
        #Setting up the dataset so that it include both sample types that are combine =
        libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
        counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
        relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
      } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
        
        if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
          f<-match(name_compare_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[f]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          #keeping only the sample type of interest
          libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
        } else { stop(paste0(paste0(name_compare_type," must be in sample type list:"),sample_types)) }
        
      }
      #combining with the negative control group
      counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
      relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
      #remove the new column "Row.names" created when merging by 0
      rownames(counts.subset_comb) <- counts.subset_comb$Row.names
      rownames(relab.subset_comb) <- relab.subset_comb$Row.names
      counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
      relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
      
      negative_comb<- append(comparison.list_select,comparison.list_negative)
      libsizes_comb <- append(libsizes_select,libsizes_negative )
      
      contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                   method="prevalence",
                                   neg=negative_comb,
                                   threshold=test_threshold)
      contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                   method="frequency",
                                   neg=negative_comb,
                                   conc=libsizes_comb,
                                   threshold=test_threshold)
      contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="combined",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
      contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                          neg=negative_comb,
                                          threshold=test_threshold)
      contam_dataframe<- contam %>% data.frame #compress the lists from the different results
      contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
      
      contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
    }
    lowerbound<-length(sample_types)+1
    upperbound<-2*length(sample_types)
    rankorder <- subset(top.relab.stats[c(lowerbound:upperbound)]) ##keeping the rank order only
    contam_result_final_t<-contam_result_final
    if (method_type=="MannWhit"){
      contam_result_final_t$contaminant <- contam_result_final_t$wilcox.contaminant  
    }
    else if (method_type=="preval") {
      contam_result_final_t$contaminant <- contam_result_final_t$prev.contaminant
    }
    else if (method_type=="freq") {
      contam_result_final_t$contaminant <- contam_result_final_t$freq.contaminant
    }
    else if (method_type=="combin") {
      contam_result_final_t$contaminant <- contam_result_final_t$combined.contaminant
    }
    contam_result_final_t<-subset(contam_result_final_t,select=c(contaminant))
    csv_output_contam <- merge(contam_result_final_t,rankorder,by="row.names")
    
    ### Come back from the DECONTAM section to add this in ###    
    # Adding contam data (saved to the end to add)
    ###figure out which taxa to highlight in red based on which test you choice 
    if (method_type == "MannWhit") {
      top.relab.stats$contaminant <- contam_result_final$wilcox.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("MannWhit")
    }
    else if (method_type == "preval") {
      top.relab.stats$contaminant <- contam_result_final$prev.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("preval")
    }
    else if (method_type == "freq") {
      top.relab.stats$contaminant <- contam_result_final$freq.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("freq")
    }
    else if (method_type == "combin") {
      top.relab.stats$contaminant <- contam_result_final$combined.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("combin")
    }
    ###########################################################################################################
    ###########################################################################################################
    top.relab.stats_final<<-top.relab.stats
    DecontKW_contaminant_list<-subset(top.relab.stats,select=c(contaminant))
    assign(paste0(paste0("Contam_list_NC_",negative_sample_type),paste0("_compare_",paste0(name_compare_type,paste0("_",output_suffix)))), DecontKW_contaminant_list,.GlobalEnv)
    print("###########################")
    print("Available contaminant list:")
    print(paste0(paste0("Contam_list_NC_",negative_sample_type),paste0("_compare_",paste0(name_compare_type,paste0("_",output_suffix)))))
    print("###########################")
    ###########################################################################################################
    ###########################################################################################################
    for (rank_sample_type in sample_types) {
      l<-match(rank_sample_type,sample_types) ###look at the element in rank_sample_type to see what their index at sample_types
      r<-l+length(sample_types)
      # Re-rank based upon the rank group decreasing from highest to lowest
      column3<-top.relab.stats[,r]
      top.relab.stats <- top.relab.stats[order(column3, decreasing = F), , drop = FALSE ]
      top.2 <- top.relab.stats
      top.3 <- head(top.2, 100)
      # Filtering ONLY the top 100 taxa 
      relab.subset.top <- relab.subset %>% 
        dplyr::filter(rownames(relab.subset) %in% rownames(top.3))
      # Re ordered and saved into a new group (previously filtered)
      # relab.subset.top.match now have only 100 taxa ordered by top background 
      relab.subset.top.match <- relab.subset.top[match(rownames(top.3),rownames(relab.subset.top)),]    
      
      top.3$color <- ifelse(top.3$contaminant, "red", "black")
      
      for (a in sample_types) {
        assign(paste0("rank_order_",a),subset(top.3, select=c(paste0("rank.order.",a), "color", paste0(paste0("rank.order.",a),".name") )))
      }      
      
      # Transposition of the figure, maybe at this point replace the column names 
      # The previous edit should have removed all character columns 
      relab.subset.top.match.transpose <- as.data.frame(t(relab.subset.top.match))
      # name a new column with categories 
      relab.subset.top.match.transpose$category <- rownames(relab.subset.top.match.transpose)
      
      for (a in sample_types) {
        relab.subset.top.match.transpose$category[grepl(a,relab.subset.top.match.transpose$category)] <- a
      }  
      
      # Reverse column order 
      relab.subset.top.match.transpose.reverse <- relab.subset.top.match.transpose[,order(ncol(relab.subset.top.match.transpose):1)]
      relab.subset.top.match.transpose.reverse_final<<-relab.subset.top.match.transpose.reverse
      
      for (a in sample_types) {
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(relab.subset.top.match.transpose.reverse, category == c(a) ))
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(get(paste0("relab.subset.top.match.transpose.reverse.",a)), select=-c(category)))
      }  
      
      #get the first plot- rank_sample_type is the sample type that is being ranked.... we want the ranked plot on the left 
      temp <- get(paste0("relab.subset.top.match.transpose.reverse.",rank_sample_type)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
      temp <- merge(temp, top.3, by.x="id", by.y="row.names")
      p <- temp %>%   mutate(id = fct_reorder(id, -get(paste0("rank.order.",rank_sample_type)))) %>%
        ggplot(., aes(x=id, y=value)) + 
        geom_boxplot(color=sample_type_color[[l]], alpha=0.2) + ylab("Relative Abundance") +
        scale_x_discrete(
          labels = function(x) {
            x_org<-x
            is_long <- nchar(x) > 35
            x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                 paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                        paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                               paste0("..",str_sub(x[is_long],-4,-1)))))
            is_short <- nchar(x) < 17
            x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                         paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                       paste0("..",str_sub(x_org[is_short],-4,-1))))))
            
            x
          }
        ) +
        coord_flip() +
        theme_bw() + 
        theme(axis.text.y=element_text(color=rev(top.3$color)),axis.title.y=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1))+
        ggtitle(rank_sample_type) 
      
      myplots <- list()  # new empty list
      myplots[[1]]<- p
      
      #remove the sample type that is being ranked. so new_sample_type_name and new_sample_type_color contain the sample types which were not ranked                
      new_sample_type_name<-sample_types[-l] #remove the lth element (lth element is the sample type that is being ranked by )
      new_sample_type_color<-sample_type_color[-l] #remove the lth element (lth element is the sample type that is being ranked by )
      for (q in 1:length(new_sample_type_name)) {
        assign(paste0("new_sample_type_color",q), new_sample_type_color[[q]])
        assign(paste0("new_sample_type_name",q), new_sample_type_name[[q]])
      }  
      #remove the axis label for the sample types which were not ranked (since the ranked sample type would be on the left most)  
      remove_y<- theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank())
      
      ###looping thro each sample type 
      number<-1
      number2<-2
      for (a in new_sample_type_name) {
        temp <- get(paste0("relab.subset.top.match.transpose.reverse.",a)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
        temp <- merge(temp, top.3, by.x="id", by.y="row.names")
        p <- temp %>%   mutate(id = fct_reorder(id, -get(paste0("rank.order.",rank_sample_type)))) %>%
          ggplot(., aes(x=id, y=value)) + 
          geom_boxplot(color=new_sample_type_color[[number]], alpha=0.2) + ylab("Relative Abundance") +
          scale_x_discrete(
            labels = function(x) {
              x_org<-x
              is_long <- nchar(x) > 35
              x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                   paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                          paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                                 paste0("..",str_sub(x[is_long],-4,-1)))))
              is_short <- nchar(x) < 17
              x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                    paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                           paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                         paste0("..",str_sub(x_org[is_short],-4,-1))))))
              x
            }
          ) +
          coord_flip() +
          theme_bw() + 
          theme(axis.text.y=element_text(color=rev(top.3$color)),axis.title.y=element_blank()) +
          theme(axis.text.x=element_text(angle=90,hjust=1))+
          ggtitle(a) 
        myplots[[number2]] <- p + remove_y
        number<-number+1
        number2<-number2+1
      }
      
      pdf_output=paste0(paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("_rank_",rank_sample_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".pdf"))
      pdf(file=pdf_output, width=14+(length(sample_types)-3)*4, height=(14+(length(sample_types)-3)*4)*10/8.5)
      show(wrap_plots(myplots, nrow=1))
      dev.off()        
      #merging in the statistics
      top.relab.stats_table$Row.names<-row.names(top.relab.stats_table)
      csv_output_contam_final<-merge(csv_output_contam, top.relab.stats_table, by="Row.names")
      
      ####export csv
      csv_output=paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".csv"))
      write.csv(csv_output_contam_final,csv_output, row.names = FALSE)
      xlsx_output=paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".xlsx"))
      write.xlsx(csv_output_contam_final, xlsx_output, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
    }    
    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################
    ##Generating boxplot rank by the corresponding sample type 
    count_table<<-as.data.frame(relab.subset_final)
    for (g in 1:length(sample_types)) {
      assign(paste0("relab.subset_final_sampletype",g), count_table[,grepl(get(paste0("sampletype",g)), colnames(count_table))])
    }
    top.relab.stats_final$color <<- ifelse(top.relab.stats_final$contaminant, "red", "black")
    
    end_bound<-3*length(sample_types)+2
    for (g in 1:length(sample_types)) {
      h<-g+length(sample_types)
      assign(paste0("rank_order_sample_type",g), as.data.frame(top.relab.stats_final[,c(h)]))
      temp_df<-get(paste0("rank_order_sample_type",g))
      rownames(temp_df)<-rownames(top.relab.stats_final)
      colnames(temp_df)<-"rank"
      temp_df3<-as.data.frame(top.relab.stats_final[,c(h,end_bound)])
      temp_df3<-temp_df3[order(temp_df3$rank, decreasing = T), ,drop = FALSE ]
      temp_df3<-tail(temp_df3, 100)
      assign(paste0("rank_order_sample_type",g), temp_df)
      assign(paste0("color_rank_order_sample_type",g), temp_df3) 
    }
    for (g in 1:length(sample_types)) {
      assign(paste0("relab.subset_final_sampletype_merge",g),merge(get(paste0("relab.subset_final_sampletype",g)),get(paste0("rank_order_sample_type",g)),by="row.names"))
      temp_df2<-get(paste0("relab.subset_final_sampletype_merge",g))
      rownames(temp_df2)=temp_df2$Row.name
      temp_df2=temp_df2[2:length(temp_df2)]
      temp_df2 <- temp_df2[order(temp_df2$rank, decreasing = F), , drop = FALSE ]
      temp_df2 <- head(temp_df2, 100)#keep top 100
      temp_df2 = subset(temp_df2, select = -c(rank))
      temp_df2 <- as.data.frame(t(temp_df2))
      temp_df2_b <- temp_df2[,order(ncol(temp_df2):1)]
      assign(paste0("relab.subset_final_sampletype_merge.reverse",g),temp_df2_b)
    }
    ###looping thro each sample type 
    number3<-1
    myplots_b<-list()
    for (g in sample_types) {
      pp <- get(paste0("relab.subset_final_sampletype_merge.reverse",number3)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE) %>%
        ggplot(., aes(x=id, y=value)) + 
        ylab("Relative Abundance") +
        scale_x_discrete(
          labels = function(x) {
            x_org<-x
            is_long <- nchar(x) > 35
            x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                 paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                        paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                               paste0("..",str_sub(x[is_long],-4,-1)))))
            is_short <- nchar(x) < 17
            x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                         paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                       paste0("..",str_sub(x_org[is_short],-4,-1))))))
            
            
            x
          }
        ) +
        coord_flip() +
        theme_bw() + 
        theme(axis.text.y=element_text(color=get(paste0("color_rank_order_sample_type",number3))$color),axis.title.y=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1))+
        ggtitle(g) 
      
      if (graph_option=="boxplot") {
        pp<-pp+geom_boxplot(color=sample_type_color[[number3]], alpha=0.2) 
      } else if(graph_option=="mean_SD") {
        pp<-pp+geom_jitter(color=sample_type_color[[number3]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x),
                       fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                       color=sample_type_color_2nd[[number3]], linewidth =1.0, size=0.4)
      } else if(graph_option=="mean_SE") {
        pp<-pp+geom_jitter(color=sample_type_color[[number3]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x)/sqrt(length(x)),
                       fun.min = function (x) ifelse( mean(x)-sd(x)/sqrt(length(x)) < 0,0, mean(x)-sd(x)/sqrt(length(x))),
                       color=sample_type_color_2nd[[number3]], linewidth =1.0, size=0.4)
      }
      myplots_b[[number3]] <- pp
      number3<-number3+1
    }
    all_plot.all <- ggarrange(plotlist=myplots_b, ncol=length(myplots), nrow=1, align="h")
    pdf_output=paste0(paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),"top"),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".pdf"))
    pdf(file=pdf_output, width=20, height=10)
    show(all_plot.all)
    dev.off()    
    #####################################################################################################
    #####################################################################################################
  }
}

##################################################################################
##################################################################################

####### generate background identification plots using the following functions (Supp Figure 3): 
decontaminant_subplot_KW (input_phyloseq= physeq.3, 
                          SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                          sample_type_var_name="Sample_Type_Involved",                                       
                          sample_types= c("BKG", "Lung.Tissue.In", "Lung.Tissue.UnIn"),
                          sample_type_color=c("grey", "red", "blue"), 
                          negative_sample_type="BKG", ###the sample type that you want to be as negative control
                          compare_type=c("Lung.Tissue.In_or_Lung.Tissue.UnIn"), 
                          method_type="prev",
                          stat_option="mean", ### the statistics to determine rank for the boxplot
                          test_threshold = 0.5) ###prevalence level for the contaminant test"  




##################################################################################
##################################################################################

#####cox PH model forest plot (Supp Figure 7)########


###from Supp table 7, filter the taxa that are significant (FDR <0.2). 
# copy these and create a new excel / csv or txt file and read it 


####Tumor results##### (Supp Figure 7A)

#only significant results with padj < 0.2 were choosen 

res <- readxl::read_excel("COX.PH.Tumor.Sig.xlsx")

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

#set color vector accroding to contam
Tumor.forest.df <- Tumor.forest.df %>% 
  mutate(contam= factor(contam)) %>% 
  mutate(color=ifelse(Tumor.forest.df$contam =="TRUE", "red", "black")) %>% 
  mutate(upper_CI= as.numeric(upper_CI)) %>% 
  mutate(lower_CI= as.numeric(lower_CI)) %>% 
  mutate(HR= as.numeric(HR)) %>% 
  mutate(nameASV.trim.colored = Tumor.forest.df$nameASV.trim) %>% 
  arrange(., desc(HR)) %>% 
  dplyr::mutate(ord = order(as.numeric(rownames(.))))

#add color of potential contaminants 
Tumor.forest.df$nameASV.trim.colored <- paste0("<span style=\"color: ", Tumor.forest.df$color, "\">", Tumor.forest.df$nameASV.trim.colored, "</span>")

#plot 
ggplot(data=Tumor.forest.df, aes(x=HR, y=fct_reorder(nameASV.trim.colored, -ord), xmin=lower_CI, xmax=upper_CI))+
  geom_point(size=5)+
  geom_errorbarh(height=0.2)+
  geom_vline(xintercept = 1, color="black", linetype="dashed", alpha=0.5)+
  ylab("")+xlab("")+
  theme_classic()+
  theme(axis.text.y= element_markdown(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20))


######the above can be repeated without contaminant data



####Lung samples#### 
#only significant results with padj < 0.2 were choosen 

res <- readxl::read_excel("COX.PH.Lung.Sig.xlsx")

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
Lung.forest.df <- inner_join(res, taxtable, by="asv")
Lung.forest.df$nameASV.trim <- paste0(Lung.forest.df$names,paste0(".."), paste0(str_sub(Lung.forest.df$asv, - 3, - 1)))

#set color vector accroding to contam
Lung.forest.df <- Lung.forest.df %>% 
  mutate(contam= factor(contam)) %>% 
  mutate(color=ifelse(Lung.forest.df$contam =="TRUE", "red", "black")) %>% 
  mutate(upper_CI= as.numeric(upper_CI)) %>% 
  mutate(lower_CI= as.numeric(lower_CI)) %>% 
  mutate(HR= as.numeric(HR)) %>% 
  mutate(nameASV.trim.colored = Lung.forest.df$nameASV.trim) %>% 
  arrange(., desc(HR)) %>% 
  dplyr::mutate(ord = order(as.numeric(rownames(.))))

#add color of potential contaminants 
Lung.forest.df$nameASV.trim.colored <- paste0("<span style=\"color: ", Lung.forest.df$color, "\">", Lung.forest.df$nameASV.trim.colored, "</span>")

#plot 
ggplot(data=Lung.forest.df, aes(x=HR, y=fct_reorder(nameASV.trim.colored, -ord), xmin=lower_CI, xmax=upper_CI))+
  geom_point(size=5)+
  geom_errorbarh(height=0.2)+
  geom_vline(xintercept = 1, color="black", linetype="dashed", alpha=0.5)+
  ylab("")+xlab("")+
  theme_classic()+
  theme(axis.text.y= element_markdown(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20))
