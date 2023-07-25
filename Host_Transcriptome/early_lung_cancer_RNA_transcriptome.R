
################################################################################
################################################################################
##################################################################################        
################################################################################        

######################Host RNA transcriptome Analysis##########################        

######Pre processing 

#Load RNA Seq file 
mycounts <-read.delim2("NovaSeq_Raw_Counts_Gene_Symbols.txt", sep="\t", row.names=1)

#Keep Only MetaData of RNASeq Samples
#Keep only samples with RNA ID that matches the RNA count data 
RNA.data = mapping.table[mapping.table$RNA_Seq_done_Final==1,]
RNA.data <- RNA.data[RNA.data$RNA_ID!= "n.a",]

#Pick up the column that have RNA ID in metadata #We want to match ID names in each element 
RNA.data$RNA_ID

#Order Meta Data by SampleId
RNA.data <- RNA.data[order(RNA.data$RNA_ID),]

#Order Count Data by SampleID and replace the X before each name (these are lung samples IDs)
colnames(mycounts) <- gsub("X", "", colnames(mycounts))
mycounts <- mycounts [, order(colnames(mycounts))]

#Confirm Sample IDs match for Count and Meta Data
table(colnames(mycounts)==as.character(RNA.data$RNA_ID))
rownames(RNA.data) <- RNA.data$RNA_ID

#confirm columns names of data and rownamnes of metadata are similar 
table(colnames(mycounts)==rownames(RNA.data))


##################################################################################        

#########################Tumor vs Lung Analysis################################

#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts))

#Convert Data to Integers to Run DESEq so you can get normalized data to be used in beta diversity analysis 
d1[] <- lapply(d1, as.integer)

#Convert Model Variable into Factor
RNA.data$Sample_Type_Involved <- as.factor(RNA.data$Sample_Type_Involved)

#create DESeq object 
dds <- DESeqDataSetFromMatrix(countData = d1, colData = RNA.data, design = ~ Sample_Type_Involved)

#Normalization Step 
dds <- estimateSizeFactors(dds)

#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]

#Transform Data 
vsd <- varianceStabilizingTransformation(dds)

#Drop Levels
dds$Sample_Type_Involved   <- droplevels(dds$Sample_Type_Involved)
vsd$Sample_Type_Involved   <- droplevels(vsd$Sample_Type_Involved)



################################################################################

#plot PCoA of Tumor vs lung (Supp Figure 8A)

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

#plot 
ggplot(newResults, aes(PC1, PC2, color= Sample_Type_Involved)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  scale_color_manual(values=c("blue", "red")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Sample_Type_Involved), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample_Type_Involved))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Lung", "Tumor")), size=10) +
  scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))

#Stats 
adonis2(vegdist ~ vsd$Sample_Type_Involved)


#####repeat for Tumor, Recurrence vs No Recurrence (Figure 3A)

#using DESeq, create a new dds object and model it against tumor progression 
dds <- DESeqDataSetFromMatrix(countData = d1, colData = RNA.data, design = ~ Progression)
vsd <- varianceStabilizingTransformation(dds)

#subset tumor table 
dds.involved <- dds [, dds$Sample_Type_Involved %in% "Lung.Tissue.In"]
vsd.involved <- vsd[, vsd$Sample_Type_Involved  %in% "Lung.Tissue.In"]

##Plot PCoA

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

ggplot(newResults, aes(PC1, PC2, color= Progression)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("orange", "#E0115F")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence")), size=10) +
  scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")+
  theme(plot.title = element_text(hjust = 0.2))

#Stats 
adonis2(vegdist ~ vsd.involved$Progression)


#######Repeat for Lung samples, Recurrence vs no Recurrence (Figure 3C)

dds.Uninvolved <- dds [, dds$Sample_Type_Involved %in% "Lung.Tissue.UnIn"]
vsd.Uninvolved <- vsd[, vsd$Sample_Type_Involved  %in% "Lung.Tissue.UnIn"]

####Plot PCoA

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

ggplot(newResults, aes(PC1, PC2, color= Progression)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("Blue", "#92A1CF")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence")), size=10) +
  scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")+
  theme(plot.title = element_text(hjust = 0.216))

#Stats 
adonis2(vegdist ~ vsd.Uninvolved$Progression)





################################################################################
################################################################################
################################################################################
################# Differential Enrichment Analysis; edgeR #########################

#Tumor vs lung # (Supp Figure 8B)

#get genes table and add 1
x= as(mycounts, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA.data, "Sample_Type_Involved")

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

#get results of edgeR
res <- resNoFilt$table
#Reverse Directionality if you need to  
res$logFC <- res$logFC*(-1)

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

# save results 
write.csv(res, file = "edgeR.results.Tumor.vs.Lung.csv")

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
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < 0.2 , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < 0.2, as.character(res$Gene.symbol),'')),size= 4, force = 25,segment.colour="grey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
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




#################################################################################

## Tumor vs Recurrence vs No recurrence### (Figure 3B)

#first select tumor samples only from metadata and counts table 
mycounts.Tumor <- mycounts %>% 
  dplyr::select(., starts_with("NYU"))

#get genes table and add 1
x= as(mycounts.Tumor, "matrix")
x= x+1 

#get tumor only table 
RNA.data.Tumor <- subset_samples(RNA.data, Sample_Type_Involved =="Lung.Tissue.In")

#overall you have 146 Tumor samples with RNA data 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA.data.Tumor, "Progression_Lab_Inv")

#create edgeR list 
dgeFull <- DGEList(counts = x, group = group, remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table

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
write.csv(res, file = "edgeR.results.Tumor_Rec_vs_noRec_.csv")

###### volcano plot (Figure 3B)
# create color variable 
res$sig <- -log10(res$FDR)
# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "#E0115F"
cols[res$logFC < 0 & res$FDR < alpha ] <- "orange"


#plot with labels 
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),size= 4, force = 25,segment.colour="grey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
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




##################################################################################


### edgeR Lung vs Recurrence vs No recurrence (Figure 3D)

#first select Lung samples only from metadata and counts table 
mycounts.Lung <- mycounts %>% 
  dplyr::select(., -starts_with("NYU"))

#get genes table and add 1
x= as(mycounts.Lung, "matrix")
x= x+1 

#get Lung only table 
RNA.data.Lung <- subset_samples(RNA.data, Sample_Type_Involved =="Lung.Tissue.UnIn")

# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA.data.Lung, "Progression_Lab_Inv")

#create edgeR list 
dgeFull <- DGEList(counts = x, group = group, remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")


res <- resNoFilt$table
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
write.csv(res, file = "edgeR.results.Lung_Rec_vs_noRec_.csv")

###### volcano plot (Figure 3D)
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "blue"
cols[res$logFC < 0 & res$FDR < alpha ] <- "#92A1CF"


#plot with labels 
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),size= 4, force = 25,segment.colour="grey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
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

##################################################################################
##################################################################################
##################################################################################


########################Random Forest of Genes (Supp Figure 10) ################################


#######Tumor samples 

#load genes dataset 
mycounts
#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
d1 = data.frame(lapply(mycounts, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

#Convert Model Variable into Factor
RNA.data$Progression_Lab_Inv <- as.factor(RNA.data$Progression_Lab_Inv)

dds <- DESeqDataSetFromMatrix(countData = d1, colData = RNA.data, design = ~ Progression_Lab_Inv)

#Normalization Step 
dds <- estimateSizeFactors(dds)
#filter genes where there are less than 3 samples with normalized counts greater than or equal to 500.
idx <- rowSums( counts(dds, normalized=TRUE) >= 500 ) >= 3
#keep only these genes 
dds <- dds[idx,]

#get genes list 
genes <- mycounts %>% 
  dplyr::filter(., idx)

#get tumor samples only 
genes <- genes %>% 
  dplyr::select(starts_with("NYU"))

#transpose to use in random forest model as features 
genes.transposed <- as.data.frame(t(genes))

#get tumor only data 
TUMOR.RNA.data <- subset_samples(RNA.data, Sample_Type_Involved == "Lung.Tissue.In")
metadata.for.RF <- sample_data(TUMOR.RNA.data) %>% 
  data.frame()

#be sure to include the samples needed 
genes.transposed<- genes.transposed %>% 
  dplyr::filter(rownames(genes.transposed) %in% rownames(metadata.for.RF))

#define outcome as recurrence 
outcome <- as.factor(metadata.for.RF$Progression_Lab)

#build a data frame of features and outcome 
rf.data_complete <- data.frame(genes.transposed, outcome)

#data preprocessing for random forest model 
#remove na 
rf.data_complete <- rf.data_complete %>% 
  drop_na()

#be sure data is numerical 
lapply(rf.data_complete, as.integer)

#replace infinite data with NA
rf.data_complete[mapply(is.infinite, rf.data_complete)] <- NA
rf.data_complete <- rf.data_complete %>% 
  mutate_if(is.character, as.factor)

# Split the data into training and testing sets
set.seed(1234)
train_index <- sample(nrow(rf.data_complete), 0.8 * nrow(rf.data_complete))
train_data <- rf.data_complete[train_index, ]
test_data <- rf.data_complete[-train_index, ]

#tune the model to get best results possible 
mtry <- sqrt(ncol(rf.data_complete))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = rf.data_complete, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get variable importance 
# Extract the mean decrease impurity (MDI) values for each feature
mdi <- randomForest::importance(fit)
mdi <- mdi %>% 
  data.frame() %>% 
  dplyr::select(MeanDecreaseGini)

# Normalize the MDI values to get Gini index values 
Gini <- mdi %>% 
  dplyr::mutate(Norm.Gini= MeanDecreaseGini/sum(MeanDecreaseGini)) %>% 
  dplyr::mutate(Gini = Norm.Gini/max(Norm.Gini)) %>% 
  dplyr::arrange(desc(Gini)) %>% 
  dplyr::select(Gini)

#add enrich group (or outcome)
gene.means <- by(t(genes), outcome, colMeans)
gene.means <- do.call(cbind, gene.means)
idx_enrich <- apply(gene.means, 1, which.max)
group_enrich <- colnames(gene.means)[idx_enrich]

#put this into your Gini index dataframe 
Gini$enrich_group <- group_enrich

#leave out genes with 0 value of Gini that don't affct the model (need to change every time you run it)
Gini.100.per<- Gini %>% 
  dplyr::slice(1:4952)
#from here subset genes of increasing percentage so you will use later to fit the model 
#top 1%
round(0.01*4952)
Genes.1.per <- Gini %>% 
  dplyr::slice(1:50)
#top 5%
round(0.05*4952)
Genes.5.per <- Gini %>% 
  dplyr::slice(1:248)
#top 10%
round(0.1*4952)
Genes.10.per <- Gini %>% 
  dplyr::slice(1:495)
#top 20%
round(0.2*4952)
Genes.20.per <- Gini %>% 
  dplyr::slice(1:990)
#top 50%
round(0.5*4952)
Genes.50.per <- Gini %>% 
  dplyr::slice(1:2476)
#top 75%
round(0.75*4952)
Genes.75.per <- Gini %>% 
  dplyr::slice(1:3714)
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
Ginidf.to.save$one_per [50:11560] <- NA
Ginidf.to.save$five_per [248:11560] <- NA
Ginidf.to.save$ten_per [495:11560] <- NA
Ginidf.to.save$twent_per [990:11560] <- NA
Ginidf.to.save$fifty_per [2476:11560] <- NA
Ginidf.to.save$sevent_fiv_per [3714:11560] <- NA

write.csv(Ginidf.to.save, file = "RF_results_Tumor_samples.csv")

######get ROC plor and AUC
#get predictions 
predictions <- predict(fit, type = "prob")
pred1 <- data.frame(predictions)

# to get pred vs observations for every row 
pred <- fit$predicted
obs <- fit$y

#now build a data frame for ROC curve
pred.df <- pred1 
pred.df$pred <- pred
pred.df$obs <- obs
pred.df <- data.frame(pred.df)
colnames(pred.df)

#to get a glimpse on sens and spec of the model:
cm <- confusionMatrix(data = pred.df$pred, reference = pred.df$obs,  mode="everything",  positive = "Recurrence" )

#now calculate sensitivity and spec. 
roc <- roc(ifelse(pred.df$obs=="Recurrence", "Recurrence", "non-Recurrence"), as.numeric(pred.df$Recurrence))

#using ggplot and pROC package to plot ROC
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

#save it 
pdf(roc.plot, file = "roc.plot.Tumor_samples_100_per.pdf", height = 8, width = 8)
roc.plot
dev.off()



####repeat model for changing numbers of top genes according to Gini index (here n=1% is provided, and it can be changed to 5%, 10%, 20%, 50%, 75%)
n <- 1

# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(Genes.1.per))

train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(rf.data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

#run the model
fit <- randomForest(outcome~., data = rf.data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions 
predictions <- predict(fit, type = "prob")
pred1 <- data.frame(predictions)

# to get pred vs observations for every row 
pred <- fit$predicted
obs <- fit$y

#now build a data frame for ROC curve
pred.df <- pred1 
pred.df$pred <- pred
pred.df$obs <- obs
pred.df <- data.frame(pred.df)
colnames(pred.df)

#to get a glimpse on sens and spec of the model:
cm <- confusionMatrix(data = pred.df$pred, reference = pred.df$obs,  mode="everything",  positive = "Recurrence" )

#now calculate sensitivity and spec. 
roc <- roc(ifelse(pred.df$obs=="Recurrence", "Recurrence", "non-Recurrence"), as.numeric(pred.df$Recurrence))

#using ggplot and pROC package to plot ROC
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

#save it 
pdf(roc.plot, file = paste0("roc.plot.Tumor_samples_", paste0(n, paste0("_per.pdf"))), height = 8, width = 8)
roc.plot
dev.off()

#this should be repeated for every top number of genes used 






######Lung Samples 


#load genes dataset 
mycounts
#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
d1 = data.frame(lapply(mycounts, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

#Convert Model Variable into Factor
RNA.data$Progression_Lab_Inv <- as.factor(RNA.data$Progression_Lab_Inv)

dds <- DESeqDataSetFromMatrix(countData = d1, colData = RNA.data, design = ~ Progression_Lab_Inv)

#Normalization Step 
dds <- estimateSizeFactors(dds)
#filter genes where there are less than 3 samples with normalized counts greater than or equal to 500.
idx <- rowSums( counts(dds, normalized=TRUE) >= 500 ) >= 3
#keep only these genes 
dds <- dds[idx,]

#get genes list 
genes <- mycounts %>% 
  dplyr::filter(., idx)

#get Lung samples only 
genes <- genes %>% 
  dplyr::select(-starts_with("NYU"))

#transpose to use in random forest model as features 
genes.transposed <- as.data.frame(t(genes))

#get Lung only data 
Lung.RNA.data <- subset_samples(RNA.data, Sample_Type_Involved == "Lung.Tissue.UnIn")
metadata.for.RF <- sample_data(Lung.RNA.data) %>% 
  data.frame()

#be sure to include the samples needed 
genes.transposed<- genes.transposed %>% 
  dplyr::filter(rownames(genes.transposed) %in% rownames(metadata.for.RF))

#define outcome as recurrence 
outcome <- as.factor(metadata.for.RF$Progression_Lab)

#build a data frame of features and outcome 
rf.data_complete <- data.frame(genes.transposed, outcome)

#data preprocessing for random forest model 
#remove na 
rf.data_complete <- rf.data_complete %>% 
  drop_na()

#be sure data is numerical 
lapply(rf.data_complete, as.integer)

#replace infinite data with NA
rf.data_complete[mapply(is.infinite, rf.data_complete)] <- NA
rf.data_complete <- rf.data_complete %>% 
  mutate_if(is.character, as.factor)

# Split the data into training and testing sets
set.seed(1234)
train_index <- sample(nrow(rf.data_complete), 0.8 * nrow(rf.data_complete))
train_data <- rf.data_complete[train_index, ]
test_data <- rf.data_complete[-train_index, ]

#tune the model to get best results possible 
mtry <- sqrt(ncol(rf.data_complete))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

######### 
#run the model
set.seed(1234)

fit <- randomForest(outcome~., data = rf.data_complete, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get variable importance 
# Extract the mean decrease impurity (MDI) values for each feature
mdi <- randomForest::importance(fit)
mdi <- mdi %>% 
  data.frame() %>% 
  dplyr::select(MeanDecreaseGini)

# Normalize the MDI values to get Gini index values 
Gini <- mdi %>% 
  dplyr::mutate(Norm.Gini= MeanDecreaseGini/sum(MeanDecreaseGini)) %>% 
  dplyr::mutate(Gini = Norm.Gini/max(Norm.Gini)) %>% 
  dplyr::arrange(desc(Gini)) %>% 
  dplyr::select(Gini)

#add enrich group (or outcome)
gene.means <- by(t(genes), outcome, colMeans)
gene.means <- do.call(cbind, gene.means)
idx_enrich <- apply(gene.means, 1, which.max)
group_enrich <- colnames(gene.means)[idx_enrich]

#put this into your Gini index dataframe 
Gini$enrich_group <- group_enrich

#leave out genes with 0 value of Gini that don't affct the model (need to change every time you run it)
Gini.100.per<- Gini %>% 
  dplyr::slice(1:4952)
#from here subset genes of increasing percentage so you will use later to fit the model 
#top 1%
round(0.01*4952)
Genes.1.per <- Gini %>% 
  dplyr::slice(1:50)
#top 5%
round(0.05*4952)
Genes.5.per <- Gini %>% 
  dplyr::slice(1:248)
#top 10%
round(0.1*4952)
Genes.10.per <- Gini %>% 
  dplyr::slice(1:495)
#top 20%
round(0.2*4952)
Genes.20.per <- Gini %>% 
  dplyr::slice(1:990)
#top 50%
round(0.5*4952)
Genes.50.per <- Gini %>% 
  dplyr::slice(1:2476)
#top 75%
round(0.75*4952)
Genes.75.per <- Gini %>% 
  dplyr::slice(1:3714)
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
Ginidf.to.save$one_per [50:11560] <- NA
Ginidf.to.save$five_per [248:11560] <- NA
Ginidf.to.save$ten_per [495:11560] <- NA
Ginidf.to.save$twent_per [990:11560] <- NA
Ginidf.to.save$fifty_per [2476:11560] <- NA
Ginidf.to.save$sevent_fiv_per [3714:11560] <- NA

write.csv(Ginidf.to.save, file = "RF_results_Lung_samples.csv")

######get ROC plor and AUC
#get predictions 
predictions <- predict(fit, type = "prob")
pred1 <- data.frame(predictions)

# to get pred vs observations for every row 
pred <- fit$predicted
obs <- fit$y

#now build a data frame for ROC curve
pred.df <- pred1 
pred.df$pred <- pred
pred.df$obs <- obs
pred.df <- data.frame(pred.df)
colnames(pred.df)

#to get a glimpse on sens and spec of the model:
cm <- confusionMatrix(data = pred.df$pred, reference = pred.df$obs,  mode="everything",  positive = "Recurrence" )

#now calculate sensitivity and spec. 
roc <- roc(ifelse(pred.df$obs=="Recurrence", "Recurrence", "non-Recurrence"), as.numeric(pred.df$Recurrence))

#using ggplot and pROC package to plot ROC
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

#save it 
pdf(roc.plot, file = "roc.plot.Lung_samples_100_per.pdf", height = 8, width = 8)
roc.plot
dev.off()

####repeat model for changing numbers of top genes according to Gini index (here n=1% is provided, and it can be changed to 5%, 10%, 20%, 50%, 75%)
n <- 1

# Split the data into training and testing sets
rf.data <- rf.data_complete %>% 
  dplyr::select(., rownames(Genes.1.per))

train_index <- sample(nrow(rf.data), 0.8 * nrow(rf.data))
train_data <- rf.data[train_index, ]
test_data <- rf.data[-train_index, ]

#tune 
mtry <- sqrt(ncol(rf.data))
tunegrid <- expand.grid(.mtry=mtry)
metric <- "ROC"
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

#run the model
fit <- randomForest(outcome~., data = rf.data, trControl=Control, metric=metric, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get predictions 
predictions <- predict(fit, type = "prob")
pred1 <- data.frame(predictions)

# to get pred vs observations for every row 
pred <- fit$predicted
obs <- fit$y

#now build a data frame for ROC curve
pred.df <- pred1 
pred.df$pred <- pred
pred.df$obs <- obs
pred.df <- data.frame(pred.df)
colnames(pred.df)

#to get a glimpse on sens and spec of the model:
cm <- confusionMatrix(data = pred.df$pred, reference = pred.df$obs,  mode="everything",  positive = "Recurrence" )

#now calculate sensitivity and spec. 
roc <- roc(ifelse(pred.df$obs=="Recurrence", "Recurrence", "non-Recurrence"), as.numeric(pred.df$Recurrence))

#using ggplot and pROC package to plot ROC
roc.df <- data.frame(roc$sensitivities, roc$specificities, roc$thresholds)

roc.plot <- roc.df %>%
  arrange(roc.thresholds) %>%
  ggplot() +
  geom_path(aes(x=1 - roc.specificities, y=roc.sensitivities)) + # connect the points in the order in which they appear in the data to form a curve
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
  coord_equal()+
  xlab("1- Specificity")+ ylab("Sensitivity")+
  theme_bw()+
  ggtitle(paste0("AUC = ", paste0(round(roc$auc,2))))

#save it 
pdf(roc.plot, file = paste0("roc.plot.Lung_samples_", paste0(n, paste0("_per.pdf"))), height = 8, width = 8)
roc.plot
dev.off()

#this should be repeated for every top number of genes used 



################################################################################
#######Plot AUC according to top genes percentage used (Supplementary Figure 10)

# tumor and lung samples AUC values generated from random forest for every top n of genes 
AUC.tumor <- c(0.82, 0.77, 0.78, 0.74, 0.67, 0.63, 0.46)
AUC.lung <- c(0.79, 0.82, 0.83, 0.81, 0.77,0.74, 0.69 )
perc.genes <- c(1, 5, 10, 20, 50, 75, 100)
AUC.df <- data.frame(perc.genes, AUC.tumor, AUC.lung)
colors <- c( "AUC.tumor" = "red", "AUC.lung" = "blue")

#Supp Figure 10A
ggplot(AUC.df, aes(x=perc.genes))+
  geom_line(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_point(aes(y=AUC.tumor, color="AUC.tumor"))+
  geom_line(aes(y=AUC.lung, color="AUC.lung"))+
  geom_point(aes(y=AUC.lung, color="AUC.lung"))+
  ylim(c(0.4, 0.85))+
  scale_x_continuous(breaks = seq(0,100,10))+
  xlab("% Top Genes")+
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


#####plot best AUC among each category of sample

## Tumor samples (Supp Fig 10B)

#build a df out of excel output that generated from random forest analysis 
RF_Tumor.samples_genes_best_auc <- read_excel("best.AUC.genes.Tumor.xlsx")

#plot bar
ggplot(RF_Tumor.samples_genes_best_auc, aes(x=one_per, y=fct_reorder(Gene, one_per, max)))+
  geom_col(position = position_dodge(0.9), fill="red", color="black", alpha = 0.7, width = 0.5)+
  ylab("")+
  xlab("Gini Index")+
  ggplot2::ggtitle("Top Gens Achieving best AUC in Tumor Samples")+
  theme_bw()+
  theme(axis.text.y = element_markdown(size =20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 26, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))

## Lung samples (Supp Fig 10C)

#build a df out of csv output (or just pull the results above if saved in specific vector)
RF_Lung.samples_genes_best_auc <- read_excel("best.AUC.genes.Lung.xlsx")

#given 450 + genes, we will choose to display the top 50 
top50.genes <- RF_Lung.samples_genes_best_auc %>% 
  slice(1:50)

#plot bar
ggplot(top50.genes, aes(x=ten_per, y=fct_reorder(Gene, ten_per, max)))+
  geom_col(position = position_dodge(0.9), fill="dodgerblue", color="black", alpha = 0.7, width = 0.5)+
  ylab("")+
  xlab("Gini Index")+
  ggplot2::ggtitle("Top Gens Achieving best AUC in Lung Samples")+
  theme_bw()+
  theme(axis.text.y = element_markdown(size =20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 26, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))
