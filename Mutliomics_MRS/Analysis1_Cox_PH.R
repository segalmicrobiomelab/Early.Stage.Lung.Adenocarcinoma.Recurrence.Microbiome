########## Pre-processing

filedir="C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data"

setwd(filedir)

pac=c("qiime2R","magrittr","tibble","tidyr","vegan","kableExtra", "ComplexHeatmap",
      'survival','OMiSA','massMap',"ape",
      "DESeq2","phyloseq","ANCOMBC","table1","ggcorrplot","reshape2","GGally","ggpubr")

sapply(pac, require, character.only=T)

# Merge to phyloseqversion

meta.data=read.delim("Surgical.Cohort.Map.txt",row.names = 1)


SVs<-read_qza(paste("16S.Microbiome","table.filtered.qza", sep="/"))
otutable=SVs$data

setdiff(rownames(meta.data),colnames(otutable))
identical(sort(rownames(meta.data)),sort(colnames(otutable)))

otutable=otutable[,rownames(meta.data)]


taxonomy<-read_qza(paste("16S.Microbiome","taxonomy.filtered.qza", sep="/"))
#taxonomy$uuid
taxtable<-taxonomy$data %>% as.tibble() %>% separate(Taxon, sep="; ",
                                                     c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
taxtable=as.matrix(taxtable)
rownames(taxtable)=taxtable[,"Feature.ID"]
taxtable=taxtable[,-c(1,9)]

setdiff(rownames(taxtable),rownames(otutable))
identical(rownames(taxtable),rownames(otutable))

otutable=otutable[rownames(taxtable),]

identical(rownames(taxtable),rownames(otutable))


tree<-read_qza(paste("16S.Microbiome","rooted-tree_quality.filtered.qza", sep="/"))
#tree$uuid
## [1] "df4a8d84-2386-421c-86e7-eebe6412a982"

OTU = otu_table(otutable, taxa_are_rows = TRUE)
TAX = tax_table(taxtable, errorIfNULL=TRUE)
SAM = sample_data(meta.data)
## Create phyloseq object using OTU, TAX, SAM, and tree
phy = phyloseq(OTU, TAX, SAM, tree$data)
save(phy,file="phy.RData")

load("phy.RData")

phy = prune_samples(sample_data(phy)$Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"), phy)

# Remove taxa with 0 abundance
phy = subset_taxa(phy, rowSums(otu_table(phy)) != 0)
phy = subset_taxa(phy, Kingdom=="k__Bacteria")

## filtering criteria
Genus.rel.table = transformSampleCounts(phy, function(x) x/sum(x))

Genus.rel.table.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.002), A = 0.01 * nsamples(Genus.rel.table))
Genus.table.com = prune_taxa(Genus.rel.table.wh1, Genus.rel.table)

phy.count=prune_taxa(Genus.rel.table.wh1, phy)

# Analysis 1: Survival analysis
## Alpha diversity

set.seed(1234)
erDF = estimate_richness(phy.count , split = TRUE)
erDF =data.frame(erDF[,c("Observed","Shannon", "Simpson") ],sample_data(phy.count))


for(j in unique(erDF$Sample_Type_Involved)) {
  
  tem=erDF[erDF$Sample_Type_Involved==j,]
  
  obstime <-as.numeric(unlist(tem[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(tem[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  cox.results=NULL
  for(ii in c("Observed","Shannon", "Simpson")) {
    
    fit.surv <- tryCatch({coxph(Surv(obstime, delta) ~ xx,data=data.frame(tem,obstime, delta, xx=scale(tem[,ii])))},
                         error = function(e) {NA})
    cox.results=rbind(cox.results,c(round(summary(fit.surv)$coefficients[2],3),
                                    paste(round(summary(fit.surv)$conf.int[3:4],3),collapse = "-"),
                                    formatC(summary(fit.surv)$coefficients[5],digits = 3)))          
    
  }
  
  rownames(cox.results)=c("Observed","Shannon", "Simpson")
  colnames(cox.results)=c("Hazard ratio","95%CI","P-value")
  
  write.csv(cox.results, file=paste(j,"COX_PH_alpha.csv", sep="_"))
  
}

## ASV level

Contam=read.csv("Contaminant.list.csv",row.names = 1)

#Genus.table.com

taxa.name=tax_table(Genus.table.com)[,2:7]
taxa.name=as.character(apply(taxa.name, 1, function(x) {xx=as.character(unlist(x));
return(paste(xx[1:6],collapse = ";")) }))
tax_table(Genus.table.com)[,7]=taxa.name

RA.com=otu_table(Genus.table.com)

meta=sample_data(Genus.table.com)

for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  
  obstime <-as.numeric(unlist(tem[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(tem[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  RA.tem=RA.com[,rownames(tem)]
  
  
  
  cox.results=NULL
  for(ii in 1:nrow(RA.tem)) {
    
    fit.surv <- tryCatch({coxph(Surv(obstime, delta) ~ xx,
                                data=data.frame(obstime,delta, xx=scale(as.numeric(RA.tem[ii,]))))},
                         error = function(e) {NA})
    
    if(length(fit.surv)==1)   cox.results=rbind(cox.results,
                                                c(rownames(RA.tem)[ii], formatC(mean(as.numeric(RA.tem[ii,delta==1])), digits = 2), formatC(mean(as.numeric(RA.tem[ii,delta==0])), digits = 2), rep(NA,3))) else   cox.results=rbind(cox.results,c(rownames(RA.tem)[ii], formatC(mean(as.numeric(RA.tem[ii,delta==1])), digits = 2), formatC(mean(as.numeric(RA.tem[ii,delta==0])), digits = 2), round(summary(fit.surv)$coefficients[2],3),
                                                                                                                                                                                                                                                   paste(round(summary(fit.surv)$conf.int[3:4],3),collapse = "-"),
                                                                                                                                                                                                                                                   summary(fit.surv)$coefficients[5]))    
                                                
  }
  
  colnames(cox.results)=c("OTU","Mean_RA_Re","Mean_RA_NoRe", "Hazard ratio","CI","pvalue")
  cox.results=data.frame(cox.results)
  cox.results$pvalue=as.numeric(cox.results$pvalue)
  cox.results$lineage=as.character(tax_table(Genus.table.com)[cox.results$OTU,7])
  cox.results$adj.p=p.adjust(cox.results$pvalue,method = "BH")
  
  rownames(cox.results)=cox.results$OTU
  contam1=Contam[Contam$asv %in%cox.results$OTU, ]
  
  cox.results=data.frame(cox.results[contam1$asv,],contaminant=contam1$potential.contaminant)
  
  cox.results=cox.results[order(cox.results$pvalue),]
  
  write.csv(cox.results, file=paste(j,"COX_PH_OTU.csv", sep="_"))
  
}