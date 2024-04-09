
pac=c("qiime2R","magrittr","tibble","tidyr","vegan","kableExtra",
      'survival','OMiSA','massMap',"ape",
      "compositions","phyloseq","table1","ggcorrplot","reshape2","GGally","ggpubr",
      "caret", "glmnet", "randomForest")

sapply(pac, require, character.only=T)

# Merge to phyloseqversion


meta.data=read.delim("Surgical.Cohort.Map.txt",row.names = 1)

SVs<-read_qza(paste("16S.Microbiome","table.filtered.qza", sep="/"))
otutable=SVs$data
otutable=otutable[,rownames(meta.data)]

taxonomy<-read_qza(paste("16S.Microbiome","taxonomy.filtered.qza", sep="/"))
#taxonomy$uuid
taxtable<-taxonomy$data %>% as.tibble() %>% separate(Taxon, sep="; ",
                                                     c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
taxtable=as.matrix(taxtable)
rownames(taxtable)=taxtable[,"Feature.ID"]
taxtable=taxtable[,-c(1,9)]

otutable=otutable[rownames(taxtable),]

tree<-read_qza(paste("16S.Microbiome","rooted-tree_quality.filtered.qza", sep="/"))

OTU = otu_table(otutable, taxa_are_rows = TRUE)
TAX = tax_table(taxtable, errorIfNULL=TRUE)
SAM = sample_data(meta.data)
## Create phyloseq object using OTU, TAX, SAM, and tree
phy = phyloseq(OTU, TAX, SAM, tree$data)




phy = prune_samples(sample_data(phy)$Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"), phy)

## 91 subjects
meta=read.csv("LNSmetafile_91sub.csv")
phy=prune_samples(sample_data(phy)$Subject_ID2 %in% meta$Subject_ID2, phy)


## filtering criteria
Genus.rel.table = transformSampleCounts(phy, function(x) x/sum(x))

Genus.rel.table.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.002), A = 0.01 * nsamples(Genus.rel.table))
Genus.table.com = prune_taxa(Genus.rel.table.wh1, Genus.rel.table)

phy.count=prune_taxa(Genus.rel.table.wh1, phy)


# Survival analysis
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
  
  if(j=="Lung.Tissue.In") title.lab="Involved: Suvival analysis for Recurrence on alpha diversity" else title.lab = "UnInvolved: Suvival analysis for Recurrence on alpha diversity" 
  
  cat('\n\n')
  
  res1=kbl(cox.results,
           row.names =TRUE, caption =title.lab) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  
  print(res1)
  
  cat('\n\n')
}


## ASV level

Contam=read.csv("Contaminant.list.csv",row.names = 1)

#Genus.table.com

taxa.name=tax_table(Genus.table.com)[,2:7]
taxa.name=as.character(apply(taxa.name, 1, function(x) {xx=as.character(unlist(x));
return(paste(xx[1:6],collapse = ";")) }))
tax_table(Genus.table.com)[,7]=taxa.name

RA.com=otu_table(Genus.table.com)
RA.com.clr=t(clr(t(RA.com)))

meta=sample_data(Genus.table.com)

for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  
  obstime <-as.numeric(unlist(tem[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(tem[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  RA.tem=RA.com.clr[,rownames(tem)]
  RA.tem.mean=RA.com[,rownames(tem)]
  
  
  cox.results=NULL
  for(ii in 1:nrow(RA.tem)) {
    
    fit.surv <- tryCatch({coxph(Surv(obstime, delta) ~ xx,
                                data=data.frame(obstime,delta, xx=scale(as.numeric(RA.tem[ii,]))))},
                         error = function(e) {NA})
    
    if(length(fit.surv)==1)   cox.results=rbind(cox.results,
                                                c(rownames(RA.tem)[ii], formatC(mean(as.numeric(RA.tem.mean[ii,delta==1])), digits = 2), formatC(mean(as.numeric(RA.tem.mean[ii,delta==0])), digits = 2), rep(NA,3))) else   cox.results=rbind(cox.results,c(rownames(RA.tem)[ii], formatC(mean(as.numeric(RA.tem.mean[ii,delta==1])), digits = 2), formatC(mean(as.numeric(RA.tem.mean[ii,delta==0])), digits = 2), round(summary(fit.surv)$coefficients[2],3),
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
  
  write.csv(cox.results, file=paste(j,"COX_PH_ASV.csv", sep="_"))
  
}

# MiRKAT-S test
## ASV level

### mean between In and UnIn samples
sampID=as.character(sample_data(phy.count)[["Subject_ID0"]])
rowind=which(!duplicated(sampID))
names(rowind)=sampID[rowind]

obstime <- as.numeric(as.character(sample_data(phy.count)[['New_time_followup.or.death']]))[rowind]
delta <- 3-unlist(lapply(strsplit(as.character(sample_data(phy.count)[['Progression_Lab_Inv']]),'\\.'),length))[rowind]
otu.tab = t(otu_table(phy.count))
tree = phy_tree(phy.count)
tax.tab = tax_table(phy.count)

############### prepare data
otu.tab.mean=t(sapply(rowind, function(rowindtmp) {
  otutmp=otu.tab[rownames(sample_data(phy))[sampID==sampID[rowindtmp]],]
  if(nrow(otutmp)==1){otutmp}else{colMeans(otutmp)}
}))


################ test
set.seed(1234)

Res.OMiSA=unlist(OMiSA(obstime, delta, otu.tab.mean, total.reads=NULL,tree, cov=NULL,
                       pow=c(1/4,1/3,1/2,1),g.unif.alpha=c(0.50), n.perm=10000))



for(j in c('Lung.Tissue.UnIn','Lung.Tissue.In')) {
  
  phy0=subset_samples(phy.count, Sample_Type_Involved == j)
  
  obstime <- as.numeric(as.character(sample_data(phy0)[['New_time_followup.or.death']]))
  delta <- 3-unlist(lapply(strsplit(as.character(sample_data(phy0)[['Progression_Lab_Inv']]),'\\.'),length))
  otu.tab = t(otu_table(phy0))
  tree = phy_tree(phy0)
  tax.tab = tax_table(phy0)
  
  set.seed(1234)
  Res.OMiSA=rbind(Res.OMiSA, unlist(OMiSA(obstime, delta, otu.tab, total.reads=NULL,tree, cov=NULL,
                                          pow=c(1/4,1/3,1/2,1),g.unif.alpha=c(0.50), n.perm=10000)))
  
}

Res.OMiSA=Res.OMiSA[,c(5,10)]
rownames(Res.OMiSA)=c("Average",'Lung.Tissue.UnIn','Lung.Tissue.In')

kbl(Res.OMiSA,
    row.names =TRUE, caption ="Associations between microbial profile and OS in community level(1000 permutations)") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))



# Two-stage microbial association mapping framework (massMap)

for(j in c('Lung.Tissue.UnIn','Lung.Tissue.In')) {
  
  phy0=subset_samples(phy.count, Sample_Type_Involved == j)
  
  obstime <- as.numeric(as.character(sample_data(phy0)[['New_time_followup.or.death']]))
  delta <- 3-unlist(lapply(strsplit(as.character(sample_data(phy0)[['Progression_Lab_Inv']]),'\\.'),length))
  otu.tab = t(otu_table(phy0))
  tree = phy_tree(phy0)
  tax.tab = tax_table(phy0)
  
  set.seed(1234)
  
  twoStageTest=massMap(X=NULL, obstime=obstime, delta=delta, otu.tab=otu.tab, tax.tab=tax.tab, tree=tree,
                                  outcome.trait = 'survival',screening.rank= 'Family',
                                  target.rank= 'OTU', alpha=0.2,n.perm=1e4)



  res=twoStageTest$res.target
  res0=twoStageTest$res.screening

  write.csv(res, file=paste(j,"massMap_target_otu.csv",sep="_"))
  write.csv(res0, file=paste(j,"massMap_screening_otu.csv",sep="_"))
  

}


# AUC analysis
## Data input and pre-processing

load("phy.RData")

phy = prune_samples(sample_data(phy)$Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"), phy)

# Remove taxa with 0 abundance
phy = subset_taxa(phy, rowSums(otu_table(phy)) != 0)
phy = subset_taxa(phy, Kingdom=="k__Bacteria")

## filtering criteria
Genus.rel.table = transformSampleCounts(phy, function(x) x/sum(x))

Genus.rel.table.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.002), A = 0.01 * nsamples(Genus.rel.table))
Genus.table.com = prune_taxa(Genus.rel.table.wh1, Genus.rel.table)

Genus.count = prune_taxa(Genus.rel.table.wh1, phy)

rm(phy); rm(Genus.rel.table); rm(Genus.rel.table.wh1)

### host transcriptome
RNA=readRDS("C:/Users/wangc22/OneDrive - NYU Langone Health/Segal/New.Early.Stage.Lung.Cancer.Data/Logistic_regression/data/rdata/matched_rm5years/tumor_fpkm_12vs0.RDS")

Genus.count <- prune_samples(sample_data(Genus.count)$Subject_ID2 %in% rownames(RNA), Genus.count)
Genus.table.com <- prune_samples(sample_data(Genus.table.com)$Subject_ID2 %in% rownames(RNA), Genus.table.com)

rm(RNA)

### For relative abundance data

for(nn in c("Normal", "Tumor")) {
  
  if(nn=="Normal") { Genus <- prune_samples(sample_data(Genus.table.com)$Sample_Type_Involved=="Lung.Tissue.UnIn", Genus.table.com)
  
  RNA=readRDS("C:/Users/wangc22/OneDrive - NYU Langone Health/Segal/New.Early.Stage.Lung.Cancer.Data/Logistic_regression/data/rdata/matched_rm5years/control_fpkm_12vs0.RDS")
  } else {
    
    Genus <- prune_samples(sample_data(Genus.table.com)$Sample_Type_Involved=="Lung.Tissue.In", Genus.table.com)
    
    RNA=readRDS("C:/Users/wangc22/OneDrive - NYU Langone Health/Segal/New.Early.Stage.Lung.Cancer.Data/Logistic_regression/data/rdata/matched_rm5years/tumor_fpkm_12vs0.RDS")
  }
  
  meta.used=sample_data(Genus)
  
  RNA=RNA[meta.used$Subject_ID2,]
  meta.used$delta <- RNA$recurrence
  Genus.used=otu_table(Genus)
  
  
  taxa.nam=rownames(Genus.used)
  Genus.used=t(Genus.used)
  colnames(Genus.used)=paste0("F", 1:length(taxa.nam))
  names(taxa.nam)=paste0("F", 1:length(taxa.nam))
  
  
  All.auc=NULL
  
  for(kk in 1:5) {
    
    folds <- createFolds(meta.used$delta, k = 10, list = FALSE)
    
    predicted.pro=NULL
    
    for(ii in 1:10) {
      
      Genus.training=Genus.used[folds!=ii,]
      RNA.training=RNA[folds!=ii,-1]
      meta.training=meta.used[folds!=ii,]
      
      
      ### host genes
      
      ## select top 200 genes according to the median absolute deviation 
      MAD=apply(RNA.training, 2, function(x) {xx=as.numeric(x); return(median(abs(xx-median(xx))))})
      RNA.training=RNA.training[,order(MAD,decreasing = TRUE)[1:200]]
      
      ## the inner CV regarding tuning parameters
      tunings=NULL
      for(alpha0 in (seq(0,20,1)/20)^2) {
        
        cvfit <- cv.glmnet(x=as.matrix(RNA.training), y=meta.training$delta, family="binomial",alpha=alpha0)
        tunings=rbind(tunings, c(alpha0, cvfit$cvm[cvfit$lambda==cvfit$lambda.min])) }  
      
      colnames(tunings)=c("alpha","cvm")
      tunings=data.frame(tunings)
      tunings=tunings[which.min(tunings$cvm),]
      
      
      ## the optimal model based on the training
      fit <- cv.glmnet(x=as.matrix(RNA.training), y=meta.training$delta, family = "binomial",
                       alpha=tunings$alpha)
      
      ## the predicted probabilities for the test
      ss.host=predict(fit, newx = as.matrix(RNA[folds==ii,colnames(RNA.training)]), s = "lambda.min",
                      type = "response") %>% as.data.frame()
      
      
      ### microbiome
      
      ## the inner CV regarding tuning parameters
      tunings=NULL
      for(alpha0 in (seq(0,20,1)/20)^2) {
        
        cvfit <- cv.glmnet(x=as.matrix(Genus.training), y=meta.training$delta, family="binomial",alpha=alpha0)
        tunings=rbind(tunings, c(alpha0, cvfit$cvm[cvfit$lambda==cvfit$lambda.min])) }  
      
      colnames(tunings)=c("alpha","cvm")
      tunings=data.frame(tunings)
      tunings=tunings[which.min(tunings$cvm),]
      
      
      ## the optimal model based on the training
      fit <- cv.glmnet(x=as.matrix(Genus.training), y=meta.training$delta, family = "binomial",
                       alpha=tunings$alpha)
      
      ## the predicted probabilities for the test
      ss.microbiome=predict(fit, newx = as.matrix(Genus.used[folds==ii,colnames(Genus.training)]), s = "lambda.min",
                            type = "response") %>% as.data.frame()
      
      
      predicted.pro=rbind(predicted.pro,cbind(ss.host,ss.microbiome))
      
    }
    
    predicted.pro=predicted.pro[meta.used$Subject_ID2,]
    
    datacom=data.frame(delta=meta.used$delta, host=as.numeric(predicted.pro[,1]), microbiome=as.numeric(predicted.pro[,2]))
    modelcom=glm(delta~host+microbiome,datacom, family="binomial")
    
    roc.host=ci(roc(datacom$delta~datacom$host))
    roc.microbiome=ci(roc(datacom$delta~datacom$microbiome))
    roc.both=ci(roc(datacom$delta~as.numeric(predict(modelcom,datacom[,-1],type = "response"))))
    
    All.auc=rbind(All.auc,c(roc.host[c(2,1,3)],roc.microbiome[c(2,1,3)],
                            roc.both[c(2,1,3)]))
    
    print(All.auc[,c(1,4,7)])
    
  }
  
  write.csv(All.auc, file=paste(c(nn,"RA.csv"),collapse = "_"))
}


#CCA for Involved samples

load("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/phy.RData")

phy = prune_samples(sample_data(phy)$Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"), phy)

meta=read.csv("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/LNSmetafile_91sub.csv")
phy=prune_samples(sample_data(phy)$Subject_ID2 %in% meta$Subject_ID2, phy)


## filtering criteria
Genus.rel.table = transformSampleCounts(phy, function(x) x/sum(x))

Genus.rel.table.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.001), A = 0.01 * nsamples(Genus.rel.table))
Genus.table.com = prune_taxa(Genus.rel.table.wh1, Genus.rel.table)

phy.count=prune_taxa(Genus.rel.table.wh1, phy)


rm(phy); rm(phy.count); rm(Genus.rel.table); rm(Genus.rel.table.wh1)


RNA=read.delim("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/RNA.host.Transcriptome/Copy of KO_Pathway_broken_Adenocarcinoma_Samples_sum_path03_JT.lns.txt", check.names = FALSE)

pathnames=RNA$path_03
names(pathnames)=paste0("f", 1:length(pathnames))

AA=sample_data(Genus.table.com)
AA.rownames=rownames(AA)

RNA=RNA[,AA.rownames]
rownames(RNA)=names(pathnames)
RNA=RNA[apply(RNA, 1, function(x) sum(x>=200)>=3),]

## AVSs having differential abundance based on edger
Taxa.sig=read_excel("Microbiome.Transcriptome.LUAD_FD_011624 JT updated.xlsx", sheet = 5, skip = 1)

Taxa.sig=Taxa.sig[!is.na(Taxa.sig$`potential.contaminant*`),]

Taxa.sig=as.character(
  sapply(Taxa.sig$`Taxa (name_ASV)`, function(x) {
    xx=unlist(strsplit(x, "_"))
    return(xx[length(xx)])
  } ))

Genus.sig.tumor = prune_taxa(rownames(tax_table(Genus.table.com)) %in% Taxa.sig, Genus.table.com)
Genus.sig.tumor= prune_samples(sample_data(Genus.sig.tumor)$Sample_Type_Involved %in% c("Lung.Tissue.In"), Genus.sig.tumor)

RNA.tumor=RNA[,rownames(sample_data(Genus.sig.tumor))]


Taxa.sig=read_excel("Microbiome.Transcriptome.LUAD_FD_011624 JT updated.xlsx", sheet = 6, skip = 1)
Taxa.sig=Taxa.sig[!is.na(Taxa.sig$`potential.contaminant*`),]

Taxa.sig=as.character(
  sapply(Taxa.sig$`Taxa (name_ASV)`, function(x) {
    xx=unlist(strsplit(x, "_"))
    return(xx[length(xx)])
  } ))

Genus.sig.lung = prune_taxa(rownames(tax_table(Genus.table.com)) %in% Taxa.sig, Genus.table.com)
Genus.sig.lung= prune_samples(sample_data(Genus.sig.lung)$Sample_Type_Involved %in% c("Lung.Tissue.UnIn"), Genus.sig.lung)

RNA.lung=RNA[,rownames(sample_data(Genus.sig.lung))]

## Genus.sig.tumor; RNA.tumor

Genus.table.com=Genus.sig.tumor; RNA=RNA.tumor

taxa.name=tax_table(Genus.table.com)[,2:7]
taxa.name=as.character(apply(taxa.name, 1, function(x) {xx=as.character(unlist(x));
return(paste(xx[1:6],collapse = ";")) }))
tax_table(Genus.table.com)[,7]=taxa.name

RA.com=otu_table(Genus.table.com)

meta=sample_data(Genus.table.com)


transcriptome=apply(RNA,1,function(x){
  if(min(x)==0) x=log2(x+1) else x=log2(x)
  return(x)})

microbiome.data=t(RA.com)

colnames(microbiome.data)=taxa.name


biomarker=vector("list", length = 2)
names(biomarker)=unlist(unique(meta[,"Progression_Lab_Inv"]))
for(j in unlist(unique(meta[,"Progression_Lab_Inv"]))) {

   print(j)

    tem.microbiome=microbiome.data[meta$Progression_Lab_Inv==j,];
    taxa.0=which(colSums(tem.microbiome)==0)
    if(length(taxa.0)>0) tem.microbiome=tem.microbiome[,-c(taxa.0)]

    tem.transcriptome=transcriptome[meta$Progression_Lab_Inv==j,]
    transcriptome.0=which(colSums(tem.transcriptome)==0)
    if(length(transcriptome.0)>0) tem.transcriptome=tem.transcriptome[,-c(transcriptome.0)]

#start_time <- Sys.time()

perm.out=CCA.permute(tem.microbiome,tem.transcriptome,trace = FALSE,nperms=100,
                     penaltyxs=seq(0,1,len=100), penaltyzs=seq(0,1,len=100))


#print(perm.out)
CCA.results=CCA(tem.microbiome,
                tem.transcriptome,
                typex ="standard",typez ="standard",
                K=3, trace=FALSE,
                penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
                v=perm.out$v.init)

  microbiome.PC=CCA.results$u
  transcriptome.PC=CCA.results$v
  estimated.cov=CCA.results$cors
  colnames(microbiome.PC)=colnames(transcriptome.PC)=names(estimated.cov)=paste0("PC",1:3)

  microbiome.PC=data.frame(taxa=colnames(tem.microbiome),microbiome.PC)
 transcriptome.PC=data.frame(gene=colnames(tem.transcriptome),transcriptome.PC)


 biomarker[[j]]=list(microbiome.PC,transcriptome.PC)

}

save(biomarker,file="CCA_biomarker_Tumor.RData")


load("CCA_biomarker_Tumor.RData")

for(j in names(biomarker)) {
  
  transcriptome.load=biomarker[[j]][[2]][,1:3]
  transcriptome.load=transcriptome.load[order(abs(as.numeric(transcriptome.load[,2])),decreasing = TRUE),]
  
  transcriptome.load$gene=pathnames[transcriptome.load$gene]
  
  write.csv(transcriptome.load, file=paste0(j,"trans_2024.csv"))
}

Contam=read.csv("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/contamlist_new.csv")


for(j in names(biomarker)) {
  
  microbiome.load=biomarker[[j]][[1]][,1:3]
  
  
  tax.name=tax_table(Genus.table.com)
  lineage=as.character(tax.name[microbiome.load$taxa,7])
  
  microbiome.load=data.frame(microbiome.load,lineage)
  
  rownames(microbiome.load)=microbiome.load$taxa
  
  contam1=Contam[Contam$asv %in%rownames(microbiome.load), ]
  
  microbiome.load=data.frame(microbiome.load[contam1$asv,],contaminant=contam1$contaminant)
  
  microbiome.load=microbiome.load[order(abs(as.numeric(microbiome.load[,2])),decreasing = TRUE),]
  
  write.csv(microbiome.load, file=paste0(j,"micro_2024.csv"))
  
}


load("CCA_biomarker_Tumor.RData")

microbiome.scale=transcriptome.scale=NULL
for (j in unlist(unique(meta[,"Progression_Lab_Inv"]))) {
  
  tem.microbiome=microbiome.data[meta$Progression_Lab_Inv==j,];
  tem.transcriptome=transcriptome[meta$Progression_Lab_Inv==j,]
  
  tem.microbiome=apply(tem.microbiome,2, function(x) if(max(x)==0) return(rep(0,length(x))) else return(scale(x)))
  tem.transcriptome=apply(tem.transcriptome,2, function(x) if(max(x)==0) return(rep(0,length(x))) else return(scale(x)))
  
  rownames(tem.microbiome)=rownames(microbiome.data)[meta$Progression_Lab_Inv==j]
  rownames(tem.transcriptome)=rownames(transcriptome)[meta$Progression_Lab_Inv==j]
  
  microbiome.scale=rbind(microbiome.scale, tem.microbiome)
  transcriptome.scale=rbind(transcriptome.scale, tem.transcriptome)
}

Res.corr=NULL
for(j in names(biomarker)) {
  
  microbiome.load=biomarker[[j]][[1]]
  transcriptome.load=biomarker[[j]][[2]]
  
  microbiome.PC=microbiome.scale[,microbiome.load$taxa] %*% apply(microbiome.load[,-1],2, as.numeric)
  transcriptome.PC=transcriptome.scale[,transcriptome.load$gene] %*% apply(transcriptome.load[,-1],2, as.numeric)
  
  for (cc in c(j, setdiff(names(biomarker),j))) {
    
    aa1=cor.test(microbiome.PC[rownames(meta)[meta$Progression_Lab_Inv==cc],1],transcriptome.PC[rownames(meta)[meta$Progression_Lab_Inv==cc],1])
    aa2=cor.test(microbiome.PC[rownames(meta)[meta$Progression_Lab_Inv==cc],2],transcriptome.PC[rownames(meta)[meta$Progression_Lab_Inv==cc],2])
    
    Res.corr=rbind(Res.corr, c(j,cc,aa1$estimate,aa1$p.value,aa2$estimate,aa2$p.value))
  }}


colnames(Res.corr)=c("Sample for identification","Sample for correlation test", "PC1_correlation","PC1_pval",
                     "PC2_correlation","PC2_pval")

Res.corr=data.frame(Res.corr)



Res.corr[,3:6]=apply(Res.corr[,3:6],2,function(x) formatC(as.numeric(x), digits = 2))

res1=kbl(Res.corr,
         row.names =FALSE, caption ="Correlation test for the identified components integrated by microbiome and transcriptome") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))


for(j in names(biomarker)) {
  
  microbiome.load=biomarker[[j]][[1]]
  transcriptome.load=biomarker[[j]][[2]]
  
  microbiome.PC=microbiome.scale[,microbiome.load$taxa] %*% apply(microbiome.load[,-1],2, as.numeric)
  transcriptome.PC=transcriptome.scale[,transcriptome.load$gene] %*% apply(transcriptome.load[,-1],2, as.numeric)
  
  obstime <-as.numeric(unlist(meta[,"New_time_followup.or.death"]))
  delta <- 3-unlist(lapply(strsplit(unlist(meta[,"Progression_Lab_Inv"]),'\\.'),length))
  
  if(unlist(strsplit(j, "\\."))[[1]]=="In") sample.plot="Involved" else sample.plot="UnInvolved" 
  
    
    data.ph=data.frame(Involved=revalue(as.character(meta$Lung_Tissue_Involved), c("0"="UnInvolved", "1"="Involved")),
                       obstime,delta, PC.cor=microbiome.PC[rownames(meta),1]*transcriptome.PC[rownames(meta),1])
    
    tem.data=data.ph[data.ph$Involved==sample.plot,]
    
    res.cat=tem.data
    
    res.cat$PC.cor=rep("low", nrow(tem.data))    
    res.cat$PC.cor[tem.data$PC.cor>mean(tem.data$PC.cor)]="high"
    
    
    # 4. Fit survival curves and visualize
    if(j== "In.Recurrence") res.cat$PC.cor=mapvalues(res.cat$PC.cor, 
                                                     from = c("high", "low"), to = c("High risk", "Low risk")) else  res.cat$PC.cor = mapvalues(res.cat$PC.cor, from = c("high", "low"), to = c("Low risk", "High risk")) 
    
    aa=ggsurvplot(
      fit = survfit(Surv(obstime, delta) ~ PC.cor,data=res.cat),
      xlab = "Days", 
      ylab = "Recurrence probability",
      legend.title = "",
      legend.labs = c("High risk", "Low risk"),
      palette = c("red","blue"),
      risk.table = TRUE,risk.table.title=NA, tables.height=0.4,
      conf.int = TRUE, pval = TRUE)+ggtitle(paste0("Identification: ", j))
    
    print(aa)    
    
    
}

#CCA for UnInvolved samples


load("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/phy.RData")

phy = prune_samples(sample_data(phy)$Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"), phy)

meta=read.csv("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/LNSmetafile_91sub.csv")
phy=prune_samples(sample_data(phy)$Subject_ID2 %in% meta$Subject_ID2, phy)


## filtering criteria
Genus.rel.table = transformSampleCounts(phy, function(x) x/sum(x))

Genus.rel.table.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.001), A = 0.01 * nsamples(Genus.rel.table))
Genus.table.com = prune_taxa(Genus.rel.table.wh1, Genus.rel.table)

phy.count=prune_taxa(Genus.rel.table.wh1, phy)

rm(phy); rm(phy.count); rm(Genus.rel.table); rm(Genus.rel.table.wh1)


RNA=read.delim("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/RNA.host.Transcriptome/Copy of KO_Pathway_broken_Adenocarcinoma_Samples_sum_path03_JT.lns.txt", check.names = FALSE)

pathnames=RNA$path_03
names(pathnames)=paste0("f", 1:length(pathnames))

AA=sample_data(Genus.table.com)
AA.rownames=rownames(AA)

RNA=RNA[,AA.rownames]
rownames(RNA)=names(pathnames)
RNA=RNA[apply(RNA, 1, function(x) sum(x>=200)>=3),]


Taxa.sig=read_excel("Microbiome.Transcriptome.LUAD_FD_011624 JT updated.xlsx", sheet = 5, skip = 1)

Taxa.sig=Taxa.sig[!is.na(Taxa.sig$`potential.contaminant*`),]

Taxa.sig=as.character(
  sapply(Taxa.sig$`Taxa (name_ASV)`, function(x) {
    xx=unlist(strsplit(x, "_"))
    return(xx[length(xx)])
  } ))

Genus.sig.tumor = prune_taxa(rownames(tax_table(Genus.table.com)) %in% Taxa.sig, Genus.table.com)
Genus.sig.tumor= prune_samples(sample_data(Genus.sig.tumor)$Sample_Type_Involved %in% c("Lung.Tissue.In"), Genus.sig.tumor)

RNA.tumor=RNA[,rownames(sample_data(Genus.sig.tumor))]



Taxa.sig=read_excel("Microbiome.Transcriptome.LUAD_FD_011624 JT updated.xlsx", sheet = 6, skip = 1)

Taxa.sig=Taxa.sig[!is.na(Taxa.sig$`potential.contaminant*`),]

Taxa.sig=as.character(
  sapply(Taxa.sig$`Taxa (name_ASV)`, function(x) {
    xx=unlist(strsplit(x, "_"))
    return(xx[length(xx)])
  } ))

Genus.sig.lung = prune_taxa(rownames(tax_table(Genus.table.com)) %in% Taxa.sig, Genus.table.com)
Genus.sig.lung= prune_samples(sample_data(Genus.sig.lung)$Sample_Type_Involved %in% c("Lung.Tissue.UnIn"), Genus.sig.lung)

RNA.lung=RNA[,rownames(sample_data(Genus.sig.lung))]

Genus.table.com=Genus.sig.lung; RNA=RNA.lung

taxa.name=tax_table(Genus.table.com)[,2:7]
taxa.name=as.character(apply(taxa.name, 1, function(x) {xx=as.character(unlist(x));
return(paste(xx[1:6],collapse = ";")) }))
tax_table(Genus.table.com)[,7]=taxa.name

RA.com=otu_table(Genus.table.com)

meta=sample_data(Genus.table.com)


transcriptome=apply(RNA,1,function(x){
  if(min(x)==0) x=log2(x+1) else x=log2(x)
  return(x)})

microbiome.data=t(RA.com)

colnames(microbiome.data)=taxa.name
biomarker=vector("list", length = 2)
names(biomarker)=unlist(unique(meta[,"Progression_Lab_Inv"]))
for(j in unlist(unique(meta[,"Progression_Lab_Inv"]))) {

   print(j)

    tem.microbiome=microbiome.data[meta$Progression_Lab_Inv==j,];
    taxa.0=which(colSums(tem.microbiome)==0)
    if(length(taxa.0)>0) tem.microbiome=tem.microbiome[,-c(taxa.0)]

    tem.transcriptome=transcriptome[meta$Progression_Lab_Inv==j,]
    transcriptome.0=which(colSums(tem.transcriptome)==0)
    if(length(transcriptome.0)>0) tem.transcriptome=tem.transcriptome[,-c(transcriptome.0)]

#start_time <- Sys.time()

perm.out=CCA.permute(tem.microbiome,tem.transcriptome,trace = FALSE,nperms=100,
                     penaltyxs=seq(0,1,len=100), penaltyzs=seq(0,1,len=100))


#print(perm.out)
CCA.results=CCA(tem.microbiome,
                tem.transcriptome,
                typex ="standard",typez ="standard",
                K=3, trace=FALSE,
                penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
                v=perm.out$v.init)

  microbiome.PC=CCA.results$u
  transcriptome.PC=CCA.results$v
  estimated.cov=CCA.results$cors
  colnames(microbiome.PC)=colnames(transcriptome.PC)=names(estimated.cov)=paste0("PC",1:3)

  microbiome.PC=data.frame(taxa=colnames(tem.microbiome),microbiome.PC)
 transcriptome.PC=data.frame(gene=colnames(tem.transcriptome),transcriptome.PC)


 biomarker[[j]]=list(microbiome.PC,transcriptome.PC)

}

save(biomarker,file="CCA_biomarker_Lung.RData")


load("CCA_biomarker_Lung.RData")

for(j in names(biomarker)) {
  
  transcriptome.load=biomarker[[j]][[2]][,1:3]
  transcriptome.load=transcriptome.load[order(abs(as.numeric(transcriptome.load[,2])),decreasing = TRUE),]
  
  transcriptome.load$gene=pathnames[transcriptome.load$gene]
  
  write.csv(transcriptome.load, file=paste0(j,"trans_2024.csv"))
}

Contam=read.csv("C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data/contamlist_new.csv")


for(j in names(biomarker)) {
  
  microbiome.load=biomarker[[j]][[1]][,1:3]
  
  
  tax.name=tax_table(Genus.table.com)
  lineage=as.character(tax.name[microbiome.load$taxa,7])
  
  microbiome.load=data.frame(microbiome.load,lineage)
  
  rownames(microbiome.load)=microbiome.load$taxa
  
  contam1=Contam[Contam$asv %in%rownames(microbiome.load), ]
  
  microbiome.load=data.frame(microbiome.load[contam1$asv,],contaminant=contam1$contaminant)
  
  microbiome.load=microbiome.load[order(abs(as.numeric(microbiome.load[,2])),decreasing = TRUE),]
  
  write.csv(microbiome.load, file=paste0(j,"micro_2024.csv"))
  
}

load("CCA_biomarker_Lung.RData")

microbiome.scale=transcriptome.scale=NULL
for (j in unlist(unique(meta[,"Progression_Lab_Inv"]))) {
  
  tem.microbiome=microbiome.data[meta$Progression_Lab_Inv==j,];
  tem.transcriptome=transcriptome[meta$Progression_Lab_Inv==j,]
  
  tem.microbiome=apply(tem.microbiome,2, function(x) if(max(x)==0) return(rep(0,length(x))) else return(scale(x)))
  tem.transcriptome=apply(tem.transcriptome,2, function(x) if(max(x)==0) return(rep(0,length(x))) else return(scale(x)))
  
  rownames(tem.microbiome)=rownames(microbiome.data)[meta$Progression_Lab_Inv==j]
  rownames(tem.transcriptome)=rownames(transcriptome)[meta$Progression_Lab_Inv==j]
  
  microbiome.scale=rbind(microbiome.scale, tem.microbiome)
  transcriptome.scale=rbind(transcriptome.scale, tem.transcriptome)
}


for(j in names(biomarker)) {
  
  microbiome.load=biomarker[[j]][[1]]
  transcriptome.load=biomarker[[j]][[2]]
  
  microbiome.PC=microbiome.scale[,microbiome.load$taxa] %*% apply(microbiome.load[,-1],2, as.numeric)
  transcriptome.PC=transcriptome.scale[,transcriptome.load$gene] %*% apply(transcriptome.load[,-1],2, as.numeric)
  
  obstime <-as.numeric(unlist(meta[,"New_time_followup.or.death"]))
  delta <- 3-unlist(lapply(strsplit(unlist(meta[,"Progression_Lab_Inv"]),'\\.'),length))
  
  if(unlist(strsplit(j, "\\."))[[1]]=="In") sample.plot="Involved" else sample.plot="UnInvolved" 

    data.ph=data.frame(Involved=revalue(as.character(meta$Lung_Tissue_Involved), c("0"="UnInvolved", "1"="Involved")),
                       obstime,delta, PC.cor=microbiome.PC[rownames(meta),1]*transcriptome.PC[rownames(meta),1])
    
    tem.data=data.ph[data.ph$Involved==sample.plot,]
    
    res.cat=tem.data
    
    res.cat$PC.cor=rep("low", nrow(tem.data))    
    res.cat$PC.cor[tem.data$PC.cor>mean(tem.data$PC.cor)]="high"
    
    
    # 4. Fit survival curves and visualize
    if(j== "UnIn.Recurrence") res.cat$PC.cor=mapvalues(res.cat$PC.cor, 
                                                       from = c("high", "low"), to = c("High risk", "Low risk")) else  res.cat$PC.cor = mapvalues(res.cat$PC.cor, from = c("high", "low"), to = c("Low risk", "High risk")) 
    
    aa=ggsurvplot(
      fit = survfit(Surv(obstime, delta) ~ PC.cor,data=res.cat),
      xlab = "Days", 
      ylab = "Recurrence probability",
      legend.title = "",
      legend.labs = c("High risk", "Low risk"),
      palette = c("red","blue"),
      risk.table = TRUE,risk.table.title=NA, tables.height=0.4,
      conf.int = TRUE, pval = TRUE)+ggtitle(paste0("Identification: ", j))
    
    
    print(aa)    
    

  }

