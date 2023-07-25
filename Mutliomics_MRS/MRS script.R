
filedir="C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data"

setwd(filedir)

pac=c("qiime2R","magrittr","tibble","tidyr","vegan","kableExtra", "ComplexHeatmap",
      'survival','OMiSA','massMap',"ape","survminer","pROC","plyr",
      "DESeq2","phyloseq","ANCOMBC","table1","ggcorrplot","reshape2","GGally","ggpubr")

sapply(pac, require, character.only=T)


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

phy.count

rm(phy); rm(phy.count); rm(Genus.rel.table); rm(Genus.rel.table.wh1)


RNA=read.delim(paste("RNA.host.Transcriptome","NovaSeq_Raw_Counts_Gene_Symbols (1).txt", sep="/"), row.names = 1,check.names = FALSE)

AA=sample_data(Genus.table.com)
AA.rownames=rownames(AA)
AA.rownames[!grepl("NYU",AA.rownames)]=gsub("\\.Lung", "", AA.rownames[!grepl("NYU",AA.rownames)])
AA.rownames[AA.rownames=="NYU1923..Lung"]="NYU1923.Lung"

both=intersect(AA.rownames,colnames(RNA))

#identical(both, AA.rownames[AA.rownames %in% both])

RNA=RNA[,both]
Genus.table.com <- prune_samples(AA.rownames %in% both, Genus.table.com)

rm(AA.rownames); rm(AA)

RNA=RNA[apply(RNA, 1, function(x) sum(x>=200)>=3),]

colnames(RNA)=rownames(sample_data(Genus.table.com))

meta=sample_data(Genus.table.com)

#####Leave one out CV
transcriptome=apply(RNA,1,function(x){
  if(min(x)==0) x=log(x+1) else x=log(x)
  return(x)})
transcriptome=t(transcriptome)

rm(RNA);

for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  RA.tem=transcriptome[,rownames(tem)]
  
  Survival_all=Survival_CV_host(meta.table=tem, RA.table=RA.tem)
  save(Survival_all, file=paste(j, "Host_Survival_CV.RData", sep="_"))
  
}


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
  RA.tem=RA.com[,rownames(tem)]
  
  Survival_all=Survival_CV(meta.table=tem, RA.table=RA.tem,Contamlist=Contam)
  save(Survival_all, file=paste(j, "Micro_Survival_CV.RData", sep="_"))
}

Survival_CV_host=function(meta.table,RA.table) {
  
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  meta.table=data.frame(meta.table,obstime,delta)
  
  PH.CV=vector("list",nrow(meta.table))
  names(PH.CV)=rownames(meta.table)
  
  for(vv in 1:nrow(meta.table)) {
    
    print(vv)
    
    cox.results=NULL
    for(ii in 1:nrow(RA.table)) {
      
      fit.surv <- tryCatch({coxph(Surv(obstime, delta) ~ xx,
                                  data=data.frame(meta.table[-vv,], xx=scale(as.numeric(RA.table[ii,-vv]))))},
                           error = function(e) {NA})
      
      if(length(fit.surv)==1)   cox.results=rbind(cox.results, c(rownames(RA.table)[ii], NA, NA)) else   cox.results=rbind(cox.results,c(rownames(RA.table)[ii], summary(fit.surv)$coefficients[2] ,summary(fit.surv)$coefficients[5]))
    }
    
    colnames(cox.results)=c("gene","HR","pvalue")
    cox.results=data.frame(cox.results)
    cox.results$pvalue=as.numeric(cox.results$pvalue)
    rownames(cox.results)=cox.results$gene
    cox.results=cox.results[!is.na(cox.results$pvalue),]
    cox.results$prank=rank(cox.results$pvalue, ties.method = "min")
    
    PH.CV[[vv]]=cox.results }
  
  return(PH.CV)
}


Survival_CV=function(meta.table,RA.table,Contamlist) {
  
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  meta.table=data.frame(meta.table,obstime,delta)
  
  PH.CV=vector("list",nrow(meta.table))
  names(PH.CV)=rownames(meta.table)
  
  for(vv in 1:nrow(meta.table)) {
    
    print(vv)
    
    cox.results=NULL
    for(ii in 1:nrow(RA.table)) {
      
      fit.surv <- tryCatch({coxph(Surv(obstime, delta) ~ xx,
                                  data=data.frame(meta.table[-vv,], xx=scale(as.numeric(RA.table[ii,-vv]))))},
                           error = function(e) {NA})
      
      if(length(fit.surv)==1)   cox.results=rbind(cox.results, c(rownames(RA.table)[ii], NA)) else   cox.results=rbind(cox.results,c(rownames(RA.table)[ii], summary(fit.surv)$coefficients[5]))
    }
    
    colnames(cox.results)=c("OTU","pvalue")
    cox.results=data.frame(cox.results)
    cox.results$pvalue=as.numeric(cox.results$pvalue)
    rownames(cox.results)=cox.results$OTU
    contam1=Contamlist[Contamlist$asv %in%cox.results$OTU, ]
    
    cox.results=data.frame(cox.results[contam1$asv,],contaminant=contam1$potential.contaminant)
    
    cox.results=cox.results[!is.na(cox.results$pvalue),]
    cox.results$prank=rank(cox.results$pvalue, ties.method = "min")
    
    PH.CV[[vv]]=cox.results }
  
  return(PH.CV)
}

###### MRSs calculation and identify optimal cutoff

for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  RA.tem=RA.com[,rownames(tem)]
  
  load(paste(j, "Micro_Survival_CV.RData", sep="_"))
  
  MRS_all=MRS_CV(meta.table=tem,RA.table=RA.tem,Survival_res=Survival_all)
  
  
  meta.table=tem;
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  meta.table=data.frame(meta.table,obstime,delta)
  
  
  Scores=NULL
  for(vv in 1:length(MRS_all)) {
    
    print(vv)
    
    AUC_data=MRS_all[[vv]]
    
    AUC=sapply(unique(AUC_data$cutoff), function(x, AUC_data0,meta.table0,IDname) {
      
      AUC_tem=AUC_data0[AUC_data0$cutoff==x,]
      AUC1=data.frame(meta.table0[AUC_tem$ID,],AUC_tem)
      AUC1=AUC1[AUC1$ID != IDname, ]
      
      return(c(x,roc(delta~MRS_all, data=AUC1)$auc,
               roc(delta~MRS_NC, data=AUC1)$auc))
      
    }, AUC_data0=AUC_data,meta.table0=meta.table,IDname=names(MRS_all)[vv])
    
    AUC=t(AUC)
    
    colnames(AUC)=c("cutoff","MRS_all","MRS_NC")
    AUC=data.frame(AUC)
    
    MRS_all_cutoff=AUC$cutoff[which.max(AUC$MRS_all)];
    MRS_NC_cutoff=AUC$cutoff[which.max(AUC$MRS_NC)]
    
    Scores=rbind(Scores, c(AUC_data$MRS_all[AUC_data$cutoff==MRS_all_cutoff & AUC_data$ID==names(MRS_all)[vv]],
                           AUC_data$MRS_NC[AUC_data$cutoff==MRS_NC_cutoff & AUC_data$ID==names(MRS_all)[vv]], MRS_all_cutoff, MRS_NC_cutoff))
    
  }
  
  colnames(Scores)=c("MRS_all","MRS_NC","MRS_all_cutoff","MRS_NC_cutoff")
  Scores=data.frame(meta.table, Scores)
  
  write.csv(Scores, file=paste(j, "MRS1.csv", sep="_"))
  
  aa1=ci(roc(delta~MRS_all, data=Scores))
  aa2=ci(roc(delta~MRS_NC, data=Scores))
  
  AUC_pred=rbind(aa1[c(2,1,3)],aa2[c(2,1,3)])
  colnames(AUC_pred)=c("AUC","LowerCI","UpperCI")
  rownames(AUC_pred)=c("MRS_all","MRS_NC")
  AUC_pred=data.frame(AUC_pred)
  
  print(AUC_pred)
  
  write.csv(AUC_pred, file=paste(j, "prediction1.csv", sep="_"))
  
}


## based on the rank, number of top taxa
MRS_CV=function(meta.table,RA.table,Survival_res) {
  
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  meta.table=data.frame(meta.table,obstime,delta)
  
  AUC.CV=vector("list",nrow(meta.table))
  names(AUC.CV)=rownames(meta.table)
  
  for(vv in 1:nrow(meta.table)) {
    
    print(vv)
    
    cox.results=Survival_res[[vv]]
    
    AUC=lapply(5:max(cox.results$prank), function(tt,cox.results,RA.table,meta.table) {
      
      Shannon1=apply(RA.table[cox.results$OTU[cox.results$prank<=tt],], 2, function(x) {if(sum(x)==0) xx=rep(0,length(x)) else xx=x/sum(x); return(xx)})
      
      Shannon1=apply(Shannon1, 2, function(x) {if(sum(x)==0) xx=0 else xx=-sum((x*log(x))[x!=0]); return(xx)})
      
      Shannon2=apply(RA.table[cox.results$OTU[cox.results$prank<=tt & !(cox.results$contaminant)],], 2, function(x) {if(sum(x)==0) xx=rep(0,length(x)) else xx=x/sum(x); return(xx)})
      
      if(is.null(dim(Shannon2))) Shannon2=rep(0,length(Shannon2)) else Shannon2=apply(Shannon2, 2, function(x) {if(sum(x)==0) xx=0 else xx=-sum((x*log(x))[x!=0]); return(xx)})
      
      AUC_data=data.frame(ID=rownames(meta.table), cutoff=rep(tt, length(Shannon1)), Shannon1=as.numeric(Shannon1), Shannon2=as.numeric(Shannon2))
      
      return(AUC_data)
      
    }, cox.results=cox.results, RA.table=RA.table, meta.table=meta.table)
    
    AUC=do.call(rbind, AUC)
    
    colnames(AUC)=c("ID","cutoff","MRS_all","MRS_NC")
    AUC=data.frame(AUC)
    
    AUC.CV[[vv]]=AUC }
  
  return(AUC.CV)
}


for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  RA.tem=RA.com[,rownames(tem)]
  
  load(paste(j, "Micro_Survival_CV.RData", sep="_"))
  
  MRS_all=MRS_CV(meta.table=tem,RA.table=RA.tem,Survival_res=Survival_all)
  
  
  meta.table=tem;
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  meta.table=data.frame(meta.table,obstime,delta)
  
  
  Scores=NULL
  for(vv in 1:length(MRS_all)) {
    
    print(vv)
    
    AUC_data=MRS_all[[vv]]
    
    AUC=sapply(unique(AUC_data$cutoff), function(x, AUC_data0,meta.table0,IDname) {
      
      AUC_tem=AUC_data0[AUC_data0$cutoff==x,]
      AUC1=data.frame(meta.table0[AUC_tem$ID,],AUC_tem)
      AUC1=AUC1[AUC1$ID != IDname, ]
      
      return(c(x,roc(delta~MRS_all, data=AUC1)$auc,
               roc(delta~MRS_NC, data=AUC1)$auc))
      
    }, AUC_data0=AUC_data,meta.table0=meta.table,IDname=names(MRS_all)[vv])
    
    AUC=t(AUC)
    
    colnames(AUC)=c("cutoff","MRS_all","MRS_NC")
    AUC=data.frame(AUC)
    
    MRS_all_cutoff=AUC$cutoff[which.max(AUC$MRS_all)];
    MRS_NC_cutoff=AUC$cutoff[which.max(AUC$MRS_NC)]
    
    Scores=rbind(Scores, c(AUC_data$MRS_all[AUC_data$cutoff==MRS_all_cutoff & AUC_data$ID==names(MRS_all)[vv]],
                           AUC_data$MRS_NC[AUC_data$cutoff==MRS_NC_cutoff & AUC_data$ID==names(MRS_all)[vv]], MRS_all_cutoff, MRS_NC_cutoff))
    
  }
  
  colnames(Scores)=c("MRS_all","MRS_NC","MRS_all_cutoff","MRS_NC_cutoff")
  Scores=data.frame(meta.table, Scores)
  
  write.csv(Scores, file=paste(j, "MRS1.csv", sep="_"))
  
  aa1=ci(roc(delta~MRS_all, data=Scores))
  aa2=ci(roc(delta~MRS_NC, data=Scores))
  
  AUC_pred=rbind(aa1[c(2,1,3)],aa2[c(2,1,3)])
  colnames(AUC_pred)=c("AUC","LowerCI","UpperCI")
  rownames(AUC_pred)=c("MRS_all","MRS_NC")
  AUC_pred=data.frame(AUC_pred)
  
  print(AUC_pred)
  
  write.csv(AUC_pred, file=paste(j, "prediction1.csv", sep="_"))
  
}



for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  RA.tem=RNA[,rownames(tem)]
  
  load(paste(j, "Host_Survival_CV.RData", sep="_"))
  
  MRS_all=Host_CV(meta.table=tem,RA.table=RA.tem,Survival_res=Survival_all)
  
  save(MRS_all, file=paste(j, "Host_Score_CV.RData", sep="_"))
  
}

## based on the rank, number of top taxa
Host_CV=function(meta.table,RA.table,Survival_res) {
  
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  
  meta.table=data.frame(meta.table,obstime,delta)
  
  AUC.CV=vector("list",nrow(meta.table))
  names(AUC.CV)=rownames(meta.table)
  
  for(vv in 1:nrow(meta.table)) {
    
    print(vv)
    
    cox.results=Survival_res[[vv]]
    
    #
    AUC=lapply(5:max(cox.results$prank), function(tt,cox.results,RA.table,meta.table, IDname) {
      
      count=RA.table[cox.results$gene[cox.results$prank<=tt],]
      Positive.Index=as.numeric(cox.results$HR[cox.results$prank<=tt])>1
      
      Host1=apply(count, 2, function(x,Positive.Index) {S1=sum(x[Positive.Index]); S2=sum(x[!Positive.Index])
      S1=ifelse(S1>0,log(S1), 0); S2=ifelse(S2>0,log(S2), 0); return(c(S1-S2))}, Positive.Index=Positive.Index)
      
      
      AUC1=data.frame(meta.table, Host1=as.numeric(Host1))
      AUC1=AUC1[rownames(AUC1) != IDname, ]
      
      auc1=roc(delta~Host1, data=AUC1)$auc;
      
      AUC_data=data.frame(ID=rownames(meta.table),cutoff=rep(tt,length(Host1)), Host1=as.numeric(Host1), 
                          auc1=rep(auc1,length(Host1)))
      
      return(AUC_data)
      
    }, cox.results=cox.results, RA.table=RA.table, meta.table=meta.table, IDname=rownames(meta.table)[vv])
    
    AUC=do.call(rbind, AUC)
    
    colnames(AUC)=c("ID","cutoff", "Host1", "auc1",)
    AUC=data.frame(AUC)
    
    AUC.CV[[vv]]=AUC }
  
  return(AUC.CV)
}


for(j in unique(meta$Sample_Type_Involved)) {
  
  tem=meta[meta$Sample_Type_Involved==j,]
  RA.tem=RNA[,rownames(tem)]
  
  load(paste(j, "Host_Score_CV.RData", sep="_"))
  
  
  meta.table=tem;
  obstime <-as.numeric(unlist(meta.table[,"New_time_followup.or.death"]))
  delta <- as.numeric(3-unlist(lapply(strsplit(unlist(meta.table[,"Progression_Lab_Inv"]),'\\.'),length)))
  meta.table=data.frame(meta.table,obstime,delta)
  
  
  Scores=NULL
  for(vv in 1:length(MRS_all)) {
    
    print(vv)
    
    AUC=MRS_all[[vv]]
    
    Host1_cutoff=AUC$cutoff[which.max(AUC$auc1)];
    
    Scores=rbind(Scores, c(AUC$Host1[AUC$cutoff==Host1_cutoff & AUC$ID==names(MRS_all)[vv]], Host1_cutoff))
    
  }
  
  colnames(Scores)=c("Host1","Host1_cutoff")
  Scores=data.frame(meta.table, Scores)
  
  write.csv(Scores, file=paste(j, "Host_score.csv", sep="_"))
  
}



## Individual features


Contam=read.csv("Contaminant.list.csv",row.names = 1)
rownames(Contam)=Contam$asv
#Genus.table.com

taxa.name=tax_table(Genus.table.com)[,2:7]
taxa.name=as.character(apply(taxa.name, 1, function(x) {xx=as.character(unlist(x));
return(paste(xx[1:6],collapse = ";")) }))
tax_table(Genus.table.com)[,7]=taxa.name

for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores=read.csv(paste(j, "MRS1.csv", sep="_"), row.names = 1)
  load(paste(j, "Micro_Survival_CV.RData", sep="_"))
  
  features_all=features_NC=NULL
  
  for(vv in rownames(Scores)) {
    
    cutoffs=Scores[vv,c("MRS_all_cutoff","MRS_NC_cutoff")]
    
    cox_res=Survival_all[[vv]]
    
    features_all=c(features_all, cox_res$OTU[cox_res$prank<=as.numeric(cutoffs[1])])
    
    features_NC=c(features_NC, cox_res$OTU[cox_res$prank<=as.numeric(cutoffs[2]) & !(cox_res$contaminant)]) }
  
  
  features_all=table(features_all)
  features_all=data.frame(ASV=names(features_all), proportion=as.numeric(features_all)/nrow(Scores)*100)
  features_all=features_all[order(features_all$proportion, decreasing = TRUE),]
  
  features_all=data.frame(features_all,lineage=tax_table(Genus.table.com)[features_all$ASV,7],
                          potential.contaminant=Contam[features_all$ASV,]$potential.contaminant)
  
  
  
  features_NC=table(features_NC)
  features_NC=data.frame(ASV=names(features_NC), proportion=as.numeric(features_NC)/nrow(Scores)*100)
  features_NC=features_NC[order(features_NC$proportion, decreasing = TRUE),]
  
  features_NC=data.frame(features_NC,lineage=tax_table(Genus.table.com)[features_NC$ASV,7])
  
  write.csv(features_all, file =paste(j, "features_all.csv", sep="_") )
  
  write.csv(features_NC, file =paste(j, "features_NC.csv", sep="_") )
  
}


for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores=read.csv(paste(j, "Host_score.csv", sep="_"), row.names = 1)
  load(paste(j, "Host_Survival_CV.RData", sep="_"))
  
  features=NULL
  
  for(vv in rownames(Scores)) {
    
    cutoffs=Scores[vv,c("Host2_cutoff")]
    
    cox_res=Survival_all[[vv]]
    
    features=c(features, cox_res$gene[cox_res$prank<=as.numeric(cutoffs)]) }
  
  
  features=table(features)
  features=data.frame(Gene=names(features), proportion=as.numeric(features)/nrow(Scores)*100)
  features=features[order(features$proportion, decreasing = TRUE),]
  
  write.csv(features, file =paste(j, "features_Host.csv", sep="_") )
  
}

## AUC

for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores1=read.csv(paste(j, "MRS1.csv", sep="_"))
  Scores=read.csv(paste(j, "Host_score.csv", sep="_"))
  # identical(rownames(Scores1),rownames(Scores))
  
  Scores=data.frame(Scores, Microbiome_all=Scores1$MRS_all,Host=Scores$Host2,
                    Microbiome_NC=Scores1$MRS_NC)
  
  Scores[,c("Microbiome_all","Microbiome_NC", "Host")]=apply(Scores[,c("Microbiome_all","Microbiome_NC", "Host")], 2, scale)
  
  AA=glm(delta~Microbiome_all+Host, data=Scores, family=binomial(link = "logit"))
  both_all=predict(AA, data=Scores)
  
  AA=glm(delta~Microbiome_NC+Host, data=Scores, family=binomial(link = "logit"))
  both_NC=predict(AA, data=Scores)
  
  Scores=data.frame(Scores,both_all, both_NC)
  
  AUC=NULL
  for(tt in c("Microbiome_all","Microbiome_NC", "Host", "both_all","both_NC")) {
    
    aa=ci(roc(Scores$delta~Scores[,tt]))
    AUC=rbind(AUC,aa[c(2,1,3)]) }
  
  colnames(AUC)=c("AUC","LowerCI", "UpperCI")
  
  AUC=data.frame(Score=c("Microbiome_all","Microbiome_NC", "Host", "Microbiome_all+Host","Microbiome_NC+Host"), AUC)
  AUC[,2:4]=apply(AUC[,2:4],2, function(x) round(as.numeric(x),3))
  
  
  
  res1=kbl(AUC,
           row.names =FALSE, caption =paste(j,"Predictive performance of MRS" ,sep=": ")) %>%
    kable_classic(full_width = T, html_font = "Cambria")
  
  print(res1)
  
  
  cat("\n")
  
  
  
  AUC$Score=factor(AUC$Score, levels = AUC$Score[5:1])
  
  p <- ggplot(AUC, aes(y = Score, x = AUC)) +
    geom_point(shape = 18, size = 5, position=position_dodge(0.5)) +  
    geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.25,position=position_dodge(0.5)) + 
    geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", cex = 1, alpha = 0.5) +
    
    xlab("AUC (95% CI)") +
    ylab(" ") + 
    theme_bw() +scale_color_discrete(name="")+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x.bottom = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))+ggtitle(j)
  
  
  
  
  
  ggsave(paste(j,"MRS_Forest_plot.pdf",sep="_"),p, width = 6,height = 4)
  
  print(p)
  
}

## Scatterplots

for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores1=read.csv(paste(j, "MRS1.csv", sep="_"))
  Scores=read.csv(paste(j, "Host_score.csv", sep="_"))
  identical(rownames(Scores1),rownames(Scores))
  
  Scores=data.frame(Scores, Microbiome_all=Scores1$MRS_all,Host=Scores$Host2,
                    Microbiome_NC=Scores1$MRS_NC)
  
  Scores[,c("Microbiome_all","Microbiome_NC", "Host")]=apply( Scores[,c("Microbiome_all","Microbiome_NC", "Host")], 2, scale)
  
  figdata= melt(Scores, measure.vars = c("Microbiome_all","Microbiome_NC", "Host"))
  
  p <- ggboxplot(figdata, x = "delta", y = "value", ylab = "Risk score", xlab="",
                 color = "delta", palette = c("blue", "red"), facet.by = "variable",
                 add = "jitter")+stat_compare_means(method = "t.test")+ggtitle(j)+
    facet_wrap(~variable, nrow=1, scales = "free")
  
  
  print(p)
  
}


for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores1=read.csv(paste(j, "MRS1.csv", sep="_"))
  Scores=read.csv(paste(j, "Host_score.csv", sep="_"))
  identical(rownames(Scores1),rownames(Scores))
  
  Scores=data.frame(Scores, Microbiome_all=-Scores1$MRS_all,Host=Scores$Host2,
                    Microbiome_NC=-Scores1$MRS_NC)
  
  
  Scores[,c("Microbiome_all","Microbiome_NC", "Host")]=apply( Scores[,c("Microbiome_all","Microbiome_NC", "Host")], 2, scale)
  Scores$delta=factor(Scores$delta)
  
  p1 <- ggplot(Scores, aes(x=Microbiome_all, y=Host, colour=delta)) +
    geom_point()+theme_bw()+color_palette(palette = c("blue", "red"))+
    geom_vline(xintercept = mean(Scores$Microbiome_all), color = "grey", linetype = "dashed", cex = 1, alpha = 0.5) +
    geom_hline(yintercept = mean(Scores$Host), color = "grey", linetype = "dashed", cex = 1, alpha = 0.5) 
  
  p2 <- ggplot(Scores, aes(x=Microbiome_NC, y=Host, colour=delta)) +
    geom_point()+theme_bw()+color_palette(palette = c("blue", "red"))+
    geom_vline(xintercept = mean(Scores$Microbiome_NC), color = "grey", linetype = "dashed", cex = 1, alpha = 0.5) +
    geom_hline(yintercept = mean(Scores$Host), color = "grey", linetype = "dashed", cex = 1, alpha = 0.5)    
  
  
  p <- ggarrange(p1,p2, common.legend = TRUE)
  
  
  print(p)
  
}

## Kaplan-Meier plots
### Cutpoint: Mean

for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores1=read.csv(paste(j, "MRS1.csv", sep="_"))
  Scores=read.csv(paste(j, "Host_score.csv", sep="_"))
  identical(rownames(Scores1),rownames(Scores))
  
  Scores=data.frame(Scores, Microbiome_all=Scores1$MRS_all,Host=Scores$Host2,
                    Microbiome_NC=Scores1$MRS_NC)
  
  AA=glm(delta~Microbiome_all+Host, data=Scores, family=binomial(link = "logit"))
  both_all=predict(AA, data=Scores)
  
  AA=glm(delta~Microbiome_NC+Host, data=Scores, family=binomial(link = "logit"))
  both_NC=predict(AA, data=Scores)
  
  Scores=data.frame(Scores,both_all, both_NC)
  
  # Scores[,c("Microbiome_all","Microbiome_NC", "Host", "both_all","both_NC")]=apply(Scores[,c("Microbiome_all","Microbiome_NC", "Host", "both_all","both_NC")], 2, scale)
  
  for(tt in c("Microbiome_all","Microbiome_NC", "Host", "both_all","both_NC") ) {
    
    Scores$tem=Scores[,tt]
    
    Scores$r=rep("Low risk", nrow(Scores))
    Scores$r[Scores$tem>=mean(Scores$tem)]="High risk"
    
    if(tt=="both_all") title="Microbiome_all+Host" else if(tt=="both_NC") title="Microbiome_NC+Host" else title=tt
    
    aa=ggsurvplot(
      fit = survfit(Surv(obstime, delta) ~ r,data=Scores),
      xlab = "Days", 
      ylab = "Overall survival probability",
      legend.title = "",
      legend.labs = c("High risk", "Low risk"),
      palette = c("red","blue"),
      risk.table = TRUE,risk.table.title=NA, tables.height=0.4,
      conf.int = TRUE, pval = TRUE)+ggtitle(paste(j, title, sep="; "))
    
    print(aa)
    
  }
  
  tt="both_NC"
  # for(tt in c("both_NC") ) {
  
  Scores$tem=Scores[,tt]
  
  Scores$r=rep("Low risk", nrow(Scores))
  Scores$r[Scores$tem>=mean(Scores$tem)]="High risk"
  
  if(tt=="both_all") title="Microbiome_all+Host" else if(tt=="both_NC") title="Microbiome_NC+Host" else title=tt
  
  aa=ggsurvplot(
    fit = survfit(Surv(obstime, delta) ~ r,data=Scores),
    xlab = "Days", 
    ylab = "Overall survival probability",
    legend.title = "",
    legend.labs = c("High risk", "Low risk"),
    palette = c("red","blue"),
    risk.table = TRUE,risk.table.title=NA, tables.height=0.4,
    conf.int = TRUE, pval = TRUE)+ggtitle(title)
  
  ## 
  
  pdf(paste(j,"KM_plot.pdf",sep="_"), width = 5, height = 6)
  print(aa, newpage = FALSE)
  dev.off()
  
}

##Another strategy, for integrated score (Microbiome_all/Microbiome_NC+Host), the sample whose Microbiome_all/Microbiome_NC is greater than the mean and Host is greater than the mean, is defined as "High risk" (The right top panel in the scatterplot), otherwise "Low risk".

for(j in c("Lung.Tissue.In", "Lung.Tissue.UnIn")) {
  
  Scores1=read.csv(paste(j, "MRS1.csv", sep="_"))
  Scores=read.csv(paste(j, "Host_score.csv", sep="_"))
  identical(rownames(Scores1),rownames(Scores))
  
  Scores=data.frame(Scores, Microbiome_all=Scores1$MRS_all,Host=Scores$Host2,
                    Microbiome_NC=Scores1$MRS_NC)
  
  Scores$r=rep("Low risk", nrow(Scores))
  Scores$r[Scores$Microbiome_NC>=mean(Scores$Microbiome_NC) & Scores$Host>=mean(Scores$Host)]="High risk"
  
  aa=ggsurvplot(
    fit = survfit(Surv(obstime, delta) ~ r,data=Scores),
    xlab = "Days", 
    ylab = "Overall survival probability",
    legend.title = "",
    legend.labs = c("High risk", "Low risk"),
    palette = c("red","blue"),
    risk.table = TRUE,risk.table.title=NA, tables.height=0.4,
    conf.int = TRUE, pval = TRUE)+ggtitle("Microbiome_NC+Host")
  
  print(aa)
  
  pdf(paste(j,"KM_plot2.pdf",sep="_"), width = 5, height = 6)
  print(aa, newpage = FALSE)
  dev.off()
  
  
}


