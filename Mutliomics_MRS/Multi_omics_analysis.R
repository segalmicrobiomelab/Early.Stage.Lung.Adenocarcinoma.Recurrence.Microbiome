
filedir="C:/Users/WANGC22/Downloads/Segal/New.Early.Stage.Lung.Cancer.Data"

setwd(filedir)

pac=c("qiime2R","magrittr","tibble","tidyr","vegan","kableExtra", "ComplexHeatmap",
      'survival','OMiSA','massMap',"ape","PMA", "plyr","survminer","pROC",
      "DESeq2","phyloseq","ANCOMBC","table1","ggcorrplot","reshape2","GGally","ggpubr")

sapply(pac, require, character.only=T)



load("phy.RData")

phy = prune_samples(sample_data(phy)$Sample_Type_Involved %in% c("Lung.Tissue.In", "Lung.Tissue.UnIn"), phy)

# Remove taxa with 0 abundance
phy = subset_taxa(phy, rowSums(otu_table(phy)) != 0)
phy = subset_taxa(phy, Kingdom=="k__Bacteria")


RNA=read.delim(paste("RNA.host.Transcriptome","Copy of KO_Pathway_broken_Adenocarcinoma_Samples_sum_path03_JT.lns.txt", sep="/"), check.names = FALSE)

pathnames=RNA$path_03
names(pathnames)=paste0("f", 1:length(pathnames))

AA=sample_data(phy)
AA.rownames=rownames(AA)
# AA.rownames[!grepl("NYU",AA.rownames)]=gsub("\\.Lung", "", AA.rownames[!grepl("NYU",AA.rownames)])
# AA.rownames[AA.rownames=="NYU1923..Lung"]="NYU1923.Lung"
#identical(both, AA.rownames[AA.rownames %in% both])

both=intersect(AA.rownames,colnames(RNA))

RNA=RNA[,both]
rownames(RNA)=names(pathnames)
phy <- prune_samples(AA.rownames %in% both, phy)
RNA=RNA[apply(RNA, 1, function(x) sum(x>=200)>=3),]

## Analyses for the components integrated by microbiome and transcriptome
### Correlation test of the components integrated by microbiome and transcriptome


## filtering criteria
Genus.rel.table = transformSampleCounts(phy, function(x) x/sum(x))

Genus.rel.table.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.002), A = 0.01 * nsamples(Genus.rel.table))
Genus.table.com = prune_taxa(Genus.rel.table.wh1, Genus.rel.table)

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

biomarker=vector("list", length = 4)
names(biomarker)=unlist(unique(meta[,"Progression_Lab_Inv"]))
for(j in unlist(unique(meta[,"Progression_Lab_Inv"]))) {

   print(j)

    tem.microbiome=microbiome.data[meta$Progression_Lab_Inv==j,];
    taxa.0=which(colSums(tem.microbiome)==0)
    if(length(taxa.0)>0) tem.microbiome=tem.microbiome[,-c(taxa.0)]

    tem.transcriptome=transcriptome[meta$Progression_Lab_Inv==j,]
    transcriptome.0=which(colSums(tem.transcriptome)==0)
    if(length(transcriptome.0)>0) tem.transcriptome=tem.transcriptome[,-c(transcriptome.0)]

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

# save(biomarker,file="biomarker_fourgroups_OTU.RData")
# 
# load("biomarker_fourgroups_OTU.RData")

#### feature contribution
for(j in names(biomarker)) {

transcriptome.load=biomarker[[j]][[2]][,1:3]
transcriptome.load=transcriptome.load[order(abs(as.numeric(transcriptome.load[,2])),decreasing = TRUE),]

transcriptome.load$gene=pathnames[transcriptome.load$gene]

write.csv(transcriptome.load, file=paste0(j,"trans.csv"))
}

Contam=read.csv("Contaminant.list.csv",row.names = 1)

for(j in names(biomarker)) {


  microbiome.load=biomarker[[j]][[1]][,1:3]
  microbiome.load=microbiome.load[order(abs(as.numeric(microbiome.load[,2])),decreasing = TRUE),]


tax.name=tax_table(Genus.table.com)
lineage=as.character(tax.name[microbiome.load$taxa,7])

microbiome.load=data.frame(microbiome.load,lineage)

rownames(microbiome.load)=microbiome.load$taxa

contam1=Contam[Contam$asv %in%rownames(microbiome.load), ]

microbiome.load=data.frame(microbiome.load[contam1$asv,],contaminant=contam1$potential.contaminant)

write.csv(microbiome.load, file=paste0(j,"micro.csv"))

}

### Involved samples
features=vector("list",4)
names(features)=names(biomarker)

for(j in names(biomarker)) {
  
  microbiome.load=biomarker[[j]][[1]][,1:3]
  transcriptome.load=biomarker[[j]][[2]][,1:3]
  
  
  microbiome.load[,2:3]=abs(microbiome.load[,2:3])
  micro.tem=microbiome.load[order(microbiome.load[,2], decreasing = TRUE),]
  micro1=micro.tem$taxa[micro.tem[,2]>=0.05]
  
  
  micro.tem=microbiome.load[order(microbiome.load[,3], decreasing = TRUE),]
  micro2=micro.tem$taxa[micro.tem[,3]>=0.05]
  
  Micro.im=unique(c(micro1[1:min(length(micro1),50)],micro2[1:min(length(micro2),50)]))
  
  transcriptome.load[,2:3]=abs(transcriptome.load[,2:3])
  trans.tem=transcriptome.load[order(transcriptome.load[,2], decreasing = TRUE),]
  trans1=trans.tem$gene[trans.tem[,2]>=0.05]
  
  trans.tem=transcriptome.load[order(transcriptome.load[,3], decreasing = TRUE),]
  trans2=trans.tem$gene[trans.tem[,3]>=0.05]
  
  Trans.im=unique(c(trans1[1:min(length(trans1),50)],trans2[1:min(length(trans2),50)]))
  
  
  
  # Micro.im=microbiome.load$taxa[apply(microbiome.load[,2:3],1, function(x) max(abs(as.numeric(x)))>=0.05)]
  # Trans.im=transcriptome.load$gene[apply(transcriptome.load[,2:3],1, function(x) max(abs(as.numeric(x)))>=0.05)]
  
  features[[j]]=list(Micro.im,Trans.im)
  
}



print(names(features))
#  "In.Recurrence"      "In.No.Recurrence"   "UnIn.No.Recurrence" "UnIn.Recurrence"   

Taxa.In.Re=features[[1]][[1]];gene.In.Re=features[[1]][[2]];
Taxa.In.No=features[[2]][[1]];gene.In.No=features[[2]][[2]];

Taxa.Re=setdiff(Taxa.In.Re, Taxa.In.No); Taxa.both=intersect(Taxa.In.Re, Taxa.In.No);
Taxa.No=setdiff(Taxa.In.No,Taxa.In.Re);
Taxa.merge=c(Taxa.No,Taxa.both,Taxa.Re)

gene.Re=setdiff(gene.In.Re, gene.In.No); gene.both=intersect(gene.In.Re, gene.In.No);
gene.No=setdiff(gene.In.No,gene.In.Re);
gene.merge=c(gene.No,gene.both,gene.Re)

delta <- 3-unlist(lapply(strsplit(unlist(meta[,"Progression_Lab_Inv"]),'\\.'),length))

data.ph=data.frame(Involved=revalue(as.character(meta$Lung_Tissue_Involved), c("0"="UnInvolved", "1"="Involved")), status=revalue(as.character(delta), c("0"="No Recurrence", "1"="Recurrence")),
                   microbiome.data[,Taxa.merge], transcriptome[,gene.merge])


tem.data=data.ph[data.ph$Involved=="Involved",]

Cov.no=cor(tem.data[tem.data$status=="No Recurrence",-c(1:2)])
Cov.no=Cov.no[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

Cov.Re=cor(tem.data[tem.data$status=="Recurrence",-c(1:2)])
Cov.Re=Cov.Re[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

Cov_merge=round(rbind(Cov.no,Cov.Re),3)

p.no = cor_pmat(tem.data[tem.data$status=="No Recurrence",-c(1:2)],method ="pearson")
p.no=p.no[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

p.Re=cor_pmat(tem.data[tem.data$status=="Recurrence",-c(1:2)],method ="pearson")
p.Re=p.Re[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

p_merge=round(rbind(p.no,p.Re),3)


tax.name=tax_table(Genus.table.com)
cor.name=sapply(as.character(tax.name[Taxa.merge,7]), function(x) paste(unlist(strsplit(x, ";"))[3:6],collapse = ";"))


Cov_merge=t(Cov_merge)
colnames(Cov_merge)=rep(cor.name,2)
rownames(Cov_merge)=pathnames[rownames(Cov_merge)]

p_merge=t(p_merge)
colnames(p_merge)=rep(cor.name,2)
rownames(p_merge)=pathnames[rownames(p_merge)]

p_merge[p_merge<=0.05]="*"
p_merge[(p_merge!="*")| is.na(p_merge)]=" "


NOR=colMeans(tem.data[tem.data$status=="No Recurrence",3:(2+length(Taxa.merge))],na.rm = TRUE)
R=colMeans(tem.data[tem.data$status=="Recurrence",3:(2+length(Taxa.merge))],na.rm = TRUE)

names(NOR)=names(R)=cor.name
#  RA.log=-log10(c(NOR,R))


ha2 = HeatmapAnnotation(barplot = anno_barplot(ifelse(c(NOR,R)==0, NA, -log10(c(NOR,R))), axis = TRUE),
                        annotation_label ="-log10(Relative Abundance)",
                        annotation_name_side ="left",
                        annotation_name_rot =00)



ha = HeatmapAnnotation(df = data.frame(Sample1 = c(rep("No Recurrence", length(Taxa.merge)), 
                                                   rep("Recurrence", length(Taxa.merge))),
                                       
                                       Sample2 = rep(c(rep("Only No Recurrence",length(Taxa.No)),
                                                       rep("Both", length(Taxa.both)),
                                                       rep("Only Recurrence", length(Taxa.Re))),2)), 
                       col = list(Sample1 = c("No Recurrence" =  "blue", "Recurrence" = "red"),
                                  Sample2 = c("Only No Recurrence" =  "blue","Both" =  "green", "Only Recurrence" = "red")),
                       annotation_label =c("Sample used for correlation calculation",
                                           "Sample used for identification"),
                       show_legend =FALSE,
                       annotation_name_side ="left")


ha3 = rowAnnotation(df = data.frame(Sample2 = c(rep("Only No Recurrence", length(gene.No)),
                                                rep("Both", length(gene.both)),
                                                rep("Only Recurrence", length(gene.Re)))), 
                    col = list(Sample2 = c("Only No Recurrence" =  "blue","Both" =  "green", 
                                           "Only Recurrence" = "red")),
                    annotation_label =c("Sample used for identification"),
                    show_legend =FALSE)



lgd_boxplot1 = Legend(labels = c("No Recurrence", "Recurrence"), title = "Sample used for correlation calculation",
                      legend_gp = gpar(fill = c("blue", "red")))

lgd_boxplot2 = Legend(labels = c("Only No Recurrence","Both", "Only Recurrence"), title = "Sample used for identification",
                      legend_gp = gpar(fill = c("blue","green", "red")))


heatmap1=Heatmap(Cov_merge,cluster_rows=FALSE,cluster_columns=FALSE, col=c("darkblue", "white", "red"),
                 row_names_side ="left",
                 heatmap_legend_param = list(title = "Correlation"),top_annotation = ha,
                 bottom_annotation = ha2,right_annotation = ha3,
                 column_names_max_height = unit(16, "cm"),
                 row_names_max_width = unit(14, "cm"),
                 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   
                   grid.text(p_merge[i, j], x, y, gp = gpar(fontsize = 10)) }
                 
)

pdf(file="Involved_Heatmap_v2.pdf",width = 10,height = 8)

draw(heatmap1, 
     heatmap_legend_list = list(lgd_boxplot1,lgd_boxplot2))

dev.off()


### UnInvolved samples

print(names(features))
#  "In.Recurrence"      "In.No.Recurrence"   "UnIn.No.Recurrence" "UnIn.Recurrence"   

Taxa.In.Re=features[[4]][[1]];gene.In.Re=features[[4]][[2]];
Taxa.In.No=features[[3]][[1]];gene.In.No=features[[3]][[2]];

Taxa.Re=setdiff(Taxa.In.Re, Taxa.In.No); Taxa.both=intersect(Taxa.In.Re, Taxa.In.No);
Taxa.No=setdiff(Taxa.In.No,Taxa.In.Re);
Taxa.merge=c(Taxa.No,Taxa.both,Taxa.Re)

gene.Re=setdiff(gene.In.Re, gene.In.No); gene.both=intersect(gene.In.Re, gene.In.No);
gene.No=setdiff(gene.In.No,gene.In.Re);
gene.merge=c(gene.No,gene.both,gene.Re)

delta <- 3-unlist(lapply(strsplit(unlist(meta[,"Progression_Lab_Inv"]),'\\.'),length))

data.ph=data.frame(Involved=revalue(as.character(meta$Lung_Tissue_Involved), c("0"="UnInvolved", "1"="Involved")), status=revalue(as.character(delta), c("0"="No Recurrence", "1"="Recurrence")),
                   microbiome.data[,Taxa.merge], transcriptome[,gene.merge])


tem.data=data.ph[data.ph$Involved!="Involved",]

Cov.no=cor(tem.data[tem.data$status=="No Recurrence",-c(1:2)])
Cov.no=Cov.no[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

Cov.Re=cor(tem.data[tem.data$status=="Recurrence",-c(1:2)])
Cov.Re=Cov.Re[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

Cov_merge=round(rbind(Cov.no,Cov.Re),3)


p.no = cor_pmat(tem.data[tem.data$status=="No Recurrence",-c(1:2)],method ="pearson")
p.no=p.no[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

p.Re=cor_pmat(tem.data[tem.data$status=="Recurrence",-c(1:2)],method ="pearson")
p.Re=p.Re[1:length(Taxa.merge),-c(1:length(Taxa.merge))]

p_merge=round(rbind(p.no,p.Re),3)


tax.name=tax_table(Genus.table.com)
cor.name=sapply(as.character(tax.name[Taxa.merge,7]), function(x) paste(unlist(strsplit(x, ";"))[3:6],collapse = ";"))


Cov_merge=t(Cov_merge)
colnames(Cov_merge)=rep(cor.name,2)
rownames(Cov_merge)=pathnames[rownames(Cov_merge)]

p_merge=t(p_merge)
colnames(p_merge)=rep(cor.name,2)
rownames(p_merge)=pathnames[rownames(p_merge)]

p_merge[p_merge<=0.05]="*"
p_merge[(p_merge!="*")| is.na(p_merge)]=" "


NOR=colMeans(tem.data[tem.data$status=="No Recurrence",3:(2+length(Taxa.merge))],na.rm = TRUE)
R=colMeans(tem.data[tem.data$status=="Recurrence",3:(2+length(Taxa.merge))],na.rm = TRUE)

ha2 = HeatmapAnnotation(barplot = anno_barplot(ifelse(c(NOR,R)==0, NA, -log10(c(NOR,R))), axis = TRUE),
                        annotation_label ="-log10(Abundance)",
                        annotation_name_side ="left",
                        annotation_name_rot =00)



ha = HeatmapAnnotation(df = data.frame(Sample1 = c(rep("No Recurrence", length(Taxa.merge)), 
                                                   rep("Recurrence", length(Taxa.merge))),
                                       
                                       Sample2 = rep(c(rep("Only No Recurrence",length(Taxa.No)),
                                                       rep("Both", length(Taxa.both)),
                                                       rep("Only Recurrence", length(Taxa.Re))),2)), 
                       col = list(Sample1 = c("No Recurrence" =  "blue", "Recurrence" = "red"),
                                  Sample2 = c("Only No Recurrence" =  "blue","Both" =  "green", "Only Recurrence" = "red")),
                       annotation_label =c("Sample used for correlation calculation",
                                           "Sample used for identification"),
                       show_legend =FALSE,
                       annotation_name_side ="left")


ha3 = rowAnnotation(df = data.frame(Sample2 = c(rep("Only No Recurrence", length(gene.No)),
                                                rep("Both", length(gene.both)),
                                                rep("Only Recurrence", length(gene.Re)))), 
                    col = list(Sample2 = c("Only No Recurrence" =  "blue","Both" =  "green", 
                                           "Only Recurrence" = "red")),
                    annotation_label =c("Sample used for identification"),
                    show_legend =FALSE)



lgd_boxplot1 = Legend(labels = c("No Recurrence", "Recurrence"), title = "Sample used for correlation calculation",
                      legend_gp = gpar(fill = c("blue", "red")))

lgd_boxplot2 = Legend(labels = c("Only No Recurrence","Both", "Only Recurrence"), title = "Sample used for identification",
                      legend_gp = gpar(fill = c("blue","green", "red")))


heatmap1=Heatmap(Cov_merge,cluster_rows=FALSE,cluster_columns=FALSE, col=c("darkblue", "white", "red"),
                 row_names_side ="left",
                 heatmap_legend_param = list(title = "Correlation"),top_annotation = ha,
                 bottom_annotation = ha2,right_annotation = ha3,
                 column_names_max_height = unit(18, "cm"),
                 row_names_max_width = unit(16, "cm"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   
                   grid.text(p_merge[i, j], x, y, gp = gpar(fontsize = 10)) }
)




pdf(file="UnInvolved_Heatmap_v2.pdf",width = 24,height = 22)

draw(heatmap1, 
     heatmap_legend_list = list(lgd_boxplot1,lgd_boxplot2))

dev.off()

