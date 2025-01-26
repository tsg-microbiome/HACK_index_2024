dir_name <- "G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\"
load(paste0(dir_name,"c56_part2AllControlStudyRun.RData"))
load(paste0(dir_name,"controlMetadata.RData"))
load(paste0(dir_name,"controlSpProfile.RData"))
load(paste0(dir_name,"controlSpProfileAll.RData"))

source(paste0(dir_name,"code_library.r"))
AllStudiesDiscoveryCohort <- read.table(paste0(dir_name,"AllStudiesDiscoveryCohort.txt"),sep="\t",row.names=1,header=TRUE)
#all_dataset_prevalent_species: 70% prevalent species of all studies
#controlStudyInfo: Information of all studies included for the CoreKeyStone identification
#r2Df: Study Ranked Envfit R2-Squared Values of Species in Each Study
#prDf: Envfit P Values of Species in Each Study
#prevalentDf: Prevalence of Different Species in Each Study
#controlSpProfile: Species Abundance
#controlMetadata: Metadata

load(paste0(dir_name,"CoreInfluencers_3Dfs.RData"))

prevalentDf <- prevalentDf_new
prDf <- prDf_new
r2Df <- r2Df_new

controlSpProfile$study_name <- controlMetadata$study_name
prevalentDf <- prevalentDf[,colnames(r2Df)]

temp <- controlSpProfileAll[rownames(controlSpProfile),rownames(r2Df)]
temp$study_name <- controlSpProfile$study_name
controlSpProfile <- temp

CoreKeyStoneDf <- apply(r2Df,2,function(x)(ifelse(x>=0.70,1,0))) * apply(prevalentDf,2,function(x)(ifelse(x>=0.65,1,0)))
MajorCoreKeyStone <- rownames(CoreKeyStoneDf[rowSums(CoreKeyStoneDf)>=12,])

CohortLifeStyle <- AllStudiesDiscoveryCohort[colnames(prevalentDf),7]
names(CohortLifeStyle) <- colnames(prevalentDf)

CohortSeqType <- AllStudiesDiscoveryCohort[colnames(prevalentDf),6]
names(CohortSeqType) <- colnames(prevalentDf)

thresholds <- seq(0,72,6)

RowSumCoreKeyStone <- rowSums(CoreKeyStoneDf)
RowSumCoreKeyStone <- as.data.frame(RowSumCoreKeyStone[order(RowSumCoreKeyStone)])
RowSumCoreKeyStone$index <- 1:nrow(RowSumCoreKeyStone)
colnames(RowSumCoreKeyStone) <- c("Total","Index")
temp_df <- RowSumCoreKeyStone[RowSumCoreKeyStone$Total>=12,]
ggplot(temp_df[temp_df$Total>=12,],aes(x=Index,y=Total))+geom_point(color=ifelse(temp_df$Total >= 12,"deepskyblue","grey"))+geom_text_repel(label=ifelse(temp_df$Total>=12,rownames(temp_df),""),max.overlaps=20,size=3.2,box.padding=0.1,min.segment.length=1)+theme_bw()

DetectionProfile <- NULL

for(i in 2:length(thresholds))
{
	range0 <- thresholds[i-1]
	range1 <- thresholds[i]
	detected <- length(which((RowSumCoreKeyStone$Total >= range0)&(RowSumCoreKeyStone$Total < range1)))
	DetectionProfile <- c(DetectionProfile,detected)
}

names(DetectionProfile) <- thresholds[2:length(thresholds)]

pcoCoreKeyStoneDf <- dudi.pco(vegdist(t(CoreKeyStoneDf[MajorCoreKeyStone,]),method="jaccard"),scannf=FALSE,nf=20)
pcoSpeciesCoreKeyStoneDf <- dudi.pco(vegdist(CoreKeyStoneDf[MajorCoreKeyStone,],method="jaccard"),scannf=FALSE,nf=20)

seqtype_PERMANOVA_CoreKeyStone <- adonis2(vegdist(pcoCoreKeyStoneDf$li,method="euclidean")~as.factor(CohortSeqType[rownames(pcoCoreKeyStoneDf$li)]))

lifestyle_PERMANOVA_CoreKeyStone <- adonis2(vegdist(pcoCoreKeyStoneDf$li,method="euclidean")~as.factor(CohortLifeStyle[rownames(pcoCoreKeyStoneDf$li)]))


pcoSchematic <- data.frame(A1=rnorm(50,0,1),A2=rnorm(50,0,1))
pcoSchematic$Color <- paste0("grey",10*as.numeric(cut(pcoSchematic$A1,breaks=quantile(pcoSchematic$A1,prob=seq(0,1,0.1)),include.lowest=TRUE)))
pcoSchematic$RandomColor <- sample(pcoSchematic$Color,50,replace=FALSE)

ggplot(pcoSchematic,aes(x=A1,y=A2))+geom_point(size=5,color=pcoSchematic$Color)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggplot(pcoSchematic,aes(x=A1,y=A2))+geom_point(size=5,color=pcoSchematic$RandomColor)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggplot(pcoSchematic,aes(x=A1,y=A2))+geom_point(size=5,color="grey50")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

mat_CoreKeyStone <- CoreKeyStoneDf[intersect(rownames(RowSumCoreKeyStone),MajorCoreKeyStone),]
#heatmap.2(mat_CoreKeyStone[,names(CohortLifeStyle)],density="none",trace="none",Rowv=FALSE,lhei=c(0.1,5),lwid=c(0.1,5),margins=c(9,30),srtRow=0.1,sepcolor="grey",sepwidth=c(0.1,0.1),rowsep=c(0:nrow(mat_CoreKeyStone)),colsep=c(0:ncol(mat_CoreKeyStone)),ColSideColor=c("Red","Blue","Green")[as.factor(CohortLifeStyle)],cexRow=1.1,col=c("white","deepskyblue"))

heatmap.2(apply(mat_KeyStone[,names(CohortLifeStyle)], 2, function(x) ifelse(x >= 0.65, 1, 0)),
          density="none",
          trace="none",
          Rowv=FALSE,
          lhei=c(0.1, 1.2),
          lwid=c(0.1, 1),
          margins=c(25,30),
          srtRow=0.1,
          sepcolor="darkgray",
          sepwidth=c(0.05, 0.05),
          rowsep=c(0:nrow(mat_KeyStone)),
          colsep=c(0:ncol(mat_KeyStone)),
          ColSideColor=c("Red", "Blue", "Green")[as.factor(CohortLifeStyle)],
          cexRow= 2,
          cexCol = 2,
          col=c("white", "deepskyblue"),
          cellnote=apply(as.matrix(prDf[intersect(rownames(RowSumCoreKeyStone), MajorCoreKeyStone), names(CohortLifeStyle)]), 2, function(x) ifelse(x <= 0.10, "*", "")),
          notecol="black", notecex=2)
mat_CoreKeyStone_Detection <- data.frame(A1=apply(mat_CoreKeyStone[,names(which(CohortLifeStyle=="IndustrializedUrban"))],1,function(x)(length(x[x==0]))),A2=apply(mat_CoreKeyStone[,names(which(CohortLifeStyle=="IndustrializedUrban"))],1,function(x)(length(x[x!=0]))),A3=apply(mat_CoreKeyStone[,names(which(CohortLifeStyle!="IndustrializedUrban"))],1,function(x)(length(x[x==0]))),A4=apply(mat_CoreKeyStone[,names(which(CohortLifeStyle!="IndustrializedUrban"))],1,function(x)(length(x[x!=0]))))

study_aggregated_mat_core <- t(aggregate(t(CoreKeyStoneDf[MajorCoreKeyStone,names(CohortLifeStyle)]),by=list(factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal"))),FUN=sum)[,-1]/as.numeric(table(CohortLifeStyle)[c(1,3,2)]))
colnames(study_aggregated_mat_core) <- rev(c("RuralTribal","UrbanRuralMixed","IndustrializedUrban"))

study_aggregated_mat_core_all <- t(aggregate(t(CoreKeyStoneDf[,names(CohortLifeStyle)]),by=list(factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal"))),FUN=sum)[,-1]/as.numeric(table(CohortLifeStyle)[c(1,3,2)]))
colnames(study_aggregated_mat_core) <- rev(c("RuralTribal","UrbanRuralMixed","IndustrializedUrban"))

#heatmap.2(study_aggregated_mat_core[rownames(mat_CoreKeyStone),],density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Greens"),lhei=c(0.1,1),margins=c(20,40),srtRow=0.5,sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:ncol(study_aggregated_mat_core),rowsep=0:nrow(study_aggregated_mat_core))

pcoPrevalenceMajorKeyStones <- dudi.pco(vegdist(t(prevalentDf[MajorCoreKeyStone,]),method="euclidean"),scannf=FALSE,nf=3)
s.class(pcoPrevalenceMajorKeyStones$li[names(CohortLifeStyle),],factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")),col=c("Red","Blue","Green"),plabels.col="black",plabels.cex=1.2)

pcoKeyStoneMajorKeyStones <- dudi.pco(vegdist(t(r2Df[MajorCoreKeyStone,]),method="euclidean"),scannf=FALSE,nf=3)
s.class(pcoKeyStoneMajorKeyStones$li[names(CohortLifeStyle),],factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")),col=c("Red","Blue","Green"),plabels.col="black",plabels.cex=1.2)

#adonis2(vegdist(pcoKeyStoneMajorKeyStones$li,method="euclidean")~factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")))
#adonis2(vegdist(t(r2Df[MajorCoreKeyStone,]),method="euclidean")~factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")))

adonis2(vegdist(t(prevalentDf[MajorCoreKeyStone,]),method="euclidean")~factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")))
adonis2(vegdist(t(r2Df[MajorCoreKeyStone,]),method="euclidean")~factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")))

study_aggregated_r_squared <- t(aggregate(t(r2Df[intersect(rownames(RowSumCoreKeyStone),MajorCoreKeyStone),names(CohortLifeStyle)]),by=list(factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal"))),FUN=median)[,-1])

study_aggregated_r_squared_all <- t(aggregate(t(r2Df[,names(CohortLifeStyle)]),by=list(factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal"))),FUN=median)[,-1])

mat_KeyStone <- as.matrix(r2Df[intersect(rownames(RowSumCoreKeyStone),MajorCoreKeyStone),])
#heatmap.2(mat_KeyStone[,names(CohortLifeStyle)],density="none",trace="none",Rowv=FALSE,lhei=c(0.1,5),lwid=c(0.1,5),margins=c(9,30),srtRow=0.1,sepcolor="grey",sepwidth=c(0.1,0.1),rowsep=c(0:nrow(mat_KeyStone)),colsep=c(0:ncol(mat_KeyStone)),ColSideColor=c("Red","Blue","Green")[as.factor(CohortLifeStyle)],cexRow=1.1,col=c("white","deepskyblue"))
#heatmap.2(apply(mat_KeyStone[,names(CohortLifeStyle)],2,function(x)(ifelse(x>=0.65,1,0))),density="none",trace="none",Rowv=FALSE,lhei=c(0.1,5),lwid=c(0.1,5),margins=c(9,30),srtRow=0.1,sepcolor="grey",sepwidth=c(0.1,0.1),rowsep=c(0:nrow(mat_KeyStone)),colsep=c(0:ncol(mat_KeyStone)),ColSideColor=c("Red","Blue","Green")[as.factor(CohortLifeStyle)],cexRow=1.1,col=c("white","deepskyblue"),cellnote=apply(as.matrix(prDf[intersect(rownames(RowSumCoreKeyStone),MajorCoreKeyStone),names(CohortLifeStyle)]),2,function(x)(ifelse(x<=0.10,"*",""))),notecol="black")

mat_Prevalence <- as.matrix(prevalentDf[intersect(rownames(RowSumCoreKeyStone),MajorCoreKeyStone),])
#heatmap.2(apply(mat_Prevalence[,names(CohortLifeStyle)],2,function(x)(ifelse(x>=0.70,1,0))),density="none",trace="none",Rowv=FALSE,lhei=c(0.1,5),lwid=c(0.1,5),margins=c(9,30),srtRow=0.1,sepcolor="grey",sepwidth=c(0.1,0.1),rowsep=c(0:nrow(mat_Prevalence)),colsep=c(0:ncol(mat_Prevalence)),ColSideColor=c("Red","Blue","Green")[as.factor(CohortLifeStyle)],cexRow=1.1,col=c("white","deepskyblue"))

r2_vec <- as.vector(apply(mat_KeyStone[,names(CohortLifeStyle)],2,function(x)(ifelse(x>=0.70,1,0))))
pr_vec <- as.vector(apply(as.matrix(prDf[intersect(rownames(RowSumCoreKeyStone),MajorCoreKeyStone),names(CohortLifeStyle)]),2,function(x)(ifelse(x<=0.05,1,0))))
(length(which((pr_vec==1)&(r2_vec==1))) + length(which((pr_vec==0)&(r2_vec==0))))/3312


SelectSpecies <- rownames(RowSumCoreKeyStone[RowSumCoreKeyStone[,1]>12,])
MeanSpAbundances <- aggregate(controlSpProfile[,rownames(CoreKeyStoneDf)],by=list(controlSpProfile$study_name),FUN=mean)[,-1]
rownames(MeanSpAbundances) <- aggregate(controlSpProfile[,rownames(CoreKeyStoneDf)],by=list(controlSpProfile$study_name),FUN=mean)[,1]
MeanSpAbundances <- t(MeanSpAbundances)
dfCorrAbundanceInfluence <- corr.test(MeanSpAbundances[SelectSpecies,],r2Df[SelectSpecies,colnames(MeanSpAbundances)])
dfCorrAbundancePrevalence <- corr.test(MeanSpAbundances[SelectSpecies,],prevalentDf[SelectSpecies,colnames(MeanSpAbundances)])

dfAllAbundanceAssociationsR <- data.frame(influenceR = diag(dfCorrAbundanceInfluence$r),prevalenceR=diag(dfCorrAbundancePrevalence$r))
dfAllAbundanceAssociationsP <- data.frame(influenceR = diag(dfCorrAbundanceInfluence$p),prevalenceR=diag(dfCorrAbundancePrevalence$p))

#heatmap.2(t(dfAllAbundanceAssociationsR),density="none",trace="none",Rowv=FALSE,Colv=FALSE,col=brewer.pal(8,"PuOr"),margins=c(20,20))
heatmap.2(t(dfAllAbundanceAssociationsR),density="none",trace="none",Rowv=FALSE,Colv=FALSE,col=brewer.pal(8,"PuOr"),margins=c(70,20),lhei=c(0.1,1),lwid=c(0.1,1),cellnote=apply(t(dfAllAbundanceAssociationsP),2,function(x)(ifelse(x<=0.05,"*",""))),notecol="black",notecex=2, cexRow=1.2, cexCol=1.1)

AllSelectSpeciesAbundanceStudy <- data.frame(SummedAbundance = rowSums(controlSpProfileAll[rownames(controlSpProfile),rownames(CoreKeyStoneDf)]),Study=controlSpProfile$study_name)
ggplot(AllSelectSpeciesAbundanceStudy,aes(x=SummedAbundance,y=Study,group=Study))+geom_boxplot()+theme(axis.text.y=element_text(size=5))


keystoneInfluence <- function(species, inputData) {
  tryCatch({
	set.seed(100)
    colIndex= which(colnames(inputData)==species)
    newData= inputData[,-colIndex]
    newData <- newData[which(rowSums(newData)!= 0), ]
    #print(dim(newData))
    newData <- newData / rowSums(newData)
	#print(unique(rowSums(newData)))
    #cat("Creating distance matrix\n")
    distanceMatrix <- vegdist(newData, method = "bray")
    #print("distance matrix done")
    #cat("Creating dudi.pco\n")
    pco <- dudi.pco(distanceMatrix, scannf = FALSE)
    #print("pco is created")
    pcoPointsDf <- pco$li
    #cat("Generating model\n")
    model <- envfit(pcoPointsDf ~ inputData[rownames(newData), species])
	r_value <- as.numeric(model$vectors$r)
	p_value <- as.numeric(model$vectors$pvals)
	return_list = list("r"=r_value,"p"=p_value)
	return(return_list)
	})
}

compute_prevalence_single_data <- function(data,threshold){
  detection_percentage <- colSums(apply(data,2,function(x)(ifelse(x>0,1,0))))/nrow(data)
  #highly_detected <- names(which(detection_percentage>=threshold))
  return(detection_percentage)
}

CoreKeyStoneProperties <- function(data,species_list,threshold1,threshold2)
{
	prevalence <- compute_prevalence_single_data(data[,species_list],threshold)
	r2_community_influence <- rep(0,196)
	names(r2_community_influence) <- species_list
	pr_community_influence <- rep(1,196)
	names(pr_community_influence) <- species_list
	for(i in 1:length(species_list))
	{
		species_name <- species_list[i]
		print(species_name)
		community_influence <- keystoneInfluence(species_name,data)
		r2_community_influence[species_name] <- community_influence$r
		pr_community_influence[species_name] <- community_influence$p
	}
	r2_community_influence <- ifelse(is.na(r2_community_influence),0,r2_community_influence)
	r2_community_influence <- rank_scale(r2_community_influence)
	pr_community_influence <- ifelse(is.na(pr_community_influence),0,pr_community_influence)
	core_keystone <- ifelse(r2_community_influence >= 0.70,1,0) * ifelse(prevalence >= 0.70,1,0)
	names(core_keystone) <- species_list
	return_list = list("prevalence"=prevalence,"r2_community_influence"=r2_community_influence,"pr_community_influence"=pr_community_influence,"core_keystone"=core_keystone)
	return(return_list)
}

meanStudyAbundance <- aggregate(controlSpProfile[,1:196],by=list(controlSpProfile$study_name),FUN=median)[,-1]
rownames(meanStudyAbundance) <- aggregate(controlSpProfile[,1:196],by=list(controlSpProfile$study_name),FUN=median)[,1]

pcoMeanStudyAbundance <- dudi.pco(vegdist(meanStudyAbundance,method="bray"),scannf=FALSE,nf=3)

TOMOverlap <- cor(as.matrix(vegdist(pcoCoreKeyStoneDf$li,method="euclidean")),as.matrix(vegdist(pcoMeanStudyAbundance$li,method="euclidean")))

df_pcoMeanStudyAbundance <- data.frame(A1=pcoMeanStudyAbundance$li[,1],A2=pcoMeanStudyAbundance$li[,2],CohortLifeStyle=CohortLifeStyle[rownames(pcoMeanStudyAbundance$li)])

df_pcoCoreKeyStoneDf <- data.frame(A1=pcoCoreKeyStoneDf$li[,1],A2=pcoCoreKeyStoneDf$li[,2],CohortLifeStyle=CohortLifeStyle[rownames(pcoCoreKeyStoneDf$li)])

#heatmap.2(TOMOverlap,density="none",trace="none",margins=c(8,15),srtRow=0.05,srtCol=89.5,lhei=c(0.5,5),cexRow=0.5,cexCol=0.7,col=brewer.pal(8,"BrBG"))

ggplot(df_pcoMeanStudyAbundance,aes(x=A1,y=A2,color=CohortLifeStyle))+geom_point(size=5)+geom_text_repel(label=rownames(df_pcoMeanStudyAbundance),size=3,box.padding=0.01,max.overlaps=20)+theme_bw()+theme(axis.text.x=element_text(size=15))

ggplot(df_pcoCoreKeyStoneDf,aes(x=A1,y=A2,color=CohortLifeStyle))+geom_point(size=5)+geom_text_repel(label=rownames(df_pcoCoreKeyStoneDf),size=3,box.padding=0.01,max.overlaps=20)+theme_bw()+theme(axis.text.x=element_text(size=15))

Procuste_MeanStudyAbundance_CoreKeyStoneDf <- procuste.randtest(df_pcoMeanStudyAbundance[, 1:2], df_pcoCoreKeyStoneDf[rownames(df_pcoMeanStudyAbundance), 1:2], nrepet = 999)
## 
## Observation: 0.6693047
## 
## Based on 999 replicates
## Simulated p-value: 0.001
## Alternative hypothesis: greater
##
## Std.Obs  Expectation     Variance
## 10.259806133  0.134364006  0.002718523

ggplot(df_pcoCoreKeyStoneDf,aes(x=A1,y=A2,color=CohortLifeStyle))+geom_point()+geom_text_repel(label=rownames(df_pcoCoreKeyStoneDf),size=2.7,box.padding=0.01,max.overlaps=20)+theme_bw()+theme(axis.text.x=element_text(size=15))

#load(paste0(dir_name,"c59_MediancorDfAndPredictedR2df.RData"))

#predictedMajorCoreKeyStoneDf <- apply(prevalentDf[MajorCoreKeyStone,rownames(predictedR2Df)],2,function(x)(ifelse(x>0,1,0))) * apply(predictedR2Df[,MajorCoreKeyStone],1,function(x)(ifelse(x>=0.70,1,0)))

#predictedMajorCoreKeyStoneDf <- t(predictedMajorCoreKeyStoneDf)

#pcoPredictedMajorCoreKeyStoneDf <- dudi.pco(vegdist(predictedMajorCoreKeyStoneDf[],method="jaccard"),scannf=FALSE,nf=3)

#df_pcoPredictedMajorCoreKeyStone <- data.frame(A1=pcoPredictedMajorCoreKeyStoneDf$li[,1],A2=pcoPredictedMajorCoreKeyStoneDf$li[,2],CohortLifeStyle=CohortLifeStyle[rownames(pcoPredictedMajorCoreKeyStoneDf$li)])

#ggplot(df_pcoPredictedMajorCoreKeyStone,aes(x=A1,y=A2,color=CohortLifeStyle))+geom_point()+geom_text_repel(label=rownames(df_pcoPredictedMajorCoreKeyStone),size=2.7,box.padding=0.01,max.overlaps=20)+theme_bw()+theme(axis.text.x=element_text(size=15))

pcoMajorCoreKeyStoneActual <- dudi.pco(vegdist(t(CoreKeyStoneDf[SelectSpecies,]),method="jaccard"),scannf=FALSE,nf=3)

df_MajorCoreKeyStoneActual <- data.frame(A1=pcoMajorCoreKeyStoneActual$li[,1],A2=pcoMajorCoreKeyStoneActual$li[,2],CohortLifeStyle=CohortLifeStyle[rownames(pcoMajorCoreKeyStoneActual$li)])

ggplot(df_MajorCoreKeyStoneActual,aes(x=A1,y=A2,color=CohortLifeStyle))+geom_point(size=5)+geom_text_repel(label=rownames(df_MajorCoreKeyStoneActual),size=3,box.padding=0.01,max.overlaps=20)+theme_bw()+theme(axis.text.x=element_text(size=15))

#pcoMajorCoreKeyStonePredicted <- dudi.pco(vegdist(predictedR2Df[,SelectSpecies],method="jaccard"),scannf=FALSE,nf=3)

#df_MajorCoreKeyStonePredicted <- data.frame(A1=pcoMajorCoreKeyStonePredicted$li[,1],A2=pcoMajorCoreKeyStonePredicted$li[,2],CohortLifeStyle=CohortLifeStyle[rownames(pcoMajorCoreKeyStonePredicted$li)])

#ggplot(df_MajorCoreKeyStonePredicted,aes(x=A1,y=A2,color=CohortLifeStyle))+geom_point(size=5)+geom_text_repel(label=rownames(df_MajorCoreKeyStonePredicted),size=3,box.padding=0.01,max.overlaps=20)+theme_bw()+theme(axis.text.x=element_text(size=15))

#s.class(pcoMeanStudyAbundance$li,as.factor(CohortSeqType[rownames(pcoMeanStudyAbundance$li)]),col=c("#f79361","#14af9d"),plabels.col="black",plabels.boxes.border="black", pgrid.draw= FALSE, plabels.boxes.alpha= 0.7, plabels.cex= 1.75)

adonis2(vegdist(pcoMeanStudyAbundance$li[,c(1,2)],method="euclidean")~as.factor(CohortSeqType[rownames(pcoMeanStudyAbundance$li)]))

#s.class(pcoCoreKeyStoneDf$li,as.factor(CohortSeqType[rownames(pcoCoreKeyStoneDf$li)]),col=c("#f79361","#14af9d"),plabels.col="black",plabels.boxes.border="black", pgrid.draw= FALSE, plabels.boxes.alpha= 0.7, plabels.cex= 1.75)

adonis2(vegdist(pcoCoreKeyStoneDf$li[,c(1,2)],method="euclidean")~as.factor(CohortSeqType[rownames(pcoCoreKeyStoneDf$li)]))


#s.class(pcoMeanStudyAbundance$li,as.factor(CohortSeqType[rownames(pcoMeanStudyAbundance$li)]),col=c("#f79361","#14af9d"),plabels.col="black",plabels.boxes.border="black", pgrid.draw= FALSE, plabels.boxes.alpha= 0.7, plabels.cex= 1.75)

adonis2(vegdist(pcoMeanStudyAbundance$li[,c(1,2)],method="euclidean")~as.factor(CohortSeqType[rownames(pcoMeanStudyAbundance$li)]))

#s.class(pcoCoreKeyStoneDf$li,as.factor(CohortSeqType[rownames(pcoCoreKeyStoneDf$li)]),col=c("#f79361","#14af9d"),plabels.col="black",plabels.boxes.border="black", pgrid.draw= FALSE, plabels.boxes.alpha= 0.7, plabels.cex= 1.75)

adonis2(vegdist(pcoCoreKeyStoneDf$li[,c(1,2)],method="euclidean")~as.factor(CohortSeqType[rownames(pcoCoreKeyStoneDf$li)]))

#s.class(pcoCoreKeyStoneDf$li,as.factor(CohortLifeStyle[rownames(pcoCoreKeyStoneDf$li)]),col=c("#e96479","#fed837","#70ad47"),plabels.col="black",plabels.boxes.border="black", pgrid.draw= FALSE, plabels.boxes.alpha= 0.7, plabels.cex= 1.75)

adonis2(vegdist(pcoCoreKeyStoneDf$li[,c(1,2)],method="euclidean")~as.factor(CohortLifeStyle[rownames(pcoCoreKeyStoneDf$li)]))

#s.class(pcoMeanStudyAbundance$li,as.factor(CohortLifeStyle[rownames(pcoMeanStudyAbundance$li)]),col=c("#e96479","#fed837","#70ad47"),plabels.col="black",plabels.boxes.border="black", pgrid.draw= FALSE, plabels.boxes.alpha= 0.7, plabels.cex= 1.75)

adonis2(vegdist(pcoMeanStudyAbundance$li[,c(1,2)],method="euclidean")~as.factor(CohortSeqType[rownames(pcoMeanStudyAbundance$li)]))

write.table(r2Df,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\r2Df_Matrix.txt",sep="\t",quote=FALSE)
write.table(prevalentDf,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\Prevalence_Matrix.txt",sep="\t",quote=FALSE)
write.table(CoreKeyStoneDf,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\CoreKeyStone_Matrix.txt",sep="\t",quote=FALSE)
write.table(prDf,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\probability_Matrix.txt",sep="\t",quote=FALSE)

RuralTribal_CoreKeyStone <- names(which(rowSums(CoreKeyStoneDf[,names(which(CohortLifeStyle=="RuralTribal"))])/length(names(which(CohortLifeStyle=="RuralTribal"))) >= 1/6))
temp_df_RuralTribal <- data.frame(DetectedStudies = rowSums(CoreKeyStoneDf[RuralTribal_CoreKeyStone,names(which(CohortLifeStyle=="RuralTribal"))])[order(rowSums(CoreKeyStoneDf[RuralTribal_CoreKeyStone,names(which(CohortLifeStyle=="RuralTribal"))]))])
temp_df_RuralTribal$Index <- 1:length(RuralTribal_CoreKeyStone)
ggplot(temp_df_RuralTribal,aes(y=DetectedStudies,x=Index))+geom_point(color="deepskyblue")+geom_text_repel(label=rownames(temp_df_RuralTribal),max.overlaps=30,size=4.2,box.padding=0.05,min.segment.length=1)+theme_bw()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))

IndustrializedUrban_CoreKeyStone <- names(which(rowSums(CoreKeyStoneDf[,names(which(CohortLifeStyle=="IndustrializedUrban"))])/length(names(which(CohortLifeStyle=="IndustrializedUrban"))) >= 1/6))
temp_df_IndustrializedUrban <- data.frame(DetectedStudies = rowSums(CoreKeyStoneDf[IndustrializedUrban_CoreKeyStone,names(which(CohortLifeStyle=="IndustrializedUrban"))])[order(rowSums(CoreKeyStoneDf[IndustrializedUrban_CoreKeyStone,names(which(CohortLifeStyle=="IndustrializedUrban"))]))])
temp_df_IndustrializedUrban$Index <- 1:length(IndustrializedUrban_CoreKeyStone)
ggplot(temp_df_IndustrializedUrban,aes(y=DetectedStudies,x=Index))+geom_point(color="deepskyblue")+geom_text_repel(label=rownames(temp_df_IndustrializedUrban),max.overlaps=30,size=4.2,box.padding=0.05,min.segment.length=1)+theme_bw()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))

UrbanRuralMixed_CoreKeyStone <- names(which(rowSums(CoreKeyStoneDf[,names(which(CohortLifeStyle=="UrbanRuralMixed"))])/length(names(which(CohortLifeStyle=="UrbanRuralMixed"))) >= 1/6))
temp_df_UrbanRuralMixed <- data.frame(DetectedStudies = rowSums(CoreKeyStoneDf[UrbanRuralMixed_CoreKeyStone,names(which(CohortLifeStyle=="UrbanRuralMixed"))])[order(rowSums(CoreKeyStoneDf[UrbanRuralMixed_CoreKeyStone,names(which(CohortLifeStyle=="UrbanRuralMixed"))]))])
temp_df_UrbanRuralMixed$Index <- 1:length(UrbanRuralMixed_CoreKeyStone)
ggplot(temp_df_UrbanRuralMixed,aes(y=DetectedStudies,x=Index))+geom_point(color="deepskyblue")+geom_text_repel(label=rownames(temp_df_UrbanRuralMixed),max.overlaps=30,size=4.2,box.padding=0.05,min.segment.length=1)+theme_bw()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))

InfluenceScoreCohortLifeStyle <- data.frame("IndustrializedUrban" = rowSums(CoreKeyStoneDf[,names(which(CohortLifeStyle=="IndustrializedUrban"))])/length(names(which(CohortLifeStyle=="IndustrializedUrban"))), "UrbanRuralMixed" = rowSums(CoreKeyStoneDf[,names(which(CohortLifeStyle=="UrbanRuralMixed"))])/length(names(which(CohortLifeStyle=="UrbanRuralMixed"))), "RuralTribal" = rowSums(CoreKeyStoneDf[,names(which(CohortLifeStyle=="RuralTribal"))])/length(names(which(CohortLifeStyle=="RuralTribal"))))

SelectSpecies <- names(which(apply(InfluenceScoreCohortLifeStyle,1,function(x)(length(x[x>0])))>=2))

AcrossCohortCoeffVariance <- apply(InfluenceScoreCohortLifeStyle[SelectSpecies,],1,function(x)(sd(x)/mean(x)))

save(AcrossCohortCoeffVariance,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\FunctionalProfiling\\AcrossCohortCoeffVariance.RData")

CoreKeyStoneScore <- rowSums(CoreKeyStoneDf)/72

InfluenceScore <- rank_scale(rowSums(CoreKeyStoneDf)/72)

VariationWithCoreKeyStoneRank <- data.frame("CV_Prevalence"=apply(prevalentDf,1,function(x)(sd(x)/mean(x))),CV_RSquared=apply(r2Df,1,function(x)(sd(x)/mean(x))),CoreKeyStoneRank=InfluenceScore)

controlSpProfile$cohort_lifestyle <- as.vector(CohortLifeStyle[controlSpProfile$study_name])

VariationWithCoreKeyStoneRank$VariationPrevCohortLifestyle <- NA
VariationWithCoreKeyStoneRank$VariationInflCohortLifestyle <- NA

for(i in 1:201)
{
	summary_lm <- summary(lm(as.numeric(prevalentDf[i,])~as.numeric(factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")))))
	VariationWithCoreKeyStoneRank[i,"VariationPrevCohortLifestyle"] <- summary_lm$coefficients[2,3]
	
	summary_lm <- summary(lm(as.numeric(r2Df[i,])~as.numeric(factor(CohortLifeStyle,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal")))))
	VariationWithCoreKeyStoneRank[i,"VariationInflCohortLifestyle"] <- summary_lm$coefficients[2,3]
	
}

CohortLifeStyleGroupedDetection <- t(apply(CoreKeyStoneDf,1,function(x)(tapply(x,CohortLifeStyle,sum))))

VariationWithCoreKeyStoneRank$DetectionVarCohortLifeStyle <- apply(t(apply(CohortLifeStyleGroupedDetection,1,function(x)(x/c(50,11,11)))),1,sd)/apply(t(apply(CohortLifeStyleGroupedDetection,1,function(x)(x/c(50,11,11)))),1,mean)

ggplot(VariationWithCoreKeyStoneRank[!is.nan(VariationWithCoreKeyStoneRank[,6]),],aes(x=CoreKeyStoneRank,y=DetectionVarCohortLifeStyle))+geom_point()+theme_bw()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+xlab("")+ylab("")

######### 15-10-2024 (OS)
######### Finding coverage of MajorCoreKeyStone as a corekeystone oin each study 
PercentageMajorSpeciesCovered <- numeric(ncol(CoreKeyStoneDf))
# Loop over each study (column) to calculate the percentage covered
for (i in seq_along(colnames(CoreKeyStoneDf))) {
  study <- colnames(CoreKeyStoneDf)[i]
  # species as a corekeystone
  speciesInStudy <- rownames(CoreKeyStoneDf)[CoreKeyStoneDf[, study] != 0]
  # number of species as a Majorcorekeystone
  numCovered <- length(intersect(speciesInStudy, MajorCoreKeyStone))
  # total number of species as a corekeystone
  totalSpecies <- length(speciesInStudy)
  # Coverage
  PercentageMajorSpeciesCovered[i] <- (numCovered / totalSpecies)
}

# Create a data frame for coverage of Majorcorekeystone as a corekeystone in  each study
PercentageMajorSpeciesCovered_Df <- data.frame(Study = colnames(CoreKeyStoneDf), PercentageCovered = PercentageMajorSpeciesCovered)
ggplot(PercentageMajorSpeciesCovered_Df, aes(x = PercentageCovered, y = Study)) +
  geom_point() + # Adjust point transparency and size
  theme(axis.text.y = element_text(size = 9)) +
  labs(x = "Coverage", y = "Study", title = "Coverage of Species Across Studies")

save(prDf,prevalentDf,r2Df,controlSpProfileAll,controlSpProfile,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\ThresholdDetermintation.RData")

save(InfluenceScore,MajorCoreKeyStone,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\InfluenceScore.RData")		

