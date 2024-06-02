#Authors: Tarini Shankar Ghosh
#The code investigates seven additional cohorts (completely distinct from the ones considered for computing the HACK Scores
#The code computes the microbiome-level HACK Scores for the samples in these cohorts and performs various comparative analyses of these with health status and longitudinal stability

library(vegan)
library(dunn.test)
library(ggplot2)
library(beanplot)
library(DescTools)

## Please change the directory name to the folder location of your data
dir_name <- "G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\"
## All input data are provided as work-spaces in github
load(paste0(dir_name,"Olssen_species.RData"))
load(paste0(dir_name,"Olssen_WGS_Species.RData"))
load(paste0(dir_name,"Pang_species.RData"))
load(paste0(dir_name,"Xu_species.RData"))
load(paste0(dir_name,"SpeciesScores.RData"))
load(paste0(dir_name,"Song_species.RData"))
load(paste0(dir_name,"Flemer_species.RData"))
load(paste0(dir_name,"Nagata_species.RData"))
load(paste0(dir_name,"Luan_species.RData"))
load(paste0(dir_name,"Saleem_species.RData"))
load(paste0(dir_name,"MicroDiab_Denmark_species.RData"))
load(paste0(dir_name,"Saleem_species.RData"))
load(paste0(dir_name,"Parbie_species.RData"))
load(paste0(dir_name,"StudyWise3RValidationData.RData"))

TopHACKs <- rownames(SpeciesScores[SpeciesScores[,4]>=0.64,])
MidHACKs <- rownames(SpeciesScores[(SpeciesScores[,4]<0.64)&(SpeciesScores[,4]>=0.33),])
LowHACKs <- rownames(SpeciesScores[SpeciesScores[,4]<0.33,])
HACKs <- rownames(tail(SpeciesScores,18))

source(paste0(dir_name,"code_library.R"))
SpeciesScores <- SpeciesScores[order(SpeciesScores[,4]),]

threshold = 0.0001

bray_followup <- function(data,follow_up)
{
	follow_up_dist <- as.data.frame(matrix(NA,nrow(data),1))
	rownames(follow_up_dist) <- rownames(follow_up)
	for(i in 1:nrow(follow_up))
	{
		T0 <- follow_up[i,1]
		T1 <- follow_up[i,2]
		if(!is.na(T1))
		{
			follow_up_dist[i,1] <- as.numeric(vegdist(data[c(T1,T0),],method="bray"))
		}
	}
	return(follow_up_dist)
}

aitchison_followup <- function(data,follow_up)
{
	data_clr <- as.data.frame(as.matrix(clr(data+0.00001)))
	data_clr <- as.data.frame(t(apply(data_clr,1,function(x)(x-min(x)))))
	follow_up_dist <- as.data.frame(matrix(NA,nrow(data),1))
	rownames(follow_up_dist) <- rownames(follow_up)
	for(i in 1:nrow(follow_up))
	{
		T0 <- follow_up[i,1]
		T1 <- follow_up[i,2]
		if(!is.na(T1))
		{
			follow_up_dist[i,1] <- as.numeric(vegdist(data_clr[c(T1,T0),],method="euclidean"))
		}
	}
	return(follow_up_dist)
}

compute_enrichment <- function(data)
{
	data_new <- apply(data,1,function(x)(mean(x)*(1-Gini(x))))
	return(data_new)
}

## Olssen 16S Validation
Olssen_data_norm <- Olssen_et_al_species[,intersect(rownames(SpeciesScores),colnames(Olssen_et_al_species))]
Olssen_data_norm <- Olssen_data_norm/rowSums(Olssen_data_norm)

#Olssen_mHACK <- compute_enrichment(apply(Olssen_data_norm[,intersect(colnames(Olssen_data_norm),TopHACKs)],2,rank_scale))

Olssen_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Olssen_data_norm[,intersect(TopHACKs,colnames(Olssen_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Olssen_data_norm[,intersect(MidHACKs,
colnames(Olssen_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Olssen_data_norm[,intersect(LowHACKs,colnames(Olssen_data_norm))],2,
rank_scale)))

Olssen_full_data_norm <- Olssen_data_norm

Olssen_follow_up <- Olssen_et_al_FollowUp

Olssen_bray_follow_up <- bray_followup(Olssen_data_norm,Olssen_follow_up)
Olssen_aitchison_follow_up <- aitchison_followup(Olssen_et_al_species,Olssen_follow_up)

Olssen_data_norm <- apply(Olssen_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Olssen_mHACK <- as.data.frame(as.matrix(Olssen_data_norm[,intersect(colnames(Olssen_data_norm),rownames(SpeciesScores))]) %*% as.matrix(SpeciesScores[intersect(colnames(Olssen_data_norm),rownames(SpeciesScores)),4]))/rowSums(Olssen_data_norm)

#Olssen_mHACK <- cor(t(Olssen_data_norm),SpeciesScores[colnames(Olssen_data_norm),4],method="spearman")

Olssen_df <- data.frame(mHACK=Olssen_mHACK,bray_follow_up=Olssen_bray_follow_up[,1],aitchison_follow_up=Olssen_aitchison_follow_up[,1],shannon=vegan::diversity(Olssen_data_norm))

Olssen_df <- Olssen_df[!is.na(Olssen_df[,2]),]

boxplot(Olssen_df[,2]~cut(Olssen_df[,1],breaks=quantile(Olssen_df[,1],prob=c(0,0.33,0.67,1)),include.lowest=TRUE),outline=FALSE,range=0.75,ylab="Bray Curtis Follow-up Distance",xlab="mHACK Index",names=c("Q1","Q2","Q3"),col=c("coral1","bisque1","darkseagreen1"),cex.axis=1.3)

dunn.test(Olssen_df[,2],cut(Olssen_df[,1],breaks=quantile(Olssen_df[,1],prob=c(0,0.33,0.67,1)),include.lowest=TRUE))

##Olssen WGS validation
Olsson_wgs_full_data_norm <- Olsson_et_al_wgs_species/rowSums(Olsson_et_al_wgs_species)
Olsson_wgs_data_norm <- Olsson_et_al_wgs_species[,intersect(rownames(SpeciesScores),colnames(Olsson_et_al_wgs_species))]

Olsson_wgs_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Olsson_wgs_data_norm[,intersect(TopHACKs,colnames(Olsson_wgs_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Olsson_wgs_data_norm[,intersect(MidHACKs,
colnames(Olsson_wgs_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Olsson_wgs_data_norm[,intersect(LowHACKs,colnames(Olsson_wgs_data_norm))],2,
rank_scale)))

Olsson_wgs_data_norm <- Olsson_wgs_data_norm/rowSums(Olsson_wgs_data_norm)
Olsson_wgs_full_data_norm <-Olsson_wgs_data_norm 
Olsson_wgs_follow_up <- Olsson_et_al_wgs_follow_up_metadata

Olsson_wgs_bray_follow_up <- bray_followup(Olsson_wgs_full_data_norm,Olsson_wgs_follow_up)
Olsson_wgs_aitchison_follow_up <- aitchison_followup(Olsson_wgs_full_data_norm,Olsson_wgs_follow_up)

#Olsson_wgs_mHACK <- as.data.frame(as.matrix(Olsson_wgs_data_norm[,intersect(colnames(Olsson_wgs_data_norm),rownames(SpeciesScores))]) %*% as.matrix(SpeciesScores[intersect(colnames(Olsson_wgs_data_norm),rownames(SpeciesScores)),4]))/(apply(Olsson_wgs_data_norm,2,function(x)(ifelse(x>0,1,0))) %*% as.matrix(SpeciesScores[intersect(colnames(Olsson_wgs_data_norm),rownames(SpeciesScores)),4]))

Olsson_wgs_data_norm <- apply(Olsson_wgs_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Olsson_wgs_mHACK <- as.data.frame(as.matrix(Olsson_wgs_data_norm) %*% as.matrix(SpeciesScores[intersect(colnames(Olsson_wgs_data_norm),rownames(SpeciesScores)),4]))/rowSums(Olsson_wgs_data_norm)

Olsson_wgs_df <- data.frame(mHACK=Olsson_wgs_mHACK,bray_follow_up=Olsson_wgs_bray_follow_up[,1],aitchison_follow_up=Olsson_wgs_aitchison_follow_up[,1],shannon=vegan::diversity(Olsson_wgs_data_norm))

boxplot(Olsson_wgs_df[,2]~cut(Olsson_wgs_df[,1],breaks=quantile(Olsson_wgs_df[,1],prob=c(0,0.33,0.67,1)),include.lowest=TRUE),outline=FALSE,range=0.75,ylab="Bray Curtis Follow-up Distance",xlab="mHACK Index",names=c("Q1","Q2","Q3"),col=c("coral1","bisque1","darkseagreen1"),cex.axis=1.3)

dunn.test(Olsson_wgs_df[,2],cut(Olsson_wgs_df[,1],breaks=quantile(Olsson_wgs_df[,1],prob=c(0,0.33,0.67,1)),include.lowest=TRUE))

## Pang et al Validation (Cross-Sectional)
Pang_full_data <- Pang_et_al_species
Pang_full_data_norm <- Pang_et_al_species/rowSums(Pang_et_al_species)
Pang_data <- Pang_et_al_species[,intersect(colnames(Pang_et_al_species),rownames(SpeciesScores))]
Pang_data_norm <- Pang_data/rowSums(Pang_data)

Pang_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Pang_data_norm[,intersect(TopHACKs,colnames(Pang_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Pang_data_norm[,intersect(MidHACKs,
colnames(Pang_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Pang_data_norm[,intersect(LowHACKs,colnames(Pang_data_norm))],2,
rank_scale)))

Pang_full_data_norm <- Pang_data_norm
Pang_data_norm <- apply(Pang_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Pang_mHACK <- as.data.frame(as.matrix(Pang_data_norm[,intersect(colnames(Pang_data_norm),rownames(SpeciesScores))]) %*% as.matrix(SpeciesScores[intersect(colnames(Pang_data_norm),rownames(SpeciesScores)),4]))/rowSums(Pang_data_norm)

#Pang_mHACK <- Pang_mHACK/rowSums(Pang_data_norm[,intersect(colnames(Pang_data_norm),rownames(SpeciesScores))])

#Pang_mHACK <- cor(t(Pang_data_norm),SpeciesScores[colnames(Pang_data_norm),4],method="spearman")

Pang_CrossSectional_df <- data.frame(mHack=Pang_HACKGrouped[intersect(rownames(Pang_CrossSectional_Metadata),rownames(Pang_data_norm)),1],Age=Pang_CrossSectional_Metadata[intersect(rownames(Pang_CrossSectional_Metadata),rownames(Pang_data_norm)),c("Age")],Health=Pang_CrossSectional_Metadata[intersect(rownames(Pang_CrossSectional_Metadata),rownames(Pang_data_norm)),c("Health_status")],Group=Pang_CrossSectional_Metadata[intersect(rownames(Pang_CrossSectional_Metadata),rownames(Pang_data_norm)),c("Group")],row.names=intersect(rownames(Pang_CrossSectional_Metadata),rownames(Pang_data_norm)))

boxplot(Pang_CrossSectional_df[,1]~factor(Pang_CrossSectional_df[,3],levels=c("F","LH","H")),names=c("Frail","Less Healthy","Healthy"),ylab="mHACK Index",xlab="Elderly Health Status",cex.axis=1.3,outline=FALSE,range=0.5,col=c("skyblue","yellowgreen","coral1"))

dunn.test(Pang_CrossSectional_df[,1],factor(Pang_CrossSectional_df[,3],levels=c("F","LH","H")))

boxplot(Pang_CrossSectional_df[,1]~factor(Pang_CrossSectional_df[,4],levels=c("20-44","45-65","66-85","90-99","100-117")),range=0.5,outline=FALSE,names=c("20-44","45-65","66-85","90-99","100-117"),xlab="Age-Group",ylab="mHACK",cex.axis=1.3,col=c("cadetblue1","aquamarine1","palegreen1","bisque1","coral1"))

dunn.test(Pang_CrossSectional_df[,1],factor(Pang_CrossSectional_df[,4],levels=c("20-44","45-65","66-85","90-99","100-117")))

#Pang et al Validation (Longitudinal Data)

temp_Pang_Followup_Metadata <- Pang_Followup_Metadata[intersect(rownames(Pang_data_norm),rownames(Pang_Followup_Metadata)),]
temp_Pang_data_norm <- Pang_data_norm[intersect(rownames(Pang_data_norm),rownames(Pang_Followup_Metadata)),]

T0_rows <- rownames(temp_Pang_Followup_Metadata)[(temp_Pang_Followup_Metadata[,2] %in% rownames(Pang_data_norm))&(temp_Pang_Followup_Metadata[,1] %in% rownames(Pang_data_norm))]
T1_rows <- paste0(T0_rows,"-2")

Pang_bray_follow_up <- bray_followup(Pang_full_data_norm[c(T0_rows,T1_rows),],Pang_Followup_Metadata[c(T0_rows,T1_rows),])
Pang_aitchison_follow_up <- aitchison_followup(Pang_full_data_norm[c(T0_rows,T1_rows),],Pang_Followup_Metadata[c(T0_rows,T1_rows),])

Pang_longitudinal_df <- data.frame(mHACK=Pang_mHACK[rownames(Pang_bray_follow_up),1],bray_follow_up=Pang_bray_follow_up[,1],shannon=vegan::diversity(Pang_full_data_norm[rownames(Pang_bray_follow_up),]),row.names=rownames(Pang_bray_follow_up))

Pang_longitudinal_df <- Pang_longitudinal_df[!is.na(Pang_longitudinal_df[,2]),]

ggplot(Pang_longitudinal_df,aes(x=mHACK,y=bray_follow_up)) + geom_point(size=2,pch=19,color="blue") + geom_point(data=subset(Pang_longitudinal_df,((Pang_longitudinal_df[,1]>=0.009)&(Pang_longitudinal_df[,2] >= 0.6))),color="red") + geom_smooth(data=subset(Pang_longitudinal_df,!(((Pang_longitudinal_df[,1]>=0.009)&(Pang_longitudinal_df[,2] >= 0.6)))),method="lm",alpha=0.5) + theme_bw() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

ggplot(Pang_longitudinal_df,aes(x=mHACK,y=bray_follow_up)) + geom_point(size=2,pch=19,color="blue") + geom_smooth(method="lm",alpha=0.5) + theme_bw() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

t <- subset(Pang_longitudinal_df,!((mHACK>=0.009)&(bray_follow_up >= 0.6)))[,3]
t1 <- subset(Pang_longitudinal_df,(mHACK>=0.009)&(bray_follow_up >= 0.6))[,3]

boxplot(t,t1,col=c("blue","red"))

cor.test(Pang_longitudinal_df[!((Pang_longitudinal_df[,1]>=0.009)&(Pang_longitudinal_df[,2]>=0.6)),1],Pang_longitudinal_df[!((Pang_longitudinal_df[,1]>=0.009)&(Pang_longitudinal_df[,2]>=0.6)),2],method="spearman")

#Xu et al Validation
Xu_data_norm <- Xu_et_al_species[,intersect(rownames(SpeciesScores),colnames(Xu_et_al_species))]/rowSums(Xu_et_al_species[,intersect(rownames(SpeciesScores),colnames(Xu_et_al_species))])

Xu_full_data_norm <- as.data.frame(Xu_data_norm)

Xu_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Xu_full_data_norm[,intersect(TopHACKs,colnames(Xu_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Xu_full_data_norm[,intersect(MidHACKs,
colnames(Xu_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Xu_full_data_norm[,intersect(LowHACKs,colnames(Xu_full_data_norm))],2,
rank_scale)))

#Xu_diversity <- vegan::diversity(Xu_data_norm)

Xu_data_norm <- apply(Xu_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

#Xu_mHACK <- Xu_diversity * (as.matrix(Xu_data_norm) %*% SpeciesScores[colnames(Xu_data_norm),4])/rowSums(Xu_data_norm)

Xu_mHACK <- (as.matrix(Xu_data_norm) %*% SpeciesScores[colnames(Xu_data_norm),4])/rowSums(Xu_data_norm)

Xu_full_df <- data.frame(mHACK=Xu_mHACK[,1],MMSE=Xu_et_al_metadata[rownames(Xu_mHACK),1],Barthel_Score=Xu_et_al_metadata[rownames(Xu_mHACK),2],Numb_Disorders=Xu_et_al_metadata[rownames(Xu_mHACK),3],Age=Xu_et_al_metadata[rownames(Xu_mHACK),"Age"])

boxplot(Xu_full_df[,"mHACK"]~cut(Xu_full_df[,"Age"],breaks=c(50,60,90,100),include.lowest=TRUE),outline=FALSE,range=0.5,ylab="mHACK",xlab="Age-Groups",cex.axis=1.3,col=c("cadetblue1","aquamarine1","bisque1"))

dunn.test(Xu_full_df[,"mHACK"],cut(Xu_full_df[,"Age"],breaks=c(50,60,90,100),include.lowest=TRUE))

#Only taking those with defined MMSE values
Xu_df <- Xu_full_df[!is.na(Xu_full_df[,2])&(Xu_full_df[,2]>0),]

boxplot(Xu_df[,"mHACK"]~cut(Xu_df[,"Numb_Disorders"],breaks=c(-1,0,1,2,5),include.lowest=TRUE),outline=FALSE,range=0.5,ylab="mHACK",xlab="Number of Diseases/Complications",cex.axis=1.3,col=c("skyblue","palegreen","thistle1","sienna1"),names=c("0","1","2",">2"))

dunn.test(Xu_df[,"mHACK"],cut(Xu_df[,"Numb_Disorders"],breaks=c(-1,0,1,2,5),include.lowest=TRUE))

Healthy_Elderly <- rownames(Xu_df[((Xu_df$MMSE >= 27)&(Xu_df$Barthel_Score == 100)),])
Not_Healthy_Elderly <- rownames(Xu_df[!((Xu_df$MMSE >= 27)&(Xu_df$Barthel_Score == 100)),])

boxplot(Xu_df[Healthy_Elderly,"mHACK"],Xu_df[Not_Healthy_Elderly,"mHACK"],names=c("Normal","Less_Healthy"),xlab="Elderly Health Status",ylab="mHACK",cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))

wilcox.test(Xu_df[Healthy_Elderly,"mHACK"],Xu_df[Not_Healthy_Elderly,"mHACK"])

## Song et al Validation (Cross-Sectional)
Song_data_norm <- Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores))]/rowSums(Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores))])

Song_full_data_norm <- Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores))]/rowSums(Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores))])

Song_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Song_full_data_norm[,intersect(TopHACKs,colnames(Song_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Song_full_data_norm[,intersect(MidHACKs,
colnames(Song_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Song_full_data_norm[,intersect(LowHACKs,colnames(Song_full_data_norm))],2,
rank_scale)))

Song_data_norm <- apply(Song_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Song_mHACK <- (Song_data_norm %*% SpeciesScores[colnames(Song_data_norm),4])/rowSums(Song_data_norm)

boxplot(Song_mHACK[which(Song_et_al_metadata[rownames(Song_mHACK),"group"]=="Normotension"),1],Song_mHACK[which(Song_et_al_metadata[rownames(Song_mHACK),"group"]=="Hypertension"),1],names=c("Normal","Hypertension"),xlab="Disease Status",ylab="mHACK",cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))

wilcox.test(Song_mHACK[which(Song_et_al_metadata[rownames(Song_mHACK),"group"]=="Normotension"),1],Song_mHACK[which(Song_et_al_metadata[rownames(Song_mHACK),"group"]=="Hypertension"),1])


## Flemer et al Validation (Cross-Sectional)
Flemer_data_norm <- Flemer_et_al_species[,intersect(colnames(Flemer_et_al_species),rownames(SpeciesScores))]/rowSums(Flemer_et_al_species[,intersect(colnames(Flemer_et_al_species),rownames(SpeciesScores))])

Flemer_full_data_norm <- as.data.frame(Flemer_data_norm)

Flemer_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Flemer_full_data_norm[,intersect(TopHACKs,colnames(Flemer_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Flemer_full_data_norm[,intersect(MidHACKs,
colnames(Flemer_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Flemer_full_data_norm[,intersect(LowHACKs,colnames(Flemer_full_data_norm))],2,
rank_scale)))

Flemer_data_norm <- apply(Flemer_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Flemer_mHACK <- (Flemer_data_norm %*% SpeciesScores[colnames(Flemer_data_norm),4])/rowSums(Flemer_data_norm)

boxplot(Flemer_mHACK[which(Flemer_et_al_metadata[rownames(Flemer_mHACK),"disease"]=="control"),1],Flemer_mHACK[which(Flemer_et_al_metadata[rownames(Flemer_mHACK),"disease"]=="CRC"),1],names=c("Normal","CRC"),xlab="Disease Status",ylab="mHACK",cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))

wilcox.test(Flemer_mHACK[which(Flemer_et_al_metadata[rownames(Flemer_mHACK),"disease"]=="control"),1],Flemer_mHACK[which(Flemer_et_al_metadata[rownames(Flemer_mHACK),"disease"]=="CRC"),1])

## Nagata et al Validation (Cross-Sectional)
Nagata_data_norm <- Nagata_et_al_species[,intersect(colnames(Nagata_et_al_species),rownames(SpeciesScores))]/rowSums(Nagata_et_al_species[,intersect(colnames(Nagata_et_al_species),rownames(SpeciesScores))])

Nagata_full_data_norm <- as.data.frame(Nagata_data_norm)

Nagata_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Nagata_full_data_norm[,intersect(TopHACKs,colnames(Nagata_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Nagata_full_data_norm[,intersect(MidHACKs,
colnames(Nagata_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Nagata_full_data_norm[,intersect(LowHACKs,colnames(Nagata_full_data_norm))],2,
rank_scale)))

Nagata_data_norm <- apply(Nagata_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Nagata_mHACK <- (Nagata_data_norm %*% SpeciesScores[colnames(Nagata_data_norm),4])/rowSums(Nagata_data_norm)

boxplot(Nagata_mHACK[which(Nagata_et_al_metadata[rownames(Nagata_mHACK),"disease_status"]=="Control"),1],Nagata_mHACK[which(Nagata_et_al_metadata[rownames(Nagata_mHACK),"disease_status"]=="PDAC"),1],names=c("Normal","PDAC"),xlab="Disease Status",ylab="mHACK",cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))

wilcox.test(Nagata_mHACK[which(Nagata_et_al_metadata[rownames(Nagata_mHACK),"disease_status"]=="Control"),1],Nagata_mHACK[which(Nagata_et_al_metadata[rownames(Nagata_mHACK),"disease_status"]=="PDAC"),1])

#Luan et al Validation
Luan_et_al_metadata$ElderlyHealthStatus <- rank_scale(Luan_et_al_metadata$ADL) - rank_scale(Luan_et_al_metadata$GDS_15)

Luan_data_norm <- Luan_et_al_species[,intersect(rownames(SpeciesScores),colnames(Luan_et_al_species))]/rowSums(Luan_et_al_species[,intersect(rownames(SpeciesScores),colnames(Luan_et_al_species))])

Luan_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Luan_data_norm[,intersect(TopHACKs,colnames(Luan_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Luan_data_norm[,intersect(MidHACKs,
colnames(Luan_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Luan_data_norm[,intersect(LowHACKs,colnames(Luan_data_norm))],2,
rank_scale)))

Luan_full_data_norm <- as.data.frame(Luan_data_norm)

Luan_data_norm <- apply(Luan_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Luan_mHACK <- (Luan_data_norm %*% SpeciesScores[colnames(Luan_data_norm),4])/rowSums(Luan_data_norm)

boxplot(Luan_mHACK[,1]~cut(Luan_et_al_metadata[rownames(Luan_mHACK),"GDS_15"],breaks=c(0,8,15),include.lowest=TRUE),cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))

wilcox.test(Luan_mHACK[,1]~cut(Luan_et_al_metadata[rownames(Luan_mHACK),"GDS_15"],breaks=c(0,8,15),include.lowest=TRUE))

#MicroDiab Denmark
MicroDiab_Denmark_data_norm <- MicroDiab_Denmark_Species[,intersect(rownames(SpeciesScores),colnames(MicroDiab_Denmark_Species))] / rowSums(MicroDiab_Denmark_Species[,intersect(rownames(SpeciesScores),colnames(MicroDiab_Denmark_Species))])

MicroDiab_Denmark_HACKGrouped <- data.frame("Top"=rowMeans(apply(
MicroDiab_Denmark_data_norm[,intersect(TopHACKs,colnames(MicroDiab_Denmark_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(MicroDiab_Denmark_data_norm[,intersect(MidHACKs,
colnames(MicroDiab_Denmark_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
MicroDiab_Denmark_data_norm[,intersect(LowHACKs,colnames(MicroDiab_Denmark_data_norm))],2,
rank_scale)))

MicroDiab_Denmark_full_data_norm <- as.data.frame(MicroDiab_Denmark_data_norm)

MicroDiab_Denmark_data_norm <- apply(MicroDiab_Denmark_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

MicroDiab_Denmark_mHACK <- as.matrix(MicroDiab_Denmark_data_norm) %*% SpeciesScores[colnames(MicroDiab_Denmark_data_norm),4]/rowSums(MicroDiab_Denmark_data_norm)

boxplot(MicroDiab_Denmark_mHACK~MicroDiab_Denmark_Metadata[rownames(MicroDiab_Denmark_mHACK),"Disease"],outline=FALSE,range=0.5,col=c("aquamarine1","bisque","coral1"),xlab="Disease Groups",ylab="mHACK",names=c("control","PreT2D","T2D"),cex.axis=1.3)

dunn.test(MicroDiab_Denmark_mHACK,MicroDiab_Denmark_Metadata[rownames(MicroDiab_Denmark_mHACK),"Disease"])

### Saleem et al Validation
Saleem_data_norm <- Saleem_et_al_species[,intersect(rownames(SpeciesScores),colnames(Saleem_et_al_species))]/rowSums(Saleem_et_al_species[,intersect(rownames(SpeciesScores),colnames(Saleem_et_al_species))])

Saleem_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Saleem_data_norm[,intersect(TopHACKs,colnames(Saleem_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Saleem_data_norm[,intersect(MidHACKs,
colnames(Saleem_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Saleem_data_norm[,intersect(LowHACKs,colnames(Saleem_data_norm))],2,
rank_scale)))

Saleem_full_data_norm <- Saleem_data_norm

Saleem_mHACK <- (as.matrix(Saleem_data_norm) %*% SpeciesScores[colnames(Saleem_data_norm),4])/rowSums(Saleem_data_norm)

boxplot(Saleem_mHACK[,1]~Saleem_et_al_metadata[rownames(Saleem_mHACK),"disease"],names=c("Normal","T2D"),xlab="Disease Status",ylab="mHACK",cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))
wilcox.test(Saleem_mHACK[,1]~Saleem_et_al_metadata[rownames(Saleem_mHACK),"disease"])

### Parbie et al Validation
Parbie_data_norm <- Parbie_et_al_species[,intersect(colnames(Parbie_et_al_species),rownames(SpeciesScores))]/rowSums(Parbie_et_al_species[,intersect(colnames(Parbie_et_al_species),rownames(SpeciesScores))])

Parbie_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Parbie_data_norm[,intersect(TopHACKs,colnames(Parbie_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Parbie_data_norm[,intersect(MidHACKs,
colnames(Parbie_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Parbie_data_norm[,intersect(LowHACKs,colnames(Parbie_data_norm))],2,
rank_scale)))

Parbie_full_data_norm <- as.data.frame(Parbie_data_norm)

Parbie_data_norm <- apply(Parbie_data_norm,2,function(x)(ifelse(x>threshold,1,0)))

Parbie_mHACK <- Parbie_data_norm %*% SpeciesScores[colnames(Parbie_data_norm),4]/rowSums(Parbie_data_norm)

boxplot(Parbie_mHACK[,1]~Parbie_et_al_metadata[rownames(Parbie_mHACK),"disease"],names=c("Normal","HIV"),xlab="Disease Status",ylab="mHACK",cex.axis=1.3,outline=FALSE,range=0.5,col=c("yellowgreen","tan1"))
wilcox.test(Parbie_mHACK[,1]~Parbie_et_al_metadata[rownames(Parbie_mHACK),"disease"])

###### Correlation with Core Influencer Score
Flemer_3R <- Flemer_3R[Flemer_3R[,2]>0,]
Flemer_3R$HACKScores <- SpeciesScores[rownames(Flemer_3R),4]
Flemer_3R$InfluenceScore <- rank_scale(Flemer_3R[,2]*Flemer_3R[,3])
ggplot(Flemer_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Flemer_C <- rownames(Flemer_3R)[Flemer_3R$CoreSp == "CoreInfluencer"]
Flemer_NC <- rownames(Flemer_3R)[Flemer_3R$CoreSp == "NonCoreInfluencer"]

Pang_3R <- Pang_3R[Pang_3R[,2]>0,]
Pang_3R$HACKScores <- SpeciesScores[rownames(Pang_3R),4]
Pang_3R$InfluenceScore <- rank_scale(Pang_3R[,2]*Pang_3R[,3])
ggplot(Pang_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Pang_C <- rownames(Pang_3R)[Pang_3R$CoreSp == "CoreInfluencer"]
Pang_NC <- rownames(Pang_3R)[Pang_3R$CoreSp == "NonCoreInfluencer"]

Xu_3R <- Xu_3R[Xu_3R[,2]>0,]
Xu_3R$HACKScores <- SpeciesScores[rownames(Xu_3R),4]
Xu_3R$InfluenceScore <- rank_scale(Xu_3R[,2]*Xu_3R[,3])
ggplot(Xu_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Xu_C <- rownames(Xu_3R)[Xu_3R$CoreSp == "CoreInfluencer"]
Xu_NC <- rownames(Xu_3R)[Xu_3R$CoreSp == "NonCoreInfluencer"]

Olssen_3R <- Olssen_3R[Olssen_3R[,2]>0,]
Olssen_3R$HACKScores <- SpeciesScores[rownames(Olssen_3R),4]
Olssen_3R$InfluenceScore <- rank_scale(Olssen_3R[,2]*Olssen_3R[,3])
ggplot(Olssen_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Olssen_C <- rownames(Olssen_3R)[Olssen_3R$CoreSp == "CoreInfluencer"]
Olssen_NC <- rownames(Olssen_3R)[Olssen_3R$CoreSp == "NonCoreInfluencer"]

Olsson_wgs_3R <- Olsson_wgs_3R[Olsson_wgs_3R[,2]>0,]
Olsson_wgs_3R$HACKScores <- SpeciesScores[rownames(Olsson_wgs_3R),4]
Olsson_wgs_3R$InfluenceScore <- rank_scale(Olsson_wgs_3R[,2]*Olsson_wgs_3R[,3])
ggplot(Olsson_wgs_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Olsson_wgs_C <- rownames(Olsson_wgs_3R)[Olsson_wgs_3R$CoreSp == "CoreInfluencer"]
Olsson_wgs_NC <- rownames(Olsson_wgs_3R)[Olsson_wgs_3R$CoreSp == "NonCoreInfluencer"]

Luan_3R <- Luan_3R[Luan_3R[,2]>0,]
Luan_3R$HACKScores <- SpeciesScores[rownames(Luan_3R),4]
Luan_3R$InfluenceScore <- rank_scale(Luan_3R[,2]*Luan_3R[,3])
ggplot(Luan_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Luan_C <- rownames(Luan_3R)[Luan_3R$CoreSp == "CoreInfluencer"]
Luan_NC <- rownames(Luan_3R)[Luan_3R$CoreSp == "NonCoreInfluencer"]

Nagata_3R <- Nagata_3R[Nagata_3R[,2]>0,]
Nagata_3R$HACKScores <- SpeciesScores[rownames(Nagata_3R),4]
Nagata_3R$InfluenceScore <- rank_scale(Nagata_3R[,2]*Nagata_3R[,3])
ggplot(Nagata_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Nagata_C <- rownames(Nagata_3R)[Nagata_3R$CoreSp == "CoreInfluencer"]
Nagata_NC <- rownames(Nagata_3R)[Nagata_3R$CoreSp == "NonCoreInfluencer"]

MicroDiab_Denmark_3R <- MicroDiab_Denmark_3R[MicroDiab_Denmark_3R[,2]>0,]
MicroDiab_Denmark_3R$HACKScores <- SpeciesScores[rownames(MicroDiab_Denmark_3R),4]
MicroDiab_Denmark_3R$InfluenceScore <- rank_scale(MicroDiab_Denmark_3R[,2]*MicroDiab_Denmark_3R[,3])
ggplot(MicroDiab_Denmark_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
MicroDiab_Denmark_C <- rownames(MicroDiab_Denmark_3R)[MicroDiab_Denmark_3R$CoreSp == "CoreInfluencer"]
MicroDiab_Denmark_NC <- rownames(MicroDiab_Denmark_3R)[MicroDiab_Denmark_3R$CoreSp == "NonCoreInfluencer"]

Saleem_3R <- Saleem_3R[Saleem_3R[,2]>0,]
Saleem_3R$HACKScores <- SpeciesScores[rownames(Saleem_3R),4]
Saleem_3R$InfluenceScore <- rank_scale(Saleem_3R[,2]*Saleem_3R[,3])
ggplot(Saleem_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Saleem_C <- rownames(Saleem_3R)[Saleem_3R$CoreSp == "CoreInfluencer"]
Saleem_NC <- rownames(Saleem_3R)[Saleem_3R$CoreSp == "NonCoreInfluencer"]

Parbie_3R <- Parbie_3R[Parbie_3R[,2]>0,]
Parbie_3R$HACKScores <- SpeciesScores[rownames(Parbie_3R),4]
Parbie_3R$InfluenceScore <- rank_scale(Parbie_3R[,2]*Parbie_3R[,3])
ggplot(Parbie_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Parbie_C <- rownames(Parbie_3R)[Parbie_3R$CoreSp == "CoreInfluencer"]
Parbie_NC <- rownames(Parbie_3R)[Parbie_3R$CoreSp == "NonCoreInfluencer"]

Song_3R <- Song_3R[Song_3R[,2]>0,]
Song_3R$HACKScores <- SpeciesScores[rownames(Song_3R),4]
Song_3R$InfluenceScore <- rank_scale(Song_3R[,2]*Song_3R[,3])
ggplot(Song_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
Song_C <- rownames(Song_3R)[Song_3R$CoreSp == "CoreInfluencer"]
Song_NC <- rownames(Song_3R)[Song_3R$CoreSp == "NonCoreInfluencer"]

beanplot(SpeciesScores[Flemer_C,4],SpeciesScores[Flemer_NC,4],SpeciesScores[Nagata_C,4],SpeciesScores[Nagata_NC,4],SpeciesScores[MicroDiab_Denmark_C,4],SpeciesScores[MicroDiab_Denmark_NC,4],SpeciesScores[Saleem_C,4],SpeciesScores[Saleem_NC,4],SpeciesScores[Parbie_C,4],SpeciesScores[Parbie_NC,4],SpeciesScores[Xu_C,4],SpeciesScores[Xu_NC,4],SpeciesScores[Pang_C,4],SpeciesScores[Pang_NC,4],SpeciesScores[Luan_C,4],SpeciesScores[Luan_NC,4],SpeciesScores[Olssen_C,4],SpeciesScores[Olssen_NC,4],SpeciesScores[Olsson_wgs_C,4],SpeciesScores[Olsson_wgs_NC,4],SpeciesScores[Song_C,4],SpeciesScores[Song_NC,4],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))


AllDetectedTaxaValidation <- names(which(table(c(rownames(Flemer_3R),rownames(Nagata_3R),rownames(MicroDiab_Denmark_3R),rownames(Saleem_3R),rownames(Parbie_3R),rownames(Xu_3R),rownames(Pang_3R),rownames(Olssen_3R),rownames(Luan_3R),rownames(Olsson_wgs_3R),rownames(Song_3R))) == 11))

CoreInfluencerDetection <- data.frame(Flemer=ifelse(Flemer_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Nagata=ifelse(Nagata_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),MicroDiab_Denmark=ifelse(MicroDiab_Denmark_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Saleem=ifelse(Saleem_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Parbie=ifelse(Parbie_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Pang=ifelse(Pang_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Xu=ifelse(Xu_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Luan=ifelse(Luan_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Olssen=ifelse(Olssen_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),Olsson_wgs=ifelse(Olsson_wgs_3R[AllDetectedTaxaValidation,"CoreSp"]=="CoreInfluencer",1,0),row.names=AllDetectedTaxaValidation)

boxplot(SpeciesScores[AllDetectedTaxaValidation,4]~cut(rowSums(CoreInfluencerDetection),breaks=c(0,1,5,11),include.lowest=TRUE),names=c("<10%","10-50%",">50%"),xlab="",ylab="",col=c("coral","bisque","darkseagreen"))


window_correlation <- function(vec1,vec2,window_size)
{
	N = length(vec1)
	M = window_size
	correlation_df <- as.data.frame(matrix(NA,(N-M+1),2))
	colnames(correlation_df) <- c("correlation","variation")
	for(i in 1:(N-M+1))
	{
		correlation_df[i,1] <- cor(vec1[i:(i+M-1)],vec2[i:(i+M-1)],method="spearman",use="pairwise.complete")
		correlation_df[i,2] <- mean(abs(vec1[i:(i+M-1)]-vec2[i:(i+M-1)]))
	}
	return(correlation_df)
}

# Stability: Bray-Curtis
corr_Olssen_bray <- corr.test(Olssen_full_data_norm,Olssen_bray_follow_up[rownames(Olssen_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")
df_corr_Olssen_bray <- data.frame("R"=corr_Olssen_bray$r,"P"=corr_Olssen_bray$p,"HACK-Index"=SpeciesScores[rownames(corr_Olssen_bray$r),4])
ggplot(df_corr_Olssen_bray,aes(x=R,y=HACK.Index))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
df_corr_Olssen_bray$dir <- ifelse(df_corr_Olssen_bray[,2]<=0.15,2*sign(df_corr_Olssen_bray[,1]),sign(df_corr_Olssen_bray[,1]))

bray_vec_Olssen_Sig_Negative <- df_corr_Olssen_bray[df_corr_Olssen_bray[,4]== -2,3]
bray_vec_Olssen_Not_Significant <- df_corr_Olssen_bray[df_corr_Olssen_bray[,3]!= -2,3]

corr_Olsson_wgs_bray <- corr.test(Olsson_wgs_full_data_norm,Olsson_wgs_bray_follow_up[rownames(Olsson_wgs_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")
df_corr_Olsson_wgs_bray <- data.frame("R"=corr_Olsson_wgs_bray$r,"P"=corr_Olsson_wgs_bray$p,"HACK-Index"=SpeciesScores[rownames(corr_Olsson_wgs_bray$r),4])
ggplot(df_corr_Olsson_wgs_bray,aes(x=R,y=HACK.Index))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
df_corr_Olsson_wgs_bray$dir <- ifelse(df_corr_Olsson_wgs_bray[,2]<=0.15,2*sign(df_corr_Olsson_wgs_bray[,1]),sign(df_corr_Olsson_wgs_bray[,1]))

bray_vec_Olsson_wgs_Sig_Negative <- df_corr_Olsson_wgs_bray[df_corr_Olsson_wgs_bray[,4]== -2,3]
bray_vec_Olsson_wgs_Not_Significant <- df_corr_Olsson_wgs_bray[df_corr_Olsson_wgs_bray[,3]!= -2,3]

corr_Pang_bray <- corr.test(Pang_full_data_norm,Pang_bray_follow_up[rownames(Pang_full_data_norm),],method="spearman",use="pairwise.complete",adjust="fdr")
df_corr_Pang_bray <- data.frame("R"=corr_Pang_bray$r[intersect(rownames(SpeciesScores),rownames(corr_Pang_bray$r)),],"P"=corr_Pang_bray$p[intersect(rownames(SpeciesScores),rownames(corr_Pang_bray$r)),],"HACK.Index"=SpeciesScores[intersect(rownames(SpeciesScores),rownames(corr_Pang_bray$r)),4])
ggplot(df_corr_Pang_bray,aes(x=R,y=HACK.Index))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
df_corr_Pang_bray$dir <- ifelse(df_corr_Pang_bray[,2]<=0.15,2*sign(df_corr_Pang_bray[,1]),sign(df_corr_Pang_bray[,1]))

bray_vec_Pang_Sig_Negative <- df_corr_Pang_bray[df_corr_Pang_bray[,4]== -2,3]
bray_vec_Pang_Not_Significant <- df_corr_Pang_bray[df_corr_Pang_bray[,3]!= -2,3]

beanplot(bray_vec_Olssen_Sig_Negative,bray_vec_Olssen_Not_Significant,bray_vec_Olsson_wgs_Sig_Negative,bray_vec_Olsson_wgs_Not_Significant,bray_vec_Pang_Sig_Negative,bray_vec_Pang_Not_Significant,side="both",what=c(0,1,1,0),overallline="median",col=list("deepskyblue","gold"))

common_rows <- intersect(rownames(df_corr_Pang_bray),intersect(rownames(df_corr_Olssen_bray),rownames(df_corr_Olsson_wgs_bray)))

bray_follow_up_directionality <- data.frame("Olssen"=df_corr_Olssen_bray[common_rows,"dir"],"Olsson_wgs"=df_corr_Olsson_wgs_bray[common_rows,"dir"],"Pang"=df_corr_Pang_bray[common_rows,"dir"],"HACK.index"=SpeciesScores[common_rows,4],row.names=common_rows)
bray_follow_up_directionality$total_sig <- apply(bray_follow_up_directionality[,c(1:3)],1,function(x)(length(x[!is.na(x)&(x== -2)])))

boxplot(bray_follow_up_directionality$HACK.index~ifelse(bray_follow_up_directionality$total_sig>=2,2,bray_follow_up_directionality$total_sig),col=c("coral","bisque","darkseagreen"),outline=FALSE,range=1,xlab="",ylab="")

## Stability: Aitchison
corr_Olssen_aitchison <- corr.test(Olssen_full_data_norm,Olssen_aitchison_follow_up[rownames(Olssen_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")
df_corr_Olssen_aitchison <- data.frame("R"=corr_Olssen_aitchison$r,"P"=corr_Olssen_aitchison$p,"HACK-Index"=SpeciesScores[rownames(corr_Olssen_aitchison$r),4])
ggplot(df_corr_Olssen_aitchison,aes(x=R,y=HACK.Index))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
df_corr_Olssen_aitchison$dir <- ifelse(df_corr_Olssen_aitchison[,2]<=0.15,2*sign(df_corr_Olssen_aitchison[,1]),sign(df_corr_Olssen_aitchison[,1]))

aitchison_vec_Olssen_Sig_Negative <- df_corr_Olssen_aitchison[df_corr_Olssen_aitchison[,4]== -2,3]
aitchison_vec_Olssen_Not_Significant <- df_corr_Olssen_aitchison[df_corr_Olssen_aitchison[,3]!= -2,3]

corr_Olsson_wgs_aitchison <- corr.test(Olsson_wgs_full_data_norm,Olsson_wgs_aitchison_follow_up[rownames(Olsson_wgs_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")
df_corr_Olsson_wgs_aitchison <- data.frame("R"=corr_Olsson_wgs_aitchison$r,"P"=corr_Olsson_wgs_aitchison$p,"HACK-Index"=SpeciesScores[rownames(corr_Olsson_wgs_aitchison$r),4])
ggplot(df_corr_Olsson_wgs_aitchison,aes(x=R,y=HACK.Index))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
df_corr_Olsson_wgs_aitchison$dir <- ifelse(df_corr_Olsson_wgs_aitchison[,2]<=0.15,2*sign(df_corr_Olsson_wgs_aitchison[,1]),sign(df_corr_Olsson_wgs_aitchison[,1]))

aitchison_vec_Olsson_wgs_Sig_Negative <- df_corr_Olsson_wgs_aitchison[df_corr_Olsson_wgs_aitchison[,4]== -2,3]
aitchison_vec_Olsson_wgs_Not_Significant <- df_corr_Olsson_wgs_aitchison[df_corr_Olsson_wgs_aitchison[,3]!= -2,3]

corr_Pang_aitchison <- corr.test(Pang_full_data_norm,Pang_aitchison_follow_up[rownames(Pang_full_data_norm),],method="spearman",use="pairwise.complete",adjust="fdr")
df_corr_Pang_aitchison <- data.frame("R"=corr_Pang_aitchison$r[intersect(rownames(SpeciesScores),rownames(corr_Pang_aitchison$r)),],"P"=corr_Pang_aitchison$p[intersect(rownames(SpeciesScores),rownames(corr_Pang_aitchison$r)),],"HACK.Index"=SpeciesScores[intersect(rownames(SpeciesScores),rownames(corr_Pang_aitchison$r)),4])
ggplot(df_corr_Pang_aitchison,aes(x=R,y=HACK.Index))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
df_corr_Pang_aitchison$dir <- ifelse(df_corr_Pang_aitchison[,2]<=0.15,2*sign(df_corr_Pang_aitchison[,1]),sign(df_corr_Pang_aitchison[,1]))

aitchison_vec_Pang_Sig_Negative <- df_corr_Pang_aitchison[df_corr_Pang_aitchison[,4]== -2,3]
aitchison_vec_Pang_Not_Significant <- df_corr_Pang_aitchison[df_corr_Pang_aitchison[,3]!= -2,3]

beanplot(aitchison_vec_Olssen_Sig_Negative,aitchison_vec_Olssen_Not_Significant,aitchison_vec_Olsson_wgs_Sig_Negative,aitchison_vec_Olsson_wgs_Not_Significant,aitchison_vec_Pang_Sig_Negative,aitchison_vec_Pang_Not_Significant,side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))

common_rows <- intersect(rownames(df_corr_Pang_aitchison),intersect(rownames(df_corr_Olssen_aitchison),rownames(df_corr_Olsson_wgs_aitchison)))

aitchison_follow_up_directionality <- data.frame("Olssen"=df_corr_Olssen_aitchison[common_rows,"dir"],"Olsson_wgs"=df_corr_Olsson_wgs_aitchison[common_rows,"dir"],"Pang"=df_corr_Pang_aitchison[common_rows,"dir"],"HACK.index"=SpeciesScores[common_rows,4],row.names=common_rows)
aitchison_follow_up_directionality$total_sig <- apply(aitchison_follow_up_directionality[,c(1:3)],1,function(x)(length(x[!is.na(x)&(x== -2)])))

## Disease
Flemer_full_data_norm$status <- ifelse(Flemer_et_al_metadata[rownames(Flemer_full_data_norm),"disease"]=="control",0,1)
rlm_Flemer <- batch_rlm(Flemer_full_data_norm,intersect(rownames(SpeciesScores),colnames(Flemer_full_data_norm)),"status")
rlm_Flemer[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Flemer),4]

Song_full_data_norm <- as.data.frame(Song_full_data_norm)
Song_full_data_norm$status <- ifelse(Song_et_al_metadata[rownames(Song_full_data_norm),"group"]=="Normotension",0,1)
rlm_Song <- batch_rlm(Song_full_data_norm,intersect(rownames(SpeciesScores),colnames(Song_full_data_norm)),"status")
rlm_Song[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Song),4]

Parbie_full_data_norm <- as.data.frame(Parbie_full_data_norm)
Parbie_full_data_norm$status <- ifelse(Parbie_et_al_metadata[rownames(Parbie_full_data_norm),"disease"] == "control",0,1)
rlm_Parbie <- batch_rlm(Parbie_full_data_norm,intersect(rownames(SpeciesScores),colnames(Parbie_full_data_norm)),"status")
rlm_Parbie[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Parbie),4]

Nagata_full_data_norm <- as.data.frame(Nagata_full_data_norm)
Nagata_full_data_norm$status <- ifelse(Nagata_et_al_metadata[rownames(Nagata_mHACK),"disease_status"]=="Control",0,1)
rlm_Nagata <- batch_rlm(Nagata_full_data_norm,intersect(rownames(SpeciesScores),colnames(Nagata_full_data_norm)),"status")
rlm_Nagata[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Nagata),4]

Saleem_full_data_norm <- as.data.frame(Saleem_full_data_norm)
Saleem_full_data_norm$status <- ifelse(Saleem_et_al_metadata[rownames(Saleem_mHACK),"disease"]=="control",0,1)
rlm_Saleem <- batch_rlm(Saleem_full_data_norm,intersect(rownames(SpeciesScores),colnames(Saleem_full_data_norm)),"status")
rlm_Saleem[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Saleem),4]

MicroDiab_Denmark_full_data_norm <- as.data.frame(MicroDiab_Denmark_full_data_norm)
MicroDiab_Denmark_full_data_norm$status <- ifelse(MicroDiab_Denmark_Metadata[rownames(MicroDiab_Denmark_mHACK),"Disease"] == "Control",0,ifelse(MicroDiab_Denmark_Metadata[rownames(MicroDiab_Denmark_mHACK),"Disease"] == "Prediabetes",1,2))
rlm_MicroDiab_Denmark <- batch_rlm(MicroDiab_Denmark_full_data_norm,intersect(rownames(SpeciesScores),colnames(MicroDiab_Denmark_full_data_norm)),"status")
rlm_MicroDiab_Denmark[,"HACK.Index"] <- SpeciesScores[rownames(rlm_MicroDiab_Denmark),4]

Xu_full_data_norm <- as.data.frame(Xu_full_data_norm)
Xu_select_rows <- rownames(Xu_full_df[!is.na(Xu_full_df[,2])&(Xu_full_df[,2]>0),])
Xu_full_data_norm <- Xu_full_data_norm[Xu_select_rows,]
Xu_full_data_norm$status1 <- ifelse(rownames(Xu_full_data_norm) %in% rownames(Xu_df[((Xu_df$MMSE >= 27)&(Xu_df$Barthel_Score == 100)),]),0,1)
Xu_full_data_norm$status2 <- ifelse(Xu_df[Xu_select_rows,"Numb_Disorders"]>0,2,Xu_df[Xu_select_rows,"Numb_Disorders"])
rlm_Xu_1 <- batch_rlm(Xu_full_data_norm,intersect(rownames(SpeciesScores),colnames(Xu_full_data_norm)),"status1")
rlm_Xu_1[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Xu_1),4]
rlm_Xu_2 <- batch_rlm(Xu_full_data_norm,intersect(rownames(SpeciesScores),colnames(Xu_full_data_norm)),"status2")
rlm_Xu_2[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Xu_2),4]

Pang_full_data_norm <- as.data.frame(Pang_full_data_norm)
Pang_filt_CrossSectional_data_norm <- Pang_full_data_norm[rownames(Pang_CrossSectional_Metadata[!is.na(Pang_CrossSectional_Metadata$Health_status),]),]
Pang_filt_CrossSectional_Metadata <- Pang_CrossSectional_Metadata[rownames(Pang_CrossSectional_Metadata[!is.na(Pang_CrossSectional_Metadata$Health_status),]),]
Pang_filt_CrossSectional_data_norm$status <- ifelse(Pang_filt_CrossSectional_Metadata$Health_status == "H",0,ifelse(Pang_filt_CrossSectional_Metadata$Health_status == "LH",1,2))
rlm_Pang_CrossSectional <- batch_rlm(Pang_filt_CrossSectional_data_norm,intersect(rownames(SpeciesScores),colnames(Pang_filt_CrossSectional_data_norm)),"status")
rlm_Pang_CrossSectional[,"HACK.Index"] <- SpeciesScores[rownames(rlm_Pang_CrossSectional),4]

beanplot(rlm_Flemer[rlm_Flemer[,4]<= -2,5],rlm_Flemer[rlm_Flemer[,4]> -2,5],rlm_Song[rlm_Song[,4]<= -2,5],rlm_Song[rlm_Song[,4]>2,5],rlm_Nagata[rlm_Nagata[,4]<= -2,5],rlm_Nagata[rlm_Nagata[,4]>2,5],rlm_Parbie[rlm_Parbie[,4]<= -2,5],rlm_Parbie[rlm_Parbie[,4]>2,5],rlm_Saleem[rlm_Saleem[,4]<= -2,5],rlm_Saleem[rlm_Saleem[,4]>2,5],rlm_MicroDiab_Denmark[rlm_MicroDiab_Denmark[,4]<= -2,5],rlm_MicroDiab_Denmark[rlm_MicroDiab_Denmark[,4]>2,5],rlm_Xu_1[rlm_Xu_1[,4]<= -2,5],rlm_Xu_1[rlm_Xu_1[,4]>2,5],rlm_Xu_2[rlm_Xu_2[,4]<= -2,5],rlm_Xu_2[rlm_Xu_2[,4]>2,5],rlm_Pang_CrossSectional[rlm_Pang_CrossSectional[,4]<= -2,5],rlm_Pang_CrossSectional[rlm_Pang_CrossSectional[,4]>2,5],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))

disease_common_rows <- names(which(table(c(rownames(rlm_Flemer),rownames(rlm_Song),rownames(rlm_Nagata),rownames(rlm_Parbie),rownames(rlm_Saleem),rownames(rlm_MicroDiab_Denmark),rownames(rlm_Xu_1),rownames(rlm_Xu_2),rownames(rlm_Pang_CrossSectional)))==9))

disease_directionality <- data.frame("Flemer"=as.numeric(rlm_Flemer[disease_common_rows,4]),"Song"=as.numeric(rlm_Song[disease_common_rows,4]),"Nagata"=as.numeric(rlm_Nagata[disease_common_rows,4]),"Saleem"=as.numeric(rlm_Saleem[disease_common_rows,4]),"Parbie"=as.numeric(rlm_Parbie[disease_common_rows,4]),"MicroDiab_Denmark"=as.numeric(rlm_MicroDiab_Denmark[disease_common_rows,4]),"Xu_1"=as.numeric(rlm_Xu_1[disease_common_rows,4]),"Xu_2"=as.numeric(rlm_Xu_2[disease_common_rows,4]),"Pang_CrossSectional"=as.numeric(rlm_Pang_CrossSectional[disease_common_rows,4]),"HACK.Index"=as.numeric(SpeciesScores[disease_common_rows,4]),row.names=disease_common_rows)

disease_directionality$SigNegative <- apply(disease_directionality[,c(1:9)],1,function(x)(length(x[x<=-2])))

boxplot(disease_directionality[,10]~cut(disease_directionality[,11],breaks=c(0,1,4,9),include.lowest=TRUE))

