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
load(paste0(dir_name,"SpeciesScores_NEW.RData"))
load(paste0(dir_name,"Song_species.RData"))
load(paste0(dir_name,"Flemer_species.RData"))
load(paste0(dir_name,"Nagata_species.RData"))
load(paste0(dir_name,"Luan_species.RData"))
load(paste0(dir_name,"Saleem_species.RData"))
load(paste0(dir_name,"MicroDiab_Denmark_species.RData"))
load(paste0(dir_name,"Parbie_species.RData"))
load(paste0(dir_name,"Hernandez_species.RData"))
load(paste0(dir_name,"Bajaj_species.RData"))
load(paste0(dir_name,"Wallen_species.RData"))
load(paste0(dir_name,"Validation_3R_Analysis.RData"))

TopHACKs <- rownames(SpeciesScores_NEW[SpeciesScores_NEW[,4]>=0.67,])
MidHACKs <- rownames(SpeciesScores_NEW[(SpeciesScores_NEW[,4]<0.67)&(SpeciesScores_NEW[,4]>=0.33),])
LowHACKs <- rownames(SpeciesScores_NEW[SpeciesScores_NEW[,4]<0.33,])
HACKs <- rownames(SpeciesScores_NEW)[SpeciesScores_NEW[,4]>=0.75]


source(paste0(dir_name,"code_library.R"))
SpeciesScores_NEW <- SpeciesScores_NEW[order(SpeciesScores_NEW[,4]),]

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
			follow_up_dist[i,1] <- as.numeric(vegdist(data[c(T1,T0),],method="bray",na.rm=TRUE))
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
			follow_up_dist[i,1] <- as.numeric(vegdist(data_clr[c(T1,T0),],method="euclidean",na.rm=TRUE))
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
Olssen_data_norm <- Olssen_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Olssen_et_al_species))]
Olssen_data_norm <- Olssen_data_norm/rowSums(Olssen_data_norm)

#Olssen_mHACK <- compute_enrichment(apply(Olssen_data_norm[,intersect(colnames(Olssen_data_norm),TopHACKs)],2,rank_scale))

Olssen_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Olssen_data_norm[,intersect(TopHACKs,colnames(Olssen_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Olssen_data_norm[,intersect(MidHACKs,
colnames(Olssen_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Olssen_data_norm[,intersect(LowHACKs,colnames(Olssen_data_norm))],2,
rank_scale)))

Olsson_full_data_norm <- Olssen_data_norm

Olsson_data_norm <- Olssen_data_norm

Olssen_follow_up <- Olssen_et_al_FollowUp

Olsson_bray_follow_up <- bray_followup(Olsson_data_norm,Olssen_follow_up)
Olsson_aitchison_follow_up <- aitchison_followup(Olsson_data_norm,Olssen_follow_up)

df_Olsson_follow_up <- data.frame(bray=Olsson_bray_follow_up,aitchison=Olsson_aitchison_follow_up)

temp_corr_Olsson <- corr.test(Olsson_data_norm[rownames(df_Olsson_follow_up),],df_Olsson_follow_up,use="pairwise.complete",method="spearman")

df_Corr_Olsson <- data.frame("r_aitchison"=temp_corr_Olsson$r[,1],"r_bray"=temp_corr_Olsson$r[,2],"p_aitchison"=temp_corr_Olsson$p[,1],"p_bray"=temp_corr_Olsson$p[,2])

df_Corr_Olsson$dir_aitchison <- ifelse((df_Corr_Olsson$r_aitchison < 0)&(df_Corr_Olsson$p_aitchison <= 0.10),1,0)

df_Corr_Olsson$dir_bray <- ifelse((df_Corr_Olsson$r_bray < 0)&(df_Corr_Olsson$p_bray <= 0.10),1,0)

df_Corr_Olsson$dir_overall <- ifelse((df_Corr_Olsson$dir_bray == 0)&(df_Corr_Olsson$dir_aitchison == 0),0,1)

df_Corr_Olsson$Influence <- SpeciesScores_NEW[rownames(df_Corr_Olsson),1]
df_Corr_Olsson$Stability <- SpeciesScores_NEW[rownames(df_Corr_Olsson),2]
df_Corr_Olsson$Health <- SpeciesScores_NEW[rownames(df_Corr_Olsson),3]
df_Corr_Olsson$HACK_Score <- SpeciesScores_NEW[rownames(df_Corr_Olsson),4]

##Olssen WGS validation
Olsson_wgs_full_data_norm <- Olsson_et_al_wgs_species/rowSums(Olsson_et_al_wgs_species)
Olsson_wgs_data_norm <- Olsson_et_al_wgs_species[,intersect(rownames(SpeciesScores_NEW),colnames(Olsson_et_al_wgs_species))]

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

df_Olsson_wgs_follow_up <- data.frame(bray=Olsson_wgs_bray_follow_up,aitchison=Olsson_wgs_aitchison_follow_up)

temp_corr_Olsson_wgs <- corr.test(Olsson_wgs_data_norm[rownames(df_Olsson_wgs_follow_up),],df_Olsson_wgs_follow_up,use="pairwise.complete",method="spearman")

df_Corr_Olsson_wgs <- data.frame("r_aitchison"=temp_corr_Olsson_wgs$r[,1],"r_bray"=temp_corr_Olsson_wgs$r[,2],"p_aitchison"=temp_corr_Olsson_wgs$p[,1],"p_bray"=temp_corr_Olsson_wgs$p[,2])

df_Corr_Olsson_wgs$dir_aitchison <- ifelse((df_Corr_Olsson_wgs$r_aitchison < 0)&(df_Corr_Olsson_wgs$p_aitchison <= 0.10),1,0)

df_Corr_Olsson_wgs$dir_bray <- ifelse((df_Corr_Olsson_wgs$r_bray < 0)&(df_Corr_Olsson_wgs$p_bray <= 0.10),1,0)

df_Corr_Olsson_wgs$dir_overall <- ifelse((df_Corr_Olsson_wgs$dir_bray == 0)&(df_Corr_Olsson_wgs$dir_aitchison == 0),0,1)

df_Corr_Olsson_wgs$Influence <- SpeciesScores_NEW[rownames(df_Corr_Olsson_wgs),1]
df_Corr_Olsson_wgs$Stability <- SpeciesScores_NEW[rownames(df_Corr_Olsson_wgs),2]
df_Corr_Olsson_wgs$Health <- SpeciesScores_NEW[rownames(df_Corr_Olsson_wgs),3]
df_Corr_Olsson_wgs$HACK_Score <- SpeciesScores_NEW[rownames(df_Corr_Olsson_wgs),4]

## Pang et al Validation (Cross-Sectional)
Pang_full_data <- Pang_et_al_species
Pang_full_data_norm <- Pang_et_al_species/rowSums(Pang_et_al_species)
Pang_data <- Pang_et_al_species[,intersect(colnames(Pang_et_al_species),rownames(SpeciesScores_NEW))]
Pang_data_norm <- Pang_data/rowSums(Pang_data)

Pang_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Pang_data_norm[,intersect(TopHACKs,colnames(Pang_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Pang_data_norm[,intersect(MidHACKs,
colnames(Pang_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Pang_data_norm[,intersect(LowHACKs,colnames(Pang_data_norm))],2,
rank_scale)))

Pang_full_data_norm <- Pang_data_norm

Pang_Control <- rownames(Pang_CrossSectional_Metadata)[(!is.na(Pang_CrossSectional_Metadata$Health_status))&(Pang_CrossSectional_Metadata$Health_status == "H")]
Pang_Unhealthy_Elderly <- rownames(Pang_CrossSectional_Metadata)[(!is.na(Pang_CrossSectional_Metadata$Health_status))&(Pang_CrossSectional_Metadata$Health_status != "H")]

#Pang et al Validation (Longitudinal Data)
temp_Pang_Followup_Metadata <- Pang_Followup_Metadata[intersect(rownames(Pang_data_norm),rownames(Pang_Followup_Metadata)),]
temp_Pang_data_norm <- Pang_data_norm[intersect(rownames(Pang_data_norm),rownames(Pang_Followup_Metadata)),]

T0_rows <- rownames(temp_Pang_Followup_Metadata)[(temp_Pang_Followup_Metadata[,2] %in% rownames(Pang_data_norm))&(temp_Pang_Followup_Metadata[,1] %in% rownames(Pang_data_norm))]
T1_rows <- paste0(T0_rows,"-2")

Pang_bray_follow_up <- bray_followup(Pang_full_data_norm[c(T0_rows,T1_rows),],Pang_Followup_Metadata[c(T0_rows,T1_rows),])
Pang_aitchison_follow_up <- aitchison_followup(Pang_full_data_norm[c(T0_rows,T1_rows),],Pang_Followup_Metadata[c(T0_rows,T1_rows),])

#Xu et al Validation
Xu_data_norm <- Xu_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Xu_et_al_species))]/rowSums(Xu_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Xu_et_al_species))])

Xu_full_data_norm <- as.data.frame(Xu_data_norm)

Xu_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Xu_full_data_norm[,intersect(TopHACKs,colnames(Xu_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Xu_full_data_norm[,intersect(MidHACKs,
colnames(Xu_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Xu_full_data_norm[,intersect(LowHACKs,colnames(Xu_full_data_norm))],2,
rank_scale)))

## Song et al Validation (Cross-Sectional)
Song_data_norm <- Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores_NEW))])

Song_full_data_norm <- Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Song_et_al_species[,intersect(colnames(Song_et_al_species),rownames(SpeciesScores_NEW))])

Song_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Song_full_data_norm[,intersect(TopHACKs,colnames(Song_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Song_full_data_norm[,intersect(MidHACKs,
colnames(Song_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Song_full_data_norm[,intersect(LowHACKs,colnames(Song_full_data_norm))],2,
rank_scale)))

## Flemer et al Validation (Cross-Sectional)
Flemer_data_norm <- Flemer_et_al_species[,intersect(colnames(Flemer_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Flemer_et_al_species[,intersect(colnames(Flemer_et_al_species),rownames(SpeciesScores_NEW))])

Flemer_full_data_norm <- as.data.frame(Flemer_data_norm)

Flemer_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Flemer_full_data_norm[,intersect(TopHACKs,colnames(Flemer_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Flemer_full_data_norm[,intersect(MidHACKs,
colnames(Flemer_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Flemer_full_data_norm[,intersect(LowHACKs,colnames(Flemer_full_data_norm))],2,
rank_scale)))

## Nagata et al Validation (Cross-Sectional)
Nagata_data_norm <- Nagata_et_al_species[,intersect(colnames(Nagata_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Nagata_et_al_species[,intersect(colnames(Nagata_et_al_species),rownames(SpeciesScores_NEW))])

Nagata_full_data_norm <- as.data.frame(Nagata_data_norm)

Nagata_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Nagata_full_data_norm[,intersect(TopHACKs,colnames(Nagata_full_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Nagata_full_data_norm[,intersect(MidHACKs,
colnames(Nagata_full_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Nagata_full_data_norm[,intersect(LowHACKs,colnames(Nagata_full_data_norm))],2,
rank_scale)))

#Luan et al Validation
Luan_data_norm <- Luan_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Luan_et_al_species))]/rowSums(Luan_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Luan_et_al_species))])

Luan_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Luan_data_norm[,intersect(TopHACKs,colnames(Luan_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Luan_data_norm[,intersect(MidHACKs,
colnames(Luan_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Luan_data_norm[,intersect(LowHACKs,colnames(Luan_data_norm))],2,
rank_scale)))

Luan_full_data_norm <- as.data.frame(Luan_data_norm)

Luan_et_al_metadata$ElderlyHealthStatus <- rank_scale(Luan_et_al_metadata$ADL) - rank_scale(Luan_et_al_metadata$GDS_15)

#MicroDiab Denmark
MicroDiab_Denmark_data_norm <- MicroDiab_Denmark_Species[,intersect(rownames(SpeciesScores_NEW),colnames(MicroDiab_Denmark_Species))] / rowSums(MicroDiab_Denmark_Species[,intersect(rownames(SpeciesScores_NEW),colnames(MicroDiab_Denmark_Species))])

MicroDiab_Denmark_HACKGrouped <- data.frame("Top"=rowMeans(apply(
MicroDiab_Denmark_data_norm[,intersect(TopHACKs,colnames(MicroDiab_Denmark_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(MicroDiab_Denmark_data_norm[,intersect(MidHACKs,
colnames(MicroDiab_Denmark_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
MicroDiab_Denmark_data_norm[,intersect(LowHACKs,colnames(MicroDiab_Denmark_data_norm))],2,
rank_scale)))

MicroDiab_Denmark_full_data_norm <- as.data.frame(MicroDiab_Denmark_data_norm)

### Saleem et al Validation
Saleem_data_norm <- Saleem_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Saleem_et_al_species))]/rowSums(Saleem_et_al_species[,intersect(rownames(SpeciesScores_NEW),colnames(Saleem_et_al_species))])

Saleem_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Saleem_data_norm[,intersect(TopHACKs,colnames(Saleem_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Saleem_data_norm[,intersect(MidHACKs,
colnames(Saleem_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Saleem_data_norm[,intersect(LowHACKs,colnames(Saleem_data_norm))],2,
rank_scale)))

Saleem_full_data_norm <- Saleem_data_norm

### Parbie et al Validation
Parbie_data_norm <- Parbie_et_al_species[,intersect(colnames(Parbie_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Parbie_et_al_species[,intersect(colnames(Parbie_et_al_species),rownames(SpeciesScores_NEW))])

Parbie_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Parbie_data_norm[,intersect(TopHACKs,colnames(Parbie_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Parbie_data_norm[,intersect(MidHACKs,
colnames(Parbie_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Parbie_data_norm[,intersect(LowHACKs,colnames(Parbie_data_norm))],2,
rank_scale)))

Parbie_full_data_norm <- as.data.frame(Parbie_data_norm)

#Hernandex et al validation (Cross Sectional)
Hernandez_data_norm <- Hernandez_et_al_species[,intersect(colnames(Hernandez_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Hernandez_et_al_species[,intersect(colnames(Hernandez_et_al_species),rownames(SpeciesScores_NEW))])

Hernandez_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Hernandez_data_norm[,intersect(TopHACKs,colnames(Hernandez_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Hernandez_data_norm[,intersect(MidHACKs,
colnames(Hernandez_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Hernandez_data_norm[,intersect(LowHACKs,colnames(Hernandez_data_norm))],2,
rank_scale)))

Hernandez_full_data_norm <- as.data.frame(Hernandez_data_norm)

#Bajaj et al validation (Cross Sectional)
Bajaj_data_norm <- Bajaj_et_al_species[,intersect(colnames(Bajaj_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Bajaj_et_al_species[,intersect(colnames(Bajaj_et_al_species),rownames(SpeciesScores_NEW))])

Bajaj_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Bajaj_data_norm[,intersect(TopHACKs,colnames(Bajaj_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Bajaj_data_norm[,intersect(MidHACKs,
colnames(Bajaj_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Bajaj_data_norm[,intersect(LowHACKs,colnames(Bajaj_data_norm))],2,
rank_scale)))

Bajaj_full_data_norm <- as.data.frame(Bajaj_data_norm)

#Wallen et al validation (Cross Sectional)
Wallen_data_norm <- Wallen_et_al_species[,intersect(colnames(Wallen_et_al_species),rownames(SpeciesScores_NEW))]/rowSums(Wallen_et_al_species[,intersect(colnames(Wallen_et_al_species),rownames(SpeciesScores_NEW))])

Wallen_HACKGrouped <- data.frame("Top"=rowMeans(apply(
Wallen_data_norm[,intersect(TopHACKs,colnames(Wallen_data_norm))],2,
rank_scale)),"Mid"=rowMeans(apply(Wallen_data_norm[,intersect(MidHACKs,
colnames(Wallen_data_norm))],2,rank_scale)),"Low"=rowMeans(apply(
Wallen_data_norm[,intersect(LowHACKs,colnames(Wallen_data_norm))],2,
rank_scale)))

Wallen_full_data_norm <- as.data.frame(Wallen_data_norm)

###### Correlation with Core Influencer Score
Flemer_3R <- Flemer_3R[Flemer_3R[,2]>0,]
Flemer_3R$HACKScores <- SpeciesScores_NEW[rownames(Flemer_3R),4]
Flemer_3R$InfluenceScore <- rank_scale(Flemer_3R[,2]*Flemer_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Flemer_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Flemer_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Flemer_C <- rownames(Flemer_3R)[Flemer_3R$CoreSp == "CoreInfluencer"]
Flemer_NC <- rownames(Flemer_3R)[Flemer_3R$CoreSp == "NonCoreInfluencer"]

##
Pang_3R <- Pang_3R[Pang_3R[,2]>0,]
Pang_3R$HACKScores <- SpeciesScores_NEW[rownames(Pang_3R),4]
Pang_3R$InfluenceScore <- rank_scale(Pang_3R[,2]*Pang_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Pang_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Pang_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Pang_C <- rownames(Pang_3R)[Pang_3R$CoreSp == "CoreInfluencer"]
Pang_NC <- rownames(Pang_3R)[Pang_3R$CoreSp == "NonCoreInfluencer"]

##
Xu_3R <- Xu_3R[Xu_3R[,2]>0,]
Xu_3R$HACKScores <- SpeciesScores_NEW[rownames(Xu_3R),4]
Xu_3R$InfluenceScore <- rank_scale(Xu_3R[,2]*Xu_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Xu_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Xu_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Xu_C <- rownames(Xu_3R)[Xu_3R$CoreSp == "CoreInfluencer"]
Xu_NC <- rownames(Xu_3R)[Xu_3R$CoreSp == "NonCoreInfluencer"]

##
Olssen_3R <- Olssen_3R[Olssen_3R[,2]>0,]
Olssen_3R$HACKScores <- SpeciesScores_NEW[rownames(Olssen_3R),4]
Olssen_3R$InfluenceScore <- rank_scale(Olssen_3R[,2]*Olssen_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Olssen_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Olssen_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19,size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Olssen_C <- rownames(Olssen_3R)[Olssen_3R$CoreSp == "CoreInfluencer"]
Olssen_NC <- rownames(Olssen_3R)[Olssen_3R$CoreSp == "NonCoreInfluencer"]

##
Olsson_wgs_3R <- Olsson_wgs_3R[Olsson_wgs_3R[,2]>0,]
Olsson_wgs_3R$HACKScores <- SpeciesScores_NEW[rownames(Olsson_wgs_3R),4]
Olsson_wgs_3R$InfluenceScore <- rank_scale(Olsson_wgs_3R[,2]*Olsson_wgs_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Olsson_wgs_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Olsson_wgs_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Olsson_wgs_C <- rownames(Olsson_wgs_3R)[Olsson_wgs_3R$CoreSp == "CoreInfluencer"]
Olsson_wgs_NC <- rownames(Olsson_wgs_3R)[Olsson_wgs_3R$CoreSp == "NonCoreInfluencer"]

##
Luan_3R <- Luan_3R[Luan_3R[,2]>0,]
Luan_3R$HACKScores <- SpeciesScores_NEW[rownames(Luan_3R),4]
Luan_3R$InfluenceScore <- rank_scale(Luan_3R[,2]*Luan_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Luan_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Luan_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Luan_C <- rownames(Luan_3R)[Luan_3R$CoreSp == "CoreInfluencer"]
Luan_NC <- rownames(Luan_3R)[Luan_3R$CoreSp == "NonCoreInfluencer"]

##
Nagata_3R <- Nagata_3R[Nagata_3R[,2]>0,]
Nagata_3R$HACKScores <- SpeciesScores_NEW[rownames(Nagata_3R),4]
Nagata_3R$InfluenceScore <- rank_scale(Nagata_3R[,2]*Nagata_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Nagata_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Nagata_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Nagata_C <- rownames(Nagata_3R)[Nagata_3R$CoreSp == "CoreInfluencer"]
Nagata_NC <- rownames(Nagata_3R)[Nagata_3R$CoreSp == "NonCoreInfluencer"]

##
MicroDiab_Denmark_3R <- MicroDiab_Denmark_3R[MicroDiab_Denmark_3R[,2]>0,]
MicroDiab_Denmark_3R$HACKScores <- SpeciesScores_NEW[rownames(MicroDiab_Denmark_3R),4]
MicroDiab_Denmark_3R$InfluenceScore <- rank_scale(MicroDiab_Denmark_3R[,2]*MicroDiab_Denmark_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_MicroDiab_Denmark_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(MicroDiab_Denmark_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

MicroDiab_Denmark_C <- rownames(MicroDiab_Denmark_3R)[MicroDiab_Denmark_3R$CoreSp == "CoreInfluencer"]
MicroDiab_Denmark_NC <- rownames(MicroDiab_Denmark_3R)[MicroDiab_Denmark_3R$CoreSp == "NonCoreInfluencer"]

##
Saleem_3R <- Saleem_3R[Saleem_3R[,2]>0,]
Saleem_3R$HACKScores <- SpeciesScores_NEW[rownames(Saleem_3R),4]
Saleem_3R$InfluenceScore <- rank_scale(Saleem_3R[,2]*Saleem_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Saleem_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Saleem_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Saleem_C <- rownames(Saleem_3R)[Saleem_3R$CoreSp == "CoreInfluencer"]
Saleem_NC <- rownames(Saleem_3R)[Saleem_3R$CoreSp == "NonCoreInfluencer"]

##
Parbie_3R <- Parbie_3R[Parbie_3R[,2]>0,]
Parbie_3R$HACKScores <- SpeciesScores_NEW[rownames(Parbie_3R),4]
Parbie_3R$InfluenceScore <- rank_scale(Parbie_3R[,2]*Parbie_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Parbie_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Parbie_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Parbie_C <- rownames(Parbie_3R)[Parbie_3R$CoreSp == "CoreInfluencer"]
Parbie_NC <- rownames(Parbie_3R)[Parbie_3R$CoreSp == "NonCoreInfluencer"]

##
Song_3R <- Song_3R[Song_3R[,2]>0,]
Song_3R$HACKScores <- SpeciesScores_NEW[rownames(Song_3R),4]
Song_3R$InfluenceScore <- rank_scale(Song_3R[,2]*Song_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Song_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Song_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Song_C <- rownames(Song_3R)[Song_3R$CoreSp == "CoreInfluencer"]
Song_NC <- rownames(Song_3R)[Song_3R$CoreSp == "NonCoreInfluencer"]

##
Hernandez_3R <- Hernandez_3R[Hernandez_3R[,2]>0,]
Hernandez_3R$HACKScores <- SpeciesScores_NEW[rownames(Hernandez_3R),4]
Hernandez_3R$InfluenceScore <- rank_scale(Hernandez_3R[,2]*Hernandez_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Hernandez_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Hernandez_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19, size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Hernandez_C <- rownames(Hernandez_3R)[Hernandez_3R$CoreSp == "CoreInfluencer"]
Hernandez_NC <- rownames(Hernandez_3R)[Hernandez_3R$CoreSp == "NonCoreInfluencer"]

##
Bajaj_3R <- Bajaj_3R[Bajaj_3R[,2]>0,]
Bajaj_3R$HACKScores <- SpeciesScores_NEW[rownames(Bajaj_3R),4]
Bajaj_3R$InfluenceScore <- rank_scale(Bajaj_3R[,2]*Bajaj_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Bajaj_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Bajaj_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19,size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Bajaj_C <- rownames(Bajaj_3R)[Bajaj_3R$CoreSp == "CoreInfluencer"]
Bajaj_NC <- rownames(Bajaj_3R)[Bajaj_3R$CoreSp == "NonCoreInfluencer"]

##
Wallen_3R <- Wallen_3R[Wallen_3R[,2]>0,]
Wallen_3R$HACKScores <- SpeciesScores_NEW[rownames(Wallen_3R),4]
Wallen_3R$InfluenceScore <- rank_scale(Wallen_3R[,2]*Wallen_3R[,3])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Corr_Wallen_3R_InfScore_HackScore.pdf", height = 4.5, width = 5)
#ggplot(Wallen_3R,aes(x=HACKScores,y=InfluenceScore))+geom_point(color="grey10",pch=19,size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

Wallen_C <- rownames(Wallen_3R)[Wallen_3R$CoreSp == "CoreInfluencer"]
Wallen_NC <- rownames(Wallen_3R)[Wallen_3R$CoreSp == "NonCoreInfluencer"]

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Beanplot_SpeciesScores_NEW_HACKscore.pdf", height = 6, width = 12)
#beanplot(SpeciesScores_NEW[Flemer_C,4],SpeciesScores_NEW[Flemer_NC,4],SpeciesScores_NEW[Nagata_C,4],SpeciesScores_NEW[Nagata_NC,4],SpeciesScores_NEW[MicroDiab_Denmark_C,4],SpeciesScores_NEW[MicroDiab_Denmark_NC,4],SpeciesScores_NEW[Saleem_C,4],SpeciesScores_NEW[Saleem_NC,4],SpeciesScores_NEW[Parbie_C,4],SpeciesScores_NEW[Parbie_NC,4],SpeciesScores_NEW[Xu_C,4],SpeciesScores_NEW[Xu_NC,4],SpeciesScores_NEW[Pang_C,4],SpeciesScores_NEW[Pang_NC,4],SpeciesScores_NEW[Luan_C,4],SpeciesScores_NEW[Luan_NC,4],SpeciesScores_NEW[Olssen_C,4],SpeciesScores_NEW[Olssen_NC,4],SpeciesScores_NEW[Olsson_wgs_C,4],SpeciesScores_NEW[Olsson_wgs_NC,4],SpeciesScores_NEW[Song_C,4],SpeciesScores_NEW[Song_NC,4],SpeciesScores_NEW[Hernandez_C,4],SpeciesScores_NEW[Hernandez_NC,4],SpeciesScores_NEW[Bajaj_C,4],SpeciesScores_NEW[Bajaj_NC,4],SpeciesScores_NEW[Wallen_C,4],SpeciesScores_NEW[Wallen_NC,4],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Beanplot_SpeciesScores_NEW_InfluenceScore.pdf", height = 6, width = 12)
#beanplot(SpeciesScores_NEW[Flemer_C,1],SpeciesScores_NEW[Flemer_NC,1],SpeciesScores_NEW[Nagata_C,1],SpeciesScores_NEW[Nagata_NC,1],SpeciesScores_NEW[MicroDiab_Denmark_C,1],SpeciesScores_NEW[MicroDiab_Denmark_NC,1],SpeciesScores_NEW[Saleem_C,1],SpeciesScores_NEW[Saleem_NC,4],SpeciesScores_NEW[Parbie_C,1],SpeciesScores_NEW[Parbie_NC,1],SpeciesScores_NEW[Xu_C,1],SpeciesScores_NEW[Xu_NC,1],SpeciesScores_NEW[Pang_C,1],SpeciesScores_NEW[Pang_NC,1],SpeciesScores_NEW[Luan_C,1],SpeciesScores_NEW[Luan_NC,1],SpeciesScores_NEW[Olssen_C,1],SpeciesScores_NEW[Olssen_NC,1],SpeciesScores_NEW[Olsson_wgs_C,1],SpeciesScores_NEW[Olsson_wgs_NC,1],SpeciesScores_NEW[Song_C,1],SpeciesScores_NEW[Song_NC,1],SpeciesScores_NEW[Hernandez_C,1],SpeciesScores_NEW[Hernandez_NC,1],SpeciesScores_NEW[Bajaj_C,1],SpeciesScores_NEW[Bajaj_NC,1],SpeciesScores_NEW[Wallen_C,1],SpeciesScores_NEW[Wallen_NC,1],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Beanplot_SpeciesScores_NEW_StabilityScore.pdf", height = 6, width = 12)
#beanplot(SpeciesScores_NEW[Flemer_C,2],SpeciesScores_NEW[Flemer_NC,2],SpeciesScores_NEW[Nagata_C,2],SpeciesScores_NEW[Nagata_NC,2],SpeciesScores_NEW[MicroDiab_Denmark_C,2],SpeciesScores_NEW[MicroDiab_Denmark_NC,2],SpeciesScores_NEW[Saleem_C,2],SpeciesScores_NEW[Saleem_NC,2],SpeciesScores_NEW[Parbie_C,2],SpeciesScores_NEW[Parbie_NC,2],SpeciesScores_NEW[Xu_C,2],SpeciesScores_NEW[Xu_NC,2],SpeciesScores_NEW[Pang_C,2],SpeciesScores_NEW[Pang_NC,2],SpeciesScores_NEW[Luan_C,2],SpeciesScores_NEW[Luan_NC,2],SpeciesScores_NEW[Olssen_C,2],SpeciesScores_NEW[Olssen_NC,2],SpeciesScores_NEW[Olsson_wgs_C,2],SpeciesScores_NEW[Olsson_wgs_NC,2],SpeciesScores_NEW[Song_C,2],SpeciesScores_NEW[Song_NC,2],SpeciesScores_NEW[Hernandez_C,2],SpeciesScores_NEW[Hernandez_NC,2],SpeciesScores_NEW[Bajaj_C,2],SpeciesScores_NEW[Bajaj_NC,2],SpeciesScores_NEW[Wallen_C,2],SpeciesScores_NEW[Wallen_NC,2],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\Beanplot_SpeciesScores_NEW_HealthScore.pdf", height = 6, width = 12)
#beanplot(SpeciesScores_NEW[Flemer_C,3],SpeciesScores_NEW[Flemer_NC,3],SpeciesScores_NEW[Nagata_C,3],SpeciesScores_NEW[Nagata_NC,3],SpeciesScores_NEW[MicroDiab_Denmark_C,3],SpeciesScores_NEW[MicroDiab_Denmark_NC,3],SpeciesScores_NEW[Saleem_C,3],SpeciesScores_NEW[Saleem_NC,3],SpeciesScores_NEW[Parbie_C,3],SpeciesScores_NEW[Parbie_NC,3],SpeciesScores_NEW[Xu_C,3],SpeciesScores_NEW[Xu_NC,3],SpeciesScores_NEW[Pang_C,3],SpeciesScores_NEW[Pang_NC,3],SpeciesScores_NEW[Luan_C,3],SpeciesScores_NEW[Luan_NC,3],SpeciesScores_NEW[Olssen_C,3],SpeciesScores_NEW[Olssen_NC,3],SpeciesScores_NEW[Olsson_wgs_C,3],SpeciesScores_NEW[Olsson_wgs_NC,3],SpeciesScores_NEW[Song_C,3],SpeciesScores_NEW[Song_NC,3],SpeciesScores_NEW[Hernandez_C,3],SpeciesScores_NEW[Hernandez_NC,3],SpeciesScores_NEW[Bajaj_C,3],SpeciesScores_NEW[Bajaj_NC,3],SpeciesScores_NEW[Wallen_C,3],SpeciesScores_NEW[Wallen_NC,3],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()


AllDetectedTaxaValidation <- names(which(table(c(rownames(Flemer_3R),rownames(Nagata_3R),rownames(MicroDiab_Denmark_3R),rownames(Saleem_3R),rownames(Parbie_3R),rownames(Xu_3R),rownames(Pang_3R),rownames(Olssen_3R),rownames(Luan_3R),rownames(Olsson_wgs_3R),rownames(Song_3R),rownames(Hernandez_3R),rownames(Wallen_3R),rownames(Bajaj_3R))) == 14))

CoreInfluencerDetection <- as.data.frame(matrix(0,nrow(SpeciesScores_NEW),14))
rownames(CoreInfluencerDetection) <- rownames(SpeciesScores_NEW)
colnames(CoreInfluencerDetection) <- c("Flemer","Nagata","MicroDiab_Denmark","Saleem","Parbie","Pang","Xu","Song","Luan","Olsson_16S","Olsson_WGS","Hernandez","Bajaj","Wallen")

CoreInfluencerDetection[Flemer_C,"Flemer"] <- 1
CoreInfluencerDetection[Nagata_C,"Nagata"] <- 1
CoreInfluencerDetection[MicroDiab_Denmark_C,"MicroDiab_Denmark"] <- 1
CoreInfluencerDetection[Saleem_C,"Saleem"] <- 1
CoreInfluencerDetection[Parbie_C,"Parbie"] <- 1
CoreInfluencerDetection[Pang_C,"Pang"] <- 1
CoreInfluencerDetection[Xu_C,"Xu"] <- 1
CoreInfluencerDetection[Song_C,"Song"] <- 1
CoreInfluencerDetection[Luan_C,"Luan"] <- 1
CoreInfluencerDetection[Olssen_C,"Olsson_16S"] <- 1
CoreInfluencerDetection[Olsson_wgs_C,"Olsson_WGS"] <- 1
CoreInfluencerDetection[Hernandez_C,"Hernandez"] <- 1
CoreInfluencerDetection[Bajaj_C,"Bajaj"] <- 1
CoreInfluencerDetection[Wallen_C,"Wallen"] <- 1

df_CoreAssociation <- data.frame(Detection=rowSums(CoreInfluencerDetection),CoreAssociation=SpeciesScores_NEW[,1],StabilityAssociation=SpeciesScores_NEW[,2],HealthAssociation=SpeciesScores_NEW[,3],HACKScore=SpeciesScores_NEW[,4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\CoreInfluenceAssociation.pdf", height = 6, width = 12)
#ggplot(df_CoreAssociation,aes(x=Detection,y=CoreAssociation))+geom_point(size=2)+geom_smooth(method='lm')+ylim(0,1)+xlim(0,13)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+ylab("") 
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\StabilityAssociation.pdf", height = 6, width = 12)
#ggplot(df_CoreAssociation,aes(x=Detection,y=StabilityAssociation))+geom_point(size=2)+geom_smooth(method='lm')+ylim(0,1)+xlim(0,13)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+ylab("") 
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\HealthAssociation.pdf", height = 6, width = 12)
#ggplot(df_CoreAssociation,aes(x=Detection,y=HealthAssociation))+geom_point(size=2)+geom_smooth(method='lm')+ylim(0,1)+xlim(0,13)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+ylab("") 
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\HACKScoreAssociation.pdf", height = 6, width = 12)
#ggplot(df_CoreAssociation,aes(x=Detection,y=HACKScore))+geom_point(size=2)+geom_smooth(method='lm')+ylim(0,1)+xlim(0,13)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+ylab("") 
#dev.off()

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

# Stability
### Olsson 16S ###
corr_Olssen_bray <- corr.test(Olssen_full_data_norm,Olssen_bray_follow_up[rownames(Olssen_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")

corr_Olssen_aitchison <- corr.test(Olssen_full_data_norm,Olssen_aitchison_follow_up[rownames(Olssen_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")

df_corr_Olssen_bray <- data.frame("R"=corr_Olssen_bray$r,"P"=corr_Olssen_bray$p,"Influence"=SpeciesScores_NEW[rownames(corr_Olssen_bray$r),1],"Stability"=SpeciesScores_NEW[rownames(corr_Olssen_bray$r),2],"Health"=SpeciesScores_NEW[rownames(corr_Olssen_bray$r),3],"HACKScore"=SpeciesScores_NEW[rownames(corr_Olssen_bray$r),4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_bray_HackScore_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_bray$dir <- ifelse(df_corr_Olssen_bray[,2]<=0.15,2*sign(df_corr_Olssen_bray[,1]),sign(df_corr_Olssen_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_bray_Influence_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_bray$dir <- ifelse(df_corr_Olssen_bray[,2]<=0.15,2*sign(df_corr_Olssen_bray[,1]),sign(df_corr_Olssen_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_bray_Stability_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_bray$dir <- ifelse(df_corr_Olssen_bray[,2]<=0.15,2*sign(df_corr_Olssen_bray[,1]),sign(df_corr_Olssen_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_bray_Health_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_bray$dir <- ifelse(df_corr_Olssen_bray[,2]<=0.15,2*sign(df_corr_Olssen_bray[,1]),sign(df_corr_Olssen_bray[,1]))

df_corr_Olssen_aitchison <- data.frame("R"=corr_Olssen_aitchison$r,"P"=corr_Olssen_aitchison$p,"Influence"=SpeciesScores_NEW[rownames(corr_Olssen_aitchison$r),1],"Stability"=SpeciesScores_NEW[rownames(corr_Olssen_aitchison$r),2],"Health"=SpeciesScores_NEW[rownames(corr_Olssen_aitchison$r),3],"HACKScore"=SpeciesScores_NEW[rownames(corr_Olssen_aitchison$r),4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_aitchison$dir <- ifelse(df_corr_Olssen_aitchison[,2]<=0.15,2*sign(df_corr_Olssen_aitchison[,1]),sign(df_corr_Olssen_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_aitchison_Influence_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_aitchison$dir <- ifelse(df_corr_Olssen_aitchison[,2]<=0.15,2*sign(df_corr_Olssen_aitchison[,1]),sign(df_corr_Olssen_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_aitchison_Stability_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_aitchison$dir <- ifelse(df_corr_Olssen_aitchison[,2]<=0.15,2*sign(df_corr_Olssen_aitchison[,1]),sign(df_corr_Olssen_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olssen_aitchison_Health_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olssen_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olssen_aitchison$dir <- ifelse(df_corr_Olssen_aitchison[,2]<=0.15,2*sign(df_corr_Olssen_aitchison[,1]),sign(df_corr_Olssen_aitchison[,1]))

sig_negative_Olsson <- union(rownames(df_corr_Olssen_aitchison[df_corr_Olssen_aitchison[,7]==-2,]),rownames(df_corr_Olssen_bray[df_corr_Olssen_bray[,7]==-2,]))

not_negative_Olsson <- setdiff(union(rownames(df_corr_Olssen_aitchison),rownames(df_corr_Olssen_bray)),sig_negative_Olsson)

### Olsson WGS ###
corr_Olsson_wgs_bray <- corr.test(Olsson_wgs_full_data_norm,Olsson_wgs_bray_follow_up[rownames(Olsson_wgs_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")

corr_Olsson_wgs_aitchison <- corr.test(Olsson_wgs_full_data_norm,Olsson_wgs_aitchison_follow_up[rownames(Olsson_wgs_full_data_norm),],use="pairwise.complete",method="spearman",adjust="fdr")

df_corr_Olsson_wgs_bray <- data.frame("R"=corr_Olsson_wgs_bray$r,"P"=corr_Olsson_wgs_bray$p,"Influence"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_bray$r),1],"Stability"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_bray$r),2],"Health"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_bray$r),3],"HACKScore"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_bray$r),4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_bray_HACKScore_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_bray,aes(x=R,y=HACKScore))+geom_point(size =2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_bray$dir <- ifelse(df_corr_Olsson_wgs_bray[,2]<=0.15,2*sign(df_corr_Olsson_wgs_bray[,1]),sign(df_corr_Olsson_wgs_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_bray_Influence_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_bray$dir <- ifelse(df_corr_Olsson_wgs_bray[,2]<=0.15,2*sign(df_corr_Olsson_wgs_bray[,1]),sign(df_corr_Olsson_wgs_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_bray_Stability_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_bray$dir <- ifelse(df_corr_Olsson_wgs_bray[,2]<=0.15,2*sign(df_corr_Olsson_wgs_bray[,1]),sign(df_corr_Olsson_wgs_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_bray_Health_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_bray$dir <- ifelse(df_corr_Olsson_wgs_bray[,2]<=0.15,2*sign(df_corr_Olsson_wgs_bray[,1]),sign(df_corr_Olsson_wgs_bray[,1]))

df_corr_Olsson_wgs_aitchison <- data.frame("R"=corr_Olsson_wgs_aitchison$r,"P"=corr_Olsson_wgs_aitchison$p,"Influence"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_aitchison$r),1],"Stability"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_aitchison$r),2],"Health"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_aitchison$r),3],"HACKScore"=SpeciesScores_NEW[rownames(corr_Olsson_wgs_aitchison$r),4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_aitchison$dir <- ifelse(df_corr_Olsson_wgs_aitchison[,2]<=0.15,2*sign(df_corr_Olsson_wgs_aitchison[,1]),sign(df_corr_Olsson_wgs_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_aitchison_Influence_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_aitchison$dir <- ifelse(df_corr_Olsson_wgs_aitchison[,2]<=0.15,2*sign(df_corr_Olsson_wgs_aitchison[,1]),sign(df_corr_Olsson_wgs_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_aitchison_Stability_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_aitchison$dir <- ifelse(df_corr_Olsson_wgs_aitchison[,2]<=0.15,2*sign(df_corr_Olsson_wgs_aitchison[,1]),sign(df_corr_Olsson_wgs_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Olsson_wgs_aitchison_Health_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Olsson_wgs_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Olsson_wgs_aitchison$dir <- ifelse(df_corr_Olsson_wgs_aitchison[,2]<=0.15,2*sign(df_corr_Olsson_wgs_aitchison[,1]),sign(df_corr_Olsson_wgs_aitchison[,1]))

sig_negative_Olsson_wgs <- union(rownames(df_corr_Olsson_wgs_aitchison[df_corr_Olsson_wgs_aitchison[,7]==-2,]),rownames(df_corr_Olsson_wgs_bray[df_corr_Olsson_wgs_bray[,7]==-2,]))

not_negative_Olsson_wgs <- setdiff(union(rownames(df_corr_Olsson_wgs_aitchison),rownames(df_corr_Olsson_wgs_bray)),sig_negative_Olsson_wgs)

### Pang ###
corr_Pang_bray <- corr.test(Pang_full_data_norm,Pang_bray_follow_up[rownames(Pang_full_data_norm),],method="spearman",use="pairwise.complete",adjust="fdr")

corr_Pang_aitchison <- corr.test(Pang_full_data_norm,Pang_aitchison_follow_up[rownames(Pang_full_data_norm),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_Pang_bray <- data.frame("R"=corr_Pang_bray$r[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_bray$r)),],"P"=corr_Pang_bray$p[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_bray$r)),],"Influence"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_bray$r)),1],"Stability"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_bray$r)),2],"Health"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_bray$r)),3],"HACKScore"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_bray$r)),4])

df_corr_Pang_aitchison <- data.frame("R"=corr_Pang_aitchison$r[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_aitchison$r)),],"P"=corr_Pang_aitchison$p[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_aitchison$r)),],"Influence"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_aitchison$r)),1],"Stability"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_aitchison$r)),2],"Health"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_aitchison$r)),3],"HACKScore"=SpeciesScores_NEW[intersect(rownames(SpeciesScores_NEW),rownames(corr_Pang_aitchison$r)),4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_bray_HACKScore_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_bray$dir <- ifelse(df_corr_Pang_bray[,2]<=0.15,2*sign(df_corr_Pang_bray[,1]),sign(df_corr_Pang_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_bray_Influence_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_bray$dir <- ifelse(df_corr_Pang_bray[,2]<=0.15,2*sign(df_corr_Pang_bray[,1]),sign(df_corr_Pang_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_bray_Stability_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_bray$dir <- ifelse(df_corr_Pang_bray[,2]<=0.15,2*sign(df_corr_Pang_bray[,1]),sign(df_corr_Pang_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_bray_Health_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_bray$dir <- ifelse(df_corr_Pang_bray[,2]<=0.15,2*sign(df_corr_Pang_bray[,1]),sign(df_corr_Pang_bray[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_aitchison$dir <- ifelse(df_corr_Pang_aitchison[,2]<=0.15,2*sign(df_corr_Pang_aitchison[,1]),sign(df_corr_Pang_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_aitchison_Influence_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_aitchison$dir <- ifelse(df_corr_Pang_aitchison[,2]<=0.15,2*sign(df_corr_Pang_aitchison[,1]),sign(df_corr_Pang_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_aitchison_Stability_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_aitchison$dir <- ifelse(df_corr_Pang_aitchison[,2]<=0.15,2*sign(df_corr_Pang_aitchison[,1]),sign(df_corr_Pang_aitchison[,1]))

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\corr_Pang_aitchison_Health_R.pdf", height = 4.5, width = 5)
#ggplot(df_corr_Pang_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
#dev.off()

df_corr_Pang_aitchison$dir <- ifelse(df_corr_Pang_aitchison[,2]<=0.15,2*sign(df_corr_Pang_aitchison[,1]),sign(df_corr_Pang_aitchison[,1]))

sig_negative_Pang <- union(rownames(df_corr_Pang_aitchison[!is.na(df_corr_Pang_aitchison[,7])&(df_corr_Pang_aitchison[,7]==-2),]),rownames(df_corr_Pang_bray[!is.na(df_corr_Pang_bray[,7])&(df_corr_Pang_bray[,7]==-2),]))

not_negative_Pang <- setdiff(union(rownames(df_corr_Pang_aitchison),rownames(df_corr_Pang_bray)),sig_negative_Pang)

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\beanplot_SpeciesScores_NEW_HACKScore_Olsson_Pang_OlssonWGS.pdf", height = 6, width = 5)
#beanplot(SpeciesScores_NEW[sig_negative_Olsson,4],SpeciesScores_NEW[not_negative_Olsson,4],SpeciesScores_NEW[sig_negative_Olsson_wgs,4],SpeciesScores_NEW[not_negative_Olsson_wgs,4],SpeciesScores_NEW[sig_negative_Pang,4],SpeciesScores_NEW[not_negative_Pang,4],side="both",what=c(0,1,1,0),overallline="median",col=list("deepskyblue","gold"))
#dev.off()

StabilityAssociationStudySpecific <- as.data.frame(matrix(0,nrow(SpeciesScores_NEW),3))
rownames(StabilityAssociationStudySpecific) <- rownames(SpeciesScores_NEW)
colnames(StabilityAssociationStudySpecific) <- c("Olsson","Olsson_WGS","Pang")

StabilityAssociationStudySpecific[sig_negative_Olsson,"Olsson"] <- 1
StabilityAssociationStudySpecific[sig_negative_Olsson_wgs,"Olsson_WGS"] <- 1
StabilityAssociationStudySpecific[sig_negative_Pang,"Pang"] <- 1

#df_StabilityDetection <- data.frame("Detection"=ifelse(rowSums(StabilityAssociationStudySpecific)>=2,2,rowSums(StabilityAssociationStudySpecific)),"Influence"=SpeciesScores_NEW[,1],"Stability"=SpeciesScores_NEW[,2],"Health"=SpeciesScores_NEW[,3],"HACKScore"=SpeciesScores_NEW[,4])

df_StabilityDetection <- data.frame("Detection"=rowSums(StabilityAssociationStudySpecific),"Influence"=SpeciesScores_NEW[,1],"Stability"=SpeciesScores_NEW[,2],"Health"=SpeciesScores_NEW[,3],"HACKScore"=SpeciesScores_NEW[,4])

## Health Association
## Flemer
Flemer_CRC <- rownames(Flemer_et_al_metadata[Flemer_et_al_metadata[,1]=="CRC",])
Flemer_Control <- rownames(Flemer_et_al_metadata[Flemer_et_al_metadata[,1]=="control",])
wilcox_Flemer <- wilcox_batch(t(Flemer_full_data_norm[Flemer_CRC,intersect(rownames(SpeciesScores_NEW),colnames(Flemer_full_data_norm))]),t(Flemer_full_data_norm[Flemer_Control,intersect(rownames(SpeciesScores_NEW),colnames(Flemer_full_data_norm))]))
df_Compare_Flemer <- data.frame("Direction"=ifelse(wilcox_Flemer[,1]<=0.10,sign(wilcox_Flemer[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Flemer),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Flemer),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Flemer),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Flemer),4])

##Song
Song_full_data_norm <- as.data.frame(Song_full_data_norm)
Song_Control <- rownames(Song_et_al_metadata[Song_et_al_metadata[,1]=="Normotension",])
Song_HT <- rownames(Song_et_al_metadata[Song_et_al_metadata[,1]!="Normotension",])
wilcox_Song <- wilcox_batch(t(Song_full_data_norm[Song_HT,intersect(rownames(SpeciesScores_NEW),colnames(Song_full_data_norm))]),t(Song_full_data_norm[Song_Control,intersect(rownames(SpeciesScores_NEW),colnames(Song_full_data_norm))]))
df_Compare_Song <- data.frame("Direction"=ifelse(wilcox_Song[,1]<=0.10,sign(wilcox_Song[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Song),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Song),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Song),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Song),4])

##Parbie
Parbie_full_data_norm <- as.data.frame(Parbie_full_data_norm)
Parbie_full_data_norm <- Parbie_full_data_norm/rowSums(Parbie_full_data_norm)
Parbie_Control <- rownames(Parbie_et_al_metadata[Parbie_et_al_metadata$disease=="control",])
Parbie_HIV <- rownames(Parbie_et_al_metadata[Parbie_et_al_metadata$disease!="control",])
wilcox_Parbie <- wilcox_batch(t(Parbie_full_data_norm[Parbie_HIV,intersect(rownames(SpeciesScores_NEW),colnames(Parbie_full_data_norm))]),t(Parbie_full_data_norm[Parbie_Control,intersect(rownames(SpeciesScores_NEW),colnames(Parbie_full_data_norm))]))
df_Compare_Parbie <- data.frame("Direction"=ifelse(wilcox_Parbie[,1]<=0.10,sign(wilcox_Parbie[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Parbie),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Parbie),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Parbie),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Parbie),4])

## Nagata
Nagata_full_data_norm <- as.data.frame(Nagata_full_data_norm)
Nagata_full_data_norm <- Nagata_full_data_norm/rowSums(Nagata_full_data_norm)
Nagata_Control <- rownames(Nagata_et_al_metadata[Nagata_et_al_metadata[,1]=="Control",])
Nagata_PDAC <- rownames(Nagata_et_al_metadata[Nagata_et_al_metadata[,1]=="PDAC",])
wilcox_Nagata <- wilcox_batch(t(Nagata_full_data_norm[Nagata_PDAC,intersect(rownames(SpeciesScores_NEW),colnames(Nagata_full_data_norm))]),t(Nagata_full_data_norm[Nagata_Control,intersect(rownames(SpeciesScores_NEW),colnames(Nagata_full_data_norm))]))
df_Compare_Nagata <- data.frame("Direction"=ifelse(wilcox_Nagata[,1]<=0.10,sign(wilcox_Nagata[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Nagata),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Nagata),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Nagata),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Nagata),4])

## Saleem
Saleem_full_data_norm <- as.data.frame(Saleem_full_data_norm)
Saleem_full_data_norm <- Saleem_full_data_norm/rowSums(Saleem_full_data_norm)
Saleem_Control <- rownames(Saleem_et_al_metadata[Saleem_et_al_metadata[,3]=="control",])
Saleem_T2D <- rownames(Saleem_et_al_metadata[Saleem_et_al_metadata[,3]!="control",])
wilcox_Saleem <- wilcox_batch(t(Saleem_full_data_norm[Saleem_T2D,intersect(rownames(SpeciesScores_NEW),colnames(Saleem_full_data_norm))]),t(Saleem_full_data_norm[Saleem_Control,intersect(rownames(SpeciesScores_NEW),colnames(Saleem_full_data_norm))]))
df_Compare_Saleem <- data.frame("Direction"=ifelse(wilcox_Saleem[,1]<=0.10,sign(wilcox_Saleem[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Saleem),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Saleem),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Saleem),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Saleem),4])

## MicroDiab Denmark
MicroDiab_Denmark_full_data_norm <- as.data.frame(MicroDiab_Denmark_full_data_norm)
MicroDiab_Denmark_full_data_norm <- MicroDiab_Denmark_full_data_norm/rowSums(MicroDiab_Denmark_full_data_norm)
MicroDiab_Denmark_Control <- rownames(MicroDiab_Denmark_Metadata[MicroDiab_Denmark_Metadata[,3] == "Control",])
MicroDiab_Denmark_T2D_Variants <- rownames(MicroDiab_Denmark_Metadata[MicroDiab_Denmark_Metadata[,3] != "Control",])
wilcox_MicroDiab_Denmark <- wilcox_batch(t(MicroDiab_Denmark_full_data_norm[MicroDiab_Denmark_T2D_Variants,intersect(rownames(SpeciesScores_NEW),colnames(MicroDiab_Denmark_full_data_norm))]),t(MicroDiab_Denmark_full_data_norm[MicroDiab_Denmark_Control,intersect(rownames(SpeciesScores_NEW),colnames(MicroDiab_Denmark_full_data_norm))]))
df_Compare_MicroDiab_Denmark <- data.frame("Direction"=ifelse(wilcox_MicroDiab_Denmark[,1]<=0.10,sign(wilcox_MicroDiab_Denmark[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_MicroDiab_Denmark),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_MicroDiab_Denmark),2],"Health"=SpeciesScores_NEW[rownames(wilcox_MicroDiab_Denmark),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_MicroDiab_Denmark),4])

## Pang
Pang_full_data_norm <- Pang_full_data_norm[,intersect(rownames(SpeciesScores_NEW),colnames(Pang_full_data_norm))]
Pang_full_data_norm <- Pang_full_data_norm/rowSums(Pang_full_data_norm)
Pang_Control <- rownames(Pang_CrossSectional_Metadata)[(!is.na(Pang_CrossSectional_Metadata$Health_status))&(Pang_CrossSectional_Metadata$Health_status == "H")]
Pang_Unhealthy_Elderly <- rownames(Pang_CrossSectional_Metadata)[(!is.na(Pang_CrossSectional_Metadata$Health_status))&(Pang_CrossSectional_Metadata$Health_status != "H")]
wilcox_Pang <- wilcox_batch(t(Pang_full_data_norm[Pang_Unhealthy_Elderly,intersect(rownames(SpeciesScores_NEW),colnames(Pang_full_data_norm))]),t(Pang_full_data_norm[Pang_Control,intersect(rownames(SpeciesScores_NEW),colnames(Pang_full_data_norm))]))
df_Compare_Pang <- data.frame("Direction"=ifelse(wilcox_Pang[,1]<=0.10,sign(wilcox_Pang[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Pang),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Pang),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Pang),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Pang),4])

## Xu
Xu_full_data_norm <- Xu_full_data_norm[,intersect(rownames(SpeciesScores_NEW),colnames(Xu_full_data_norm))]
Xu_full_data_norm <- Xu_full_data_norm/rowSums(Xu_full_data_norm)
Xu_Control <- grep("NA",rownames(Xu_et_al_metadata[((Xu_et_al_metadata$MMSE >= 27)&(Xu_et_al_metadata$MMSE >= 27)),]),value=TRUE,invert=TRUE)
Xu_Unhealthy_Elderly <- grep("NA",rownames(Xu_et_al_metadata[!((Xu_et_al_metadata$MMSE >= 27)&(Xu_et_al_metadata$Barthel_Score == 100)),]),value=TRUE,invert=TRUE)
wilcox_Xu <- wilcox_batch(t(Xu_full_data_norm[Xu_Unhealthy_Elderly,intersect(rownames(SpeciesScores_NEW),colnames(Xu_full_data_norm))]),t(Xu_full_data_norm[Xu_Control,intersect(rownames(SpeciesScores_NEW),colnames(Xu_full_data_norm))]))
df_Compare_Xu <- data.frame("Direction"=ifelse(wilcox_Xu[,1]<=0.10,sign(wilcox_Xu[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Xu),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Xu),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Xu),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Xu),4])

## Hernandez
Hernandez_full_data_norm <- Hernandez_full_data_norm[,intersect(rownames(SpeciesScores_NEW),colnames(Hernandez_full_data_norm))]
Hernandez_full_data_norm <- Hernandez_full_data_norm/rowSums(Hernandez_full_data_norm)
Hernandez_Control <- rownames(Hernandez_et_al_metadata[Hernandez_et_al_metadata$diabetes_status != "T2D",])
Hernandez_T2D <- rownames(Hernandez_et_al_metadata[Hernandez_et_al_metadata$diabetes_status == "T2D",])
wilcox_Hernandez <- wilcox_batch(t(Hernandez_full_data_norm[Hernandez_T2D,intersect(rownames(SpeciesScores_NEW),colnames(Hernandez_full_data_norm))]),t(Hernandez_full_data_norm[Hernandez_Control,intersect(rownames(SpeciesScores_NEW),colnames(Hernandez_full_data_norm))]))
df_Compare_Hernandez <- data.frame("Direction"=ifelse(wilcox_Hernandez[,1]<=0.10,sign(wilcox_Hernandez[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Hernandez),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Hernandez),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Hernandez),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Hernandez),4])

## Bajaj
Bajaj_full_data_norm <- Bajaj_full_data_norm/rowSums(Bajaj_full_data_norm)
Bajaj_Control <- rownames(Bajaj_et_al_metadata)[Bajaj_et_al_metadata[,1]=="Controls"]
Bajaj_CD_ITB <- rownames(Bajaj_et_al_metadata)[Bajaj_et_al_metadata[,1]=="ITB"]
wilcox_Bajaj <- wilcox_batch(t(Bajaj_full_data_norm[Bajaj_CD_ITB,intersect(rownames(SpeciesScores_NEW),colnames(Bajaj_full_data_norm))]),t(Bajaj_full_data_norm[Bajaj_Control,intersect(rownames(SpeciesScores_NEW),colnames(Bajaj_full_data_norm))]))
df_Compare_Bajaj <- data.frame("Direction"=ifelse(wilcox_Bajaj[,1]<=0.10,sign(wilcox_Bajaj[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Bajaj),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Bajaj),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Bajaj),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Bajaj),4])

## Wallen
Wallen_full_data_norm <- Wallen_full_data_norm[,intersect(rownames(SpeciesScores_NEW),colnames(Wallen_full_data_norm))]
Wallen_full_data_norm <- Wallen_full_data_norm/rowSums(Wallen_full_data_norm)
Wallen_Parkinsons <- rownames(Wallen_et_al_metadata[Wallen_et_al_metadata$health_state=="Parkinsons",])
Wallen_Control <- rownames(Wallen_et_al_metadata[Wallen_et_al_metadata$health_state!="Parkinsons",])
wilcox_Wallen <- wilcox_batch(t(Wallen_full_data_norm[Wallen_Parkinsons,intersect(rownames(SpeciesScores_NEW),colnames(Wallen_full_data_norm))]),t(Wallen_full_data_norm[Wallen_Control,intersect(rownames(SpeciesScores_NEW),colnames(Wallen_full_data_norm))]))
df_Compare_Wallen <- data.frame("Direction"=ifelse(wilcox_Wallen[,1]<=0.10,sign(wilcox_Wallen[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Wallen),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Wallen),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Wallen),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Wallen),4])

## Luan
Luan_full_data_norm <- Luan_full_data_norm[,intersect(rownames(SpeciesScores_NEW),colnames(Luan_full_data_norm))]
Luan_full_data_norm <- Luan_full_data_norm/rowSums(Luan_full_data_norm)
Luan_GDS <- rownames(Luan_et_al_metadata)[Luan_et_al_metadata$GDS_15>=8]
Luan_Control <- rownames(Luan_et_al_metadata)[Luan_et_al_metadata$GDS_15<8]
wilcox_Luan <- wilcox_batch(t(Luan_full_data_norm[Luan_GDS,intersect(rownames(SpeciesScores_NEW),colnames(Luan_full_data_norm))]),t(Luan_full_data_norm[Luan_Control,intersect(rownames(SpeciesScores_NEW),colnames(Luan_full_data_norm))]))
df_Compare_Luan <- data.frame("Direction"=ifelse(wilcox_Luan[,1]<=0.10,sign(wilcox_Luan[,2]),0),"Influence"=SpeciesScores_NEW[rownames(wilcox_Luan),1],"Stability"=SpeciesScores_NEW[rownames(wilcox_Luan),2],"Health"=SpeciesScores_NEW[rownames(wilcox_Luan),3],"HACKScore"=SpeciesScores_NEW[rownames(wilcox_Luan),4])

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\beanplot_Compared_directionality_NegOne_PosFour.pdf", height = 6, width = 12)
#beanplot(df_Compare_Flemer[df_Compare_Flemer[,1]==-1,4],df_Compare_Flemer[df_Compare_Flemer[,1]!=-1,4],df_Compare_Song[df_Compare_Song[,1]==-1,4],df_Compare_Song[df_Compare_Song[,1]!=-1,4],df_Compare_Parbie[df_Compare_Parbie[,1]==-1,4],df_Compare_Parbie[df_Compare_Parbie[,1]!=-1,4],df_Compare_Nagata[df_Compare_Nagata[,1]==-1,4],df_Compare_Nagata[df_Compare_Nagata[,1]!=-1,4],df_Compare_Saleem[df_Compare_Saleem[,1]==-1,4],df_Compare_Saleem[df_Compare_Saleem[,1]!=-1,4],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]==-1,4],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]!=-1,4],df_Compare_Pang[df_Compare_Pang[,1]==-1,4],df_Compare_Pang[df_Compare_Pang[,1]!=-1,4],df_Compare_Xu[df_Compare_Xu[,1]==-1,4],df_Compare_Xu[df_Compare_Xu[,1]!=-1,4],df_Compare_Hernandez[df_Compare_Hernandez[,1]==-1,4],df_Compare_Hernandez[df_Compare_Hernandez[,1]!=-1,4],df_Compare_Bajaj[df_Compare_Bajaj[,1]==-1,4],df_Compare_Bajaj[df_Compare_Bajaj[,1]!=-1,4],df_Compare_Wallen[df_Compare_Wallen[,1]==-1,4],df_Compare_Wallen[df_Compare_Wallen[,1]!=-1,4],df_Compare_Luan[df_Compare_Luan[,1]==-1,4],df_Compare_Luan[df_Compare_Luan[,1]!=-1,4],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\beanplot_Compared_directionality_NegOne_PosFive.pdf", height = 6, width = 12)
#beanplot(df_Compare_Flemer[df_Compare_Flemer[,1]==-1,5],df_Compare_Flemer[df_Compare_Flemer[,1]!=-1,5],df_Compare_Song[df_Compare_Song[,1]==-1,5],df_Compare_Song[df_Compare_Song[,1]!=-1,5],df_Compare_Parbie[df_Compare_Parbie[,1]==-1,5],df_Compare_Parbie[df_Compare_Parbie[,1]!=-1,5],df_Compare_Nagata[df_Compare_Nagata[,1]==-1,5],df_Compare_Nagata[df_Compare_Nagata[,1]!=-1,5],df_Compare_Saleem[df_Compare_Saleem[,1]==-1,5],df_Compare_Saleem[df_Compare_Saleem[,1]!=-1,5],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]==-1,5],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]!=-1,5],df_Compare_Pang[df_Compare_Pang[,1]==-1,5],df_Compare_Pang[df_Compare_Pang[,1]!=-1,5],df_Compare_Xu[df_Compare_Xu[,1]==-1,5],df_Compare_Xu[df_Compare_Xu[,1]!=-1,5],df_Compare_Hernandez[df_Compare_Hernandez[,1]==-1,5],df_Compare_Hernandez[df_Compare_Hernandez[,1]!=-1,5],df_Compare_Bajaj[df_Compare_Bajaj[,1]==-1,5],df_Compare_Bajaj[df_Compare_Bajaj[,1]!=-1,5],df_Compare_Wallen[df_Compare_Wallen[,1]==-1,5],df_Compare_Wallen[df_Compare_Wallen[,1]!=-1,5],df_Compare_Luan[df_Compare_Luan[,1]==-1,5],df_Compare_Luan[df_Compare_Luan[,1]!=-1,5],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\beanplot_Compared_directionality_NegOne_PosTwo.pdf", height = 6, width = 12)
#beanplot(df_Compare_Flemer[df_Compare_Flemer[,1]==-1,2],df_Compare_Flemer[df_Compare_Flemer[,1]!=-1,2],df_Compare_Song[df_Compare_Song[,1]==-1,2],df_Compare_Song[df_Compare_Song[,1]!=-1,2],df_Compare_Parbie[df_Compare_Parbie[,1]==-1,2],df_Compare_Parbie[df_Compare_Parbie[,1]!=-1,2],df_Compare_Nagata[df_Compare_Nagata[,1]==-1,2],df_Compare_Nagata[df_Compare_Nagata[,1]!=-1,2],df_Compare_Saleem[df_Compare_Saleem[,1]==-1,2],df_Compare_Saleem[df_Compare_Saleem[,1]!=-1,2],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]==-1,2],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]!=-1,2],df_Compare_Pang[df_Compare_Pang[,1]==-1,2],df_Compare_Pang[df_Compare_Pang[,1]!=-1,2],df_Compare_Xu[df_Compare_Xu[,1]==-1,2],df_Compare_Xu[df_Compare_Xu[,1]!=-1,2],df_Compare_Hernandez[df_Compare_Hernandez[,1]==-1,2],df_Compare_Hernandez[df_Compare_Hernandez[,1]!=-1,2],df_Compare_Bajaj[df_Compare_Bajaj[,1]==-1,2],df_Compare_Bajaj[df_Compare_Bajaj[,1]!=-1,2],df_Compare_Wallen[df_Compare_Wallen[,1]==-1,2],df_Compare_Wallen[df_Compare_Wallen[,1]!=-1,2],df_Compare_Luan[df_Compare_Luan[,1]==-1,2],df_Compare_Luan[df_Compare_Luan[,1]!=-1,2],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

#pdf("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\Revised_additional_validation\\beanplot_Compared_directionality_NegOne_PosThree.pdf", height = 6, width = 12)
#beanplot(df_Compare_Flemer[df_Compare_Flemer[,1]==-1,3],df_Compare_Flemer[df_Compare_Flemer[,1]!=-1,3],df_Compare_Song[df_Compare_Song[,1]==-1,3],df_Compare_Song[df_Compare_Song[,1]!=-1,3],df_Compare_Parbie[df_Compare_Parbie[,1]==-1,3],df_Compare_Parbie[df_Compare_Parbie[,1]!=-1,3],df_Compare_Nagata[df_Compare_Nagata[,1]==-1,3],df_Compare_Nagata[df_Compare_Nagata[,1]!=-1,3],df_Compare_Saleem[df_Compare_Saleem[,1]==-1,3],df_Compare_Saleem[df_Compare_Saleem[,1]!=-1,3],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]==-1,3],df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]!=-1,3],df_Compare_Pang[df_Compare_Pang[,1]==-1,3],df_Compare_Pang[df_Compare_Pang[,1]!=-1,3],df_Compare_Xu[df_Compare_Xu[,1]==-1,3],df_Compare_Xu[df_Compare_Xu[,1]!=-1,3],df_Compare_Hernandez[df_Compare_Hernandez[,1]==-1,3],df_Compare_Hernandez[df_Compare_Hernandez[,1]!=-1,3],df_Compare_Bajaj[df_Compare_Bajaj[,1]==-1,3],df_Compare_Bajaj[df_Compare_Bajaj[,1]!=-1,3],df_Compare_Wallen[df_Compare_Wallen[,1]==-1,3],df_Compare_Wallen[df_Compare_Wallen[,1]!=-1,3],df_Compare_Luan[df_Compare_Luan[,1]==-1,3],df_Compare_Luan[df_Compare_Luan[,1]!=-1,3],side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"))
#dev.off()

DirectionalityHealthAssociation <- as.data.frame(matrix(0,nrow(SpeciesScores_NEW),12))
rownames(DirectionalityHealthAssociation) <- rownames(SpeciesScores_NEW)
colnames(DirectionalityHealthAssociation) <- c("Flemer","Song","Parbie","Nagata","Saleem","MicroDiab_Denmark","Pang","Xu","Hernandez","Bajaj","Wallen","Luan")

DirectionalityHealthAssociation[rownames(df_Compare_Flemer[df_Compare_Flemer[,1]==-1,]),"Flemer"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Song[df_Compare_Song[,1]==-1,]),"Song"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Parbie[df_Compare_Parbie[,1]==-1,]),"Parbie"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Nagata[df_Compare_Nagata[,1]==-1,]),"Nagata"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Saleem[df_Compare_Saleem[,1]==-1,]),"Saleem"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_MicroDiab_Denmark[df_Compare_MicroDiab_Denmark[,1]==-1,]),"MicroDiab_Denmark"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Pang[df_Compare_Pang[,1]==-1,]),"Pang"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Xu[df_Compare_Xu[,1]==-1,]),"Xu"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Hernandez[df_Compare_Hernandez[,1]==-1,]),"Hernandez"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Bajaj[df_Compare_Bajaj[,1]==-1,]),"Bajaj"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Wallen[df_Compare_Wallen[,1]==-1,]),"Wallen"] = 1
DirectionalityHealthAssociation[rownames(df_Compare_Luan[df_Compare_Luan[,1]==-1,]),"Luan"] = 1

df_HealthAssociation <- data.frame("Detection"=rowSums(DirectionalityHealthAssociation),"Influence"=SpeciesScores_NEW[rownames(DirectionalityHealthAssociation),1],"Stability"=SpeciesScores_NEW[rownames(DirectionalityHealthAssociation),2],"Health"=SpeciesScores_NEW[rownames(DirectionalityHealthAssociation),3],"HACKScore"=SpeciesScores_NEW[rownames(DirectionalityHealthAssociation),4])


save.image("241108_AdditionalValidation.RData")