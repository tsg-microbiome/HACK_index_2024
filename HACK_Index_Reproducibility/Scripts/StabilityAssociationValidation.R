library(beanplot)

StudyList_IndustrializedUrban <- c("Baxter_2019", "ClooneyA_2020", "Dahl_2016", "Deehan_2020", "HalfvarsonJ_2017", "IaniroG_2022", "KangJ_2022", "Kovatcheva_2015", "LloydPriceJ_2019", "LouisS_2016", "Mars_2020", "MetwalyA_2020", "NUAGE", "PallejaA_2018", "PoyetM_2019", "RaymondF_2016", "Tap_2015", "VincentC_2016", "WedenojaS_2022")

StudyList_Others <- c("DavidLA_2015", "Morales_2016", "TeeM_2022")

StudyList_WGS <- c("IaniroG_2022", "KangJ_2022", "LloydPriceJ_2019", "LouisS_2016", "Mars_2020", "PallejaA_2018", "RaymondF_2016", "TeeM_2022", "VincentC_2016")

StudyList_16S <- c("Baxter_2019", "ClooneyA_2020", "Dahl_2016", "DavidLA_2015", "Deehan_2020", "HalfvarsonJ_2017", "Healey_2018", "Kovatcheva_2015", "MetwalyA_2020", "Morales_2016", "NUAGE", "PoyetM_2019", "Tap_2015", "WedenojaS_2022")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\StabiltyScoresSubjectGroup.RData")


load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AMIT\\Discovery_and_validation\\Discovery_and_validation_abund_distance_196_species.RData")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\CoreInfluencers_3Dfs.RData")

source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\code_library.R")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\SpeciesScores_NEW.RData")

#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\241025_StabilityAssociation.RData")

hack_top_17 <- rownames(SpeciesScores_NEW[SpeciesScores_NEW[,4]>=0.75,])

Only16S_taxa <- c("Blautia_faecis","Blautia_luti","Eubacterium_desmolans","Oscillibacter_valericigenes","Oscillospira_guilliermondii","Clostridium_lactatifermentans","Clostridium_ruminantium","Acetanaerobacterium_elongatum","Clostridium_thermocellum","Sporobacter_termitidis","Clostridium_glycolicum","Clostridium_sporosphaeroides")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\stability_analysis_spProfile.RData")

combined_df <- as.data.frame(cbind(stability_analysis_spProfile[rownames(combined_rel_abund_followup_dist_all_201sp_df),],combined_rel_abund_followup_dist_all_201sp_df[,colnames(combined_rel_abund_followup_dist_all_201sp_df)[202:213]]))

detected_taxa_16S <- c(detected_taxa,Only16S_taxa)

print("SequencingType")

meta_16S_BrayCurtis <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",StudyList_16S)

meta_16S_BrayCurtis_new <- compute_meta_lm_group(combined_df,detected_taxa_16S,"BrayCurtis_dist","Study_Name",StudyList_16S)

meta_16S_Aitchison <- compute_meta_lm_group(combined_df,detected_taxa,"Aitchison_dist","Study_Name",StudyList_16S)

meta_16S_Aitchison_new <- compute_meta_lm_group(combined_df,detected_taxa_16S,"Aitchison_dist","Study_Name",StudyList_16S)

meta_WGS_BrayCurtis <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",StudyList_WGS)

meta_WGS_Aitchison <- compute_meta_lm_group(combined_df,detected_taxa,"Aitchison_dist","Study_Name",StudyList_WGS)

df_StabilityScore_SequenceType = data.frame("16S"=((rank_scale1((-1) * sign(meta_16S_Aitchison$beta) * (-log(meta_16S_Aitchison$qval,10)) * meta_16S_Aitchison$consistency)) + (rank_scale1((-1) * sign(meta_16S_BrayCurtis$beta) * (-log(meta_16S_BrayCurtis$qval,10)) * meta_16S_BrayCurtis$consistency)))/2,"WGS"=((rank_scale1((-1) * sign(meta_WGS_Aitchison$beta) * (-log(meta_WGS_Aitchison$qval,10)) * meta_WGS_Aitchison$consistency)) + (rank_scale1((-1) * sign(meta_WGS_BrayCurtis$beta) * (-log(meta_WGS_BrayCurtis$qval,10)) * meta_WGS_BrayCurtis$consistency)))/2,row.names=detected_taxa)

StabilityScores_16S_only <- data.frame("16S"=((rank_scale1((-1) * sign(meta_16S_Aitchison_new$beta) * (-log(meta_16S_Aitchison_new$qval,10)) * meta_16S_Aitchison_new$consistency)) + (rank_scale1((-1) * sign(meta_16S_BrayCurtis_new$beta) * (-log(meta_16S_BrayCurtis_new$qval,10)) * meta_16S_BrayCurtis_new$consistency)))/2,row.names=detected_taxa_16S)



df_StabilityScore_SequenceType <- df_StabilityScore_SequenceType[!is.na(df_StabilityScore_SequenceType[,1]),]

df_StabilityScore_SequenceType$CV = apply(df_StabilityScore_SequenceType[,1:2],1,function(x)(sd(x)/mean(x)))

df_StabilityScore_SequenceType$StabilityScore <- SpeciesScores_NEW[rownames(df_StabilityScore_SequenceType),2]
df_StabilityScore_SequenceType$HACKScore <- SpeciesScores_NEW[rownames(df_StabilityScore_SequenceType),4]

ggplot(df_StabilityScore_SequenceType,aes(x=StabilityScore,y=CV))+geom_point(col=ifelse(rownames(df_StabilityScore_SequenceType)%in%hack_top_17,"blue","grey50"),size=5)+theme_bw()+geom_smooth(method="lm",color="firebrick4",linewidth=2)+xlab("")+ylab("")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

boxplot(df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),1],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),1],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

boxplot(df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),2],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),2],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

beanplot(df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),2],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),2],df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),1],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),1],col=list("blue","grey50"),side="both",what=c(1,1,1,0),overallline="median")

print("CohortType")

meta_IndustrializedUrban_BrayCurtis <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",StudyList_IndustrializedUrban)

meta_IndustrializedUrban_Aitchison <- compute_meta_lm_group(combined_df,detected_taxa,"Aitchison_dist","Study_Name",StudyList_IndustrializedUrban)

meta_Others_BrayCurtis <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",StudyList_Others)

meta_Others_Aitchison <- compute_meta_lm_group(combined_df,detected_taxa,"Aitchison_dist","Study_Name",StudyList_Others)

df_StabilityScore_CohortType = data.frame("IndustrializedUrban"=((rank_scale1((-1) * sign(meta_IndustrializedUrban_Aitchison$beta) * (-log(meta_IndustrializedUrban_Aitchison$qval,10)) * meta_IndustrializedUrban_Aitchison$consistency)) + (rank_scale1((-1) * sign(meta_IndustrializedUrban_BrayCurtis$beta) * (-log(meta_IndustrializedUrban_BrayCurtis$qval,10)) * meta_IndustrializedUrban_BrayCurtis$consistency)))/2,"Others"=((rank_scale1((-1) * sign(meta_Others_Aitchison$beta) * (-log(meta_Others_Aitchison$qval,10)) * meta_Others_Aitchison$consistency)) + (rank_scale1((-1) * sign(meta_Others_BrayCurtis$beta) * (-log(meta_Others_BrayCurtis$qval,10)) * meta_Others_BrayCurtis$consistency)))/2,row.names=detected_taxa)


df_StabilityScore_CohortType <- df_StabilityScore_CohortType[!is.na(df_StabilityScore_CohortType[,1]),]

df_StabilityScore_CohortType$CV = apply(df_StabilityScore_CohortType[,1:2],1,function(x)(sd(x)/mean(x)))

df_StabilityScore_CohortType$StabilityScore <- SpeciesScores_NEW[rownames(df_StabilityScore_CohortType),2]
df_StabilityScore_CohortType$HACKScore <- SpeciesScores_NEW[rownames(df_StabilityScore_CohortType),4]

ggplot(df_StabilityScore_CohortType,aes(x=StabilityScore,y=CV))+geom_point(col=ifelse(rownames(df_StabilityScore_CohortType)%in%hack_top_17,"blue","grey50"),size=5)+theme_bw()+geom_smooth(method="lm",color="firebrick4",linewidth=2)+xlab("")+ylab("")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

boxplot(df_StabilityScore_CohortType[intersect(rownames(df_StabilityScore_CohortType),hack_top_17),1],df_StabilityScore_CohortType[setdiff(rownames(df_StabilityScore_CohortType),hack_top_17),1],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

boxplot(df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),2],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),2],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

beanplot(df_StabilityScore_CohortType[intersect(rownames(df_StabilityScore_CohortType),hack_top_17),2],df_StabilityScore_CohortType[setdiff(rownames(df_StabilityScore_CohortType),hack_top_17),2],df_StabilityScore_CohortType[intersect(rownames(df_StabilityScore_CohortType),hack_top_17),1],df_StabilityScore_CohortType[setdiff(rownames(df_StabilityScore_CohortType),hack_top_17),1],col=list("blue","grey50"),side="both",what=c(1,1,1,0),overallline="median")


print("Subject Type")

df_StabilityScore_SubjectType <- as.data.frame(matrix(NA,nrow(df_StabilityScore_CohortType),5))
rownames(df_StabilityScore_SubjectType) <- rownames(df_StabilityScore_CohortType)
colnames(df_StabilityScore_SubjectType) <- c("UC","CD","Control","Others","Fiber")

for(i in 1:nrow(df_StabilityScore_CohortType))
{
	Species <- rownames(df_StabilityScore_CohortType)[i]
	x <- c(MetaLM_BrayCurtis_UC[Species,"Stability_Score"],MetaLM_Aitchison_UC[Species,"Stability_Score"])
	df_StabilityScore_SubjectType[Species,"UC"] <- mean(x[!is.na(x)])
	
	x <- c(MetaLM_BrayCurtis_CD[Species,"Stability_Score"],MetaLM_Aitchison_CD[Species,"Stability_Score"])
	df_StabilityScore_SubjectType[Species,"CD"] <- mean(x[!is.na(x)])
	
	x <- c(MetaLM_BrayCurtis_Others[Species,"Stability_Score"],MetaLM_Aitchison_Others[Species,"Stability_Score"])
	df_StabilityScore_SubjectType[Species,"Others"] <- mean(x[!is.na(x)])
	
	x <- c(MetaLM_BrayCurtis_Control[Species,"Stability_Score"],MetaLM_Aitchison_Control[Species,"Stability_Score"])
	df_StabilityScore_SubjectType[Species,"Control"] <- mean(x[!is.na(x)])
	
	x <- c(MetaLM_BrayCurtis_Fiber[Species,"Stability_Score"],MetaLM_Aitchison_Fiber[Species,"Stability_Score"])
	df_StabilityScore_SubjectType[Species,"Fiber"] <- mean(x[!is.na(x)])
}
	
df_StabilityScore_SubjectType$CV <- apply(df_StabilityScore_SubjectType,1,function(x)(sd(x[!is.na(x)])/mean(x[!is.na(x)])))

df_StabilityScore_SubjectType$Overall <- df_StabilityScore_CohortType$StabilityScore

ggplot(df_StabilityScore_SubjectType,aes(x=Overall,y=CV))+geom_point(col=ifelse(rownames(df_StabilityScore_SubjectType)%in%hack_top_17,"blue","grey50"),size=5)+theme_bw()+geom_smooth(method="lm",color="firebrick4",linewidth=2)+xlab("")+ylab("")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

beanplot(df_StabilityScore_SubjectType[intersect(rownames(df_StabilityScore_SubjectType),hack_top_17),1],df_StabilityScore_SubjectType[setdiff(rownames(df_StabilityScore_SubjectType),hack_top_17),1],df_StabilityScore_SubjectType[intersect(rownames(df_StabilityScore_SubjectType),hack_top_17),2],df_StabilityScore_SubjectType[setdiff(rownames(df_StabilityScore_SubjectType),hack_top_17),2],df_StabilityScore_SubjectType[intersect(rownames(df_StabilityScore_SubjectType),hack_top_17),3],df_StabilityScore_SubjectType[setdiff(rownames(df_StabilityScore_SubjectType),hack_top_17),3],df_StabilityScore_SubjectType[intersect(rownames(df_StabilityScore_SubjectType),hack_top_17),4],df_StabilityScore_SubjectType[setdiff(rownames(df_StabilityScore_SubjectType),hack_top_17),4],df_StabilityScore_SubjectType[intersect(rownames(df_StabilityScore_SubjectType),hack_top_17),5],df_StabilityScore_SubjectType[setdiff(rownames(df_StabilityScore_SubjectType),hack_top_17),5],col=list("blue","grey50"),side="both",what=c(1,1,1,0),overallline="median")

beanplot(df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),1],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),1],df_StabilityScore_SequenceType[intersect(rownames(df_StabilityScore_SequenceType),hack_top_17),2],df_StabilityScore_SequenceType[setdiff(rownames(df_StabilityScore_SequenceType),hack_top_17),2],df_StabilityScore_CohortType[intersect(rownames(df_StabilityScore_CohortType),hack_top_17),1],df_StabilityScore_CohortType[setdiff(rownames(df_StabilityScore_CohortType),hack_top_17),1],df_StabilityScore_CohortType[intersect(rownames(df_StabilityScore_CohortType),hack_top_17),2],df_StabilityScore_CohortType[setdiff(rownames(df_StabilityScore_CohortType),hack_top_17),2],col=list("blue","grey50"),side="both",what=c(1,1,0,0))