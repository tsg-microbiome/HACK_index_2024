load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AMIT\\Discovery_and_validation\\Discovery_and_validation_abund_distance_196_species.RData")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\InfluenceScore.RData")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AMIT\\Discovery_and_validation\\stability_analysis_spProfile.RData")

source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\code_library.R")

#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\CoreKeyStones.RData")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\FunctionalProfiling\\AllSpecies.RData")

combined_rel_abund_followup_dist_full_df <- as.data.frame(rbind(combined_rel_abund_followup_dist_discovery_df,combined_rel_abund_followup_dist_validation_df))

combined_rel_abund_followup_dist_all_201sp_df <- stability_analysis_spProfile[rownames(combined_rel_abund_followup_dist_full_df),names(InfluenceScore)]

MetaDataColumns <- c("Subject_ID","Sample_ID","Time_Point","Follow_up","Time_point_type","Treatment_Condition","Timepoint_difference","Study_Name","Jaccard_dist","Aitchison_dist","Kendall_dist","BrayCurtis_dist")

combined_rel_abund_followup_dist_all_201sp_df[,MetaDataColumns] <- combined_rel_abund_followup_dist_full_df[,MetaDataColumns]

combined_df <- combined_rel_abund_followup_dist_all_201sp_df

detected_taxa <- colnames(combined_df)[1:201]

combined_StudyNames <- unique(combined_df$Study_Name)

combined_BrayCurtisDist <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",combined_StudyNames)
combined_BrayCurtisDist <- combined_BrayCurtisDist[!is.na(combined_BrayCurtisDist$consistency),]
combined_BrayCurtisDist$SpeciesType <- ifelse(rownames(combined_BrayCurtisDist) %in% MajorCoreKeyStone,"CoreConnected","Other")

combined_StudyNames <- unique(combined_df$Study_Name)

combined_BrayCurtisDist <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",combined_StudyNames)
combined_BrayCurtisDist <- combined_BrayCurtisDist[!is.na(combined_BrayCurtisDist$consistency),]
combined_BrayCurtisDist$SpeciesType <- ifelse(rownames(combined_BrayCurtisDist) %in% MajorCoreKeyStone,"CoreKeyStone","Other")

combined_AitchisonDist <- compute_meta_lm_group(combined_df,detected_taxa,"Aitchison_dist","Study_Name",combined_StudyNames)
combined_AitchisonDist <- combined_AitchisonDist[!is.na(combined_AitchisonDist$consistency),]
combined_AitchisonDist$SpeciesType <- ifelse(rownames(combined_AitchisonDist) %in% MajorCoreKeyStone,"CoreKeyStone","Other")

all_study_names <- unique(combined_df$Study_Name)

for(i in 1:10)
{
	iter <- i
	select_studies <- sample(all_study_names,12,replace=FALSE)
	assign(paste0("study_list_iter",i),select_studies)
	meta_BrayCurtis <- compute_meta_lm_group(combined_df,detected_taxa,"BrayCurtis_dist","Study_Name",select_studies)
	assign(paste0("iter",i,"_metaBrayCurtisDist"),meta_BrayCurtis)
	meta_Aitchison <- compute_meta_lm_group(combined_df,detected_taxa,"Aitchison_dist","Study_Name",select_studies)
	assign(paste0("iter",i,"_metaAitchisonDist"),meta_Aitchison)
}

iteration_beta_Aitchison <- data.frame(iter1=iter1_metaAitchisonDist$beta,iter2=iter2_metaAitchisonDist$beta,iter3=iter3_metaAitchisonDist$beta,iter4=iter4_metaAitchisonDist$beta,iter5=iter5_metaAitchisonDist$beta,iter6=iter6_metaAitchisonDist$beta,iter7=iter7_metaAitchisonDist$beta,iter8=iter8_metaAitchisonDist$beta,iter9=iter9_metaAitchisonDist$beta,iter10=iter10_metaAitchisonDist$beta)

iteration_beta_Aitchison <- (-1) * iteration_beta_Aitchison

iteration_qval_Aitchison <- data.frame(iter1=iter1_metaAitchisonDist$qval,iter2=iter2_metaAitchisonDist$qval,iter3=iter3_metaAitchisonDist$qval,iter4=iter4_metaAitchisonDist$qval,iter5=iter5_metaAitchisonDist$qval,iter6=iter6_metaAitchisonDist$qval,iter7=iter7_metaAitchisonDist$qval,iter8=iter8_metaAitchisonDist$qval,iter9=iter9_metaAitchisonDist$qval,iter10=iter10_metaAitchisonDist$qval)

iteration_qval_Aitchison <- (-1) * log(iteration_qval_Aitchison,10)

iteration_stability_Aitchison <- data.frame(iter1_stability = rank_scale1((-1) * iter1_metaAitchisonDist$beta * (-1) * log(iter1_metaAitchisonDist$qval,10) * iter1_metaAitchisonDist$consistency),
iter2_stability = rank_scale1((-1) * iter2_metaAitchisonDist$beta * (-1) * log(iter2_metaAitchisonDist$qval,10) * iter2_metaAitchisonDist$consistency),
iter3_stability = rank_scale1((-1) * iter3_metaAitchisonDist$beta * (-1) * log(iter3_metaAitchisonDist$qval,10) * iter3_metaAitchisonDist$consistency),
iter4_stability = rank_scale1((-1) * iter4_metaAitchisonDist$beta * (-1) * log(iter4_metaAitchisonDist$qval,10) * iter4_metaAitchisonDist$consistency),
iter5_stability = rank_scale1((-1) * iter5_metaAitchisonDist$beta * (-1) * log(iter5_metaAitchisonDist$qval,10) * iter5_metaAitchisonDist$consistency),
iter6_stability = rank_scale1((-1) * iter6_metaAitchisonDist$beta * (-1) * log(iter6_metaAitchisonDist$qval,10) * iter6_metaAitchisonDist$consistency),
iter7_stability = rank_scale1((-1) * iter7_metaAitchisonDist$beta * (-1) * log(iter7_metaAitchisonDist$qval,10) * iter7_metaAitchisonDist$consistency),
iter8_stability = rank_scale1((-1) * iter8_metaAitchisonDist$beta * (-1) * log(iter8_metaAitchisonDist$qval,10) * iter8_metaAitchisonDist$consistency),
iter9_stability = rank_scale1((-1) * iter9_metaAitchisonDist$beta * (-1) * log(iter9_metaAitchisonDist$qval,10) * iter9_metaAitchisonDist$consistency),
iter10_stability = rank_scale1((-1) * iter10_metaAitchisonDist$beta * (-1) * log(iter10_metaAitchisonDist$qval,10) * iter10_metaAitchisonDist$consistency),row.names=detected_taxa)

mean_stability_Aitchison <- apply(iteration_stability_Aitchison,1,function(x)(mean(x[!is.na(x)])))
IQR_stability_Aitchison <- apply(iteration_stability_Aitchison,1,function(x)(IQR(x[!is.na(x)])))
stability_aitchison <- data.frame(mean=mean_stability_Aitchison,IQR=IQR_stability_Aitchison)

iteration_meta_Aitchison <- data.frame(median_beta = apply(iteration_beta_Aitchison,1,function(x)(median(x[!is.na(x)]))), iqr_beta = apply(iteration_beta_Aitchison,1,function(x)(IQR(x[!is.na(x)]))), median_qval = apply(iteration_qval_Aitchison,1,function(x)(median(x[!is.na(x)]))), iqr_qval = apply(iteration_qval_Aitchison,1,function(x)(IQR(x[!is.na(x)]))),row.names=detected_taxa)

ggplot(iteration_meta_Aitchison,aes(x=median_beta,y=median_qval))+geom_point(color=ifelse(rownames(iteration_meta_Aitchison) %in% MajorCoreKeyStone, "cornflowerblue", "red"))+geom_errorbar(aes(ymin=median_qval-iqr_qval,ymax=median_qval+iqr_qval),alpha=0.2)+geom_errorbar(aes(xmin=median_beta-iqr_beta,xmax=median_beta+iqr_beta),alpha=0.2)+theme_bw()+geom_text_repel(label=ifelse(rownames(iteration_meta_Aitchison) %in% MajorCoreKeyStone, rownames(iteration_meta_Aitchison), ""),box.padding=0.0001,max.overlaps=50,point.padding=0.001,size=3) + ylim(c(0,9.5)) + geom_hline(yintercept=0)+geom_vline(xintercept=0)

iteration_beta_BrayCurtis <- data.frame(iter1=iter1_metaBrayCurtisDist$beta,iter2=iter2_metaBrayCurtisDist$beta,iter3=iter3_metaBrayCurtisDist$beta,iter4=iter4_metaBrayCurtisDist$beta,iter5=iter5_metaBrayCurtisDist$beta,iter6=iter6_metaBrayCurtisDist$beta,iter7=iter7_metaBrayCurtisDist$beta,iter8=iter8_metaBrayCurtisDist$beta,iter9=iter9_metaBrayCurtisDist$beta,iter10=iter10_metaBrayCurtisDist$beta)

iteration_beta_BrayCurtis <- (-1) * iteration_beta_BrayCurtis

iteration_qval_BrayCurtis <- data.frame(iter1=iter1_metaBrayCurtisDist$qval,iter2=iter2_metaBrayCurtisDist$qval,iter3=iter3_metaBrayCurtisDist$qval,iter4=iter4_metaBrayCurtisDist$qval,iter5=iter5_metaBrayCurtisDist$qval,iter6=iter6_metaBrayCurtisDist$qval,iter7=iter7_metaBrayCurtisDist$qval,iter8=iter8_metaBrayCurtisDist$qval,iter9=iter9_metaBrayCurtisDist$qval,iter10=iter10_metaBrayCurtisDist$qval)

iteration_qval_BrayCurtis <- (-1) * log(iteration_qval_BrayCurtis,10)

iteration_meta_BrayCurtis <- data.frame(median_beta = apply(iteration_beta_BrayCurtis,1,function(x)(median(x[!is.na(x)]))), iqr_beta = apply(iteration_beta_BrayCurtis,1,function(x)(IQR(x[!is.na(x)]))), median_qval = apply(iteration_qval_BrayCurtis,1,function(x)(median(x[!is.na(x)]))), iqr_qval = apply(iteration_qval_BrayCurtis,1,function(x)(IQR(x[!is.na(x)]))),row.names=detected_taxa)

ggplot(iteration_meta_BrayCurtis,aes(x=median_beta,y=median_qval))+geom_point(color=ifelse(rownames(iteration_meta_BrayCurtis) %in% MajorCoreKeyStone, "cornflowerblue", "red"))+geom_errorbar(aes(ymin=median_qval-iqr_qval,ymax=median_qval+iqr_qval),alpha=0.2)+geom_errorbar(aes(xmin=median_beta-iqr_beta,xmax=median_beta+iqr_beta),alpha=0.2)+theme_bw()+geom_text_repel(label=ifelse(rownames(iteration_meta_BrayCurtis) %in% MajorCoreKeyStone, rownames(iteration_meta_BrayCurtis), ""),box.padding=0.0001,max.overlaps=50,point.padding=0.001,size=3) + ylim(c(0,9.5)) + geom_hline(yintercept=0)+geom_vline(xintercept=0)

stability_AitchisonDist <- data.frame(mean=apply(iteration_stability_Aitchison,1,function(x)(mean(x[!is.na(x)]))),IQR=apply(iteration_stability_Aitchison,1,function(x)(IQR(x[!is.na(x)]))))

iteration_stability_BrayCurtis <- data.frame(iter1_stability = rank_scale1((-1) * iter1_metaBrayCurtisDist$beta * (-1) * log(iter1_metaBrayCurtisDist$qval,10) * iter1_metaBrayCurtisDist$consistency),
iter2_stability = rank_scale1((-1) * iter2_metaBrayCurtisDist$beta * (-1) * log(iter2_metaBrayCurtisDist$qval,10) * iter2_metaBrayCurtisDist$consistency),
iter3_stability = rank_scale1((-1) * iter3_metaBrayCurtisDist$beta * (-1) * log(iter3_metaBrayCurtisDist$qval,10) * iter3_metaBrayCurtisDist$consistency),
iter4_stability = rank_scale1((-1) * iter4_metaBrayCurtisDist$beta * (-1) * log(iter4_metaBrayCurtisDist$qval,10) * iter4_metaBrayCurtisDist$consistency),
iter5_stability = rank_scale1((-1) * iter5_metaBrayCurtisDist$beta * (-1) * log(iter5_metaBrayCurtisDist$qval,10) * iter5_metaBrayCurtisDist$consistency),
iter6_stability = rank_scale1((-1) * iter6_metaBrayCurtisDist$beta * (-1) * log(iter6_metaBrayCurtisDist$qval,10) * iter6_metaBrayCurtisDist$consistency),
iter7_stability = rank_scale1((-1) * iter7_metaBrayCurtisDist$beta * (-1) * log(iter7_metaBrayCurtisDist$qval,10) * iter7_metaBrayCurtisDist$consistency),
iter8_stability = rank_scale1((-1) * iter8_metaBrayCurtisDist$beta * (-1) * log(iter8_metaBrayCurtisDist$qval,10) * iter8_metaBrayCurtisDist$consistency),
iter9_stability = rank_scale1((-1) * iter9_metaBrayCurtisDist$beta * (-1) * log(iter9_metaBrayCurtisDist$qval,10) * iter9_metaBrayCurtisDist$consistency),
iter10_stability = rank_scale1((-1) * iter10_metaBrayCurtisDist$beta * (-1) * log(iter10_metaBrayCurtisDist$qval,10) * iter10_metaBrayCurtisDist$consistency),row.names=detected_taxa)

stability_BrayCurtisDist <- data.frame(mean=apply(iteration_stability_BrayCurtis,1,function(x)(mean(x[!is.na(x)]))),IQR=apply(iteration_stability_BrayCurtis,1,function(x)(IQR(x[!is.na(x)]))))

stabilityRank <- as.data.frame(cbind(stability_AitchisonDist,stability_BrayCurtisDist))

stabilityRank <- stabilityRank[!is.na(stabilityRank[,1]),]

colnames(stabilityRank) <- c("mean_Aitchison","IQR_Aitchison","mean_BrayCurtis","IQR_BrayCurtis")

stabilityRank$MeanStabilityScore <- apply(stabilityRank[,c(1,3)],1,mean)

stabilityRank <- stabilityRank[order(stabilityRank$MeanStabilityScore),]
stabilityRank$Rank <- c(1:nrow(stabilityRank))

stabilityRank <- stabilityRank[rev(rownames(stabilityRank)),]

stabilityRank$Rank <- rev(stabilityRank$Rank)

StabilityScore <- stabilityRank$MeanStabilityScore
names(StabilityScore) <- rownames(stabilityRank)

Treatment_Conditions <- rownames(table(combined_df[,c("Treatment_Condition","Study_Name")]))

## Condition: CD
data_CD <- combined_df[combined_df$Treatment_Condition == "CD",]

MetaLM_BrayCurtis_CD <- compute_meta_lm_group(data_CD,detected_taxa,"BrayCurtis_dist","Study_Name",unique(data_CD$Study_Name))

MetaLM_BrayCurtis_CD <- MetaLM_BrayCurtis_CD[!is.na(MetaLM_BrayCurtis_CD[,1]),]

MetaLM_BrayCurtis_CD$SpeciesType <- ifelse(rownames(MetaLM_BrayCurtis_CD) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_BrayCurtis_CD$Stability_Score <- rank_scale1((-1) * MetaLM_BrayCurtis_CD$beta * (-1) * log(MetaLM_BrayCurtis_CD$qval,10) * MetaLM_BrayCurtis_CD$consistency)

MetaLM_Aitchison_CD <- compute_meta_lm_group(data_CD,detected_taxa,"Aitchison_dist","Study_Name",unique(data_CD$Study_Name))

MetaLM_Aitchison_CD <- MetaLM_Aitchison_CD[!is.na(MetaLM_Aitchison_CD[,1]),]

MetaLM_Aitchison_CD$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_CD) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_CD$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_CD$beta * (-1) * log(MetaLM_Aitchison_CD$qval,10) * MetaLM_Aitchison_CD$consistency)


## Condition: Control
data_Control <- combined_df[combined_df$Treatment_Condition == "Control",]

MetaLM_BrayCurtis_Control <- compute_meta_lm_group(data_Control,detected_taxa,"BrayCurtis_dist","Study_Name",unique(data_Control$Study_Name))

MetaLM_BrayCurtis_Control <- MetaLM_BrayCurtis_Control[!is.na(MetaLM_BrayCurtis_Control[,1]),]

MetaLM_BrayCurtis_Control$SpeciesType <- ifelse(rownames(MetaLM_BrayCurtis_Control) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_BrayCurtis_Control$Stability_Score <- rank_scale1((-1) * MetaLM_BrayCurtis_Control$beta * (-1) * log(MetaLM_BrayCurtis_Control$qval,10) * MetaLM_BrayCurtis_Control$consistency)

MetaLM_Aitchison_Control <- compute_meta_lm_group(data_Control,detected_taxa,"Aitchison_dist","Study_Name",unique(data_Control$Study_Name))

MetaLM_Aitchison_Control <- MetaLM_Aitchison_Control[!is.na(MetaLM_Aitchison_Control[,1]),]

MetaLM_Aitchison_Control$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_Control) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_Control$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_Control$beta * (-1) * log(MetaLM_Aitchison_Control$qval,10) * MetaLM_Aitchison_Control$consistency)

## Condition: Fiber
data_Fiber <- combined_df[combined_df$Treatment_Condition == "Fiber",]

MetaLM_Aitchison_Fiber <- compute_meta_lm_group(data_Fiber,detected_taxa,"Aitchison_dist","Study_Name",unique(data_Fiber$Study_Name))

MetaLM_Aitchison_Fiber <- MetaLM_Aitchison_Fiber[!is.na(MetaLM_Aitchison_Fiber[,1]),]

MetaLM_Aitchison_Fiber$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_Fiber) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_Fiber$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_Fiber$beta * (-1) * log(MetaLM_Aitchison_Fiber$qval,10) * MetaLM_Aitchison_Fiber$consistency)

#MetaLM_BrayCurtis_Fiber <- compute_meta_lm_group(data_Fiber,detected_taxa,"BrayCurtis_dist","Study_Name",unique(data_Fiber$Study_Name))

MetaLM_BrayCurtis_Fiber <- MetaLM_Aitchison_Fiber

MetaLM_BrayCurtis_Fiber <- MetaLM_BrayCurtis_Fiber[!is.na(MetaLM_BrayCurtis_Fiber[,1]),]

MetaLM_BrayCurtis_Fiber$SpeciesType <- ifelse(rownames(MetaLM_BrayCurtis_Fiber) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_BrayCurtis_Fiber$Stability_Score <- rank_scale1((-1) * MetaLM_BrayCurtis_Fiber$beta * (-1) * log(MetaLM_BrayCurtis_Fiber$qval,10) * MetaLM_BrayCurtis_Fiber$consistency)

MetaLM_Aitchison_Fiber <- compute_meta_lm_group(data_Fiber,detected_taxa,"Aitchison_dist","Study_Name",unique(data_Fiber$Study_Name))

MetaLM_Aitchison_Fiber <- MetaLM_Aitchison_Fiber[!is.na(MetaLM_Aitchison_Fiber[,1]),]

MetaLM_Aitchison_Fiber$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_Fiber) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_Fiber$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_Fiber$beta * (-1) * log(MetaLM_Aitchison_Fiber$qval,10) * MetaLM_Aitchison_Fiber$consistency)

## Condition: UC
data_UC <- combined_df[combined_df$Treatment_Condition == "UC",]

MetaLM_BrayCurtis_UC <- compute_meta_lm_group(data_UC,detected_taxa,"BrayCurtis_dist","Study_Name",unique(data_UC$Study_Name))

MetaLM_BrayCurtis_UC <- MetaLM_BrayCurtis_Fiber[!is.na(MetaLM_BrayCurtis_UC[,1]),]

MetaLM_BrayCurtis_UC$SpeciesType <- ifelse(rownames(MetaLM_BrayCurtis_UC) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_BrayCurtis_UC$Stability_Score <- rank_scale1((-1) * MetaLM_BrayCurtis_UC$beta * (-1) * log(MetaLM_BrayCurtis_UC$qval,10) * MetaLM_BrayCurtis_UC$consistency)

MetaLM_Aitchison_UC <- compute_meta_lm_group(data_UC,detected_taxa,"Aitchison_dist","Study_Name",unique(data_UC$Study_Name))

MetaLM_Aitchison_UC <- MetaLM_Aitchison_UC[!is.na(MetaLM_Aitchison_UC[,1]),]

MetaLM_Aitchison_UC$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_UC) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_UC$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_UC$beta * (-1) * log(MetaLM_Aitchison_UC$qval,10) * MetaLM_Aitchison_UC$consistency)

## Condition: Other
data_Others <- combined_df[combined_df$Treatment_Condition == "Others",]

MetaLM_BrayCurtis_Others <- compute_meta_lm_group(data_Others,detected_taxa,"BrayCurtis_dist","Study_Name",unique(data_Others$Study_Name))

MetaLM_BrayCurtis_Others <- MetaLM_BrayCurtis_Others[!is.na(MetaLM_BrayCurtis_Others[,1]),]

MetaLM_BrayCurtis_Others$SpeciesType <- ifelse(rownames(MetaLM_BrayCurtis_Others) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_BrayCurtis_Others$Stability_Score <- rank_scale1((-1) * MetaLM_BrayCurtis_Others$beta * (-1) * log(MetaLM_BrayCurtis_Others$qval,10) * MetaLM_BrayCurtis_Others$consistency)

MetaLM_Aitchison_Others <- compute_meta_lm_group(data_Others,detected_taxa,"Aitchison_dist","Study_Name",unique(data_Others$Study_Name))

MetaLM_Aitchison_Others <- MetaLM_Aitchison_Others[!is.na(MetaLM_Aitchison_Others[,1]),]

MetaLM_Aitchison_Others$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_Others) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_Others$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_Others$beta * (-1) * log(MetaLM_Aitchison_Others$qval,10) * MetaLM_Aitchison_Others$consistency)

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AMIT\\Discovery_and_validation\\Discovery_validation_196sp_with_abund_difference_196_species.RData")

combined_diff_df <- as.data.frame(rbind(combined_rel_abund_followup_dist_discovery_196sp_df_diff[,c(1:196)],combined_rel_abund_followup_dist_validation_196sp_df_diff[,c(1:196)]))

combined_full_df <- as.data.frame(cbind(combined_df,combined_diff_df))

compute_meta_lm_single_adjust <- function(data,var1,var2,var3,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	grouping_list_new <- grouping_list
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		tryCatch(               
							expr = {                     
										f <- as.formula(paste0(var1,"~",var3,"+",var2))
										temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
										summary_temp_rlm <- summary(temp_rlm)
										temp_meta[i,2] <- summary_temp_rlm$coefficients[3,3]
										temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
										temp_meta[i,4] <- 1
										temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
										temp_meta[i,6] <- sign(temp_meta[i,2])
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
		
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		#print(temp_meta)
	}
	temp_meta <- temp_meta[temp_meta[,2]!=0,]
	#temp_meta <- mutate(temp_meta,study_id=grouping_list_new)
	grouping_list_new <- unique(temp_meta$dataset)
	print(dim(temp_meta))
	print(grouping_list_new)
	print(temp_meta)
	rownames(temp_meta) <- grouping_list_new
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

compute_meta_lm_group_diff_adjust <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
	return_out <- as.data.frame(matrix(NA,length(feature_list),9))
	rownames(return_out) <- feature_list
	colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir")
	for(i in 1:length(feature_list))
	{
		species_name <- feature_list[i]
		diff_species_name <- paste0("diff_",species_name)
		tryCatch(               
							expr = {
									temp_data <- data[,c(species_name,diff_species_name,metadata_var,grouping_var)]
									grouping_list <- unique(temp_data$Study_Name)
									temp_res <- compute_meta_lm_single_adjust(combined_full_df,metadata_var,species_name,diff_species_name,grouping_var,grouping_list)
									return_out[i,"beta"] <- temp_res$model$beta
									return_out[i,"pval"] <- temp_res$model$pval
									return_out[i,"ci.ub"] <- temp_res$model$ci.ub
									return_out[i,"ci.lb"] <- temp_res$model$ci.lb
									return_out[i,"tau2"] <- temp_res$model$tau2
									return_out[i,"QE"] <- temp_res$model$QE
									return_out[i,"QEp"] <- temp_res$model$QEp	
									return_out[i,"consistency"] <- length(which(sign(temp_res$df_studies$di) == sign(as.numeric(temp_res$model$beta))))/nrow(temp_res$df_studies)	
		
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
		
		
		
		
	}
	return_out$qval <- p.adjust(return_out$pval,method="fdr")
	return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta)))
	return(return_out)
}

MetaLM_BrayCurtis_AbundanceDiff_Adjusted <- compute_meta_lm_group_diff_adjust(combined_full_df,detected_taxa,"BrayCurtis_dist","Study_Name",unique(combined_full_df$Study_Name))

MetaLM_BrayCurtis_AbundanceDiff_Adjusted <- MetaLM_BrayCurtis_AbundanceDiff_Adjusted[!is.na(MetaLM_Aitchison_UC[,1]),]

MetaLM_BrayCurtis_AbundanceDiff_Adjusted$SpeciesType <- ifelse(rownames(MetaLM_BrayCurtis_AbundanceDiff_Adjusted) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_BrayCurtis_AbundanceDiff_Adjusted$Stability_Score <- rank_scale1((-1) * MetaLM_BrayCurtis_AbundanceDiff_Adjusted$beta * (-1) * log(MetaLM_BrayCurtis_AbundanceDiff_Adjusted$qval,10) * MetaLM_BrayCurtis_AbundanceDiff_Adjusted$consistency)

MetaLM_Aitchison_AbundanceDiff_Adjusted <- compute_meta_lm_group_diff_adjust(combined_full_df,detected_taxa,"Aitchison_dist","Study_Name",unique(combined_full_df$Study_Name))

MetaLM_Aitchison_AbundanceDiff_Adjusted <- MetaLM_Aitchison_AbundanceDiff_Adjusted[!is.na(MetaLM_Aitchison_UC[,1]),]

MetaLM_Aitchison_AbundanceDiff_Adjusted$SpeciesType <- ifelse(rownames(MetaLM_Aitchison_AbundanceDiff_Adjusted) %in% MajorCoreKeyStone,"CoreKeyStone","Others")

MetaLM_Aitchison_AbundanceDiff_Adjusted$Stability_Score <- rank_scale1((-1) * MetaLM_Aitchison_AbundanceDiff_Adjusted$beta * (-1) * log(MetaLM_Aitchison_AbundanceDiff_Adjusted$qval,10) * MetaLM_Aitchison_AbundanceDiff_Adjusted$consistency)

StabilityScoreVariants <- as.data.frame(matrix(0,length(detected_taxa),7))
rownames(StabilityScoreVariants) <- detected_taxa
colnames(StabilityScoreVariants) <- c("OverallStabilityScore","AbundanceDiffAdjustedStability","StabilityUC","StabilityCD","StabilityControl","StabilityFiber","StabilityOthers")

for(i in 1:length(detected_taxa))
{
	species <- detected_taxa[i]
	StabilityScoreVariants[species,1] <- stabilityRank[species,"MeanStabilityScore"]
	StabilityScoreVariants[species,2] <- (MetaLM_Aitchison_AbundanceDiff_Adjusted[species,"Stability_Score"]+MetaLM_BrayCurtis_AbundanceDiff_Adjusted[species,"Stability_Score"])/2
	StabilityScoreVariants[species,3] <- (MetaLM_Aitchison_UC[species,"Stability_Score"]+MetaLM_BrayCurtis_UC[species,"Stability_Score"])/2
	StabilityScoreVariants[species,4] <- (MetaLM_Aitchison_CD[species,"Stability_Score"]+MetaLM_BrayCurtis_CD[species,"Stability_Score"])/2
	StabilityScoreVariants[species,5] <- (MetaLM_Aitchison_Control[species,"Stability_Score"]+MetaLM_BrayCurtis_Control[species,"Stability_Score"])/2
	StabilityScoreVariants[species,6] <- (MetaLM_Aitchison_Fiber[species,"Stability_Score"]+MetaLM_BrayCurtis_Fiber[species,"Stability_Score"])/2
	StabilityScoreVariants[species,7] <- (MetaLM_Aitchison_Others[species,"Stability_Score"]+MetaLM_BrayCurtis_Others[species,"Stability_Score"])/2
}

InfluenceStabilityScore <- data.frame(StabilityScore=StabilityScore,InfluenceScore=InfluenceScore[names(StabilityScore)],row.names=names(StabilityScore))

ggplot(InfluenceStabilityScore,aes(x=StabilityScore,y=InfluenceScore))+geom_point()+geom_text_repel(label=ifelse(InfluenceStabilityScore$StabilityScore>=0.70&InfluenceStabilityScore$InfluenceScore>=0.70,sub("_"," ",rownames(InfluenceStabilityScore)),""),box.padding=0.1,point.padding=0.0001,color="blue2")+ylim(c(0.25,1.0))+xlim(c(0.25,1.05))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),text=element_text(family="TT Times New Roman"))+geom_vline(xintercept=0.70)+geom_hline(yintercept=0.7)

ggplot(InfluenceStabilityScore,aes(x=StabilityScore,y=InfluenceScore))+geom_point()+geom_text_repel(label=sub("_"," ",rownames(InfluenceStabilityScore)),box.padding=0.1,point.padding=0.0001,color=ifelse(InfluenceStabilityScore$StabilityScore>=0.70&InfluenceStabilityScore$InfluenceScore>=0.70,"darkblue",ifelse(InfluenceStabilityScore$StabilityScore>=0.70&InfluenceStabilityScore$InfluenceScore<0.70,"mediumorchid4",ifelse(InfluenceStabilityScore$StabilityScore<0.70&InfluenceStabilityScore$InfluenceScore>=0.70,"mediumorchid4",ifelse(InfluenceStabilityScore$StabilityScore<0.70&InfluenceStabilityScore$InfluenceScore<0.70,"red2","white")))))+ylim(c(0.25,1.0))+xlim(c(0.25,1.05))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),text=element_text(family="TT Times New Roman"))+geom_vline(xintercept=0.70)+geom_hline(yintercept=0.7)

ggplot(iteration_meta_Aitchison,aes(x=(-1)*median_beta,y=median_qval))+geom_point()+theme_bw()+geom_text_repel(label=ifelse(iteration_meta_Aitchison$median_qval>=1,sub("_"," ",rownames(iteration_meta_Aitchison)),""),max.overlaps=30,box.padding=0.01,point.padding=0.00001,color="blue2",size=3) + geom_hline(yintercept=0) + geom_vline(xintercept=0)

ggplot(iteration_meta_BrayCurtis,aes(x=(-1)*median_beta,y=median_qval))+geom_point()+theme_bw()+geom_text_repel(label=ifelse(iteration_meta_BrayCurtis$median_qval>=1,sub("_"," ",rownames(iteration_meta_BrayCurtis)),""),max.overlaps=30,box.padding=0.01,point.padding=0.00001,color="blue2",size=3) + geom_hline(yintercept=0) + geom_vline(xintercept=0)

ggplot(MetaLM_Aitchison_AbundanceDiff_Adjusted,aes(x=beta,y=ifelse(-log(qval,10)>=6,6,-log(qval,10))))+geom_point()+theme_bw()+geom_text_repel(label=ifelse(MetaLM_Aitchison_AbundanceDiff_Adjusted$qval<=0.10,sub("_"," ",rownames(MetaLM_Aitchison_AbundanceDiff_Adjusted)),""),max.overlaps=30,box.padding=0.01,point.padding=0.00001,color="blue2",size=3) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + xlab("Estimate") + ylab("-log(Q,10)")

ggplot(MetaLM_BrayCurtis_AbundanceDiff_Adjusted,aes(x=beta,y=ifelse(-log(qval,10)>=6,6,-log(qval,10))))+geom_point()+theme_bw()+geom_text_repel(label=ifelse(MetaLM_BrayCurtis_AbundanceDiff_Adjusted$qval<=0.10,sub("_"," ",rownames(MetaLM_BrayCurtis_AbundanceDiff_Adjusted)),""),max.overlaps=30,box.padding=0.01,point.padding=0.00001,color="blue2",size=3) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + xlab("Estimate") + ylab("-log(Q,10)")

iteration_meta_BrayCurtis$dir <- ifelse(iteration_meta_BrayCurtis$median_qval>=1,2*sign(iteration_meta_BrayCurtis$median_beta),sign(iteration_meta_BrayCurtis$median_beta))

iteration_meta_Aitchison$dir <- ifelse(iteration_meta_Aitchison$median_qval>=1,2*sign(iteration_meta_Aitchison$median_beta),sign(iteration_meta_Aitchison$median_beta))

Directionality_MajorCoreKeyStone <- data.frame(Overall_Aitchison=combined_AitchisonDist[MajorCoreKeyStone,"dir"],Aitchison_Adjusted=MetaLM_Aitchison_AbundanceDiff_Adjusted[MajorCoreKeyStone,"dir"],Aitchison_UC=MetaLM_Aitchison_UC[MajorCoreKeyStone,"dir"],Aitchison_CD=MetaLM_Aitchison_CD[MajorCoreKeyStone,"dir"],Aitchison_Controls=MetaLM_Aitchison_Control[MajorCoreKeyStone,"dir"],Aitchison_Fiber=MetaLM_Aitchison_Fiber[MajorCoreKeyStone,"dir"],Aitchison_Others=MetaLM_Aitchison_Others[MajorCoreKeyStone,"dir"],Overall_BrayCurtis=combined_BrayCurtisDist[MajorCoreKeyStone,"dir"],BrayCurtis_Adjusted=MetaLM_BrayCurtis_AbundanceDiff_Adjusted[MajorCoreKeyStone,"dir"],BrayCurtis_UC=MetaLM_BrayCurtis_UC[MajorCoreKeyStone,"dir"],BrayCurtis_CD=MetaLM_BrayCurtis_CD[MajorCoreKeyStone,"dir"],BrayCurtis_Controls=MetaLM_BrayCurtis_Control[MajorCoreKeyStone,"dir"],BrayCurtis_Fiber=MetaLM_BrayCurtis_Fiber[MajorCoreKeyStone,"dir"],BrayCurtis_Others=MetaLM_BrayCurtis_Others[MajorCoreKeyStone,"dir"],row.names=MajorCoreKeyStone)

Directionality_MajorCoreKeyStone <- apply(Directionality_MajorCoreKeyStone,2,function(x)(ifelse(is.na(x),0,x)))

StabilityScoreVariants <- data.frame(OverallStabilityScore=stabilityRank$MeanStabilityScore,AbundanceDiffAdjustedStability=(MetaLM_Aitchison_AbundanceDiff_Adjusted[rownames(stabilityRank),"Stability_Score"]+MetaLM_BrayCurtis_AbundanceDiff_Adjusted[rownames(stabilityRank),"Stability_Score"])/2,StabilityScoreUC = (MetaLM_Aitchison_UC[rownames(stabilityRank),"Stability_Score"]+MetaLM_BrayCurtis_UC[rownames(stabilityRank),"Stability_Score"])/2,StabilityScoreCD = (MetaLM_Aitchison_CD[rownames(stabilityRank),"Stability_Score"]+MetaLM_BrayCurtis_CD[rownames(stabilityRank),"Stability_Score"])/2,StabilityScoreControl = (MetaLM_Aitchison_Control[rownames(stabilityRank),"Stability_Score"]+MetaLM_BrayCurtis_Control[rownames(stabilityRank),"Stability_Score"])/2,StabilityScoreFiber = (MetaLM_Aitchison_Fiber[rownames(stabilityRank),"Stability_Score"]+MetaLM_BrayCurtis_Fiber[rownames(stabilityRank),"Stability_Score"])/2,StabilityScoreOther = (MetaLM_Aitchison_Others[rownames(stabilityRank),"Stability_Score"]+MetaLM_BrayCurtis_Others[rownames(stabilityRank),"Stability_Score"])/2,row.names=rownames(stabilityRank))

batch_rlm_grouped_ab_adjust <- function(data,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	df_est <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_est) <- grouping_list
	colnames(df_est) <- variable_group
	df_p_val <- as.data.frame(matrix(1,length(grouping_list),length(variable_group)))
	rownames(df_p_val) <- grouping_list
	colnames(df_p_val) <- variable_group
	for(i in 1:length(grouping_list))
	{
		study_name <- grouping_list[i]
		print(study_name)
		study_samples <- rownames(metadata[metadata[,grouping_feature] == study_name,])
		for(j in 1:length(variable_group))
		{
			tryCatch(               
						expr = {                     
									species_name <- variable_group[j]
									diff_species_name <- paste0("diff_",species_name)
									print(species_name)
									print(length(species_val[species_val>0]))
									if(length(species_val[species_val>0])>0)
									{
										temp_rlm <- rlm(as.formula(paste0(metadata_feature,"~",diff_species_name,"+",species_name)),data=data)
										summary_temp_rlm <- summary(temp_rlm)
										print(summary_temp_rlm)
										df_est[i,j] <- summary_temp_rlm$coefficients[3,3]
										df_p_val[i,j] <- f.robftest(temp_rlm,var=species_name)$p.value
									}
								},
								error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
								},
								finally = {            
										print("finally Executed")
								}
						)
				

			
		}
		
	}
	df_q_val <- apply(df_p_val,2,p.adjust)
	l_fisher <- p.adjust(apply(df_q_val,2,function(x)(sumlog(x)$p)),method="fdr")
	print("l_fisher generated")
	df_dir <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_dir) <- grouping_list
	colnames(df_dir) <- variable_group
	for(i in 1:length(grouping_list))
	{
		for(j in 1:length(variable_group))
		{
			df_dir[i,j] <- ifelse(df_q_val[i,j]<=0.10,3*sign(df_est[i,j]),ifelse(df_p_val[i,j]<=0.05,2*sign(df_est[i,j]),1*sign(df_est[i,j])))
		}
	}
	return_list <- list("est"=df_est,"p.value"=df_p_val,"q.value"=df_q_val,"fisher"=l_fisher,"dir"=df_dir)
	return(return_list)
}

StabilityScore <- round(StabilityScore,2)

save(StabilityScore,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\StabilityScore.RData")