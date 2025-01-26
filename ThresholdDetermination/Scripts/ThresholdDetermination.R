library("reshape2")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\For_Upload\\ThresholdDetermination\\controlSpProfileAll.RData")
load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\For_Upload\\ThresholdDetermination\\controlSpProfile.RData")
load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\For_Upload\\ThresholdDetermination\\3R.RData")

source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\code_library.r")

print("Determination of Thresholds")
controlSpProfileAll <- controlSpProfileAll[rownames(controlSpProfile),]
controlSpProfileAll$study_name <- controlSpProfile[,"study_name"]

study_list <- colnames(prevalentDf)

AllSpeciesDetectionPattern <- compute_detection(controlSpProfileAll,setdiff(colnames(controlSpProfileAll),"study_name"),"study_name",study_list)

AllSpeciesMeanAbundance <- aggregate(controlSpProfileAll[,setdiff(colnames(controlSpProfileAll),"study_name")],by=list(controlSpProfileAll$study_name),FUN=mean)[,-1]
rownames(AllSpeciesMeanAbundance) <- aggregate(controlSpProfileAll[,setdiff(colnames(controlSpProfileAll),"study_name")],by=list(controlSpProfileAll$study_name),FUN=mean)[,1]
AllSpeciesMeanAbundance <- t(AllSpeciesMeanAbundance)

study_thresholds <- seq(0.05,0.95,by=0.05)
detection_thresholds <- seq(0.05,0.95,by=0.05)

df_numb_species <- matrix(NA,length(study_thresholds),length(detection_thresholds))
rownames(df_numb_species) <- study_thresholds
colnames(df_numb_species) <- detection_thresholds

df_representation_90_plus <- matrix(NA,length(study_thresholds),length(detection_thresholds))
rownames(df_representation_90_plus) <- study_thresholds
colnames(df_representation_90_plus) <- detection_thresholds

df_representation_70_minus <- matrix(NA,length(study_thresholds),length(detection_thresholds))
rownames(df_representation_70_minus) <- study_thresholds
colnames(df_representation_70_minus) <- detection_thresholds

for(i in 1:length(detection_thresholds))
{
	detect <- detection_thresholds[i]
	for(j in 1:length(study_thresholds))
	{
		print(j)
		study_perc <- study_thresholds[j]
		df_numb_species[i,j] <- length(which(apply(AllSpeciesDetectionPattern,1,function(x)(length(x[x>=detect])))/ncol(AllSpeciesDetectionPattern)>=study_perc))
		temp_species <- names(which(apply(AllSpeciesDetectionPattern,1,function(x)(length(x[x>=detect])))/ncol(AllSpeciesDetectionPattern)>=study_perc))
		if(length(temp_species) > 1)
		{
			df_representation_90_plus[i,j] <- length(which(colSums(AllSpeciesMeanAbundance[temp_species,])>=0.90))/72
			
			df_representation_70_minus[i,j] <- length(which(colSums(AllSpeciesMeanAbundance[temp_species,])<0.70))/72
		}
		else
		{
			df_representation_90_plus[i,j] <- length(which(sum(AllSpeciesMeanAbundance[temp_species,])>=0.90))/72
			
			df_representation_70_minus[i,j] <- length(which(sum(AllSpeciesMeanAbundance[temp_species,])<0.70))/72
		}
	}
}

df_gut_associated_identification <- data.frame("number_of_species"=as.numeric(df_numb_species),"representation_90_plus"=as.numeric(df_representation_90_plus),"representation_70_minus"=as.numeric(df_representation_70_minus))

ggplot(df_gut_associated_identification)+geom_point(aes(x=number_of_species,y=representation_90_plus),col="chartreuse2",size=4)+geom_point(aes(x=number_of_species,y=representation_70_minus),col="pink2",size=4)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(size=1,color="grey30"))

AllGutAssociatedSpecies <- names(which(apply(AllSpeciesDetectionPattern,1,function(x)(length(x[x>=0.05])))/ncol(AllSpeciesDetectionPattern)>=0.50))[1:201]

selectControlSpProfile <- controlSpProfileAll[,AllGutAssociatedSpecies]

save(selectControlSpProfile,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\selectControlSpProfile.RData")

study_list <- colnames(prevalentDf)

detection_with_diff_thresholds <- matrix(NA,length(detection_thresholds),ncol(prevalentDf))
rownames(detection_with_diff_thresholds) <- detection_thresholds
colnames(detection_with_diff_thresholds) <- colnames(prevalentDf)

for(i in 1:length(detection_thresholds))
{
	thres <- detection_thresholds[i]
	detection_with_diff_thresholds[i,] <- apply(AllSpeciesDetectionPattern[AllGutAssociatedSpecies,],2,function(x)(length(x[x>=thres])))
}



df_jaccard <- matrix(NA,length(detection_thresholds),length(study_list))
rownames(df_jaccard) <- detection_thresholds
colnames(df_jaccard) <- study_list

df_representation <- matrix(NA,length(detection_thresholds),length(study_list))
rownames(df_representation) <- detection_thresholds
colnames(df_representation) <- study_list

for(i in 1:length(detection_thresholds))
{
	threshold <- detection_thresholds[i]
	print(threshold)
	for(j in 1:length(study_list))
	{
		study_name <- study_list[j]
		temp_core <- rownames(AllSpeciesDetectionPattern)[which(AllSpeciesDetectionPattern[,study_name]>=threshold)]
		if(length(temp_core)>1)
		{
			temp_sp_profile <- controlSpProfileAll[controlSpProfileAll$study_name == study_name,temp_core]
			temp_sp_profile <- temp_sp_profile[rowSums(temp_sp_profile)>0,colSums(temp_sp_profile)>0]
			temp_jaccard <- as.matrix(vegdist(temp_sp_profile,method="jaccard",binary=TRUE))
			diag(temp_jaccard) <- NA
			df_jaccard[i,j] <- median(apply(temp_jaccard,1,function(x)(x[!is.na(x)])))
			df_representation[i,j] <- median(rowSums(temp_sp_profile))
		}
		else
		{
			df_jaccard[i,j] <- 0
			df_representation[i,j] <- median(sum(temp_sp_profile))
		}
	}
}

df_patterns <- melt(1-df_jaccard)
colnames(df_patterns) <- c("Threshold","Study","Jaccard_Similarity")
df_patterns$Representation <- melt(df_representation)$value
colnames(df_patterns)[4] <- "Representation"
df_patterns <- df_patterns[df_patterns$Representation<=1,]

ggplot(df_patterns)+geom_boxplot(aes(y=Jaccard_Similarity,x=Threshold,group=Threshold),fill="green",alpha=0.5)+geom_boxplot(aes(y=Representation,x=Threshold,group=Threshold),fill="blue",alpha=0.5)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(size=0.5,color="grey70"))


rank_threshold <- seq(0.05,0.95,by=0.05)

df_sig_overlap <- matrix(0,length(rank_threshold),length(study_list))
rownames(df_sig_overlap) <- rank_threshold
colnames(df_sig_overlap) <- study_list

for(i in 1:length(rank_threshold))
{
	thres <- rank_threshold[i]
	temp_core_influence <- (apply(prevalentDf,2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df,2,function(x)ifelse(x>=thres,1,0)))
	temp_sig <- apply(prDf,2,function(x)(ifelse(x<=0.05,1,0)))
	for(j in 1:length(study_list))
	{
		#df_sig_overlap[i,j] <- 1-vegdist(rbind(temp_core_influence[,j],temp_sig[,j]),method="jaccard")[1]
		temp_r2_vec <- as.numeric(temp_core_influence[,j])
		temp_pr_vec <- as.numeric(temp_sig[,j])
		df_sig_overlap[i,j] <- (length(which((temp_pr_vec==1)&(temp_r2_vec==1))) + length(which((temp_pr_vec==0)&(temp_r2_vec==0))))/length(temp_r2_vec)
	}
}

df_shortlisted_thresholds <- as.data.frame(cbind(melt(df_numb_species)[which((melt(df_representation_90_plus)[,3]>=0.90)&(melt(df_representation_70_minus)[,3]==0)),],melt(df_representation_90_plus)[which((melt(df_representation_90_plus)[,3]>=0.90)&(melt(df_representation_70_minus)[,3]==0)),3],melt(df_representation_70_minus)[which((melt(df_representation_90_plus)[,3]>=0.90)&(melt(df_representation_70_minus)[,3]==0)),3]))

colnames(df_shortlisted_thresholds) <- c("X","Y","Number","90_plus","70_minus")

df_shortlisted_thresholds <- df_shortlisted_thresholds[order(df_shortlisted_thresholds[,3]),]
