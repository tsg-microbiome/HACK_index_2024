load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Meslier_V_2019_Metadata_SpProfile.RData")

subjectIDs <- unique(sapply(strsplit(rownames(Meslier_V_2019.species.matrix),"_"),function(x)(x[[1]])))

MED_subjectID <- unique(sapply(strsplit(grep("MED",rownames(Meslier_V_2019.species.matrix),value=TRUE),"_"),function(x)(x[[1]])))

CONTROL_subjectID <- unique(sapply(strsplit(grep("CONTROL",rownames(Meslier_V_2019.species.matrix),value=TRUE),"_"),function(x)(x[[1]])))

CONTROL_baseline <- paste0(CONTROL_subjectID,"_CONTROLbaseline")
CONTROL_weeks4 <- paste0(CONTROL_subjectID,"_CONTROLweeks4")
CONTROL_weeks8 <- paste0(CONTROL_subjectID,"_CONTROLweeks8")

MED_baseline <- paste0(MED_subjectID,"_MEDbaseline")
MED_weeks4 <- paste0(MED_subjectID,"_MEDweeks4")
MED_weeks8 <- paste0(MED_subjectID,"_MEDweeks8")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\SpeciesScores.RData")
source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\code_library.R")

SpeciesScores <- SpeciesScores[order(SpeciesScores[,4]),]
common_species <- intersect(rownames(SpeciesScores),colnames(Meslier_V_2019.species.matrix))

Meslier_V_2019_species <- Meslier_V_2019.species.matrix[,common_species]/rowSums(Meslier_V_2019.species.matrix[,common_species])

Meslier_V_2019_species_non_ranked <- Meslier_V_2019_species

Meslier_V_2019_species <- apply(Meslier_V_2019_species,2,rank_scale)

IncreasedMedDiet <- sub("weeks4","",c(MED_weeks4,CONTROL_weeks4)[which(Meslier_V_2019_diet_metadata[c(MED_weeks4,CONTROL_weeks4),"ItMedIndex"] - Meslier_V_2019_diet_metadata[c(MED_baseline,CONTROL_baseline),"ItMedIndex"]>0)])

DecreasedMedDiet <- sub("weeks4","",c(MED_weeks4,CONTROL_weeks4)[which(Meslier_V_2019_diet_metadata[c(MED_weeks4,CONTROL_weeks4),"ItMedIndex"] - Meslier_V_2019_diet_metadata[c(MED_baseline,CONTROL_baseline),"ItMedIndex"]<=0)])

diff_species_IncreasedMedDiet <- Meslier_V_2019_species[paste0(IncreasedMedDiet,"weeks4"),common_species] - Meslier_V_2019_species[paste0(IncreasedMedDiet,"baseline"),common_species]

diff_species_DecreasedMedDiet <- Meslier_V_2019_species[paste0(DecreasedMedDiet,"weeks4"),common_species] - Meslier_V_2019_species[paste0(DecreasedMedDiet,"baseline"),common_species]

t_wilcox_Increased_Decreased <- wilcox_batch(t(diff_species_IncreasedMedDiet),t(diff_species_DecreasedMedDiet))

diff_CONTROL_species <- Meslier_V_2019_species[paste0(CONTROL_subjectID,"_CONTROLweeks4"),]-Meslier_V_2019_species[paste0(CONTROL_subjectID,"_CONTROLbaseline"),]

diff_MED_species <- Meslier_V_2019_species[paste0(MED_subjectID,"_MEDweeks4"),]-Meslier_V_2019_species[paste0(MED_subjectID,"_MEDbaseline"),]

t_wilcox_CONTROL_MED <- wilcox_batch(t(diff_MED_species),t(diff_CONTROL_species))

diff_NonHACKs_MED_weeks4_baseline <- apply(Meslier_V_2019_species[MED_weeks4,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)-apply(Meslier_V_2019_species[MED_baseline,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)
diff_HACKs_MED_weeks4_baseline <- apply(Meslier_V_2019_species[MED_weeks4,rownames(tail(SpeciesScores),18)],1,mean)-apply(Meslier_V_2019_species[MED_baseline,rownames(tail(SpeciesScores),18)],1,mean)
diff_NonHACKs_CONTROL_weeks4_baseline <- apply(Meslier_V_2019_species[CONTROL_weeks4,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)-apply(Meslier_V_2019_species[CONTROL_baseline,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)
diff_HACKs_CONTROL_weeks4_baseline <- apply(Meslier_V_2019_species[CONTROL_weeks4,rownames(tail(SpeciesScores),18)],1,mean)-apply(Meslier_V_2019_species[CONTROL_baseline,rownames(tail(SpeciesScores),18)],1,mean)

diff_NonHACKs_MED_weeks8_weeks4 <- apply(Meslier_V_2019_species[MED_weeks8,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)-apply(Meslier_V_2019_species[MED_weeks4,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)
diff_HACKs_MED_weeks8_weeks4 <- apply(Meslier_V_2019_species[MED_weeks8,rownames(tail(SpeciesScores),18)],1,mean)-apply(Meslier_V_2019_species[MED_weeks4,rownames(tail(SpeciesScores),18)],1,mean)
diff_NonHACKs_CONTROL_weeks8_weeks4 <- apply(Meslier_V_2019_species[CONTROL_weeks8,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)-apply(Meslier_V_2019_species[CONTROL_weeks4,setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)
diff_HACKs_CONTROL_weeks8_weeks4 <- apply(Meslier_V_2019_species[CONTROL_weeks8,rownames(tail(SpeciesScores),18)],1,mean)-apply(Meslier_V_2019_species[CONTROL_weeks4,rownames(tail(SpeciesScores),18)],1,mean)

diff_NonHACKs_DecreasedMedDiet_weeks4_baseline <- apply(Meslier_V_2019_species[paste0(DecreasedMedDiet,"weeks4"),setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)-apply(Meslier_V_2019_species[paste0(DecreasedMedDiet,"baseline"),setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)
diff_HACKs_DecreasedMedDiet_weeks4_baseline <- apply(Meslier_V_2019_species[paste0(DecreasedMedDiet,"weeks4"),rownames(tail(SpeciesScores),18)],1,mean)-apply(Meslier_V_2019_species[paste0(DecreasedMedDiet,"baseline"),rownames(tail(SpeciesScores),18)],1,mean)
diff_NonHACKs_IncreasedMedDiet_weeks4_baseline <- apply(Meslier_V_2019_species[paste0(IncreasedMedDiet,"weeks4"),setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)-apply(Meslier_V_2019_species[paste0(IncreasedMedDiet,"baseline"),setdiff(common_species,rownames(tail(SpeciesScores),18))],1,mean)
diff_HACKs_IncreasedMedDiet_weeks4_baseline <- apply(Meslier_V_2019_species[paste0(IncreasedMedDiet,"weeks4"),rownames(tail(SpeciesScores),18)],1,mean)-apply(Meslier_V_2019_species[paste0(IncreasedMedDiet,"baseline"),rownames(tail(SpeciesScores),18)],1,mean)

HACK_Species <- rownames(tail(SpeciesScores,18))
NonHACK_Species <- setdiff(common_species,HACK_Species)

df_SpeciesGroups <- data.frame("HACK"=apply(Meslier_V_2019_species[,HACK_Species],1,mean),"NonHACK"=apply(Meslier_V_2019_species[,NonHACK_Species],1,mean))

Meslier_V_2019_diet_rankscaled <- as.data.frame(apply(Meslier_V_2019_diet_metadata[,c(6,8,10,12,14,15,16,17,19,22:25)]/Meslier_V_2019_diet_metadata[,"energy"],2,rank_scale))

#Meslier_V_2019_diet_rankscaled$sfa_pufa <- (Meslier_V_2019_diet_metadata$sfa_g+0.00001)/(Meslier_V_2019_diet_metadata$pufa_g+0.00001)

print("Level1")

corr_diet_species <- corr.test(Meslier_V_2019_species,Meslier_V_2019_diet_rankscaled[rownames(Meslier_V_2019_species),],method="spearman",use="pairwise.complete")

CorrDietSpecies <- as.data.frame(corr_diet_species$r)

CorrDietSpecies$HACK <- SpeciesScores[rownames(corr_diet_species$r),4]

CorrDietHACK <- corr.test(CorrDietSpecies,adjust="fdr")

df_CorrDietHACK <- data.frame("r"=CorrDietHACK$r[1:13,14],"p"=CorrDietHACK$p[1:13,14])

df_CorrDietHACK$Tags <- c("Carbohydrates_Overall","Sugars","Proteins","Lipids","Saturated Fatty Acid","Mono-Unsaturated Fatty Acids","Poly-Unsaturated Fatty Acids","Dietary Fiber","Alcohol","Vegetable Proteins","Animal Proteins","Veg Protein:Animal Protein","Fruits Vegetable") 

ggplot(df_CorrDietHACK,aes(x=r,y=-log(p,10)))+geom_point()+theme_bw()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_text_repel(label=ifelse(df_CorrDietHACK$p<=0.10,df_CorrDietHACK$Tags,""),color=ifelse(df_CorrDietHACK$r<0,"red","blue"))+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

diff_diet_MED_weeks4_baseline <- Meslier_V_2019_diet_rankscaled[c(MED_weeks4),]-Meslier_V_2019_diet_rankscaled[c(MED_baseline),]

diff_species_MED_weeks4_baseline <- Meslier_V_2019_species[c(MED_weeks4),]-Meslier_V_2019_species[c(MED_baseline),]

corr_MED_diff_diet_species <- corr.test(diff_species_MED_weeks4_baseline,diff_diet_MED_weeks4_baseline,use="pairwise.complete")

CorrMEDDiffDietSpecies <- as.data.frame(corr_MED_diff_diet_species$r)

CorrMEDDiffDietSpecies$HACK <- SpeciesScores[rownames(corr_MED_diff_diet_species$r),4]

CorrMEDDiffDietHACK <- corr.test(CorrMEDDiffDietSpecies,method="spearman",adjust="fdr")

df_CorrMEDDiffDietHACK <- data.frame("r"=CorrMEDDiffDietHACK$r[1:13,14],"p"=CorrMEDDiffDietHACK$p[1:13,14])

df_CorrMEDDiffDietHACK$Tags <- c("Carbohydrates_Overall","Sugars","Proteins","Lipids","Saturated Fatty Acid","Mono-Unsaturated Fatty Acids","Poly-Unsaturated Fatty Acids","Dietary Fiber","Alcohol","Vegetable Proteins","Animal Proteins","Veg Protein:Animal Protein","Fruits Vegetable") 

ggplot(df_CorrMEDDiffDietHACK,aes(x=r,y=-log(p,10)))+geom_point()+theme_bw()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_text_repel(label=ifelse(df_CorrMEDDiffDietHACK$p<=0.10,df_CorrMEDDiffDietHACK$Tags,""),color=ifelse(df_CorrMEDDiffDietHACK$r<0,"red","blue"))+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

#
diff_diet_CONTROL_weeks4_baseline <- Meslier_V_2019_diet_rankscaled[c(CONTROL_weeks4),]-Meslier_V_2019_diet_rankscaled[c(CONTROL_baseline),]

diff_species_CONTROL_weeks4_baseline <- Meslier_V_2019_species[c(CONTROL_weeks4),]-Meslier_V_2019_species[c(CONTROL_baseline),]

corr_CONTROL_diff_diet_species <- corr.test(diff_species_CONTROL_weeks4_baseline,diff_diet_CONTROL_weeks4_baseline,use="pairwise.complete")

CorrCONTROLDiffDietSpecies <- as.data.frame(corr_CONTROL_diff_diet_species$r)

CorrCONTROLDiffDietSpecies$HACK <- SpeciesScores[rownames(corr_CONTROL_diff_diet_species$r),4]

CorrCONTROLDiffDietHACK <- corr.test(CorrCONTROLDiffDietSpecies,method="spearman",adjust="fdr")

df_CorrCONTROLDiffDietHACK <- data.frame("r"=CorrCONTROLDiffDietHACK$r[1:13,14],"p"=CorrCONTROLDiffDietHACK$p[1:13,14])

df_CorrCONTROLDiffDietHACK$Tags <- c("Carbohydrates_Overall","Sugars","Proteins","Lipids","Saturated Fatty Acid","Mono-Unsaturated Fatty Acids","Poly-Unsaturated Fatty Acids","Dietary Fiber","Alcohol","Vegetable Proteins","Animal Proteins","Veg Protein:Animal Protein","Fruits Vegetable") 

ggplot(df_CorrCONTROLDiffDietHACK,aes(x=r,y=-log(p,10)))+geom_point()+theme_bw()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_text_repel(label=ifelse(df_CorrCONTROLDiffDietHACK$p<=0.10,df_CorrCONTROLDiffDietHACK$Tags,""),color=ifelse(df_CorrCONTROLDiffDietHACK$r<0,"red","blue"))+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))


CorrMedDietScoreSpecies <- corr.test(Meslier_V_2019_species,Meslier_V_2019_diet_metadata[rownames(Meslier_V_2019_species),"ItMedIndex"])

df_CorrMedDietScoreSpecies <- data.frame("r"=CorrMedDietScoreSpecies$r,"p"=CorrMedDietScoreSpecies$p)

df_CorrMedDietScoreSpecies$HACK <- SpeciesScores[rownames(df_CorrMedDietScoreSpecies),4]

ggplot(df_CorrMedDietScoreSpecies,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

diff_MEDAdherence <- Meslier_V_2019_diet_metadata[c(MED_weeks4),"ItMedIndex"]-Meslier_V_2019_diet_metadata[c(MED_baseline),"ItMedIndex"]

diff_CONTROLAdherence <- Meslier_V_2019_diet_metadata[c(CONTROL_weeks4),"ItMedIndex"]-Meslier_V_2019_diet_metadata[c(CONTROL_baseline),"ItMedIndex"]

corr_MED_diff_species_med_adherence <- corr.test(diff_species_MED_weeks4_baseline,diff_MEDAdherence,method="spearman",adjust="fdr")

df_CorrMEDDiffSpeciesMedAdherence <- data.frame("r"=corr_MED_diff_species_med_adherence$r,"p"=corr_MED_diff_species_med_adherence$p)

df_CorrMEDDiffSpeciesMedAdherence$HACK <- SpeciesScores[rownames(df_CorrMEDDiffSpeciesMedAdherence),4]

ggplot(df_CorrMEDDiffSpeciesMedAdherence,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

corr_CONTROL_diff_species_med_adherence <- corr.test(diff_species_CONTROL_weeks4_baseline,diff_CONTROLAdherence,method="spearman",adjust="fdr")

df_CorrCONTROLDiffSpeciesMedAdherence <- data.frame("r"=corr_CONTROL_diff_species_med_adherence$r,"p"=corr_CONTROL_diff_species_med_adherence$p)

df_CorrCONTROLDiffSpeciesMedAdherence$HACK <- SpeciesScores[rownames(df_CorrCONTROLDiffSpeciesMedAdherence),4]

ggplot(df_CorrCONTROLDiffSpeciesMedAdherence,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

boxplot(CorrDietSpecies[,"dietary_fiber_g"]~cut(CorrDietSpecies[,"HACK"],breaks=c(0,0.4,0.7,1),include.lowest=TRUE),ylab="Correlation Dietary Fiber",xlab="HACK Index",col=c("coral","cornsilk","cyan"),outline=FALSE)


HACKs <- rownames(tail(SpeciesScores,18))
NonHACKs <- setdiff(colnames(Meslier_V_2019_species),HACKs)

pcoDietChange <- dudi.pco(vegdist(as.matrix(rbind(diff_diet_MED_weeks4_baseline,diff_diet_CONTROL_weeks4_baseline)),method="euclidean"),scannf=FALSE,nf=3)

pcoSpeciesChange <- dudi.pco(vegdist(as.matrix(rbind(diff_species_MED_weeks4_baseline,diff_species_CONTROL_weeks4_baseline)),method="euclidean"),scannf=FALSE,nf=3)

pcoHACKSpeciesChange <- dudi.pco(vegdist(as.matrix(rbind(diff_species_MED_weeks4_baseline[,HACKs],diff_species_CONTROL_weeks4_baseline[,HACKs])),method="euclidean"),scannf=FALSE,nf=3)

pcoNonHACKSpeciesChange <- dudi.pco(vegdist(as.matrix(rbind(diff_species_MED_weeks4_baseline[,NonHACKs],diff_species_CONTROL_weeks4_baseline[,NonHACKs])),method="euclidean"),scannf=FALSE,nf=3)

diff_diet_weeks4_baseline <- as.data.frame(rbind(diff_diet_MED_weeks4_baseline,diff_diet_CONTROL_weeks4_baseline))

diff_species_weeks4_baseline <- as.data.frame(rbind(diff_species_MED_weeks4_baseline,diff_species_CONTROL_weeks4_baseline))

Meslier_V_2019_diet_metadata$Groups <- Meslier_V_2019_metadata[rownames(Meslier_V_2019_diet_metadata),"timepoint"]

boxplot(Meslier_V_2019_diet_metadata[,"ItMedIndex"]~Meslier_V_2019_diet_metadata[,"Groups"],col="coral2",cex.axis=1.5)

corr_Species_Med_CONTROL_baseline <- corr.test(Meslier_V_2019_species[CONTROL_baseline,intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species))],Meslier_V_2019_diet_metadata[CONTROL_baseline,"ItMedIndex"])

dfCorrSpeciesMED_Control_baseline <- data.frame(r=apply(corr_Species_Med_CONTROL_baseline$r,2,function(x)(ifelse(is.na(x),0,x))),HACK=SpeciesScores[intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species)),4])

ggplot(dfCorrSpeciesMED_Control_baseline,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm") + theme_bw() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

corr_Species_Med_CONTROL_weeks4 <- corr.test(Meslier_V_2019_species[CONTROL_weeks4,intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species))],Meslier_V_2019_diet_metadata[CONTROL_weeks4,"ItMedIndex"])

dfCorrSpeciesMED_Control_weeks4 <- data.frame(r=apply(corr_Species_Med_CONTROL_weeks4$r,2,function(x)(ifelse(is.na(x),0,x))),HACK=SpeciesScores[intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species)),4])

ggplot(dfCorrSpeciesMED_Control_weeks4,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm") + theme_bw() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

corr_Species_Med_CONTROL_weeks8 <- corr.test(Meslier_V_2019_species[CONTROL_weeks8,intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species))],Meslier_V_2019_diet_metadata[CONTROL_weeks8,"ItMedIndex"])

dfCorrSpeciesMED_Control_weeks8 <- data.frame(r=apply(corr_Species_Med_CONTROL_weeks8$r,2,function(x)(ifelse(is.na(x),0,x))),HACK=SpeciesScores[intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species)),4])

ggplot(dfCorrSpeciesMED_Control_weeks8,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm") + theme_bw() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

corr_Species_Med_MED_baseline <- corr.test(Meslier_V_2019_species[MED_baseline,intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species))],Meslier_V_2019_diet_metadata[MED_baseline,"ItMedIndex"])

dfCorrSpeciesMED_MED_baseline <- data.frame(r=apply(corr_Species_Med_MED_baseline$r,2,function(x)(ifelse(is.na(x),0,x))),HACK=SpeciesScores[intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species)),4])

ggplot(dfCorrSpeciesMED_MED_baseline,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm") + theme_bw() + theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))

corr_Species_Med_MED_weeks4 <- corr.test(Meslier_V_2019_species[MED_weeks4,intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species))],Meslier_V_2019_diet_metadata[MED_weeks4,"ItMedIndex"])

dfCorrSpeciesMED_MED_weeks4 <- data.frame(r=apply(corr_Species_Med_MED_weeks4$r,2,function(x)(ifelse(is.na(x),0,x))),HACK=SpeciesScores[intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species)),4])

ggplot(dfCorrSpeciesMED_MED_weeks4,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm") + theme_bw() + theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))

corr_Species_Med_MED_weeks8 <- corr.test(Meslier_V_2019_species[MED_weeks8,intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species))],Meslier_V_2019_diet_metadata[MED_weeks8,"ItMedIndex"])

dfCorrSpeciesMED_MED_weeks8 <- data.frame(r=apply(corr_Species_Med_MED_weeks8$r,2,function(x)(ifelse(is.na(x),0,x))),HACK=SpeciesScores[intersect(rownames(SpeciesScores),colnames(Meslier_V_2019_species)),4])

ggplot(dfCorrSpeciesMED_MED_weeks8,aes(x=r,y=HACK))+geom_point()+geom_smooth(method="lm") + theme_bw() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

CorrHACKDietChange <- corr.test(pcoDietChange$li[,c(1,2)],diff_species_weeks4_baseline[rownames(pcoDietChange$li),HACKs])

df_CorrHACKDietChange <- as.data.frame(t(CorrHACKDietChange$r))

ggplot(df_CorrHACKDietChange,aes(x=A1,y=A2))+geom_point()+geom_segment(x=0,y=0,xend=df_CorrHACKDietChange$A1,yend=df_CorrHACKDietChange$A2)+ theme_bw() + geom_hline(yintercept=0,color="lightslateblue",size=1) + geom_vline(xintercept=0,color="lightslateblue",size=1) + geom_text_repel(label=rownames(df_CorrHACKDietChange),color=ifelse((df_CorrHACKDietChange$A1<0)&(df_CorrHACKDietChange$A2<0),"blue",ifelse((df_CorrHACKDietChange$A1>0)&(df_CorrHACKDietChange$A2>0),"firebrick4","grey30"))) 

df_CorrHACKDietChange$PointType <- "Taxa"

Centroids <- as.data.frame(rbind(apply(pcoDietChange$li[MED_weeks4,c(1:2)],2,mean),apply(pcoDietChange$li[CONTROL_weeks4,c(1:2)],2,mean)))

rownames(Centroids) <- c("MED DIET","CONTROL DIET")

Centroids$PointType <- "Centroid"

CombinedHACKDietChange <- as.data.frame(rbind(df_CorrHACKDietChange,Centroids))

ggplot(CombinedHACKDietChange,aes(x=A1,y=A2))+geom_point(size=ifelse(CombinedHACKDietChange$PointType == "HACK",1,2))+geom_segment(x=0,y=0,xend=CombinedHACKDietChange$A1,yend=CombinedHACKDietChange$A2)+ theme_bw() + geom_hline(yintercept=0,color="lightslateblue",size=1) + geom_vline(xintercept=0,color="lightslateblue",size=1) + geom_text_repel(label=rownames(CombinedHACKDietChange),color=ifelse((CombinedHACKDietChange$A1<0)&(CombinedHACKDietChange$A2<0),"blue",ifelse((CombinedHACKDietChange$A1>0)&(CombinedHACKDietChange$A2>0),"firebrick4","grey30")),size=ifelse(CombinedHACKDietChange$PointType == "HACK",1,2)) 

#MedDiet-Associated Food Intake