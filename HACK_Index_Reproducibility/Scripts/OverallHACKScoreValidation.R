load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\SpeciesScores_NEW.RData")
load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\HealthScores.RData")
load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\StabilityScores.RData")
load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\CoreScores.RData")

df_Combined_SequencingType_WGS <- as.data.frame(matrix(0,201,3))
rownames(df_Combined_SequencingType_WGS) <- rownames(SpeciesScores_NEW)
colnames(df_Combined_SequencingType_WGS) <- c("CoreAssociation","StabilityAssociation","HealthAssociation")

df_Combined_SequencingType_WGS[rownames(df_CoreScore_SequencingType),"CoreAssociation"] <- df_CoreScore_SequencingType[,"WGS"]
df_Combined_SequencingType_WGS[rownames(df_StabilityScore_SequencingType),"StabilityAssociation"] <- df_StabilityScore_SequencingType[,"WGS"]
df_Combined_SequencingType_WGS[rownames(df_HealthScore_SequencingType),"HealthAssociation"] <- df_HealthScore_SequencingType[,"WGS"]
df_Combined_SequencingType_WGS$HACKScore <- ifelse(is.nan(apply(apply(df_Combined_SequencingType_WGS,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x))))),0,apply(apply(df_Combined_SequencingType_WGS,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x)))))


df_Combined_SequencingType_16S <- as.data.frame(matrix(0,201,3))
rownames(df_Combined_SequencingType_16S) <- rownames(SpeciesScores_NEW)
colnames(df_Combined_SequencingType_16S) <- c("CoreAssociation","StabilityAssociation","HealthAssociation")

df_Combined_SequencingType_16S[rownames(df_CoreScore_SequencingType),"CoreAssociation"] <- df_CoreScore_SequencingType[,"X16S"]
df_Combined_SequencingType_16S[rownames(df_StabilityScore_SequencingType),"StabilityAssociation"] <- df_StabilityScore_SequencingType[,"X16S"]
df_Combined_SequencingType_16S[rownames(df_HealthScore_SequencingType),"HealthAssociation"] <- df_HealthScore_SequencingType[,"X16S"]
df_Combined_SequencingType_16S$HACKScore <- ifelse(is.nan(apply(apply(df_Combined_SequencingType_16S,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x))))),0,apply(apply(df_Combined_SequencingType_16S,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x)))))

df_Combined_CohortType_IndustrializedUrban <- as.data.frame(matrix(0,201,3))
rownames(df_Combined_CohortType_IndustrializedUrban) <- rownames(SpeciesScores_NEW)
colnames(df_Combined_CohortType_IndustrializedUrban) <- c("CoreAssociation","StabilityAssociation","HealthAssociation")

df_Combined_CohortType_IndustrializedUrban[rownames(df_CoreScore_CohortType),"CoreAssociation"] <- df_CoreScore_CohortType[,"IndustrializedUrban"]
df_Combined_CohortType_IndustrializedUrban[rownames(df_StabilityScore_CohortType),"StabilityAssociation"] <- df_StabilityScore_CohortType[,"IndustrializedUrban"]
df_Combined_CohortType_IndustrializedUrban[rownames(df_HealthScore_CohortType),"HealthAssociation"] <- df_HealthScore_CohortType[,"IndustrializedUrban"]
df_Combined_CohortType_IndustrializedUrban$HACKScore <- ifelse(is.nan(apply(apply(df_Combined_CohortType_IndustrializedUrban,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x))))),0,apply(apply(df_Combined_CohortType_IndustrializedUrban,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x)))))


df_Combined_CohortType_Others <- as.data.frame(matrix(0,201,3))
rownames(df_Combined_CohortType_Others) <- rownames(SpeciesScores_NEW)
colnames(df_Combined_CohortType_Others) <- c("CoreAssociation","StabilityAssociation","HealthAssociation")

df_Combined_CohortType_Others[rownames(df_CoreScore_CohortType),"CoreAssociation"] <- (df_CoreScore_CohortType[,"UrbanRuralMixed"]+df_CoreScore_CohortType[,"RuralTribal"])/2
df_Combined_CohortType_Others[rownames(df_StabilityScore_CohortType),"StabilityAssociation"] <- df_StabilityScore_CohortType[,"Others"]
df_Combined_CohortType_Others[rownames(df_HealthScore_CohortType),"HealthAssociation"] <- df_HealthScore_CohortType[,"UrbanRuralMixed"]
df_Combined_CohortType_Others$HACKScore <- ifelse(is.nan(apply(apply(df_Combined_CohortType_Others,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x))))),0,apply(apply(df_Combined_CohortType_Others,2,rank_scale),1,function(x)(mean(x)*(1-Gini(x)))))

#df_HACKScore <- data.frame("WGS"=
