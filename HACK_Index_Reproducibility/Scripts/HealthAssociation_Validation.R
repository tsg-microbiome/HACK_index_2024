load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\SpeciesScores_NEW.RData")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\CellReportsRevision\\Section3_DiseaseAssociationInvestigationFigures\\HealthScore_SubGroups\\HealthScore_16s_WGS.RData")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\CellReportsRevision\\Section3_DiseaseAssociationInvestigationFigures\\HealthScore_SubGroups\\HealthScore_UrbanRural_IndustrializedUrban.RData")

hack_top_17 <- rownames(SpeciesScores_NEW[SpeciesScores_NEW[,4]>=0.75,])

df_HealthScore_SequencingType = data.frame("WGS"=HealthScore_WGS[intersect(rownames(HealthScore_WGS),rownames(HealthScore_16s)),1],"16S"=HealthScore_16s[intersect(rownames(HealthScore_WGS),rownames(HealthScore_16s)),1],row.names=intersect(rownames(HealthScore_WGS),rownames(HealthScore_16s)))

df_HealthScore_SequencingType$CV <- apply(df_HealthScore_SequencingType,1,function(x)(sd(x)/mean(x)))

df_HealthScore_SequencingType[rownames(df_HealthScore_SequencingType),"Overall"] <- SpeciesScores_NEW[rownames(df_HealthScore_SequencingType),3]

ggplot(df_HealthScore_SequencingType,aes(x=Overall,y=CV))+geom_point(col=ifelse(rownames(df_HealthScore_SequencingType)%in%hack_top_17,"blue","grey50"),size=5)+theme_bw()+geom_smooth(method="lm",color="firebrick4",linewidth=2)+xlab("")+ylab("")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

df_HealthScore_CohortType = data.frame("IndustrializedUrban"=HealthScore_IndustrializedUrban[intersect(rownames(HealthScore_IndustrializedUrban),rownames(HealthScore_UrbanRuralMixed)),1],"UrbanRuralMixed"=HealthScore_UrbanRuralMixed[intersect(rownames(HealthScore_IndustrializedUrban),rownames(HealthScore_UrbanRuralMixed)),1],row.names=intersect(rownames(HealthScore_IndustrializedUrban),rownames(HealthScore_UrbanRuralMixed)))

df_HealthScore_SequencingType$CV <- apply(df_HealthScore_SequencingType,1,function(x)(sd(x)/mean(x)))

df_HealthScore_SequencingType[rownames(df_HealthScore_SequencingType),"Overall"] <- SpeciesScores_NEW[rownames(df_HealthScore_SequencingType),3]