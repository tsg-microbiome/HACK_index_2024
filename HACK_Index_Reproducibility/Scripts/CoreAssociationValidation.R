source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\code_library.R")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\HACK_Index_Reproducibility\\CoreInfluencers_3Dfs.RData")

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\SpeciesScores_NEW.RData")

StudyList_IndustrializedUrban <- c("AG", "Arivale", "AsnicarF_2021", "BackhedF_2015", "Barton_et_al_2018", "BedarfJR_2017", "Bengtsson-PalmeJ_2015", "CosteaPI_2017", "Cronin_et_al_2018", "DeFilippisF_2019", "FengQ_2015", "FerrettiP_2018", "Ghosh_et_al_2020", "HMP_2012", "HMP_2019_t2d", "HallAB_2017", "HanniganGD_2017", "HansenLBS_2018", "IjazUZ_2017", "Jeffery_et_al_2020", "JieZ_2017", "KarlssonFH_2013", "LiJ_2014", "LiJ_2017", "LifeLinesDeep_2016", "MehtaRS_2018", "MetaCardis_2020_a", "NagySzakalD_2017", "NielsenHB_2014", "Odamaki", "QinJ_2012", "QinN_2014", "SankaranarayananK_2015", "SchirmerM_2016", "ShaoY_2019", "ThomasAM_2018a", "ThomasAM_2018b", "ThomasAM_2019_c", "VogtmannE_2016", "WampachL_2018", "WirbelJ_2018", "XieH_2016", "YachidaS_2019", "YassourM_2016", "YassourM_2018", "YeZ_2018", "YuJ_2015", "ZeeviD_2015", "ZellerG_2014", "ZhuF_2020")

StudyList_UrbanRuralMixed <- c("BritoIL_2016", "CuestaZuluaga_et_al_2018", "DhakanDB_2019", "He", "Kedia_et_al", "KeohaneDM_2020", "LogMPie", "MicroDiab_India", "MobegiF_2020", "PehrssonE_2016", "RenallN_2023")

StudyList_RuralTribal <- c("KaurK_2020", "LiuW_2016", "LokmerA_2019", "Obregon-TitoAJ_2015", "PasolliE_2019", "RampelliS_2015", "RubelMA_2020", "SmitsSA_2017", "TettAJ_2019_a", "TettAJ_2019_b", "TettAJ_2019_c")

StudyList_WGS <- c("AsnicarF_2021", "BackhedF_2015", "Barton_et_al_2018", "BedarfJR_2017", "Bengtsson-PalmeJ_2015", "BritoIL_2016", "CosteaPI_2017", "Cronin_et_al_2018", "DeFilippisF_2019", "DhakanDB_2019", "FengQ_2015", "FerrettiP_2018", "Ghosh_et_al_2020", "HMP_2012", "HMP_2019_t2d", "HallAB_2017", "HanniganGD_2017", "HansenLBS_2018", "IjazUZ_2017", "Jeffery_et_al_2020", "JieZ_2017", "KarlssonFH_2013", "KaurK_2020", "KeohaneDM_2020", "LiJ_2014", "LiJ_2017", "LifeLinesDeep_2016", "LiuW_2016", "LokmerA_2019", "MehtaRS_2018", "MetaCardis_2020_a", "MobegiF_2020", "NagySzakalD_2017", "NielsenHB_2014", "Obregon-TitoAJ_2015", "PasolliE_2019", "PehrssonE_2016", "QinJ_2012", "QinN_2014", "RampelliS_2015", "RenallN_2023", "RubelMA_2020", "SankaranarayananK_2015", "SchirmerM_2016", "ShaoY_2019", "SmitsSA_2017", "TettAJ_2019_a", "TettAJ_2019_b", "TettAJ_2019_c", "ThomasAM_2018a", "ThomasAM_2018b", "ThomasAM_2019_c", "VogtmannE_2016", "WampachL_2018", "WirbelJ_2018", "XieH_2016", "YachidaS_2019", "YassourM_2016", "YassourM_2018", "YeZ_2018", "YuJ_2015", "ZeeviD_2015", "ZellerG_2014", "ZhuF_2020")

StudyList_16S <- c("AG", "Arivale", "CuestaZuluaga_et_al_2018", "He", "Kedia_et_al", "LogMPie", "MicroDiab_India", "Odamaki")

hack_top_17 <- rownames(SpeciesScores_NEW[SpeciesScores_NEW[,4]>=0.75,])

print("SequenceType Analysis")

SelectSpecies_16S <- names(which(apply(prevalentDf_new[,StudyList_16S],1,function(x)(length(x[x>0.05])))>=0.50))

SelectSpecies_WGS <- names(which(apply(prevalentDf_new[,StudyList_WGS],1,function(x)(length(x[x>0.05])))>=0.50))

SelectSpecies <- names(which(table(c(SelectSpecies_16S,SelectSpecies_WGS))==2))

df_SequencingType = data.frame("WGS"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,StudyList_WGS],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,StudyList_WGS],2,function(x)(ifelse(x>=0.70,1,0))))),"16S"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,StudyList_16S],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,StudyList_16S],2,function(x)(ifelse(x>=0.70,1,0))))),"Overall"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,],2,function(x)(ifelse(x>=0.70,1,0))))))

df_SequencingType$CV = apply(df_SequencingType[,1:2],1,function(x)(sd(x)/mean(x)))

df_SequencingType <- df_SequencingType[((df_SequencingType[,1]>0)&(df_SequencingType[,2]>0)),]

ggplot(df_SequencingType,aes(x=Overall,y=CV))+geom_point(col=ifelse(rownames(df_SequencingType)%in%hack_top_17,"blue","grey50"),size=5)+theme_bw()+geom_smooth(method="lm",color="firebrick4",linewidth=2)+xlab("")+ylab("")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

boxplot(df_SequencingType[intersect(rownames(df_SequencingType),hack_top_17),1],df_SequencingType[setdiff(rownames(df_SequencingType),hack_top_17),1],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

boxplot(df_SequencingType[intersect(rownames(df_SequencingType),hack_top_17),2],df_SequencingType[setdiff(rownames(df_SequencingType),hack_top_17),2],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

print("CohortTypeAnalysis")

SelectSpecies_IndustrializedUrban <- names(which(apply(prevalentDf_new[,StudyList_IndustrializedUrban],1,function(x)(length(x[x>0.05])))>=0.50))

SelectSpecies_UrbanRuralMixed <- names(which(apply(prevalentDf_new[,StudyList_UrbanRuralMixed],1,function(x)(length(x[x>0.05])))>=0.50))

SelectSpecies_RuralTribal <- names(which(apply(prevalentDf_new[,StudyList_RuralTribal],1,function(x)(length(x[x>0.05])))>=0.50))

SelectSpecies <- names(which(table(c(SelectSpecies_RuralTribal,SelectSpecies_UrbanRuralMixed,SelectSpecies_IndustrializedUrban))==3))

df_CohortType = data.frame("IndustrializedUrban"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,StudyList_IndustrializedUrban],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,StudyList_IndustrializedUrban],2,function(x)(ifelse(x>=0.70,1,0))))),"UrbanRuralMixed"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,StudyList_UrbanRuralMixed],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,StudyList_UrbanRuralMixed],2,function(x)(ifelse(x>=0.70,1,0))))),"RuralTribal"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,StudyList_RuralTribal],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,StudyList_RuralTribal],2,function(x)(ifelse(x>=0.70,1,0))))),"Overall"=rank_scale(rowSums(apply(prevalentDf_new[SelectSpecies,],2,function(x)(ifelse(x>=0.65,1,0))) * apply(r2Df_new[SelectSpecies,],2,function(x)(ifelse(x>=0.70,1,0))))))

df_CohortType$CV = apply(df_CohortType[,1:3],1,function(x)(sd(x)/mean(x)))

#df_CohortType <- df_CohortType[names(which(apply(df_CohortType[,1:3],1,function(x)(length(x[x>0])))>2)),]

ggplot(df_CohortType,aes(x=Overall,y=CV))+geom_point(col=ifelse(rownames(df_CohortType)%in%hack_top_17,"blue","grey50"),size=5)+theme_bw()+geom_smooth(method="lm",color="firebrick4",linewidth=2)+xlab("")+ylab("")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

boxplot(df_CohortType[intersect(rownames(df_CohortType),hack_top_17),1],df_CohortType[setdiff(rownames(df_CohortType),hack_top_17),1],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

boxplot(df_CohortType[intersect(rownames(df_CohortType),hack_top_17),2],df_CohortType[setdiff(rownames(df_CohortType),hack_top_17),2],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)

boxplot(df_CohortType[intersect(rownames(df_CohortType),hack_top_17),3],df_CohortType[setdiff(rownames(df_CohortType),hack_top_17),3],range=1,col=c("cadetblue1","antiquewhite1"),outline=FALSE,cex.axis=1.2)