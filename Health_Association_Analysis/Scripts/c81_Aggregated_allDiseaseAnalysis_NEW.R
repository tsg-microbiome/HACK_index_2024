### 29-09-2024 (Sourav Goswami) Disease Analysis ###

library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gplots)
library(vegan)
library(ade4)
library(adegraphics)
library(writexl)

load("G:/My Drive/HACK paper/CellReports_revision/disease_analysis_revision.RData")
source("G:/My Drive/HACK paper/CellReports_revision/disease_analysis_overall_output/c73_diseaseAnalysisFunction_NEW.R")
load("G:/My Drive/HACK paper/CellReports_revision/disease_analysis_overall_output/MajorCoreKeyStone.RData")
load("G:/My Drive/HACK paper/CellReports_revision/disease_analysis_overall_output/DiseaseAnalysis_overall_iterative_analysis_SG_OS.RData")

species_to_check <- NewGutAssociatedSpecies

outputFolder= "G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output"

intermediateOutput_NEW = targetIntermediateOutput(controlSampleList,diseaseSampleList,disease_analysis_spProfile, unique(combinedCohorts$disease), species_to_check, outputFolder, "mannWhitneyIntermediateOutput.RData", TRUE, speciesGroupList= list())

load("G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\mannWhitneyIntermediateOutput.RData")


keyStoneType= c()
for(species in species_to_check)
{
  keyStoneType= c(keyStoneType,ifelse(species %in% MajorCoreKeyStone,"MajorCore","NonMajorCore"))
}

names(keyStoneType) <- species_to_check
speciesGroupingVariable= keyStoneType
speciesGroupingColour= c("#20B2AA","#DC143C")


associationOutput(intermediateOutput_NEW, unique(combinedCohorts$disease), outputFolder, speciesGroupingVariable, speciesGroupingColour, TRUE, TRUE)


MainOutputFolder= "G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c80_iterativeAnalysis"

iterativeOutput= diseaseAnalysisIterationPipeline(10, allDiseases, speciesGroupingVariable, speciesGroupingColour, controlSampleList, diseaseSampleList, MainOutputFolder,AllCombinedSpProfile, species_to_check, speciesGroupList= list())

save(iterativeOutput, file= paste0(MainOutputFolder,"\\iterativeOutput.RData"))
iterationHealthScore_rankedScaledDifference= as.data.frame(iterativeOutput[["iterationOutputDf_ranked_scaledDifference"]])

# taking the mean of HACKScores from all the 10 iterations.
HealthScore= apply(iterationHealthScore_rankedScaledDifference,1,mean)

# save(HealthScore, file= "G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c81_allDiseaseOutput\\healthAssociationScore.RData")
save(HealthScore, file= "G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c80_iterativeAnalysis\\healthAssociationScore.RData")


iterationHealthScore_rankedScaledDifference$species= rownames(iterationHealthScore_rankedScaledDifference)

iterationHealthScore_rankedScaledDifference= iterationHealthScore_rankedScaledDifference %>% 
  arrange(iterationHealthScore_rankedScaledDifference$species)

write_xlsx(iterationHealthScore_rankedScaledDifference, paste0(MainOutputFolder,"\\iterationHealthScore_rankedScaledDifference.xlsx"))



diseaseAssociationSummary= read_excel("G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\heatmap2Df.xlsx")

dotPlotDf= diseaseAssociationSummary[,c("scaledDifference", "ranked_scaledDifference", "species","marker","taxaGroup")]

dotPlotDf$category= NA
dotPlotDf[which(dotPlotDf$scaledDifference<0),"category"]=-1
dotPlotDf[which(dotPlotDf$scaledDifference>0),"category"]= 1
dotPlotDf[is.na(dotPlotDf$category),"category"]=0

dotPlotDf= dotPlotDf %>% 
  arrange(ranked_scaledDifference)

dotPlotDf$s_no= c(1:201)

dotPlotDf$species <- factor(dotPlotDf$species, levels = unique(dotPlotDf$species))
dotPlotDf$category= as.factor(dotPlotDf$category)

pdf("G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\ranked_sacledDifference_dotPlot.pdf", height = 10, width= 15)
ggplot(dotPlotDf, aes(x=s_no,
                      y= ranked_scaledDifference, color= category))+
  geom_point(aes(alpha= 0.1, size= 0.2))+
  scale_color_manual(values = c("#c90101","#705500","#018d00","#126560","#c90101"))+
  theme_bw()+
  xlab("")+
  ylab("Health Score")+
  theme(legend.position = "none")+
  theme(axis.text = element_text(color="black"))+
  geom_label_repel(data= filter(filter(dotPlotDf,marker=="yes")),aes(label= species, color= taxaGroup),
                   size= 3.5, segment.size= 0.1, force= 80, max.overlaps = 83)
dev.off()

heatmap2_df_significant <- heatmap2_df[!is.na(heatmap2_df$diseaseAssociation), ]

Association_heatmap_colnames$indicator <- ifelse(Association_heatmap_colnames$x %in% heatmap2_df_significant$species, "YES", "NO")
colnames(Association_heatmap_colnames)[1] <- "species"
Association_heatmap_colnames$association <- heatmap2_df_significant[match(Association_heatmap_colnames$species, heatmap2_df_significant$species),18]
Association_heatmap_colnames_significant <- Association_heatmap_colnames[!is.na(Association_heatmap_colnames$association), ]

Association_heatmap_rownames <- read.table("AssociationHeatmap_rownames.txt", header = T)

write.xlsx(Association_heatmap_colnames,"Association_heatmap_colnames_modified.xlsx",row.names = F)
write.xlsx(Association_heatmap_colnames_significant,"Association_heatmap_colnames_significant.xlsx",row.names = F)



