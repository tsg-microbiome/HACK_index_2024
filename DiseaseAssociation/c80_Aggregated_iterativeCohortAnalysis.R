# Dr. Tarini Shankar Ghosh, Abhishek Goel
# 19-11-2023

# Input Files:
# 1) c70_discoveryCohortSampleInfo.xlsx
# 2) c70_validationCohortSampleInfo.xlsx
# 3) c10_species_to_check.RData
# 4) c71_finalData.RData
# 5) c56_imageCode_MajorCoreKeyStone.RData
# 6) c73_diseaseAnalysisFunction.R

# Output Files:
# 1) c80_iterativeAnalysis

# the code will perform the iterative analysis on the of the studies.
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

combinedCohorts= read_excel("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c70_discoveryCohortSampleInfo.xlsx")
validationCohort= read_excel("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c70_validationCohortSampleInfo.xlsx")

combinedCohorts= rbind(combinedCohorts,validationCohort)

load("G:/My Drive/IIITD/project/control_runs/input_files/c10_species_to_check.RData")
load("G:\\My Drive\\IIITD\\project\\control_runs\\datasets\\c71_finalData.RData")
load("G:/My Drive/IIITD/project/control_runs/output_files/c56_imageCode_MajorCoreKeyStone.RData")
source("G:/My Drive/IIITD/project/control_runs/codes/c73_diseaseAnalysisFunction.R")


# filtering out the disease__study combinations 
combinedCohorts= combinedCohorts %>% 
  filter(`80_20`==1) %>% 
  filter(controlcount>=15 & diseaseCount>=15)

# getting the control samples of the shortlisted studies for the discovery cohort.
controlSampleList= list()
uniqueStudies= unique(combinedCohorts$study)

# storing the control samples for all the studies in the list.
for(study in uniqueStudies)
{
  currStudy= AllCombinedMetadata %>% 
    filter(!metaType%in% c("batch4","batch3")) %>% 
    filter(study_name==study) %>% 
    filter(study_condition=="control")
  
  controlSampleList[[study]]= currStudy$sample_id
}

# expanding the metadata based on disease category
expandedMetadata= AllCombinedMetadata %>% 
  filter(!metaType %in% c("batch4","batch3")) %>% 
  separate_rows(diseaseCat, sep=";")

allDiseases= unique(combinedCohorts$disease)  

# it will contain the sample_ids for the disease__study combinations
diseaseSampleList= list()

for(dis in allDiseases)
{
  temp= combinedCohorts %>% 
    filter(disease==dis)
  
  uniqueStudies= unique(temp$study)
  
  for(study in uniqueStudies)
  {
    tempMeta= expandedMetadata %>% 
      filter(study_name==study) %>% 
      filter(diseaseCat==dis)
    
    diseaseSampleList[[dis]][[study]]= tempMeta$sample_id
  }
  
}

keyStoneType= c()
for(species in species_to_check)
{
  keyStoneType= c(keyStoneType,ifelse(species %in% MajorCoreKeyStone,"MajorCore","NonMajorCore"))
}

names(keyStoneType) <- species_to_check
speciesGroupingVariable= keyStoneType
speciesGroupingColour= c("#20B2AA","#DC143C")

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


# creating the correlation heatmap for the healt scores obtained in all the iterations.
iterativeScores= as.data.frame(read_excel("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c80_iterativeAnalysis\\iterationHealthScore_rankedScaledDifference.xlsx"))

rownames(iterativeScores)= iterativeScores$species

# removing the species name column from the df.
iterativeScores= iterativeScores[,-ncol(iterativeScores)]

# arranging the columns as per the iteration number in ascending order.
columns= sort(as.numeric(colnames(iterativeScores)))

iterativeScores= iterativeScores[,as.character(columns)]

corMatrix= as.data.frame(cor(iterativeScores, method = "spearman"))

library(ggplot2)
heatmap_union_core_df= corMatrix
heatmap_union_core_df$iteration1= c(rownames(corMatrix))

library(tidyverse)
union_df_2=heatmap_union_core_df %>%
  pivot_longer(!iteration1, names_to = "iteration2", values_to = "correlation")


union_df_2$iteration1 <- factor(union_df_2$iteration1, levels = unique(union_df_2$iteration1))
union_df_2$iteration2 <- factor(union_df_2$iteration2, levels = unique(union_df_2$iteration2))

library(ggplot2)
# color_vector <- c("#eddfb3", "#ca955c","#132b43","#56b1f7")
# color_vector <- c("white","#56b1f7")
ggplot(union_df_2, aes(x = iteration1, y = iteration2, fill = correlation))+
  geom_tile()+
  scale_fill_gradient(low= "yellow", high = "#d24531",limits= c(0,1)) +
  theme(axis.text.y = element_text(size = 10,color = "black"))+
  theme(axis.text.x= element_text(size=10, color="black"))+
  geom_tile(color = "black")

ggsave("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c81_allDiseaseOutput\\ranked_sacledDifference_score_correlation.pdf", height = 8, width= 12)

# getting the cor matrix in excel table as per the heatmap
corMatrix= corMatrix[sort(as.numeric(rownames(corMatrix)), decreasing = TRUE),]
corMatrix= round(corMatrix,2)
library(writexl)
write_xlsx(corMatrix, "G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c81_allDiseaseOutput\\corMatrix_iterationAnalysis.xlsx")
