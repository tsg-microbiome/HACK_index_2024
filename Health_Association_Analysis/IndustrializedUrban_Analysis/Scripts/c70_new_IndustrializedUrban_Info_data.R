## Replicating c70_ValidationDiseaseSampleInfo.R for IndustrializedUrban cohort

# Omprakash
# 23-09-2024

# the code will get the information of all the diseases and their control and disease samples for IndustrializedUrban cohort


## Load the IndustrializedUrban_metadata and IndustrializedUrban_spProfile 
c1_new_cohort_wise_master_data <- new.env()
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/c1_new_cohort_wise_master_data.RData",envir = c1_new_cohort_wise_master_data)
attach(c1_new_cohort_wise_master_data)
IndustrializedUrban_metadata <- IndustrializedUrban_metadata
detach(c1_new_cohort_wise_master_data)
rm(c1_new_cohort_wise_master_data)

library(dplyr)
library(tidyverse)
library(writexl)

expandedMetadata= IndustrializedUrban_metadata %>% 
  filter(study_condition=="disease") %>% 
  separate_rows(diseaseCat, sep=";") %>% 
  filter(diseaseCat!="other")

uniqueDiseases = unique(expandedMetadata$diseaseCat)

# getting total rows for the disease study combinations
rowCount=0

for(disease in uniqueDiseases)
{
  # getting all the samples for the given disease.
  temp= expandedMetadata[expandedMetadata$diseaseCat==disease,]
  totalStudies= length(unique(temp$study_name))
  rowCount= rowCount+totalStudies
}

# creating a empty dataframe that will store control and disease sample count for each of the disease__study combination.
IndustrializedUrbanCohortInfo= as.data.frame(matrix(nrow= rowCount, ncol= 4))
colnames(IndustrializedUrbanCohortInfo)= c("disease", "study", "controlcount", "diseaseCount")

cnt= 1
for(disease in uniqueDiseases)
{
  # getting disease specific samples
  temp= expandedMetadata[expandedMetadata$diseaseCat==disease,]
  
  # getting the studies that have the given disease
  uniqueStudies= unique(temp$study_name)
  
  for(study in uniqueStudies)
  {
    # creating a vector that will contain the information about the current disease__study combination.
    currRow= c(disease, study)
    
    # getting the number of non-diseases (control) samples for the given study.
    controlMeta= IndustrializedUrban_metadata %>% 
      filter(study_name==study) %>% 
      filter(study_condition== "control")
    
    controlCount= nrow(controlMeta)
    
    # disease sample count in the specific study
    diseaseCount= nrow(temp[temp$study_name==study,])
    
    currRow= c(currRow,controlCount)
    currRow= c(currRow, diseaseCount)
    
    IndustrializedUrbanCohortInfo[cnt,]= currRow
    cnt= cnt+1
  }
  
}

IndustrializedUrbanCohortInfo$controlcount= as.numeric(IndustrializedUrbanCohortInfo$controlcount)
IndustrializedUrbanCohortInfo$diseaseCount= as.numeric(IndustrializedUrbanCohortInfo$diseaseCount)

# getting only those disease_study combinations that has the corresponding control match samples.
IndustrializedUrbanCohortInfo= IndustrializedUrbanCohortInfo %>% 
  filter(controlcount>0) %>% 
  arrange(controlcount)

IndustrializedUrbanCohortInfo= IndustrializedUrbanCohortInfo %>% 
  arrange(disease)

# creating columns that will store the infromation about the percentage of control and disease samples for the given combination. 
# It is done on the disease__study combinations instead of the all the samples for the given disease condition from all the studies having that disease becase
# we are considering control matched samples for the disease in each study. If the given study has significnatly large number control and disease samples then 
# that may generate the bias. Therefore, for each study there should be enough disease and control matched samples independently.

# here the 60,70,80 mean that the control sample percent or disease sample percent for the particular disease__study combination doesn't excede these numbers.
# simply neither of them can constitute 61%,71%, 81% of the total number of samples.
IndustrializedUrbanCohortInfo$`60_40`= NA
IndustrializedUrbanCohortInfo$`70_30`= NA
IndustrializedUrbanCohortInfo$`80_20`= NA

for(index in 1:nrow(IndustrializedUrbanCohortInfo))
{
  total= IndustrializedUrbanCohortInfo[index,"controlcount"]  + IndustrializedUrbanCohortInfo[index, "diseaseCount"]
  value= IndustrializedUrbanCohortInfo[index,"controlcount"]/total
  
  if(value<0.4 | value>0.6)
  {
    IndustrializedUrbanCohortInfo[index,"60_40"]=0
  }
  else
  {
    IndustrializedUrbanCohortInfo[index,"60_40"]=1
  }
  
  
  if(value<0.3 | value>0.7)
  {
    IndustrializedUrbanCohortInfo[index,"70_30"]=0
  }
  
  else
  {
    IndustrializedUrbanCohortInfo[index,"70_30"]=1
  }
  
  if(value<0.2 | value>0.8)
  {
    IndustrializedUrbanCohortInfo[index,"80_20"]=0
  }
  
  else
  {
    IndustrializedUrbanCohortInfo[index,"80_20"]=1
  }
  
  
}

write_xlsx(IndustrializedUrbanCohortInfo, "C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/IndustrializedUrban/c70_new_IndustrializedUrban_Info.xlsx")


save.image("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/IndustrializedUrban/c70_new_IndustrializedUrban_Info_data.RData")
