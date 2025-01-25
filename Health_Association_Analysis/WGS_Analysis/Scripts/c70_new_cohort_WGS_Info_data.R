## Replicating c70_ValidationDiseaseSampleInfo.R for cohort_WGS cohort

# Omprakash
# 25-09-2024

# the code will get the information of all the diseases and their control and disease samples for cohort_WGS cohort


## Load the metadata_WGS 
c1_new_cohort_wise_master_data <- new.env()
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/c1_new_cohort_wise_master_data.RData",envir = c1_new_cohort_wise_master_data)
attach(c1_new_cohort_wise_master_data)
metadata_WGS <- metadata_WGS
detach(c1_new_cohort_wise_master_data)
rm(c1_new_cohort_wise_master_data)

library(dplyr)
library(tidyverse)
library(writexl)

expandedMetadata= metadata_WGS %>% 
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
cohort_WGSCohortInfo= as.data.frame(matrix(nrow= rowCount, ncol= 4))
colnames(cohort_WGSCohortInfo)= c("disease", "study", "controlcount", "diseaseCount")

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
    controlMeta= metadata_WGS %>% 
      filter(study_name==study) %>% 
      filter(study_condition== "control")
    
    controlCount= nrow(controlMeta)
    
    # disease sample count in the specific study
    diseaseCount= nrow(temp[temp$study_name==study,])
    
    currRow= c(currRow,controlCount)
    currRow= c(currRow, diseaseCount)
    
    cohort_WGSCohortInfo[cnt,]= currRow
    cnt= cnt+1
  }
  
}

cohort_WGSCohortInfo$controlcount= as.numeric(cohort_WGSCohortInfo$controlcount)
cohort_WGSCohortInfo$diseaseCount= as.numeric(cohort_WGSCohortInfo$diseaseCount)

# getting only those disease_study combinations that has the corresponding control match samples.
cohort_WGSCohortInfo= cohort_WGSCohortInfo %>% 
  filter(controlcount>0) %>% 
  arrange(controlcount)

cohort_WGSCohortInfo= cohort_WGSCohortInfo %>% 
  arrange(disease)

# creating columns that will store the infromation about the percentage of control and disease samples for the given combination. 
# It is done on the disease__study combinations instead of the all the samples for the given disease condition from all the studies having that disease becase
# we are considering control matched samples for the disease in each study. If the given study has significnatly large number control and disease samples then 
# that may generate the bias. Therefore, for each study there should be enough disease and control matched samples independently.

# here the 60,70,80 mean that the control sample percent or disease sample percent for the particular disease__study combination doesn't excede these numbers.
# simply neither of them can constitute 61%,71%, 81% of the total number of samples.
cohort_WGSCohortInfo$`60_40`= NA
cohort_WGSCohortInfo$`70_30`= NA
cohort_WGSCohortInfo$`80_20`= NA

for(index in 1:nrow(cohort_WGSCohortInfo))
{
  total= cohort_WGSCohortInfo[index,"controlcount"]  + cohort_WGSCohortInfo[index, "diseaseCount"]
  value= cohort_WGSCohortInfo[index,"controlcount"]/total
  
  if(value<0.4 | value>0.6)
  {
    cohort_WGSCohortInfo[index,"60_40"]=0
  }
  else
  {
    cohort_WGSCohortInfo[index,"60_40"]=1
  }
  
  
  if(value<0.3 | value>0.7)
  {
    cohort_WGSCohortInfo[index,"70_30"]=0
  }
  
  else
  {
    cohort_WGSCohortInfo[index,"70_30"]=1
  }
  
  if(value<0.2 | value>0.8)
  {
    cohort_WGSCohortInfo[index,"80_20"]=0
  }
  
  else
  {
    cohort_WGSCohortInfo[index,"80_20"]=1
  }
  
  
}

write_xlsx(cohort_WGSCohortInfo, "C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/Shotgun_sequenced/c70_new_cohort_WGS_Info.xlsx")


save.image("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/Shotgun_sequenced/c70_new_cohort_WGS_Info_data.RData")
