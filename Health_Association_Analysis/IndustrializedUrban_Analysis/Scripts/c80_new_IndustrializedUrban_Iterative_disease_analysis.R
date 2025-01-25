######### This code is a replication of c80_Aggregated_iterativeCohortAnalysis.R with some modifications for IndustrializedUrban
## Omprakash 
## Dt.25-09-2024


## This code will use c73 for IndustializedUrban cohort and run disease analysis for this cohort

############################################################################## Load the data
## Load the IndustrializedUrban_metadata and IndustrializedUrban_spProfile 
c1_new_cohort_wise_master_data <- new.env()
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/c1_new_cohort_wise_master_data.RData",envir = c1_new_cohort_wise_master_data)
attach(c1_new_cohort_wise_master_data)
IndustrializedUrban_metadata <- IndustrializedUrban_metadata
IndustrializedUrban_spProfile <- IndustrializedUrban_spProfile
detach(c1_new_cohort_wise_master_data)
rm(c1_new_cohort_wise_master_data)

## Load the functions that are used in disease analysis (c73_new_IndustrializedUrban_disease_analysis_functions.R)
source("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/IndustrializedUrban/c73_new_IndustrializedUrban_disease_analysis_functions.R")
## Load the species to check in the data 
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/NewGutAssociatedSpecies.RData")

## Load the summarised information for IndustrializedUrban cohort from c70 (c70_new_IndustrializedUrban_Info_data.RData)
c70_new_IndustrializedUrban_Info_data <- new.env()
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/IndustrializedUrban/c70_new_IndustrializedUrban_Info_data.RData",envir = c70_new_IndustrializedUrban_Info_data)
attach(c70_new_IndustrializedUrban_Info_data)
IndustrializedUrbanCohortInfo <- IndustrializedUrbanCohortInfo
detach(c70_new_IndustrializedUrban_Info_data)
rm(c70_new_IndustrializedUrban_Info_data)


################################################################################ Now extract the data for input for functions
## filtering out the disease__study combinations 
IndustrializedUrbanCohortInfo= IndustrializedUrbanCohortInfo %>% 
  filter(`80_20`==1) %>% 
  filter(controlcount>=15 & diseaseCount>=15)

# getting the control samples of the shortlisted studies for the discovery cohort.
controlSampleList= list()
uniqueStudies= unique(IndustrializedUrbanCohortInfo$study)

# storing the control samples for all the studies in the list.
for(study in uniqueStudies)
{
  currStudy= IndustrializedUrban_metadata %>% 
    filter(study_name==study) %>% 
    filter(study_condition=="control")
  
  controlSampleList[[study]]= currStudy$sample_id
}

# expanding the metadata based on disease category
expandedMetadata= IndustrializedUrban_metadata %>% 
  separate_rows(diseaseCat, sep=";")

allDiseases= unique(IndustrializedUrbanCohortInfo$disease)  

# it will contain the sample_ids for the disease__study combinations
diseaseSampleList= list()

for(dis in allDiseases)
{
  temp= IndustrializedUrbanCohortInfo %>% 
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

################################################################################ Now Run the functions for Disease Analysis

MainOutputFolder= "C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\disease_analysis\\new_files\\IndustrializedUrban\\disease_analysis_results"

IndustrializedUrban_iterative_disease<- diseaseAnalysisIterationPipeline(10, unique(IndustrializedUrbanCohortInfo$disease), controlSampleList, diseaseSampleList, MainOutputFolder, IndustrializedUrban_spProfile, NewGutAssociatedSpecies)
  
################################################################################# Get the Mean Health score for each of the species

iterationHealthScore_rankedScaledDifference= as.data.frame(IndustrializedUrban_iterative_disease[["iterationOutputDf_ranked_scaledDifference"]])

IndustrializedUrban_HealthScore = as.data.frame(apply(iterationHealthScore_rankedScaledDifference,1,mean))
colnames(IndustrializedUrban_HealthScore) <- "health_score"

################################################################################# Save the workspace and individual Health Score

save(IndustrializedUrban_HealthScore, file = "C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\disease_analysis\\new_files\\IndustrializedUrban\\disease_analysis_results\\IndustrializedUrban_HealthScore.RData")

save.image("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/IndustrializedUrban/c80_new_IndustrializedUrban_Iterative_disease_analysis_workspace.RData")

################################################################################  Do the correlation of health score for IndustrializedUrban with the health score of overall Data (Previously done by Abhishek)
