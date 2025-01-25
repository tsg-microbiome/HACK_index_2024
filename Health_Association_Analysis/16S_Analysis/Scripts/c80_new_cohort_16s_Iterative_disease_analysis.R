######### This code is a replication of c80_Aggregated_iterativeCohortAnalysis.R with some modifications for cohort_16s
## Omprakash 
## Dt.25-09-2024


## This code will use c73 for cohort_16s cohort and run disease analysis for this cohort

############################################################################## Load the data
## Load the metadata_16s and spProfile_16s 
c1_new_cohort_wise_master_data <- new.env()
load("G:\\.shortcut-targets-by-id\\1-VXvgWg44QVMpipxNRFjdw2Oe4gs2t6O\\CoreFinder\\Manuscript\\CellReportsRevision\\Section3_DiseaseAssociationInvestigationFigures\\Amplicon_sequenced\\c1_new_cohort_wise_master_data.RData",envir = c1_new_cohort_wise_master_data)
attach(c1_new_cohort_wise_master_data)
metadata_16s <- metadata_16s
spProfile_16s <- spProfile_16s
detach(c1_new_cohort_wise_master_data)
rm(c1_new_cohort_wise_master_data)

## Load the functions that are used in disease analysis (c73_new_cohort_16s_disease_analysis_functions.R)
source("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/Amplicon_sequenced/c73_new_cohort_16s_disease_analysis_functions.R")
## Load the species to check in the data 
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/NewGutAssociatedSpecies.RData")

## Load the summarised information for cohort_16s cohort from c70 (c70_new_cohort_16s_Info_data.RData)
c70_new_cohort_16s_Info_data <- new.env()
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/Amplicon_sequenced/c70_new_cohort_16s_Info_data.RData",envir = c70_new_cohort_16s_Info_data)
attach(c70_new_cohort_16s_Info_data)
cohort_16sCohortInfo <- cohort_16sCohortInfo
detach(c70_new_cohort_16s_Info_data)
rm(c70_new_cohort_16s_Info_data)


################################################################################ Now extract the data for input for functions
## filtering out the disease__study combinations 
cohort_16sCohortInfo= cohort_16sCohortInfo %>% 
  filter(`80_20`==1) %>% 
  filter(controlcount>=15 & diseaseCount>=15)

# getting the control samples of the shortlisted studies for the discovery cohort.
controlSampleList= list()
uniqueStudies= unique(cohort_16sCohortInfo$study)

# storing the control samples for all the studies in the list.
for(study in uniqueStudies)
{
  currStudy= metadata_16s %>% 
    filter(study_name==study) %>% 
    filter(study_condition=="control")
  
  controlSampleList[[study]]= currStudy$sample_id
}

# expanding the metadata based on disease category
expandedMetadata= metadata_16s %>% 
  separate_rows(diseaseCat, sep=";")

allDiseases= unique(cohort_16sCohortInfo$disease)  

# it will contain the sample_ids for the disease__study combinations
diseaseSampleList= list()

for(dis in allDiseases)
{
  temp= cohort_16sCohortInfo %>% 
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

MainOutputFolder= "C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\disease_analysis\\new_files\\Amplicon_sequenced\\disease_analysis_results"

#Extrcat the intersect species to check
filtered_NewGutAssociatedSpecies <- sort(intersect(colnames(spProfile_16s),NewGutAssociatedSpecies))

cohort_16s_iterative_disease<- diseaseAnalysisIterationPipeline(10, unique(cohort_16sCohortInfo$disease), controlSampleList, diseaseSampleList, MainOutputFolder, spProfile_16s, filtered_NewGutAssociatedSpecies)

################################################################################# Get the Mean Health score for each of the species

iterationHealthScore_rankedScaledDifference= as.data.frame(cohort_16s_iterative_disease[["iterationOutputDf_ranked_scaledDifference"]])

cohort_16s_HealthScore = as.data.frame(apply(iterationHealthScore_rankedScaledDifference,1,mean))
colnames(cohort_16s_HealthScore) <- "health_score"

################################################################################# Save the workspace and individual Health Score

save(cohort_16s_HealthScore, file = "C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\disease_analysis\\new_files\\Amplicon_sequenced\\disease_analysis_results\\cohort_16s_HealthScore.RData")

save.image("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/Amplicon_sequenced/c80_new_cohort_16s_Iterative_disease_analysis_workspace.RData")

################################################################################  Do the correlation of health score for cohort_16s with the health score of overall Data (Previously done by Abhishek)
