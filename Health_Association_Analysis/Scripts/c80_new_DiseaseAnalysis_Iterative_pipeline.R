
############################################################################## Load the data
## Load the disease_analysis_metadata and disease_analysis_spProfile 

c1_new_cohort_wise_master_data <- new.env()
load("disease_analysis_input1.RData",envir = c1_new_cohort_wise_master_data)
load("disease_analysis_input2.RData",envir = c1_new_cohort_wise_master_data)
load("disease_analysis_input3.RData",envir = c1_new_cohort_wise_master_data)
attach(c1_new_cohort_wise_master_data)
disease_analysis_metadata <- disease_analysis_metadata
disease_analysis_spProfile <- disease_analysis_spProfile
detach(c1_new_cohort_wise_master_data)
rm(c1_new_cohort_wise_master_data)

## Load the functions that are used in disease analysis
source("c73_new_DiseaseAnalysis_cohort_functions.R")
## Load the species to check in the data 
load("NewGutAssociatedSpecies.RData")

## Load the summarised information
c70_new_DiseaseAnalysis_CohortInfo_data <- new.env()
load("C:/Users/ompra/Downloads/Cell_reports_revision/disease_analysis/new_files/DiseaseAnalysis_combinedCohort/c70_new_DiseaseAnalysis_CohortInfo_data.RData",envir = c70_new_DiseaseAnalysis_CohortInfo_data)
attach(c70_new_DiseaseAnalysis_CohortInfo_data)
DiseaseAnalysis_cohortInfo <- DiseaseAnalysis_cohortInfo
detach(c70_new_DiseaseAnalysis_CohortInfo_data)
rm(c70_new_DiseaseAnalysis_CohortInfo_data)


################################################################################ Now extract the data for input for functions
## filtering out the disease_study combinations 
DiseaseAnalysis_cohortInfo= DiseaseAnalysis_cohortInfo %>% 
  filter(`80_20`==1) %>% 
  filter(controlcount>=15 & diseaseCount>=15)

# getting the control samples of the shortlisted studies for the discovery cohort.
controlSampleList= list()
uniqueStudies= unique(DiseaseAnalysis_cohortInfo$study)

# storing the control samples for all the studies in the list.
for(study in uniqueStudies)
{
  currStudy= disease_analysis_metadata %>% 
    filter(study_name==study) %>% 
    filter(study_condition=="control")
  
  controlSampleList[[study]]= currStudy$sample_id
}

# expanding the metadata based on disease category
expandedMetadata= disease_analysis_metadata %>% 
  separate_rows(diseaseCat, sep=";")

allDiseases= unique(DiseaseAnalysis_cohortInfo$disease)  

# it will contain the sample_ids for the disease__study combinations
diseaseSampleList= list()

for(dis in allDiseases)
{
  temp= DiseaseAnalysis_cohortInfo %>% 
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

#### Now Run the functions for Disease Analysis

MainOutputFolder= "G:\\My Drive\\HACK paper\\CellReports_revision\\Final_Submission\\Final_v3\\Scripts_GitHub\\Section3_DiseaseAssociationInvestigation\\DATA"

#Extract the intersecting species to check
filtered_NewGutAssociatedSpecies <- sort(intersect(colnames(disease_analysis_spProfile),NewGutAssociatedSpecies))

DiseaseAnalysis_Cohort_iterative_disease <- diseaseAnalysisIterationPipeline(10, unique(DiseaseAnalysis_cohortInfo$disease), controlSampleList, diseaseSampleList, MainOutputFolder, disease_analysis_spProfile, filtered_NewGutAssociatedSpecies)

################################################################################# Get the Mean Health score for each of the species

iterationHealthScore_rankedScaledDifference= as.data.frame(DiseaseAnalysis_Cohort_iterative_disease[["iterationOutputDf_ranked_scaledDifference"]])

Overall_DiseaseAnalysis_HealthScore = as.data.frame(apply(iterationHealthScore_rankedScaledDifference,1,mean))
colnames(Overall_DiseaseAnalysis_HealthScore) <- "health_score"

################################################################################# Save the workspace and individual Health Score

save(Overall_DiseaseAnalysis_HealthScore, file = "G:\\My Drive\\HACK paper\\CellReports_revision\\Final_Submission\\Final_v3\\Scripts_GitHub\\Section3_DiseaseAssociationInvestigation\\DATA\\Overall_DiseaseAnalysis_HealthScore.RData")

save.image("G:\\My Drive\\HACK paper\\CellReports_revision\\Final_Submission\\Final_v3\\Scripts_GitHub\\Section3_DiseaseAssociationInvestigation\\DATA\\c80_new_DiseaseAnalysis_Iterative_pipeline_workspace.RData")

