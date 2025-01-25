### This script is the 2nd part of c56 (3R Analysis). Here I have did this for Validation datasets

# code name- c2_part2_AllValidationControlStudyRun.R

# the code is written to put all the r2 value into 1 df 

# to get the core species of each study based on intersection of prevalent species and species with r2 value
# within top75 percentile

# finally creating a core microbiome df out of them
library(dplyr)
library(writexl)
library(ggplot2)
library(ggrepel)

#load ValidationDatasetsFor3R
#load("C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/c1_All_Validation_studies_SpProfile.RData")
load("c1_All_Validation_studies_SpProfile.RData")

# Store all the studies data in one list
object_names <- ls()
AllStudy_list <- list()
for (i in 1:14) {
  AllStudy_list[[object_names[i]]] <- get(object_names[i])
}

#load output from 1st part of c2
#load("C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/c2_part1_all_validation_dataset_prevalent_species.RData")
load("c2_part1_all_validation_dataset_prevalent_species.RData")


folderPath= "C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\validation_3R_all_studies\\output_3R_Part1\\"
allStudy= list.files(folderPath)

# total rows for the r2 df
rowSize= length(list.files(paste0(folderPath,allStudy[1],"\\")))

r2Df= as.data.frame(matrix(nrow= rowSize, ncol= length(allStudy)))
rownames(r2Df)= list.files(paste0(folderPath,allStudy[1],"\\"))
rownames(r2Df)= sub(".RData","",rownames(r2Df))

colnames(r2Df)= allStudy

prDf= as.data.frame(matrix(nrow= rowSize, ncol= length(allStudy)))
rownames(prDf)= list.files(paste0(folderPath,allStudy[1],"\\"))
rownames(prDf)= sub(".RData","",rownames(prDf))

colnames(prDf)= allStudy

rank_scale=function(x)
{
  # x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

# filling the values in r2Df and prDf
sC=1
for(study in allStudy)
{
  printed= paste0(sC," ",study)
  print(printed)
  sC= sC+1
  studyPath= paste0(folderPath,study,"//")
  allSpecies= list.files(studyPath)
  
  spC=1
  for(species in allSpecies)
  {
    printed= paste0(spC," ",species)
    spC= spC+1
    # print(printed)
    
    load(paste0(studyPath,species))
    
    # if(check$`p-value`<=0.05)
    # {
    #   r2Df[sub(".RData","",species),study]= check$r2 
    # }
    r2Df[sub(".RData","",species),study]= check$r2 
    prDf[sub(".RData","",species),study]= check$`p-value` 
    
  }
}

r2DfNonRankScaled= r2Df
r2Df= as.data.frame(apply(r2Df,2,rank_scale))


prevalentDf= as.data.frame(matrix(nrow= nrow(r2Df), ncol= ncol(r2Df)))
rownames(prevalentDf)= rownames(r2Df)
colnames(prevalentDf)= colnames(r2Df)

stCount=1
for(study in colnames(prevalentDf))
{
  print(stCount)
  stCount= stCount+1
  
  currStudy=AllStudy_list[[study]] %>% 
    replace(is.na(.),0)
  
  currStudy= as.data.frame(t(apply(currStudy,1,function(x)(ifelse(x>0,1,0)))))
  
  for(species in rownames(prevalentDf))
  {
    if(species %in% colnames(currStudy))
    {
      prevalentDf[species,study]= sum(currStudy[,species])/nrow(currStudy)
    }
    else
    {
      prevalentDf[species,study]=0
    }
    
  }
}


###################################################################################################
##################################################################################################
## Create new df for CoreInfluencer and Non-coreInfluencer species in each study.
## Separate the study specific 3R data

#### Create one more Df "csDf" CoreSpecies Df for all the 10 studies
csDf <- data.frame(matrix(ncol = length(r2Df), nrow = nrow(r2Df)))
colnames(csDf) <- colnames(r2Df)
rownames(csDf) <- rownames(r2Df)

# Loop through each row and column to set values in csDf based on conditions
for (i in 1:nrow(r2Df)) {
  for (col in colnames(r2Df)) {
    if (r2Df[i, col] >= 0.7 & prevalentDf[i, col] >= 0.65) {
      csDf[i, col] <- "CoreInfluencer"
    } else {
      csDf[i, col] <- "NonCoreInfluencer"
    }
  }
}


#### Assigning new colnames to the prDf, prevalentDf, r2Df, r2DfNonRankScaled, csDf
col_names <- colnames(csDf)
# Replace the pattern
new_col_names <- sub("_select_data_norm2$", "_3R", col_names)
# Rename the columns in the data frame
colnames(csDf) <- new_col_names
colnames(r2Df) <- new_col_names
colnames(prDf) <- new_col_names
colnames(prevalentDf) <- new_col_names
colnames(r2DfNonRankScaled) <- new_col_names

########### Grouping these dfs according to study
all_dfs <- list(prDf, prevalentDf, r2Df, r2DfNonRankScaled,csDf)

# Create empty list to store data frames
StudyWise_dfs_list <- list()

# Loop through each column name and create data frames dynamically
for (col in colnames(r2Df)) {
  # Extract columns from each data frame and combine into a single data frame
  temp_df <- do.call(cbind, lapply(all_dfs, function(df) df[col]))
  
  # Assign column names
  colnames(temp_df) <- c("pval", "prevalent", "r2val", "r2val_NonRankScaled","CoreSp")
  
  # Assign the data frame to the list
  StudyWise_dfs_list[[col]] <- temp_df
}

##################################################################################################
#################################################################################################

# save only few dfs and lists to upload on drive (sir)

#save(csDf,r2Df,r2DfNonRankScaled, prDf,prevalentDf,StudyWise_dfs_list,AllStudy_list, file= "C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/c2_part2_3R_Analysis_Validation_14Stdy.RData")
save(csDf,r2Df,r2DfNonRankScaled, prDf,prevalentDf,StudyWise_dfs_list,AllStudy_list, file= "c2_part2_3R_Analysis_Validation_14Stdy.RData")

###### Code to separate dfs from StudyWise_dfs_list to global environment
for (col_name in names(StudyWise_dfs_list)) {
  assign(col_name, StudyWise_dfs_list[[col_name]], envir = .GlobalEnv)
}

all_3R_dfs <- ls()[grep("3R",ls())]
save(list = all_3R_dfs, file = "Validation_3R_Analysis.RData")



# save the whole environment
save.image(file= "c2_part2_workspace_all_validation_3R_analysis.RData")
