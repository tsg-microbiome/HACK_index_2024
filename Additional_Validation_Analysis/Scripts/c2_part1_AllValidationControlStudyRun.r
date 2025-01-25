### This script is the 1st part of c56 (3R Analysis). For Validation datasets

# code name- c1_AllValidationControlStudyRun.R


# In this startegy, all 201 species are taken one by one, and  distance matrix of the rest of the species profile is created
# after renormalizing it. Then applying the envfit of the vegan package on the distance matrix and the selected species.


#### Load all the species profile for 14 different studies.
#load("C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/c1_All_Validation_studies_SpProfile.RData")
load("c2_All_Validation_studies_SpProfile.RData")

#created the list of all the study profile objects
object_names <- ls()

AllStudy_list <- list()
for (i in 1:14) {
  AllStudy_list[[object_names[i]]] <- get(object_names[i])
}

#### Load the libraries
library(vegan)
library(dplyr)
library(ade4)

# inpudata
# 1) list of control studies
# 2) species that has to be used for the model

#load("C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/NewGutAssociatedSpecies.RData")
load("NewGutAssociatedSpecies.RData")

species_to_check <- NewGutAssociatedSpecies
rm(NewGutAssociatedSpecies)

# function to find out the prevalent species in a study

# data= input study
# threshold= it will make sure that the given species is detected in at least this many samples
compute_prevalent_single_data <- function(data,threshold){		
  detection_percentage <- colSums(apply(data,2,function(x)(ifelse(x>0,1,0))))/nrow(data)
  highly_detected <- names(which(detection_percentage>=threshold))
  return(highly_detected)
}

# function for the envFit model

# species= the species on which the model has to applied against all the other species.
# inputData= it will take the study as input
# outputFolder= this will take the name of the ouptutfolder where the output of the input speices will get stored
keystoneInfluence <- function(species, inputData, outputFolder) {
  tryCatch({
    
    # removing selected species
    colIndex= which(colnames(inputData)==species)
    newData= inputData[,-colIndex]
    
    # removing empty rows
    newData <- newData[which(rowSums(newData)!= 0), ]
    print(dim(newData))
    
    # normalizing the data 
    newData <- newData / rowSums(newData)
    
    # checking if the data got renormalized
    print(unique(rowSums(newData)))
    
    cat("Creating distance matrix\n")
    distanceMatrix <- vegdist(newData, method = "bray")
    
    print("distance matrix done")
    cat("Creating dudi.pco\n")
    
    pco <- dudi.pco(distanceMatrix, scannf = FALSE)
    
    print("pco is created")
    pcoPointsDf <- pco$li
    
    cat("Generating model\n")
    # since, empty rows are removed it rownames(newData) is specifically mentioned to make sure the sample sequence is same while running the model
    model <- envfit(pcoPointsDf ~ inputData[rownames(newData), species])
    
    check= as.data.frame(t(c(model$vector$r, model$vector$pvals)))
    colnames(check)= c("r2","p-value")
    
    # dir.create(paste0("/home/abhishek/microbiome_project/output_files/c56_AllStudyOutput/",outputFolder,"/"))
    save(check, file= paste0("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\validation_3R_all_studies\\output_3R_Part1\\",outputFolder,"\\",species,".RData"))
    # save(check, file= "D:\\downloads\\check.RData")
    cat(species, "done\n")
    return()
  }, error = function(e) {
    cat("Error in keystoneInfluence for species:", species, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    return("abhs")  # Return an error indicator if needed
  })
}

all_dataset_prevalent_species= list()

# if some species are absent in dataset then adding them with zero abundance for each sample
for(study in names(AllStudy_list))		#iterate over each study
{
  #View(data_list[[index]])
  #print(index)
  #print("##############################")
  #cnt=0
  for(species in species_to_check)		# iterate over each test species
  {
    #species_name= paste0("",species)
    if (species %in% colnames(AllStudy_list[[study]]))	# if test species present as column name in study then it will move to next test species
    {
      next
      #cnt=cnt+1
    }
    else											# if it is not present then it will create a new column with test species name and all values will be zero
    {	
      AllStudy_list[[study]][[species]]= with(AllStudy_list[[study]],0)		
      #print("false")
    }
  }
  #print(cnt)
}

studyCount=1
for(study in names(AllStudy_list))
{
  printed= paste0(studyCount, "-----------------------------------------------------------------------------------------", study)
  studyCount= studyCount+1
  print(printed)
  AllStudy_list[[study]]= AllStudy_list[[study]] %>% 
    replace(is.na(.),0)	
  
  print("computing the prevalent species from studies")
  all_dataset_prevalent_species[[study]]= compute_prevalent_single_data(AllStudy_list[[study]], 0.70)
  
  # normalizing the dataset
  print("normalizing the data")
  train_data_norm= AllStudy_list[[study]]/rowSums(AllStudy_list[[study]])		
  train_data_norm= train_data_norm %>% 
    replace(is.na(.),0)
  
  train_data_norm= train_data_norm[rowSums(train_data_norm)>0,]
  
  inputData= train_data_norm[,species_to_check]
  
  # printed= paste0("----------------------------------------------",study,"-------------------------------------------")
  # print(printed)
  count=1
  dir.create(paste0("C:\\Users\\ompra\\Downloads\\Cell_reports_revision\\validation_3R_all_studies\\output_3R_Part1\\",study,"\\"))
  
  for(species in species_to_check)
  {
    printed= paste0("#####################",count," ",species,"###########################")
    count= count+1
    print(printed)
    keystoneInfluence(species, inputData, study)
  }
}


save(all_dataset_prevalent_species,file="C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/c2_part1_all_validation_dataset_prevalent_species.RData")
save.image("C:/Users/ompra/Downloads/Cell_reports_revision/validation_3R_all_studies/c2_part1_workspace_all_validation_dataset.RData")

