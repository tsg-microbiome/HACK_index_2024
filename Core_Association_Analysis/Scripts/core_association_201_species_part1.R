
# output: 
# 1) prDf_new
# 2) r2Df_new_ranked
# 3) prevalentDf_new
# 4) r2Df_new (not needed but save for future use)


library(vegan)
library(dplyr)
library(ade4)
library(readxl)

# Importing metadata and species profile
#load("C:/Users/ompra/Downloads/Cell_reports_revision/core_influencers2/CoreInf_SpProfile_Metadata.RData")
load("CoreInf_SpProfile_Metadata.RData")

# Importing the 201 filtered species names
#load("C:/Users/ompra/Downloads/Cell_reports_revision/core_influencers2/NewGutAssociatedSpecies.RData")
load("NewGutAssociatedSpecies.RData")


# Remove the columns with colSums == 0
CoreInf_SpProfile <- CoreInf_SpProfile[,which(colSums(CoreInf_SpProfile)!=0),]

# Add study_name column to species profile to run the study specififc analysis using below functions
CoreInf_SpProfile$sample_id <- rownames(CoreInf_SpProfile)
CoreInf_SpProfile <- left_join(CoreInf_SpProfile,CoreInf_Metadata[,c("sample_id","study_name")],by = "sample_id")
CoreInf_SpProfile <- data.frame(CoreInf_SpProfile)
rownames(CoreInf_SpProfile) <- CoreInf_SpProfile$sample_id
CoreInf_SpProfile$sample_id <- NULL


# some of the sample names are not maching due to extension of study name in sample name, so add manually the stduy name
CoreInf_SpProfile$study_name <- ifelse(is.na(CoreInf_SpProfile$study_name),"RampelliS_2015", CoreInf_SpProfile$study_name)


###############

## Create the function for finding species in one study that have prevalance more or equal to 70%

compute_prevalent_single_data <- function(data,threshold){		
  detection_percentage <- colSums(apply(data,2,function(x)(ifelse(x>0,1,0))))/nrow(data)
  highly_detected <- names(which(detection_percentage>=threshold))
  return(highly_detected)
}

##############

# create the function for finding core influencers
keystoneInfluence <- function(All_species_Profile, species_for_core_finding) {
  
  ## Create an empty list to store the dfs, one df for one study which will be generated after the analysis
  output_all_stdy_dfs <- list()
  
  ## Create an empty list to store the species which have abundance above 70%
  all_dataset_prevalent_species <- list()
  
  ## Extract the unique study_names from the species profile
  unique_study_names <- unique(All_species_Profile$study_name)
  
  countn <- 1
  ## Add the for loop on study_names
  for (stdy in unique_study_names) {
    tryCatch({
      print(paste0("############## ",countn," ",stdy, " #################"))
      countn <- countn + 1
      
      # separate species_profile for one study
      one_species_prof <- All_species_Profile[which(All_species_Profile$study_name == stdy),]
      
      # remove study_name column from the species_profile
      one_species_prof <- one_species_prof[, !(colnames(one_species_prof) %in% c("study_name"))]
      
      # Replace NA value with 0 if any present in species_profile of current study
      one_species_prof[is.na(one_species_prof)] <- 0
      
      ## Storing the species which are detected in at least 70% samples
      all_dataset_prevalent_species[[stdy]] <- compute_prevalent_single_data(one_species_prof, 0.70)
      
      # Normalize the species_profile of current study
      one_species_prof_norm <- one_species_prof / rowSums(one_species_prof)
      
      # Replace NA value with 0 if any present in species_profile of current study
      one_species_prof_norm[is.na(one_species_prof_norm)] <- 0
      
      # remove rows whose row sums is 0
      one_species_prof_norm <- one_species_prof_norm[rowSums(one_species_prof_norm) > 0, ]
      
      # Filter only species that are present in species_for_core_finding (New species associated with gut)
      inputData <- one_species_prof_norm[, species_for_core_finding]
      
      ## Run the ENV-Fit on each species 
      # count <- 1
      
      # Create an empty dataframe to store the p-value and r value for each species for current study
      summary_envfit_onestudy <- as.data.frame(matrix(NA, nrow = length(species_for_core_finding), ncol = 2))
      colnames(summary_envfit_onestudy) <- c("r2", "p-value")
      rownames(summary_envfit_onestudy) <- species_for_core_finding
      
      ## Loop over the species in species_for_core_finding to find the ENV-Fit results.
      for (species in species_for_core_finding) {
        tryCatch({
          # printed <- paste0("#####################", count, " ", stdy, " ", species, "###########################")
          # count <- count + 1
          # print(printed)
          
          # removing selected species (1R)
          colIndex <- which(colnames(inputData) == species)
          newData <- inputData[, -colIndex]
          
          # removing empty rows
          newData <- newData[which(rowSums(newData) != 0), ]
          #print(dim(newData))
          
          # normalizing the data (2R)
          newData <- newData / rowSums(newData)
          
          # checking if the data got renormalized
          cat(length(which(rowSums(newData) == 0)))
          
          #cat("Creating distance matrix\n")
          distanceMatrix <- vegdist(newData, method = "bray")
          
          #print("distance matrix done")
          #cat("Creating dudi.pco\n")
          
          pco <- dudi.pco(distanceMatrix, scannf = FALSE)
          
          #print("pco is created")
          pcoPointsDf <- pco$li
          
          #cat("Generating model\n") (3R)
          # since, empty rows are removed the rownames(newData) are specifically mentioned to make sure the sample sequence is same while running the model
          model <- envfit(pcoPointsDf ~ inputData[rownames(newData), species])
          
          check <- as.data.frame(t(c(model$vector$r, model$vector$pvals)))
          colnames(check) <- c("r2", "p-value")
          
          summary_envfit_onestudy[species, "r2"] <- check[1, "r2"]
          summary_envfit_onestudy[species, "p-value"] <- check[1, "p-value"]
          
          #cat(species, "done\n")
        }, error = function(e) {
          cat("Error in processing species:", species, "in study:", stdy, "\n")
          print(e)
        })
      }
      cat(" done", "\n")
      output_all_stdy_dfs[[stdy]] <- summary_envfit_onestudy
      gc()
    }, error = function(e) {
      cat("Error in processing study:", stdy, "\n")
      print(e)
    })
  }
  final_output_list <- list(ENV_fit_summary = output_all_stdy_dfs, all_dataset_prevalent_species = all_dataset_prevalent_species)
  return(final_output_list)
}



#############
# see if the species in NewGutAssociations is same as that in core influencers.
selected_species <-  intersect(colnames(CoreInf_SpProfile),NewGutAssociatedSpecies)

# Run the above code for all the studies previously used in finding Core Influencers
new_keystoneInfl <- keystoneInfluence(CoreInf_SpProfile,selected_species)


###########
all_study_wise_keystone <- new_keystoneInfl[["ENV_fit_summary"]]


#Build the function to extract pvalues and r2values of all studies in two different dataframes. (generalized function where column name has to be given which we want to extract from study specific dataframes)

extract_same_column_from_multiple_dfs_in_list <- function(df_list,x) {
  combined_df <- data.frame(matrix(NA, nrow = length(selected_species),ncol = 1))
  colnames(combined_df) <- "species"
  combined_df$species <- selected_species
  
  for(nam in names(df_list)){
    one_df <- df_list[[nam]]
    one_df$species <- rownames(one_df)
    
    combined_df <- left_join(combined_df,one_df[,c(x,"species"), drop =FALSE], by = "species")
    combined_df <- data.frame(combined_df)
    
  }
  rownames(combined_df) <- combined_df$species
  combined_df$species <- NULL
  colnames(combined_df) <- names(df_list)
  return(combined_df)
}

# r2Df_new 
r2Df_new <- extract_same_column_from_multiple_dfs_in_list(all_study_wise_keystone,"r2")

prDf_new <- extract_same_column_from_multiple_dfs_in_list(all_study_wise_keystone,"p-value")



############
# Makr the r2Df_new as ranked df
rank_scale=function(x)
{
  # x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

r2Df_new_ranked = as.data.frame(apply(r2Df_new,2,rank_scale))


#############



#############
## Function for prevalance
prevalence_function <- function(All_species_profile,species_to_check,studies){
  
  prevalentDf= as.data.frame(matrix(NA, nrow= length(species_to_check), ncol= length(studies)))
  rownames(prevalentDf)= species_to_check
  colnames(prevalentDf)= studies
  
  stCount=1
  for(stdy in colnames(prevalentDf))
  {
    print(paste0("##### ",stCount," ",stdy," #####"))
    stCount= stCount+1
    
    # separate species_profile for one study
    one_species_prof <- All_species_profile[which(All_species_profile$study_name == stdy),]
    
    # remove study_name column from the species_profile
    one_species_prof <- one_species_prof[, !(colnames(one_species_prof) %in% c("study_name"))]
    
    # Replace NA value with 0 if any present in species_profile of current study
    one_species_prof[is.na(one_species_prof)] <- 0
    
    # Remove the rows that have zero rowsums
    one_species_prof <- one_species_prof[which(rowSums(one_species_prof)!=0),]
    
    # Remove the columns that have zero colsums
    one_species_prof <- one_species_prof[,which(colSums(one_species_prof)!=0)]
    
    currStudy= as.data.frame(t(apply(one_species_prof,1,function(x)(ifelse(x>0,1,0)))))
    
    for(species in rownames(prevalentDf))
    {
      if(species %in% colnames(currStudy))
      {
        prevalentDf[species,stdy]= sum(currStudy[,species])/nrow(currStudy)
      }
      else
      {
        prevalentDf[species,stdy]=0
      }
      
    }
  }
  return(prevalentDf)
}

prevalentDf_new <- prevalence_function(CoreInf_SpProfile,selected_species,unique(CoreInf_SpProfile$study_name))


################
# Filtering species following detection in r2Df_new_rankscaled and prevalentDf_new
CoreKeyStoneDf <- data.frame(apply(r2Df_new_ranked,2,function(x)(ifelse(x>=0.70,1,0))) * apply(prevalentDf_new,2,function(x)(ifelse(x>=0.65,1,0))))

# add the detection column in CoreKeyStone_df
CoreKeyStoneDf$detection <- rowSums(CoreKeyStoneDf)
CoreKeyStoneDf <- CoreKeyStoneDf[order(CoreKeyStoneDf$detection, decreasing = TRUE), ]
CoreKeyStoneDf$species <- rownames(CoreKeyStoneDf)

# Adding Influence_score using Rank_Scale function
CoreKeyStoneDf$Influence <- (CoreKeyStoneDf$detection)/72
CoreKeyStoneDf$Influence <- rank_scale(CoreKeyStoneDf$Influence)

# Note: Don't use study wise prevalant specie sgenerated uisng core key stone function from new_keystoneInfl dataframe (Bcos it is based on 70% prevalance but we have used 65% now)

# Extract the species and their Influence score in new df
Influence_Score_df <- CoreKeyStoneDf[,74:75]


# Save the Important dfs.
#save(r2Df_new_ranked,r2Df_new, prDf_new,
#     prevalentDf_new,Influence_Score_df,
#     file= "C:/Users/ompra/Downloads/Cell_reports_revision/core_influencers2/c56_core_key_stone_new.RData")


#save(Influence_Score_df, file = "C:/Users/ompra/Downloads/Cell_reports_revision/core_influencers2/CoreInfluencerScore_new.RData")

#save.image("C:/Users/ompra/Downloads/Cell_reports_revision/core_influencers2/core_association_201_species.RData")
save.image("core_association_201_species.RData")



