
setwd("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output")
library(xlsx)
library(dplyr)

###########################################################################

load("G:/My Drive/HACK paper/CellReports_revision/Functional_Analysis_Overall_Output/union_features_func_prof_all_genomes.RData")
load("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c87_hackAssociated_functinonal_Ids.RData")

##### CONSERVATION PROFILE GENERATION HACKS ######

df <- union_features_func_prof_all_genomes

features <- names(df)[2:ncol(df)]
hacks <- hackSpecies

df[,2:ncol(df)] <- ifelse(df[,2:ncol(df)] > 0, 1,0)

conservation_df <- as.data.frame(matrix(0, nrow = length(features), ncol = length(hacks)))
rownames(conservation_df) <- features
colnames(conservation_df) <- hacks
  
for (i in 1:length(hacks)) {
  species_df <- df[df$Species == hacks[i],]
  species_df <- species_df[,-1]
  for (j in 1:length(features)){
    present_in_genomes <- sum(species_df[,j] > 0)
    total_genomes <- nrow(species_df)
      
    conservation_df[j, i] <- as.numeric((present_in_genomes / total_genomes) * 100)
  }
  print(paste0("done for:",hacks[i]))
}

write.xlsx(conservation_df,"G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\Conservation_df_HACKS.xlsx",row.names = T)

########### PERCENTAGE DETECTION ############

union_features_HACKS <- union_hackAssociatedGenes_Df[rownames(union_hackAssociatedGenes_Df) %in% hackSpecies,]
union_features_nonHACKS <- union_hackAssociatedGenes_Df[rownames(union_hackAssociatedGenes_Df) %in% nonHackSpecies,]

calculate_detection_percentage <- function(df) {
  
  percentage_greater_than_zero <- colSums(df > 0) / nrow(df) * 100
  
  percentage_detection_df <- data.frame(
    Column = colnames(df),
    Percentage = percentage_greater_than_zero
  )
  
  return(percentage_detection_df)
}

percentage_detection_HACKS_df <- calculate_detection_percentage(union_features_HACKS)
percentage_detection_nonHACKS_df <- calculate_detection_percentage(union_features_nonHACKS)

percentage_detection_overall_df <- data.frame(features = percentage_detection_HACKS_df$Column, 
                                              percentage_detection_HACKS = percentage_detection_HACKS_df$Percentage,
                                              percentage_detection_nonHACKS = percentage_detection_nonHACKS_df$Percentage)

####################################################################

suppl_table22 <- as.data.frame(cbind(percentage_detection_overall_df,conservation_df))
suppl_table22$detection_ratio <- suppl_table22$percentage_detection_HACKS/suppl_table22$percentage_detection_nonHACKS
suppl_table22$conservation <- rowSums(suppl_table22[, 4:19] > 90)
suppl_table22$score <- ifelse(suppl_table22$detection_ratio >= 2 & suppl_table22$conservation >= 2, 1 ,0)
suppl_table22 <- suppl_table22[order(suppl_table22$score, decreasing = T),]
score1_df <- suppl_table22[1:150,]
score1_df_sorted <- score1_df[order(score1_df$conservation, decreasing = T),]

SupplementaryTable22_final <- as.data.frame(rbind(score1_df_sorted,suppl_table22[151:1019,]))

write.xlsx(SupplementaryTable22_final,"SupplementaryTable22_final.xlsx",row.names = F)

