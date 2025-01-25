
library(dplyr)

load("c93_logisticRegression_FeatureSelection.RData")
load("FunctionalClassificationAnalysis_stage2_part1.RData")
load("FunctionalClassificationAnalysis_stage2_part2.RData")

mannWhitney_batch = function(x,y)
{
  print("in man funciton")
  p_array <- NULL;
  medianDirection <- NULL;
  meanDirection <- NULL;
  
  z <- intersect(rownames(x),rownames(y));
  
  for(i in 1:length(z))
  {
    printed= paste0(i," out of ", length(z))
    
    printed= paste0(i," ",z[i])
   
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
    
    medianDirection[i] <- ifelse(median(as.numeric(x[z[i],]),na.rm=TRUE) > median(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(median(as.numeric(x[z[i],]),na.rm=TRUE) < median(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    
    meanDirection[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
   
    i <- i + 1;
  }
  out <- as.data.frame(cbind(medianDirection, meanDirection, p_array,p.adjust(p_array, method="fdr")));
  colnames(out)[4]= "q_array"
  
  
  rownames(out) <- z;
  # NA values will be either in p-value column or q-value column, if we take them as 1 it will still be insignificant and won't be considered.
  out <- as.data.frame(apply(out,1,function(x)(ifelse(is.nan(x),1,x))));
  # print(head(out))
  return(t(out));
}

overlapping_HACKS <- SpeciesScores[rownames(SpeciesScores) %in% rownames(df_select_features),]
overlapping_HACKS <- overlapping_HACKS[rownames(df_select_features),]

hackGenes= function(functionalProfile)
{
  cat("profile dimensitons before removing empty rows and columns ",dim(functionalProfile),"\n")
  # removing colums with zero sum
  functionalProfile= functionalProfile[,colSums(functionalProfile)>0]
  
  # removing rows with zero sum
  functionalProfile= functionalProfile[rowSums(functionalProfile)>0,]
  
  cat("profile dimensitons after removing empty rows and columns ",dim(functionalProfile),"\n")
  
  functionalProfile$HACKScore <- overlapping_HACKS$HACKScore
    
  # dividing the data frame based on hack and non hack species, the threshold for hackScore is 0.75
  hackDf <- as.data.frame(t(functionalProfile[which(functionalProfile$HACKScore>=0.75),]))
  nonHackDf <- as.data.frame(t(functionalProfile[which(functionalProfile$HACKScore<0.75),]))
  
  
  # removing the hack score form the dataframe
  hackDf= hackDf[-which(rownames(hackDf)=="HACKScore"),]
  nonHackDf= nonHackDf[-which(rownames(nonHackDf)=="HACKScore"),]
  
  cat("dimensions of hackDf are: ",dim(hackDf),"\n")
  cat("dimensions of nonhackDf are: ",dim(nonHackDf),"\n")
  
  wilcox_output= mannWhitney_batch(nonHackDf,hackDf)
  
  
  cnt=1
  # getting the categorical directions as per q_value, p_value and meanDirection
  vector= c()
  for(rows in 1:nrow(wilcox_output))
  {
    if(wilcox_output[rows,"meanDirection"]==0)
    {
      vector= c(vector,0)
    }
    else if(wilcox_output[rows,"q_array"]<= 0.15)
    {
      vector= c(vector, wilcox_output[rows,"meanDirection"]*3)
    }
    else if(wilcox_output[rows,"p_array"]<= 0.05)
    {
      vector= c(vector, wilcox_output[rows,"meanDirection"]*2)
    }
    else
    {
      vector= c(vector, wilcox_output[rows,"meanDirection"]*1)
    }
    
  }
  wilcox_output= as.data.frame(wilcox_output)
  
  wilcox_output$association= vector
  
  # applying 2nd approach
  
  hack_MeanDf= apply(hackDf,1,mean)
  
  nonHack_MeanDf= apply(nonHackDf,1,mean)
  
  meanDf= data.frame(hack_MeanDf,nonHack_MeanDf)
  meanDf$ratio= (meanDf$hack_MeanDf)/(meanDf$nonHack_MeanDf)
  rownames(meanDf)= rownames(hackDf)
  
  meanDf= meanDf %>% 
    arrange(-ratio)
  
  outputList= list()
  
  outputList[["wilcox_output"]]= wilcox_output
  outputList[["meanDf"]]= meanDf
  outputList[["hackSpecies"]]= colnames(hackDf)
  outputList[["nonHackSpecies"]]= colnames(nonHackDf)
  
  return(outputList)
  
}

# running the above function for all the above functions.
output= hackGenes(df_select_features)

# taking only those function Ids that are singificantly more in abundance in hack species
check= output$wilcox_output %>% 
  filter(association %in% c(-3,-2))

table(check$association)


# intersect_hackAssociatedGenes_Df= df_select_features[,intersect(rownames(check), rownames(logisticModel_Df))]
union_hackAssociatedGenes_Df= df_select_features[,union(rownames(check), rownames(filtered_logisticModel_Df))]

# adding function name to the features
featureFunctionType= data.frame(colnames(union_hackAssociatedGenes_Df))
colnames(featureFunctionType)= "Feature"
featureFunctionType$functionalProfile= NA


# getting the funcitnoal features of all the function type.
allFunctionalProfiles= ls()[grepl("df_select",ls())]
allFunctionalProfiles= setdiff(allFunctionalProfiles,"df_select_features")

functionalProfileList= list()

for(functions in allFunctionalProfiles)
{
  functionalProfileList[[functions]]= get(functions)
}

for(functions in names(functionalProfileList))
{
  index= which(featureFunctionType$Feature %in% colnames(functionalProfileList[[functions]])) 
  featureFunctionType[index,"functionalProfile"]= sub("df_select_","",functions)
}
rownames(featureFunctionType)= featureFunctionType$Feature

hackSpecies= output$hackSpecies
nonHackSpecies= output$nonHackSpecies


#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------

# heatmap for the wilcox based selected functions.


library(ggplot2)


SpeciesScores= SpeciesScores[rownames(df_select_features),]

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% hackSpecies) %>% 
  arrange(HACKScore)

hackFunctionalProfile= df_select_features[rownames(TempOrder),rownames(check)]

# creating non hack functional profile

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% nonHackSpecies) %>% 
  arrange(HACKScore)

nonHackFunctionalProfile= df_select_features[rownames(TempOrder),rownames(check)]

# combining both the hack and nonHack profile
FunctionalProfile= rbind(hackFunctionalProfile, nonHackFunctionalProfile)
FunctionalProfile[FunctionalProfile>0]=1


# creating the heatmap for the wiclox based functions.

heatmap_union_core_df= FunctionalProfile
heatmap_union_core_df$species= rownames(FunctionalProfile)

library(tidyverse)
union_df_2=heatmap_union_core_df %>%
  pivot_longer(!species, names_to = "Features", values_to = "presence")

catVector= c()
for(rows in 1:nrow(union_df_2))
{
  catVector= c(catVector,ifelse(unlist(union_df_2[rows,"species"]) %in% hackSpecies,"HACK", "non HACK"))
}

# adding category column in uniion_df_2
union_df_2= cbind(union_df_2,catVector)
union_df_2$function_type= NA


# getting the funcitnoal features of all the functin type.
allFunctionalProfiles= ls()[grepl("df_select",ls())]
allFunctionalProfiles= setdiff(allFunctionalProfiles,"df_select_features")

functionalProfileList= list()

for(functions in allFunctionalProfiles)
{
  functionalProfileList[[functions]]= get(functions)
}

for(functions in names(functionalProfileList))
{
  index= which(union_df_2$Features %in% colnames(functionalProfileList[[functions]])) 
  union_df_2[index,"function_type"]= sub("df_select_","",functions)
}


union_df_2$function_type <- factor(union_df_2$function_type, levels = unique(union_df_2$function_type))

union_df_2$catVector <- factor(union_df_2$catVector, levels = unique(union_df_2$catVector))

union_df_2$species <- factor(union_df_2$species, levels = unique(union_df_2$species))

library(ggplot2)
# color_vector <- c("#eddfb3", "#ca955c","#132b43","#56b1f7")
color_vector <- c("#132b43","#56b1f7")

pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c87_wilcoxFeatures.pdf", width = 20, height = 10);
ggplot(union_df_2, aes(x = Features, y = species, fill = factor(presence)))+
  geom_tile()+
  scale_fill_manual(values = color_vector) +
  theme(axis.text.y = element_text(size = 5,color = "black"))+
  theme(axis.text.x= element_text(size=3, color="black", angle= 90))+
  geom_tile(color = "black")+
  facet_grid( catVector ~ function_type , scales="free",space="free")+
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8))
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
# crating a heatmap by adding one more filter on intersect features of the wilcox and logistic outputs. We will consider only
# those features which are at least 2 folds higher in the hack species.

SpeciesScores= SpeciesScores[rownames(df_select_features),]

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% hackSpecies) %>% 
  arrange(HACKScore)

hackFunctionalProfile= df_select_features[rownames(TempOrder),colnames(union_hackAssociatedGenes_Df)]
hackFunctionalProfile[hackFunctionalProfile>0]=1

# creating non hack functional profile

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% nonHackSpecies) %>% 
  arrange(HACKScore)

nonHackFunctionalProfile= df_select_features[rownames(TempOrder),colnames(union_hackAssociatedGenes_Df)]
nonHackFunctionalProfile[nonHackFunctionalProfile>0]=1

# creating a data frame that contains the information about the prevalence of the features in hack and non hack species
FeatureDetectionRate= data.frame(matrix(nrow=ncol(union_hackAssociatedGenes_Df), ncol= 2))
rownames(FeatureDetectionRate)= colnames(union_hackAssociatedGenes_Df)
colnames(FeatureDetectionRate)= c("Hack_detectionRate", "nonHack_detectionRate")

# making sure that both hack and non hack profile are in same order as that of the FeatureDetectionRate
hackFunctionalProfile= hackFunctionalProfile[,colnames(union_hackAssociatedGenes_Df)]
nonHackFunctionalProfile= nonHackFunctionalProfile[,colnames(union_hackAssociatedGenes_Df)]


FeatureDetectionRate$Hack_detectionRate= round(colSums(hackFunctionalProfile)/nrow(hackFunctionalProfile),2)
FeatureDetectionRate$nonHack_detectionRate= round(colSums(nonHackFunctionalProfile)/nrow(nonHackFunctionalProfile),2)

FeatureDetectionRate$foldChange= FeatureDetectionRate$Hack_detectionRate/FeatureDetectionRate$nonHack_detectionRate


#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------

# getting the frequency of features cross the different fold change threshold. It will help to decide what fold change to take
# to select the final features.

thresholds= seq(1,5,0.1)

freqDb= as.data.frame(matrix(nrow= length(thresholds), ncol= 2))
rownames(freqDb)= thresholds
colnames(freqDb)= c("Thresholds", "Frequency")
freqDb$Thresholds= rownames(freqDb)

for(value in thresholds)
{
  featureCnt= length(which(FeatureDetectionRate$foldChange>=value))
  freqDb[as.character(value),"Frequency"]= featureCnt
}

ggplot(freqDb, aes(x= Thresholds,y = Frequency)) +
  geom_bar(stat= "identity")+
  geom_text(aes(label= Frequency, vjust= -1))

Fold_featureList= rownames(FeatureDetectionRate)[which(FeatureDetectionRate$foldChange>=2)]

# checking the count for each function type for the selected features.
table(featureFunctionType[Fold_featureList,"functionalProfile"])

#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
# library(gplots)
# creating the heatmap based on the selected threshold

SpeciesScores= SpeciesScores[rownames(df_select_features),]

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% hackSpecies) %>% 
  arrange(HACKScore)

hackFunctionalProfile= df_select_features[rownames(TempOrder),Fold_featureList]
hackFunctionalProfile[hackFunctionalProfile>0]=1


# the feature must be identified at least in 25% of hack species
index= which(colSums(hackFunctionalProfile)/nrow(hackFunctionalProfile)>=0.25) # .25 is taken instead of .25 coz it should get detected in at least 6 species
hackFunctionalProfile= hackFunctionalProfile[,index]

#---------------------------
# this is the final hack associated funcitonal profile for which the heatmap will also be made.
final_HackFunctionalProfile= df_select_features[,colnames(hackFunctionalProfile)]
#---------------------------

# creating non hack functional profile

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% nonHackSpecies) %>% 
  arrange(HACKScore)

nonHackFunctionalProfile= df_select_features[rownames(TempOrder),Fold_featureList]
nonHackFunctionalProfile[nonHackFunctionalProfile>0]=1
nonHackFunctionalProfile= nonHackFunctionalProfile[,colnames(hackFunctionalProfile)]


# combining both the hack and nonHack profile
FunctionalProfile= rbind(hackFunctionalProfile, nonHackFunctionalProfile)
FunctionalProfile[FunctionalProfile>0]=1


# creating the heatmap for the selected functional ids.
heatmap_union_core_df= FunctionalProfile
heatmap_union_core_df$species= rownames(FunctionalProfile)

library(tidyverse)
union_df_2=heatmap_union_core_df %>%
  pivot_longer(!species, names_to = "Features", values_to = "presence")

catVector= c()
for(rows in 1:nrow(union_df_2))
{
  catVector= c(catVector,ifelse(unlist(union_df_2[rows,"species"]) %in% hackSpecies,"HACK", "non HACK"))
}

# adding category column in uniion_df_2
union_df_2= cbind(union_df_2,catVector)
union_df_2$function_type= NA


# getting the funcitnoal features of all the functin type.
allFunctionalProfiles= ls()[grepl("df_select",ls())]
allFunctionalProfiles= setdiff(allFunctionalProfiles,"df_select_features")

functionalProfileList= list()

for(functions in allFunctionalProfiles)
{
  functionalProfileList[[functions]]= get(functions)
}

for(functions in names(functionalProfileList))
{
  index= which(union_df_2$Features %in% colnames(functionalProfileList[[functions]])) 
  union_df_2[index,"function_type"]= sub("df_select_","",functions)
}


union_df_2$function_type <- factor(union_df_2$function_type, levels = unique(union_df_2$function_type))

union_df_2$catVector <- factor(union_df_2$catVector, levels = unique(union_df_2$catVector))

union_df_2$species <- factor(union_df_2$species, levels = unique(union_df_2$species))

library(ggplot2)
# color_vector <- c("#eddfb3", "#ca955c","#132b43","#56b1f7")
color_vector <- c("#132b43","#56b1f7")    

pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c87_HackFinal_Features.pdf", width = 20, height = 10);
ggplot(union_df_2, aes(x = Features, y = species, fill = factor(presence)))+
  geom_tile()+
  scale_fill_manual(values = color_vector) +
  theme(axis.text.y = element_text(size = 5,color = "black"))+
  theme(axis.text.x= element_text(size=3, color="black", angle= 90))+
  geom_tile(color = "black")+
  facet_grid( catVector ~ function_type , scales="free",space="free")+
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8))
dev.off()

#------------------------------------------------------------------------------------------------------------------------------------
# getting the plot that shows that the number of features detected are getting reduced with reduce in the HACK Score

copySpeciesScore= SpeciesScores


SpeciesScores= SpeciesScores[rownames(final_HackFunctionalProfile),]


binaryData= as.data.frame(apply(final_HackFunctionalProfile,2,function(x)(ifelse(x>0,1,0))))

# getting the total features detected for each species.
FeatureCountData= as.data.frame(apply(binaryData,1,sum))

colnames(FeatureCountData)= "FeatureCount"                  


FeatureCountData$HACKScore= SpeciesScores[rownames(FeatureCountData),4]

FeatureCountData= FeatureCountData %>% 
  arrange(HACKScore)

FeatureCountData$rank= c(1:122)


cor(FeatureCountData$FeatureCount, FeatureCountData$HACKScore) # 0.54


pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c87_HackFinal_Features_dotPlot.pdf");
ggplot(FeatureCountData, aes(x=rank,
                             y= FeatureCount))+
  geom_point(aes(color= "#6c0000", size= 0.1))+
  theme_bw()+
  xlab("Ranked Species based on HACK Score")+
  ylab("Functional Feature Count")+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size= 25,color="black"))+
  geom_smooth(method= "lm", color= "black")
dev.off()



save(final_HackFunctionalProfile,
     union_hackAssociatedGenes_Df,
     featureFunctionType,
     hackSpecies,
     nonHackSpecies,
     file= "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c87_hackAssociated_functinonal_Ids.RData")
