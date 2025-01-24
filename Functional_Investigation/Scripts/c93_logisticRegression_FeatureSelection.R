
# applying logistic regression on the entire selected features to filter them out further.


library(dplyr)

load("G:/My Drive/HACK paper/CellReports_revision/Final_Submission/Final_v3/Scripts_GitHub/Section10_FunctionalInvestigation/DATA/FunctionalClassificationAnalysis_stage2_part1.RData")
load("G:/My Drive/HACK paper/CellReports_revision/Final_Submission/Final_v3/Scripts_GitHub/Section10_FunctionalInvestigation/DATA/FunctionalClassificationAnalysis_stage2_part2.RData")
load("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\SpeciesScores_FunctionalAnalysis.RData")

length(intersect(rownames(SpeciesScores), rownames(df_select_features)))

#save(SpeciesScores,file='G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\SpeciesScores_FunctionalAnalysis.RData')

copySpeciesScore= SpeciesScores

SpeciesScores= SpeciesScores[intersect(rownames(SpeciesScores), rownames(df_select_features)),]
df_select_features= df_select_features[rownames(SpeciesScores),]


HACKScore= SpeciesScores$HACKScore

# converting data into binary data based on function detection in species as it is required for logistic regression.
# df_select_features[df_select_features>0]=1


rank_scale=function(x)
{
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

df_select_features= as.data.frame(apply(df_select_features,2,rank_scale))

logisticModel_Df= as.data.frame(matrix(nrow= ncol(df_select_features), ncol= 2))
colnames(logisticModel_Df)= c("Estimate","p_value")
rownames(logisticModel_Df)= colnames(df_select_features)


featureCnt=1

for(features in rownames(logisticModel_Df))
{
  print(featureCnt)
  featureCnt= featureCnt+1
  
  model <- glm(df_select_features[,features] ~ HACKScore,
               family = gaussian) #binomial(link = "logit")) 
  
  model_sum= summary(model)
  
  p_value= model_sum$coefficients[2,4]
  estimate= model_sum$coefficients[2,1]
  
  logisticModel_Df[features,"p_value"]= p_value
  logisticModel_Df[features,"Estimate"]= estimate
}

logisticModel_Df$q_value= p.adjust(logisticModel_Df$p_value, method="fdr")

filtered_logisticModel_Df= logisticModel_Df %>% 
  filter(q_value <=0.05) %>% 
  filter(Estimate>0)


logisticModel_Df$association= ifelse(logisticModel_Df$Estimate>0,"Positive","Negative")

table(logisticModel_Df$association)

save(filtered_logisticModel_Df, file= "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c93_logisticRegression_FeatureSelection.RData")


#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------


# creating the functional profile based on the logistic regression features.
SpeciesScores= SpeciesScores[rownames(df_select_features),]

# creating hack functional profile
hackSpecies= rownames(SpeciesScores)[SpeciesScores$HACKScore>=0.75]

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% hackSpecies) %>% 
  arrange(HACKScore)

hackFunctionalProfile= df_select_features[rownames(TempOrder),rownames(filtered_logisticModel_Df)]

# creating non hack functional profile
nonHackSpecies= rownames(SpeciesScores)[SpeciesScores$HACKScore<0.75]

TempOrder= SpeciesScores %>% 
  filter(rownames(SpeciesScores) %in% nonHackSpecies) %>% 
  arrange(HACKScore)

nonHackFunctionalProfile= df_select_features[rownames(TempOrder),rownames(filtered_logisticModel_Df)]

# combining both the hack and nonHack profile
FunctionalProfile= rbind(hackFunctionalProfile, nonHackFunctionalProfile)

copyProfile= FunctionalProfile
FunctionalProfile[FunctionalProfile>0]=1


# creating the heatmap of logistic regression based features.
library(ggplot2)

heatmap_union_core_df= FunctionalProfile
heatmap_union_core_df$species= rownames(FunctionalProfile)

library(tidyverse)
union_df_2=heatmap_union_core_df %>%
  pivot_longer(!species, names_to = "Features", values_to = "presence")


hackSpecies= rownames(SpeciesScores)[which(SpeciesScores$HACKScore>=0.75)]


catVector= c()
for(rows in 1:nrow(union_df_2))
{
  catVector= c(catVector,ifelse(unlist(union_df_2[rows,"species"]) %in% hackSpecies,"HACK", "Non HACK"))
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

pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c93_logisticFeatures.pdf", width = 20, height = 10);
ggplot(union_df_2, aes(x = Features, y = species, fill = factor(presence)))+
  geom_tile()+
  scale_fill_manual(values = color_vector) +
  theme(axis.text.y = element_text(size = 5,color = "black"))+
  theme(axis.text.x= element_text(size=3, color="black", angle= 90))+
  geom_tile(color = "black")+
  facet_grid( catVector ~ function_type , scales="free",space="free")+
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8))+
  labs(x = "Functional Features", y = "Species")
dev.off()



