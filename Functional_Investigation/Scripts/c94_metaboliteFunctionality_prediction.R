# the code will check if the hack score can predict the metabolity activity using the logisitic regression.

library(ggplot2)
library(ggrepel)
library(randomForest)
library(dplyr)

load("FunctionalProfiling\\combined_metabolite_map.RData")
load("MetaboliteMap_Output.RData")
load("SpeciesScores_FunctionalAnalysis.RData")


# selecting the top features from the rf_metabolite model
rf_metabolite_map= rf_metabolite
select_metabolite_map <- names(which(cumsum(rf_metabolite_map$importance[order(rf_metabolite_map$importance),])/cumsum(rf_metabolite_map$importance[order(rf_metabolite_map$importance),])[length(rf_metabolite_map$importance)]>=0.05))


# converting the names of features back into the orignial form.
# -------------------------------------------------------------------------------------------------------------------
select_metabolite_map= sub("compatible__","",select_metabolite_map)

# replace . with 4 underscores
select_metabolite_map= gsub("____","\\.",select_metabolite_map)

# replace - with 3 underscores
select_metabolite_map= gsub("___","-",select_metabolite_map)
# -------------------------------------------------------------------------------------------------------------------

metaboliteMap_functionalProfile= combined_metabolite_map[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),select_metabolite_map]

HACKScore= unlist(SpeciesScores[rownames(metaboliteMap_functionalProfile), "HACKScore"])

# creatingt the dataframe that will store the estimate and p-value for each of the metabolite function deduced using logistic regression.
logisticModel_Df= as.data.frame(matrix(nrow= ncol(metaboliteMap_functionalProfile), ncol= 2))
colnames(logisticModel_Df)= c("Estimate","p_value")
rownames(logisticModel_Df)= colnames(metaboliteMap_functionalProfile)


featureCnt=1

for(features in rownames(logisticModel_Df))
{
  print(featureCnt)
  featureCnt= featureCnt+1
  
  model <- glm(metaboliteMap_functionalProfile[,features] ~ HACKScore,
               family = binomial(link = "logit")) 
  
  model_sum= summary(model)
  
  p_value= model_sum$coefficients[2,4]
  estimate= model_sum$coefficients[2,1]
  
  logisticModel_Df[features,"p_value"]= p_value
  logisticModel_Df[features,"Estimate"]= estimate
}


# -------------------------------------------------------------------------------------------------------------------

# creating the volcano plot for the logistic features
logisticModel_Df$association= ifelse(logisticModel_Df$Estimate>0,"Positive","Negative")

logisticModel_Df$log_pValue= -log10(logisticModel_Df$p_value)

rownames(logisticModel_Df)= gsub("_"," ",rownames(logisticModel_Df))

# getting only those labels that we want on the plot
labels= rownames(logisticModel_Df[which(logisticModel_Df$p_value<=0.1),])

# making category so that species with different p-value range can be coloured differently
logisticModel_Df$category= NA

logisticModel_Df[which(logisticModel_Df$p_value<=0.05),"category"]= "main"
logisticModel_Df[which(logisticModel_Df$p_value>0.05 & logisticModel_Df$p_value<=0.1),"category"]= "marginal"
logisticModel_Df[which(logisticModel_Df$p_value>0.1),"category"]= "insignificant"

logisticModel_Df$category= paste0(logisticModel_Df$category, "__", logisticModel_Df$association)
logisticModel_Df$category= as.factor(logisticModel_Df$category)

pdf("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c94_logisticFeatures_volcanoPlot.pdf", width = 16, height = 7)
ggplot(logisticModel_Df, aes(x=Estimate,
                             y= log_pValue, color= category))+
  scale_color_manual(values = c("gray", "gray", "red4","blue4","red3","#5d98a8"))+
  geom_point(aes(size= 0.1))+
  scale_x_continuous(limits = c(-10, 10))+
  theme_bw()+
  xlab("Estimate")+
  ylab("-log(p-value)base10")+
  theme(legend.position = "none")+
  theme(legend.position="none")+
  theme(axis.text = element_text(size= 17,color="black"))+
  geom_label_repel(data= subset(logisticModel_Df, p_value<=0.1),aes(label= labels),
                   size= 5.5, segment.size= 0.2, force= 8.2, label.size = NA, max.overlaps = 20)+
  geom_hline(yintercept = 1, linetype = "dotted", size = 1)+
  geom_hline(yintercept = 1.301, linetype = "dotted", size = 1)+
  geom_vline(xintercept = 0, size = 1)
  
dev.off()


# -------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------
# predicting the HACKScore using the rf_selected random forest features
rownames(logisticModel_Df)= gsub(" ","_",rownames(logisticModel_Df))

rf_MetaboliteMap <- randomForest(SpeciesScores[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),4]~.,combined_metabolite_map[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),rownames(logisticModel_Df)])

cor(rf_MetaboliteMap$y, rf_MetaboliteMap$predicted, method= "spearman") # 0.53

# pearson correlation
cor(rf_MetaboliteMap$y, rf_MetaboliteMap$predicted) # 0.51

# -------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------

filtered_logisticModel_df= logisticModel_Df %>% 
  filter(p_value<=0.05)

# taking only the logistic selected model as they are better to predic the hack score
metaboliteMap_functionalProfile= combined_metabolite_map[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),rownames(filtered_logisticModel_df),]

# training random forest model to check if the metabolite can predict the hack socre
rf_logistic_selectedFeatures <- randomForest(SpeciesScores[intersect(rownames(SpeciesScores),rownames(metaboliteMap_functionalProfile)),4]~.,metaboliteMap_functionalProfile[intersect(rownames(SpeciesScores),rownames(metaboliteMap_functionalProfile)),])

# spearman correlation
cor.test(rf_logistic_selectedFeatures$y, rf_logistic_selectedFeatures$predicted, method= "spearman") #p-value = 1.666e-11, rho = 0.56

# pearson correlation
cor(rf_logistic_selectedFeatures$y, rf_logistic_selectedFeatures$predicted) # 0.54

# getting the mean square error for the predicted values
mean(((rf_logistic_selectedFeatures$y)^2) - ((rf_logistic_selectedFeatures$predicted)^2))

df= data.frame(rf_logistic_selectedFeatures$y, rf_logistic_selectedFeatures$predicted)
colnames(df)= c("actual", "predicted")
df$HACKScore= SpeciesScores[rownames(df),"HACKScore"]

df$category= ifelse(df$HACKScore>=0.5,"HACK", "Non HACK")

pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c94_MetabolitePrediction_plot.pdf")
ggplot(df, aes(x=actual,
               y= predicted))+
  geom_point(aes(color= "#6c0000", size= 0.1))+
  theme_bw()+
  xlab("Actual")+
  ylab("Predicted")+
  labs(x= "Actual HACK Score", y= "Predicted HACK Score")+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size= 16,color="black"))+
  geom_smooth(method= "lm", color= "black")
dev.off()

save(logisticModel_Df, filtered_logisticModel_df
     , file= "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c94_metabolite_selectedFeatures.RData")

supplementary_table20 <- combined_metabolite_map[rownames(combined_metabolite_map) %in% rownames(df_rf_metaboliteMap),]
write.xlsx(supplementary_table20,"supplementary_table20.xlsx",row.names = T)

tableS21_B <- logisticModel_Df[,c(1,2)]
tableS21_B <- tableS21_B[order(tableS21_B$p_value, decreasing = F),]
write.xlsx(tableS21_B,"tableS21_B.xlsx")

feature_imp <- as.data.frame(rf_metabolite_map$importance)
feature_imp <- as.data.frame(feature_imp[order(feature_imp$IncNodePurity,decreasing = T),,drop = F])
rownames(feature_imp) <- gsub("compatible__","",rownames(feature_imp))
write.xlsx(feature_imp,"feature_imp.xlsx",row.names = T)

