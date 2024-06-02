# Dr. Tarini Shankar Ghosh, Abhishek Goel
# 10-01-2024

# Input Files:
# 1) FunctionalClassificationAnalysis_stage2.RData

# Output Files:
# 1) c105_images

library(ggplot2)

load("G:\\My Drive\\IIITD\\HACK_project\\control_runs\\input_files\\FunctionalProfiling/FunctionalClassificationAnalysis_stage2.RData")

outOfBox_variables= ls()[grepl("df_rf",ls())]

outOfBox_variables= setdiff(outOfBox_variables, "df_rf_metaboliteMap")

for(dataframe in outOfBox_variables)
{
  df= get(dataframe)
  # pdf(file = paste0("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\images\\c105_images\\c94_",dataframe,".pdf"))
  my_plot= ggplot(df, aes(x=Actual,
                          y= Predicted))+
    geom_point(aes(color= "#6c0000", size= 0.1))+
    theme_bw()+
    xlab("Actual")+
    ylab("Predicted")+
    labs(x= "Actual HACK Score", y= "Predicted HACK Score")+
    theme(legend.position = "none")+
    theme(axis.text = element_text(size= 24,color="black"))+
    geom_smooth(method= "lm", color= "black")+
    theme(axis.title.x = element_text(size = 30),  # Size for x-axis title
          axis.title.y = element_text(size = 30))
  ggsave(paste0("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\images\\c105_images\\c105_",dataframe,".pdf"), plot = my_plot)
  # dev.off()
}

# gettin the violin plot for the bootstrap results for all the files
# creating dataframe for the correaltion values for each funcitonal profile

corDf= as.data.frame(matrix(nrow= 800, ncol= 2))
colnames(corDf)= c("functionalProfile","correlationValues")

library(tidyverse)

medianCor_df= corDf %>% 
  group_by(functionalProfile) %>% 
  summarise_all(median, na.rm= TRUE)

write.csv(medianCor_df, "G:\\My Drive\\IIITD\\HACK_project\\control_runs\\Knowlede_Transfer\\FunctionalAnalysis\\outputFiles\\c105_medianCorValues.csv")

bootstrap_variables= ls()[grepl("bootstrap_rf", ls())]
bootstrap_variables= setdiff(bootstrap_variables, "bootstrap_rf_metaboliteMap")

index= 1
for(dfs in bootstrap_variables)
{
  currList= get(dfs)
  corValues= currList[["cor"]]
  corDf[index:(index+99),"functionalProfile"]= sub("bootstrap_rf_","",dfs)
  corDf[index:(index+99), "correlationValues"]= corValues
  
  index= index+100
}

my_plot=ggplot(corDf, aes(x = functionalProfile, y = correlationValues, fill= "red")) +
  geom_violin(alpha=0.6) +  # Use geom_violin instead of geom_boxplot
  labs(x = "Functional Profile", y = "Pearson Correlation Values")+
  stat_summary(fun.y="median", geom="point", size=2, color="black")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size= 24,color="black"))+
  theme(axis.text.x  = element_text(size= 24,color="black", angle= 90))+
  theme(axis.title.x = element_text(size = 30),  # Size for x-axis title
        axis.title.y = element_text(size = 30))+
  scale_y_continuous(limits = c(0, 1))

ggsave(paste0("G:\\My Drive\\IIITD\\HACK_project\\control_runs\\output_files\\images\\c105_images\\c105_bootstrap_violinPlot",".pdf"), plot = my_plot, height = 15, width= 25)


