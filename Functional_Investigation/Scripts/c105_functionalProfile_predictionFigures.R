
library(ggplot2)

load("FunctionalClassificationAnalysis_stage2_part1.RData")
load("FunctionalClassificationAnalysis_stage2_part2.RData")

outOfBox_variables= ls()[grepl("df_rf",ls())]

outOfBox_variables= setdiff(outOfBox_variables, "df_rf_metaboliteMap")

for(dataframe in outOfBox_variables)
{
  df= get(dataframe)
  my_plot= ggplot(df, aes(x=Actual,
                          y= Predicted))+
    geom_point(aes(color= "#6c0000", size= 0.07))+
    theme_bw()+
    xlab("Actual")+
    ylab("Predicted")+
    labs(x= "Actual HACK Score", y= "Predicted HACK Score")+
    theme(legend.position = "none")+
    theme(axis.text = element_text(size= 14,color="black"))+
    geom_smooth(method= "lm", color= "black")+
    theme(axis.title.x = element_text(size = 14),  # Size for x-axis title
          axis.title.y = element_text(size = 14))
  ggsave(paste0("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\FuntionalOutput_Plots\\c105_",dataframe,".pdf"), plot = my_plot)
}

# gettin the violin plot for the bootstrap results for all the files
# creating dataframe for the correaltion values for each funcitonal profile

corDf= as.data.frame(matrix(nrow= 800, ncol= 2))
colnames(corDf)= c("functionalProfile","correlationValues")
corDf$functionalProfile <- c(rep("BiGG", 100),
                             rep("EC_core",100),rep("cazy",100),
                             rep("cog",100),rep("kegg_ko_core",100),
                             rep("kegg_module",100),
                             rep("kegg_reaction",100),rep("pfam",100))
corDf$correlationValues <- c(bootstrap_rf_BiGG$cor,
                             bootstrap_rf_EC_core$cor,
                             bootstrap_rf_cazy$cor,bootstrap_rf_cog$cor,
                             bootstrap_rf_kegg_ko_core$cor,bootstrap_rf_kegg_module$cor,
                             bootstrap_rf_kegg_reaction$cor,bootstrap_rf_pfam$cor)

library(tidyverse)

medianCor_df= corDf %>% 
  group_by(functionalProfile) %>% 
  summarise_all(median, na.rm= TRUE)

write.csv(medianCor_df, "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c105_medianCorValues.csv")

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
  scale_y_continuous(limits = c(0, 1)) + 
  theme(plot.margin = margin(t = 50, r = 10, b = 10, l = 10, unit = "pt"))

ggsave(paste0("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\FuntionalOutput_Plots\\c105_bootstrap_violinPlot",".pdf"), plot = my_plot, height = 10, width= 25)


