
# the code will select the important features from each functional profile and combine them into a single funcitonal profile.

allFunctionOutput_Folder= "G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\c102_funtionalProfileOutput\\"

allOutputs= list.files(allFunctionOutput_Folder)

for(outputs in allOutputs)
{
  filePath= paste0(allFunctionOutput_Folder,outputs)
  load(filePath)
}


load("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\modified_functional_prof_core.RData")

load("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\SpeciesScores_FunctionalAnalysis.RData")

load("G:\\My Drive\\HACK paper\\CellReports_revision\\Functional_Analysis_Overall_Output\\combined_metabolite_map.RData")

#selecting the top features for each function using random forest function feature importance.
#-------------------------------------------------------------------------------------------------------------------------------------

select_cazy <- names(which(cumsum(rf_cazy$importance[order(rf_cazy$importance),])/cumsum(rf_cazy$importance[order(rf_cazy$importance),])[length(rf_cazy$importance)]>=0.05))
# making the column names compatible with the random forest model.
select_cazy= sub("compatible__","",select_cazy)

# replace . with 4 underscores
select_cazy= gsub("____","\\.",select_cazy)

# replace - with 3 underscores
select_cazy= gsub("___","-",select_cazy)

#------------------------------------------------
select_kegg_module <- names(which(cumsum(rf_kegg_module$importance[order(rf_kegg_module$importance),])/cumsum(rf_kegg_module$importance[order(rf_kegg_module$importance),])[length(rf_kegg_module$importance)]>=0.05))
# making the column names compatible with the random forest model.
select_kegg_module= sub("compatible__","",select_kegg_module)

# replace . with 4 underscores
select_kegg_module= gsub("____","\\.",select_kegg_module)

# replace - with 3 underscores
select_kegg_module= gsub("___","-",select_kegg_module)

#------------------------------------------------
select_kegg_reaction <- names(which(cumsum(rf_kegg_reaction$importance[order(rf_kegg_reaction$importance),])/cumsum(rf_kegg_reaction$importance[order(rf_kegg_reaction$importance),])[length(rf_kegg_reaction$importance)]>=0.05))

# making the column names compatible with the random forest model.
select_kegg_reaction= sub("compatible__","",select_kegg_reaction)

# replace . with 4 underscores
select_kegg_reaction= gsub("____","\\.",select_kegg_reaction)

# replace - with 3 underscores
select_kegg_reaction= gsub("___","-",select_kegg_reaction)

#------------------------------------------------
select_EC_core <- names(which(cumsum(rf_EC_core$importance[order(rf_EC_core$importance),])/cumsum(rf_EC_core$importance[order(rf_EC_core$importance),])[length(rf_EC_core$importance)]>=0.05))
# making the column names compatible with the random forest model.
select_EC_core= sub("compatible__","",select_EC_core)

# replace . with 4 underscores
select_EC_core= gsub("____","\\.",select_EC_core)

# replace - with 3 underscores
select_EC_core= gsub("___","-",select_EC_core)


#------------------------------------------------
select_metabolite_map <- names(which(cumsum(rf_metabolite$importance[order(rf_metabolite$importance),])/cumsum(rf_metabolite$importance[order(rf_metabolite$importance),])[length(rf_metabolite$importance)]>=0.05))

# making the column names compatible with the random forest model.
select_metabolite_map= sub("compatible__","",select_metabolite_map)

# replace . with 4 underscores
select_metabolite_map= gsub("____","\\.",select_metabolite_map)

# replace - with 3 underscores
select_metabolite_map= gsub("___","-",select_metabolite_map)

#------------------------------------------------
select_pfam <- names(which(cumsum(rf_pfam$importance[order(rf_pfam$importance),])/cumsum(rf_pfam$importance[order(rf_pfam$importance),])[length(rf_pfam$importance)]>=0.05))

# making the column names compatible with the random forest model.
select_pfam= sub("compatible__","",select_pfam)

# replace . with 4 underscores
select_pfam= gsub("____","\\.",select_pfam)

# replace - with 3 underscores
select_pfam= gsub("___","-",select_pfam)

#------------------------------------------------
select_BiGG <- names(which(cumsum(rf_BiGG$importance[order(rf_BiGG$importance),])/cumsum(rf_BiGG$importance[order(rf_BiGG$importance),])[length(rf_BiGG$importance)]>=0.05))

# making the column names compatible with the random forest model.
select_BiGG= sub("compatible__","",select_BiGG)

# replace . with 4 underscores
select_BiGG= gsub("____","\\.",select_BiGG)

# replace - with 3 underscores
select_BiGG= gsub("___","-",select_BiGG)

#------------------------------------------------
select_cog <- names(which(cumsum(rf_cog$importance[order(rf_cog$importance),])/cumsum(rf_cog$importance[order(rf_cog$importance),])[length(rf_cog$importance)]>=0.05))

# making the column names compatible with the random forest model.
select_cog= sub("compatible__","",select_cog)

# replace . with 4 underscores
select_cog= gsub("____","\\.",select_cog)

# replace - with 3 underscores
select_cog= gsub("___","-",select_cog)
#------------------------------------------------
select_kegg_ko_core <- names(which(cumsum(rf_kegg_ko_core$importance[order(rf_kegg_ko_core$importance),])/cumsum(rf_kegg_ko_core$importance[order(rf_kegg_ko_core$importance),])[length(rf_kegg_ko_core$importance)]>=0.05))

# making the column names compatible with the random forest model.
select_kegg_ko_core= sub("compatible__","",select_kegg_ko_core)

# replace . with 4 underscores
select_kegg_ko_core= gsub("____","\\.",select_kegg_ko_core)

# replace - with 3 underscores
select_kegg_ko_core= gsub("___","-",select_kegg_ko_core)

#------------------------------------------------

df_select_cazy <- 
  as.data.frame(cbind(cazy[intersect(rownames(cazy),rownames(SpeciesScores)),select_cazy],SpeciesScores[intersect(rownames(cazy),rownames(SpeciesScores)),4]))
colnames(df_select_cazy)[ncol(df_select_cazy)] <- "HACKScore"

df_select_kegg_module <- 
  as.data.frame(cbind(kegg_module[intersect(rownames(kegg_module),rownames(SpeciesScores)),select_kegg_module],SpeciesScores[intersect(rownames(kegg_module),rownames(SpeciesScores)),4]))
colnames(df_select_kegg_module)[ncol(df_select_kegg_module)] <- "HACKScore"

df_select_kegg_reaction <- 
  as.data.frame(cbind(kegg_reaction[intersect(rownames(kegg_reaction),rownames(SpeciesScores)),select_kegg_reaction],SpeciesScores[intersect(rownames(kegg_reaction),rownames(SpeciesScores)),4]))
colnames(df_select_kegg_reaction)[ncol(df_select_kegg_reaction)] <- "HACKScore"

df_select_EC_core <- 
  as.data.frame(cbind(EC_core[intersect(rownames(EC_core),rownames(SpeciesScores)),select_EC_core],SpeciesScores[intersect(rownames(EC_core),rownames(SpeciesScores)),4]))
colnames(df_select_EC_core)[ncol(df_select_EC_core)] <- "HACKScore"

df_select_metabolite_map <- 
  as.data.frame(cbind(combined_metabolite_map[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),select_metabolite_map],SpeciesScores[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),4]))
colnames(df_select_metabolite_map)[ncol(df_select_metabolite_map)] <- "HACKScore"

df_select_pfam <- 
  as.data.frame(cbind(pfam[intersect(rownames(pfam),rownames(SpeciesScores)),select_pfam],SpeciesScores[intersect(rownames(pfam),rownames(SpeciesScores)),4]))
colnames(df_select_pfam)[ncol(df_select_pfam)] <- "HACKScore"

df_select_BiGG <- 
  as.data.frame(cbind(BiGG[intersect(rownames(BiGG),rownames(SpeciesScores)),select_BiGG],SpeciesScores[intersect(rownames(BiGG),rownames(SpeciesScores)),4]))
colnames(df_select_BiGG)[ncol(df_select_BiGG)] <- "HACKScore"

df_select_cog <- 
  as.data.frame(cbind(cog[intersect(rownames(cog),rownames(SpeciesScores)),select_cog],SpeciesScores[intersect(rownames(cog),rownames(SpeciesScores)),4]))
colnames(df_select_cog)[ncol(df_select_cog)] <- "HACKScore"

df_select_kegg_ko_core <- 
  as.data.frame(cbind(kegg_ko_core[intersect(rownames(kegg_ko_core),rownames(SpeciesScores)),select_kegg_ko_core],SpeciesScores[intersect(rownames(kegg_ko_core),rownames(SpeciesScores)),4]))
colnames(df_select_kegg_ko_core)[ncol(df_select_kegg_ko_core)] <- "HACKScore"

df_select_features <- as.data.frame(cbind(cazy[intersect(rownames(cazy),rownames(SpeciesScores)),select_cazy],kegg_module[intersect(rownames(kegg_module),rownames(SpeciesScores)),select_kegg_module],kegg_reaction[intersect(rownames(kegg_reaction),rownames(SpeciesScores)),select_kegg_reaction],EC_core[intersect(rownames(EC_core),rownames(SpeciesScores)),select_EC_core],pfam[intersect(rownames(pfam),rownames(SpeciesScores)),select_pfam],BiGG[intersect(rownames(BiGG),rownames(SpeciesScores)),select_BiGG],cog[intersect(rownames(cog),rownames(SpeciesScores)),select_cog],kegg_ko_core[intersect(rownames(kegg_ko_core),rownames(SpeciesScores)),select_kegg_ko_core]))

save(list = ls()[1:30],file="G:/My Drive/HACK paper/CellReports_revision/Final_Submission/Final_v3/Scripts_GitHub/Section10_FunctionalInvestigation/DATA/FunctionalClassificationAnalysis_stage2_part1.RData")
save(list = ls()[31:69],file="G:/My Drive/HACK paper/CellReports_revision/Final_Submission/Final_v3/Scripts_GitHub/Section10_FunctionalInvestigation/DATA/FunctionalClassificationAnalysis_stage2_part2.RData")

































