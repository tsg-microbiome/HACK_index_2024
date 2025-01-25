#### Its the revised code of c73_diseaseAnalysisFunction.R (This is for UrbanRuralMixed only.)

library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gplots)
library(vegan)
library(ade4)
library(adegraphics)
library(writexl)

####################################################################

mannWhitney_batch = function(x,y){
  print("in man funciton")
  p_array <- NULL;
  medianDirection <- NULL;
  meanDirection <- NULL;
  
  z <- intersect(rownames(x),rownames(y));
  # print(rownames(x))
  for(i in 1:length(z))
  {
    printed= paste0(i," out of ", length(z))
    # print(printed)
    printed= paste0(i," ",z[i])
    # print(printed)
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
    
    medianDirection[i] <- ifelse(median(as.numeric(x[z[i],]),na.rm=TRUE) > median(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(median(as.numeric(x[z[i],]),na.rm=TRUE) < median(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    # print(medianDirection[i])
    meanDirection[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    # print(meanDirection[i])
    i <- i + 1;
  }
  out <- as.data.frame(cbind(medianDirection, meanDirection, p_array,p.adjust(p_array, method="fdr")));
  colnames(out)[4]= "q_array"
  
  
  rownames(out) <- z;
  # NA values will be either in p-value column or q-value column, if we take them as 1 it will still be insignificant and won't be considered.
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  # print(head(out))
  return(t(out));
}


#####################################################################

targetIntermediateOutput= function(controlSampleList, diseaseSampleList, AllCombinedSpProfile, uniqueDiseases, species_to_check, saveLocationString, fileName, normalizationRequired)
{
  # making sure that all the diseases from uniqueDiseases are present in diseaseSampleList
  if(length(setdiff(uniqueDiseases,names(diseaseSampleList)))>0)
  {
    print("The sample ids are not available for all the disesases you took please make sure diseaseSampleList contains all the uniqueDisease.")
    return()
  }
  
  # making sure that for all the studies taken for the disease profile their control abundance profile counterpart is available.
  
  allStudies= c()
  for(dis in names(diseaseSampleList))
  {
    allStudies= union(allStudies, names(diseaseSampleList[[dis]]))
  }
  
  if(length(setdiff(allStudies,names(controlSampleList)))>0)
  {
    print("all the studies that are taken for disease abundance profile, their control abundance profile counter part is not available")
    return()
  }
  
  
  # the list will contain the output for all the selected diseases.
  intermediateOutput= list() 
  
  cnt=1
  for(dis in uniqueDiseases)
  {
    printed= paste0(cnt," out of ", length(uniqueDiseases))
    print(printed)
    cnt= cnt+1
    
    # for the given disease taking all the disease and control samples from all the studies in which the disease was present.
    diseaseSpecificStudies= names(diseaseSampleList[[dis]])
    
    diseaseSpecificControlSamples= c()
    diseaseSpecificDiseaseSamples= c()
    
    for(study in diseaseSpecificStudies)
    {
      diseaseSpecificControlSamples= c(diseaseSpecificControlSamples, controlSampleList[[study]])
    }
    
    for(study in diseaseSpecificStudies)
    {
      diseaseSpecificDiseaseSamples= c(diseaseSpecificDiseaseSamples, diseaseSampleList[[dis]][[study]])
    }
    
    # removing empty rows from the control and disease species profile after taking 196 species.
    # for our disease analysis there are no empty rows, all the samples were considered.
    
    # global_allControlSamples<<- union(global_allControlSamples, diseaseSpecificControlSamples)
    # global_allDiseaseSamples<<- union(global_allDiseaseSamples, diseaseSpecificDiseaseSamples)
    
    controlSpProfile= data.frame()
    diseaseSpProfile= data.frame()
    
    print("working on specific species vector")
    
    # checking if there any species that are absent in the AllCombinedSpProfile, if yes then they have to be added with 0 values
    absentSpecies= setdiff(species_to_check,colnames(AllCombinedSpProfile))
    
    if(length(absentSpecies)>0)
    {
      print("some of the species_to_check species are absent in the species profile matrix")
      # print(absentSpecies)
      for(species in absentSpecies)
      {
        AllCombinedSpProfile[[species]]= with(AllCombinedSpProfile,0)
        printed= paste0("absent species is ", species)
        print(printed)
        
        # making sure that the species got added with the 0 values.
        print(unique(AllCombinedSpProfile[,species]))
      }
      
    }
    
    controlSpProfile= AllCombinedSpProfile[diseaseSpecificControlSamples, species_to_check]
    controlSpProfile= controlSpProfile[rowSums(controlSpProfile)>0,]
    # print(head(controlSpProfile))
    # normalizing the sp profile
    if(normalizationRequired==TRUE)
    {
      print("normalizing the control data")
      controlSpProfile= controlSpProfile/rowSums(controlSpProfile)  
      print(unique(rowSums(controlSpProfile)))
    }
    
    # transposing the data to make it compatible with the mannWhitney_batch funciton.
    controlSpProfile= as.data.frame(t(controlSpProfile))
    
    diseaseSpProfile= AllCombinedSpProfile[diseaseSpecificDiseaseSamples, species_to_check]
    diseaseSpProfile= diseaseSpProfile[rowSums(diseaseSpProfile)>0,]
    # normalizing the sp profile
    if(normalizationRequired==TRUE)
    {
      print("normalizing the disease data")
      diseaseSpProfile= diseaseSpProfile/rowSums(diseaseSpProfile)  
      print(unique(rowSums(diseaseSpProfile)))
    }
    
    # transposing the data to make it compatible with the mannWhitney_batch funciton.
    diseaseSpProfile= as.data.frame(t(diseaseSpProfile))
    
    printed= paste0("disease df dimensions =  ", dim(diseaseSpProfile))
    print(printed)
    
    printed= paste0("control df dimensions =  ", dim(controlSpProfile))
    print(printed)
    intermediateOutput[[dis]]= mannWhitney_batch(diseaseSpProfile,controlSpProfile)
    
  }
  
  save(intermediateOutput, file= paste0(saveLocationString,"\\",fileName))
  return(intermediateOutput)
}


#####################################################################

associationOutput= function(intermediateOutput, uniqueDiseases, outputFolder, rowDendogramed, colDendogramed)
{
  # creating the empty dataframe that contains information about the selected species(speceies_to_check) and the total number of diseases. It will contain association
  # value for all the species for each disease.
  meanBasedOutput= as.data.frame(matrix(nrow= nrow(intermediateOutput[[1]]), ncol= length(intermediateOutput)))
  rownames(meanBasedOutput)= rownames(intermediateOutput[[1]])
  colnames(meanBasedOutput)= names(intermediateOutput)
  
  cnt=1
  for(output in colnames(meanBasedOutput))
  {
    
    # getting the categorical directions as per q_value, p_value and meanDirection
    vector= c()
    for(rows in 1:nrow(meanBasedOutput))
    {
      if(intermediateOutput[[output]][rows,"meanDirection"]==0)
      {
        vector= c(vector,0)
      }
      else if(intermediateOutput[[output]][rows,"q_array"]<= 0.15)
      {
        vector= c(vector, intermediateOutput[[output]][rows,"meanDirection"]*3)
      }
      else if(intermediateOutput[[output]][rows,"p_array"]<= 0.05)
      {
        vector= c(vector, intermediateOutput[[output]][rows,"meanDirection"]*2)
      }
      else
      {
        vector= c(vector, intermediateOutput[[output]][rows,"meanDirection"]*1)
      }
      
    }
    meanBasedOutput[,output]= vector
    
  }
  
  print("getting the heatmap carpet")
  pdf(file = paste0(outputFolder, "\\heatmap.pdf"), width = 20, height = 15)
  heatmap2 <- heatmap.2(as.matrix(t(meanBasedOutput)) , density= "none", trace= "none", Rowv=rowDendogramed, Colv= colDendogramed)
  # print(names(heatmap2))
  # env$heatmap= heatmap2
  
  heatmapValues= heatmap2$carpet
  save(heatmapValues, file= paste0(outputFolder,"\\heatmapCarpet.RData"))
  write.table(rownames(heatmapValues), paste0(outputFolder,"\\AssociationHeatmap_colnames.txt"))
  write.table(colnames(heatmapValues), paste0(outputFolder,"\\AssociationHeatmap_rownames.txt"))
  
  dev.off()
  
  # getting the association count
  associationCount= data.frame("a")
  # checkMeanOutput <<- meanBasedOutput
  # it will get associations data for all 196 species instead of filtered one 
  # tempMeanOuptut<<- meanBasedOutput
  for(species in rownames(meanBasedOutput))
  {
    # print(species)
    currDf= as.data.frame(t(data.frame(table(unlist(meanBasedOutput[species,])))))
    columnNames= currDf[1,]
    currDf= as.data.frame(currDf[-1,])
    colnames(currDf)= columnNames
    currDf$species= species
    
    associationCount= merge(associationCount,currDf, all.x= TRUE, all.y= TRUE)
  }
  
  index = which(colnames(associationCount)=="X.a.")
  print(colnames(associationCount))
  associationCount= associationCount[,-index]
  
  species= associationCount$species
  
  index = which(colnames(associationCount)=="species")
  associationCount= associationCount[,-index]
  
  associationCount= as.data.frame(apply(associationCount,2,as.numeric))
  rownames(associationCount)= species
  
  associationCount[is.na(associationCount)]=0
  print("checkig if all the columns are present or not")
  for(columns in c("-3","-2","2","3","-1","0","1"))
  {
    if(!columns %in% colnames(associationCount))
    {
      associationCount[[columns]]= with(associationCount,0)  
      print(columns)
    }
    
    # making sure that the species got added with the 0 values.
    # print(unique(associationCount[,columns]))
  }
  # tempAssociation <<- associationCount
  
  associationCount$totalSignificant= rowSums(associationCount[,c("-3","3")])
  associationCount$totalNegativeSignificant= associationCount[,"-3"]
  associationCount$totalPositiveSignificant= associationCount[,"3"]
  
  associationCount$totalNonSignificant= rowSums(associationCount[,c("-1","0","1")])
  associationCount= associationCount[,c("-3","-2","-1","0","1","2","3","totalSignificant","totalNegativeSignificant",
                                        "totalPositiveSignificant","totalNonSignificant")]
  
  associationCount$total= length(uniqueDiseases)
  
  associationCount$scaledDifference= (associationCount$totalNegativeSignificant - associationCount$totalPositiveSignificant)/associationCount$total
  
  rank_scale=function(x)
  {
    # x <- rank(x);
    y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
    y <- ifelse(is.nan(y),0,y)
    return(y);
  }
  
  
  associationCount$ranked_scaledDifference= rank_scale(associationCount$scaledDifference)
  
  write_xlsx(associationCount, paste0(outputFolder,"\\associationCount.xlsx"))
  
  # for heatmap2
  heatmap2Df= associationCount
  heatmap2Df= heatmap2Df[rownames(heatmapValues),]
  
  heatmap2Df$species= rownames(heatmap2Df)
  heatmap2Df$species <- factor(heatmap2Df$species, levels = unique(heatmap2Df$species))
  
  # getting the names of only those species that will be labeled, if significant or non significant 
  # contribute to >=70% then it will be labeled
  heatmap2Df$marker= NA
  View(heatmap2Df)
  counter=1
  
  for(species in rownames(heatmap2Df))
  {
    printed= paste0(counter," ", species)
    print(printed)
    counter= counter+1
    
    if(heatmap2Df[species,"totalSignificant"]==0)
    {
      heatmap2Df[species,"marker"]= "no"
    }
    
    else if(heatmap2Df[species,"totalSignificant"]/heatmap2Df[species,"total"]<0.30)
    {
      heatmap2Df[species,"marker"]= "no"
    }
    
    else if(heatmap2Df[species,"totalNegativeSignificant"]/heatmap2Df[species,"totalSignificant"]>=0.70)
    {
      heatmap2Df[species,"marker"]= "yes"
    }
    else if(heatmap2Df[species,"totalPositiveSignificant"]/heatmap2Df[species,"totalSignificant"]>=0.70)
    {
      heatmap2Df[species,"marker"]= "yes"
    }
    else
    {
      heatmap2Df[species,"marker"]= "no"
    }
    
  }
  
  heatmap2Df$diseaseAssociation= NA
  index= which(heatmap2Df$marker=="yes" & heatmap2Df$totalNegativeSignificant> heatmap2Df$totalPositiveSignificant)
  
  heatmap2Df[index, "diseaseAssociation"]= "Negative Association With Disease"
  
  index= which(heatmap2Df$marker=="yes" & heatmap2Df$totalNegativeSignificant < heatmap2Df$totalPositiveSignificant)
  heatmap2Df[index, "diseaseAssociation"]= "Positive Association With Disease"
  
  # saving the selected major taxa that are have yes in marker column and have more negatively
  # associated studies
  
  
  write_xlsx(heatmap2Df,paste0(outputFolder,"\\heatmap2Df.xlsx"))
  
  # melting the data
  df_long= heatmap2Df[,c("totalNegativeSignificant", "totalPositiveSignificant","species")] %>%
    pivot_longer(!species, names_to = "variable", values_to = "value")
  
  df_long2= heatmap2Df[,c("totalNegativeSignificant", "totalPositiveSignificant","marker")] %>%
    pivot_longer(!marker, names_to = "variable", values_to = "value")
  
  df_long$marker= df_long2$marker
  pdf(file = paste0(outputFolder, "/linePlot.pdf"), width = 20, height = 10)
  plot= ggplot(df_long,
               aes(x= species,
                   y= value, group=variable))+
    geom_line(aes(color= variable),alpha = 0.5,size=0.5)+
    geom_point() +
    scale_color_manual(values=c("green", "red"))+
    xlab("Species")+
    ylab("sutdy count")+
    theme(axis.text = element_text(color= "black"))+
    theme_bw()+
    theme(axis.text.x= element_text(angle= 90, vjust = 0.5, size=6))
  
  dev.off()
  return() 
}

######################################################################

diseaseAnalysisIterationPipeline= function(iteration, allDiseases, controlSampleList, diseaseSampleList, MainOutputFolder, AllCombinedSpProfile, species_to_check)
{
  # it will store the information about what diseases are used in each iteration.
  iterationSpecificDiseaseInfo= list()
  
  for(repeation in 1:iteration)
  {
    printed= paste0("######################################################################################  ", repeation,
                    "  ######################################################################################")
    print(printed)
    
    # randomly picking 65% of diseases out of all 28 disesaes. 18 in this case.
    totalSelected_diseases= round((6*0.65),0)
    
    selectedIndex= sample(seq(1, 6), totalSelected_diseases, replace = FALSE)
    uniqueDiseases= allDiseases[selectedIndex]
    iterationSpecificDiseaseInfo[[repeation]]= uniqueDiseases
    
    # getting all the disease samples for the selected diseases.
    AllDiseaseSamples= c()
    iterationSpecificStudies= c()
    
    for(dis in uniqueDiseases)
    {
      for(study in names(diseaseSampleList[[dis]]))
      {
        AllDiseaseSamples= union(AllDiseaseSamples, diseaseSampleList[[dis]][[study]]) 
      }
      iterationSpecificStudies= union(iterationSpecificStudies, names(diseaseSampleList[[dis]]))
    }
    
    # getting all the control samples for the selected diseases.
    AllControlSamples= c()
    for(study in iterationSpecificStudies)
    {
      AllControlSamples= c(AllControlSamples, controlSampleList[[study]])
    }
    
    # taking subset for control and disease samples as per the current iteration.
    subControlSampleList= controlSampleList[iterationSpecificStudies]
    subDiseaseSampleList= diseaseSampleList[uniqueDiseases]
    
    # removing the empty rows after filtering 196 species
    tempSpProfile= AllCombinedSpProfile[AllControlSamples, species_to_check]
    tempSpProfile= tempSpProfile[rowSums(tempSpProfile)>0,]
    AllControlSamples= rownames(tempSpProfile)
    
    
    tempSpProfile= AllCombinedSpProfile[AllDiseaseSamples,species_to_check]
    tempSpProfile= tempSpProfile[rowSums(tempSpProfile)>0,]
    AllDiseaseSamples= rownames(tempSpProfile)
    
    
    
    dir.create(paste0(MainOutputFolder,"\\",repeation))
    outputFolder= paste0(MainOutputFolder,"\\",repeation)
    
    # there will be difference in the count of the samples in the list and the sample vectors as we removed some of the samples from control and disease after filtering 196 species.
    save(subControlSampleList,
         subDiseaseSampleList,
         AllControlSamples,
         AllDiseaseSamples,
         file= paste0(outputFolder,"\\SampleInfo.RData"))
    
    
    intermediateOutput= targetIntermediateOutput(subControlSampleList, subDiseaseSampleList, AllCombinedSpProfile, uniqueDiseases, species_to_check, outputFolder, "mannWhitneyIntermediateOutput.RData", TRUE)
    
    associationOutput(intermediateOutput, uniqueDiseases, outputFolder, TRUE, TRUE)
  }
  
  # stroing the results of all iterations in the list.
  allIterationFolders= list.files(MainOutputFolder)
  
  iterationOutputDf_direction= as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_direction)= species_to_check
  colnames(iterationOutputDf_direction)= allIterationFolders
  
  iterationOutputDf_divisionScore= as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_divisionScore)= species_to_check
  colnames(iterationOutputDf_divisionScore)= allIterationFolders
  
  iterationOutputDf_subtractionScore= as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_subtractionScore)= species_to_check
  colnames(iterationOutputDf_subtractionScore)= allIterationFolders
  
  iterationOutputDf_ranked_scaledDifference = as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_ranked_scaledDifference)= species_to_check
  colnames(iterationOutputDf_ranked_scaledDifference)= allIterationFolders
  
  iterationOutputDf_scaledDifference = as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_scaledDifference)= species_to_check
  colnames(iterationOutputDf_scaledDifference)= allIterationFolders
  
  
  for(folder in allIterationFolders)
  {
    heatmapInfo= read_excel(paste0(MainOutputFolder,"\\",folder,"\\heatmap2Df.xlsx"))
    heatmapInfo= as.data.frame(heatmapInfo)
    rownames(heatmapInfo)= heatmapInfo$species
    
    iterationOutputDf_ranked_scaledDifference[species_to_check,folder]= heatmapInfo[species_to_check,"ranked_scaledDifference"]
    iterationOutputDf_scaledDifference[species_to_check,folder]= heatmapInfo[species_to_check,"scaledDifference"]
  }
  
  
  outputList= list()
  outputList[["iterationSpecificDiseaseNames"]]= iterationSpecificDiseaseInfo
  outputList[["iterationOutputDf_ranked_scaledDifference"]]= iterationOutputDf_ranked_scaledDifference
  outputList[["iterationOutputDf_scaledDifference"]]= iterationOutputDf_scaledDifference
  
  return(outputList) 
  dev.off()
}

##################################################################################################


