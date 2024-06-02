# Dr. Tarnini Shankar Ghosh, Abhishek Goel
# 20-11-2023

# Writing the function to find the association between the abundances between control and disease samples.

mannWhitney_batch = function(x,y)
{
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








#----------------------------------------------------------------------------------------------------------------------------


# input for function
# 1) controlSampleList- It is the list of studies, each study contains the control sample_ids.
# 2) diseaseSampleList- It is the list of diseases, each element of the list represents a disease. Each element[disease] is a list that contains studies which contains the sample_ids corresponding to the disease.
# 3) AllCombinedSpProfile- It is the matrix having abundance information of all the samples. [row= sample ids, column= species]
# 4) uniqueDiseases= It is a vector of diseases that contains the name of the diseases on which you want to work.
# 5) species_to_check- These are the target species on which the analysis will be performed.
# 6) saveLocationString- This is the location where your list of intermediate output of all the diseases will be saved. The folder that you stating  must be there else error will be prodcued as the function doesn't create the folder while processing.
# 7) fileName: Intermediate output will be saved with this name.
# 8) normalizationRequired- This is the logical argument that takes TRUE or FALSE values, if TRUE is taken then your data will be normalized (total sum scailing) else not.
# 9) speciesGroupList- This is the list of the groups, each group contains the name of the species. This is optional, by default it will take the empty list.


targetIntermediateOutput= function(controlSampleList, diseaseSampleList, AllCombinedSpProfile, uniqueDiseases, species_to_check, saveLocationString, fileName, normalizationRequired, speciesGroupList)
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
	print(length(speciesGroupList))
	if(length(speciesGroupList)==0)
	{
		tempSpProfile= AllCombinedSpProfile[diseaseSpecificControlSamples, species_to_check]
		tempSpProfile= tempSpProfile[rowSums(tempSpProfile)>0,]
		diseaseSpecificControlSamples= rownames(tempSpProfile)
		
		
		tempSpProfile= AllCombinedSpProfile[diseaseSpecificDiseaseSamples,species_to_check]
		tempSpProfile= tempSpProfile[rowSums(tempSpProfile)>0,]
		diseaseSpecificDiseaseSamples= rownames(tempSpProfile)
	}
	
	# global_allControlSamples<<- union(global_allControlSamples, diseaseSpecificControlSamples)
	# global_allDiseaseSamples<<- union(global_allDiseaseSamples, diseaseSpecificDiseaseSamples)
	
    controlSpProfile= data.frame()
    diseaseSpProfile= data.frame()
    
    # if group of species is not considered then
    if(length(speciesGroupList)==0)
    {
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
    }
    
    else # if group species list is given then this will work.
    {
      print("working on groups")
      
      # working on control profile
      controlSpProfile= as.data.frame(matrix(nrow= length(speciesGroupList), ncol=length(diseaseSpecificControlSamples)))
      rownames(controlSpProfile)= names(speciesGroupList)
      colnames(controlSpProfile)= diseaseSpecificControlSamples
	  
	  print("initial control dimesions")
	  print(dim(controlSpProfile))
	  
      # print(length(diseaseSpecificControlSamples))
	  # print(length(diseaseSpecificDiseaseSamples))
	  
      for(group in names(speciesGroupList))
      {
        availableSpecies= intersect(speciesGroupList[[group]], colnames(AllCombinedSpProfile))
        # printed= paste0(group," total species= ", length(speciesGroupList[[group]]), " available species= ", length(availableSpecies))
        # print(printed)
        tempProfile= AllCombinedSpProfile[diseaseSpecificControlSamples,availableSpecies]
		tempProfile= tempProfile[rowSums(tempProfile)>0,]
        # tempProfile= tempProfile[rowSums(tempProfile>0),]
		
		printed= paste0("control_group dim after removing empty rows= ",dim(tempProfile))
		print(printed)
		
		
        if(normalizationRequired==TRUE)
        {
		  print("performing normalization")
          tempProfile= tempProfile/rowSums(tempProfile)  
          tempProfile[is.na(tempProfile)]=0
        }
        
        tempProfile= as.data.frame(t(tempProfile))
		printed= paste0("control_group dim after transposing the data= ",dim(tempProfile))
		print(printed)
		
		
        # if(ncol(tempProfile)==0)
		# {
			# print("there are no samples for this group of species")
			# controlSpProfile[group,]= 0
			# next
		# }
		
        cumulativeVector= apply(tempProfile,2,sum)
		print("cumulative sample count: ")
		print(length(cumulativeVector))
		
		print("checking the colnames of the tempProfile")
        print(head(colnames(tempProfile)))
		print(head(colnames(controlSpProfile)))
		
        controlSpProfile[group,]= cumulativeVector
		print("control dimensions after adding values of group specific samples")
		print(dim(controlSpProfile))
		
		print("next group----------------------------------------------------------------------------------"	)
      }
	  
      controlSpProfile[is.na(controlSpProfile)]=0
	  controlSpProfile= controlSpProfile[,colSums(controlSpProfile)>0]
	  
	  printed= paste0("control_group dim after adding samples of all the groups= ",dim(controlSpProfile))
	  print(printed)
		
      # working on disease profile
	  
	  print("################################   working on disease profile   #############################################")
	  
      diseaseSpProfile= as.data.frame(matrix(nrow= length(speciesGroupList), ncol=length(diseaseSpecificDiseaseSamples)))
      rownames(diseaseSpProfile)= names(speciesGroupList)
      colnames(diseaseSpProfile)= diseaseSpecificDiseaseSamples
      
	  print("initial disease dimesions")
	  print(dim(diseaseSpProfile))
	  
      for(group in names(speciesGroupList))
      {
        availableSpecies= intersect(speciesGroupList[[group]], colnames(AllCombinedSpProfile))
        # printed= paste0(group," total species= ", length(speciesGroupList[[group]]), " available species= ", length(availableSpecies))
        # print(printed)
		
        tempProfile= AllCombinedSpProfile[diseaseSpecificDiseaseSamples,availableSpecies]
        # tempProfile= tempProfile[rowSums(tempProfile>0),]
		tempProfile= tempProfile[rowSums(tempProfile)>0,]
		printed= paste0("control_group dim after removing empty rows= ",dim(tempProfile))
		print(printed)
		
        if(normalizationRequired==TRUE)
        {
          tempProfile= tempProfile/rowSums(tempProfile)  
          tempProfile[is.na(tempProfile)]=0
        }
        
        tempProfile= as.data.frame(t(tempProfile))
		
		printed= paste0("control_group dim after transposing the data= ",dim(tempProfile))
		print(printed)
		
        # checkTemp <<- tempProfile
		# if(ncol(tempProfile)==0)
		# {
			# diseaseSpProfile[group,]= 0
			# next
		# }
        cumulativeVector= apply(tempProfile,2,sum)
		print("cumulative sample count: ")
		print(length(cumulativeVector))
		
		print("checking the colnames of the tempProfile")
        print(head(colnames(tempProfile)))
		print(head(colnames(diseaseSpProfile)))
		
        # TempCumulativeVector<<- cumulativeVector
        # checkTempProfile <<- tempProfile
        diseaseSpProfile[group,]= cumulativeVector
		
		print("control dimensions after adding values of group specific samples")
		print(dim(diseaseSpProfile))
		print("next group----------------------------------------------------------------------------------"	)
		
      }
	  diseaseSpProfile[is.na(diseaseSpProfile)]=0
	  diseaseSpProfile= diseaseSpProfile[,colSums(diseaseSpProfile)>0]
	  
	  # these are the global variables, so make sure to create their empty vector outside the function.
	  AllControlSamples= union(AllControlSamples, colnames(controlSpProfile))
	  AllDiseaseSamples= union(AllDiseaseSamples, colnames(diseaseSpProfile))
	  
	  printed= paste0("disease_group dim after adding samples of all the groups= ",dim(diseaseSpProfile))
	  print(printed)
      
    }
    # print((controlSpProfile[1:5,1:5]))
    # print(diseaseSpProfile[1:5,1:5])
    # checkDisease <<-diseaseSpProfile
    # checkControl <<- controlSpProfile
	printed= paste0("disease df dimensions =  ", dim(diseaseSpProfile))
	print(printed)
	
	printed= paste0("control df dimensions =  ", dim(controlSpProfile))
	print(printed)
    intermediateOutput[[dis]]= mannWhitney_batch(diseaseSpProfile,controlSpProfile)
	
  }
  
  save(intermediateOutput, file= paste0(saveLocationString,"\\",fileName))
  return(intermediateOutput)
}
#----------------------------------------------------------------------------------------------------------------------------

# intermediateOuptut: generated from "targetIntermediateOutput" function.
# speciesGroupingVariable: 
# speciesGroupingColour
# rowDendogramed
# colDendogramed
# speciesGroupingVariable- A named vector with names as species and values as the group/category
# speciesGroupingColour- Colour for each group you selected.

associationOutput= function(intermediateOutput, uniqueDiseases, outputFolder, speciesGroupingVariable, speciesGroupingColour, rowDendogramed, colDendogramed)
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
  heatmap2 <- heatmap.2(as.matrix(t(meanBasedOutput)) , density= "none", trace= "none", Rowv=rowDendogramed, Colv= colDendogramed)
  # print(names(heatmap2))
  # env$heatmap= heatmap2
  
  heatmapValues= heatmap2$carpet
  save(heatmapValues, file= paste0(outputFolder,"\\heatmapCarpet.RData"))
  write.table(rownames(heatmapValues), paste0(outputFolder,"\\AssociationHeatmap_colnames.txt"))
  write.table(colnames(heatmapValues), paste0(outputFolder,"\\AssociationHeatmap_rownames.txt"))
  
  dev.off()
  
  print("getting the association heatmap")
  if(length(speciesGroupingVariable)!=0)
  {
    tryCatch({
      pdf(paste0(outputFolder,"\\AssociationHeatmap.pdf"), width = 20, height = 10)
      # par(mar=c(1,1,1,1))
      heatmap.2(as.matrix(t(meanBasedOutput)) , density="none", trace="none", Rowv=rowDendogramed, Colv= colDendogramed, lhei=c(0.1,10), lwid=c(0.1,10), margins=c(10,10), sepcolor="black", sepwidth=c(0,0), rowsep=c(0:ncol(meanBasedOutput)), colsep=c(0:nrow(meanBasedOutput)), ColSideColor= c(speciesGroupingColour)[as.factor(speciesGroupingVariable)],cexRow=1.5,srtCol= 89.5,col=c("#607d3b","#98bf64","white", "white","white","#b38b6d","#7b403b"))
      dev.off()
    }, error = function(e) {
      cat("Error in keystoneInfluence for species:", species, "\n")
      cat("Error message:", conditionMessage(e), "\n")
      # return("abhs")  # Return an error indicator if needed
    })
  }
  else
  {
    tryCatch({
      pdf(paste0(outputFolder,"\\AssociationHeatmap.pdf"), width = 20, height = 10)
      # par(mar=c(1,1,1,1))
      heatmap.2(as.matrix(t(meanBasedOutput)) , density="none", trace="none", Rowv=rowDendogramed, Colv= colDendogramed, lhei=c(0.1,10), lwid=c(0.1,10), margins=c(10,10), sepcolor="black", sepwidth=c(0,0), rowsep=c(0:ncol(meanBasedOutput)), colsep=c(0:nrow(meanBasedOutput)),cexRow=1.5,,srtCol= 89.5,col=c("#607d3b","#98bf64","white", "white","white","#b38b6d","#7b403b"))
      dev.off()
    }, error = function(e) {
      cat("Error in keystoneInfluence for species:", "\n")
      cat("Error message:", conditionMessage(e), "\n")
      # return("abhs")  # Return an error indicator if needed
    })
  }
  
  
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
  
  
  
  if(length(speciesGroupingVariable)!=0)
  {
    associationCount$taxaGroup= NA
    associationCount[names(speciesGroupingVariable),"taxaGroup"]= speciesGroupingVariable
  }
  
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
  
  # creating a box plot of associations
  if(length(speciesGroupingVariable)!=0)
  {
    plot= ggplot(associationCount, aes(x = taxaGroup, y = totalNegativeSignificant, fill = taxaGroup)) +
      geom_violin(alpha=0.6) +  # Use geom_violin instead of geom_boxplot
      labs(x = "taxaGroup", y = "totalNegativeSignificant") +
      scale_fill_manual(values = speciesGroupingColour, name = "taxaGroup") +
      stat_summary(fun.y="median", geom="point", size=2, color="black")+
      theme_bw()
    
    ggsave("totalNegativeSignificant.pdf", plot, path = outputFolder)
    
    print("for totalNegativeSignificant")
    # print(wilcox.test(associationCount[associationCount$taxaType=="MajorTaxa","totalNegativeSignificant"],
    #                   associationCount[associationCount$taxaType=="nonMajorTaxa","totalNegativeSignificant"]))
    
    # output= capture.output(wilcox.test(associationCount[associationCount$taxaGroup=="MajorTaxa","totalNegativeSignificant"],
                                        # associationCount[associationCount$taxaGroup=="nonMajorTaxa","totalNegativeSignificant"]))
    
    # writeLines(output, paste0(outputFolder,"\\Wilcox_totalNegativeSignificant.txt"))
    
    plot= ggplot(associationCount, aes(x = taxaGroup, y = totalPositiveSignificant, fill = taxaGroup)) +
      geom_violin(alpha=0.6) +  # Use geom_violin instead of geom_boxplot
      labs(x = "taxaGroup", y = "totalPositiveSignificant") +
      scale_fill_manual(values = speciesGroupingColour, name = "taxaGroup") +
      stat_summary(fun.y="median", geom="point", size=2, color="black")+
      theme_bw()
    
    ggsave("totalPositiveSignificant.pdf", plot, path = outputFolder)
    
    print("for totalPositiveSignificant")
    # output= capture.output(wilcox.test(associationCount[associationCount$taxaGroup=="MajorTaxa","totalPositiveSignificant"],
                                       # associationCount[associationCount$taxaGroup=="nonMajorTaxa","totalPositiveSignificant"]))
    
    # writeLines(output, paste0(outputFolder,"\\Wilcox_totalPositiveSignificant.txt"))
    
    plot= ggplot(associationCount, aes(x = taxaGroup, y = totalSignificant, fill = taxaGroup)) +
      geom_violin(alpha=0.6) +  # Use geom_violin instead of geom_boxplot
      labs(x = "taxaGroup", y = "totalSignificant") +
      scale_fill_manual(values = speciesGroupingColour, name = "taxaGroup") +
      stat_summary(fun.y="median", geom="point", size=2, color="black")+
      theme_bw()
    
    ggsave("totalSignificant.pdf", plot, path = outputFolder)
    
    
    print("for totalSignificant")
    # output= capture.output(wilcox.test(associationCount[associationCount$taxaGroup=="MajorTaxa","totalSignificant"],
                                       # associationCount[associationCount$taxaGroup=="nonMajorTaxa","totalSignificant"]))
    
    # writeLines(output, paste0(outputFolder,"\\Wilcox_totalSignificant.txt"))
  }
  # generating the line plots for heatmap2 based on total significant and non significnat study count.
  
  # for heatmap2
  heatmap2Df= associationCount
  heatmap2Df= heatmap2Df[rownames(heatmapValues),]
  
  heatmap2Df$species= rownames(heatmap2Df)
  heatmap2Df$species <- factor(heatmap2Df$species, levels = unique(heatmap2Df$species))
  
  # getting the names of only those species that will be labeled, if significant or non significnat 
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
  
  if("MajorCore" %in% speciesGroupingVariable)
  {
    check= heatmap2Df %>% 
      filter(marker=="yes") %>% 
      filter(taxaGroup=="MajorTaxa") %>% 
      filter(totalNegativeSignificant>totalPositiveSignificant)
    
    healthAssociatedMajorKeystone= sort(as.character(check$species))
    
    save(healthAssociatedMajorKeystone, file = paste0(outputFolder,"\\healthAssociatedMajorKeystone.RData"))
  }
  
  write_xlsx(heatmap2Df,paste0(outputFolder,"\\heatmap2Df.xlsx"))
  
  if("MajorCore" %in% speciesGroupingVariable)
  {
    table(heatmap2Df[heatmap2Df$marker=="yes","taxaType"])  
  }
  
  # melting the data
  df_long= heatmap2Df[,c("totalNegativeSignificant", "totalPositiveSignificant","species")] %>%
    pivot_longer(!species, names_to = "variable", values_to = "value")
  
  df_long2= heatmap2Df[,c("totalNegativeSignificant", "totalPositiveSignificant","marker")] %>%
    pivot_longer(!marker, names_to = "variable", values_to = "value")
  
  df_long$marker= df_long2$marker
  
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
  
  ggsave("linePlot.pdf", plot, path = outputFolder, width = 20, height= 10)
  dev.off()
  return() 
}


# input arguments
# 1) iteration- Number of iterations that you want to perform on the disease analysis]
# 2) allDiseases- Vector of all the diseases.
# 3) speciesGroupingVariable- A named vector with names as species and values as the group/category
# 4) speciesGroupingColour- Colour for each group you selected.
# 5) controlSampleList- This is the list containing studies as element, each study is a vector of sample_ids 
# 6) diseaseSampleList- This is the list containing sub lists with thier name as disease. Each such sublist contains the elements as studies which is a vector of sample ids.
# 7) MainOutputFolder- This is the main folder in which output of all the iterations will be saved. You have to create it manually.

diseaseAnalysisIterationPipeline= function(iteration, allDiseases, speciesGroupingVariable, speciesGroupingColour, controlSampleList, diseaseSampleList, MainOutputFolder, AllCombinedSpProfile, species_to_check, speciesGroupList= list())
{
  # it will store the information about what diseases are used in each iteration.
  iterationSpecificDiseaseInfo= list()
  
  for(repeation in 1:iteration)
  {
    printed= paste0("######################################################################################  ", repeation,
                    "  ######################################################################################")
    print(printed)
	
	# randomly picking 65% of diseases out of all 28 disesaes. 18 in this case.
	totalSelected_diseases= round((28*0.65),0)
	
    selectedIndex= sample(seq(1, 28), totalSelected_diseases, replace = FALSE)
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
         AllCombinedMetadata,
         file= paste0(outputFolder,"\\SampleInfo.RData"))
    

    intermediateOutput= targetIntermediateOutput(subControlSampleList, subDiseaseSampleList, AllCombinedSpProfile, uniqueDiseases, species_to_check, outputFolder, "mannWhitneyIntermediateOutput.RData", TRUE, speciesGroupList)
    # load("G:\\My Drive\\IIITD\\project\\control_runs\\output_files\\c74_disease_discoveryCohortOutput\\mannWhitneyIntermediateOutput.RData")
    
    # sometimes this part gives error on Rstudio, in such case move to R console to run it
    associationOutput(intermediateOutput, uniqueDiseases, outputFolder, speciesGroupingVariable, speciesGroupingColour, TRUE, TRUE)
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
}



