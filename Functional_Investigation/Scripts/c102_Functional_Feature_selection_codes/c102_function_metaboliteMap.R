# Abhishek Goel
# 16-12-2023
# code name- function_metaboliteMap_abhishek.R

# functional profiling for MetaboliteMap


library(randomForest)
library(vegan)
library(ade4)
# library(randomForestExplainer)
library(ggplot2)
library(ggrepel)

dir_name <- "/storage32Tb/sourav/HACK_REVISION/"

load(paste0(dir_name,"SpeciesScores_FunctionalAnalysis.RData"))

load(paste0(dir_name,"combined_metabolite_map.RData"))

bootstrap_rf = function(data,metadata,metadata_column,iter,size)
{
	error_array <- NULL;
	cor_array <- NULL;
	common_rows <- intersect(rownames(data),rownames(metadata))
	data <- data[common_rows,]
	metadata <- metadata[common_rows,]
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	predictionPattern <- as.data.frame(matrix(NA,length(common_rows),iter));
	rownames(predictionPattern) <- common_rows
	colnames(predictionPattern) <- 1:iter
	for(i in 1:iter)
	{
		print(i);
		trainRows <- sample(rownames(data),as.integer(size*nrow(data)),replace=FALSE);
		testRows <- setdiff(rownames(data),trainRows);
		#print(data[trainRows,])
		rfTrain <- randomForest(metadata[trainRows,metadata_column]~.,data[trainRows,]);
		print(rfTrain)
		featureProfile[i,] <- sapply(colnames(data),function(x)(ifelse(x %in% rownames(rfTrain$importance),rfTrain$importance[x,],0)));
		predictionPattern[testRows,i] <- (predict(rfTrain,data[testRows,])-metadata[testRows,metadata_column])**2
		error_array[i] <- mean((predict(rfTrain,data[testRows,])-metadata[testRows,metadata_column])**2);
		cor_array[i] <- cor(predict(rfTrain,data[testRows,]),metadata[testRows,metadata_column]);
		print(error_array[i]);
		print(cor_array[i]);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(data);
	returnList <- list("err"=error_array,"cor"=cor_array,"featureProfile"=featureProfile,"predictionPattern"=predictionPattern);
}

# making the column names compatible with the random forest model.
colnames(combined_metabolite_map)= paste0("compatible__",colnames(combined_metabolite_map))

# replace . with 4 underscores
colnames(combined_metabolite_map)= gsub("\\.","____",colnames(combined_metabolite_map))

# replace - with 3 underscores
colnames(combined_metabolite_map)= gsub("-","___",colnames(combined_metabolite_map))


print("generating random forest model")
rf_metabolite <- randomForest(SpeciesScores[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),4]~.,combined_metabolite_map[intersect(rownames(combined_metabolite_map),rownames(SpeciesScores)),])

save.image("/storage32Tb/sourav/HACK_REVISION/metabolite/metabolite_map_temp.RData")

df_rf_metaboliteMap <- data.frame(Actual=rf_metabolite$y,Predicted=rf_metabolite$predicted)
save.image("/storage32Tb/sourav/HACK_REVISION/metabolite/metabolite_map_temp.RData")

pdf(file = "/storage32Tb/sourav/HACK_REVISION/metabolite/MetaboliteMap_FunctionCorPlot.pdf", width = 15, height = 15);
ggplot(df_rf_metaboliteMap,aes(x=Actual,y=Predicted))+geom_point(color="coral2",size=3)+geom_smooth(method="lm",color="blue")+theme_bw()+theme(axis.text.x=element_text(size=40),axis.text.y=element_text(size=40),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))
dev.off()

save.image("/storage32Tb/sourav/HACK_REVISION/metabolite/metabolite_map_temp.RData")


print("generating the bootstrap rf model")
bootstrap_rf_metaboliteMap <- bootstrap_rf(combined_metabolite_map,SpeciesScores,"HACKScore",100,0.5)
save.image("/storage32Tb/sourav/HACK_REVISION/metabolite/metabolite_map_temp.RData")

print("generating error rates")
err_rates_metaboliteMap <- apply(bootstrap_rf_metaboliteMap$predictionPattern,1,function(x)(mean(x[!is.na(x)])))

save.image("/storage32Tb/sourav/HACK_REVISION/metabolite/metabolite_map_temp.RData")

print("everything done")
save(rf_metabolite, df_rf_metaboliteMap, bootstrap_rf_metaboliteMap, err_rates_metaboliteMap, file= "/storage32Tb/sourav/HACK_REVISION/metabolite/MetaboliteMap_Output.RData")