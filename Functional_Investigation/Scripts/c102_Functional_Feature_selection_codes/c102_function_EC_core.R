# Abhishek Goel
# 16-12-2023

# code name- c102_function_EC_core_abhishek.R

# -------------------------------------------------------------------------------------------------------------------
library(randomForest)
library(vegan)
library(ade4)
# library(randomForestExplainer)
library(ggplot2)
library(ggrepel)

dir_name <- "/storage32Tb/sourav/HACK_REVISION/"

load(paste0(dir_name,"SpeciesScores_FunctionalAnalysis.RData"))

load(paste0(dir_name,"modified_functional_prof_core.RData"))

# -------------------------------------------------------------------------------------------------------------------
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
# -------------------------------------------------------------------------------------------------------------------

df= EC_core
functionName= "EC_core"

# -------------------------------------------------------------------------------------------------------------------
# making the column names compatible with the random forest model.
colnames(df)= paste0("compatible__",colnames(df))

# replace . with 4 underscores
colnames(df)= gsub("\\.","____",colnames(df))

# replace - with 3 underscores
colnames(df)= gsub("-","___",colnames(df))

# -------------------------------------------------------------------------------------------------------------------
print("generating random forest model")
rf_model <- randomForest(SpeciesScores[intersect(rownames(df),rownames(SpeciesScores)),4]~.,df[intersect(rownames(df),rownames(SpeciesScores)),])

save.image(paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_temp.RData"))

rf_df <- data.frame(Actual=rf_model$y,Predicted=rf_model$predicted)
save.image(paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_temp.RData"))

pdf(file = paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_FunctionCorPlot.pdf"), width = 15, height = 15);
ggplot(rf_df,aes(x=Actual,y=Predicted))+geom_point(color="coral2",size=3)+geom_smooth(method="lm",color="blue")+theme_bw()+theme(axis.text.x=element_text(size=40),axis.text.y=element_text(size=40),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))
dev.off()

save.image(paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_temp.RData"))
# -------------------------------------------------------------------------------------------------------------------

print("generating the bootstrap rf model")
bootstrap_model <- bootstrap_rf(df,SpeciesScores,"HACKScore",100,0.5)
save.image(paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_temp.RData"))

print("generating error rates")
err_rates <- apply(bootstrap_model$predictionPattern,1,function(x)(mean(x[!is.na(x)])))

save.image(paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_temp.RData"))

print("everything done")
# -------------------------------------------------------------------------------------------------------------------

rf_EC_core= rf_model
df_rf_EC_core= rf_df
bootstrap_rf_EC_core= bootstrap_model
err_rates_EC_core= err_rates

save(rf_EC_core, df_rf_EC_core, bootstrap_rf_EC_core, err_rates_EC_core, file= paste0("/storage32Tb/sourav/HACK_REVISION/ec/",functionName,"_Output.RData"))