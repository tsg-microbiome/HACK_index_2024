print("Loading Libraries")
library(compositions)
library(igraph)
library(metap)
library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(gplots)
library(RColorBrewer)
library(sfsmisc)
library(pcaPP)
library(psych)
library(dendextend)
library(gplots)
library(RColorBrewer)
library(pROC)
library(randomForest)
library(gtools)
library(vegan)
library(ade4)
library(effsize)
library(ggplot2)
library(pcaPP)
library(ggrepel)
library(adegraphics)
library(ccrepe)

compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
	colnames(temp_meta) <- c("dataset","ri","ni")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		dat1 <- data[data[,grouping_variable]==group,var1]
		dat2 <- data[data[,grouping_variable]==group,var2]
		temp_meta[i,2] <- cor.fk(dat1,dat2)
		temp_meta[i,3] <- length(dat1)
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	return(res)
}

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
	detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
	rownames(detection_matrix) <- var1_list
	colnames(detection_matrix) <- grouping_list
	for(i in 1:length(var1_list))
	{
		var1 <- var1_list[i]
		for(j in 1:length(grouping_list))
		{
			group <- grouping_list[j]
			detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]>0))/length(data[data[,grouping_variable]==group,var1])
		}
	}
	return(detection_matrix)
}

compute_hedges_g <- function(data,var1_list,var2,grouping_variable,grouping_list)
{
	hedges_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
	colnames(hedges_matrix) <- grouping_list
	rownames(hedges_matrix) <- var1_list
	p_value_matrix <- data.frame(matrix(1,length(var1_list),length(grouping_list)))
	colnames(p_value_matrix) <- grouping_list
	rownames(p_value_matrix) <- var1_list
	for(i in 1:length(var1_list))
	{
		var1 <- var1_list[i]
		for(j in 1:length(grouping_list))
		{
			group <- grouping_list[j]
			data_group_Case <- data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1]
			data_group_Control <- data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1]
			hedges_matrix[i,j] <- as.numeric(effsize::cohen.d(data_group_Case,data_group_Control,hedges.correction=TRUE)$estimate)
			p_value_matrix[i,j] <- as.numeric(wilcox.test(data_group_Case,data_group_Control)$p.value)
		}
	}
	return_list = list("hedges"=hedges_matrix,"p_value"=p_value_matrix)
	return(return_list)
}

compute_meta_effsize <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),7))
	colnames(temp_meta) <- c("dataset","m1i","m2i","sd1i","sd2i","n1i","n2i")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		data_group_Case <- data[(data[,var2]==1)&(data[,grouping_variable]==group),var1]
		data_group_Control <- data[(data[,var2]==0)&(data[,grouping_variable]==group),var1]
		temp_meta[i,2] <- mean(data_group_Case)
		temp_meta[i,3] <- mean(data_group_Control)
		temp_meta[i,4] <- sd(data_group_Case)
		temp_meta[i,5] <- sd(data_group_Control)
		temp_meta[i,6] <- length(data_group_Case)
		temp_meta[i,7] <- length(data_group_Control)
		print(temp_meta)
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="SMD",m1i=m1i,m2i=m2i,sd1i=sd1i,sd2i=sd2i,n1i=n1i,n2i=n2i,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	#res$slabs <- rownames(temp_meta)
	res_list <- list("df_studies"=temp_meta,"model"=res)
	return(res_list)
}

compute_meta_effsize2 <- function(data,var1,list1,list2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),7))
	colnames(temp_meta) <- c("dataset","m1i","m2i","sd1i","sd2i","n1i","n2i")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		data_group_Case <- data[list1,var1]
		data_group_Control <- data[list2,var1]
		temp_meta[i,2] <- mean(data_group_Case)
		temp_meta[i,3] <- mean(data_group_Control)
		temp_meta[i,4] <- sd(data_group_Case)
		temp_meta[i,5] <- sd(data_group_Control)
		temp_meta[i,6] <- length(data_group_Case)
		temp_meta[i,7] <- length(data_group_Control)
		print(temp_meta)
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="SMD",m1i=m1i,m2i=m2i,sd1i=sd1i,sd2i=sd2i,n1i=n1i,n2i=n2i,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	#res$slabs <- rownames(temp_meta)
	res_list <- list("df_studies"=temp_meta,"model"=res)
	return(res_list)
}


compute_meta_lm <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		vec1 <- data[data[,grouping_variable]==group,var1]
		vec2 <- data[data[,grouping_variable]==group,var2]
		
		#print(paste0(group,",",length(vec1[abs(vec1)>0]),",",length(vec2[abs(vec2)>0])))
		if((length(vec1[abs(vec1)>0]) > 0)&&(length(vec2[abs(vec2)>0]) > 0))
		{
			#print(data[data[,grouping_variable]==group,c(var1,var2)])
			print(group)
			tryCatch(               
							expr = {    
										f <- as.formula(paste0(var1,"~",var2))
										temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
										summary_temp_rlm <- summary(temp_rlm)
										#print("Enter")
										print(summary_temp_rlm)
										temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
										temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
										temp_meta[i,4] <- 1
										temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
										temp_meta[i,6] <- sign(temp_meta[i,2])
										#levels <- unique(data[data[,grouping_variable]==group,var2])
										#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
										#temp <- ifelse(is.nan(temp),0,temp)
										#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
										#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
										#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
										#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
										#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
									
									},
							error = function(e){ 
										temp_meta[i,2] <- 0
										temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
										temp_meta[i,4] <- 1
										temp_meta[i,5] <- 1
										temp_meta[i,6] <- 1
										print(e)
										print("Error observed. Moving to next")
									},
							finally = {            
										print("finally Executed")
							}
						)
			
			
		}
		else
		{
			#temp_meta[i,2] <- 0
			#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			#temp_meta[i,4] <- 1
			#temp_meta[i,5] <- 1
			#temp_meta[i,6] <- 1
		}
	}
	temp_meta <- temp_meta[!is.na(temp_meta[,"ti"])&(temp_meta[,"ti"] != 0),]
	grouping_list <- temp_meta$dataset
	print(temp_meta)
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}


compute_meta_lm_single_adjust <- function(data,var1,var2,var3,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		f <- as.formula(paste0(var1,"~",var3,"+",var2))
		temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
		summary_temp_rlm <- summary(temp_rlm)
		temp_meta[i,2] <- summary_temp_rlm$coefficients[3,3]
		temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		temp_meta[i,4] <- 1
		temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
		temp_meta[i,6] <- sign(temp_meta[i,2])
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

compute_meta_lm_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
	return_out <- as.data.frame(matrix(NA,length(feature_list),10))
	rownames(return_out) <- feature_list
	colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","consistency","qval","dir")
	for(i in 1:length(feature_list))
	{
		species_name <- feature_list[i]
		tryCatch(               
							expr = {      
									print(species_name)
									temp_res <- compute_meta_lm(data,species_name,metadata_var,grouping_var,grouping_list)
									return_out[i,"beta"] <- temp_res$model$beta
									return_out[i,"pval"] <- temp_res$model$pval
									return_out[i,"ci.ub"] <- temp_res$model$ci.ub
									return_out[i,"ci.lb"] <- temp_res$model$ci.lb
									return_out[i,"tau2"] <- temp_res$model$tau2
									return_out[i,"QE"] <- temp_res$model$QE
									return_out[i,"QEp"] <- temp_res$model$QEp
									return_out[i,"consistency"] <- length(which(sign(temp_res$df_studies$di) == sign(as.numeric(temp_res$model$beta))))/nrow(temp_res$df_studies)		
									},
									error = function(e){ 
										return_out[i,"beta"] <- 0
										return_out[i,"pval"] <- 1
										return_out[i,"ci.ub"] <- 0
										return_out[i,"ci.lb"] <- 0
										return_out[i,"tau2"] <- 1
										return_out[i,"QE"] <- 1
										return_out[i,"QEp"] <- 1
										return_out[i,"consistency"] <- 0
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
		
	}
	return_out$qval <- p.adjust(return_out$pval,method="fdr")
	return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta)))
	return(return_out)
}

batch_rlm_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	df_est <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_est) <- grouping_list
	colnames(df_est) <- variable_group
	df_p_val <- as.data.frame(matrix(1,length(grouping_list),length(variable_group)))
	rownames(df_p_val) <- grouping_list
	colnames(df_p_val) <- variable_group
	for(i in 1:length(grouping_list))
	{
		study_name <- grouping_list[i]
		print(study_name)
		study_samples <- rownames(metadata[metadata[,grouping_feature] == study_name,])
		for(j in 1:length(variable_group))
		{
			tryCatch(               
						expr = {                     
									species_name <- variable_group[j]
									print(species_name)
									species_val <- data[study_samples,species_name]
									print(length(species_val[species_val>0]))
									if(length(species_val[species_val>0])>0)
									{
										temp_rlm <- rlm(data[study_samples,species_name]~metadata[study_samples,metadata_feature])
										summary_temp_rlm <- summary(temp_rlm)
										print(summary_temp_rlm)
										df_est[i,j] <- summary_temp_rlm$coefficients[2,3]
										df_p_val[i,j] <- f.robftest(temp_rlm)$p.value
									}
								},
								error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
								},
								finally = {            
										print("finally Executed")
								}
						)
				

			
		}
		
	}
	df_q_val <- apply(df_p_val,2,p.adjust)
	l_fisher <- p.adjust(apply(df_q_val,2,function(x)(sumlog(x)$p)),method="fdr")
	print("l_fisher generated")
	df_dir <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_dir) <- grouping_list
	colnames(df_dir) <- variable_group
	for(i in 1:length(grouping_list))
	{
		for(j in 1:length(variable_group))
		{
			df_dir[i,j] <- ifelse(df_q_val[i,j]<=0.10,3*sign(df_est[i,j]),ifelse(df_p_val[i,j]<=0.05,2*sign(df_est[i,j]),1*sign(df_est[i,j])))
		}
	}
	return_list <- list("est"=df_est,"p.value"=df_p_val,"q.value"=df_q_val,"fisher"=l_fisher,"dir"=df_dir)
	return(return_list)
}

batch_rlm_grouped_ab_adjust <- function(data,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	df_est <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_est) <- grouping_list
	colnames(df_est) <- variable_group
	df_p_val <- as.data.frame(matrix(1,length(grouping_list),length(variable_group)))
	rownames(df_p_val) <- grouping_list
	colnames(df_p_val) <- variable_group
	for(i in 1:length(grouping_list))
	{
		study_name <- grouping_list[i]
		print(study_name)
		study_samples <- rownames(data[data[,grouping_feature] == study_name,])
		for(j in 1:length(variable_group))
		{
			tryCatch(               
						expr = {                     
									species_name <- variable_group[j]
									diff_species_name <- paste0("diff_",species_name)
									print(species_name)
									species_val <- data[study_samples,species_name]
									print(length(species_val[species_val>0]))
									if(length(species_val[species_val>0])>0)
									{
										temp_rlm <- rlm(as.formula(paste0(metadata_feature,"~",diff_species_name,"+",species_name)),data=data[study_samples,])
										summary_temp_rlm <- summary(temp_rlm)
										print(summary_temp_rlm)
										df_est[i,j] <- summary_temp_rlm$coefficients[3,3]
										df_p_val[i,j] <- f.robftest(temp_rlm,var=species_name)$p.value
									}
								},
								error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
								},
								finally = {            
										print("finally Executed")
								}
						)
				

			
		}
		
	}
	df_q_val <- t(apply(df_p_val,1,function(x)(p.adjust(x,method="fdr"))))
	l_fisher <- p.adjust(apply(df_q_val,2,function(x)(sumlog(x)$p)),method="fdr")
	print("l_fisher generated")
	df_dir <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_dir) <- grouping_list
	colnames(df_dir) <- variable_group
	for(i in 1:length(grouping_list))
	{
		for(j in 1:length(variable_group))
		{
			df_dir[i,j] <- ifelse(df_q_val[i,j]<=0.15,3*sign(df_est[i,j]),ifelse(df_p_val[i,j]<=0.05,2*sign(df_est[i,j]),1*sign(df_est[i,j])))
		}
	}
	return_list <- list("est"=df_est,"p.value"=df_p_val,"q.value"=df_q_val,"fisher"=l_fisher,"dir"=df_dir)
	return(return_list)
}


batch_rem_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	common_rows <- intersect(rownames(data),rownames(metadata))
	df_est <- as.data.frame(matrix(0,length(variable_group),5))
	rownames(df_est) <- variable_group
	colnames(df_est) <- c("est","pval","qval","dir",metadata_feature)
	df_est$est <- 0
	df_est$pval <- 1
	df_est$qval <- 1
	for(j in 1:length(variable_group))
	{
			species_name <- variable_group[j]
			print(species_name)
			species_val <- data[common_rows,species_name]
			metadata_val <- metadata[common_rows,metadata_feature]
			metadata_grouping <- metadata[common_rows,grouping_feature]
			print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				tryCatch(               
							expr = {                     
									df_temp <- data.frame(sp=species_val,meta_val=metadata_val,meta_grp=metadata_grouping,row.names=common_rows)
									temp_res <- compute_meta_lm(df_temp,"sp","meta_val","meta_grp",grouping_list)
									print(temp_res$model)
									df_est[j,1] <- as.numeric(temp_res$model$beta)
									df_est[j,2] <- temp_res$model$pval	
		
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
								
			}
		
	}
	df_est[,3] <- p.adjust(df_est[,2],method="fdr")
	df_est[,4] <- ifelse(df_est[,3]<=0.10,3*sign(df_est[,1]),ifelse(df_est[,2]<=0.05,2*sign(df_est[,1]),1*sign(df_est[,1])))
	df_est[,5] <- metadata_feature
	return(df_est)
}

batch_rem2_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	common_rows <- intersect(rownames(data),rownames(metadata))
	df_est <- as.data.frame(matrix(0,length(variable_group),5))
	rownames(df_est) <- variable_group
	colnames(df_est) <- c("est","pval","qval","dir",metadata_feature)
	df_est$est <- 0
	df_est$pval <- 1
	df_est$qval <- 1
	for(j in 1:length(variable_group))
	{
			species_name <- variable_group[j]
			print(species_name)
			species_val <- data[common_rows,species_name]
			metadata_val <- metadata[common_rows,metadata_feature]
			metadata_grouping <- metadata[common_rows,grouping_feature]
			print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				tryCatch(               
							expr = {                     
									df_temp <- data.frame(sp=species_val,meta_val=as.factor(metadata_val),meta_grp=metadata_grouping,row.names=common_rows)
									temp_res <- compute_meta_effsize(df_temp,"sp","meta_val","meta_grp",grouping_list)
									print(temp_res$model)
									df_est[j,1] <- as.numeric(temp_res$model$beta)
									df_est[j,2] <- temp_res$model$pval	
		
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
								
			}
		
	}
	df_est[,3] <- p.adjust(df_est[,2],method="fdr")
	df_est[,4] <- ifelse(df_est[,3]<=0.10,3*sign(df_est[,1]),ifelse(df_est[,2]<=0.05,2*sign(df_est[,1]),1*sign(df_est[,1])))
	df_est[,5] <- metadata_feature
	return(df_est)
}

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

rem_network1 <- function(data,species_list,group_name,study_list)
{
        species_data <- data[,species_list]
        species_data$group <- data[,group_name]
        est_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
		rownames(est_matrix) <- species_list
		colnames(est_matrix) <- species_list
		pval_matrix <- as.data.frame(matrix(1,length(species_list),length(species_list)))
		rownames(pval_matrix) <- species_list
		colnames(pval_matrix) <- species_list
		consistency_matrix <- as.data.frame(matrix(1,length(species_list),length(feature_list)))
		rownames(consistency_matrix) <- species_list
		colnames(consistency_matrix) <- feature_list
        for(i in 1:length(species_list))
        {
                species1 <- species_list[i]
                #print(species1)
                for(j in 1:length(species_list))
                {
                        species2 <- species_list[j]
                        print(paste0(species1,",",species2))
                        if(species1 != species2)
                        {
								tryCatch(               
								expr = {                     
										temp_rem <- compute_meta_corr(data,species1,species2,group_name,study_list)
										#print(temp_rem)
										consistency_matrix[species1,species2] <- length(which(sign(temp_rem$df_studies$di)==sign(as.numeric(temp_rem$model$beta))))/length(temp_rem$df_studies$di)
										est_matrix[species1,species2] <- as.numeric(temp_rem$beta)
										pval_matrix[species1,species2] <- temp_rem$pval
		
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
								)
                               
                        }
                }
        }
        qval_matrix <- apply(pval_matrix,2,function(x)(p.adjust(x,method="fdr")))
		dir_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
		rownames(dir_matrix) <- species_list
		colnames(dir_matrix) <- species_list
		for(i in 1:length(species_list))
        {
          for(j in 1:length(species_list))
          {
			dir_matrix[i,j] <- ifelse(qval_matrix[i,j]<=0.001,sign(est_matrix[i,j]),0)
		  }
		}
		return_list <- list("est"=est_matrix,"pval"=pval_matrix,"consistency"=consistency_matrix,"qval"=qval_matrix,"dir"=dir_matrix)
		return(return_list)
}

rem_network2 <- function(data,species_list,feature_list,group_name,study_list)
{
        species_data <- data[,species_list]
        species_data$group <- data[,group_name]
        est_matrix <- as.data.frame(matrix(0,length(species_list),length(feature_list)))
		rownames(est_matrix) <- species_list
		colnames(est_matrix) <- feature_list
		pval_matrix <- as.data.frame(matrix(1,length(species_list),length(feature_list)))
		rownames(pval_matrix) <- species_list
		colnames(pval_matrix) <- feature_list
		consistency_matrix <- as.data.frame(matrix(1,length(species_list),length(feature_list)))
		rownames(consistency_matrix) <- species_list
		colnames(consistency_matrix) <- feature_list
        for(i in 1:length(species_list))
        {
                species1 <- species_list[i]
                #print(species1)
                for(j in 1:length(feature_list))
                {
                        species2 <- feature_list[j]
                        print(paste0(species1,",",species2))
                        if(species1 != species2)
                        {
								tryCatch(               
								expr = {                     
										temp_rem <- compute_meta_lm(data,species1,species2,group_name,study_list)
										#print(temp_rem)
										consistency_matrix[species1,species2] <- length(which(sign(temp_rem$df_studies$di)==sign(as.numeric(temp_rem$model$beta))))/length(temp_rem$df_studies$di)
										est_matrix[species1,species2] <- as.numeric(temp_rem$model$beta)
										pval_matrix[species1,species2] <- temp_rem$model$pval
		
									},
									error = function(e){         
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
								)
                               
                        }
                }
        }
        qval_matrix <- apply(pval_matrix,2,function(x)(p.adjust(x,method="fdr")))
		dir_matrix <- as.data.frame(matrix(0,length(species_list),length(feature_list)))
		rownames(dir_matrix) <- species_list
		colnames(dir_matrix) <- feature_list
		for(i in 1:length(species_list))
        {
          for(j in 1:length(feature_list))
          {
			dir_matrix[i,j] <- ifelse((qval_matrix[i,j]<=0.01)&&(consistency_matrix[i,j]>=0.70),2*sign(est_matrix[i,j]),ifelse((pval_matrix[i,j]<=0.05)&&(consistency_matrix[i,j]>=0.70),sign(est_matrix[i,j]),0))
		  }
		}
		return_list <- list("est"=est_matrix,"pval"=pval_matrix,"consistency"=consistency_matrix,"qval"=qval_matrix,"dir"=dir_matrix)
		return(return_list)
}

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
        detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
        rownames(detection_matrix) <- var1_list
        colnames(detection_matrix) <- grouping_list
        for(i in 1:length(var1_list))
        {
                var1 <- var1_list[i]
                for(j in 1:length(grouping_list))
                {
                        group <- grouping_list[j]
                        detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]!=0))/length(data[data[,grouping_variable]==group,var1])
                }
        }
        return(detection_matrix)
}

wilcox_batch = function(x,y)
{
	p_array <- NULL;
	type_array <- NULL;
	mean1_array <- NULL;
	mean2_array <- NULL;
	x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
	y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
	z <- intersect(rownames(x),rownames(y));
	for(i in 1:length(z))
	{
		p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
		type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
		mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
		mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
		i <- i + 1;
	}
	out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array,method="fdr"),mean1_array,mean2_array));
	rownames(out) <- z;
	out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
	return(t(out));
}

wilcox_batch_paired = function(x,y)
{
	p_array <- NULL;
	type_array <- NULL;
	mean1_array <- NULL;
	mean2_array <- NULL;
	x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
	y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
	z <- intersect(rownames(x),rownames(y));
	for(i in 1:length(z))
	{
		p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]),paired=TRUE)$p.value;
		type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
		mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
		mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
		i <- i + 1;
	}
	out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array),mean1_array,mean2_array));
	rownames(out) <- z;
	out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
	return(t(out));
}

wilcox_batch1 = function(x,y)
{
	p_array <- NULL;
	type_array <- NULL;
	mean1_array <- NULL;
	mean2_array <- NULL;
	#x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
	#y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
	z <- intersect(rownames(x),rownames(y));
	for(i in 1:length(z))
	{
		p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
		type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
		mean1_array[i] <- median(as.numeric(x[z[i],]),na.rm=TRUE);
		mean2_array[i] <- median(as.numeric(y[z[i],]),na.rm=TRUE);
		i <- i + 1;
	}
	out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array),mean1_array,mean2_array));
	rownames(out) <- z;
	out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
	return(t(out));
}



iterative_compare = function(data,window,disease,control,iter,fractiondisease,fractioncontrol)
{
	set.seed(100);
	#featureProfile <- as.data.frame(matrix(NA,200,ncol(data)));
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	totalDiffFeatures <- NULL;
	totalLostFeatures <- NULL;
	totalGainedFeatures <- NULL;
	lengthFeatures <- NULL;
	featureAssociation <- NULL;
	nominalPValues <- NULL;
	fdrPValues <- NULL;
	medianGroup1 <- NULL;
	medianGroup2 <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	logChangeGainLoss <- NULL;
	for(i in 1:iter)
	{
		tempDiseaseSamples <- sample(intersect(window,disease),as.integer(fractiondisease*length(intersect(window,disease))),replace=TRUE);
		tempControlSamples <- sample(intersect(window,control),as.integer(fractioncontrol*length(intersect(window,control))),replace=TRUE);	
		tempDisease <- data[tempDiseaseSamples,];
		tempControl <- data[tempControlSamples,];
		numberDiseasedSamples[i] <- length(tempDiseaseSamples);
		numberControlSamples[i] <- length(tempControlSamples);
		compareFeatures <- wilcox_batch(t(tempDisease),t(tempControl));
		if(i == 1)
		{
			featureAssociation <- as.data.frame(compareFeatures[compareFeatures[,3] < 0.15,2]);
			nominalPValues <- as.data.frame(compareFeatures[,1]);
			fdrPValues <- as.data.frame(compareFeatures[,3]);
			medianGroup1 <- as.data.frame(compareFeatures[,4]);
			medianGroup2 <- as.data.frame(compareFeatures[,5]);
			
		}
		else
		{
			if(nrow(featureAssociation) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[compareFeatures[,3] < 0.15,2]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(featureAssociation,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(featureAssociation,temp1,all=TRUE,by="row.names")[,1];
					featureAssociation <- apply(temp,2,function(x)(ifelse(is.na(x),0,x)));
				}
			}
			else
			{
				featureAssociation <- as.data.frame(compareFeatures[compareFeatures[,3] < 0.15,2]);
			}

			if(nrow(nominalPValues) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,1]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(nominalPValues,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(nominalPValues,temp1,all=TRUE,by="row.names")[,1];
					nominalPValues <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				nominalPValues <- as.data.frame(compareFeatures[,1]);
			}

			if(nrow(fdrPValues) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[compareFeatures[,3] < 1,3]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(fdrPValues,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(fdrPValues,temp1,all=TRUE,by="row.names")[,1];
					fdrPValues <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				fdrPValues <- as.data.frame(compareFeatures[,3]);
			}

			if(nrow(medianGroup1) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,4]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(medianGroup1,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(medianGroup1,temp1,all=TRUE,by="row.names")[,1];
					medianGroup1 <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				medianGroup1 <- as.data.frame(compareFeatures[,4]);
			}

			if(nrow(medianGroup2) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,5]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(medianGroup2,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(medianGroup2,temp1,all=TRUE,by="row.names")[,1];
					medianGroup2 <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				medianGroup2 <- as.data.frame(compareFeatures[,5]);
			}
		}
		totalDiffFeatures[i] <- ifelse(is.null(nrow(compareFeatures[(compareFeatures[,3] < 0.15),])),0,nrow(compareFeatures[(compareFeatures[,3] < 0.15),]));
		totalLostFeatures[i] <- ifelse(is.null(nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (-1)),])),0,nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (-1)),]));
		totalGainedFeatures[i] <- ifelse(is.null(nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (1)),])),0,nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (1)),]));
		lengthFeatures[i] <- length(as.vector(ifelse(compareFeatures[,3] < 0.2,compareFeatures[,2],0)));
		logChangeGainLoss[i] <- log((totalGainedFeatures[i] + 0.0001)/(totalLostFeatures[i] + 0.0001))
		i <- i + 1;
	}
	#colnames(featureProfile) <- colnames(data);
	returnList <- list("totalDiffFeatures" = totalDiffFeatures,"totalLostFeatures" = totalLostFeatures,"totalGainedFeatures" = totalGainedFeatures,"featureAssociation"=featureAssociation,"nominalP"=nominalPValues,"fdrP"=fdrPValues,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"logChangeGainLoss"=logChangeGainLoss,"medianDisease"=medianGroup1,"medianControl"=medianGroup2);
	return(returnList);
}

iterative_compare_fixed_size = function(data,window,disease,control,iter,diseaseSize,controlSize)
{
	set.seed(100);
	#featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	totalDiffFeatures <- NULL;
	totalLostFeatures <- NULL;
	totalGainedFeatures <- NULL;
	lengthFeatures <- NULL;
	tempDiseaseSamples <- NULL;
	tempControlSamples <- NULL;
	nominalPValues <- NULL;
	fdrPValues <- NULL;
	medianGroup1 <- NULL;
	medianGroup2 <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	logChangeGainLoss <- NULL;
	for(i in 1:iter)
	{
		if(diseaseSize < length(intersect(window,disease)))
		{
			tempDiseaseSamples <- sample(intersect(window,disease),diseaseSize,replace=FALSE);
		}
		else
		{
			tempDiseaseSamples <- sample(intersect(window,disease),diseaseSize,replace=TRUE);
		}
		if(controlSize >= length(intersect(window,disease)))
		{
			tempControlSamples <- sample(intersect(window,control),controlSize,replace=FALSE);	
		}
		else
		{
			tempControlSamples <- sample(intersect(window,control),controlSize,replace=TRUE);	
		}
		
		tempDisease <- data[tempDiseaseSamples,];
		tempControl <- data[tempControlSamples,];
		numberDiseasedSamples[i] <- length(tempDiseaseSamples);
		numberControlSamples[i] <- length(tempControlSamples);
		compareFeatures <- wilcox_batch(t(tempDisease),t(tempControl));
		if(i == 1)
		{
			nominalPValues <- as.data.frame(compareFeatures[,1]);
			fdrPValues <- as.data.frame(compareFeatures[,3]);
			medianGroup1 <- as.data.frame(compareFeatures[,4]);
			medianGroup2 <- as.data.frame(compareFeatures[,5]);
			
		}
		else
		{
			
			if(nrow(nominalPValues) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,1]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(nominalPValues,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(nominalPValues,temp1,all=TRUE,by="row.names")[,1];
					nominalPValues <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				nominalPValues <- as.data.frame(compareFeatures[,1]);
			}

			if(nrow(fdrPValues) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,3]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(fdrPValues,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(fdrPValues,temp1,all=TRUE,by="row.names")[,1];
					fdrPValues <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}

			}
			else
			{
				fdrPValues <- as.data.frame(compareFeatures[,3]);
			}

			if(nrow(medianGroup1) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,4]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(medianGroup1,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(medianGroup1,temp1,all=TRUE,by="row.names")[,1];
					medianGroup1 <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				medianGroup1 <- as.data.frame(compareFeatures[,4]);
			}

			if(nrow(medianGroup2) > 0)
			{
				temp1 <- as.data.frame(compareFeatures[,5]);
				if(nrow(temp1) > 0)
				{
					temp <- merge(medianGroup2,temp1,all=TRUE,by="row.names")[,-1];
					rownames(temp) <- merge(medianGroup2,temp1,all=TRUE,by="row.names")[,1];
					medianGroup2 <- apply(temp,2,function(x)(ifelse(is.na(x),1,x)));
				}
			}
			else
			{
				medianGroup2 <- as.data.frame(compareFeatures[,5]);
			}
		}
		totalDiffFeatures[i] <- ifelse(is.null(nrow(compareFeatures[(compareFeatures[,3] < 0.15),])),0,nrow(compareFeatures[(compareFeatures[,3] < 0.15),]));
		totalLostFeatures[i] <- ifelse(is.null(nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (-1)),])),0,nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (-1)),]));
		totalGainedFeatures[i] <- ifelse(is.null(nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (1)),])),0,nrow(compareFeatures[(compareFeatures[,3] < 0.15)&(compareFeatures[,2] == (1)),]));
		lengthFeatures[i] <- length(as.vector(ifelse(compareFeatures[,3] < 0.2,compareFeatures[,2],0)));
		logChangeGainLoss[i] <- log((totalGainedFeatures[i] + 0.0001)/(totalLostFeatures[i] + 0.0001))
		i <- i + 1;
	}
	#colnames(featureProfile) <- colnames(data);
	returnList <- list("totalDiffFeatures" = totalDiffFeatures,"totalLostFeatures" = totalLostFeatures,"totalGainedFeatures" = totalGainedFeatures,"nominalP"=nominalPValues,"fdrP"=fdrPValues,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"logChangeGainLoss"=logChangeGainLoss,"medianDisease"=medianGroup1,"medianControl"=medianGroup2,"tempDiseaseSamples"=tempDiseaseSamples,"tempControlSamples"=tempControlSamples);
	return(returnList);
}

mod_feature_profile = function(medianDisease,medianControl,PValue,Pval)
{
		selectSpecies <- union(rownames(medianDisease)[which(apply(medianDisease,1,median)>0)],rownames(medianControl)[which(apply(medianControl,1,median)>0)]);
		t <- ifelse(medianDisease >  medianControl,1,ifelse(medianDisease <  medianControl,-1,0)) * PValue;
		
		totalGainedFeatures <- as.numeric(apply(t,2,function(x)(length(which((x > 0)&(abs(x) <= Pval))))));
		totalLostFeatures <- as.numeric(apply(t,2,function(x)(length(which((x < 0)&(abs(x) <= Pval))))));
		logChangeGainLoss <- log((totalGainedFeatures + 0.0001)/(totalLostFeatures + 0.0001));
		returnList <- list("totalLostFeatures" = totalLostFeatures,"totalGainedFeatures" = totalGainedFeatures,"logChangeGainLoss"=logChangeGainLoss,"featureProfile"=t,"SelectSpecies"=selectSpecies);
		return(returnList);
}

mod_feature_profile2 = function(medianDisease,medianControl,PValue,Pval)
{
		selectSpecies <- union(rownames(medianDisease)[which(apply(medianDisease,1,median)>0)],rownames(medianControl)[which(apply(medianControl,1,median)>0)]);
		t <- ifelse(medianDisease[selectSpecies,] >  medianControl[selectSpecies,],1,ifelse(medianDisease[selectSpecies,] <  medianControl[selectSpecies,],-1,0)) * apply(PValue[selectSpecies,],2,function(x)(p.adjust(x)));
		
		totalGainedFeatures <- as.numeric(apply(t,2,function(x)(length(which((x > 0)&(abs(x) <= Pval))))));
		totalLostFeatures <- as.numeric(apply(t,2,function(x)(length(which((x < 0)&(abs(x) <= Pval))))));
		logChangeGainLoss <- log((totalGainedFeatures + 0.0001)/(totalLostFeatures + 0.0001));
		returnList <- list("totalLostFeatures" = totalLostFeatures,"totalGainedFeatures" = totalGainedFeatures,"logChangeGainLoss"=logChangeGainLoss,"featureProfile"=t,"SelectSpecies"=selectSpecies);
		return(returnList);
}


iterative_rf = function(data,window,disease,control,iter)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(as.integer(length(intersect(window,disease))/2) < length(intersect(window,control)))
		{
			trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		}
		else
		{
			trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=TRUE);
		}
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		testDisease <- setdiff(intersect(window,disease),trainDisease);
		testControl <- setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray);
	return(returnList);
		
}

iterative_rf_variable_track = function(data,window,disease,control,iter,variable_list)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	min_level <- as.data.frame(matrix(NA,iter,length(variable_list)));
	average_occurrence <- as.data.frame(matrix(NA,iter,length(variable_list)));
	root_occurrence <- as.data.frame(matrix(NA,iter,length(variable_list)));
	p_value_array <- as.data.frame(matrix(NA,iter,length(variable_list)));
	control_samples <- matrix(NA,iter,as.integer(length(intersect(window,disease))/2));
	disease_samples <- matrix(NA,iter,as.integer(length(intersect(window,disease))/2));
	interaction_edges <- as.data.frame(matrix(NA,iter,length(variable_list)));
	colnames(min_level) <- variable_list;
	colnames(average_occurrence) <- variable_list;
	colnames(p_value_array) <- variable_list;
	colnames(root_occurrence) <- variable_list;
	colnames(interaction_edges) <- variable_list;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(as.integer(length(intersect(window,disease))/2) < length(intersect(window,control)))
		{
			trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		}
		else
		{
			trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=TRUE);
		}
		control_samples[i,] <- trainControl;
		disease_samples[i,] <- trainDisease;
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain,localImp = TRUE);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		testDisease <- setdiff(intersect(window,disease),trainDisease);
		testControl <- setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		imp_matrix <- measure_importance(rfTempComp);
		temp_matrix <- min_depth_interactions(rfTempComp,variable_list);
		for(j in 1:length(variable_list))
		{
			min_level[i,variable_list[j]] <- imp_matrix[imp_matrix[,1]==variable_list[j],2];
			average_occurrence[i,variable_list[j]] <- (imp_matrix[imp_matrix[,1]==variable_list[j],3]+1)/(imp_matrix[imp_matrix[,1]==variable_list[j],6]+1);
			root_occurrence[i,variable_list[j]] <- imp_matrix[imp_matrix[,1]==variable_list[j],7];
			p_value_array[i,variable_list[j]] <- imp_matrix[imp_matrix[,1]==variable_list[j],8];
			interaction_edges[i,variable_list[j]] <- length(which(temp_matrix[(temp_matrix[,2]==variable_list[j]),4]!=0))
		}
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"MinLevel"=min_level,"AverageOccurrence"=average_occurrence,"isRoot"=root_occurrence,"pValue"=p_value_array,"interactionEdges"=interaction_edges,"control_samples"=control_samples,"disease_samples"=disease_samples);
	return(returnList);
		
}

iterative_rf_variable_track_1 = function(data,window,disease,control,iter,variable_list,Young,Middle,Elderly)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	length_young <- as.integer(length(intersect(Young,intersect(window,disease)))/2)
	length_middle <- as.integer(length(intersect(Middle,intersect(window,disease)))/2)
	length_elderly <- as.integer(length(intersect(Elderly,intersect(window,disease)))/2)
	total_length <- length_young + length_middle + length_elderly;
	min_level <- as.data.frame(matrix(NA,iter,length(variable_list)));
	average_occurrence <- as.data.frame(matrix(NA,iter,length(variable_list)));
	root_occurrence <- as.data.frame(matrix(NA,iter,length(variable_list)));
	p_value_array <- as.data.frame(matrix(NA,iter,length(variable_list)));
	control_samples <- matrix(NA,iter,total_length);
	disease_samples <- matrix(NA,iter,total_length);
	interaction_edges <- as.data.frame(matrix(NA,iter,length(variable_list)));
	colnames(min_level) <- variable_list;
	colnames(average_occurrence) <- variable_list;
	colnames(p_value_array) <- variable_list;
	colnames(root_occurrence) <- variable_list;
	colnames(interaction_edges) <- variable_list;
	
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		
		trainDisease <- c(sample(intersect(Young,intersect(window,disease)),length_young,replace=FALSE),sample(intersect(Middle,intersect(window,disease)),length_middle,replace=FALSE),sample(intersect(Elderly,intersect(window,disease)),length_elderly,replace=FALSE))		
		#trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(length_young < length(intersect(Young,intersect(window,control))))
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=FALSE);
		}
		else
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=TRUE);
		}
		if(length_middle < length(intersect(Middle,intersect(window,control))))
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=FALSE);
		}
		else
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=TRUE);
		}
		if(length_elderly < length(intersect(Elderly,intersect(window,control))))
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=FALSE);
		}
		else
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=TRUE);
		}
		trainControl <- c(trainControlYoung,trainControlMiddle,trainControlElderly)
		control_samples[i,] <- trainControl;
		disease_samples[i,] <- trainDisease;
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain,localImp = TRUE);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		testDisease <- setdiff(intersect(window,disease),trainDisease);
		testControl <- setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		imp_matrix <- measure_importance(rfTempComp);
		temp_matrix <- min_depth_interactions(rfTempComp,mean_sample="top_trees",vars=variable_list);
		for(j in 1:length(variable_list))
		{
			min_level[i,variable_list[j]] <- imp_matrix[imp_matrix[,1]==variable_list[j],2];
			average_occurrence[i,variable_list[j]] <- (imp_matrix[imp_matrix[,1]==variable_list[j],3]+1)/(imp_matrix[imp_matrix[,1]==variable_list[j],6]+1);
			root_occurrence[i,variable_list[j]] <- imp_matrix[imp_matrix[,1]==variable_list[j],7];
			p_value_array[i,variable_list[j]] <- imp_matrix[imp_matrix[,1]==variable_list[j],8];
			interaction_edges[i,variable_list[j]] <- length(which(temp_matrix[(temp_matrix[,2]==variable_list[j]),4]!=0))
		}
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"MinLevel"=min_level,"AverageOccurrence"=average_occurrence,"isRoot"=root_occurrence,"pValue"=p_value_array,"interactionEdges"=interaction_edges,"control_samples"=control_samples,"disease_samples"=disease_samples);
	return(returnList);
		
}

iterative_rf_age_match = function(data,window,disease,control,iter,Young,Middle,Elderly,length_young,length_middle,length_elderly)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	length_young <- as.integer(length(intersect(Young,intersect(window,disease)))/2)
	length_middle <- as.integer(length(intersect(Middle,intersect(window,disease)))/2)
	length_elderly <- as.integer(length(intersect(Elderly,intersect(window,disease)))/2)
	total_length <- length_young + length_middle + length_elderly;
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	control_samples <- matrix(NA,iter,total_length);
	disease_samples <- matrix(NA,iter,total_length);
	print(total_length)
	
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		trainDisease <- c(sample(intersect(Young,intersect(window,disease)),length_young,replace=FALSE),sample(intersect(Middle,intersect(window,disease)),length_middle,replace=FALSE),sample(intersect(Elderly,intersect(window,disease)),length_elderly,replace=FALSE))		
		#trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(length_young < length(intersect(Young,intersect(window,control))))
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=FALSE);
		}
		else
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=TRUE);
		}
		if(length_middle < length(intersect(Middle,intersect(window,control))))
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=FALSE);
		}
		else
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=TRUE);
		}
		if(length_elderly < length(intersect(Elderly,intersect(window,control))))
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=FALSE);
		}
		else
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=TRUE);
		}
		trainControl <- c(trainControlYoung,trainControlMiddle,trainControlElderly)
		control_samples[i,] <- trainControl;
		disease_samples[i,] <- trainDisease;
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		#testDisease <- c(sample(intersect(Young,setdiff(intersect(window,disease),trainDisease)),length_young,replace=FALSE),sample(intersect(Middle,setdiff(intersect(window,disease),trainDisease)),length_middle,replace=FALSE),sample(intersect(Elderly,setdiff(intersect(window,disease),trainDisease)),length_elderly,replace=FALSE))
		#testControl <- c(sample(intersect(Young,setdiff(intersect(window,control),trainControl)),length_young,replace=FALSE),sample(intersect(Middle,setdiff(intersect(window,control),trainControl)),length_middle,replace=FALSE),sample(intersect(Elderly,setdiff(intersect(window,control),trainControl)),length_elderly,replace=FALSE))
		testDisease = setdiff(intersect(window,disease),trainDisease);
		testControl = setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		print(median(AUCArray));
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"control_samples"=control_samples,"disease_samples"=disease_samples);
	return(returnList);
		
}

iterative_rf_age_match_2 = function(data,data1,window,disease,control,iter,Young,Middle,Elderly,length_young,length_middle,length_elderly)
{
	set.seed(100);
	#featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	length_young <- as.integer(length(intersect(Young,intersect(window,disease)))/2)
	length_middle <- as.integer(length(intersect(Middle,intersect(window,disease)))/2)
	length_elderly <- as.integer(length(intersect(Elderly,intersect(window,disease)))/2)
	total_length <- length_young + length_middle + length_elderly;
	AUCArray <- NULL;
	AUCArray1 <- NULL;
	SensitivityArray <- NULL;
	SensitivityArray1 <- NULL;
	SpecificityArray <- NULL;
	SpecificityArray1 <- NULL;
	control_samples <- matrix(NA,iter,total_length);
	disease_samples <- matrix(NA,iter,total_length);
	print(total_length)
	
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		trainDisease <- c(sample(intersect(Young,intersect(window,disease)),length_young,replace=FALSE),sample(intersect(Middle,intersect(window,disease)),length_middle,replace=FALSE),sample(intersect(Elderly,intersect(window,disease)),length_elderly,replace=FALSE))		
		#trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(length_young < length(intersect(Young,intersect(window,control))))
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=FALSE);
		}
		else
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=TRUE);
		}
		if(length_middle < length(intersect(Middle,intersect(window,control))))
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=FALSE);
		}
		else
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=TRUE);
		}
		if(length_elderly < length(intersect(Elderly,intersect(window,control))))
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=FALSE);
		}
		else
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=TRUE);
		}
		trainControl <- c(trainControlYoung,trainControlMiddle,trainControlElderly)
		control_samples[i,] <- trainControl;
		disease_samples[i,] <- trainDisease;
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
		tempTrain1 <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain1) <- make.names(rownames(tempTrain),unique=TRUE);
		
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		rfTempComp1 <- randomForest(as.factor(TrainTags)~.,tempTrain1);
		#featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		#testDisease <- c(sample(intersect(Young,setdiff(intersect(window,disease),trainDisease)),length_young,replace=FALSE),sample(intersect(Middle,setdiff(intersect(window,disease),trainDisease)),length_middle,replace=FALSE),sample(intersect(Elderly,setdiff(intersect(window,disease),trainDisease)),length_elderly,replace=FALSE))
		#testControl <- c(sample(intersect(Young,setdiff(intersect(window,control),trainControl)),length_young,replace=FALSE),sample(intersect(Middle,setdiff(intersect(window,control),trainControl)),length_middle,replace=FALSE),sample(intersect(Elderly,setdiff(intersect(window,control),trainControl)),length_elderly,replace=FALSE))
		testDisease = setdiff(intersect(window,disease),trainDisease);
		testControl = setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		rfTempPredict1 <- predict(rfTempComp1,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		AUCArray1[i] <- auc(TestTags,rfTempPredict1[,2])[1];
		print(AUCArray[i]);
		print(AUCArray1[i]);
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SensitivityArray1[i] <- length(which(predict(rfTempComp1,data[testDisease,])=="Diseased"))/length(predict(rfTempComp1,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		SpecificityArray1[i] <- length(which(predict(rfTempComp1,data[testControl,])=="Control"))/length(predict(rfTempComp1,data[testControl,]));
		SpecificityArray1[i] <- length(which(predict(rfTempComp1,data[testControl,])=="Control"))/length(predict(rfTempComp1,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"AUC1"=AUCArray1,"Sensitivity1"=SensitivityArray1,"Specificity1"=SpecificityArray1,"control_samples"=control_samples,"disease_samples"=disease_samples);
	return(returnList);
		
}

iterative_rf_age_match1 = function(data,window,disease,control,iter,Young,Middle,Elderly,length_young,length_middle,length_elderly)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	total_length <- length_young + length_middle + length_elderly;
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	control_samples <- matrix(NA,iter,total_length);
	disease_samples <- matrix(NA,iter,total_length);
	print(total_length)
	
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		trainDisease <- c(sample(intersect(Young,intersect(window,disease)),length_young,replace=FALSE),sample(intersect(Middle,intersect(window,disease)),length_middle,replace=FALSE),sample(intersect(Elderly,intersect(window,disease)),length_elderly,replace=FALSE))		
		#trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(length_young < length(intersect(Young,intersect(window,control))))
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=FALSE);
		}
		else
		{
			trainControlYoung <- sample(intersect(Young,intersect(window,control)),length_young,replace=TRUE);
		}
		if(length_middle < length(intersect(Middle,intersect(window,control))))
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=FALSE);
		}
		else
		{
			trainControlMiddle <- sample(intersect(Middle,intersect(window,control)),length_middle,replace=TRUE);
		}
		if(length_elderly < length(intersect(Elderly,intersect(window,control))))
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=FALSE);
		}
		else
		{
			trainControlElderly <- sample(intersect(Elderly,intersect(window,control)),length_elderly,replace=TRUE);
		}
		trainControl <- c(trainControlYoung,trainControlMiddle,trainControlElderly)
		control_samples[i,] <- trainControl;
		disease_samples[i,] <- trainDisease;
		
		print(length(trainDisease));
		print(length(trainControl));
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		print(dim(tempTrain));
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		testDisease <- setdiff(intersect(window,disease),trainDisease);
		testControl <- setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		print(AUCArray[i]);
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"control_samples"=control_samples,"disease_samples"=disease_samples);
	return(returnList);
		
}

iterative_rf_two = function(data,window1,window2,disease,control,fold)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,100,ncol(data)));
	AUCArray <- NULL;
	AccuracyArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:50)
	{
		print(i);
		trainDisease <- sample(intersect(window1,disease),as.integer(length(intersect(window1,disease))/fold),replace=FALSE);
		if(as.integer(length(intersect(window1,disease))/fold) <= length(intersect(window1,control)))
		{
			trainControl <- sample(intersect(window1,control),as.integer(length(intersect(window1,disease))/fold),replace=FALSE);
		}
		else
		{
			trainControl <- sample(intersect(window1,control),as.integer(length(intersect(window1,disease))/fold),replace=TRUE);
		}
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		numberDiseasedSamples[i] <- nrow(trainDisease);
		numberControlSamples[i] <- nrow(trainControl);
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		if(length(setdiff(window1,window2)) == 0)
		{
			testDisease <- setdiff(intersect(window1,disease),trainDisease);
			testControl <- setdiff(intersect(window1,control),trainControl);
		}
		else
		{
			testDisease <- sample(intersect(window2,disease),as.integer(length(intersect(window2,disease))),replace=FALSE);
			testControl <- sample(intersect(window2,control),as.integer(length(intersect(window2,control))),replace=FALSE);
		}
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		AccuracyArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	#returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
	returnList <- list("AUC"=AUCArray,"Accuracy"=AccuracyArray,"numberDiseasedSamples"=numberDiseasedSamples,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"numberControlSamples"=numberControlSamples,"featureProfile"=featureProfile);
	return(returnList);
		
}

iterative_rf_three = function(data,window1,window2,disease,control,trainsize,testsize,iter)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,100,ncol(data)));
	AUCArray <- NULL;
	AccuracyArray <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		print(i);
		#trainDisease <- sample(intersect(window1,disease),trainsize,replace=FALSE);
		if(trainsize <= length(intersect(window1,disease)))
		{
			trainDisease <- sample(intersect(window1,disease),trainsize,replace=FALSE);
		}
		else
		{
			trainDisease <- sample(intersect(window1,disease),trainsize,replace=TRUE);
		}
		if(trainsize <= length(intersect(window1,control)))
		{
			trainControl <- sample(intersect(window1,control),trainsize,replace=FALSE);
		}
		else
		{
			trainControl <- sample(intersect(window1,control),trainsize,replace=TRUE);
		}
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		numberDiseasedSamples[i] <- nrow(trainDisease);
		numberControlSamples[i] <- nrow(trainControl);
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		if(length(setdiff(window1,window2)) == 0)
		{
			if(length(setdiff(intersect(window1,disease),trainDisease)) > testsize)
			{
				testDisease <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=FALSE);
			}
			else
			{
				testDisease <- setdiff(intersect(window1,disease),trainDisease);
			}
			
			if(length(setdiff(intersect(window1,control),trainControl)) > testsize)
			{
				testControl <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=FALSE);
			}
			else
			{
				testControl <- setdiff(intersect(window1,control),trainControl);
			}
			
		}
		else
		{
			if(length(intersect(window2,disease)) > testsize)
			{
				testDisease <- sample(intersect(window2,disease),testsize,replace=FALSE);
			}
			else
			{
				testDisease <- intersect(window2,disease);
			}
			
			if(length(intersect(window2,control)) > testsize)
			{
				testControl <- sample(intersect(window2,control),testsize,replace=FALSE);
			}
			else
			{
				testControl <- intersect(window2,control);
			}
			
			
		}
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rownames(tempTest) <- make.names(rownames(tempTest),unique=TRUE);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		print(AUCArray[i]);
		AccuracyArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	#returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
	returnList <- list("AUC"=AUCArray,"Accuracy"=AccuracyArray,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"featureProfile"=featureProfile);
	return(returnList);
		
}

iterative_rf_two_profiles = function(data1,data2,window,iter,disease,control)
{

	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data2)));
	AUCArray1 <- NULL;
	AUCArray2 <- NULL;
	
	for(i in 1:iter)
	{

		trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE)
		trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		#trainDisease <- sample(intersect(window,disease),trainSize,replace=FALSE)








		#trainControl <- sample(intersect(window,control),trainSize,replace=FALSE);
		
		tempTrain1 <- rbind(data1[trainDisease,],data1[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp1 <- randomForest(as.factor(TrainTags)~.,tempTrain1);
		
		tempTrain2 <- rbind(data2[trainDisease,],data2[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp2 <- randomForest(as.factor(TrainTags)~.,tempTrain2);
		
		
		
		featureProfile[i,] <- sapply(colnames(tempTrain2),function(x)(ifelse(x %in% rownames(rfTempComp2$importance),rfTempComp2$importance[x,],0)));

		testDisease1 <- setdiff(intersect(window,disease),trainDisease);
		testControl1 <- setdiff(intersect(window,control),trainControl);

		tempTest1 <- rbind(data1[testDisease1,],data1[testControl1,]);
		tempTest2 <- rbind(data2[testDisease1,],data2[testControl1,]);
		
		
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease1)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl1)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict1 <- predict(rfTempComp1,tempTest1,type="vote",norm.votes=TRUE)
		AUCArray1[i] <- auc(TestTags,rfTempPredict1[,2])[1];
		rfTempPredict2 <- predict(rfTempComp2,tempTest2,type="vote",norm.votes=TRUE)
		AUCArray2[i] <- auc(TestTags,rfTempPredict2[,2])[1];
		
		i <- i + 1;

	}

	colnames(featureProfile) <- colnames(tempTrain2);

	returnList <- list("featureProfile"=featureProfile,"AUC1"=AUCArray1,"AUC2"=AUCArray2);

	return(returnList);

		

}

iterative_rf_two_profiles_uniform_windows = function(data1,data2,window1,window2,window3,iter,disease,control)
{

	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data2)));
	AUCArray1 <- NULL;
	AUCArray2 <- NULL;
	AUCArray1a <- NULL;
	AUCArray2a <- NULL;
	AUCArray3 <- NULL;
	for(i in 1:iter)
	{
		trainDisease <- c(sample(intersect(window1,disease),10,replace=FALSE),sample(intersect(window2,disease),10,replace=FALSE),sample(intersect(window3,disease),10,replace=FALSE));
		trainControl <- c(sample(intersect(window1,control),10,replace=FALSE),sample(intersect(window2,control),10,replace=FALSE),sample(intersect(window3,control),10,replace=FALSE));
		#trainDisease <- c(sample(intersect(window1,disease),as.integer(length(intersect(c(window1,window2,window3),disease))/6),replace=TRUE),sample(intersect(window2,disease),as.integer(length(intersect(c(window1,window2,window3),disease))/6),replace=TRUE),sample(intersect(window3,disease),as.integer(length(intersect(c(window1,window2,window3),disease))/6),replace=TRUE));
		#trainControl <- c(sample(intersect(window1,control),as.integer(length(intersect(c(window1,window2,window3),disease))/6),replace=FALSE),sample(intersect(window2,control),as.integer(length(intersect(c(window1,window2,window3),disease))/6),replace=FALSE),sample(intersect(window3,control),as.integer(length(intersect(c(window1,window2,window3),disease))/6),replace=FALSE));
		#trainDisease <- sample(intersect(window,disease),trainSize,replace=FALSE)
		#trainControl <- sample(intersect(window,control),trainSize,replace=FALSE);
		print(i);
		
		tempTrain1 <- rbind(data1[trainDisease,],data1[trainControl,]);
		rownames(tempTrain1) <- make.names(rownames(tempTrain1),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp1 <- randomForest(as.factor(TrainTags)~.,tempTrain1);
		
		tempTrain2 <- rbind(data2[trainDisease,],data2[trainControl,]);
		rownames(tempTrain2) <- make.names(rownames(tempTrain2),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp2 <- randomForest(as.factor(TrainTags)~.,tempTrain2);
		
		
		permutedRows <- permute(c(trainDisease,trainControl));
		tempTrain1a <- data1[permutedRows,];
		rownames(tempTrain1a) <- make.names(rownames(tempTrain1a),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp1a <- randomForest(as.factor(TrainTags)~.,tempTrain1a);
		
		
		tempTrain2a <- data2[permutedRows,];
		rownames(tempTrain2a) <- make.names(rownames(tempTrain2a),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp2a <- randomForest(as.factor(TrainTags)~.,tempTrain2a);
		
		
		tempTrain3 <- as.data.frame(data2[c(trainDisease,trainControl),"Age"]);
		rownames(tempTrain3) <- make.names(rownames(tempTrain3),unique=TRUE);
		colnames(tempTrain3) <- "Age";
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp3 <- randomForest(as.factor(TrainTags)~.,tempTrain3);
		
		
		featureProfile[i,] <- sapply(colnames(tempTrain2),function(x)(ifelse(x %in% rownames(rfTempComp2$importance),rfTempComp2$importance[x,],0)));

		testDisease <- setdiff(intersect(c(window1,window2,window3),disease),trainDisease);
		testControl <- setdiff(intersect(c(window1,window2,window3),control),trainControl);

		tempTest1 <- rbind(data1[testDisease,],data1[testControl,]);
		tempTest2 <- rbind(data2[testDisease,],data2[testControl,]);
		tempTest3 <- as.data.frame(data2[c(testDisease,testControl),"Age"]);
		colnames(tempTest3) <- "Age";
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict1 <- predict(rfTempComp1,tempTest1,type="vote",norm.votes=TRUE)
		AUCArray1[i] <- auc(TestTags,rfTempPredict1[,2])[1];
		rfTempPredict1a <- predict(rfTempComp1a,tempTest1,type="vote",norm.votes=TRUE)
		AUCArray1a[i] <- auc(TestTags,rfTempPredict1a[,2])[1];
		rfTempPredict2 <- predict(rfTempComp2,tempTest2,type="vote",norm.votes=TRUE)
		AUCArray2[i] <- auc(TestTags,rfTempPredict2[,2])[1];
		rfTempPredict2a <- predict(rfTempComp2a,tempTest2,type="vote",norm.votes=TRUE)
		AUCArray2a[i] <- auc(TestTags,rfTempPredict2a[,2])[1];
		rfTempPredict3 <- predict(rfTempComp3,tempTest3,type="vote",norm.votes=TRUE)
		AUCArray3[i] <- auc(TestTags,rfTempPredict3[,2])[1];
		i <- i + 1;

	}

	colnames(featureProfile) <- colnames(tempTrain2);

	returnList <- list("featureProfile"=featureProfile,"AUC1"=AUCArray1,"AUC1a"=AUCArray1a,"AUC2"=AUCArray2,"AUC2a"=AUCArray2a,"AUC3"=AUCArray3);

	return(returnList);

		

}



iterative_rf_fixed_size = function(data,window1,window2,disease,control,iter,trainsize,testsize)
{
	set.seed(150);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	trainDiseaseSamples <- NULL;
	testDiseaseSamples <- NULL;
	trainControlSamples <- NULL;
	testControlSamples <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		if(length(intersect(window1,disease)) > trainsize)
		{
			trainDisease <- sample(intersect(window1,disease),trainsize,replace=FALSE);
		}
		else
		{
			trainDisease <- sample(intersect(window1,disease),trainsize,replace=TRUE);
		}

		if(length(intersect(window1,control)) > trainsize)
		{
			trainControl <- sample(intersect(window1,control),trainsize,replace=FALSE);
		}
		else
		{
			trainControl <- sample(intersect(window1,control),trainsize,replace=TRUE);
		}
		
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		numberDiseasedSamples[i] <- nrow(trainDisease);
		numberControlSamples[i] <- nrow(trainControl);
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain,importance=TRUE);
		#print("model created");
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,3],0)));
		if(length(setdiff(window1,window2)) == 0)
		{
			
			if(length(setdiff(intersect(window1,disease),trainDisease)) > testsize)
			{
				#print("LOOP1a");
				testDisease <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=FALSE);
				#print(testDisease);
			}
			else
			{	
				#print(intersect(window1,disease));
				testDisease <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=TRUE);
			}
			
			if(length(setdiff(intersect(window1,control),trainControl)) > testsize)
			{
				testControl <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=FALSE);
			}
			else
			{
				testControl <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=TRUE);
			}
			#print("test created");
			
		}
		else
		{
			if(length(intersect(window2,disease)) > testsize)
			{
				testDisease <- sample(intersect(window2,disease),testsize,replace=FALSE);
			}
			else
			{
				testDisease <- sample(intersect(window2,disease),testsize,replace=TRUE);
			}
			
			if(length(intersect(window2,control)) > testsize)
			{
				testControl <- sample(intersect(window2,control),testsize,replace=FALSE);
			}
			else
			{
				testControl <- sample(intersect(window2,control),testsize,replace=TRUE);
			}
			
			
		}
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rownames(tempTest) <- make.names(rownames(tempTest),unique=TRUE);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE);
		print(i);
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		print(median(AUCArray));
		SensitivityArray[i] <- length(which(predict(rfTempComp,tempTest[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTest[1:testsize,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,tempTest[(testsize+1):nrow(tempTest),])=="Control"))/length(predict(rfTempComp,tempTest[(testsize+1):nrow(tempTest),]));
		if(i == 100)
		{
			trainDiseaseSamples <- trainDisease;
			testDiseaseSamples <- testDisease;
			trainControlSamples <- trainControl;
			testControlSamples <- testControl;
		}
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	#returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
	returnList <- list("AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"featureProfile"=featureProfile,"trainDisease100"=trainDiseaseSamples,"trainControl100"=trainControlSamples,"testDisease100"=testDiseaseSamples,"testControl100"=testControlSamples);
	return(returnList);
		
}

iterative_rf_generic_disease= function(data,window1,window2,disease1,disease2,disease3,disease4,disease5,control,iter,subsample,testsize)
{
	set.seed(150);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	trainDiseaseSamples <- NULL;
	testDiseaseSamples <- NULL;
	trainControlSamples <- NULL;
	testControlSamples <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{

		trainDisease <- c(sample(intersect(disease1,window1),subsample,replace=FALSE),sample(intersect(disease2,window1),subsample,replace=FALSE),sample(intersect(disease3,window1),subsample,replace=FALSE),sample(intersect(disease4,window1),subsample,replace=FALSE),sample(intersect(disease5,window1),subsample,replace=FALSE));
		trainControl <- sample(intersect(control,window1),5*subsample,replace=FALSE);
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		numberDiseasedSamples[i] <- nrow(trainDisease);
		numberControlSamples[i] <- nrow(trainControl);
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		#print("model created");
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		if(length(setdiff(window1,window2)) == 0)
		{
			testDisease <- c(sample(setdiff(intersect(window1,disease1),trainDisease),subsample,replace=FALSE),sample(setdiff(intersect(window1,disease2),trainDisease),subsample,replace=FALSE),sample(setdiff(intersect(window1,disease3),trainDisease),subsample,replace=FALSE),sample(setdiff(intersect(window1,disease4),trainDisease),subsample,replace=FALSE),sample(setdiff(intersect(window1,disease5),trainDisease),subsample,replace=FALSE));
			testControl <- sample(setdiff(intersect(window1,control),trainControl),5*subsample,replace=FALSE);
					
			
		}
		else
		{
			
			testDisease <- c(sample(intersect(disease1,window2),subsample,replace=FALSE),sample(intersect(disease2,window2),subsample,replace=FALSE),sample(intersect(disease3,window2),subsample,replace=FALSE),sample(intersect(disease4,window2),subsample,replace=FALSE),sample(intersect(disease5,window2),subsample,replace=FALSE));
			testControl <- sample(intersect(control,window2),5*subsample,replace=FALSE);
			
		}
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rownames(tempTest) <- make.names(rownames(tempTest),unique=TRUE);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE);
		print(i);
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		print(AUCArray[i]);
		SensitivityArray[i] <- length(which(predict(rfTempComp,tempTest[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTest[1:testsize,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,tempTest[(testsize+1):nrow(tempTest),])=="Control"))/length(predict(rfTempComp,tempTest[(testsize+1):nrow(tempTest),]));
		if(i == 100)
		{
			trainDiseaseSamples <- trainDisease;
			testDiseaseSamples <- testDisease;
			trainControlSamples <- trainControl;
			testControlSamples <- testControl;
		}
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	#returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
	returnList <- list("AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"featureProfile"=featureProfile,"trainDisease100"=trainDiseaseSamples,"trainControl100"=trainControlSamples,"testDisease100"=testDiseaseSamples,"testControl100"=testControlSamples);
	return(returnList);
		
}





iterative_rf_low_size = function(data,window1,window2,disease,control)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,100,ncol(data)));
	AUCArray <- NULL;
	AccuracyArray <- NULL;
	numberDiseasedSamples <- NULL;
	numberControlSamples <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:100)
	{
		trainDisease <- sample(intersect(window1,disease),10,replace=FALSE);
		trainControl <- sample(intersect(window1,control),10,replace=FALSE);
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainControl)] <- "Control";
		numberDiseasedSamples[i] <- nrow(trainDisease);
		numberControlSamples[i] <- nrow(trainControl);
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		if(length(setdiff(window1,window2)) == 0)
		{
			testDisease <- setdiff(intersect(window1,disease),trainDisease);
			testControl <- setdiff(intersect(window1,control),trainControl);
		}
		else
		{
			testDisease <- sample(intersect(window2,disease),as.integer(length(intersect(window2,disease))),replace=TRUE);
			testControl <- sample(intersect(window2,control),as.integer(length(intersect(window2,control))),replace=TRUE);
		}
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		AccuracyArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	#returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
	returnList <- list("AUC"=AUCArray,"Accuracy"=AccuracyArray,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"featureProfile"=featureProfile);
	return(returnList);
		
}	

bootstrap_rf = function(data,metadata,metadata_column,iter,size)
{
	error_array <- NULL;
	cor_array <- NULL;
	common_rows <- intersect(rownames(data),rownames(metadata))
	data <- data[common_rows,]
	metadata <- metadata[common_rows,]
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	for(i in 1:iter)
	{
		print(i);
		trainRows <- sample(rownames(data),as.integer(size*nrow(data)),replace=FALSE);
		testRows <- setdiff(rownames(data),trainRows);
		rfTrain <- randomForest(metadata[trainRows,metadata_column]~.,data[trainRows,]);
		featureProfile[i,] <- sapply(colnames(rfTrain),function(x)(ifelse(x %in% rownames(rfTrain$importance),rfTrain$importance[x,],0)));
		error_array[i] <- mean((predict(rfTrain,data[testRows,])-metadata[testRows,metadata_column])**2);
		cor_array[i] <- cor(predict(rfTrain,data[testRows,]),metadata[testRows,metadata_column]);
		print(error_array[i]);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(data);
	returnList <- list("cor"=cor_array,"featureProfile"=featureProfile);
}

bootstrap_rf_cor = function(data,fim,iter,size)
{
	cor_array <- NULL;
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	for(i in 1:iter)
	{
		print(i);
		trainRows <- sample(rownames(data),as.integer(size*nrow(data)),replace=FALSE);
		testRows <- setdiff(rownames(data),trainRows);
		rfTrain <- randomForest(fim[trainRows]~.,data[trainRows,]);
		#featureProfile[i,] <- sapply(colnames(rfTrain),function(x)(ifelse(x %in% rownames(rfTrain$importance),rfTrain$importance[x,],0)));
		cor_array[i] <- cor(predict(rfTrain,data[testRows,]),fim[testRows]);
		print(cor_array[i]);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(data);
	returnList <- list("cor"=cor_array,"featureProfile"=featureProfile);
	return(cor_array);
}

aggregate_profile = function(compareProfile,pval)
{
	temp <- mod_feature_profile(compareProfile$medianDisease,compareProfile$medianControl,compareProfile$nominalP,pval);
	modDigitalProfile <- apply(temp$featureProfile[temp$SelectSpecies,],2,function(x)(ifelse((x < 0)&(abs(x) < pval),-1,ifelse((x > 0)&(abs(x) < pval),1,0))));
	return(modDigitalProfile);
}

combined_boxplot = function(data,features,group1,group2,horizontal)
{
	library(ggplot2)
	library(reshape)
	tempCombinedData <- rbind(data[group1,species],data[group2,species]);
	meltFrame <- melt(tempCombinedData);
	groupTags <- as.data.frame(c(rep("Group1",length(group1)),rep("Group2",length(group2))));
	rownames(groupTags) <- c(group1,group2);
	meltFrame[,4] <- apply(meltFrame,1,function(x){groupTags[x[1],1]});
	colnames(meltFrame) <- c("Samples","Features","Value","Group");
	if(!horizontal)
	{
		p <- ggplot(meltFrame,aes(x=factor(Features),Group,y=Value,fill=factor(Group)))+ geom_boxplot(outlier.size=-1) + coord_cartesian(ylim=c(0,200))+theme_bw();
	}
	else
	{
		p <- ggplot(meltFrame,aes(x=factor(Features),Group,y=Value,fill=factor(Group)))+ geom_boxplot(outlier.size=-1) + coord_cartesian(xlim=c(0,200))+theme_bw()+coord_flip();
	}
	return(p);
}

bootstrap_envfit=function(data,meta,size,iter)
{
	r2 <- matrix(NA,iter,ncol(meta));
	for(i in 1:iter)
	{
		tempRows <- sample(rownames(data),size,replace=FALSE);
		pcoTemp <- dudi.pco(as.dist(1-cor(t(data[tempRows,]),method="spearman")/2),scannf=FALSE,nf=3);
		for(j in 1:ncol(meta))
		{
			if(is.factor(meta[,j]))
			{
				tempEnvfit <- envfit(pcoTemp,as.data.frame(meta[tempRows,j]),na.rm=TRUE,500);
				r2[i,j] <- as.numeric(tempEnvfit$vectors$r[1]);
			}
			else
			{
				tempEnvfit <- envfit(pcoTemp,as.data.frame(meta[tempRows,j]),na.rm=TRUE,500);
				r2[i,j] <- as.numeric(tempEnvfit$vectors$r);
			}
		}
		
		
	}
	colnames(r2) <- colnames(meta);
	return(r2);
}
	
	
do_fisher = function(t1,pThreshold)
{
	#forFisher <- cbind(apply(t1,1,function(x)(length(x[((abs(x) < 0.05)&(x > 0))]))),apply(t1,1,function(x)(length(x[((abs(x) < 0.05)&(x < 0))]))),sum(apply(t1,1,function(x)(length(x[((abs(x) < 0.05)&(x > 0))]))))-apply(t1,1,function(x)(length(x[((abs(x) < 0.05)&(x > 0))]))),sum(apply(t1,1,function(x)(length(x[((abs(x) < 0.05)&(x < 0))]))))-apply(t1,1,function(x)(length(x[((abs(x) < 0.05)&(x < 0))]))));
	forFisher <- cbind(apply(t1,1,function(x)(length(x[((abs(x) < pThreshold)&(x > 0))]))),apply(t1,1,function(x)(length(x[((abs(x) < pThreshold)&(x < 0))]))),sum(apply(t1,1,function(x)(length(x[(x > 0)])))),sum(apply(t1,1,function(x)(length(x[((x < 0))])))));
	fisherOut <- cbind(forFisher,apply(forFisher,1,function(x)(fisher.test(matrix(x,2,2))$p)),p.adjust(apply(forFisher,1,function(x)(fisher.test(matrix(x,2,2))$p))),ifelse(forFisher[,1] > forFisher[,2],1,-1));
	colnames(fisherOut) <- c("SpGain","SpLoss","NonSpGain","NonSpLoss","PVal","FDRPVal","Direction");
	fisherOut <- fisherOut[fisherOut[,6] < 0.1,];
	GainedSpecies <- rownames(fisherOut[fisherOut[,7] == 1,]);
	LostSpecies <- rownames(fisherOut[fisherOut[,7] == -1,]);
	
	returnList <- list("fisherOut"=fisherOut,"GainedSpecies"=GainedSpecies,"LostSpecies"=LostSpecies);
	return(returnList);
}
	
iterative_coinertia = function(data1,data2,iter,size)
{
	coinArray <- NULL;
	for(i in 1:iter)
	{
		tempRows <- sample(intersect(rownames(data1),rownames(data2)),size,replace=FALSE);
		tempdata1 <- data1[tempRows,];
		tempdata2 <- data2[tempRows,];
		distdata1 <- as.dist(1-cor(t(tempdata1),method="spearman")/2);
		distdata2 <- as.dist(1-cor(t(tempdata2),method="spearman")/2);
		pcodata1 <- dudi.pco(distdata1,scannf=FALSE,nf=5);
		pcodata2 <- dudi.pco(distdata2,scannf=FALSE,nf=5);
		coinArray[i] <- coinertia(pcodata1,pcodata2,scan=FALSE,nf=3)$RV;
		i <- i + 1;
	}
	return(coinArray);
}

transform_venn =function(t,VennRows)
{
	
	tempFrame <- cbind(aggregate(t,by=list(apply(t,1,function(x)(paste0(x[1],x[2],x[3],x[4])))),FUN=sum)[,1],as.numeric(apply(aggregate(t,by=list(apply(t,1,function(x)(paste0(x[1],x[2],x[3],x[4])))),FUN=sum),1,function(x)(max(x[2],x[3],x[4],x[5])))));
	rownames(tempFrame) <- tempFrame[,1];
	VennCount <- NULL;
	for(i in 1:length(VennRows))
	{
		VennCount[i] <- as.numeric(ifelse(VennRows[i] %in% rownames(tempFrame),tempFrame[VennRows[i],2],0));
		i <- i + 1;
	}
	return(VennCount);
}

bootstrap_correlation=function(vector1,vector2,sample_list,iter,sub_sample)
{
	corr_array <- NULL;
	for(i in 1:iter)
	{
			tempRows <- sample(sample_list,sub_sample,replace=FALSE);
			corr_array[i] <- cor(vector1[tempRows],vector2[tempRows],method="kendall");
			i <- i + 1;
	}
	return(corr_array);
}

bootstrap_correlation_1=function(data,var1,var2,sub_sample,iter)
{
	corr_array <- NULL;
	for(i in 1:iter)
	{
			tempRows <- sample(rownames(data),sub_sample,replace=FALSE);
			corr_array[i] <- cor(data[tempRows,c(var1,var2)],method="kendall")[1,2];
			i <- i + 1;
	}
	return(corr_array);
}

bootstrap_envfit_factors=function(data,meta,index,size,iter)
{
	set.seed(100);
	r2 <- matrix(NA,iter,2);
	for(i in 1:iter)
	{
		tempRows <- sample(rownames(data),size,replace=FALSE);
		pcoTemp <- dudi.pco(vegdist(data[tempRows,],method="euclidean"),scannf=FALSE,nf=3);
		tempEnvfit <- envfit(pcoTemp,as.data.frame(meta[tempRows,index]),na.rm=TRUE,500);
		#permute_tempEnvfit <- envfit(pcoTemp,as.data.frame(meta[permute(tempRows),index]),na.rm=TRUE,500);
		r2[i,1] <- as.numeric(tempEnvfit$factors$r[1]);
		#r2[i,2] <- as.numeric(permute_tempEnvfit$factors$r[1]);	
		r2[i,2] <- as.numeric(tempEnvfit$factors$pvals);
		#r2[i,4] <- as.numeric(permute_tempEnvfit$factors$pvals);	
		
	}
	colnames(r2) <- c("ActualR2","ActualP");
	return(r2);
}

bootstrap_adonis_factors=function(data,meta,index,size,iter)
{
	set.seed(100);
	r2 <- matrix(NA,iter,3);
	for(i in 1:iter)
	{
		tempRows <- sample(rownames(data),size,replace=FALSE);
		distTemp <- as.dist(1-cor(t(data[tempRows,]),method="spearman")/2);
		t <- adonis(distTemp~meta[tempRows,index],method="euclidean")
		print(t)
		#permute_tempEnvfit <- envfit(pcoTemp,as.data.frame(meta[permute(tempRows),index]),na.rm=TRUE,500);
		r2[i,1] <- as.numeric(t$aov.tab$R2[1]);
		print(r2[i,1]);
		#r2[i,2] <- as.numeric(permute_tempEnvfit$factors$r[1]);	
		r2[i,2] <- as.numeric(t$aov.tab$Pr[1]);
		r2[i,3] <- as.numeric(t$aov.tab$F.Model)[1]
		#r2[i,4] <- as.numeric(permute_tempEnvfit$factors$pvals);	
		
	}
	colnames(r2) <- c("ActualR2","ActualP","FModel");
	return(r2);
}

anova_confounding=function(data,meta_country,meta_age,species_list)
{
	return_matrix <- matrix(NA,length(species_list),1);
	colnames(return_matrix) <- c("Association");
	rownames(return_matrix) <- species_list;
	p_value <- NULL;
	f_value <- NULL;
	for(i in 1:length(species_list))
	{
		print(species_list[i]);
		fit1 <- lm(data[,species_list[i]]~meta_country);
		fit2 <- lm(data[,species_list[i]]~meta_country+meta_age);
		fit3 <- lm(data[,species_list[i]]~meta_age);
		arr1 <- sign(as.numeric(fit2$coefficients[length(fit2$coefficients)]));
		tanova <- anova(fit1,fit2);
		print(tanova);
		arr2 <- tanova$Pr[2];
		return_matrix[i,] <- arr1*p.adjust(arr2);
		i <- i + 1;
	}
	return(return_matrix);
}



compare_row_wise=function(data,frame1,frame2)
{
	pProfile <- as.data.frame(matrix(NA,nrow(frame1),ncol(data)));
	colnames(pProfile) <- colnames(data);
	for(i in 1:nrow(frame1))
	{
		pProfile[i,] <- -log(wilcox_batch1(t(data[frame1[i,],]),t(data[frame2[i,],]))[,1],10);
		i <- i + 1;
	}
	return(pProfile);
}

filter_list_with_region_confounding=function(frame,confounding)
{
	filtered_list <- as.data.frame(matrix(NA,nrow(frame),ncol(frame)));
	rownames(filtered_list) <-  rownames(frame);

	for(i in 1:nrow(frame))
	{
		if(ncol(frame)==2)
		{
			filtered_list[rownames(frame)[i],1] <- ifelse(abs(confounding[i,1])<0.15,frame[i,1],0);
			filtered_list[rownames(frame)[i],2] <- ifelse(abs(confounding[i,2])<0.15,frame[i,2],0);
		}
		if(ncol(frame)==3)
		{
			filtered_list[rownames(frame)[i],1] <- ifelse(abs(confounding[i,1])<0.15,frame[i,1],0);
			filtered_list[rownames(frame)[i],2] <- ifelse(abs(confounding[i,2])<0.15,frame[i,2],0);
			filtered_list[rownames(frame)[i],3] <- ifelse(abs(confounding[i,3])<0.15,frame[i,3],0);
		}
		i <- i + 1;
	}
	return(filtered_list);
}

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

rank_scale1=function(x)
{
	x <- ifelse(is.na(x),NA,rank(x))
	min_x <- min(x[!is.na(x)])
	max_x <- max(x[!is.na(x)])
	y <- ifelse(is.na(x),NA,(x-min_x)/(max_x-min_x))
	return(y);
}

rank_scale2=function(x,range_min,range_max)
{
	x <- rank(x);
	y <- range_min + (range_max - range_min)*(rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	return(y);
}

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

range_scale2=function(x,range_min,range_max)
{
	y <- range_min + (range_max - range_min)*(x-min(x))/(max(x)-min(x));
	return(y);
}

effect_size_calculator=function(x,y)
{
	library(effsize);
	species_list <- intersect(colnames(x),colnames(y));
	effect_size = matrix(NA,length(species_list),5);
	rownames(effect_size) <- species_list;
	for(i in 1:length(species_list))
	{
		print(species_list[i]);
		gg <- effsize::cohen.d(x[,species_list[i]],y[,species_list[i]])
		#print(gg);
		gt <- wilcox.test(x[,species_list[i]],y[,species_list[i]])
		effect_size[i,1] <- as.numeric(gg$estimate[1]);
		effect_size[i,2] <- as.factor(gg$magnitude[1]);
		effect_size[i,3] <- as.numeric(gt$p.value[1]);
		effect_size[i,4] <- mean(x[,species_list[i]]);
		effect_size[i,5] <- mean(y[,species_list[i]]);
	}
	colnames(effect_size) <- c("Estimate","Level","NominalP","MedianGroup1","MedianGroup2");
	effect_size <- apply(effect_size,2,function(x)(ifelse(is.nan(x),0,x)))
	#effect_size <- effect_size[names(which(!is.nan(rowSums(effect_size)))),];
	return(effect_size);
}
	
iterative_effsize=function(data1,data2,iter,size)
{
	print(dim(data1));
	print(dim(data2));
	nspecies <- intersect(colnames(data1),colnames(data2));
	
	effsize_matrix <- as.data.frame(matrix(NA,iter,length(nspecies)));
	pvalue_matrix <- as.data.frame(matrix(NA,iter,length(nspecies)));
	if(length(nspecies) > 1)
	{
		data1 <- data1[,nspecies];
		data2 <- data2[,nspecies];
		for(i in 1:iter)
		{
			tempRow1 <- sample(rownames(data1),size,replace=FALSE);
			tempRow2 <- sample(rownames(data2),size,replace=FALSE);
			gg <- effect_size_calculator(data1[tempRow1,],data2[tempRow2,]);
			#print(nrow(gg));
			effsize_matrix[i,] <- gg[,1];
			pvalue_matrix[i,] <- -log(gg[,3],10);
		}
	}
	else
	{
		data1_rownames <- rownames(data1);
		data2_rownames <- rownames(data2);
		data1 <- as.data.frame(data1[,nspecies]);
		data2 <- as.data.frame(data2[,nspecies]);
		rownames(data1) <- data1_rownames;
		rownames(data2) <- data2_rownames;
		for(i in 1:iter)
		{
			tempRow1 <- rownames(data1)[sample(1:nrow(data1),size,replace=FALSE)];
			tempRow2 <- rownames(data2)[sample(1:nrow(data2),size,replace=FALSE)];
			gg <- effect_size_calculator(data1[tempRow1,],data2[tempRow2,]);
			#print(nrow(gg));
			effsize_matrix[i,] <- gg[,1];
			pvalue_matrix[i,] <- -log(gg[,3],10);
		}
	}
	colnames(effsize_matrix) <- nspecies;
	colnames(pvalue_matrix) <- nspecies;
	pvalue_matrix <- apply(pvalue_matrix,2,function(x)(ifelse(is.finite(x),x,0)));
	returnlist <- list("EffSize"=effsize_matrix,"PValue"=pvalue_matrix);
	return(returnlist);
}

combined_boxplot = function(data,features,group1,group2,horizontal)
{
	library(ggplot2)
	library(reshape)
	
	tempCombinedData <- rbind(data[group1,species],data[group2,species]);
	meltFrame <- melt(tempCombinedData);
	groupTags <- as.data.frame(c(rep("Group1",length(group1)),rep("Group2",length(group2))));
	rownames(groupTags) <- c(group1,group2);
	meltFrame[,4] <- apply(meltFrame,1,function(x){groupTags[x[1],1]});
	colnames(meltFrame) <- c("Samples","Features","Value","Group");
	if(!horizontal)
	{
		p <- ggplot(meltFrame,aes(x=factor(Features),Group,y=Value,fill=factor(Group)))+ geom_boxplot(outlier.size=-1) + coord_cartesian(ylim=c(0,200))+theme_bw();
	}
	else
	{
		p <- ggplot(meltFrame,aes(x=factor(Features),Group,y=Value,fill=factor(Group)))+ geom_boxplot(outlier.size=-1) + coord_cartesian(xlim=c(0,200))+theme_bw()+coord_flip();
	}
	return(p);
}

batch_dunns = function(data,factor)
{
	species_list <- colnames(data);
	length_factor <- length(unique(factor));
	pair_factor <- length_factor*(length_factor-1)/2;
	comparisons <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	comparisons_adjusted <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	z_changes <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	temp_comp <- NULL;
	k_test <- NULL;
	for(i in 1:length(species_list))
	{
		#print(i)
		temp_dunn <- dunn.test(data[,species_list[i]],factor,method="bh");
		comparisons[i,] <- temp_dunn$P.adjusted;
		z_changes[i,] <- temp_dunn$Z;
		k_test[i] <- kruskal.test(data[,species_list[i]]~factor)$p.value
		temp_comp <- temp_dunn$comparisons;
		i <- i + 1;
	}
	names(k_test) <- species_list
	rownames(comparisons) <- species_list;
	colnames(comparisons) <- temp_comp;
	rownames(z_changes) <- species_list;
	colnames(z_changes) <- temp_comp;
	comparisons_adjusted <- apply(comparisons,2,p.adjust);
	final_trends <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	for(i in 1:length(species_list))
	{
		for(j in 1:length(temp_comp))
		{
			final_trends[i,j] <- ifelse(comparisons_adjusted[i,j] <= 0.10,2*sign(z_changes[i,j]),ifelse(comparisons[i,j] <= 0.001,1*sign(z_changes[i,j]),0));
			j <- j + 1;
		}
		i <- i + 1;
	}
	rownames(final_trends) <- species_list;
	colnames(final_trends) <- temp_comp;
	return_list <- list("kruskal"=k_test,"PValue"=comparisons,"CorrectedP"=comparisons_adjusted,"Z"=z_changes,"final_trends"=final_trends);
	return(return_list);
}

tag_age_classification = function(age_association)
{
	class <- NULL;
	for(i in 1:nrow(age_association))
	{
		
		if((age_association[i,1] < age_association[i,2])&(age_association[i,2] < age_association[i,3]))
		{
			if((abs(age_association[i,5]) > 1.3)&(abs(age_association[i,6]) > 1.3))
			{
				class[i] <- 1;
			}
			else if((abs(age_association[i,4]) > 1.3)&(abs(age_association[i,6]) > 1.3))
			{
				class[i] <- 1;
			}
			else if(abs(age_association[i,6]) > 1.3)
			{
				class[i] <- 2;
			}
			else if(abs(age_association[i,5]) > 1.3)
			{
				class[i] <- 2;
			}
			else
			{
				class[i] <- 3;
			}
		}
		else if((age_association[i,1] > age_association[i,2])&&(age_association[i,2] > age_association[i,3]))
		{
			if((abs(age_association[i,5]) > 1.3)&&(abs(age_association[i,6]) > 1.3))
			{
				class[i] <- 5;
			}
			else if((abs(age_association[i,4]) > 1.3)&(abs(age_association[i,6]) > 1.3))
			{
				class[i] <- 5;
			}
			else if(abs(age_association[i,6]) > 1.3)
			{
				class[i] <- 4;
			}
			else if(abs(age_association[i,5]) > 1.3)
			{
				class[i] <- 2;
			}
			else
			{
				class[i] <- 3;
			}
		}
		else
		{
			class[i] <- 3;
		}
		
	}
	return(class);
}

merge_three_vectors = function(array1,array2,array3)
{
	frame1 <- as.data.frame(matrix(NA,length(array1),1));
	rownames(frame1) <- array1;
	frame1[,1] <- 1;

	frame2 <- as.data.frame(matrix(NA,length(array2),1));
	rownames(frame2) <- array2;
	frame2[,1] <- 1;

	frame3 <- as.data.frame(matrix(NA,length(array3),1));
	rownames(frame3) <- array3;
	frame3[,1] <- 1;

	temp0 <- merge(frame1,frame2,by="row.names",all=TRUE)[,-1]
	rownames(temp0) <- merge(frame1,frame2,by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,2,function(x)(ifelse(is.na(x),0,x)));

	temp1 <- merge(temp0,frame3,by="row.names",all=TRUE)[,-1]
	rownames(temp1) <- merge(temp0,frame3,by="row.names",all=TRUE)[,1]
	temp1 <- apply(temp1,2,function(x)(ifelse(is.na(x),0,x)));
	
	return(temp1);
}

convert_to_weights = function(feature_weights,pcoData,fac1,fac2)
{
	refined_mat <- matrix(NA,nrow(feature_weights),ncol(pcoData))
	rownames(refined_mat) <- rownames(feature_weights)
	for(i in 1:ncol(feature_weights))
	{
		print(i)
		refined_mat[,i] <- ifelse(feature_weights[,i] >= 0, abs(feature_weights[,i])*fac1*max(pcoData[,i]),abs(feature_weights[,i])*fac2*min(pcoData[,i]))
		i <- i + 1;
	}
	return(refined_mat)
}

lm_batch = function(vector,data,cols)
{
			mat <- matrix(NA,length(cols),2)
			rownames(mat) <- cols
			for(i in 1:length(cols))
			{
				t <-summary(lm(rank_scale(vector)~rank_scale(data[,cols[i]])))
				mat[i,1] <- t$coefficients[2,3]
				mat[i,2] <- t$coefficients[2,4]
				i <- i + 1
			}
			return(mat)
}

lm_batch = function(vector,data,cols)
{
			mat <- matrix(NA,length(cols),2)
			rownames(mat) <- cols
			for(i in 1:length(cols))
			{
				t <-summary(lm(rank_scale(vector)~rank_scale(data[,cols[i]])))
				mat[i,1] <- t$coefficients[2,3]
				mat[i,2] <- t$coefficients[2,4]
				i <- i + 1
			}
			return(mat)
}

iterative_lm = function(data,col1,col2,sample,iter)
{
			e <- NULL
			for(i in 1:iter)
			{
				tempRows <- sample(nrow(data),sample,replace=FALSE)
				t <-summary(lm(as.formula(paste0(col1,"~",col2)),data=data[tempRows,]))
				e[i] <- t$coefficients[nrow(t$coefficients),1]
			}
			return(e)
}

iterative_lme = function(data,disease,control,features,region,iter,size)
{
	AICTab <- as.data.frame(matrix(NA,length(features),iter))
	TValueTab <- as.data.frame(matrix(NA,length(features),iter))
	rownames(AICTab) <- features
	rownames(TValueTab) <- features
	for(j in 1:iter)
	{
		print(features[j])
		trainDisease <- sample(disease,20,replace=FALSE)
		trainControl <- sample(control,20,replace=FALSE)
		#print(trainDisease)
		#print(trainControl)
		for(i in 1:length(features))
		{
			#print (length(droplevels(combined_species_profile_with_age_country_final[c(trainDisease,trainControl),region])))
			t <- t <- summary(glmer(as.numeric(c(rep(1,length(trainDisease)),rep(0,length(trainControl))))~data[c(trainDisease,trainControl),features[i]]+(1|droplevels(data[c(trainDisease,trainControl),region])),family="binomial"))
			print(j)
			#print(summary(t))
			#AICTab[i,j] <- as.numeric(t$AICtab[1])
			#AICTab[i,j] <- t$AICtab
			#TValueTab[i,j] <- t$coefficients[2,3]
			i <- i + 1
		}
		j <- j + 1
	}
	returnList = list("AIC"=AICTab,tValue=TValueTab)
	return(returnList)
}
	
# Function to calculate AICc for PERMANOVA. Requires input from adonis or adonis2 {vegan}


AICc.PERMANOVA <- function(adonis.model) {
    
    # check to see if object is an adonis model...
    
    if (!(adonis.model$aov.tab[1,1] >= 1))
        stop("object not output of adonis {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model
    # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
    
    RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
    MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
    
    k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
    
    nn <- nrow(adonis.model$model.matrix)
    
    # AIC : 2*k + n*ln(RSS)
    # AICc: AIC + [2k(k+1)]/(n-k-1)

    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
    
    # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
    # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
    
    AIC <- 2*k + nn*log(RSS)
    AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
    AIC.MSE <- 2*k + nn * log(MSE)
    AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
    AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
    
    output <- list("AIC" = AIC, "AIC.g" = AIC.g, "AICc" = AICc,
                   "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                   "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k)
    
    return(output)   

}

AICc.PERMANOVA2 <- function(adonis2.model) {
    
    # check to see if object is an adonis2 model...
    
    if (is.na(adonis2.model$SumOfSqs[1]))
        stop("object not output of adonis2 {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model Calculating AICc
    # using residual sum of squares (RSS or SSE) since I don't think that adonis
    # returns something I can use as a likelihood function... maximum likelihood
    # and MSE estimates are the same when distribution is gaussian See e.g.
    # https://www.jessicayung.com/mse-as-maximum-likelihood/;
    # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
    # So using RSS or MSE estimates is fine as long as the residuals are
    # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
    # conditional likelihoods then AIC is not valid. However, comparing models
    # with different error distributions is ok (above link).
    
    
    RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
    MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
    
    k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    
    # AIC : 2*k + n*ln(RSS/n)
    # AICc: AIC + [2k(k+1)]/(n-k-1)
    
    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
    # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
    
    # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
    # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
    
    AIC <- 2*k + nn*log(RSS/nn)
    AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
    AIC.MSE <- 2*k + nn * log(MSE)
    AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
    AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
    
    output <- list("AIC" = AIC, "AICc" = AICc, "AIC.g" = AIC.g, 
                   "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                   "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k, "N" = nn)
    
    return(output)   
    
}

adjust_by_metadata <- function(data,metadata,meta_field)
{
	out_matrix <- as.data.frame(matrix(NA,nrow(data),ncol(data)))
	rownames(out_matrix) <- rownames(data)
	colnames(out_matrix) <- colnames(data)
	for(i in 1:ncol(data))
	{
		Taxa <- colnames(data)[i]
		temp_df <- data.frame(Taxa=data[,Taxa],Meta1=as.factor(metadata[,meta_field[1]]))
		out_matrix[,i] <- lm(rank_scale(Taxa)~Meta1,data=temp_df)$residuals
	}
	return(out_matrix)
}

core_central_taxa <- function(data)
{
	print("Calculating prevalence")
	prevalence <- colSums(apply(data,2,function(x)(ifelse(x>0,1,0))))/nrow(data)
	print("Computing ccrepe")
	ccrepe_data <- ccrepe(data)
	ccrepe_data$p.values = apply(ccrepe_data$p.values,2,function(x)(ifelse(is.na(x),1,x)))
	print("Computing Centrality Measures")
	betweenness <- rank_scale(betweenness(graph_from_adjacency_matrix(apply(apply(ccrepe_data$p.values,2,function(x)(p.adjust(x,method="fdr"))),2,function(x)(ifelse(x<=0.15,1,0))))))
	degree <- rank_scale(degree(graph_from_adjacency_matrix(apply(apply(ccrepe_data$p.values,2,function(x)(p.adjust(x,method="fdr"))),2,function(x)(ifelse(x<=0.15,1,0))))))
	degree[setdiff(names(prevalence),names(degree))] <- 0
	betweenness[setdiff(names(prevalence),names(betweenness))] <- 0
	df_species_properties <- data.frame(prevalence=prevalence,degree=degree[names(prevalence)],betweenness=betweenness[names(prevalence)])
	core_degree_prevalence <- names(which(apply(df_species_properties[,c(1,2)],1,function(x)(length(x[x>=0.65])))==2))
	core_betweenness_prevalence <- names(which(apply(df_species_properties[,c(1,3)],1,function(x)(length(x[x>=0.65])))==2))
	return_list <- list("species_properties"=df_species_properties,"core_degree_prevalence"=core_degree_prevalence,"core_betweenness_prevalence"=core_betweenness_prevalence)
	return(return_list)
}

compute_prevalent_single_data <- function(data,threshold){
  detection_percentage <- colSums(apply(data,2,function(x)(ifelse(x>0,1,0))))/nrow(data)
  highly_detected <- names(which(detection_percentage>=threshold))
  return(highly_detected)
}

keystoneInfluence <- function(species, inputData) {
  tryCatch({
	set.seed(100)
    colIndex= which(colnames(inputData)==species)
    newData= inputData[,-colIndex]
    newData <- newData[which(rowSums(newData)!= 0), ]
    print(dim(newData))
    newData <- newData / rowSums(newData)
	print(unique(rowSums(newData)))
    cat("Creating distance matrix\n")
    distanceMatrix <- vegdist(newData, method = "bray")
    print("distance matrix done")
    cat("Creating dudi.pco\n")
    pco <- dudi.pco(distanceMatrix, scannf = FALSE)
    print("pco is created")
    pcoPointsDf <- pco$li
    cat("Generating model\n")
    model <- envfit(pcoPointsDf ~ inputData[rownames(newData), species])
	print(model)
	r_value <- as.numeric(model$vectors$r)
	p_value <- as.numeric(model$vectors$pvals)
	return_list = list("r"=r_value,"p"=p_value)
	return(return_list)
	})
}

compute_meta_lm_double_adjust <- function(data,var1,var2,var3,var4,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		#print(group)
		temp_data <- data[data[,grouping_variable]==group,]
		length_var3 <- length(unique(temp_data[,var3]))
		#print(length_var3)
		length_var4 <- length(unique(temp_data[,var4]))
		#print(length_var4)
		vec1 <- temp_data[,var1]
		vec2 <- temp_data[,var2]
		#print(length(vec1[vec1 > 0]))
		#print(length(vec2[vec2 > 0]))
		if((length(vec1[vec1 > 0]) > 0)&(length(vec2[vec2 > 0]) > 0))
		{
			if((length_var3 == 1) & (length_var4 == 1))
			{
				#print("Loop1")
				f <- as.formula(paste0(var1,"~",var2))
				temp_rlm <- rlm(f,data=temp_data)
				summary_temp_rlm <- summary(temp_rlm)
				#print("Enter")
				#print(summary_temp_rlm)
				temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
				temp_meta[i,3] <- nrow(temp_data)
				temp_meta[i,4] <- 1
				temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
				temp_meta[i,6] <- sign(temp_meta[i,2])
			}
			else if((length_var3 > 1) & (length_var4 == 1))
			{
				#print("Loop2")
				f <- as.formula(paste0(var1,"~",var3,"+",var2))
				temp_rlm <- rlm(f,data=temp_data)
				summary_temp_rlm <- summary(temp_rlm)
				temp_meta[i,2] <- summary_temp_rlm$coefficients[3,3]
				temp_meta[i,3] <- nrow(temp_data)
				temp_meta[i,4] <- 1
				temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
				temp_meta[i,6] <- sign(temp_meta[i,2])
			}
			else if((length_var3 == 1) & (length_var4 > 1))
			{
				#print("Loop3")
				f <- as.formula(paste0(var1,"~",var4,"+",var2))
				temp_rlm <- rlm(f,data=temp_data)
				summary_temp_rlm <- summary(temp_rlm)
				temp_meta[i,2] <- summary_temp_rlm$coefficients[3,3]
				temp_meta[i,3] <- nrow(temp_data)
				temp_meta[i,4] <- 1
				temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
				temp_meta[i,6] <- sign(temp_meta[i,2])
			}
			else if((length_var3 > 1) & (length_var4 > 1))
			{
				#print("Loop4")
				f <- as.formula(paste0(var1,"~",var4,"+",var3,"+",var2))
				temp_rlm <- rlm(f,data=temp_data)
				summary_temp_rlm <- summary(temp_rlm)
				temp_meta[i,2] <- summary_temp_rlm$coefficients[4,3]
				temp_meta[i,3] <- nrow(temp_data)
				temp_meta[i,4] <- 1
				temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
				temp_meta[i,6] <- sign(temp_meta[i,2])
			}
			else
			{
				temp_meta[i,2] <- 0
				temp_meta[i,3] <- nrow(temp_data)
				temp_meta[i,4] <- 1
				temp_meta[i,5] <- 1
				temp_meta[i,6] <- 1
			}
		}
		else
		{
			temp_meta[i,2] <- 0
			temp_meta[i,3] <- nrow(temp_data)
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- 1
			temp_meta[i,6] <- 1
		}
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

compute_meta_lm_group_adjust <- function(data,feature_list,metadata_var,confound_var1,grouping_var,grouping_list)
{
	return_out <- as.data.frame(matrix(NA,length(feature_list),9))
	rownames(return_out) <- feature_list
	colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir")
	for(i in 1:length(feature_list))
	{
		species_name <- feature_list[i]
		tryCatch(               
							expr = {                     
									temp_res <- compute_meta_lm_single_adjust(data,metadata_var,species_name,confound_var1,grouping_var,grouping_list)
									return_out[i,"beta"] <- temp_res$model$beta
									return_out[i,"pval"] <- temp_res$model$pval
									return_out[i,"ci.ub"] <- temp_res$model$ci.ub
									return_out[i,"ci.lb"] <- temp_res$model$ci.lb
									return_out[i,"tau2"] <- temp_res$model$tau2
									return_out[i,"QE"] <- temp_res$model$QE
									return_out[i,"QEp"] <- temp_res$model$QEp	
		
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
		
		
		
		
	}
	return_out$qval <- p.adjust(return_out$pval,method="fdr")
	return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta)))
	return(return_out)
}

batch_rlm <- function(data,feature_list,metadata_var)
{
	return_mat <- as.data.frame(matrix(NA,length(feature_list),4))
	rownames(return_mat) <- feature_list
	colnames(return_mat) <- c("beta","p_val","q_val","dir")
	return_mat[,1] <- 0
	return_mat[,2] <- 1
	for(i in 1:length(feature_list))
	{
		feature_name <- feature_list[i]
		
		tryCatch(               
							expr = {                     
									f <- as.formula(paste0(feature_name,"~",metadata_var))
									temp_rlm <- rlm(f,data=data)
									summary_temp_rlm <- summary(temp_rlm)
									return_mat[i,1] <- summary_temp_rlm$coefficients[2,3]
									return_mat[i,2] <- f.robftest(temp_rlm,var=metadata_var)$p.value
		
									},
									error = function(e){ 
										return_mat[i,1] <- 0
										return_mat[i,2] <- 1
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
						)
	}
	return_mat[,1] <- ifelse(is.nan(return_mat[,1]),0,return_mat[,1])
	return_mat[,3] <- p.adjust(return_mat[,2],method="fdr")
	return_mat[,4] <- ifelse(return_mat[3] <= 0.15,3*sign(return_mat[,1]),ifelse(return_mat[,2] <= 0.05,2*sign(return_mat[,1]),sign(return_mat[,1])))
	return(return_mat)
}

oob_validation <- function(data,metadata,features_to_check,grouping_column,metadata_column,groups_list)
{
	df_Perform <- as.data.frame(matrix(NA,length(groups_list),2))
	rownames(df_Perform) <- groups_list
	colnames(df_Perform) <- c("oob_corr","oob_p")
	
	pairwise_matrix_r <- as.data.frame(matrix(NA,length(groups_list),length(groups_list)))
	rownames(pairwise_matrix_r) <- groups_list
	colnames(pairwise_matrix_r) <- groups_list
	
	pairwise_matrix_p <- as.data.frame(matrix(NA,length(groups_list),length(groups_list)))
	rownames(pairwise_matrix_p) <- groups_list
	colnames(pairwise_matrix_p) <- groups_list
	
	mat_importance <- as.data.frame(matrix(0,length(groups_list),length(features_to_check)))
	rownames(mat_importance) <- groups_list
	colnames(mat_importance) <- features_to_check
	
	group_list_rows <- rownames(data[data[,grouping_column] %in% groups_list,])
	rows_with_metadata <- rownames(metadata[!is.na(metadata[,metadata_column]),])
	df_pred <- data[intersect(rows_with_metadata,group_list_rows),features_to_check]
	
	for(i in 1:length(groups_list))
	{
		print(i)
		group_name <- groups_list[i]
		train_rows <- intersect(rownames(df_pred),rownames(metadata)[metadata[,grouping_column] == group_name])
		print(train_rows)
		df_train_pred <- df_pred[train_rows,]
		df_train_pred <- apply(df_train_pred,2,function(x)(ifelse(is.nan(x),0,x)))
		df_train_pred <- apply(df_train_pred,2,function(x)(ifelse(is.na(x),0,x)))
		df_train_pred <- df_train_pred[rowSums(df_train_pred)>0,colSums(df_train_pred)>0]
		train_rows <- rownames(df_train_pred)
		
		df_train_response <- metadata[train_rows,metadata_column]
		
		rf_oob_temp <- randomForest(df_train_response~.,df_train_pred)
		
		df_oob_corr <- cor.test(rf_oob_temp$predicted,rf_oob_temp$y,method="spearman",use="pairwise.complete")
		
		print(df_oob_corr)
		
		mat_importance[group_name,rownames(rf_oob_temp$importance)] <- rf_oob_temp$importance
		
		for(j in 1:length(groups_list))
		{
			test_study <- groups_list[j]
			if(i == j)
			{
				pairwise_matrix_r[i,j] <- df_test_study_corr$estimate
				pairwise_matrix_p[i,j] <- df_test_study_corr$p.value
			}
			else
			{
				test_study_rows <- intersect(rownames(df_pred),rownames(metadata)[metadata[,grouping_column] == test_study])
				test_study_rows <- intersect(rows_with_metadata,test_study_rows)
				df_pred_test_study <- data[test_study_rows,colnames(df_train_pred)]
				response_test_study <- metadata[test_study_rows,metadata_column]
				df_test_study_corr <- cor.test(predict(rf_oob_temp,df_pred_test_study),response_test_study,method="spearman",use="pairwise.complete")
				pairwise_matrix_r[i,j] <- df_test_study_corr$estimate
				pairwise_matrix_p[i,j] <- df_test_study_corr$p.value
			}
		}
		
		df_Perform[i,1] <- df_oob_corr$estimate
		df_Perform[i,2] <- df_oob_corr$p.value
		
		i <- i + 1
	}
	
	df_Correlation <- list("r"=pairwise_matrix_r,"p"=pairwise_matrix_p)
	
	returnlist <- list("Perform"=df_Perform,"Correlation"=df_Correlation,"FeatureImportance"=mat_importance)
	
	return(returnlist)
}
