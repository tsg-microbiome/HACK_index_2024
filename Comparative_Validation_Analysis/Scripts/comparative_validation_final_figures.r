#### Plots for Comparative Evaluation
## Dt. 07.11.2024

## OS

library(pheatmap)
library(ggplot2)
library(vegan)
library(compositions)
library(pcaPP)
##Load the workspace
#load("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/241030_Comparative_Evaluation.RData")
load("241030_Comparative_Evaluation.RData")

####################### Panel 2 (Part1)
colnames(Directionality_Disease_Validation) <- c("Shannon","Kendall_Uniqueness","Dysbiosis_Score","GMWI","HACK_Top_17","GMWI_hack_top_17")
Directionality_Disease_Validation$GMWI_hack_top_17 <- NULL

### Directionality_Disease_Validation
Directionality_Disease_Validation[Directionality_Disease_Validation == -1 | Directionality_Disease_Validation == 1] <- 0
Directionality_Disease_Validation[Directionality_Disease_Validation == -2] <- -1
Directionality_Disease_Validation[Directionality_Disease_Validation == 2] <- 1

# pdf("C:/Users/ompra/Downloads/Cell_reports_revision/Comparative_validation/Directionality_Disease_Validation.pdf", width = 15, height = 10)
# # Generate the heatmap
# my_palette <- c("cyan3", "white", "palevioletred1")
# 
# pheatmap(t(Directionality_Disease_Validation),
#          color = my_palette,
#          #breaks = breaks,
#          fontsize_row = 13,
#          fontsize_col =13,
#          cellheight = 20,
#          cellwidth = 32,
#          cluster_rows = FALSE,
#          cluster_cols = T,
#          border_color = "black",
#          treeheight_col = 0
# )
# dev.off()
### before plotting arrange the names of columns in the following manner
Directionality_Disease_Validation <- Directionality_Disease_Validation[,c("Shannon","Kendall_Uniqueness","GMWI","Dysbiosis_Score","HACK_Top_17")]

rnames_vector <- c("Luan","Xu","Pang","Wallen","Saleem","Parbie","Nagata","MicroDiab_Denmark","Bajaj","Hernandez","Flemer","Song")                                               
Directionality_Disease_Validation <- Directionality_Disease_Validation[rnames_vector,]

#pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel2/Directionality_Disease_Validation_panel2.pdf", width = 15, height = 10)
# Generate the heatmap
#my_palette <- c("cornflowerblue", "white", "darkgoldenrod1")

#pheatmap(t(Directionality_Disease_Validation),
#         color = my_palette,
#         #breaks = breaks,
#         fontsize_row = 13,
#         fontsize_col =13,
#         cellheight = 20,
#         cellwidth = 32,
#         cluster_rows = FALSE,
#         cluster_cols = F,
#         border_color = "black",
#         treeheight_col = 0
#)
#dev.off()




####################### Panel 2 Part 2

library(readxl)

rownames(df_Confounding) <- colnames(Directionality_Disease_Validation)

#pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel2/df_Confounding_Panel2_Part2.pdf", width = 15, height = 10)
# Generate the heatmap
#my_palette <- c("white", "#008080")

#pheatmap(df_Confounding,
#         color = my_palette,
#         #breaks = breaks,
#         fontsize_row = 13,
#         fontsize_col =13,
#         cellheight = 20,
#         cellwidth = 32,
#         cluster_rows = F,
#         cluster_cols = T,
#         border_color = "black",
#         treeheight_col = 0
#)
#dev.off()


####################### Panel 3, pre-required dfs
#### bray_followup Function 
bray_followup <- function(data,follow_up)
{
  follow_up_dist <- as.data.frame(matrix(NA,nrow(data),1))
  rownames(follow_up_dist) <- rownames(follow_up)
  for(i in 1:nrow(follow_up))
  {
    T0 <- follow_up[i,1]
    T1 <- follow_up[i,2]
    if((!is.na(T1))&(T1 %in% rownames(data))&(T0 %in% rownames(data)))
    {
      follow_up_dist[i,1] <- as.numeric(vegdist(data[c(T1,T0),],method="bray"))
    }
  }
  return(follow_up_dist)
}

#### aitchison_followup Function
aitchison_followup <- function(data,follow_up)
{
  data_clr <- as.data.frame(as.matrix(clr(data+0.00001)))
  data_clr <- as.data.frame(t(apply(data_clr,1,function(x)(x-min(x)))))
  follow_up_dist <- as.data.frame(matrix(NA,nrow(data),1))
  rownames(follow_up_dist) <- rownames(follow_up)
  for(i in 1:nrow(follow_up))
  {
    T0 <- follow_up[i,1]
    T1 <- follow_up[i,2]
    if((!is.na(T1))&(T1 %in% rownames(data))&(T0 %in% rownames(data)))
    {
      follow_up_dist[i,1] <- as.numeric(vegdist(data_clr[c(T1,T0),],method="euclidean"))
    }
  }
  return(follow_up_dist)
}



#### Pang 
Pang_et_al_species <- Pang_et_al_species/rowSums(Pang_et_al_species)

temp_Pang_Followup_Metadata <- Pang_Followup_Metadata[intersect(rownames(Pang_et_al_species),rownames(Pang_Followup_Metadata)),]

T0_rows <- rownames(temp_Pang_Followup_Metadata)[(temp_Pang_Followup_Metadata[,2] %in% rownames(Pang_et_al_species))&(temp_Pang_Followup_Metadata[,1] %in% rownames(Pang_et_al_species))]
T1_rows <- paste0(T0_rows,"-2")

Pang_bray_follow_up <- bray_followup(Pang_et_al_species[c(T0_rows,T1_rows),],Pang_Followup_Metadata[c(T0_rows,T1_rows),])
Pang_aitchison_follow_up <- aitchison_followup(Pang_et_al_species[c(T0_rows,T1_rows),],Pang_Followup_Metadata[c(T0_rows,T1_rows),])

df_PangBrayDist <- data.frame("BrayFollowUpDist"=Pang_bray_follow_up[,1],"Shannon"=diversity(Pang_et_al_species[rownames(Pang_bray_follow_up),]),"KendallUniqueness"=kendall_uniqueness(Pang_et_al_species[rownames(Pang_bray_follow_up),]),"DysbiosisScore"=dysbiosis_score(Pang_et_al_species[rownames(Pang_bray_follow_up),],rownames(Pang_bray_follow_up)),"GMWI"=GMHI(Pang_et_al_species[rownames(Pang_bray_follow_up),]),"hack_top_17"=hack_top_17(Pang_et_al_species[rownames(Pang_bray_follow_up),]),row.names=rownames(Pang_bray_follow_up))

#### Olsson 16s 
#load("H:/.shortcut-targets-by-id/1-VXvgWg44QVMpipxNRFjdw2Oe4gs2t6O/CoreFinder/AdditionalValidation/Olssen_species.RData")
load("Olssen_species.RData")

Olssen_et_al_species <- Olssen_et_al_species/rowSums(Olssen_et_al_species)

Olssen_16S_bray_followup_dist <- bray_followup(Olssen_et_al_species,Olssen_et_al_FollowUp)

df_OlssenBrayDist <- data.frame("BrayFollowUpDist"=Olssen_16S_bray_followup_dist[,1],"Shannon"=diversity(Olssen_et_al_species[rownames(Olssen_16S_bray_followup_dist),]),"KendallUniqueness"=kendall_uniqueness(Olssen_et_al_species[rownames(Olssen_16S_bray_followup_dist),]),"DysbiosisScore"=dysbiosis_score(Olssen_et_al_species[rownames(Olssen_16S_bray_followup_dist),],rownames(Olssen_16S_bray_followup_dist)),"GMWI"=GMHI(Olssen_et_al_species[rownames(Olssen_16S_bray_followup_dist),]),"hack_top_17"=hack_top_17(Olssen_et_al_species[rownames(Olssen_16S_bray_followup_dist),]),row.names=rownames(Olssen_16S_bray_followup_dist))

df_OlssenBrayDist <- df_OlssenBrayDist[!is.na(df_OlssenBrayDist[,1]),]


#### Olsson WGS
#load("H:/.shortcut-targets-by-id/1-VXvgWg44QVMpipxNRFjdw2Oe4gs2t6O/CoreFinder/AdditionalValidation/Olssen_WGS_Species.RData")
load("Olssen_WGS_Species.RData")

Olssen_et_al_wgs_species <- Olsson_et_al_wgs_species
Olssen_et_al_wgs_species <- Olssen_et_al_wgs_species/rowSums(Olssen_et_al_wgs_species)
Olssen_et_al_wgs_FollowUp <- Olsson_et_al_wgs_follow_up_metadata

Olssen_wgs_bray_followup_dist <- bray_followup(Olssen_et_al_wgs_species,Olssen_et_al_wgs_FollowUp)

df_wgs_OlssenBrayDist <- data.frame("BrayFollowUpDist"=Olssen_wgs_bray_followup_dist[,1],"Shannon"=diversity(Olssen_et_al_wgs_species[rownames(Olssen_wgs_bray_followup_dist),]),"KendallUniqueness"=kendall_uniqueness(Olssen_et_al_wgs_species[rownames(Olssen_wgs_bray_followup_dist),]),"DysbiosisScore"=dysbiosis_score(Olssen_et_al_wgs_species[rownames(Olssen_wgs_bray_followup_dist),],rownames(Olssen_wgs_bray_followup_dist)),"GMWI"=GMHI(Olssen_et_al_wgs_species[rownames(Olssen_wgs_bray_followup_dist),]),"hack_top_17"=hack_top_17(Olssen_et_al_wgs_species[rownames(Olssen_wgs_bray_followup_dist),]),row.names=rownames(Olssen_wgs_bray_followup_dist))

df_wgs_OlssenBrayDist <- df_wgs_OlssenBrayDist[!is.na(df_wgs_OlssenBrayDist[,1]),]


####################### Panel 3
### Stability Box Plots and Correlation Plots

####### df_OlssenBrayDist (16s)
colnames(df_OlssenBrayDist) <- c("BrayFollowUpDist","Shannon","KendallUniqueness","DysbiosisScore","GMWI","HACK_Top_17")
df_OlssenBrayDist <- df_OlssenBrayDist[,c("BrayFollowUpDist","Shannon","KendallUniqueness","GMWI","DysbiosisScore","HACK_Top_17")]

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_16s/HACKTop17_Olssen_16s.pdf", width = 4, height = 5)
boxplot(df_OlssenBrayDist$HACK_Top_17~cut(df_OlssenBrayDist[,1],breaks=c(0,0.3,0.67,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray HACK_Top_17",
        xlab = "breaks",
        ylab = "Values")
dev.off()


pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_16s/GMWI_Olssen_16s.pdf", width = 4, height = 5)
boxplot(df_OlssenBrayDist$GMWI~cut(df_OlssenBrayDist[,1],breaks=c(0,0.3,0.67,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray GMWI",
        xlab = "Breaks",
        ylab = "Values")
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_16s/DysbiosisScore_Olssen_16s.pdf", width = 4, height = 5)
boxplot(df_OlssenBrayDist$DysbiosisScore~cut(df_OlssenBrayDist[,1],breaks=c(0,0.3,0.67,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray DysbiosisScore",
        xlab = "Breaks",
        ylab = "Values")
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_16s/KendallUniqueness_Olssen_16s.pdf", width = 4, height = 5)
boxplot(df_OlssenBrayDist$KendallUniqueness~cut(df_OlssenBrayDist[,1],breaks=c(0,0.3,0.67,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray KendallUniqueness",
        xlab = "Breaks",
        ylab = "Values")
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_16s/Shannon_Olssen_16s.pdf", width = 4, height = 5)
boxplot(df_OlssenBrayDist$Shannon~cut(df_OlssenBrayDist[,1],breaks=c(0,0.3,0.67,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray Shannon",
        xlab = "Breaks",
        ylab = "Values")
dev.off()





####### df_wgs_OlssenBrayDist
colnames(df_wgs_OlssenBrayDist) <- c("BrayFollowUpDist","Shannon","KendallUniqueness","DysbiosisScore","GMWI","HACK_Top_17")
df_wgs_OlssenBrayDist <- df_wgs_OlssenBrayDist[,c("BrayFollowUpDist","Shannon","KendallUniqueness","GMWI","DysbiosisScore","HACK_Top_17")]
                                       
pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_WGS/HACK_Top_17_Olssen_WGS.pdf", width = 4, height = 5)
boxplot(df_wgs_OlssenBrayDist$HACK_Top_17~cut(df_wgs_OlssenBrayDist[,1],breaks=c(0,0.3,0.50,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray Hack_Top_18",
        xlab = "breaks",
        ylab = "Values")
dev.off()


pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_WGS/GMWI_Olssen_WGS.pdf", width = 4, height = 5)
boxplot(df_wgs_OlssenBrayDist$GMWI~cut(df_wgs_OlssenBrayDist[,1],breaks=c(0,0.3,0.50,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray GMWI",
        xlab = "Breaks",
        ylab = "Values")
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_WGS/DysbiosisScore_Olssen_WGS.pdf", width = 4, height = 5)
boxplot(df_wgs_OlssenBrayDist$DysbiosisScore~cut(df_wgs_OlssenBrayDist[,1],breaks=c(0,0.3,0.50,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray DysbiosisScore",
        xlab = "Breaks",
        ylab = "Values")
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_WGS/KendallUniqueness_Olssen_WGS.pdf", width = 4, height = 5)
boxplot(df_wgs_OlssenBrayDist$KendallUniqueness~cut(df_wgs_OlssenBrayDist[,1],breaks=c(0,0.3,0.50,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray KendallUniqueness",
        xlab = "Breaks",
        ylab = "Values")
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Olssen_WGS/Shannon_Olssen_WGS.pdf", width = 4, height = 5)
boxplot(df_wgs_OlssenBrayDist$Shannon~cut(df_wgs_OlssenBrayDist[,1],breaks=c(0,0.3,0.50,1),include.lowest=TRUE),outline=FALSE,
        col = c("lightseagreen", "gold2", "indianred3"),  # Fill colors for the boxes
        main = "OlssenBray Shannon",
        xlab = "Breaks",
        ylab = "Values")
dev.off()




####### df_PangBrayDist
colnames(df_PangBrayDist) <- c("BrayFollowUpDist","Shannon","KendallUniqueness","DysbiosisScore","GMWI","HACK_Top_17")
df_PangBrayDist <- df_PangBrayDist[,c("BrayFollowUpDist","Shannon","KendallUniqueness","GMWI","DysbiosisScore","HACK_Top_17")]

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Pang/HACK_Top_17_Pang.pdf", width = 4.5, height = 5)
ggplot(df_PangBrayDist,aes(x=BrayFollowUpDist,y=HACK_Top_17))+geom_point()+geom_smooth(method='lm')+theme_bw()
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Pang/GMWI_Pang.pdf", width = 4.5, height = 5)
ggplot(df_PangBrayDist,aes(x=BrayFollowUpDist,y=GMWI))+geom_point()+geom_smooth(method='lm')+theme_bw()
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Pang/DysbiosisScore_Pang.pdf", width = 4.5, height = 5)
ggplot(df_PangBrayDist,aes(x=BrayFollowUpDist,y=DysbiosisScore))+geom_point()+geom_smooth(method='lm')+theme_bw()
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Pang/KendallUniqueness_Pang.pdf", width = 4.5, height = 5)
ggplot(df_PangBrayDist,aes(x=BrayFollowUpDist,y=KendallUniqueness))+geom_point()+geom_smooth(method='lm')+theme_bw()
dev.off()

pdf("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/panel3/Pang/Shannon_Pang.pdf", width = 4.5, height = 5)
ggplot(df_PangBrayDist,aes(x=BrayFollowUpDist,y=Shannon))+geom_point()+geom_smooth(method='lm')+theme_bw()
dev.off()



all_objects <- ls()
temp_objects <- grep("^temp", all_objects, value = TRUE)
rm(list = temp_objects)


#save.image("C:/Users/ompra/Downloads/Cell_reports_revision/comparative_validation2/revised_comparative_validation/20241107_comparative_validation_final.RData")
save.image("20241107_comparative_validation_final.RData")
