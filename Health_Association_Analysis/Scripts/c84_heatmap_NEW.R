
library("scatterplot3d") 
library(writexl)
library(dplyr)
library(xlsx)
load("SpeciesScores_NEW.RData")

SpeciesScores= SpeciesScores_NEW %>% 
  arrange(-HACKScore)


healthySpecies= rownames(SpeciesScores)[which(SpeciesScores$Health >= 0.70 & SpeciesScores$Influence >= 0.70 & SpeciesScores$Stability >=0.70)]

index= which(SpeciesScores$Health >= 0.70 & SpeciesScores$Influence >= 0.70 & SpeciesScores$Stability >=0.70)

heatmapDf= SpeciesScores[,c(1,2,3)]

heatmapDf[heatmapDf>=0.70]= 1
heatmapDf[heatmapDf<0.70]= 0



# this has to be run on R-console, using this to cluster the species.
library(gplots)
pdf("G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\SpeciesScores_heatmap_modified.pdf", height = 15, width = 15)
heatmapObject= heatmap.2(as.matrix(heatmapDf) , density= "none", trace= "none", Rowv=TRUE, Colv= FALSE)
dev.off()
save(heatmapObject, file= "G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\c84_heamap2_object_new.RData")

load("G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\c84_heamap2_object_new.RData")

heatmapDendogramed_df= as.data.frame(t(heatmapObject$carpet))

# getting all those species that have at least one of the 3 scores >=0.70
filteredHetmap_df= heatmapDendogramed_df[which(rowSums(heatmapDendogramed_df)>0),]

# creating the heatmap
heatmap_union_core_df= filteredHetmap_df

heatmap_union_core_df$species= rownames(filteredHetmap_df)

#################################################################
## 17-10-24 (SG) Continuous heatmap

SpeciesScores_heatmap <- SpeciesScores_NEW[rownames(SpeciesScores_NEW) %in% rownames(heatmap_union_core_df),]
SpeciesScores_heatmap <- SpeciesScores_heatmap[,1:3]


library(pheatmap)

SpeciesScores_heatmap_transposed <- t(SpeciesScores_heatmap)

annotation_matrix <- ifelse(SpeciesScores_heatmap_transposed >= 0.70, "*", "")
light_palette <- colorRampPalette(c("deeppink","white","royalblue"))

colnames(SpeciesScores_heatmap_transposed) <- gsub("_", " ", colnames(SpeciesScores_heatmap_transposed))

pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\c84_heatmap_SpeciesScores.pdf", width = 20, height = 8.5)
pheatmap(SpeciesScores_heatmap_transposed,
         cluster_rows = FALSE,  
         cluster_cols = FALSE, 
         display_numbers = annotation_matrix, 
         number_color = "black",  
         fontsize_number = 18,  
         color = light_palette(100),
         border_color = "black",
         fontsize_row = 12,
         angle_col = 90,
         fontsize_col = 14)  
dev.off()

write.xlsx(as.data.frame(t(SpeciesScores_heatmap_transposed)),"G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\SpeciesScores_heatmap_transposed.xlsx",row.names = T)

All3Associations_columns <- as.data.frame(SpeciesScores_heatmap_transposed[, colSums(SpeciesScores_heatmap_transposed >= 0.70) == nrow(SpeciesScores_heatmap_transposed)])

###################################################################







