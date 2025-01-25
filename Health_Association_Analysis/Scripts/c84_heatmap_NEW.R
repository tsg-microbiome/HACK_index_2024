

library("scatterplot3d") 
library(writexl)
library(dplyr)
library(xlsx)
load("G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\SpeciesScores_NEW.RData")

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

# common_sp_temp_df <- SpeciesScores_NEW[rownames(SpeciesScores) %in% rownames(heatmap_union_core_df),]
# common_sp_temp <- rownames(SpeciesScores_NEW[rownames(SpeciesScores) %in% rownames(heatmap_union_core_df),])
# test <- heatmap_union_core_df[common_sp_temp,]
# heatmap_union_core_df <- test

library(tidyverse)
union_df_2=heatmap_union_core_df %>%
  pivot_longer(!species, names_to = "Features", values_to = "presence")

catVector= c()
for(rows in 1:nrow(union_df_2))
{
  catVector= c(catVector,ifelse(unlist(union_df_2[rows,"species"]) %in% healthySpecies,"HACK Species", "non HACK Species"))
}

# adding category column in union_df_2
union_df_2= cbind(union_df_2,catVector)

union_df_2$catVector <- factor(union_df_2$catVector, levels = unique(union_df_2$catVector))

union_df_2$species <- factor(union_df_2$species, levels = unique(union_df_2$species))

union_df_2$Features <- factor(union_df_2$Features, levels = c("Health", "Stability", "Influence"))

union_df_2$species_mod <- gsub("_", " ", union_df_2$species)
union_df_2$species_mod <- factor(union_df_2$species_mod, levels = unique(union_df_2$species_mod))

library(ggplot2)

color_vector <- c("white","#56b1f7")    

pdf(file = "G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\c84_heatmap_modified.pdf", width = 20, height = 8.5)
ggplot(union_df_2, aes(x = species_mod, y = Features, fill = factor(presence))) +
  geom_tile() +
  scale_fill_manual(values = color_vector) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5,
                               color = "black")
  ) +
  geom_tile(color = "black") +
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8)) +
  labs(x = "Species", y = "Scores")
dev.off()

filteredHetmap_df$detection= ""

index= which(filteredHetmap_df$Influence==1)

filteredHetmap_df[index,"detection"]= paste0("Influence", filteredHetmap_df[index,"detection"])

index= which(filteredHetmap_df$Stability==1)
filteredHetmap_df[index,"detection"]= paste0(filteredHetmap_df[index,"detection"], ", Stability")

index= which(filteredHetmap_df$Health==1)
filteredHetmap_df[index,"detection"]= paste0(filteredHetmap_df[index,"detection"], ", Health")

filteredHetmap_df$detection= sub("^, ","", filteredHetmap_df$detection)

write.xlsx(filteredHetmap_df, "G:\\My Drive\\HACK paper\\CellReports_revision\\disease_analysis_overall_output\\c84_heatmapDf_modified.xlsx", row.names = T)

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







