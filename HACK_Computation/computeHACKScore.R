library(ggplot2)
library(LaplacesDemon)
library(DescTools)

dir_name <- "G:\\My Drive\\Lab\\Projects\\CoreFinder\\FunctionalProfiling\\"
load(paste0(dir_name,"InfluenceScore.RData"))
load(paste0(dir_name,"StabilityScore.RData"))
load(paste0(dir_name,"healthAssociationScore.RData"))

CommonSpecies <- intersect(intersect(rownames(InfluenceScore),names(StabilityScore)),names(HealthScore))

SpeciesScores <- data.frame(Influence = InfluenceScore[CommonSpecies,],Stability = StabilityScore[CommonSpecies],Health = HealthScore[CommonSpecies])

SpeciesScores$HACKScore <- apply(SpeciesScores,1,mean) * (1-apply(SpeciesScores,1,Gini))

Modes_HACKScore <- Modes(SpeciesScores$HACKScore)$modes
is.trimodal(SpeciesScores$HACKScore)
Size_HACKScore <- Modes(SpeciesScores$HACKScore)$size

save(SpeciesScores,file="G:\\My Drive\\Lab\\Projects\\CoreFinder\\FunctionalProfiling\\SpeciesScores.RData")



