rm(list=ls())

library(DescTools)

load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\SpeciesScores_NEW.RData")
source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\code_library.R")

compute_enrichment <- function(data)
{
	data_new <- apply(data,1,function(x)(mean(x)*(1-Gini(x))))
	return(data_new)
}

SpeciesScores_NEW <- SpeciesScores_NEW[order(SpeciesScores_NEW[,4]),]
##
TopHACKs <- rownames(SpeciesScores_NEW)[SpeciesScores_NEW[,4]>=0.67]
MidHACKs <- rownames(SpeciesScores_NEW)[(SpeciesScores_NEW[,4]<0.)&(SpeciesScores_NEW[,4]>=0.33)]
LowHACKs <- rownames(SpeciesScores_NEW)[SpeciesScores_NEW[,4]<0.33]

Top17 <- rownames(tail(SpeciesScores_NEW,17))
G1 <- rownames(tail(SpeciesScores_NEW,30))
G2 <- rownames(head(tail(SpeciesScores_NEW,60),30))
G3 <- rownames(head(tail(SpeciesScores_NEW,90),30))
G4 <- rownames(head(tail(SpeciesScores_NEW,120),30))
G5 <- rownames(head(tail(SpeciesScores_NEW,150),30))
G6 <- rownames(head(tail(SpeciesScores_NEW,180),30))
G7 <- rownames(head(SpeciesScores_NEW,21))


#FrankelAE_2017
FrankelAE_2017_species <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\FrankelAE_2017.csv",sep=",",row.names=1,header=TRUE)

FrankelAE_2017_metadata <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\FrankelAE_2017.txt",sep="\t",row.names=1,header=TRUE)

FrankelAE_2017_common_rows <- intersect(rownames(FrankelAE_2017_species),rownames(FrankelAE_2017_metadata))

FrankelAE_2017_species <- FrankelAE_2017_species[FrankelAE_2017_common_rows,]
FrankelAE_2017_metadata <- FrankelAE_2017_metadata[FrankelAE_2017_common_rows,]

FrankelAE_2017_species <- as.data.frame(FrankelAE_2017_species[,intersect(rownames(SpeciesScores_NEW),colnames(FrankelAE_2017_species))]/rowSums(FrankelAE_2017_species[,intersect(rownames(SpeciesScores_NEW),colnames(FrankelAE_2017_species))]))

#Temporary rank scaling
FrankelAE_2017_species <- as.data.frame(apply(FrankelAE_2017_species,2,rank_scale))


FrankelAE_2017_sp_responders <- FrankelAE_2017_species[rownames(FrankelAE_2017_metadata)[FrankelAE_2017_metadata$ORR=="yes"],intersect(rownames(SpeciesScores_NEW),colnames(FrankelAE_2017_species))]

FrankelAE_2017_sp_non_responders <- FrankelAE_2017_species[rownames(FrankelAE_2017_metadata)[FrankelAE_2017_metadata$ORR=="no"],intersect(rownames(SpeciesScores_NEW),colnames(FrankelAE_2017_species))]

wilcox_FrankelAE_2017 <- as.data.frame(wilcox_batch(t(FrankelAE_2017_sp_responders),t(FrankelAE_2017_sp_non_responders)))

wilcox_FrankelAE_2017$MeanLogFoldChange <- log((wilcox_FrankelAE_2017[,4]+1e-6)/(wilcox_FrankelAE_2017[,5]+1e-6),10)

wilcox_FrankelAE_2017$HACK <- SpeciesScores_NEW[rownames(wilcox_FrankelAE_2017),4]

df_FrankelAE_2017_HACKGrouped <- data.frame("Top17"=rowMeans(apply(FrankelAE_2017_species[,intersect(Top17,colnames(FrankelAE_2017_species))],2,rank_scale)),"G1"=rowMeans(apply(FrankelAE_2017_species[,intersect(G1,colnames(FrankelAE_2017_species))],2,rank_scale)),"G2"=rowMeans(apply(FrankelAE_2017_species[,intersect(G2,colnames(FrankelAE_2017_species))],2,rank_scale)),"G3"=rowMeans(apply(FrankelAE_2017_species[,intersect(G3,colnames(FrankelAE_2017_species))],2,rank_scale)),"G4"=rowMeans(apply(FrankelAE_2017_species[,intersect(G4,colnames(FrankelAE_2017_species))],2,rank_scale)),"G5"=rowMeans(apply(FrankelAE_2017_species[,intersect(G5,colnames(FrankelAE_2017_species))],2,rank_scale)),"G6"=rowMeans(apply(FrankelAE_2017_species[,intersect(G6,colnames(FrankelAE_2017_species))],2,rank_scale)),"G7"=rowMeans(apply(FrankelAE_2017_species[,intersect(G7,colnames(FrankelAE_2017_species))],2,rank_scale)))

#df_FrankelAE_2017_HACKGrouped <- data.frame("Top"=rowMeans(FrankelAE_2017_species[,intersect(TopHACKs,colnames(FrankelAE_2017_species))]),"Mid"=rowMeans(FrankelAE_2017_species[,intersect(MidHACKs,colnames(FrankelAE_2017_species))]),"Low"=rowMeans(FrankelAE_2017_species[,intersect(LowHACKs,colnames(FrankelAE_2017_species))]))

df_FrankelAE_2017_HACKGrouped$response <- NA
df_FrankelAE_2017_HACKGrouped[rownames(FrankelAE_2017_metadata)[FrankelAE_2017_metadata$ORR=="yes"],"response"] <- "yes"
df_FrankelAE_2017_HACKGrouped[rownames(FrankelAE_2017_metadata)[FrankelAE_2017_metadata$ORR=="no"],"response"] <- "no"
df_FrankelAE_2017_HACKGrouped$study_name <- "FrankelAE_2017"

FrankelAE_2017_species$response <- NA
FrankelAE_2017_species[rownames(FrankelAE_2017_metadata)[FrankelAE_2017_metadata$ORR=="yes"],"response"] <- "yes"
FrankelAE_2017_species[rownames(FrankelAE_2017_metadata)[FrankelAE_2017_metadata$ORR=="no"],"response"] <- "no"
FrankelAE_2017_species$study_name <- "FrankelAE_2017"

#GopalakrishnanV_2018
GopalakrishnanV_2018_species <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\GopalakrishnanV_2018.csv",sep=",",row.names=1,header=TRUE)

GopalakrishnanV_2018_metadata <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\GopalakrishnanV_2018.txt",sep="\t",row.names=1,header=TRUE)

GopalakrishnanV_2018_common_rows <- intersect(rownames(GopalakrishnanV_2018_species),rownames(GopalakrishnanV_2018_metadata))

GopalakrishnanV_2018_species <- GopalakrishnanV_2018_species[GopalakrishnanV_2018_common_rows,]
GopalakrishnanV_2018_metadata <- GopalakrishnanV_2018_metadata[GopalakrishnanV_2018_common_rows,]

GopalakrishnanV_2018_species <- GopalakrishnanV_2018_species[,intersect(rownames(SpeciesScores_NEW),colnames(GopalakrishnanV_2018_species))]/rowSums(GopalakrishnanV_2018_species[,intersect(rownames(SpeciesScores_NEW),colnames(GopalakrishnanV_2018_species))])

#LeeKA_2022
LeeKA_2022_species <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\LeeKA_2022.csv",sep=",",row.names=1,header=TRUE)

LeeKA_2022_metadata <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\LeeKA_2022.txt",sep="\t",row.names=1,header=TRUE)

LeeKA_2022_common_rows <- intersect(rownames(LeeKA_2022_species),rownames(LeeKA_2022_metadata))

LeeKA_2022_species <- LeeKA_2022_species[LeeKA_2022_common_rows,]
LeeKA_2022_metadata <- LeeKA_2022_metadata[LeeKA_2022_common_rows,]

LeeKA_2022_species <- LeeKA_2022_species[,intersect(rownames(SpeciesScores_NEW),colnames(LeeKA_2022_species))]/rowSums(LeeKA_2022_species[,intersect(rownames(SpeciesScores_NEW),colnames(LeeKA_2022_species))])

#Temporary rank scaling
LeeKA_2022_species <- as.data.frame(apply(LeeKA_2022_species,2,rank_scale))


LeeKA_2022_sp_responders <- LeeKA_2022_species[rownames(LeeKA_2022_metadata)[LeeKA_2022_metadata$ORR=="yes"],intersect(rownames(SpeciesScores_NEW),colnames(LeeKA_2022_species))]

LeeKA_2022_sp_non_responders <- LeeKA_2022_species[rownames(LeeKA_2022_metadata)[LeeKA_2022_metadata$ORR=="no"],intersect(rownames(SpeciesScores_NEW),colnames(LeeKA_2022_species))]

wilcox_LeeKA_2022 <- as.data.frame(wilcox_batch(t(LeeKA_2022_sp_responders),t(LeeKA_2022_sp_non_responders)))

wilcox_LeeKA_2022$MeanLogFoldChange <- log((wilcox_LeeKA_2022[,4]+1e-6)/(wilcox_LeeKA_2022[,5]+1e-6),10)

wilcox_LeeKA_2022$HACK <- SpeciesScores_NEW[rownames(wilcox_LeeKA_2022),4]

df_LeeKA_2022_HACKGrouped <- data.frame("Top17"=rowMeans(apply(LeeKA_2022_species[,intersect(Top17,colnames(LeeKA_2022_species))],2,rank_scale)),"G1"=rowMeans(apply(LeeKA_2022_species[,intersect(G1,colnames(LeeKA_2022_species))],2,rank_scale)),"G2"=rowMeans(apply(LeeKA_2022_species[,intersect(G2,colnames(LeeKA_2022_species))],2,rank_scale)),"G3"=rowMeans(apply(LeeKA_2022_species[,intersect(G3,colnames(LeeKA_2022_species))],2,rank_scale)),"G4"=rowMeans(apply(LeeKA_2022_species[,intersect(G4,colnames(LeeKA_2022_species))],2,rank_scale)),"G5"=rowMeans(apply(LeeKA_2022_species[,intersect(G5,colnames(LeeKA_2022_species))],2,rank_scale)),"G6"=rowMeans(apply(LeeKA_2022_species[,intersect(G6,colnames(LeeKA_2022_species))],2,rank_scale)),"G7"=rowMeans(apply(LeeKA_2022_species[,intersect(G7,colnames(LeeKA_2022_species))],2,rank_scale)))

#df_LeeKA_2022_HACKGrouped <- data.frame("Top"=rowMeans(LeeKA_2022_species[,intersect(TopHACKs,colnames(LeeKA_2022_species))]),"Mid"=rowMeans(LeeKA_2022_species[,intersect(MidHACKs,colnames(LeeKA_2022_species))]),"Low"=rowMeans(LeeKA_2022_species[,intersect(LowHACKs,colnames(LeeKA_2022_species))]))


df_LeeKA_2022_HACKGrouped$response <- NA
df_LeeKA_2022_HACKGrouped[rownames(LeeKA_2022_metadata)[LeeKA_2022_metadata$ORR=="yes"],"response"] <- "yes"
df_LeeKA_2022_HACKGrouped[rownames(LeeKA_2022_metadata)[LeeKA_2022_metadata$ORR=="no"],"response"] <- "no"
df_LeeKA_2022_HACKGrouped$study_name <- "LeeKA_2022"

LeeKA_2022_species$response <- NA
LeeKA_2022_species[rownames(LeeKA_2022_metadata)[LeeKA_2022_metadata$ORR=="yes"],"response"] <- "yes"
LeeKA_2022_species[rownames(LeeKA_2022_metadata)[LeeKA_2022_metadata$ORR=="no"],"response"] <- "no"
LeeKA_2022_species$study_name <- "LeeKA_2022"


#MatsonV_2018
MatsonV_2018_species <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\MatsonV_2018.csv",sep=",",row.names=1,header=TRUE)

MatsonV_2018_metadata <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\MatsonV_2018.txt",sep="\t",row.names=1,header=TRUE)

MatsonV_2018_common_rows <- intersect(rownames(MatsonV_2018_species),rownames(MatsonV_2018_metadata))

MatsonV_2018_species <- MatsonV_2018_species[MatsonV_2018_common_rows,]
MatsonV_2018_metadata <- MatsonV_2018_metadata[MatsonV_2018_common_rows,]

MatsonV_2018_species <- MatsonV_2018_species[,intersect(rownames(SpeciesScores_NEW),colnames(MatsonV_2018_species))]/rowSums(MatsonV_2018_species[,intersect(rownames(SpeciesScores_NEW),colnames(MatsonV_2018_species))])

#Temporary rank scaling
MatsonV_2018_species <- as.data.frame(apply(MatsonV_2018_species,2,rank_scale))


MatsonV_2018_sp_responders <- MatsonV_2018_species[rownames(MatsonV_2018_metadata)[MatsonV_2018_metadata$anti_PD_1 =="responder"],intersect(rownames(SpeciesScores_NEW),colnames(MatsonV_2018_species))]

MatsonV_2018_sp_non_responders <- MatsonV_2018_species[rownames(MatsonV_2018_metadata)[MatsonV_2018_metadata$anti_PD_1 =="non_responder"],intersect(rownames(SpeciesScores_NEW),colnames(MatsonV_2018_species))]

wilcox_MatsonV_2018 <- as.data.frame(wilcox_batch(t(MatsonV_2018_sp_responders),t(MatsonV_2018_sp_non_responders)))

wilcox_MatsonV_2018$MeanLogFoldChange <- log((wilcox_MatsonV_2018[,4]+1e-6)/(wilcox_MatsonV_2018[,5]+1e-6),10)

wilcox_MatsonV_2018$HACK <- SpeciesScores_NEW[rownames(wilcox_MatsonV_2018),4]

df_MatsonV_2018_HACKGrouped <- data.frame("Top17"=rowMeans(apply(MatsonV_2018_species[,intersect(Top17,colnames(MatsonV_2018_species))],2,rank_scale)),"G1"=rowMeans(apply(MatsonV_2018_species[,intersect(G1,colnames(MatsonV_2018_species))],2,rank_scale)),"G2"=rowMeans(apply(MatsonV_2018_species[,intersect(G2,colnames(MatsonV_2018_species))],2,rank_scale)),"G3"=rowMeans(apply(MatsonV_2018_species[,intersect(G3,colnames(MatsonV_2018_species))],2,rank_scale)),"G4"=rowMeans(apply(MatsonV_2018_species[,intersect(G4,colnames(MatsonV_2018_species))],2,rank_scale)),"G5"=rowMeans(apply(MatsonV_2018_species[,intersect(G5,colnames(MatsonV_2018_species))],2,rank_scale)),"G6"=rowMeans(apply(MatsonV_2018_species[,intersect(G6,colnames(MatsonV_2018_species))],2,rank_scale)),"G7"=rowMeans(apply(MatsonV_2018_species[,intersect(G7,colnames(MatsonV_2018_species))],2,rank_scale)))

#df_MatsonV_2018_HACKGrouped <- data.frame("Top"=rowMeans(MatsonV_2018_species[,intersect(TopHACKs,colnames(MatsonV_2018_species))]),"Mid"=rowMeans(MatsonV_2018_species[,intersect(MidHACKs,colnames(MatsonV_2018_species))]),"Low"=rowMeans(MatsonV_2018_species[,intersect(LowHACKs,colnames(MatsonV_2018_species))]))

df_MatsonV_2018_HACKGrouped$response <- NA
df_MatsonV_2018_HACKGrouped[rownames(MatsonV_2018_metadata)[MatsonV_2018_metadata$anti_PD_1=="responder"],"response"] <- "yes"
df_MatsonV_2018_HACKGrouped[rownames(MatsonV_2018_metadata)[MatsonV_2018_metadata$anti_PD_1=="non_responder"],"response"] <- "no"
df_MatsonV_2018_HACKGrouped$study_name <- "MatsonV_2018"

MatsonV_2018_species$response <- NA
MatsonV_2018_species[rownames(MatsonV_2018_metadata)[MatsonV_2018_metadata$anti_PD_1=="responder"],"response"] <- "yes"
MatsonV_2018_species[rownames(MatsonV_2018_metadata)[MatsonV_2018_metadata$anti_PD_1=="non_responder"],"response"] <- "no"
MatsonV_2018_species$study_name <- "MatsonV_2018"

#WindTT_2020
WindTT_2020_species <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\WindTT_2020.csv",sep=",",row.names=1,header=TRUE)

WindTT_2020_metadata <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\WindTT_2020.txt",sep="\t",row.names=1,header=TRUE)

WindTT_2020_common_rows <- intersect(rownames(WindTT_2020_species),rownames(WindTT_2020_metadata))

WindTT_2020_species <- WindTT_2020_species[WindTT_2020_common_rows,]
WindTT_2020_metadata <- WindTT_2020_metadata[WindTT_2020_common_rows,]

WindTT_2020_species <- WindTT_2020_species[,intersect(rownames(SpeciesScores_NEW),colnames(WindTT_2020_species))]/rowSums(WindTT_2020_species[,intersect(rownames(SpeciesScores_NEW),colnames(WindTT_2020_species))])

#Temporary rank scaling
WindTT_2020_species <- as.data.frame(apply(WindTT_2020_species,2,rank_scale))

WindTT_2020_sp_responders <- WindTT_2020_species[rownames(WindTT_2020_metadata)[WindTT_2020_metadata$ORR=="yes"],intersect(rownames(SpeciesScores_NEW),colnames(WindTT_2020_species))]

WindTT_2020_sp_non_responders <- WindTT_2020_species[rownames(WindTT_2020_metadata)[WindTT_2020_metadata$ORR=="no"],intersect(rownames(SpeciesScores_NEW),colnames(WindTT_2020_species))]

wilcox_WindTT_2020 <- as.data.frame(wilcox_batch(t(WindTT_2020_sp_responders),t(WindTT_2020_sp_non_responders)))

wilcox_WindTT_2020$MeanLogFoldChange <- log((wilcox_WindTT_2020[,4]+1e-6)/(wilcox_WindTT_2020[,5]+1e-6),10)

wilcox_WindTT_2020$HACK <- SpeciesScores_NEW[rownames(wilcox_WindTT_2020),4]

df_WindTT_2020_HACKGrouped <- data.frame("Top17"=rowMeans(apply(WindTT_2020_species[,intersect(Top17,colnames(WindTT_2020_species))],2,rank_scale)),"G1"=rowMeans(apply(WindTT_2020_species[,intersect(G1,colnames(WindTT_2020_species))],2,rank_scale)),"G2"=rowMeans(apply(WindTT_2020_species[,intersect(G2,colnames(WindTT_2020_species))],2,rank_scale)),"G3"=rowMeans(apply(WindTT_2020_species[,intersect(G3,colnames(WindTT_2020_species))],2,rank_scale)),"G4"=rowMeans(apply(WindTT_2020_species[,intersect(G4,colnames(WindTT_2020_species))],2,rank_scale)),"G5"=rowMeans(apply(WindTT_2020_species[,intersect(G5,colnames(WindTT_2020_species))],2,rank_scale)),"G6"=rowMeans(apply(WindTT_2020_species[,intersect(G6,colnames(WindTT_2020_species))],2,rank_scale)),"G7"=rowMeans(apply(WindTT_2020_species[,intersect(G7,colnames(WindTT_2020_species))],2,rank_scale)))

#df_WindTT_2020_HACKGrouped <- data.frame("Top"=rowMeans(WindTT_2020_species[,intersect(TopHACKs,colnames(WindTT_2020_species))]),"Mid"=rowMeans(WindTT_2020_species[,intersect(MidHACKs,colnames(WindTT_2020_species))]),"Low"=rowMeans(WindTT_2020_species[,intersect(LowHACKs,colnames(WindTT_2020_species))]))

df_WindTT_2020_HACKGrouped$response <- NA
df_WindTT_2020_HACKGrouped[rownames(WindTT_2020_metadata)[WindTT_2020_metadata$ORR=="yes"],"response"] <- "yes"
df_WindTT_2020_HACKGrouped[rownames(WindTT_2020_metadata)[WindTT_2020_metadata$ORR=="no"],"response"] <- "no"
df_WindTT_2020_HACKGrouped$study_name <- "WindTT_2020"

WindTT_2020_species$response <- NA
WindTT_2020_species[rownames(WindTT_2020_metadata)[WindTT_2020_metadata$ORR=="yes"],"response"] <- "yes"
WindTT_2020_species[rownames(WindTT_2020_metadata)[WindTT_2020_metadata$ORR=="no"],"response"] <- "no"
WindTT_2020_species$study_name <- "WindTT_2020"


#RoutyB_2018
RoutyB_2018_species <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\Routy_Species.txt",sep="\t",row.names=1,header=TRUE)

RoutyB_2018_metadata <- read.table("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\AntiCancer\\Routy_Metadata.txt",sep="\t",row.names=1,header=TRUE)

RoutyB_2018_common_rows <- intersect(rownames(RoutyB_2018_species),rownames(RoutyB_2018_metadata))

RoutyB_2018_species <- RoutyB_2018_species[RoutyB_2018_common_rows,]
RoutyB_2018_metadata <- RoutyB_2018_metadata[RoutyB_2018_common_rows,]

RoutyB_2018_species <- RoutyB_2018_species[,intersect(rownames(SpeciesScores_NEW),colnames(RoutyB_2018_species))]/rowSums(RoutyB_2018_species[,intersect(rownames(SpeciesScores_NEW),colnames(RoutyB_2018_species))])

#Temporary rank scaling
RoutyB_2018_species <- as.data.frame(apply(RoutyB_2018_species,2,rank_scale))

RoutyB_2018_sp_responders <- RoutyB_2018_species[rownames(RoutyB_2018_metadata)[RoutyB_2018_metadata$anti_PD_1 =="responder"],intersect(rownames(SpeciesScores_NEW),colnames(RoutyB_2018_species))]

RoutyB_2018_sp_non_responders <- RoutyB_2018_species[rownames(RoutyB_2018_metadata)[RoutyB_2018_metadata$anti_PD_1 =="non_responder"],intersect(rownames(SpeciesScores_NEW),colnames(RoutyB_2018_species))]

wilcox_RoutyB_2018 <- as.data.frame(wilcox_batch(t(RoutyB_2018_sp_responders),t(RoutyB_2018_sp_non_responders)))

wilcox_RoutyB_2018$MeanLogFoldChange <- log((wilcox_RoutyB_2018[,4]+1e-6)/(wilcox_RoutyB_2018[,5]+1e-6),10)

wilcox_RoutyB_2018$HACK <- SpeciesScores_NEW[rownames(wilcox_RoutyB_2018),4]

df_RoutyB_2018_HACKGrouped <- data.frame("Top17"=rowMeans(apply(RoutyB_2018_species[,intersect(Top17,colnames(RoutyB_2018_species))],2,rank_scale)),"G1"=rowMeans(apply(RoutyB_2018_species[,intersect(G1,colnames(RoutyB_2018_species))],2,rank_scale)),"G2"=rowMeans(apply(RoutyB_2018_species[,intersect(G2,colnames(RoutyB_2018_species))],2,rank_scale)),"G3"=rowMeans(apply(RoutyB_2018_species[,intersect(G3,colnames(RoutyB_2018_species))],2,rank_scale)),"G4"=rowMeans(apply(RoutyB_2018_species[,intersect(G4,colnames(RoutyB_2018_species))],2,rank_scale)),"G5"=rowMeans(apply(RoutyB_2018_species[,intersect(G5,colnames(RoutyB_2018_species))],2,rank_scale)),"G6"=rowMeans(apply(RoutyB_2018_species[,intersect(G6,colnames(RoutyB_2018_species))],2,rank_scale)),"G7"=rowMeans(apply(RoutyB_2018_species[,intersect(G7,colnames(RoutyB_2018_species))],2,rank_scale)))

#df_RoutyB_2018_HACKGrouped <- data.frame("Top"=rowMeans(RoutyB_2018_species[,intersect(TopHACKs,colnames(RoutyB_2018_species))]),"Mid"=rowMeans(RoutyB_2018_species[,intersect(MidHACKs,colnames(RoutyB_2018_species))]),"Low"=rowMeans(RoutyB_2018_species[,intersect(LowHACKs,colnames(RoutyB_2018_species))]))

df_RoutyB_2018_HACKGrouped$response <- NA
df_RoutyB_2018_HACKGrouped[rownames(RoutyB_2018_metadata)[RoutyB_2018_metadata$anti_PD_1=="responder"],"response"] <- "yes"
df_RoutyB_2018_HACKGrouped[rownames(RoutyB_2018_metadata)[RoutyB_2018_metadata$anti_PD_1=="non_responder"],"response"] <- "no"
df_RoutyB_2018_HACKGrouped$study_name <- "RoutyB_2018"

RoutyB_2018_species$response <- NA
RoutyB_2018_species[rownames(RoutyB_2018_metadata)[RoutyB_2018_metadata$anti_PD_1=="responder"],"response"] <- "yes"
RoutyB_2018_species[rownames(RoutyB_2018_metadata)[RoutyB_2018_metadata$anti_PD_1=="non_responder"],"response"] <- "no"
RoutyB_2018_species$study_name <- "RoutyB_2018"

#Merging of Grouped Data and Meta-Analysis
temp0 <- merge(t(df_RoutyB_2018_HACKGrouped[,c(1:8)]),t(df_FrankelAE_2017_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(t(df_RoutyB_2018_HACKGrouped[,c(1:8)]),t(df_FrankelAE_2017_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(t(temp0),t(df_LeeKA_2022_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(t(temp0),t(df_LeeKA_2022_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,1]
temp1 <- apply(temp1,1,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(t(temp1),t(df_MatsonV_2018_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(t(temp1),t(df_MatsonV_2018_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,1,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(t(temp2),t(df_WindTT_2020_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(t(temp2),t(df_WindTT_2020_HACKGrouped[,c(1:8)]),by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,1,function(x)(ifelse(is.na(x),0,x)))

df_AllStudies_HACKGrouped <- as.data.frame(temp3)

df_AllStudies_HACKGrouped$response <- NA
df_AllStudies_HACKGrouped[rownames(df_FrankelAE_2017_HACKGrouped),"response"] <- df_FrankelAE_2017_HACKGrouped$response
df_AllStudies_HACKGrouped[rownames(df_LeeKA_2022_HACKGrouped),"response"] <- df_LeeKA_2022_HACKGrouped$response
df_AllStudies_HACKGrouped[rownames(df_MatsonV_2018_HACKGrouped),"response"] <- df_MatsonV_2018_HACKGrouped$response
df_AllStudies_HACKGrouped[rownames(df_WindTT_2020_HACKGrouped),"response"] <- df_WindTT_2020_HACKGrouped$response
df_AllStudies_HACKGrouped[rownames(df_RoutyB_2018_HACKGrouped),"response"] <- df_RoutyB_2018_HACKGrouped$response

df_AllStudies_HACKGrouped$study_name <- NA
df_AllStudies_HACKGrouped[rownames(df_FrankelAE_2017_HACKGrouped),"study_name"] <- df_FrankelAE_2017_HACKGrouped$study_name
df_AllStudies_HACKGrouped[rownames(df_LeeKA_2022_HACKGrouped),"study_name"] <- df_LeeKA_2022_HACKGrouped$study_name
df_AllStudies_HACKGrouped[rownames(df_MatsonV_2018_HACKGrouped),"study_name"] <- df_MatsonV_2018_HACKGrouped$study_name
df_AllStudies_HACKGrouped[rownames(df_WindTT_2020_HACKGrouped),"study_name"] <- df_WindTT_2020_HACKGrouped$study_name
df_AllStudies_HACKGrouped[rownames(df_RoutyB_2018_HACKGrouped),"study_name"] <- df_RoutyB_2018_HACKGrouped$study_name

df_AllStudies_HACKGrouped$binary_response <- ifelse(df_AllStudies_HACKGrouped$response == "yes",1,0)

Groups_Batch_REM_Response <- batch_rem2_grouped(df_AllStudies_HACKGrouped,df_AllStudies_HACKGrouped,c("G1","G2","G3","G4","G5","G6","G7"),"binary_response","study_name",unique(df_AllStudies_HACKGrouped$study_name))

#Merging of Grouped Data and Meta-Analysis
temp0 <- merge(t(RoutyB_2018_species[,setdiff(colnames(RoutyB_2018_species),c("study_name","study_name"))]),t(FrankelAE_2017_species[,setdiff(colnames(FrankelAE_2017_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(t(RoutyB_2018_species[,setdiff(colnames(RoutyB_2018_species),c("study_name","study_name"))]),t(FrankelAE_2017_species[,setdiff(colnames(FrankelAE_2017_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(t(temp0),t(LeeKA_2022_species[,setdiff(colnames(LeeKA_2022_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(t(temp0),t(LeeKA_2022_species[,setdiff(colnames(LeeKA_2022_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,1]
temp1 <- apply(temp1,1,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(t(temp1),t(MatsonV_2018_species[,setdiff(colnames(MatsonV_2018_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(t(temp1),t(MatsonV_2018_species[,setdiff(colnames(MatsonV_2018_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,1,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(t(temp2),t(WindTT_2020_species[,setdiff(colnames(WindTT_2020_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(t(temp2),t(WindTT_2020_species[,setdiff(colnames(WindTT_2020_species),c("study_name","study_name"))]),by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,1,function(x)(ifelse(is.na(x),0,x)))

AllStudies_species <- as.data.frame(temp3)

AllStudies_species$response <- NA
AllStudies_species[rownames(FrankelAE_2017_species),"response"] <- FrankelAE_2017_species$response
AllStudies_species[rownames(LeeKA_2022_species),"response"] <- LeeKA_2022_species$response
AllStudies_species[rownames(MatsonV_2018_species),"response"] <- MatsonV_2018_species$response
AllStudies_species[rownames(WindTT_2020_species),"response"] <- WindTT_2020_species$response
AllStudies_species[rownames(RoutyB_2018_species),"response"] <- RoutyB_2018_species$response

AllStudies_species$study_name <- NA
AllStudies_species[rownames(FrankelAE_2017_species),"study_name"] <- FrankelAE_2017_species$study_name
AllStudies_species[rownames(LeeKA_2022_species),"study_name"] <- LeeKA_2022_species$study_name
AllStudies_species[rownames(MatsonV_2018_species),"study_name"] <- MatsonV_2018_species$study_name
AllStudies_species[rownames(WindTT_2020_species),"study_name"] <- WindTT_2020_species$study_name
AllStudies_species[rownames(RoutyB_2018_species),"study_name"] <- RoutyB_2018_species$study_name

AllStudies_species$binary_response <- ifelse(AllStudies_species$response == "yes",1,0)

Species_REM_Response <- as.data.frame(batch_rem2_grouped(AllStudies_species,AllStudies_species,intersect(rownames(SpeciesScores_NEW),colnames(AllStudies_species)),"binary_response","study_name",unique(AllStudies_species$study_name)))

Species_REM_Response$TaxaType <- ifelse(rownames(Species_REM_Response) %in% Top17,"Darkolivegreen2","firebrick3")

Species_REM_Response$HACK <- SpeciesScores_NEW[rownames(Species_REM_Response),4]
Species_REM_Response$Influence <- SpeciesScores_NEW[rownames(Species_REM_Response),1]
Species_REM_Response$Stability <- SpeciesScores_NEW[rownames(Species_REM_Response),2]
Species_REM_Response$Health <- SpeciesScores_NEW[rownames(Species_REM_Response),3]

ggplot(Species_REM_Response[Species_REM_Response$consistency>=0.8,],aes(x=est,y=-log(pval,10)))+geom_point(col="blue4",size=2)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

ggplot(Species_REM_Response[Species_REM_Response$consistency>=0.8,],aes(x=est,y=Influence))+geom_point(col="blue4",size=2)+geom_smooth(method='lm')+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

ggplot(Groups_Batch_REM_Response,aes(x=est,y=-log(pval,10)))+geom_point(size=2)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+geom_text_repel(label=rownames(Groups_Batch_REM_Response)) + geom_hline(yintercept=1.3010,col="deepskyblue")

Study_Names <- unique(AllStudies_species$study_name)

Top17 <- rownames(tail(SpeciesScores_NEW,17))
G1 <- rownames(tail(SpeciesScores_NEW,30))
G2 <- rownames(head(tail(SpeciesScores_NEW,60),30))
G3 <- rownames(head(tail(SpeciesScores_NEW,90),30))
G4 <- rownames(head(tail(SpeciesScores_NEW,120),30))
G5 <- rownames(head(tail(SpeciesScores_NEW,150),30))
G6 <- rownames(head(tail(SpeciesScores_NEW,180),30))
G7 <- rownames(head(SpeciesScores_NEW,21))

FrankelAE_2017_species_grouped <- data.frame("Top17"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),Top17)],2,rank_scale)),"G1"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G1)],2,rank_scale)),"G2"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G2)],2,rank_scale)),"G3"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G3)],2,rank_scale)),"G4"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G4)],2,rank_scale)),"G5"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G5)],2,rank_scale)),"G6"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G6)],2,rank_scale)),"G7"=compute_enrichment(apply(FrankelAE_2017_species[,intersect(colnames(FrankelAE_2017_species),G7)],2,rank_scale)))
FrankelAE_2017_species_grouped$study_name <- "FrankelAE_2017"
FrankelAE_2017_species_grouped$binary_response <- ifelse(rownames(FrankelAE_2017_species_grouped) %in% rownames(FrankelAE_2017_sp_responders),1,0)

LeeKA_2022_species_grouped <- data.frame("Top17"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),Top17)],2,rank_scale)),"G1"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G1)],2,rank_scale)),"G2"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G2)],2,rank_scale)),"G3"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G3)],2,rank_scale)),"G4"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G4)],2,rank_scale)),"G5"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G5)],2,rank_scale)),"G6"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G6)],2,rank_scale)),"G7"=compute_enrichment(apply(LeeKA_2022_species[,intersect(colnames(LeeKA_2022_species),G7)],2,rank_scale)))
LeeKA_2022_species_grouped$study_name <- "LeeKA_2022"
LeeKA_2022_species_grouped$binary_response <- ifelse(rownames(LeeKA_2022_species_grouped) %in% rownames(LeeKA_2022_sp_responders),1,0)

RoutyB_2018_species_grouped <- data.frame("Top17"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),Top17)],2,rank_scale)),"G1"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G1)],2,rank_scale)),"G2"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G2)],2,rank_scale)),"G3"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G3)],2,rank_scale)),"G4"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G4)],2,rank_scale)),"G5"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G5)],2,rank_scale)),"G6"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G6)],2,rank_scale)),"G7"=compute_enrichment(apply(RoutyB_2018_species[,intersect(colnames(RoutyB_2018_species),G7)],2,rank_scale)))
RoutyB_2018_species_grouped$study_name <- "RoutyB_2018"
RoutyB_2018_species_grouped$binary_response <- ifelse(rownames(RoutyB_2018_species_grouped) %in% rownames(RoutyB_2018_sp_responders),1,0)

WindTT_2020_species_grouped <- data.frame("Top17"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),Top17)],2,rank_scale)),"G1"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G1)],2,rank_scale)),"G2"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G2)],2,rank_scale)),"G3"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G3)],2,rank_scale)),"G4"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G4)],2,rank_scale)),"G5"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G5)],2,rank_scale)),"G6"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G6)],2,rank_scale)),"G7"=compute_enrichment(apply(WindTT_2020_species[,intersect(colnames(WindTT_2020_species),G7)],2,rank_scale)))
WindTT_2020_species_grouped$study_name <- "WindTT_2020"
WindTT_2020_species_grouped$binary_response <- ifelse(rownames(WindTT_2020_species_grouped) %in% rownames(WindTT_2020_sp_responders),1,0)

MatsonV_2018_species_grouped <- data.frame("Top17"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),Top17)],2,rank_scale)),"G1"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G1)],2,rank_scale)),"G2"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G2)],2,rank_scale)),"G3"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G3)],2,rank_scale)),"G4"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G4)],2,rank_scale)),"G5"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G5)],2,rank_scale)),"G6"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G6)],2,rank_scale)),"G7"=compute_enrichment(apply(MatsonV_2018_species[,intersect(colnames(MatsonV_2018_species),G7)],2,rank_scale)))
MatsonV_2018_species_grouped$study_name <- "MatsonV_2018"
MatsonV_2018_species_grouped$binary_response <- ifelse(rownames(MatsonV_2018_species_grouped) %in% rownames(MatsonV_2018_sp_responders),1,0)

AllStudies_species_grouped <- as.data.frame(rbind(FrankelAE_2017_species_grouped,LeeKA_2022_species_grouped,RoutyB_2018_species_grouped,WindTT_2020_species_grouped,MatsonV_2018_species_grouped))

SpeciesGroup_REM_Response <- as.data.frame(batch_rem2_grouped(AllStudies_species_grouped,AllStudies_species_grouped,c("G1","G2","G3","G4","G5","G6","G7"),"binary_response","study_name",unique(AllStudies_species_grouped$study_name)))

Directions_Species <- as.data.frame(matrix(0,length(intersect(rownames(SpeciesScores_NEW),colnames(AllStudies_species))),length(Study_Names)))
rownames(Directions_Species) <- intersect(rownames(SpeciesScores_NEW),colnames(AllStudies_species))
colnames(Directions_Species) <- Study_Names

for(i in 1:nrow(Directions_Species))
{
	Species <- rownames(Directions_Species)[i]
	for(j in 1:ncol(Directions_Species))
	{
		Study <- Study_Names[j]
		Wilcox <- get(paste0("wilcox_",Study))
		if(Species %in% rownames(Wilcox))
		{
			Directions_Species[Species,Study] <- ifelse(Wilcox[Species,3]<=0.15,3*sign(Wilcox[Species,2]),ifelse(Wilcox[Species,1]<=0.05,2*sign(Wilcox[Species,2]),sign(Wilcox[Species,2])))
		}
	}
}

Group_colors <- c("Violet","Slateblue4","Blue","Green","Yellow","Orange","Red")
names(Group_colors) <- c("G1","G2","G3","G4","G5","G6","G7")

Species_REM_Response[intersect(G1,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G1"]
Species_REM_Response[intersect(G2,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G2"]
Species_REM_Response[intersect(G3,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G3"]
Species_REM_Response[intersect(G4,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G4"]
Species_REM_Response[intersect(G5,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G5"]
Species_REM_Response[intersect(G6,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G6"]
Species_REM_Response[intersect(G7,rownames(Species_REM_Response)),"TaxaType"] <- Group_colors["G7"]

ggplot(Species_REM_Response[,],aes(x=est,y=-log(pval,10)))+geom_point(col=Species_REM_Response[,"TaxaType"],size=2)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15)) + geom_hline(yintercept=-log(0.05,10))+geom_text_repel(label=ifelse(Species_REM_Response[,"pval"]<=0.05,rownames(Species_REM_Response),""),max.overlaps=20)

ggplot(Species_REM_Response[Species_REM_Response$consistency>=0.6,],aes(x=HACK,y=est))+geom_point(color=Species_REM_Response[Species_REM_Response$consistency>=0.6,"TaxaType"],size=2)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

Species_REM_Response$Group <- NA
Species_REM_Response[intersect(G1,rownames(Species_REM_Response)),"Group"] <- "G1"
Species_REM_Response[intersect(G2,rownames(Species_REM_Response)),"Group"] <- "G2"
Species_REM_Response[intersect(G3,rownames(Species_REM_Response)),"Group"] <- "G3"
Species_REM_Response[intersect(G4,rownames(Species_REM_Response)),"Group"] <- "G4"
Species_REM_Response[intersect(G5,rownames(Species_REM_Response)),"Group"] <- "G5"
Species_REM_Response[intersect(G6,rownames(Species_REM_Response)),"Group"] <- "G6"
Species_REM_Response[intersect(G7,rownames(Species_REM_Response)),"Group"] <- "G7"

Filtered_Species_REM_Response <- Species_REM_Response[!is.na(Species_REM_Response$consistency)&(Species_REM_Response$consistency>=0.6),c("est","pval","Group","TaxaType")]

Filtered_Species_REM_Response$PointType <- "Taxa"

Groups_Batch_REM_Response$TaxaType <- NA
Groups_Batch_REM_Response["G1","TaxaType"] <- Group_colors["G1"]
Groups_Batch_REM_Response["G2","TaxaType"] <- Group_colors["G2"]
Groups_Batch_REM_Response["G3","TaxaType"] <- Group_colors["G3"]
Groups_Batch_REM_Response["G4","TaxaType"] <- Group_colors["G4"]
Groups_Batch_REM_Response["G5","TaxaType"] <- Group_colors["G5"]
Groups_Batch_REM_Response["G6","TaxaType"] <- Group_colors["G6"]
Groups_Batch_REM_Response["G7","TaxaType"] <- Group_colors["G7"]

Groups_Batch_REM_Response$Group <- NA
Groups_Batch_REM_Response["G1","Group"] <- "G1"
Groups_Batch_REM_Response["G2","Group"] <- "G2"
Groups_Batch_REM_Response["G3","Group"] <- "G3"
Groups_Batch_REM_Response["G4","Group"] <- "G4"
Groups_Batch_REM_Response["G5","Group"] <- "G5"
Groups_Batch_REM_Response["G6","Group"] <- "G6"
Groups_Batch_REM_Response["G7","Group"] <- "G7"

Groups_Batch_REM_Response$PointType <- "Group"

ggplot(Groups_Batch_REM_Response[Groups_Batch_REM_Response$TaxaType != "Top17",],aes(x=est,y=-log(pval,10)))+geom_point(col=Groups_Batch_REM_Response[Groups_Batch_REM_Response$TaxaType != "Top17",],size=2)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15)) + geom_hline(yintercept=-log(0.05,10))+geom_text_repel(labels=rownames(Groups_Batch_REM_Response[Groups_Batch_REM_Response$TaxaType != "Top17",]))

Cumulate_Directions <- data.frame("Responder_Positive"=apply(Directions_Species,1,function(x)(length(x[x>=1]))),"Responder_Negative"=apply(Directions_Species,1,function(x)(length(x[x<=-1]))))

boxplot(SpeciesScores_NEW[rownames(Cumulate_Directions),4]~cut(Cumulate_Directions[,1]-Cumulate_Directions[,2],breaks=c(-5,-2,0,2,5),include.lowest=TRUE),xlab="",ylab="",col=c("Violetred2","Yellow2","Slateblue2","Darkolivegreen"))












