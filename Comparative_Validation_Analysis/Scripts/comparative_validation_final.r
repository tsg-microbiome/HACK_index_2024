library(DescTools)
library(dunn.test)
#source("G:\\My Drive\\Lab\\Projects\\CoreFinder\\Manuscript\\Analysis_Data\\code_library.R")
source("code_library.R")

#library(vegan)
library(pcaPP)
library(dplyr)

diversity <- vegan::diversity
kendall_uniqueness <- function(data)
{
	dist_mat <- 1-cor.fk(t(data))
	diag(dist_mat) <- NA
	kendall_uniqueness <- apply(dist_mat,1,function(x)(min(x[!is.na(x)])))
	return(kendall_uniqueness)
}

dysbiosis_score <- function(data,control_group)
{
	dist_mat <- as.matrix(vegdist(data,method="bray"))
	diag(dist_mat) <- NA
	dysbiosis_score <- apply(dist_mat[,control_group],1,function(x)(median(x[!is.na(x)])))
	return(dysbiosis_score)
}

GMHI <- function(data) {
  
  library(dplyr)
  
  tmp1 <- data.frame(t(data),check.rows = F,check.names = F) 
  
  MH_species <- c("Alistipes_senegalensis","Bacteroidales_bacterium_ph8","Bifidobacterium_adolescentis","Bifidobacterium_angulatum","Bifidobacterium_catenulatum","Lachnospiraceae_bacterium_8_1_57FAA","Sutterella_wadsworthensis")
  
  MN_species <- c("Anaerotruncus_colihominis","Atopobium_parvulum","Bifidobacterium_dentium","Blautia_producta","candidate_division_TM7_single_cell_isolate_TM7c","Clostridiales_bacterium_1_7_47FAA","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_clostridioforme","Clostridium_hathewayi","Clostridium_nexile","Clostridium_ramosum","Clostridium_symbiosum","Eggerthella_lenta","Erysipelotrichaceae_bacterium_2_2_44A","Flavonifractor_plautii","Fusobacterium_nucleatum","Gemella_morbillorum","Gemella_sanguinis","Granulicatella_adiacens","Holdemania_filiformis","Klebsiella_pneumoniae","Lachnospiraceae_bacterium_1_4_56FAA","Lachnospiraceae_bacterium_2_1_58FAA","Lachnospiraceae_bacterium_3_1_57FAA_CT1","Lachnospiraceae_bacterium_5_1_57FAA","Lachnospiraceae_bacterium_9_1_43BFAA","Lactobacillus_salivarius","Peptostreptococcus_stomatis","Ruminococcaceae_bacterium_D16","Ruminococcus_gnavus","Solobacterium_moorei","Streptococcus_anginosus","Streptococcus_australis","Streptococcus_gordonii","Streptococcus_infantis","Streptococcus_mitis","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_sp_4_3_54A2FAA","Subdoligranulum_variabile","Veillonella_atypica")
  
  
  MH_species_metagenome <- tmp1[row.names(tmp1) %in% MH_species, ]
  MN_species_metagenome <- tmp1[row.names(tmp1) %in% MN_species, ]
  
  
  alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
  MH_shannon <- apply((MH_species_metagenome), 2, alpha) 
  MN_shannon <- apply((MN_species_metagenome), 2, alpha) 
  
  R_MH <- apply(MH_species_metagenome, 2, function(i) (sum(i > 0))) 
  R_MN <- apply(MN_species_metagenome, 2, function(i) (sum(i > 0)))
  
  
  #By default from methods
  MH_prime <- 7
  MN_prime <- 31
  
  psi_MH <- ((R_MH/MH_prime)*MH_shannon) 
  psi_MN <- ((R_MN/MN_prime)*MN_shannon)
  
  GMHI <- suppressWarnings(data.frame(log10((psi_MH+0.00001)/(psi_MN+0.00001))))
  colnames(GMHI) <- c("GMHI")
  GMHI$GMHI[is.nan(GMHI$GMHI)] <- 0
  return(GMHI)
}



hack_top_17 <- function(data)
{
	hack_top_17 <- c("Faecalibacterium_prausnitzii","Bacteroides_uniformis","Odoribacter_splanchnicus","Fusicatenibacter_saccharivorans","Coprococcus_catus","Eubacterium_rectale","Oscillibacter_sp_57_20","Eubacterium_ramulus","Agathobaculum_butyriciproducens","Roseburia_inulinivorans","Roseburia_intestinalis","Ruminococcus_bromii","Roseburia_hominis","Alistipes_shahii","Eubacterium_eligens","Roseburia_faecis","Alistipes_putredinis")
	
	#hack_top_17 <- c("Faecalibacterium_prausnitzii","Bacteroides_uniformis","Odoribacter_splanchnicus","Fusicatenibacter_saccharivorans","Coprococcus_catus","Eubacterium_rectale","Oscillibacter_sp_57_20","Eubacterium_ramulus","Agathobaculum_butyriciproducens","Roseburia_inulinivorans","Roseburia_intestinalis","Ruminococcus_bromii")
	
	hack_top_17_score <- rowMeans(apply(data[,intersect(colnames(data),hack_top_17)],2,rank_scale))
	
	return(hack_top_17_score)
}

GMHI_hack_top_17 <- function(data)
{
	MH_Species <- c("Alistipes_senegalensis","Bacteroidales_bacterium_ph8","Bifidobacterium_adolescentis","Bifidobacterium_angulatum","Bifidobacterium_catenulatum","Lachnospiraceae_bacterium_8_1_57FAA","Sutterella_wadsworthensis")
	
	hack_top_17 <- c("Faecalibacterium_prausnitzii","Bacteroides_uniformis","Odoribacter_splanchnicus","Fusicatenibacter_saccharivorans","Coprococcus_catus","Eubacterium_rectale","Oscillibacter_sp_57_20","Eubacterium_ramulus","Agathobaculum_butyriciproducens","Roseburia_inulinivorans","Roseburia_intestinalis","Ruminococcus_bromii","Roseburia_hominis","Alistipes_shahii","Eubacterium_eligens","Roseburia_faecis","Alistipes_putredinis","Parabacteroides_merdae","Eubacterium_ventriosum","Gemmiger_formicilis")
	
	MH_Species <- unique(c(hack_top_17,MH_Species))

	MN_Species <- c("Anaerotruncus_colihominis","Atopobium_parvulum","Bifidobacterium_dentium","Blautia_producta","candidate_division_TM7_single_cell_isolate_TM7c","Clostridiales_bacterium_1_7_47FAA","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_clostridioforme","Clostridium_hathewayi","Clostridium_nexile","Clostridium_ramosum","Clostridium_symbiosum","Eggerthella_lenta","Erysipelotrichaceae_bacterium_2_2_44A","Flavonifractor_plautii","Fusobacterium_nucleatum","Gemella_morbillorum","Gemella_sanguinis","Granulicatella_adiacens","Holdemania_filiformis","Klebsiella_pneumoniae","Lachnospiraceae_bacterium_1_4_56FAA","Lachnospiraceae_bacterium_2_1_58FAA","Lachnospiraceae_bacterium_3_1_57FAA_CT1","Lachnospiraceae_bacterium_5_1_57FAA","Lachnospiraceae_bacterium_9_1_43BFAA","Lactobacillus_salivarius","Peptostreptococcus_stomatis","Ruminococcaceae_bacterium_D16","Ruminococcus_gnavus","Solobacterium_moorei","Streptococcus_anginosus","Streptococcus_australis","Streptococcus_gordonii","Streptococcus_infantis","Streptococcus_mitis","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_sp_4_3_54A2FAA","Subdoligranulum_variabile","Veillonella_atypica")

	MH_Diversity <- diversity(data[,intersect(MH_Species,colnames(data))])
	MN_Diversity <- diversity(data[,intersect(MN_Species,colnames(data))])

	Weighted_MH_Diversity <- MH_Diversity * (length(intersect(MH_Species,colnames(data)))/length(MH_Species))

	Weighted_MN_Diversity <- MN_Diversity * (length(intersect(MN_Species,colnames(data)))/length(MN_Species))

	GMHI <- log(((Weighted_MH_Diversity+0.00001)/(Weighted_MN_Diversity+0.00001)),10)

	return(GMHI)
}

print("Loading Datasets")	
print("Flemer") 
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Flemer_species.RData")
load("Flemer_species.RData")
Flemer_CRC <- rownames(Flemer_et_al_metadata[Flemer_et_al_metadata[,1]=="CRC",])
Flemer_Control <- rownames(Flemer_et_al_metadata[Flemer_et_al_metadata[,1]=="control",])
Flemer_Summary_Statistics <- data.frame("Shannon"=diversity(Flemer_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Flemer_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Flemer_et_al_species,Flemer_Control),"GMHI"=GMHI(Flemer_et_al_species),"hack_top_17_score"=hack_top_17(Flemer_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Flemer_et_al_species),"Seq_Type"="16S","Cohort_Type"="I","Study_Name"="Flemer")

wilcox_Flemer <- wilcox_batch(t(Flemer_Summary_Statistics[Flemer_CRC,1:6]),t(Flemer_Summary_Statistics[Flemer_Control,1:6]))

print("Nagata")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Nagata_species.RData")
load("Nagata_species.RData")
Nagata_Control <- rownames(Nagata_et_al_metadata[Nagata_et_al_metadata[,1]=="Control",])
Nagata_PDAC <- rownames(Nagata_et_al_metadata[Nagata_et_al_metadata[,1]=="PDAC",])
Nagata_Summary_Statistics <- data.frame("Shannon"=diversity(Nagata_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Nagata_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Nagata_et_al_species,Nagata_Control),"GMHI"=GMHI(Nagata_et_al_species),"hack_top_17_score"=hack_top_17(Nagata_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Nagata_et_al_species),"Seq_Type"="WGS","Cohort_Type"="I","Study_Name"="Nagata")

wilcox_Nagata <- wilcox_batch(t(Nagata_Summary_Statistics[Nagata_PDAC,1:6]),t(Nagata_Summary_Statistics[Nagata_Control,1:6]))

print("MicroDiab")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\MicroDiab_Denmark_species.RData")
load("MicroDiab_Denmark_species.RData")
MicroDiab_Denmark_et_al_species <- MicroDiab_Denmark_Species/rowSums(MicroDiab_Denmark_Species)
MicroDiab_Denmark_Control <- rownames(MicroDiab_Denmark_Metadata[MicroDiab_Denmark_Metadata[,3] == "Control",])
MicroDiab_Denmark_T2D_Variants <- rownames(MicroDiab_Denmark_Metadata[MicroDiab_Denmark_Metadata[,3] != "Control",])
MicroDiab_Denmark_Summary_Statistics <- data.frame("Shannon"=diversity(MicroDiab_Denmark_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(MicroDiab_Denmark_et_al_species),"Dysbiosis_Score"=dysbiosis_score(MicroDiab_Denmark_et_al_species,MicroDiab_Denmark_Control),"GMHI"=GMHI(MicroDiab_Denmark_et_al_species),"hack_top_17_score"=hack_top_17(MicroDiab_Denmark_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(MicroDiab_Denmark_et_al_species),"Seq_Type"="16S","Cohort_Type"="I","Study_Name"="MicroDiab_Denmark")

wilcox_MicroDiab_Denmark <- wilcox_batch(t(MicroDiab_Denmark_Summary_Statistics[MicroDiab_Denmark_T2D_Variants,1:6]),t(MicroDiab_Denmark_Summary_Statistics[MicroDiab_Denmark_Control,1:6]))

print("Saleem")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Saleem_species.RData")
load("Saleem_species.RData")
Saleem_Control <- rownames(Saleem_et_al_metadata[Saleem_et_al_metadata[,3]=="control",])
Saleem_T2D <- rownames(Saleem_et_al_metadata[Saleem_et_al_metadata[,3]!="control",])
Saleem_Summary_Statistics <- data.frame("Shannon"=diversity(Saleem_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Saleem_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Saleem_et_al_species,Saleem_Control),"GMHI"=GMHI(Saleem_et_al_species),"hack_top_17_score"=hack_top_17(Saleem_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Saleem_et_al_species),"Seq_Type"="16S","Cohort_Type"="RM","Study_Name"="Saleem")

wilcox_Saleem <- wilcox_batch(t(Saleem_Summary_Statistics[Saleem_T2D,1:6]),t(Saleem_Summary_Statistics[Saleem_Control,1:6]))

print("Parbie")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Parbie_species.RData")
load("Parbie_species.RData")
Parbie_et_al_species <- Parbie_et_al_species/rowSums(Parbie_et_al_species)
Parbie_Control <- rownames(Parbie_et_al_metadata[Parbie_et_al_metadata$disease=="control",])
Parbie_HIV <- rownames(Parbie_et_al_metadata[Parbie_et_al_metadata$disease!="control",])
Parbie_Summary_Statistics <- data.frame("Shannon"=diversity(Parbie_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Parbie_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Parbie_et_al_species,Parbie_Control),"GMHI"=GMHI(Parbie_et_al_species),"hack_top_17_score"=hack_top_17(Parbie_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Parbie_et_al_species),"Seq_Type"="16S","Cohort_Type"="RM","Study_Name"="Parbie")

wilcox_Parbie <- wilcox_batch(t(Parbie_Summary_Statistics[Parbie_HIV,1:6]),t(Parbie_Summary_Statistics[Parbie_Control,1:6]))

print("Song")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Song_species.RData")
load("Song_species.RData")
Song_Control <- rownames(Song_et_al_metadata[Song_et_al_metadata[,1]=="Normotension",])
Song_HT <- rownames(Song_et_al_metadata[Song_et_al_metadata[,1]!="Normotension",])
Song_Summary_Statistics <- data.frame("Shannon"=diversity(Song_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Song_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Song_et_al_species,Song_Control),"GMHI"=GMHI(Song_et_al_species),"hack_top_17_score"=hack_top_17(Song_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Song_et_al_species),"Seq_Type"="16S","Cohort_Type"="I","Study_Name"="Song")

wilcox_Song <- wilcox_batch(t(Song_Summary_Statistics[Song_HT,1:6]),t(Song_Summary_Statistics[Song_Control,1:6]))

print("Xu")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Xu_species.RData")
load("Xu_species.RData")
Xu_Control <- grep("NA",rownames(Xu_et_al_metadata[((Xu_et_al_metadata$MMSE >= 27)&(Xu_et_al_metadata$MMSE >= 27)),]),value=TRUE,invert=TRUE)
Xu_Unhealthy_Elderly <- grep("NA",rownames(Xu_et_al_metadata[!((Xu_et_al_metadata$MMSE >= 27)&(Xu_et_al_metadata$Barthel_Score == 100)),]),value=TRUE,invert=TRUE)

Xu_Summary_Statistics <- data.frame("Shannon"=diversity(Xu_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Xu_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Xu_et_al_species,Xu_Control),"GMHI"=GMHI(Xu_et_al_species),"hack_top_17_score"=hack_top_17(Xu_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Xu_et_al_species),"Seq_Type"="WGS","Cohort_Type"="I","Study_Name"="Xu")

wilcox_Xu <- wilcox_batch(t(Xu_Summary_Statistics[Xu_Unhealthy_Elderly,1:6]),t(Xu_Summary_Statistics[Xu_Control,1:6]))

print("Hernandez")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Hernandez_species.RData")
load("Hernandez_species.RData")
Hernandez_Control <- rownames(Hernandez_et_al_metadata[Hernandez_et_al_metadata$diabetes_status != "T2D",])
Hernandez_T2D <- rownames(Hernandez_et_al_metadata[Hernandez_et_al_metadata$diabetes_status == "T2D",])
Hernandez_Summary_Statistics <- data.frame("Shannon"=diversity(Hernandez_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Hernandez_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Hernandez_et_al_species,Hernandez_Control),"GMHI"=GMHI(Hernandez_et_al_species),"hack_top_17_score"=hack_top_17(Hernandez_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Hernandez_et_al_species),"Seq_Type"="16S","Cohort_Type"="RM","Study_Name"="Hernandez")

wilcox_Hernandez <- wilcox_batch(t(Hernandez_Summary_Statistics[Hernandez_T2D,1:6]),t(Hernandez_Summary_Statistics[Hernandez_Control,1:6]))

print("Wallen")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Wallen_species.RData")
load("Wallen_species.RData")
Wallen_Parkinsons <- rownames(Wallen_et_al_metadata[Wallen_et_al_metadata$health_state=="Parkinsons",])
Wallen_Control <- rownames(Wallen_et_al_metadata[Wallen_et_al_metadata$health_state!="Parkinsons",])

Wallen_Summary_Statistics <- data.frame("Shannon"=diversity(Wallen_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Wallen_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Wallen_et_al_species,Wallen_Control),"GMHI"=GMHI(Wallen_et_al_species),"hack_top_17_score"=hack_top_17(Wallen_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Wallen_et_al_species),"Seq_Type"="16S","Cohort_Type"="I","Study_Name"="Wallen")

wilcox_Wallen <- wilcox_batch(t(Wallen_Summary_Statistics[Wallen_Parkinsons,1:6]),t(Wallen_Summary_Statistics[Wallen_Control,1:6]))

print("Bajaj")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Bajaj_species.RData")
load("Bajaj_species.RData")
Bajaj_Control <- rownames(Bajaj_et_al_metadata)[Bajaj_et_al_metadata[,1]=="Controls"]
Bajaj_CD_ITB <- rownames(Bajaj_et_al_metadata)[Bajaj_et_al_metadata[,1]=="ITB"]

Bajaj_Summary_Statistics <- data.frame("Shannon"=diversity(Bajaj_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Bajaj_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Bajaj_et_al_species,Bajaj_Control),"GMHI"=GMHI(Bajaj_et_al_species),"hack_top_17_score"=hack_top_17(Bajaj_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Bajaj_et_al_species),"Seq_Type"="16S","Cohort_Type"="RM","Study_Name"="Bajaj")

wilcox_Bajaj <- wilcox_batch(t(Bajaj_Summary_Statistics[Bajaj_CD_ITB,1:6]),t(Bajaj_Summary_Statistics[Bajaj_Control,1:6]))

print("Pang")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Pang_species.RData")
load("Pang_species.RData")
Pang_et_al_species <- Pang_et_al_species/rowSums(Pang_et_al_species)

Pang_Control <- rownames(Pang_CrossSectional_Metadata)[(!is.na(Pang_CrossSectional_Metadata$Health_status))&(Pang_CrossSectional_Metadata$Health_status == "H")]
Pang_Unhealthy_Elderly <- rownames(Pang_CrossSectional_Metadata)[(!is.na(Pang_CrossSectional_Metadata$Health_status))&(Pang_CrossSectional_Metadata$Health_status != "H")]

Pang_Summary_Statistics <- data.frame("Shannon"=diversity(Pang_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Pang_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Pang_et_al_species,Pang_Control),"GMHI"=GMHI(Pang_et_al_species),"hack_top_17_score"=hack_top_17(Pang_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Pang_et_al_species),"Seq_Type"="16S","Cohort_Type"="I","Study_Name"="Pang")

wilcox_Pang <- wilcox_batch(t(Pang_Summary_Statistics[Pang_Unhealthy_Elderly,1:6]),t(Pang_Summary_Statistics[Pang_Control,1:6]))

print("Luan")
#load("G:\\My Drive\\Lab\\Projects\\CoreFinder\\AdditionalValidation\\Luan_species.RData")
load("Luan_species.RData")
Luan_GDS <- rownames(Luan_et_al_metadata)[Luan_et_al_metadata$GDS_15>=8]
Luan_Control <- rownames(Luan_et_al_metadata)[Luan_et_al_metadata$GDS_15<8]

Luan_Summary_Statistics <- data.frame("Shannon"=diversity(Luan_et_al_species),"Kendall_Uniqueness"=kendall_uniqueness(Luan_et_al_species),"Dysbiosis_Score"=dysbiosis_score(Luan_et_al_species,Luan_Control),"GMHI"=GMHI(Luan_et_al_species),"hack_top_17_score"=hack_top_17(Luan_et_al_species),"GMHI_hack_top_17"=GMHI_hack_top_17(Luan_et_al_species),"Seq_Type"="WGS","Cohort_Type"="I","Study_Name"="Luan")

wilcox_Luan <-  wilcox_batch(t(Luan_Summary_Statistics[Luan_GDS,1:6]),t(Luan_Summary_Statistics[Luan_Control,1:6]))


Directionality_Disease_Validation <- as.data.frame(matrix(NA,12,6))
rownames(Directionality_Disease_Validation) <- c("Flemer","Nagata","MicroDiab_Denmark","Saleem","Parbie","Song","Xu","Hernandez","Wallen","Bajaj","Pang","Luan")
colnames(Directionality_Disease_Validation) <- c("Shannon","Kendall_Uniqueness","Dysbiosis_Score","GMWI","hack_top_17","GMWI_hack_top_17")

Directionality_Disease_Validation["Flemer",] <- ifelse(wilcox_Flemer[,1]<=0.05,2*wilcox_Flemer[,2],wilcox_Flemer[,2])
Directionality_Disease_Validation["Nagata",] <- ifelse(wilcox_Nagata[,1]<=0.05,2*wilcox_Nagata[,2],wilcox_Nagata[,2])
Directionality_Disease_Validation["MicroDiab_Denmark",] <- ifelse(wilcox_MicroDiab_Denmark[,1]<=0.05,2*wilcox_MicroDiab_Denmark[,2],wilcox_MicroDiab_Denmark[,2])
Directionality_Disease_Validation["Saleem",] <- ifelse(wilcox_Saleem[,1]<=0.05,2*wilcox_Saleem[,2],wilcox_Saleem[,2])
Directionality_Disease_Validation["Parbie",] <- ifelse(wilcox_Parbie[,1]<=0.05,2*wilcox_Parbie[,2],wilcox_Parbie[,2])
Directionality_Disease_Validation["Song",] <- ifelse(wilcox_Song[,1]<=0.05,2*wilcox_Song[,2],wilcox_Song[,2])
Directionality_Disease_Validation["Xu",] <- ifelse(wilcox_Xu[,1]<=0.05,2*wilcox_Xu[,2],wilcox_Xu[,2])
Directionality_Disease_Validation["Hernandez",] <- ifelse(wilcox_Hernandez[,1]<=0.05,2*wilcox_Hernandez[,2],wilcox_Hernandez[,2])
Directionality_Disease_Validation["Wallen",] <- ifelse(wilcox_Wallen[,1]<=0.05,2*wilcox_Wallen[,2],wilcox_Wallen[,2])
Directionality_Disease_Validation["Bajaj",] <- ifelse(wilcox_Bajaj[,1]<=0.05,2*wilcox_Bajaj[,2],wilcox_Bajaj[,2])
Directionality_Disease_Validation["Pang",] <- ifelse(wilcox_Pang[,1]<=0.05,2*wilcox_Pang[,2],wilcox_Pang[,2])
Directionality_Disease_Validation["Luan",] <- ifelse(wilcox_Luan[,1]<=0.05,2*wilcox_Luan[,2],wilcox_Luan[,2])

Combined_Summary_Statistics <- as.data.frame(rbind(Flemer_Summary_Statistics[,c(1:8)],Nagata_Summary_Statistics[,c(1:8)],MicroDiab_Denmark_Summary_Statistics[,c(1:8)],Saleem_Summary_Statistics[,c(1:8)],Parbie_Summary_Statistics[,c(1:8)],Song_Summary_Statistics[,c(1:8)],Xu_Summary_Statistics[,c(1:8)],Hernandez_Summary_Statistics[,c(1:8)],Wallen_Summary_Statistics[,c(1:8)],Bajaj_Summary_Statistics[,c(1:8)],Pang_Summary_Statistics[,c(1:8)],Luan_Summary_Statistics[,c(1:8)]))

Combined_Control <- c(Flemer_Control,Nagata_Control,MicroDiab_Denmark_Control,Saleem_Control,Parbie_Control,Song_Control,Xu_Control,Hernandez_Control,Wallen_Control,Bajaj_Control,Pang_Control,Luan_Control)

Combined_Summary_Statistics$Subject_Type <- "Diseased"
Combined_Summary_Statistics[Combined_Control,"Subject_Type"] <- "Control"

Combined_Summary_Statistics$Study_Name <- NA
Combined_Summary_Statistics[rownames(Flemer_Summary_Statistics),"Study_Name"] <- "Flemer"
Combined_Summary_Statistics[rownames(Nagata_Summary_Statistics),"Study_Name"] <- "Nagata"
Combined_Summary_Statistics[rownames(MicroDiab_Denmark_Summary_Statistics),"Study_Name"] <- "MicroDiab_Denmark"
Combined_Summary_Statistics[rownames(Saleem_Summary_Statistics),"Study_Name"] <- "Saleem"
Combined_Summary_Statistics[rownames(Parbie_Summary_Statistics),"Study_Name"] <- "Parbie"
Combined_Summary_Statistics[rownames(Song_Summary_Statistics),"Study_Name"] <- "Song"
Combined_Summary_Statistics[rownames(Xu_Summary_Statistics),"Study_Name"] <- "Xu"
Combined_Summary_Statistics[rownames(Xu_Summary_Statistics),"Study_Name"] <- "Xu"
Combined_Summary_Statistics[rownames(Hernandez_Summary_Statistics),"Study_Name"] <- "Hernandez"
Combined_Summary_Statistics[rownames(Wallen_Summary_Statistics),"Study_Name"] <- "Wallen"
Combined_Summary_Statistics[rownames(Bajaj_Summary_Statistics),"Study_Name"] <- "Bajaj"
Combined_Summary_Statistics[rownames(Pang_Summary_Statistics),"Study_Name"] <- "Pang"
Combined_Summary_Statistics[rownames(Luan_Summary_Statistics),"Study_Name"] <- "Luan"

df_Confounding <- as.data.frame(matrix(NA,5,3))
rownames(df_Confounding) <- c("Shannon","Kendall_Uniqueness","GMHI","Dysbiosis_Score","hack_top_17_score")
colnames(df_Confounding) <- c("Study_Name","Cohort_Type","Seq_Type")

df_Confounding["Shannon","Study_Name"] <- ifelse(summary(lm(Shannon~as.numeric(as.factor(Study_Name)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Shannon","Cohort_Type"] <- ifelse(summary(lm(Shannon~as.numeric(as.factor(Cohort_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Shannon","Seq_Type"] <- ifelse(summary(lm(Shannon~as.numeric(as.factor(Seq_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Kendall_Uniqueness","Study_Name"] <- ifelse(summary(lm(Kendall_Uniqueness~as.numeric(as.factor(Study_Name)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Kendall_Uniqueness","Cohort_Type"] <- ifelse(summary(lm(Kendall_Uniqueness~as.numeric(as.factor(Cohort_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Kendall_Uniqueness","Seq_Type"] <- ifelse(summary(lm(Kendall_Uniqueness~as.numeric(as.factor(Seq_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["GMHI","Study_Name"] <- ifelse(summary(lm(GMHI~as.numeric(as.factor(Study_Name)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["GMHI","Cohort_Type"] <- ifelse(summary(lm(GMHI~as.numeric(as.factor(Cohort_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["GMHI","Seq_Type"] <- ifelse(summary(lm(GMHI~as.numeric(as.factor(Seq_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Dysbiosis_Score","Study_Name"] <- ifelse(summary(lm(Dysbiosis_Score~as.numeric(as.factor(Study_Name)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Dysbiosis_Score","Cohort_Type"] <- ifelse(summary(lm(Dysbiosis_Score~as.numeric(as.factor(Cohort_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["Dysbiosis_Score","Seq_Type"] <- ifelse(summary(lm(Dysbiosis_Score~as.numeric(as.factor(Seq_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["hack_top_17_score","Study_Name"] <- ifelse(summary(lm(hack_top_17_score~as.numeric(as.factor(Study_Name)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["hack_top_17_score","Cohort_Type"] <- ifelse(summary(lm(hack_top_17_score~as.numeric(as.factor(Cohort_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

df_Confounding["hack_top_17_score","Seq_Type"] <- ifelse(summary(lm(hack_top_17_score~as.numeric(as.factor(Seq_Type)),data=Combined_Summary_Statistics[Combined_Summary_Statistics$Subject_Type == "Control",]))$coefficients[2,4]<=0.05,1,0)

temp0 <- merge(t(Flemer_et_al_species),t(Nagata_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(t(Flemer_et_al_species),t(Nagata_et_al_species),by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(t(temp0),t(MicroDiab_Denmark_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(t(temp0),t(MicroDiab_Denmark_et_al_species),by="row.names",all=TRUE)[,1]
temp1 <- apply(temp1,1,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(t(temp1),t(Saleem_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(t(temp1),t(Saleem_et_al_species),by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,1,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(t(temp2),t(Parbie_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(t(temp2),t(Parbie_et_al_species),by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,1,function(x)(ifelse(is.na(x),0,x)))

temp4 <- merge(t(temp3),t(Song_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp4) <- merge(t(temp3),t(Song_et_al_species),by="row.names",all=TRUE)[,1]
temp4 <- apply(temp4,1,function(x)(ifelse(is.na(x),0,x)))

temp5 <- merge(t(temp4),t(Xu_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp5) <- merge(t(temp4),t(Xu_et_al_species),by="row.names",all=TRUE)[,1]
temp5 <- apply(temp5,1,function(x)(ifelse(is.na(x),0,x)))

temp6 <- merge(t(temp5),t(Hernandez_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp6) <- merge(t(temp5),t(Hernandez_et_al_species),by="row.names",all=TRUE)[,1]
temp6 <- apply(temp6,1,function(x)(ifelse(is.na(x),0,x)))

temp7 <- merge(t(temp6),t(Wallen_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp7) <- merge(t(temp6),t(Wallen_et_al_species),by="row.names",all=TRUE)[,1]
temp7 <- apply(temp7,1,function(x)(ifelse(is.na(x),0,x)))

temp8 <- merge(t(temp7),t(Bajaj_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp8) <- merge(t(temp7),t(Bajaj_et_al_species),by="row.names",all=TRUE)[,1]
temp8 <- apply(temp8,1,function(x)(ifelse(is.na(x),0,x)))

temp9 <- merge(t(temp8),t(Pang_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp9) <- merge(t(temp8),t(Pang_et_al_species),by="row.names",all=TRUE)[,1]
temp9 <- apply(temp9,1,function(x)(ifelse(is.na(x),0,x)))

temp10 <- merge(t(temp9),t(Luan_et_al_species),by="row.names",all=TRUE)[,-1]
rownames(temp10) <- merge(t(temp9),t(Luan_et_al_species),by="row.names",all=TRUE)[,1]
temp10 <- apply(temp10,1,function(x)(ifelse(is.na(x),0,x)))

Combined_Species <- as.data.frame(temp10)

Study_Names <- unique(Combined_Summary_Statistics$Study_Name)

df_LODO <- as.data.frame(matrix(NA,length(Study_Names),5))
rownames(df_LODO) <- Study_Names
colnames(df_LODO) <- c("Shannon","Kendall_Uniqueness","GMHI","Dysbiosis_Score","hack_top_17_score")

for(i in 1:length(Study_Names))
{
	Study_Name <- Study_Names[i]
	print(Study_Name)
	print("Generating LODO Summary Statistics")
	reference <- rownames(Combined_Summary_Statistics[(Combined_Summary_Statistics$Study_Name != Study_Name)&(Combined_Summary_Statistics$Subject_Type == "Control"),])
	diseased <- rownames(Combined_Summary_Statistics[(Combined_Summary_Statistics$Study_Name == Study_Name)&(Combined_Summary_Statistics$Subject_Type != "Control"),])
	temp_species <- Combined_Species[c(reference,diseased),]
	temp_species <- temp_species[,colSums(temp_species)>0]
	temp_summary_statistics <- data.frame("Shannon"=diversity(temp_species),"Kendall_Uniqueness"=kendall_uniqueness(temp_species),"Dysbiosis_Score"=dysbiosis_score(temp_species,reference),"GMHI"=GMHI(temp_species),"hack_top_17"=hack_top_17(temp_species))
}


all_objects <- ls()
temp_objects <- grep("^temp", all_objects, value = TRUE)
rm(list = temp_objects)

save.image("241030_Comparative_Evaluation.RData")