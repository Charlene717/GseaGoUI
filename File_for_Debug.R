##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)


##### Load Packages  #####
source("FUN_Package_InstLoad.R")
FUN_Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT","shiny","shinyFiles")
FUN_BiocManager.set <- c("clusterProfiler","enrichplot","pathview")
## Set the desired organism
# organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
# c(organism,"fgsea")

FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
rm(FUN_Basic.set, FUN_BiocManager.set)

##### Function setting #####
## Call function
source("FUN_GSEA_ANAL.R"); source("FUN_GSEA_LargeGeneSet.R"); source("FUN_GSEA_ggplot.R"); source("FUN_GSEA_ForOFFL.R");
source("FUN_Find_Markers.R"); source("FUN_DEG_Analysis.R"); source("FUN_VolcanoPlot.R");
source("FUN_ggPlot_vline.R"); source("FUN_DistrPlot.R"); source("FUN_Group_GE.R"); source("FUN_Beautify_ggplot.R");



#### GeneExp.df ####
GeneExp.df <- read.csv("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/Xena_TCGA_LGG_GE.tsv",
                       sep = "\t", row.names = 1, check.names = F)
colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
GeneExp.df <- data.frame(Gene = row.names(GeneExp.df),GeneExp.df)

# GeneExpS.df <- read.csv("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/Xena_TCGA_LGG_GE_S.tsv",
#                        sep = "\t", row.names = 1, check.names = F)

GeneExpS.df <- GeneExp.df[sample(1:nrow(GeneExp.df),5000),]


#### Anno.df ####
Anno.df <- read.csv("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/TCGA.LGG.sampleMap_LGG_clinicalMatrix.tsv",
                    sep = "\t", row.names = 1, check.names = F)

Anno.Set <- c("X_INTEGRATION", "X_PATIENT","sample_type","gender","X_primary_disease","X_primary_site",
              "eczema_history","family_history_of_cancer","first_presenting_symptom",
              "followup_treatment_success","hay_fever_history","headache_history","histological_type","history_ionizing_rt_to_head",
              "history_of_neoadjuvant_treatment", "neoplasm_histologic_grade",
              "seizure_history", "supratentorial_localization","tissue_source_site","age_at_initial_pathologic_diagnosis")

Anno.df <- Anno.df[,Anno.Set]
Anno.df <- data.frame(ID = row.names(Anno.df),Anno.df)


###################################################
#### VolcanoPlot ####
  DEG_ANAL.lt <-
    FUN_DEG_Analysis(GeneExp.df, Anno.df,
                     GroupType = "sample_type" ,GroupCompare = c("Primary Tumor","Recurrent Tumor"),
                     ThrSet = list(LogFC = c("logFC", 1),
                                   pVal = c("PValue",0.05)),
                     TarGeneName = "ATRX", GroupMode = Mode_Group, SampleID = "X_INTEGRATION",
                     Save.Path = Save.Path, SampleName = SampleName, AnnoName = "AvB")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]
  rownames(DE_Extract.df) <- seq(1:nrow(DE_Extract.df))
  Pos.List <- DEG_ANAL.lt[["DE_Extract_FltH.set"]]
  Neg.List <- DEG_ANAL.lt[["DE_Extract_FltL.set"]]
  Plot.VolcanoPlot <- VolcanoPlot(DE_Extract.df ,
                                  Pos.List,
                                  Neg.List,
                                  log2FC = 1, PValueSet = 0.05,  ShowGeneNum = 5)
  Plot.VolcanoPlot

#### GSEA analysis ####
  GeneSet.df <- read.csv("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/2022-11-03_GSEA_Geneset_Pathway_3Database_WithoutFilter.txt",
                         sep = "\t", header = F, check.names = F) %>% as.data.frame()
  FUN_GSEA_ANAL(DE_Extract.df, pathwayGeneSet = GeneSet.df,
                TarGeneName = "ATRX", GroupMode = Mode_Group,
                ThrSet = list(LogFC = c("logFC", 1),
                              pVal = c("PValue",0.05)),
                Species = "Homo sapiens", # Speices type can check by msigdbr_species()
                Save.Path = Save.Path, SampleName = SampleName, AnnoName = "Path")
