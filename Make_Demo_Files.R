##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)


#### GeneExp.df ####
GeneExp.df <- read.csv("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/Xena_TCGA_LGG_GE.tsv",
                       sep = "\t", row.names = 1, check.names = F)
colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
GeneExp.df <- data.frame(Gene = row.names(GeneExp.df),GeneExp.df)

# GeneExpS.df <- read.csv("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/Xena_TCGA_LGG_GE_S.tsv",
#                        sep = "\t", row.names = 1, check.names = F)

GeneExpS.df <- GeneExp.df[sample(1:nrow(GeneExp.df),2000),]
write.table(GeneExpS.df, paste0("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/",Sys.Date(),"_","Xena_TCGA_LGG_GE_S.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

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

write.table(Anno.df, paste0("D:/Dropbox/##_GitHub/##_Charlene/GseaGoUI/Input_TCGA/",Sys.Date(),"_","Xena_TCGA_LGG_Anno_S.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)
