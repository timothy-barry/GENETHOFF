###  Format onco gene list from oncoKG

# The output must be a ";" separated table with 2 columns only and named as below: 
# - ensembl.transcriptID
# - Onco_annotation


library(tidyverse)

update <- format.Date(x = Sys.Date(),format = "%Y-%m-%d")
genome <- "GRCh38"


## Example for GRCh38 using oncoKB database
onco_list <- read.delim("https://www.oncokb.org/api/v1/utils/cancerGeneList.txt")

onco_list_output <- onco_list %>% 
  filter(GRCh38.Isoform!="",X..of.occurrence.within.resources..Column.J.P.>1) %>%
  rename("GRCh38.Isoform"="ensembl.transcriptID",
         "Is.Tumor.Suppressor.Gene"= "TSG",
         "Is.Oncogene" = "Oncogene") %>% 
  select(ensembl.transcriptID,TSG,Oncogene) %>% 
  pivot_longer(cols = -ensembl.transcriptID) %>% 
  filter(value == "Yes") %>% 
  group_by(ensembl.transcriptID) %>% 
  summarise(Onco_annotation = paste(unique(name),collapse = "|"))

write.table(onco_list_output, file = paste("02-ressources/OncoList_OncoKB_",genome,"_",update,".tsv",sep = ""), sep=";")

