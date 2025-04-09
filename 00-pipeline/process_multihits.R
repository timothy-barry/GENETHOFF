library(tidyverse)
library(data.table)


bed <- fread("../debug/epe10.5_Cas_NGG_K562.UMI.ODN.trimmed.filtered.multi.realigned.bed",header=F)
bed <- bed %>%  mutate(V4 = str_remove(V4,"/2$"))


info <- fread("../debug/epe10.5_Cas_NGG_K562.UMI.ODN.trimmed.filtered.multi.realigned.info",header=F)
info <- info %>% mutate(AS = as.numeric(str_remove(V5,"AS:i:")),
                        FragLength = abs(V4)) %>% 
  select(AS,FragLength)


# bind bed and info -----------------------------------------------------------------

bed <- bed %>% bind_cols(info)


# Filter the bed file ---------------------------------------------------------------
# keep alignment with higher AS

bed <- bed %>% 
  group_by(V4) %>%
  slice_max(n = 1,order_by = AS,with_ties = T) 


bed2 <- bed %>% 
  ungroup %>% 
  mutate(cut_site = case_when(V6=="+"~ V2,
                                                       TRUE ~ V3)) %>% 
  separate(V4, into = c("RNAME","UMI"),sep="_") %>% 
  filter(AS> -20) %>% 
  group_by(RNAME,UMI) %>% 
  mutate(n_match = n())

bed_multi <- bed2 %>% filter(n_match>1) 

x <- bed_multi%>%  group_by(V1,cut_site,V6) %>% summarise(reads = n_distinct(RNAME), umi = n_distinct(UMI),matches = median(n_match)) %>% mutate(UMI_norm = umi/matches)


bed_single_rescued <- bed2 %>% filter(n_match == 1) 






