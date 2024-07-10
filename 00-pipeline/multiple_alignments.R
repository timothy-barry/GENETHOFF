####-------------------------------------------------------###
####-------------------------------------------------------###
####-------------------------------------------------------###
####-------------------------------------------------------###
####-------------------------------------------------------###
####-------------------------------------------------------###
####-------------------------------------------------------###


args <- commandArgs(trailingOnly = T)


# Load libraries --------------------------------------------------------------------

library(tidyverse,quietly = T, verbose = F)
library(Biostrings,quietly = T,verbose = F)

# for debugging
#args <- c("epe10.5_Cas_NGG_K562.cluster_slop.fa","epe10.5_Cas_NGG_K562.cluster_slop.bed","epe10.5_Cas_NGG_K562.collapsefragPerISCluster.bed","GTTTGCCTTGTCAAGGCTAT")

# Load sample results ---------------------------------------------------------------


fasta <- readDNAStringSet(args[1],
                          use.names = T)

clusters <- read.delim(args[2],
                       header=F, 
                       col.names = c("chromosome",
                                     "start_cluster",
                                     "end_cluster",
                                     "clusterID",
                                     "N_orientations_cluster",
                                     "medianMAPQ_cluster",
                                     "N_IS_cluster",
                                     "N_UMI_cluster",
                                     "N_reads_cluster"))

bed_collapsedUMI <- read.delim(args[3],
                               header=F, 
                               col.names = c("chromosome",
                                             "start_IS",
                                             "end_IS",
                                             "IS_ID",
                                             "MedianMAPQ_IS",
                                             "strand_IS",
                                             "N_UMI_IS",
                                             "Nreads_IS",
                                             "UMI_list",
                                             "ReadPerUMI",
                                             "clusterID"))
  
grna <- DNAString(args[4])
  
## Calculate pairwise alignments between insertion site sequence and gRNA sequence in forward and reverse orientation
watson= pairwiseAlignment(pattern = fasta,subject = grna,type = "local-global",gapOpening=-100,)
crick = pairwiseAlignment(pattern = reverseComplement(fasta),subject = grna,type = "local-global",gapOpening=-100,)


df <- data.frame(Whatson = nmismatch(watson), 
                 Crick = nmismatch(crick)) %>%
  rownames_to_column("rank")

## keep orientation with the smallest number of mismatches
df  <- df %>% filter((Whatson<=6 | Crick <=6) & Whatson!=Crick)


watson_subset <- watson[as.numeric(df$rank[which(df$Whatson<df$Crick)])]
crick_subseted <- crick[as.numeric(df$rank[which(df$Whatson>df$Crick)])]

df <- data.frame(watson_subset@pattern@range) %>% 
  bind_cols(data.frame(sequence=alignedPattern(watson_subset),
                       mismatches = nmismatch(watson_subset)) %>% 
              rownames_to_column('position')) %>% 
  mutate(strand_guide = "+") %>% 
  bind_rows(data.frame(crick_subseted@pattern@range) %>% 
              bind_cols(data.frame(sequence=alignedPattern(crick_subseted),
                                   mismatches = nmismatch(crick_subseted)) %>% 
                          rownames_to_column('position')) %>% 
              mutate(strand_guide = "-"))

df <- df %>% 
  separate(position, into = c("clusterID", "chromosome","start_chr","end_chr"), sep = "[:]+|-",convert = T)

df <- df %>%
  mutate(cut = case_when(strand_guide == "+" ~ start_chr + end -3 ,
                              TRUE ~ end_chr -width -start+ 3 +2))

## annotate clusters with identified gRNA sequence
cluster_annnotated <- clusters %>% 
  left_join(df %>% 
              select(-start,-end,-width,-start_chr,-end_chr) ,by = c("chromosome","clusterID"))


IS_in_clusters = bed_collapsedUMI %>% 
  left_join(cluster_annnotated %>%  
              select(clusterID,cut), 
            by = c("clusterID")) %>%
  filter(!is.na(cut))

IS_in_clusters = IS_in_clusters %>% 
  mutate(relative_distance = start_IS - cut)
  


save(cluster_annnotated, fasta,grna, bed_collapsedUMI, clusters, watson,crick,IS_in_clusters, file = args[6])

#
# nmismatch(x)
# score(x)
# nedit(x)
# nmatch(x)
# pid(x)
# pattern(x[c(1,2,3)])
# subject(x)
# aligned(x)
# coverage(x)
# start(x)
# width(x)
