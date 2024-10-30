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
#args <- c("debug/simulated3.cluster_slop.fa","debug/simulated3.cluster_slop.bed","debug/simulated3.collapsefragPerISCluster.bed","CTAGGTATGTTATCTGAAAGT","")

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

pam = DNAString("NNN")
pam_length = nchar(pam)

offset = 4


## Calculate pairwise alignments between insertion site sequence and gRNA sequence in forward and reverse orientation
watson= pwalign::pairwiseAlignment(pattern = fasta,subject = grna,type = "local-global",gapOpening=0,)
crick = pwalign::pairwiseAlignment(pattern = reverseComplement(fasta),subject = grna,type = "local-global",gapOpening=0,)


# find strand with best score

## NEW 
align_stat <- data.frame(position = names(fasta),
                         watson_score = pwalign::score(watson),
                         crick_score =  pwalign::score(crick),
                         watson_edit = pwalign::nedit(watson),
                         crick_edit = pwalign::nedit(crick),
                         watson_mm = pwalign::nmismatch(watson),
                         crick_mm = pwalign::nmismatch(crick),
                         watson_pid = pwalign::pid(watson,),
                         crick_pid = pwalign::pid(crick)) %>% 
  separate("position", into = c("clusterID","cluster"),sep = "::",convert = T)

## use edit instead of mismatches (it includes indels and mismatches)
align_stat <- align_stat %>%
  filter(watson_edit <= 6 | crick_edit <= 6) %>%
  mutate(grna_orientation = case_when(watson_score > crick_score ~ "watson",
                                      watson_score < crick_score  ~"crick",
                                      TRUE ~NA))

## WATSON strand analysis

watson_best <- align_stat %>% 
  filter(grna_orientation == "watson") %>% 
  select(clusterID,cluster,grna_orientation,starts_with("watson")) %>%
  rename_all(~str_remove(.,"watson_"))

fasta_temp <- fasta[watson_best$clusterID]

watson_sub <- watson[watson_best$clusterID]
watson_best$seq_gDNA <- pwalign::alignedPattern(watson_sub) %>% as.character
watson_best$seq_gRNA <- pwalign::alignedSubject(watson_sub)%>% as.character
watson_best$mismatches_position_gRNA <- lapply(watson_sub@subject@mismatch,toString) %>% unlist 
watson_best <- watson_best %>% bind_cols(data.frame(watson_sub@pattern@range))
watson_best <- watson_best %>% left_join(fasta_temp %>% data.frame %>% select(sequence_window = ".") %>% rownames_to_column('id') %>% separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))





pams <- lapply(seq_along(watson_sub), function(x) {
  start = watson_best$end[x]+1
  end = watson_best$end[x] + pam_length
  length = seqlengths(fasta[watson_best$clusterID[x]])
  if((end + pam_length)> length){
    "..."
  } else {
    subseq(fasta[watson_best$clusterID[x]],start  , end )
    
  }
}
)

pams <- do.call(c,pams)

pams_align  <- pwalign::pairwiseAlignment(pattern = pams,subject = pam,type = "local-global",gapOpening=10)

watson_best$pam_gDNA <- pwalign::alignedPattern(pams_align) %>% as.character
watson_best$pam_gRNA <- pwalign::alignedSubject(pams_align) %>% as.character
watson_best <-watson_best %>% mutate(rank=row_number())
watson_best <- watson_best %>% bind_cols(pwalign::nindel(watson_sub)@insertion %>% data.frame %>% select(Insertion_count = Length, Insertion_cum_width = WidthSum))
watson_best <- watson_best %>% bind_cols(pwalign::nindel(watson_sub)@deletion %>% data.frame %>% select(Deletion_count = Length, Deletion_cum_width = WidthSum))


indels_list <- pwalign::indel(watson_sub)

indels_table <- bind_rows(
  bind_rows(
    "insertions" = pwalign::insertion(indels_list)%>%
      data.frame,
    "deletions"= pwalign::deletion(indels_list) %>%
      data.frame, .id = "indel") %>% 
    group_by(indel,group = group) %>% 
    summarise(start = toString(start),
              end = toString(end),
              width = toString(width))) %>% 
  pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
  select(group,ends_with("deletions"),ends_with("insertions")) %>% 
  arrange(group)


watson_best <- watson_best %>% left_join(indels_table, by = c("rank"="group"))







## CRICK strand analysis

crick_best <- align_stat %>% 
  filter(grna_orientation == "crick") %>% 
  select(clusterID,cluster,grna_orientation,starts_with("crick")) %>%
  rename_all(~str_remove(.,"crick_"))

fasta_temp <- fasta[crick_best$clusterID]

crick_sub <- crick[crick_best$clusterID]
crick_best$seq_gDNA <- pwalign::alignedPattern(crick_sub) %>% as.character
crick_best$seq_gRNA <- pwalign::alignedSubject(crick_sub)%>% as.character
crick_best$mismatches_position_gRNA <- lapply(crick_sub@subject@mismatch,toString) %>% unlist 
crick_best <- crick_best %>% bind_cols(data.frame(crick_sub@pattern@range))
crick_best <- crick_best %>% left_join(fasta_temp %>% data.frame %>% select(sequence_window = ".") %>% rownames_to_column('id') %>% separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))


pams <- lapply(seq_along(crick_sub), function(x) {
  start = crick_best$end[x]+1
  end = crick_best$end[x] + pam_length
  length = seqlengths(fasta[crick_best$clusterID[x]])
  if((end + pam_length)> length){
    "..."
  } else {
    subseq(fasta[crick_best$clusterID[x]],start  , end )
    
  }
}
)

pams <- do.call(c,pams)

pams_align  <- pwalign::pairwiseAlignment(pattern = pams,subject = pam,type = "local-global",gapOpening=10)
crick_best <-crick_best %>% mutate(rank=row_number())
crick_best$pam_gDNA <- pwalign::alignedPattern(pams_align) %>% as.character
crick_best$pam_gRNA <- pwalign::alignedSubject(pams_align) %>% as.character
crick_best <- crick_best %>% bind_cols(pwalign::nindel(crick_sub)@insertion %>% data.frame %>% select(Insertion_count = Length, Insertion_cum_width = WidthSum))
crick_best <- crick_best %>% bind_cols(pwalign::nindel(crick_sub)@deletion %>% data.frame %>% select(Deletion_count = Length, Deletion_cum_width = WidthSum))


indels_list <- pwalign::indel(crick_sub)

indels_table <- bind_rows(
  bind_rows(
    "insertions" = pwalign::insertion(indels_list)%>%
      data.frame,
    "deletions"= pwalign::deletion(indels_list) %>%
      data.frame, .id = "indel") %>% 
    group_by(indel,group = group) %>% 
    summarise(start = toString(start),
              end = toString(end),
              width = toString(width))) %>% 
  pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
  select(group,ends_with("deletions"),ends_with("insertions")) %>% 
  arrange(group)


crick_best <- crick_best %>% left_join(indels_table, by = c("rank"="group"))

# collapse watson and crick
best <- bind_rows(watson_best,crick_best)

# reorder columns : 
best <- best %>% 
  select(clusterID, cluster, sequence_window,
         grna_orientation,
         seq_gDNA, seq_gRNA,
         alignment_start_gDNA = start, alignment_end_gDNA = end, alignment_width_gDNA = width,
         alignment_score=score,
         Identity_pct = pid,
         N_edits = edit,
         N_mismatches = mm, mismatches_position_gRNA, 
         Insertion_count,Insertion_cum_width,start_insertions,end_insertions,width_insertions,
         Deletion_count,Deletion_cum_width,start_deletions,end_deletions,width_deletions,
         pam_gDNA,pam_gRNA )



## get cuting site based on CAS offset and strand orientation
  
  
  if(nrow(best)>0){
    best <- best %>% 
      separate(cluster, into = c("chromosome","start_chr","end_chr"), sep = "[:]+|-",convert = T)
    
    best <- best %>%
      mutate(cut_gRNa_alignment = case_when(grna_orientation == "watson" ~ start_chr + alignment_end_gDNA - offset -1 ,
                             TRUE ~ end_chr - alignment_width_gDNA - alignment_start_gDNA + offset + 2))
  } else{
  
  
  col.idx <- c("start" ,"end", "width", "clusterID", "chromosome", "start_chr",  "end_chr","sequence", "mismatches", "strand_guide","cut_gRNa_alignment")
  best <- data.frame(matrix(nrow = 0, ncol = length(col.idx) ))
  names(best) <- col.idx 
  }



## annotate clusters with identified gRNA sequence
cluster_annnotated <- clusters %>% 
  left_join(best %>% 
              select(-start_chr,-end_chr,-alignment_width_gDNA,-alignment_start_gDNA,-alignment_end_gDNA) ,by = c("chromosome","clusterID"))


IS_in_clusters = bed_collapsedUMI %>% 
  left_join(cluster_annnotated %>%  
              select(clusterID,cut_gRNa_alignment), 
            by = c("clusterID")) %>%
  filter(!is.na(cut_gRNa_alignment))

IS_in_clusters = IS_in_clusters %>% 
  mutate(relative_distance = start_IS - cut_gRNa_alignment)



modal_cut_position <- IS_in_clusters %>% 
  group_by(clusterID,start_IS) %>%
  summarise(count_UMI = sum(N_UMI_IS)) %>% 
  mutate(UMI_proportion = round(count_UMI / sum(count_UMI)*100,digits = 1)) %>% 
  group_by(clusterID) %>% 
  slice_max(n = 1, order_by = UMI_proportion) %>%
  group_by(clusterID,UMI_proportion) %>% 
  summarise(cut_modal_position = toString(start_IS))

cluster_annnotated <- cluster_annnotated %>%  left_join(modal_cut_position, by = c("clusterID"))


save(cluster_annnotated, fasta,grna, bed_collapsedUMI, clusters, watson,crick,IS_in_clusters, file = args[6])

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



# test 
# library(ggmsa)
# library(msa)
# 
# x <- DNAStringSet(c(as.character(grna),best$seq_gDNA),use.names = T)
# names(x) <- c("gRNA",best$clusterID)
# 
# 
# msa_x <- msa(x,order = "input")
# tidy_msa_data <- tidy_msa(msa_x@unmasked)
# 
# ggplot(tidy_msa_data, aes(position,name,fill=character,)) + geom_tile(col = "black") + theme_bw()
# 
# 
# 
# 
# 
# 
# ggmsa(x, char_width = 0.5, seq_name = T,color= "Zappo_NT")

