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
#args <- c("EBS_Cas_NGG_hDMD.cluster_slop.fa","EBS_Cas_NGG_hDMD.cluster_slop.bed","EBS_Cas_NGG_hDMD.collapsefragPerISCluster.bed","TCTTCCGGAACAAAGTTGCT","ebs","NNN",4)

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
gRNA_name <- args[5]

grna@metadata$name <- gRNA_name

pam <- DNAString(args[6])
pam_length <- nchar(pam)

offset <- as.numeric(args[7])

## Calculate pairwise alignments between insertion site sequence and gRNA sequence in forward and reverse orientation
watson= pairwiseAlignment(pattern = fasta,subject = grna,type = "local-global",gapOpening=0,)
crick = pairwiseAlignment(pattern = reverseComplement(fasta),subject = grna,type = "local-global",gapOpening=0,)


# find strand with best score

## NEW 
align_stat <- data.frame(position = names(fasta),
                         watson_score = score(watson),
                         crick_score =  score(crick),
                         watson_edit = nedit(watson),
                         crick_edit = nedit(crick),
                         watson_mm = nmismatch(watson),
                         crick_mm = nmismatch(crick),
                         watson_pid = pid(watson,),
                         crick_pid = pid(crick)) %>% 
  separate("position", into = c("clusterID","cluster"),sep = "::",convert = T)

## use edit instead of mismatches (it includes indels and mismatches)
align_stat <- align_stat %>%
  filter(watson_edit <= 6 | crick_edit <= 6) %>%
  mutate(grna_orientation = case_when(watson_score > crick_score ~ "watson",
                                      watson_score < crick_score  ~ "crick",
                                      TRUE ~ NA))

if(nrow(align_stat)>0){
  ## WATSON strand analysis
  
  watson_best <- align_stat %>% 
    filter(grna_orientation == "watson") %>% 
    select(clusterID,cluster,grna_orientation,starts_with("watson")) %>%
    rename_all(~str_remove(.,"watson_"))
  
  if(nrow(watson_best)>0){
    fasta_temp <- fasta[watson_best$clusterID]
    
    watson_sub <- watson[watson_best$clusterID]
    watson_best$seq_gDNA <- alignedPattern(watson_sub) %>% as.character
    watson_best$seq_gRNA <- alignedSubject(watson_sub)%>% as.character
    watson_best$mismatches_position_gRNA <- lapply(watson_sub@subject@mismatch,toString) %>% unlist 
    watson_best <- watson_best %>% bind_cols(data.frame(watson_sub@pattern@range))
    watson_best <- watson_best %>% 
      left_join(fasta_temp %>%
                  data.frame %>% select(sequence_window = ".") %>%
                  rownames_to_column('id') %>% separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))
    
    
    
    
    
    pams <- lapply(seq_along(watson_sub), function(x) {
      start = watson_best$end[x]+1
      end = watson_best$end[x] + pam_length
      length = seqlengths(fasta[watson_best$clusterID[x]])
      if((end + pam_length)> length){
        DNAStringSet("...")
      } else {
        subseq(fasta[watson_best$clusterID[x]],start  , end )
        
      }
    }
    )
    
    pams <- do.call(c,pams)
    
    pams_align  <- pairwiseAlignment(pattern = pams,subject = pam,type = "local-global",gapOpening=10)
    
    watson_best$pam_gDNA <- alignedPattern(pams_align) %>% as.character
    watson_best$pam_gRNA <- alignedSubject(pams_align) %>% as.character
    watson_best <-watson_best %>% mutate(rank=row_number())
    watson_best <- watson_best %>% bind_cols(nindel(watson_sub)@insertion %>% data.frame %>% select(Insertion_count = Length, Insertion_cum_width = WidthSum))
    watson_best <- watson_best %>% bind_cols(nindel(watson_sub)@deletion %>% data.frame %>% select(Deletion_count = Length, Deletion_cum_width = WidthSum))
    
    
    indels_list <- indel(watson_sub)
    
    indels_table <- bind_rows(
      bind_rows(
        "insertions" = insertion(indels_list)%>%
          data.frame,
        "deletions"= deletion(indels_list) %>%
          data.frame, .id = "indel") %>% 
        group_by(indel,group = group) %>% 
        summarise(start = toString(start),
                  end = toString(end),
                  width = toString(width))) %>% 
      pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
      select(group,ends_with("deletions"),ends_with("insertions")) %>% 
      arrange(group)
    
    
    watson_best <- watson_best %>% left_join(indels_table, by = c("rank"="group"))
    
  } 
  
  
  
  
  
  ## CRICK strand analysis
  
  crick_best <- align_stat %>% 
    filter(grna_orientation == "crick") %>% 
    select(clusterID,cluster,grna_orientation,starts_with("crick")) %>%
    rename_all(~str_remove(.,"crick_"))
  
  if(nrow(crick_best)>0) {
    fasta_temp <- fasta[crick_best$clusterID]
    
    crick_sub <- crick[crick_best$clusterID]
    crick_best$seq_gDNA <- alignedPattern(crick_sub) %>% as.character
    crick_best$seq_gRNA <- alignedSubject(crick_sub)%>% as.character
    crick_best$mismatches_position_gRNA <- lapply(crick_sub@subject@mismatch,toString) %>% unlist 
    crick_best <- crick_best %>% bind_cols(data.frame(crick_sub@pattern@range))
    crick_best <- crick_best %>% left_join(fasta_temp %>% data.frame %>% select(sequence_window = ".") %>% rownames_to_column('id') %>% separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))
    
    
    pams <- lapply(seq_along(crick_sub), function(x) {
      start = crick_best$end[x]+1
      end = crick_best$end[x] + pam_length
      length = seqlengths(fasta[crick_best$clusterID[x]])
      if((end + pam_length)> length){
        DNAStringSet("...")
      } else {
        subseq(fasta[crick_best$clusterID[x]],start  , end )
        
      }
    }
    )
    
    pams <- do.call(c,pams)
    
    pams_align  <- pairwiseAlignment(pattern = pam,subject = pam,type = "local-global",gapOpening=10)
    crick_best <-crick_best %>% mutate(rank=row_number())
    crick_best$pam_gDNA <- alignedPattern(pams_align) %>% as.character
    crick_best$pam_gRNA <- alignedSubject(pams_align) %>% as.character
    crick_best <- crick_best %>% bind_cols(nindel(crick_sub)@insertion %>% data.frame %>% select(Insertion_count = Length, Insertion_cum_width = WidthSum))
    crick_best <- crick_best %>% bind_cols(nindel(crick_sub)@deletion %>% data.frame %>% select(Deletion_count = Length, Deletion_cum_width = WidthSum))
    
    
    indels_list <- indel(crick_sub)
    
    indels_table <- bind_rows(
      bind_rows(
        "insertions" = insertion(indels_list)%>%
          data.frame,
        "deletions"= deletion(indels_list) %>%
          data.frame, .id = "indel") %>% 
        group_by(indel,group = group) %>% 
        summarise(start = toString(start),
                  end = toString(end),
                  width = toString(width))) %>% 
      pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
      select(group,ends_with("deletions"),ends_with("insertions")) %>% 
      arrange(group)
    
    
    crick_best <- crick_best %>% left_join(indels_table, by = c("rank"="group"))
  }
  if(nrow(watson_best)>0 & nrow(crick_best)>0){
    # collapse watson and crick
    best <- bind_rows(watson_best,crick_best)
  } else {
    if(nrow(watson_best)>0){
      best <- watson_best
    } else {
      best <- crick_best
    }
  }
  
  
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
  
}

save(cluster_annnotated, fasta,grna, bed_collapsedUMI, clusters, watson,crick,IS_in_clusters, file = args[8])


