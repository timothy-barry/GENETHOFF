

# Set R options
options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))


# load necessary libraries
library(tidyverse,quietly = T, verbose = F,warn.conflicts = F)
library(GenomicRanges,quietly = T, verbose = F,warn.conflicts = F)
library(annotatr,quietly = T, verbose = F,warn.conflicts = F)
#library(biomaRt,quietly = T, verbose = F,warn.conflicts = F)
library(writexl,quietly = T, verbose = F,warn.conflicts = F)
library(rtracklayer,quietly = T, verbose = F,warn.conflicts = F)
library(pwalign,quietly = T, verbose = F,warn.conflicts = F)


# get argument from command line 
args <- commandArgs(trailingOnly = T)

# Parse arguments to objects
annotation <- args[1]
files <- args[2]
output<- args[3]


###########################################################
## debug
# files = "../debug/05-Report/epe1"
# annotation <- "GRCh38"
###########################################################

# load sample rdata file
names(files) <- str_remove(basename(files),pattern = ".rdata")
load(files)



## if file is not empty ...
if(exists("cluster_annotated") && nrow(cluster_annotated)>0){
  
  ## sort & convert some features type
  results_df <- cluster_annotated %>% 
    arrange(desc(N_UMI_cluster)) %>% 
    mutate(gRNA = as.character(grna),
           gRNA_name = grna@metadata$name,
           library = names(files)) 
  rm(cluster_annotated);
  
  
  # load pre-prepared annotation GRange object
  gtf <- readRDS(paste("../02-ressources/",annotation,".rds",sep=""))
  
  # convert OTS to gRanges -----------------------------------------------------
  results_granges <- makeGRangesFromDataFrame(results_df, 
                                              ignore.strand = T,
                                              keep.extra.columns = T,
                                              strand.field = "strand_guide",
                                              start.field = "start_cluster",
                                              end.field = "end_cluster",
                                              seqnames.field = "chromosome",
                                              na.rm = T)
  
  
  
  
  # annotate gRanges -------------------------------------------------------------------
  
  results_granges_annot <- annotate_regions(regions = results_granges,
                                            annotations = gtf,
                                            ignore.strand = T,
                                            minoverlap = 1)
  
  
  # convert to data.frame
  results_granges_df <- data.frame(results_granges_annot) %>%
    dplyr::select(clusterID,starts_with('anno') ) 
  
  
  
  
  
  # collapse annotation per cluster -----------------------------------------------------
  results_granges_df_annot <- results_granges_df %>%
    distinct(clusterID,annot.type,annot.gene_id,annot.gene_name,annot.gene_biotype)%>%
    group_by(clusterID,annot.gene_id,annot.gene_name,annot.gene_biotype) %>%
    summarise(pos =toString(annot.type)) %>% 
    mutate(position=case_when(str_detect(pos,"exon") ~ "exon", TRUE ~ "intron")) %>%
    dplyr::select(-pos) %>% 
    group_by(clusterID) %>% 
    summarise(gene_ensemblID = toString(annot.gene_id),
              Symbol = toString(annot.gene_name),
              gene_type = toString(annot.gene_biotype),
              position = toString(position))
  
  
  
  
 
  
  
  # Annotate ------------------------------------------------------------------
  
  
  
  results_granges_df_annot <- results_df %>%  left_join(results_granges_df_annot, by = "clusterID")
  
  
  
  # Annotate clusters that have a similar alignment pattern --> multihits -------------
  
  
  multihits <- results_granges_df_annot %>% 
    filter(!is.na(Alignment))  %>%    #keep clusters with gRNa match
    distinct(clusterID,Alignment)%>% 
    arrange(Alignment) %>% 
    mutate(multi_cluster = as.numeric(factor(Alignment))) %>%     # assign ID to clusters with identical alignment pattern
    group_by(multi_cluster) %>% 
    mutate(SimilarAlignmentCount = n_distinct(clusterID)) %>%     # count how many clusters share the same alignment pattern
    distinct(clusterID,multi_cluster,SimilarAlignmentCount)       # select columns for annotation
  
  
  
  results_granges_df_annot <- results_granges_df_annot %>% 
    left_join(multihits, by = "clusterID") %>% 
    group_by(multi_cluster) %>% 
    mutate(n_UMI_multiSum = sum(N_UMI_cluster),                 # calculate sum of UMI counts per ID. As we don't know which site is the correct one, we sum up all UMI from all site with same pattern
           multiHit = SimilarAlignmentCount>1)                  # if there are more than 1 site per ID, it's a multihit
  
  

  # Save to csv file  -----------------------------------------------------------------

  
  write.table(results_granges_df_annot, str_replace(output,pattern = "xlsx$","csv"),sep=",",row.names = F,quote=F)
  
  
  # save to excel file -----------------------------------------------------
  
  sp <- split(results_granges_df_annot,results_granges_df_annot$library)
  writexl::write_xlsx(sp,path = output,col_names = T,format_headers = T)
  
}
