

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
onco_list <- args[4]


###########################################################
## debug
# files = "05-Report/Cpf1.rdata"
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
  
  
  
  ## Add oncogene / tumors suppressing gene annotation
  
  # human list examples : https://bioinfo.uth.edu/TSGene/download.cgi & https://bioinfo-minzhao.org/ongene/download.html 
  
  if(onco_list!=""){
    
    onco_list_df <- read.delim(onco_list,sep="\t") 
    
    if(all(colnames(onco_list_df) %in% c("ensembl.transcriptID","Is.Oncogene","Is.Tumor.Suppressor.Gene") & nrow(onco_list_df)>0)){
      
      ## this file must contain "ensembl.transcriptID", column "Is.Oncogene" and "Is.Tumor.Suppressor.Gene" as in https://www.oncokb.org/cancer-genes
      
      onco_list_df <- onco_list_df %>%
        distinct(ensembl.transcriptID,Is.Oncogene,Is.Tumor.Suppressor.Gene)
      
      results_granges_df <- results_granges_df %>%
        left_join(onco_list_df, by = c("annot.transcript_id" = "ensembl.transcriptID") ) %>% 
        select(-annot.transcript_id) %>%  distinct() %>% 
        replace_na(replace = list(Is.Oncogene="no",Is.Tumor.Suppressor.Gene="no"))
    } else{
      
      errorCondition("Oncogene list is empty or is not correctly formatted. 
                     \nColnames must be : 'ensembl.transcriptID','Is.Oncogene','Is.Tumor.Suppressor.Gene'")
      
    }
    } else {
    warning("Oncogene list was empty. Annotating with NAs")
    results_granges_df <- results_granges_df %>% mutate(Is.Oncogene = NA,Is.Tumor.Suppressor.Gene=NA )
  }
  

  
  # collapse annotation per cluster -----------------------------------------------------
  results_granges_df_annot <- results_granges_df %>%
    distinct(clusterID,annot.type,annot.gene_id,annot.gene_name,annot.gene_biotype,Is.Oncogene,Is.Tumor.Suppressor.Gene)%>%
    group_by(clusterID,annot.gene_id,annot.gene_name,annot.gene_biotype,,Is.Oncogene,Is.Tumor.Suppressor.Gene) %>%
    summarise(pos =toString(annot.type)) %>% 
    mutate(position=case_when(str_detect(pos,"exon") ~ "exon", TRUE ~ "intron")) %>%
    dplyr::select(-pos) %>% 
    group_by(clusterID) %>% 
    summarise(gene_ensemblID = toString(annot.gene_id),
              Symbol = toString(annot.gene_name),
              gene_type = toString(annot.gene_biotype),
              position = toString(position),
              Is.Oncogene = toString(Is.Oncogene),
              Is.Tumor.Suppressor.Gene = toString(Is.Tumor.Suppressor.Gene))
  
  

  # Annotate ------------------------------------------------------------------
  
  results_granges_df_annot <- results_df %>%  left_join(results_granges_df_annot, by = "clusterID")
  
  # Save to csv file  -----------------------------------------------------------------

  
  write.table(results_granges_df_annot, str_replace(output,pattern = "xlsx$","csv"),sep=",",row.names = F,quote=F)
  
  
  # save to excel file -----------------------------------------------------
  
  sp <- split(results_granges_df_annot,results_granges_df_annot$library)
  writexl::write_xlsx(sp,path = output,col_names = T,format_headers = T)
  
}
