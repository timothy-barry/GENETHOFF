options(tidyverse.quiet = TRUE,warn = -1,verbose = F)





library(tidyverse,quietly = T, verbose = F)
library(GenomicRanges,quietly = T, verbose = F)
library(annotatr,quietly = T, verbose = F)
library(biomaRt,quietly = T, verbose = F)
library(writexl,quietly = T, verbose = F)
library(rtracklayer,quietly = T, verbose = F)

args <- commandArgs(trailingOnly = T)
gtf <- args[1]
files <- args[2]
output<- args[3]



names(files) <- str_remove(basename(files),pattern = ".rdata")


load(x)

if(exists("cluster_annotated") && nrow(cluster_annotated)>0){
  results_df <- cluster_annotated %>% 
    arrange(desc(N_UMI_cluster)) %>% 
    mutate(gRNA = as.character(grna),
           gRNA_name = grna@metadata$name,
           library = names(files)) %>%
    filter(str_starts(chromosome,"chr"))   # remove scaffolds and unplaced sequences
  rm(cluster_annnotated);rm(y)
}



# build annotation from gtf file

mart_gr <- import(gtf)
mart_gr <- mart_gr[mart_gr$type =="gene"]



# convert insertions to gRanges -----------------------------------------------------

results_granges <-   makeGRangesFromDataFrame(results_df, ignore.strand = T,keep.extra.columns = T, strand.field = "strand_guide",start.field = "start_cluster",end.field = "end_cluster",seqnames.field = "chromosome", na.rm = T)



# annotate gRanges -------------------------------------------------------------------

results_granges_annot <- annotate_regions(regions = results_granges,annotations = mart_gr,ignore.strand = T,minoverlap = 1)
results_granges_df <- data.frame(results_granges_annot)

# collapse annotation per cluster -----------------------------------------------------

results_granges_df_annot <- results_granges_df %>% 
  distinct() %>% 
  group_by(across(c(library,seqnames:gRNA ))) %>%
  summarise(gene_ensemblID = toString(annot.gene_id),
            Symbol = toString(annot.gene_name),
            gene_type = toString(annot.gene_type)) %>%
  arrange(desc(N_UMI_cluster))


results_granges_df_annot <- results_df %>%  left_join(results_granges_df_annot)

write.table(results_granges_df_annot, str_replace(output,pattern = "xlsx$","tsv"),sep="\t",row.names = F,quote=F)


# save to excel file -----------------------------------------------------

sp <- split(results_granges_df_annot,results_granges_df_annot$library)
writexl::write_xlsx(sp,path = output,col_names = T,format_headers = T)


