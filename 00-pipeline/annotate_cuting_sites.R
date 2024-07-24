
library(tidyverse,quietly = T, verbose = F)
library(GenomicRanges,quietly = T, verbose = F)
library(annotatr,quietly = T, verbose = F)
library(biomaRt,quietly = T, verbose = F)
library(writexl,quietly = T, verbose = F)


args <- commandArgs(trailingOnly = T)

files <- list.files(path = "05-Report/",pattern = "rdata", recursive = T, full.names = T)

names(files) <- str_remove(basename(files),pattern = ".rdata")

results <- lapply(files, function(x){
  load(x)
  cluster_annnotated %>% arrange(desc(N_UMI_cluster)) %>% 
    mutate(gRNA = as.character(grna)) %>%
    filter(str_starts(chromosome,"chr"))
})


results_df <- bind_rows(results,.id = "library")


# build annotation from bioMart -----------------------------------------------------------------

ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'hsapiens_gene_ensembl',
                      version = 112)


# listEnsemblArchives()
 attributes = listAttributes(ensembl)
# attributes


features = attributes[c(204,205,213,214,215,220,29,62,78),1]

mart <- getBM(attributes = features,
              mart = ensembl)



mart = mart %>% 
  mutate(strand = case_when(strand == "1" ~ "+",
                            TRUE ~ "-"),
         chromosome_name = case_when(chromosome_name %in% c(1:100,"X","Y")~ paste("chr",chromosome_name,sep=""),
                                     chromosome_name == "MT" ~ "chrM",
                                     TRUE ~ chromosome_name)) %>% 
  filter(str_starts(chromosome_name, "chr"))



mart_gr <- makeGRangesFromDataFrame(df = mart , seqnames.field = "chromosome_name",start.field = "start_position",end.field = "end_position",keep.extra.columns = T,ignore.strand = F,strand.field = "strand")


# convert insertions to gRanges -----------------------------------------------------

results_granges <-   makeGRangesFromDataFrame(results_df, ignore.strand = T,keep.extra.columns = T, strand.field = "strand_guide",start.field = "start_cluster",end.field = "end_cluster",seqnames.field = "chromosome", na.rm = T)



# results_granges <- lapply(results, function(x) {
#   makeGRangesFromDataFrame(x, ignore.strand = T,keep.extra.columns = T, strand.field = "strand_guide",start.field = "start_cluster",end.field = "end_cluster",seqnames.field = "chromosome", na.rm = T)
# })


# annotate gRanges -------------------------------------------------------------------

results_granges_annot <- annotate_regions(regions = results_granges,annotations = mart_gr,ignore.strand = T,minoverlap = 1)
results_granges_df <- data.frame(results_granges_annot)

# 
# results_granges_annot <- lapply(results_granges, function(x) {
#   annotate_regions(regions = x,annotations = mart_gr,ignore.strand = T,minoverlap = 1)}
#   )
# 
# 
# results_granges_df <- lapply(results_granges_annot, function(x){
#   x %>%  data.frame
# })


# collapse annotation per cluster -----------------------------------------------------
results_granges_df_annot <- results_granges_df %>% 
  distinct() %>% 
  group_by(across(c(library,seqnames:gRNA ))) %>%
  summarise(geneID=toString(annot.entrezgene_id),
            gene_ensemblID = toString(annot.ensembl_gene_id),
            Symbol = toString(annot.hgnc_symbol),
            gene_type = toString(annot.gene_biotype)) %>%
  arrange(desc(N_UMI_cluster))


# results_granges_df_annot <- lapply(results_granges_df, function(x) {
#   x %>% 
#     
#     distinct() %>% 
#     group_by(across(seqnames:gRNA )) %>%
#     summarise(geneID=toString(annot.entrezgene_id),
#               gene_ensemblID = toString(annot.ensembl_gene_id),
#               Symbol = toString(annot.hgnc_symbol),
#               gene_type = toString(annot.gene_biotype)) %>%
#     arrange(desc(N_UMI_cluster))}
# )
results_granges_df_annot <- results_df %>%  left_join(results_granges_df_annot)



write.table(results_granges_df_annot, "05-Report/summary.tsv",sep="\t",row.names = F,quote=F)

# save to excel file -----------------------------------------------------

sp <- split(results_granges_df_annot,results_granges_df_annot$library)
writexl::write_xlsx(sp,path = "05-Report/summary.xlsx",col_names = T,format_headers = T)



#make graph
# 
# lapply(results_granges_df_annot, function(x){
#   
#   y = x %>% ungroup %>%  filter(!is.na(sequence)) %>% 
#    slice_max(order_by = N_UMI_cluster,n=10)
#   
#   z <- DNAStringSet(x = c(unique(y$gRNA),y$sequence))
#   names(z) <- c("gRNA",y$clusterID)
#   
#   ggmsa::ggmsa(z,color = "Shapely_NT",
#                char_width = 0.5,
#                disagreement = F,
#                use_dot = T,
#                consensus_views = T, 
#                font = "DroidSansMono",
#                seq_name = T) + 
#     coord_fixed(ratio = 1) +
#     labs(title = "test")+
#     scale_y_discrete(limits = factor(rev(names(z))))
# 
# })

