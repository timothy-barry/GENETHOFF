options(tidyverse.quiet = TRUE,conflicts.policy = list(warn = FALSE),verbose = F)
options(knitr.kable.NA = '')

library(tidyverse,quietly = T,warn.conflicts = F,verbose = F)
library(readxl,quietly = T,warn.conflicts = F,verbose = F)
library(kableExtra,quietly = T,warn.conflicts = F,verbose = F)
library(GenomicRanges,quietly = T,warn.conflicts = F,verbose = F)
library(yaml,quietly = T,warn.conflicts = F,verbose = F)
library(rmdformats,quietly = T,warn.conflicts = F,verbose = F)
library(ggrepel)
library(data.table,quietly = T,warn.conflicts = F,verbose = F)
library(DT, quietly = T,warn.conflicts = F,verbose = F)


args <- commandArgs(trailingOnly = T)



summary_files = args[1] 
sampleInfo <- read.delim(args[2], sep=";")
config <- read_yaml(args[3])
predicted_files = args[4] 
max_clusters <- as.numeric(args[5])
minUMI_alignments_figure <- as.numeric(args[6])
min_predicted_distance <- as.numeric(args[7])

# debug
# summary_files = "05-Report/g53_a_Cas_NNN_K562_summary.xlsx 05-Report/g20_a_Cas_NNN_K562_summary.xlsx 05-Report/g16_Cas_NNN_K562_summary.xlsx 05-Report/g12_Cas_NNN_K562_summary.xlsx 05-Report/epe18.2_Cas_NNN_K562_summary.xlsx 05-Report/epe10.5_Cas_NGG_K562_summary.xlsx 05-Report/g20_b_Cas_NGG_K562_summary.xlsx 05-Report/g53_b_Cas_NGG_K562_summary.xlsx 05-Report/EBS_Cas_NGG_hDMD_summary.xlsx"
# sampleInfo <- read.delim("sampleInfo.csv",sep=";")
# config=read_yaml("guideSeq_GNT.yml")
# predicted_files = "06-offPredict/human_GCATCATCCTGGTACCAGGANNN.csv 06-offPredict/human_TTTATCACAGGCTCCAGGAANNN.csv 06-offPredict/human_CAGATAACTGGGCCAACCATNNN.csv 06-offPredict/human_TTGTAATCAGCAGTACCATTNNN.csv 06-offPredict/human_GTATCCTCTTGGGGGCCCCTNNN.csv 06-offPredict/human_GTTTGCCTTGTCAAGGCTATNGG.csv 06-offPredict/human_TTTATCACAGGCTCCAGGAANGG.csv 06-offPredict/human_GCATCATCCTGGTACCAGGANGG.csv 06-offPredict/human_TCTTCCGGAACAAAGTTGCTNGG.csv"
# max_clusters = 100
# minUMI_alignments_figure = 1
# min_predicted_distance = 100





# Load demux stat -------------------------------------------------------------------------

files <- list.files("05-Report/", pattern = "stat$", full.names = T)
names(files) <- str_remove(basename(files),"\\.stat")



stats <- lapply(files, read.delim)
stats <- stats %>% 
  bind_rows() %>% 
  distinct() %>% 
  select(file,num_seqs) %>% 
  mutate(library = str_match(file,"/(.+)_R1")[,2]) %>% 
  mutate(step = case_when(str_detect(file, "_R1.UMI.fastq.gz")~"demultiplexed",
                          str_ends(file, "_R1.UMI.ODN.fastq.gz")~"ODN checked",
                          str_ends(file, "R1.UMI.ODN.trimmed.fastq.gz")~NA,
                          str_ends(file, "_R1.UMI.ODN.trimmed.filtered.fastq.gz")~"trimmed-filtered",
         TRUE ~ NA)) %>% 
  filter(!is.na(step)) %>% 
  select(-file) %>% 
  group_by(library) %>% 
  mutate(prop = round(num_seqs/dplyr::first(num_seqs)*100,digits = 2),
         value = paste(format(num_seqs,big.mark=",")," (",prop,"%)",sep="")) %>% 
  pivot_wider(id_cols = "library",names_from = "step", values_from = "value")



## Load summary excel file

 

summary_files <- unlist(str_split(summary_files," "))
summary <- lapply(summary_files,read_excel)

names(summary) <- str_remove(basename(summary_files),"\\_summary.xlsx")


libraries_count <- length(summary)


# Get probable ON-targets -----------------------------------------------------------
# based on smallest EDITS in crRNA & PAM

summary <- lapply(summary,function(x){
  
  best <- x %>% 
    slice_min(n = 1,with_ties = T, order_by = N_edits) %>% 
    slice_min(n = 1,with_ties = T, order_by = PAM_indel_count) %>% 
    #filter(N_edits==0) %>% 
    select(library,clusterID) %>% 
    mutate(best=TRUE)

  abundance <- x  %>% 
    filter(!is.na(Alignment)) %>% 
    left_join(best) %>%
    mutate(Relative_abundance = round(N_UMI_cluster / sum(N_UMI_cluster) *100,digits = 2)) %>% 
    select(library,clusterID, best,Relative_abundance )
  
  x <-x %>% left_join(abundance)
  return(x)
  
  
})


best_aligns <- summary %>%  
  bind_rows() %>% 
  filter(best==T) %>%
  arrange(library) %>% 
  select(library, clusterID,chromosome,cut_modal_position,cut_gRNa_alignment, Edits_gRNA = N_edits, Edits_PAM=PAM_indel_count,N_IS_cluster, N_UMI_cluster,Relative_abundance)





# figure 7 ------------------------------

files <- list.files("04-IScalling/", pattern="UMIs_per_IS_in_Cluster.bed", full.names = T)
names(files) <- str_remove(basename(files),".UMIs_per_IS_in_Cluster.bed")

IS <- lapply(files, read.delim,header=F)

IS <- IS %>% bind_rows(.id="library")

IS <- IS %>% inner_join(best_aligns, by = c("library","V11"="clusterID"))

IS <- IS %>% mutate(rel_dist_gRNA = V2 - cut_gRNa_alignment,
                    rel_dist_mod = V2 - as.numeric(cut_modal_position))

IS_stat <- IS %>%
  group_by(library,V6,chromosome,cut_gRNa_alignment ,rel_dist_gRNA,clusterID=V11 ) %>%
  summarise(count=sum(V7),
            UMI = log10(sum(V7)+1)) %>%
  mutate(UMI = case_when(V6 == "+"~UMI,
                         TRUE ~ -UMI))

fig7 <- ggplot(IS_stat, aes(rel_dist_gRNA,UMI,fill = factor(sign(UMI)))) +
  geom_vline(xintercept = 0, lty = 2,col = "black") +
  geom_col(show.legend = F) + 
  facet_wrap(~library+paste(chromosome,cut_gRNa_alignment ),scales= "free_y",ncol = 3) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  labs(x = "Distance to cut site (bp)", y = "log10(UMI)")+
  scale_fill_discrete()


# Get stats -------------------------------------------------------------------------

stats_summary <- lapply(seq_along(summary), function(x){
  
  df <- summary[[x]]
  
  df %>% summarise(
                   Reads_aligned = sum(N_reads_cluster),
                   UMIs = sum(N_UMI_cluster),
                   Insertions = sum(N_IS_cluster),
                   total = n(),
                   "dual Orientation" = length(which(N_orientations_cluster==2)),
                   "crRNA matched"= length(which(!is.na(Alignment))),
                   "multiple Cuts" = length(which(N_IS_cluster>3)))
  }
  )

names(stats_summary) <- names(summary)

stats_summary <- stats_summary %>% bind_rows(.id="library")


# Load predictions ------------------------------------------------------------------
cat("Reading predicted cutting sites\n")
predicted_files <- unlist(str_split(predicted_files," "))
predict_gRNA <- lapply(predicted_files,data.table::fread  )
names(predict_gRNA) <- str_remove(basename(predicted_files),"\\.csv")



# Perform analysis for each library ------------------------------------------------------------------

tables_html <- lapply(seq_along(summary), function(x){
  lib_name <- names(summary)[x]
  cat("processing",lib_name,"\n")
  
  
  
  df <- summary[[x]] 
  
  df_grna <- df %>% 
    dplyr::filter(!is.na(Alignment))
  
  if(nrow(df_grna)>0){

    ### load prediction
    if(nrow(df_grna)>0){
      grna_name <- sampleInfo %>% 
        filter(sampleName == names(summary)[x]) %>% 
        mutate("grna_name" = paste(Genome,"_",gRNA_sequence,PAM_sequence,sep="")) %>% pull(grna_name)
      
      if(length(grna_name)==1){
        cat("Using .....", grna_name,"\n")
        
        if(nrow(predict_gRNA[[grna_name]])>0){

          predict <- predict_gRNA[[grna_name]] %>% 
            filter(str_starts(Chromosome,"chr")) %>% 
            mutate(Chromosome = str_remove_all(Chromosome,"_.+$"))


          gr_predict <- makeGRangesFromDataFrame(predict,start.field = "Position", end.field = "Position",seqnames.field = "chromosome",keep.extra.columns = F)
          gr_data <- makeGRangesFromDataFrame(df_grna, start.field = "cut_gRNa_alignment", end.field = "cut_gRNa_alignment",seqnames.field = "chromosome",keep.extra.columns = F)
          near_df <- distanceToNearest(gr_data,gr_predict)
          near_df <- near_df[near_df@elementMetadata$distance <= min_predicted_distance,]
          
          predict <- predict %>%
            distinct(chromosome=Chromosome,DNA=toupper(AlignedText),crRNA=toupper(AlignedTarget))
          
          x3 <- df_grna %>% 
            unite("DNA", seq_gDNA,pam_gDNA,sep="",remove = F ) %>% 
            unite("crRNA",seq_gRNA,pam_gRNA,sep="",remove=F) %>% 
            select(chromosome, clusterID,DNA,crRNA)
          
          x4 <- x3 %>% semi_join(predict) %>% 
            select(chromosome,clusterID) %>% 
            mutate(predicted_alignment="yes")

          x3 <- df_grna %>% 
            left_join(x4) %>% 
            replace_na(list(predicted_alignment="no")) #%>% 
            
          x3$predicted_position <- "no"
          x3$predicted_position[near_df@from] <- "yes"
          
          x3 <- x3 %>% mutate(predicted = case_when(predicted_position=="yes" | predicted_alignment =="yes" ~ 'yes',
                                                    TRUE ~ "no"),
                              bulge = case_when(str_detect(seq_gDNA,"-") & str_detect(seq_gRNA,"-") ~ "both",
                                                str_detect(seq_gRNA,"-") ~ "gDNA",
                                                str_detect(seq_gDNA,"-") ~ "gRNA",
                                                TRUE ~ "none")) %>% 
            select(clusterID,predicted,bulge)
          
          
          df <- df %>% left_join(x3)
          return(df)
          
          
        }
      }
    }
  }
}
)

names(tables_html) <- names(summary)


# figure 3 --------------------------------------------------------------------------

data_fig3 <- lapply(tables_html, function(x){
  
  x %>% group_by(chromosome,"has_gRNA?"=!is.na(Alignment)) %>% 
    summarise(clusters = n(), reads = sum(N_reads_cluster), UMI = sum(N_UMI_cluster))
}
)

data_fig3 <- data_fig3 %>% bind_rows(.id="library") %>% mutate(chromosome = factor(chromosome, levels = paste("chr",c(1:22,"X","Y","M"),sep=""),),
                                                                       "has_gRNA?"= factor(`has_gRNA?`,levels = c("TRUE","FALSE")))

fig3 <- ggplot(data_fig3, aes( chromosome, clusters,fill = `has_gRNA?`)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))





# figure 3b --------------------------------------------------------------------------
fig3b <- ggplot(data_fig3, aes( chromosome, UMI,fill = `has_gRNA?`)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))






# figure 4 --------------------------------------------------------------------------

data_fig4 <- lapply(tables_html, function(x){
  
  x %>% filter(!is.na(Alignment)) %>% group_by(chromosome,predicted) %>%  summarise(clusters = n(), reads = sum(N_reads_cluster), UMI = sum(N_UMI_cluster))
}
)

data_fig4 <- data_fig4 %>% bind_rows(.id="library") %>% mutate(chromosome = factor(chromosome, levels = paste("chr",c(1:22,"X","Y","M"),sep="")),
                                                               predicted= factor(predicted,levels = c("yes","no")))

fig4 <- ggplot(data_fig4, aes( chromosome, clusters,fill=predicted)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))

# figure 4b --------------------------------------------------------------------------


fig4b <- ggplot(data_fig4, aes( chromosome, UMI,fill=predicted)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))



# figure 5 --------------------------------------------------------------------------

fig5 <- lapply(tables_html, function(x){
  
  x %>% filter(N_UMI_cluster>minUMI_alignments_figure, !is.na(Alignment))
}
)

fig5_data <- fig5 %>% bind_rows(.id="library") %>% group_by(library) %>%  mutate(rank_desc = row_number(dplyr::desc(N_UMI_cluster)))



fig5 <- ggplot(fig5_data, aes(rank_desc,N_UMI_cluster)) +
  geom_step(col = "black", direction = "hv") +
  geom_point(pch=19)+
  facet_wrap(~library, ncol = 3, scales = "free") + 
  theme_bw(base_size = 12)+
  scale_y_log10() + 
  scale_x_log10() + 
  ggrepel::geom_text_repel(data = . %>% group_by(library) %>% slice_head(n=3), 
                           aes(rank_desc,N_UMI_cluster,label=paste(chromosome,cut_gRNa_alignment)),inherit.aes = F,col = "purple",cex=3,nudge_x = 1,nudge_y = 0.5,force = 5,direction= "x")+
  geom_point(data = . %>% filter(N_edits==0,PAM_indel_count==0), aes(rank_desc,N_UMI_cluster),col = "red",cex= 2,inherit.aes = F)+
  labs(x = "Ranked clusters by decreasing abundance", y = "UMI count per cluster",caption = paste("Only clusters with >",minUMI_alignments_figure, "UMIs",sep=" "))+
  scale_color_manual(values = c("black","green3"))





# figure 6 --------------------------------------------------------------------------


data_fig6 <- lapply(tables_html, function(x){
  
  x %>% filter(!is.na(Alignment))
}
)

data_fig6 <- data_fig6 %>% bind_rows(.id="library") 

fig6 <- ggplot(data_fig6, aes(N_edits,fill = predicted)) +
  geom_bar(col = "black", stat = ) +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  theme_bw(base_size = 12)+
  scale_x_continuous(breaks = c(0:10))+
  scale_fill_manual(values = c("black","green3"))



# figure 6b --------------------------------------------------------------------------


fig6b <- ggplot(data_fig6, aes(PAM_indel_count,fill = predicted)) +
  geom_bar(col = "black", stat = ) +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  theme_bw(base_size = 12)+
  scale_x_continuous(breaks = c(0:10))+
  scale_fill_manual(values = c("black","green3"))







# tables off targets ----------------------------------------------------------------


tables_off <- lapply(seq_along(tables_html),function(x){
  
  cat("Generating html table ",names(tables_html)[x],"\n")
  
  df <- tables_html[[x]] %>% filter(!is.na(Alignment))  ## keep only clusters with gRNA match
  
  x3 <-  df  %>% ungroup %>% 
    slice_max(n = max_clusters,order_by = N_UMI_cluster,with_ties = F) %>% 
    filter(N_UMI_cluster > minUMI_alignments_figure) %>% 
    rowwise() %>%   # this step is very slow, try to find an alternative ASAP (maybe use data.table package ?)
    mutate(seq_gRNA_html = 
             paste(text_spec(background_as_tile=F,monospace = T,
                             strsplit(
                               seq_gRNA,split="")[[1]],
                             background = recode(strsplit(
                               seq_gRNA,split="")[[1]], A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = ""),
           seq_gDNA_html = 
             paste(text_spec(background_as_tile=F,monospace = T,
                             strsplit(
                               seq_gDNA,split="")[[1]],
                             background = recode(strsplit(
                               seq_gDNA,split="")[[1]], A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = ""),
           pam_gRNA_html = 
             paste(text_spec(background_as_tile=F,monospace = T,
                             strsplit(
                               pam_gRNA,split="")[[1]],
                             background = recode(strsplit(
                               pam_gRNA,split="")[[1]],N="grey", A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = ""),
           pam_gDNA_html = 
             paste(text_spec(background_as_tile=F,monospace = T,
                             strsplit(
                               pam_gDNA,split="")[[1]],
                             background = recode(strsplit(
                               pam_gDNA,split="")[[1]], A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = "")) %>% 
    
    mutate(alignment_html = paste("gRNA: ",seq_gRNA_html," ",pam_gRNA_html," <br>gDNA: ",seq_gDNA_html," ",pam_gDNA_html,sep=""),
           Symbol_html = str_replace_all(Symbol,", "," <br>"),
           position_html = str_replace_all(position,", "," <br>"))
  
  x3$predicted_alignment_html <- cell_spec(x3$predicted, background = ifelse(x3$predicted == "yes", "#129749", "white"))

  
  x3 <- x3 %>% 
    mutate("cut offset" = as.numeric(cut_modal_position) - as.numeric(cut_gRNa_alignment)) %>% 
    select(library,chromosome,cut_gRNa_alignment, `cut offset`,
                      alignment=alignment_html,
                      UMIs=N_UMI_cluster,
                      "Edits crRNA"=N_edits, N_mismatches, n_indels, soft_trim,
                      "Edits pam"= PAM_indel_count,
                      Symbol=Symbol_html,
                      Position=position_html,
                      predicted= predicted_alignment_html,
                      bulge) %>% 
    unite(col = "gRNA cut<br>position",chromosome,cut_gRNa_alignment,sep = ":") %>% 
    unite("MM_indels_softClip",N_mismatches, n_indels, soft_trim,sep = "_",remove = T) 
  
  kb <- kbl(x3 %>% select(-library),
            escape = F,
            align=c("l","r","r",rep('c', 9))) %>%
    kable_classic_2(full_width = F,html_font = "helvetica") %>%
    kable_styling(bootstrap_options = c("condensed","hover","stripped"),
                  font_size = 12,
                  fixed_thead = T) %>%

    column_spec(1:11,extra_css = "vertical-align:middle;")

    save_kable(kb,file = paste("05-Report/report-files/",names(summary)[x],"_offtargets.html",sep=""),self_contained=T)


    # dt <- DT::datatable(x3 %>% select(-library),escape = F, extensions = 'Buttons', rownames = FALSE,filter = 'top',
    #                     options = list(
    #                       pageLength = 20,
    #                       autoWidth = F,
    #                       lengthMenu = list(20,50,-1),
    #                       fixedHeader = TRUE,
    #                       dom = 'Blcfrtip',
    #                       buttons = c('copy', 'csv', 'excel'),
    #                       columnDefs = list(list(className = 'dt-center', targets = c(3,5,6,7,8))),
    #                       initComplete = DT::JS(
    #                         "function(settings, json) {",
    #                         "$('body').css({'font-family': 'Calibri', 'font-size': '10px'});",
    #                         "}"
    #                       )
    #                     )
    # )
    # 
    # return(paste("05-Report/report-files/",names(summary)[x],"_offtargets.html",sep=""))
    # 
     #htmlwidgets::saveWidget(widget = dt, file = paste(names(summary)[x],"_offtargets.html",sep=""),selfcontained = T)

})

names(tables_off) <- names(tables_html)



# Save report data ------------------------------------------------------------------

save(list = c("tables_off", "tables_html","sampleInfo", "stats", "stats_summary","best_aligns", grep(x=ls(),"fig",value = T)), file = "05-Report/report.rdata")


# Make report -----------------------------------------------------------------------






