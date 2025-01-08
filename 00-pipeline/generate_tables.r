options(tidyverse.quiet = TRUE,warn = -1,verbose = F)


library(tidyverse,quietly = T,warn.conflicts = F,verbose = F)
library(readxl,quietly = T,warn.conflicts = F,verbose = F)
library(kableExtra,quietly = T,warn.conflicts = F,verbose = F)
library(GenomicRanges,quietly = T,warn.conflicts = F,verbose = F)
library(yaml,quietly = T,warn.conflicts = F,verbose = F)
library(rmdformats,quietly = T,warn.conflicts = F,verbose = F)
library(ggrepel)
library(data.table)


args <- commandArgs(trailingOnly = T)





config <- read_yaml(args[3])
sampleInfo <- read.delim(args[2], sep=";")


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
                          str_ends(file, "R1.UMI.ODN.trimmed.fastq.gz")~"trimmed",
                          str_ends(file, "_R1.UMI.ODN.trimmed.filtered.fastq.gz")~"filtered",
         TRUE ~ NA)) %>% 
  filter(!is.na(step)) %>% 
  select(-file) %>% 
  group_by(library) %>% 
  mutate(prop = round(num_seqs/dplyr::first(num_seqs)*100,digits = 2),
         value = paste(num_seqs," (",prop,"%)",sep="")) %>% 
  pivot_wider(id_cols = "library",names_from = "step", values_from = "value")



## Load summary excel file

path = args[1]  # summary file xls
# summary <- path %>% excel_sheets() %>% 
#   set_names() %>% 
#   map(read_excel, path = path)

path=unlist(str_split(path," "))
summary <- lapply(path,read_excel)

names(summary) <- str_remove(basename(path),"\\_summary.xlsx")


libraries_count <- length(summary)


# Get probable ON-targets -----------------------------------------------------------
# based on smallest EDITS in crRNA & PAM

best_aligns <- lapply(summary,function(x){
  
  x %>% 
    slice_min(n = 1,with_ties = T, order_by = N_edits) %>% 
    slice_min(n = 1,with_ties = T, order_by = PAM_indel_count) %>% 
    select(library, clusterID,N_IS_cluster, N_UMI_cluster,N_orientations_cluster,N_edits, PAM_indel_count,cut_modal_position,cut_gRNa_alignment)

})

best_aligns <- best_aligns %>%  bind_rows()

# figure 7 ------------------------------

files <- list.files("04-IScalling/", pattern="ISCluster.bed", full.names = T)
names(files) <- str_remove(basename(files),".collapsefragPerISCluster.bed")

IS <- lapply(files, read.delim,header=F)

IS <- IS %>% bind_rows(.id="library")

IS <- IS %>% inner_join(best_aligns, by = c("library","V11"="clusterID"))

IS <- IS %>% mutate(rel_dist_gRNA = V2 - cut_gRNa_alignment,
                    rel_dist_mod = V2 - as.numeric(cut_modal_position))

IS_stat <- IS %>%
  group_by(library,V6,rel_dist_gRNA,clusterID=V11 ) %>%
  summarise(UMI = log10(sum(V7)+1)) %>%
  mutate(UMI = case_when(V6 == "+"~UMI,
                         TRUE ~ -UMI))

fig7 <- ggplot(IS_stat, aes(rel_dist_gRNA,UMI,fill = factor(sign(UMI)))) +
  geom_vline(xintercept = 0, lty = 2,col = "black") +
  geom_col(show.legend = F) + 
  facet_wrap(~library,scales= "free_y",ncol = 3) +
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
})

names(stats_summary) <- names(summary)

stats_summary <- stats_summary %>% bind_rows(.id="library")


# Load predictions ------------------------------------------------------------------


path = args[4]  # summary file xls
# summary <- path %>% excel_sheets() %>% 
#   set_names() %>% 
#   map(read_excel, path = path)

path=unlist(str_split(path," "))
predict_gRNA <- lapply(path,data.table::fread  )

names(predict_gRNA) <- str_remove(basename(path),"\\.csv")



# Perform analysis for each library ------------------------------------------------------------------

tables_html <- lapply(seq_along(summary), function(x){
  lib_name <- names(summary)[x]
  cat("processing",lib_name,"\n")
  
  df <- summary[[x]]
  if(nrow(df)>0){
    df <- df %>% 
      mutate(score = case_when(N_orientations_cluster==2 ~ 1,
                               TRUE ~ 0))  %>% 
      mutate(score = case_when(N_IS_cluster>1 ~ score +1 ,
                               TRUE ~ score +0)) %>% 
      mutate(score = case_when(N_UMI_cluster>3 ~ score + 1,
                               TRUE ~ score +0)) %>% 
      mutate(score = case_when(N_edits<=6 ~ score + 1,
                               TRUE ~ score +0))%>% 
      mutate(score = case_when(!is.na(grna_orientation)~score + 1,
                               TRUE ~ score +0))
    
    ### filter cluster with high score
    df <- df %>% 
      dplyr::filter(!is.na(Alignment))
    
    ### load prediction
    if(nrow(df)>0){
      grna_name <- sampleInfo %>% filter(sampleName == names(summary)[x]) %>% mutate("grna_name" = paste(Genome,"_",gRNA_sequence,PAM_sequence,sep="")) %>% pull(grna_name)
      
      if(length(grna_name)==1){
        cat("Using .....", grna_name,"\n")
        
        if(nrow(predict_gRNA[[grna_name]])>0){
          predict <- predict_gRNA[[grna_name]] %>% 
            mutate(Chromosome = str_remove_all(Chromosome,"_.+$")) %>% 
            filter(str_starts(Chromosome,"chr"))
          
          gr_predict <- makeGRangesFromDataFrame(predict,start.field = "Position", end.field = "Position",seqnames.field = "chromosome",keep.extra.columns = F)
          gr_data <- makeGRangesFromDataFrame(df, start.field = "cut_gRNa_alignment", end.field = "cut_gRNa_alignment",seqnames.field = "chromosome",keep.extra.columns = F)
          near_df <- distanceToNearest(gr_data,gr_predict)
          near_df <- near_df[near_df@elementMetadata$distance<100,]
          
          
          
          predict <- predict %>%
            distinct(chromosome=Chromosome,DNA=toupper(AlignedText),crRNA=toupper(AlignedTarget))
          
          
          x3 <- df %>% 
            unite("DNA", seq_gDNA,pam_gDNA,sep="",remove = F ) %>% 
            left_join(predict,by=c("chromosome","DNA")) %>% 
            rowwise() %>% 
            mutate(seq_gRNA_html = 
                     paste(text_spec(background_as_tile=F,
                                     strsplit(
                                       seq_gRNA,split="")[[1]],
                                     background = recode(strsplit(
                                       seq_gRNA,split="")[[1]], A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = ""),
                   seq_gDNA_html = 
                     paste(text_spec(background_as_tile=F,
                                     strsplit(
                                       seq_gDNA,split="")[[1]],
                                     background = recode(strsplit(
                                       seq_gDNA,split="")[[1]], A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = ""),
                   pam_gRNA_html = 
                     paste(text_spec(background_as_tile=F,
                                     strsplit(
                                       pam_gRNA,split="")[[1]],
                                     background = recode(strsplit(
                                       pam_gRNA,split="")[[1]],N="grey", A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = ""),
                   pam_gDNA_html = 
                     paste(text_spec(background_as_tile=F,
                                     strsplit(
                                       pam_gDNA,split="")[[1]],
                                     background = recode(strsplit(
                                       pam_gDNA,split="")[[1]], A="#129749",T="#d62839",C="#255c99",G="#f7b32b")),collapse = "")) %>% 
            
            mutate(alignment_html = paste("gRNA: ",seq_gRNA_html," ",pam_gRNA_html," <br>gDNA: ",seq_gDNA_html," ",pam_gDNA_html,sep=""),
                   Symbol_html = str_replace_all(Symbol,", "," <br>"),
                   predicted_alignment = case_when(!is.na(crRNA)~ "yes",
                                                   TRUE ~ "no"))
          x3$predicted_position <- "no"
          x3$predicted_position[near_df@from] <- "yes"
          
          x3$predicted_alignment_html <- cell_spec(x3$predicted_alignment, background = ifelse(x3$predicted_alignment == "yes", "#129749", "white"))
          x3$predicted_position_html <- cell_spec(x3$predicted_position, background = ifelse(x3$predicted_position == "yes", "#129749", "white"))
          
          return(x3)
          
          
        }
      }
    }
  }
}
)

names(tables_html) <- names(summary)


# figure 3 --------------------------------------------------------------------------

data_fig3 <- lapply(summary, function(x){
  
  x %>% group_by(chromosome,"has_crRNA?"=!is.na(Alignment)) %>% summarise(clusters = n(), reads = sum(N_reads_cluster), UMI = sum(N_UMI_cluster))
}
)

data_fig3 <- data_fig3 %>% bind_rows(.id="library") %>% mutate(chromosome = factor(chromosome, levels = paste("chr",c(1:22,"X","Y","M"),sep=""),),
                                                                       "has_crRNA?"= factor(`has_crRNA?`,levels = c("TRUE","FALSE")))

fig3 <- ggplot(data_fig3, aes( chromosome, clusters,fill = `has_crRNA?`)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))





# figure 3b --------------------------------------------------------------------------
fig3b <- ggplot(data_fig3, aes( chromosome, UMI,fill = `has_crRNA?`)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))






# figure 4 --------------------------------------------------------------------------

data_fig4 <- lapply(tables_html, function(x){
  
  x %>% filter(!is.na(Alignment),N_UMI_cluster>3) %>% group_by(chromosome,predicted_position) %>%  summarise(clusters = n(), reads = sum(N_reads_cluster), UMI = sum(N_UMI_cluster))
}
)

data_fig4 <- data_fig4 %>% bind_rows(.id="library") %>% mutate(chromosome = factor(chromosome, levels = paste("chr",c(1:22,"X","Y","M"),sep="")),
                                                     predicted_position= factor(predicted_position,levels = c("yes","no")))

fig4 <- ggplot(data_fig4, aes( chromosome, clusters,fill=predicted_position)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))

# figure 4b --------------------------------------------------------------------------


fig4b <- ggplot(data_fig4, aes( chromosome, UMI,fill=predicted_position)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))



# figure 5 --------------------------------------------------------------------------

fig5 <- lapply(tables_html, function(x){
  
  x %>% filter(!is.na(Alignment),N_UMI_cluster>3)
}
)

fig5 <- fig5 %>% bind_rows(.id="library") %>% group_by(library) %>%  mutate(rank_desc = row_number(dplyr::desc(N_UMI_cluster)))



fig5 <- ggplot(fig5, aes(rank_desc,N_UMI_cluster)) +
  geom_step(col = "black", direction = "hv") +
  geom_point()+
  facet_wrap(~library, ncol = 3, scales = "free") + 
  theme_bw(base_size = 12)+
  scale_y_log10() + 
  scale_x_log10() + 
  ggrepel::geom_text_repel(data = . %>% group_by(library) %>% slice_head(n=3), 
                           aes(rank_desc,N_UMI_cluster,label=paste(chromosome,cut_gRNa_alignment)),inherit.aes = F,col = "purple",cex=3,nudge_x = 1,nudge_y = 0.5,force = 5,direction= "x")+
  geom_point(data = . %>% filter(N_edits==0,PAM_indel_count==0), aes(rank_desc,N_UMI_cluster),col = "red",cex= 2,inherit.aes = F)+
  labs(x = "Ranked clusters by decreasing abundance", y = "UMI count per cluster")





# figure 6 --------------------------------------------------------------------------


data_fig6 <- lapply(tables_html, function(x){
  
  x %>% filter(!is.na(Alignment),N_UMI_cluster>3)
}
)

data_fig6 <- data_fig6 %>% bind_rows(.id="library") 

fig6 <- ggplot(data_fig6, aes(N_edits,fill = predicted_position)) +
  geom_bar(col = "black", stat = ) +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  theme_bw(base_size = 12)+
  scale_x_continuous(breaks = c(0:10))



# figure 6b --------------------------------------------------------------------------


fig6b <- ggplot(data_fig6, aes(PAM_indel_count,fill = predicted_position)) +
  geom_bar(col = "black", stat = ) +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  theme_bw(base_size = 12)+
  scale_x_continuous(breaks = c(0:10))







# tables off targets ----------------------------------------------------------------


tables_off <- lapply(seq_along(tables_html),function(x){
  
  cat("processing",names(tables_html)[x],"\n")
  
  df <- tables_html[[x]]
  
  x3 <-  df  %>%
    select(library,chromosome,cut_gRNa_alignment,
           alignment=alignment_html,
           UMI=N_UMI_cluster, "#Edits crRNA"=N_edits,"#Edits pam"= PAM_indel_count,
           Symbol=Symbol_html,
           pred.align= predicted_alignment_html,
           pred.pos = predicted_position_html,
           score) %>% 
    unite(col = "position",chromosome,cut_gRNa_alignment,sep = ":") %>% 
    filter(score>config$minScore) %>% 
    select(-score)
  
  kb <- kbl(x3 %>% select(-library),
            escape = F,
            align=rep('c', 5)) %>%
    kable_classic_2(full_width = F,html_font = "helvetica") %>% 
    kable_styling(bootstrap_options = c("condensed","hover","stripped"),
                  font_size = 12,
                  fixed_thead = T) %>% 
    
    column_spec(1:8,extra_css = "vertical-align:middle;") %>% 
    column_spec(2, monospace = T) 
  
    save_kable(kb,file = paste("05-Report/report-files/",names(summary)[x],"_offtargets.html",sep=""),self_contained=T)
  
})

names(tables_off) <- names(tables_html)



# Save report data ------------------------------------------------------------------

save(list = c("tables_off", "tables_html","sampleInfo", "stats", "stats_summary", grep(x=ls(),"fig",value = T)), file = "05-Report/report.rdata")


# Make report -----------------------------------------------------------------------






