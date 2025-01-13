library(rtracklayer,quietly = T, verbose = F,warn.conflicts = F)

## convert gtf file to R object

args <- commandArgs(trailingOnly = T)


# read yaml config file ------------------------------------------------------------------
config <- yaml::read_yaml(args[1])


# Get list of species ---------------------------------------------------------------
species <- names(config$genome)


# build annotation objects ----------------------------------------------------------

# for each species:

for ( my_specie in species){
  cat("Preparing annotation for", my_specie,"\n")
  gtf <- rtracklayer::import(config$genome[[my_specie]]$annotation)
  gtf <- gtf[gtf$type %in% c("gene","exon")]
  saveRDS(object = gtf, file = paste("../02-ressources/",my_specie,".rds",sep = "" ))
  
}







