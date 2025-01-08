options(tidyverse.quiet = TRUE,warn = -1,verbose = F)


library(rmdformats,quietly = T,warn.conflicts = F,verbose = F)


args <- commandArgs(trailingOnly = T)


load(args[1])


rmarkdown::render(input = "../00-pipeline/guideSeq_GNT_report.rmd", 
                  output_format = "readthedown", 
                  output_dir = dirname(args[2]),
                  output_file = basename(args[2]),
                  clean = T)
