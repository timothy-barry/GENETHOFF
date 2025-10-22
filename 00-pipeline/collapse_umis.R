# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

tab <- data.table::fread(file = input_file) 
tab_summary <- tab |> 
  dplyr::group_by(start, chr, strand) |>
  dplyr::summarize(n_umis = dplyr::n(), read_count = sum(count)) |>
  dplyr::ungroup()

saveRDS(my_results, file = output_file)
