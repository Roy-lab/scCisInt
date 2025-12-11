#scCisInt_postprocessing_combine_all_predictions
#Set Defaults ----
options(scipen = 999)

source('scCisInt_aux_functions.R')

#Load Libraries ----
quiet_library("data.table")
quiet_library("dplyr")
quiet_library("tidyr")
quiet_library("caret")  # data manipulatio
quiet_library("cluster")    # clustering algorithms
quiet_library("GenomicRanges")


args = commandArgs(trailingOnly = T)

if (length(args) == 0)
{
  stop("Usage: Rscript scCisInt_postprocessing_combine_all_predictions.R [indir (string)] [outfile (string)] [file suffix, optional (default = consensus_EP_predictions.txt)] [filter kmmeans cluster, optional (default = 1)]", call. = FALSE)
} else if (length(args) > 0)
{
  print("Combining scCisInt predictions.")
}

## Input Arguments -----
indir = args[1]
outfile = args[2]


if(length(args) == 3){
  file_suffix <- args[3]
}else{
  file_suffix <- 'consensus_EP_predictions.txt'
}

if(length(args) == 4)
{
  filter_kmeans_cluster = as.numeric(args[4])
}else{
  filter_kmeans_cluster = 1
}

files <- list.files(path = indir, pattern = file_suffix)

full_dat <- data.table()
for (f in files)
{
  dat_in <- read_head_detect(paste(indir, f, sep ="/"))
  dat <- dat_in[[1]]
  is_header <- dat_in[[2]]
  
  if(!is_header)
  {
    names(dat) <- c("interact_region", "promoter_region", "pred")
  }
  
  dat_keep <- filter(dat, kmeans <= filter_kmeans_cluster)
  full_dat <- bind_rows(full_dat, dat_keep)
}
  
fwrite(full_dat, file=outfile, quote=F, sep="\t", row.names = F, col.names = T)
  
