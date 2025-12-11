#scCisInt_postprocessing_enhancer_scores

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
quiet_library("purrr")


args = commandArgs(trailingOnly = T)

if (length(args) != 2 && length(args) != 3)
{
  stop("Usage: Rscript scCisInt_postprocessing_combine_all_predictions.R [indir (string)] [outfile (string)] [config_file, optional (string)]", call. = FALSE)
} else if (length(args) > 0)
{
  print("summarizing enhancer scores, scCisInt predictions.")
}

## Input Arguments -----
infile = args[1]
outfile = args[2]

if(length(args) == 3){
  rules_file <- args[3]
}else{
  rules_file <- NULL
}

## Read input file 
dat_in <- read_head_detect(infile)
dat <- dat_in[[1]]
is_header <- dat_in[[2]]

## Rules read ----- 
rules <- NULL
if(!is.null(rules_file)) {
  rules <- read.table(rules_file)
  names(rules) <- c("column", "rules")
  
  if(any(!(rules$column %in% names(dat)))){
    stop("Rules reference missing columns")
  }
}


## Summarize with rules -----
if(!is.null(rules)) {
  
  summarized <- summarize(dat, .by = "interact_region", 
                          number_of_promoter_bins = n(), 
                          impscore_avg = mean(pred), 
                          impscore_sum = sum(pred),
                          !!!set_names(
                            map2(rules$column, rules$rule, ~{
                              rlang::expr(summary_fun(!!rlang::sym(.x), !!.y))
                            }),
                            paste(rules$column, rules$rules, sep = "_")
                          )
                        )
  
} else {
  ## Summarize base  -----
  summarized <- summarize(dat, .by = "interact_region", 
                          number_of_promoter_bins = n(), 
                          impscore_avg = mean(pred), 
                          impscore_sum = sum(pred)
  )
}


fwrite(summarized, file=outfile, quote=F, sep="\t", row.names = F, col.names = T)
