#aux_functions.R

quiet_library <- function(pkg) {
  suppressPackageStartupMessages(
    suppressWarnings(
      library(pkg, character.only = TRUE)
    )
  )
}


#Function for save loading with variable header:
read_head_detect <- function(infile, sep = "\t"){
  first_line <- readLines(infile, n = 1) ## Read first line. 
  fields <- unlist(strsplit(first_line, split = sep))
  is_header <- all(is.na(suppressWarnings(as.numeric(fields)))) ## Check if any are numeric
  dat <- read.table(infile, header = is_header)
  return(list(dat, is_header))
}
