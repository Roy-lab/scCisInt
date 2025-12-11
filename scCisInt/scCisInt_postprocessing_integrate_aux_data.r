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
  stop("At least one argument must be supplied(inputfile).n", call. = FALSE)
} else if (length(args) > 0)
{
  print("Apply Kmeans to scCisInt results and add ranking and distance metrics.")
}

## Input Arguments -----
infile = args[1]
auxfile = args[2]
auxfile_type = args[3]
outfile = args[4]

in_list <- read_head_detect(infile)
pred <- in_list[[1]]
is_header <- in_list[[2]]

## Set reasonable defaults for scCisInt Inputs.
if (!is_header)
{
  names(pred)[1:3] = c("interact_region", "promoter_region", "pred")
}



if (auxfile_type == "enhancer") {
  pred_enhancers <- pred$interact_region
  pred_enhancers <- strsplit(pred_enhancers, split = '_')
  pred_enhancers_gr <- GRanges(
    seqnames = sapply(pred_enhancers, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(pred_enhancers, function(x) {as.numeric(x[[2]])}),
      end = sapply(pred_enhancers, function(x) {as.numeric(x[[3]])}-1)
    )
  )
  
  ## load aux data
  aux_data_in <- read_head_detect(auxfile)
  aux_data <- aux_data_in[[1]]
  is_header <- aux_data_in[[2]]
  
  if (!is_header ||
      !all(names(aux_data)[1:3] ==  c("chr", "start", "end"))) {
    names(aux_data)[1:3] = c("chr", "start", "end")
  }
  
  aux_data_gr <- makeGRangesFromDataFrame(
    aux_data,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  
  ## Identify all pairwise intersects
  hits <- findOverlaps(pred_enhancers_gr, aux_data_gr)
  
  ## Create annotation object from intersects
  gr_res <- pred_enhancers_gr[queryHits(hits)]
  mcols(gr_res) <- cbind(mcols(gr_res), mcols(aux_data_gr)[subjectHits(hits), , drop = FALSE])
  
  ## Add to predictions cleanly
  gr_res_df <- select(as.data.frame(gr_res), !any_of(c("width", "strand")))
  gr_res_df <- mutate(gr_res_df, interact_region = paste(seqnames, start, end, sep = "_"))
  gr_res_df <- select(gr_res_df, !any_of(c("seqnames", "start", "end")))
  pred <- left_join(pred, gr_res_df)
  
} else if (auxfile_type == "promoter") {
  pred_promoters <- pred$promoter_region
  pred_promoters <- strsplit(pred_promoters, split = '_')
  pred_promoters_gr <- GRanges(seqnames = sapply(pred_promoters, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(pred_promoters, function(x) {as.numeric(x[[2]])}),
      end = sapply(pred_promoters, function(x) {as.numeric(x[[3]])}-1)
    )
  )
  
  ## load aux data
  aux_data_in <- read_head_detect(auxfile)
  aux_data <- aux_data_in[[1]]
  is_header <- aux_data_in[[2]]
  
  if (!is_header ||
      !all(names(aux_data)[1:3] ==  c("chr", "start", "end"))) {
    names(aux_data)[1:3] = c("chr", "start", "end")
  }
  
  aux_data_gr <- makeGRangesFromDataFrame(
    aux_data,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  
  ## Identify all pairwise intersects
  hits <- findOverlaps(pred_promoters_gr, aux_data_gr)
  
  ## Create annotation object from intersects
  gr_res <- pred_promoters_gr[queryHits(hits)]
  mcols(gr_res) <- cbind(mcols(gr_res), mcols(aux_data_gr)[subjectHits(hits), , drop = FALSE])
  
  ## Add to predictions cleanly
  gr_res_df <- select(as.data.frame(gr_res), !any_of(c("width", "strand")))
  gr_res_df <- mutate(gr_res_df, interact_region = paste(seqnames, start, end, sep = "_"))
  gr_res_df <- select(gr_res_df, !any_of(c("seqnames", "start", "end")))
  pred <- left_join(pred, gr_res_df)
  
  
} else if (auxfile_type == "edge_undirected") {
  ## We assume HiC-interaction calling format e.g. Bin1 \t Bin2 \t extra columns of any length
  ## Prep enhancer
  pred_enhancers <- pred$interact_region
  pred_enhancers <- strsplit(pred_enhancers, split = '_')
  pred_enhancers_gr <- GRanges(seqnames = sapply(pred_enhancers, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(pred_enhancers, function(x) {as.numeric(x[[2]])}),
                     end = sapply(pred_enhancers, function(x) {as.numeric(x[[3]])} - 1)  ### This is stupid but needs to be here for open ending on scCIsInt bins
    )
  )
  
  ## Prep promoter
  pred_promoters <- pred$promoter_region
  pred_promoters <- strsplit(pred_promoters, split = '_')
  pred_promoters_gr <- GRanges(
    seqnames = sapply(pred_promoters, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(pred_promoters, function(x) {as.numeric(x[[2]])}),
                     end = sapply(pred_promoters, function(x) {as.numeric(x[[3]])} - 1)
    )
  )
  
  ## Aux data load
  aux_data_in <- read_head_detect(auxfile)
  aux_data <- aux_data_in[[1]]
  is_header <- aux_data_in[[2]]
  
  if (!is_header) {
    names(aux_data)[1:2] = c("bin1", "bin2")
  }
  
  bin1 <- strsplit(aux_data$bin1, split = '_')
  bin1_gr <- GRanges(
    seqnames = sapply(bin1, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(bin1, function(x) {as.numeric(x[[2]])}),
                     end = sapply(bin1, function(x) {as.numeric(x[[3]])}) - 1)
  )
  
  bin2 <- strsplit(aux_data$bin2, split = '_')
  bin2_gr <- GRanges(seqnames = sapply(bin2, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(bin2, function(x) {as.numeric(x[[2]])}),
                     end = sapply(bin2, function(x) {as.numeric(x[[3]])}) - 1
    )
  )
  
  ## Find possible intersects (e in bin1, p in bin2)
  
  pred_enhancers_to_bin1  <- findOverlaps(pred_enhancers_gr, bin1_gr)
  pred_promoters_to_bin2  <- findOverlaps(pred_promoters_gr, bin2_gr)
  
  edges_keep <- intersect(pred_enhancers_to_bin1, pred_promoters_to_bin2)
  aux_data_filtered <- aux_data[subjectHits(edges_keep), ]
  pred_match <- pred[queryHits(edges_keep), ]
  
  combine_1 <- cbind(pred_match, aux_data_filtered)
  
  
  ## Find possible intersects (e in bin2, p in bin1)
  
  pred_enhancers_to_bin2  <- findOverlaps(pred_enhancers_gr, bin2_gr)
  pred_promoters_to_bin1  <- findOverlaps(pred_promoters_gr, bin1_gr)
  
  edges_keep <- intersect(pred_enhancers_to_bin2, pred_promoters_to_bin1)
  aux_data_filtered <- aux_data[subjectHits(edges_keep), ]
  names(aux_data_filtered) <- c("bin2", "bin1") 
  pred_match <- pred[queryHits(edges_keep), ]
  
  combine_2 <- cbind(pred_match, aux_data_filtered)
  
  combine <- rbind(combine_1, combine_2)
  combine <- distinct(combine)
  
  pred <- left_join(pred, combine)
  
} else if (auxfile_type == "edge_directed") {
  ## Prep enhancer
  pred_enhancers <- pred$interact_region
  pred_enhancers <- strsplit(pred_enhancers, split = '_')
  pred_enhancers_gr <- GRanges(
    seqnames = sapply(pred_enhancers, function(x) {
      x[[1]]
    }),
    ranges = IRanges(
      start = sapply(pred_enhancers, function(x) {
        as.numeric(x[[2]])
      }),
      end = sapply(pred_enhancers, function(x) {
        as.numeric(x[[3]])
      } - 1)  ### This is stupid but needs to be here for open ending on scCIsInt bins
    )
  )
  
  ## Prep promoter
  pred_promoters <- pred$promoter_region
  pred_promoters <- strsplit(pred_promoters, split = '_')
  pred_promoters_gr <- GRanges(seqnames = sapply(pred_promoters, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(pred_promoters, function(x) {as.numeric(x[[2]])}),
                      end = sapply(pred_promoters, function(x) {as.numeric(x[[3]])} - 1)
    )
  )
  
  ## Aux data load
  aux_data_in <- read_head_detect(auxfile)
  aux_data <- aux_data_in[[1]]
  is_header <- aux_data_in[[2]]
  
  if (!is_header) {
    names(aux_data)[1:2] = c("bin1", "bin2")
  }
  
  bin1 <- strsplit(aux_data$bin1, split = '_')
  bin1_gr <- GRanges(
    seqnames = sapply(bin1, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(bin1, function(x) {as.numeric(x[[2]])}), 
                     end = sapply(bin1, function(x) {as.numeric(x[[3]])}) - 1
    )
  )
  
  bin2 <- strsplit(aux_data$bin2, split = '_')
  bin2_gr <- GRanges(
    seqnames = sapply(bin2, function(x) {x[[1]]}),
    ranges = IRanges(start = sapply(bin2, function(x) {as.numeric(x[[2]])}),
                     end = sapply(bin2, function(x) {as.numeric(x[[3]])}) - 1
    )
  )
  
  ## Find possible intersects (e in bin1, p in bin2)
  
  pred_enhancers_to_bin1  <- findOverlaps(pred_enhancers_gr, bin1_gr)
  pred_promoters_to_bin2  <- findOverlaps(pred_promoters_gr, bin2_gr)
  
  edges_keep <- intersect(pred_enhancers_to_bin1, pred_promoters_to_bin2)
  aux_data_filtered <- aux_data[subjectHits(edges_keep), ]
  pred_match <- pred[queryHits(edges_keep), ]
  
  combine_1 <- cbind(pred_match, aux_data_filtered)
  pred <- left_join(pred, combine_1)
  
} else{
  stop(
    "Invalid auxfile_type, supported values are enhancer, promoter, edge_undirected, edge_directed"
  )
}

fwrite(pred, file=outfile, quote=F, sep="\t", row.names = F, col.names = T)

