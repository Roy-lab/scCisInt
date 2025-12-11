#Set Defaults ---- 
options(scipen=999)

source('scCisInt_aux_functions.R')

#Load Libraries ----  
quiet_library("data.table")
quiet_library("dplyr")
quiet_library("tidyr")
quiet_library("caret")
quiet_library("cluster")    # clustering algorithms
#library(factoextra) # clustering algorithms & visualization


args=commandArgs(trailingOnly = T);
if(length(args)==0)
{
  stop("Usage: Rscript scCisInt_postprocessing_kmeans.R [infile (str)] [outfile (str)]",call. = FALSE)
}else if(length(args)>0)
{
  print("Apply Kmeans to scCisInt results and add ranking and distance metrics.")
}



## Input Arguments -----  
infile = args[1]
outfile = args[2]

pred_in=read_head_detect(infile) ## Read all prediction for TSS bins
pred <- pred_in[[1]]

## Set reasonable defaults if name do not already exist
if(!pred_in[[2]])
{
  names(pred)[1:3]=c("interact_region","promoter_region","pred")
}



names(pred)[1:3]=c("interact_region","promoter_region","pred") ## rename values of predicted interactions
pred$RF_rank=rank(desc(pred$pred)) ## Add a ranking 
bin1=do.call(rbind,strsplit(as.character(pred$interact_region),"_")) ## get enhancer starts 
mid1kb=(as.numeric(bin1[,2])+as.numeric(bin1[,3]))/2 ## Midpoint of enhancer start
bin2=do.call(rbind,strsplit(as.character(pred$promoter_region), "_"))
geneloc=(as.numeric(bin2[,2])+as.numeric(bin2[,3]))/2 ## Midpoint of enhancer start 
pred$distance=mid1kb-geneloc ## distance to gene tss 
id0=which(pred$distance>0) ## Find right directional

## All of the above code was for prepping for HiC. We no longer need this. SO we skip to 
k3 <- kmeans(pred$pred, centers = 3, nstart = 25)
new_order <- order(k3$centers[,1], decreasing = TRUE)
map <- setNames(1:3, new_order)
pred$kmeans=map[as.character(k3$cluster)]

pred$clustmean=0.0
for(cl in 1:3){
  pred$clustmean[which(pred$kmeans==cl)]=mean(pred$pred[pred$kmeans == cl])
}

fwrite(pred, file=outfile, quote=F, sep="\t", row.names = F, col.names = T)



