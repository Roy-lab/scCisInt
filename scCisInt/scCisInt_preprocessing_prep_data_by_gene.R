#!/user/bin/env Rscript
## prepare enhancer and promoter data for each gene:
args=commandArgs(trailingOnly = T);
if(length(args)<3)
{
  stop("Usage: Rscript scATAC_EPpredict_preData_forgenes.R [gene name (str)] [pseudobulk_file (str)] [gene_bin_file (str)] [radius: optional default 1e6 (int)]  [outdir: optional, default 'Results' (str)] [verbose: optional, default 'False' (bool)]" ,call. = FALSE)
}else if(length(args)>0)
{
  print("Start!!")
}
options(scipen=999)

# Args management
genename=args[1]
nmfdata=args[2] 
genebinfile=args[3]

if(length(args) >= 4)
{
  radius=args[4]
}else
{
  radius=1e6
}

if(length(args) >= 5)
{
  outdir=args[5]
}else
{
  outdir="Results"
}

if(length(args) >= 6)
{
  verbose = args[6];
}else
{
  verbose = F
}

# Load accessibility pseudobulk data
library(data.table)
d0=fread(nmfdata)
d0=d0[,-1]
df=t(d0)
d=data.table(Var1=row.names(df),df)
genes=c(genename)


g0 = fread(genebinfile)
gene_parts=strsplit(g0$gene, '_')
g=g0

colnames(g)[1] <- "bin"
g$gene=sapply(gene_parts, function(x){x[1]})
g$type=sapply(gene_parts, function(x){x[4]})



for(i in 1:length(genes)){
  gene=genes[i]
  
  if(verbose){
    print(gene)
  }
  
  idx <- tryCatch(
    {
      which(grepl(gene, g$gene))
    },
    error = function(e) {
      message("Unaable to find gene match in gene id list", e$message)
      return(NULL)   # or integer(0), NA, etc.
    }
  )
  if(length(idx) == 0 || any(is.na(idx))){
    message("Unable to find gene match in gene id list.")
    return(NULL)
  }

  
  for(j in idx)
  {
    gene_region = g$bin[[j]]
    
    if(file.exists(paste0(outdir, "/", gene_region, "_enhancer_data.txt")))
    {
      if(verbose){
        print(paste0("Input data files already prepared. Skipping ", gene_region, '.'))
      }
      next
    }
    
    gene_region_split= unlist(strsplit(gene_region, '_'))
    chr = gene_region_split[1]
    geneloc = (as.numeric(gene_region_split[2]) + as.numeric(gene_region_split[3]))/2 
    
    ## Identify region/bin set for each gene 
    d1 = d[grep(paste0(chr, "_"), d$Var1), ]
    d1$Var1=as.character(d1$Var1)
    region=do.call(rbind,strsplit(d1$Var1,"_"))
    startp=as.numeric(region[,2])
    endp=as.numeric(region[,3])
    midp=(startp+endp)/2
    distance=abs(midp-geneloc)
    ide = which(distance <= radius)
    d2=d1[ide,]
    
    ## Prepare enhancer and promoter matrix files
    ### Check psuedobulk accessibility
    csum=colSums(d2[,-1])
    idc=which(csum>0)
    kpcols=c("Var1",names(idc))
    
    
    ### Check promoter presence 
    id=grep(gene_region,d2$Var1)
    if(length(id)==0){
      next 
    }
    
    if(verbose)
    {
      print(gene_region)
    }
    write.table(d2$Var1, file = paste0(outdir, '/', gene_region, "_bins.txt"), quote = F, sep ="\n", row.names = F, col.names = F)
    
    ## atac profiles of candidate enhancers within 1MB:
    enhancer=d2[-id,..kpcols]
    ## atac profiles of promoter:
    promoter=d2[id,..kpcols]
    
    if(verbose)
    {
      print(enhancer[1:3,1:3])
      print(promoter[1,1:3])
    }
    
    enhancer_region=enhancer$Var1
    enhancer_X=data.frame(t(enhancer[,-1]))
    names(enhancer_X)=enhancer_region
    promoter_Y=data.frame(t(promoter[,-1]))
    names(promoter_Y)=gene_region
    
    if(verbose){
      print(paste0("After transformation enhancer cells=",nrow(enhancer_X)," enhancers=",ncol(enhancer_X)))
    }
    
    fwrite(enhancer_X,file=paste0(outdir, "/", gene_region, "_enhancer_data.txt"), quote=F, sep="\t", row.names = T, col.names = T)
    fwrite(promoter_Y,file=paste0(outdir, "/", gene_region, "_promoter_data.txt"), quote=F, sep="\t", row.names = T, col.names = T)
  }
}
  
