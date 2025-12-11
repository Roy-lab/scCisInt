#!/user/bin/env Rscript
## prepare enhancer and promoter data for each gene:
args=commandArgs(trailingOnly = T);
if(length(args)<5)
{
  stop("Usage: Rscript scATAC_EPpredict_preData_forgenes.R [chromosome (string)] [lower_bound (int)] [upper_bound (int)] [pseudobulk_file (str)] [gene_bin_file (str)] [radius: optional default 1e6 (int)]  [outdir: optional, default 'Results' (str)] [verbose: optional, default 'False' (bool)]" ,call. = FALSE)
}else if(length(args)>0)
{
  print("Start!!")
}
options(scipen=999)

# Args management
chromosome=args[1]
lower_bound=args[2]
upper_bound=args[3]
nmfdata=args[4] 
genebinfile=args[5]

if(length(args) >= 6)
{
  radius=args[6]
}else
{
  radius=1e6
}

if(length(args) >= 7)
{
  outdir=args[7]
}else
{
  outdir="Results"
}

if(length(args) >= 8)
{
  verbose = args[8];
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
g0$gene = strsplit(g0$gene, ";")
g0 <- g0[, .(gene = unlist(gene)), by = setdiff(names(g0), "gene")]
gene_parts=strsplit(g0$gene, '_')
g=g0

colnames(g)[1] <- "bin"
g$gene=sapply(gene_parts, function(x){x[1]})
g$type=sapply(gene_parts, function(x){x[4]})

bin_parts=strsplit(g0$bin1kb, '_')
g$bin_start = sapply(bin_parts, function(x){as.numeric(x[2])})
g$bin_end = sapply(bin_parts, function(x){as.numeric(x[3])})

chr_idx = which(g$chr == chromosome)
bin_start_idx = which(g$bin_start >= lower_bound)
bin_end_idx = which(g$bin_end <= upper_bound)
filtered_regions = Reduce(intersect, list(chr_idx, bin_start_idx, bin_end_idx))

g_filt = g[filtered_regions, ]
genes <- unique(g_filt$gene)



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
