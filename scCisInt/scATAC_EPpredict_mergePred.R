#!/user/bin/env Rscript
args=commandArgs(trailingOnly = T);
if(length(args)==0)
{ 
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
}else if(length(args)>0)
{ 
  print("Start merging predictions of multiple bins.")
}
options(scipen=999)
genename=args[1]
chr=args[2]
indir="Results/"
outdir="Results/"
library(dplyr)
library(data.table)
#indir="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/multitask_matfact_scatac/data/NMF_sparse_bp_K70/"
#outdir="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/multitask_matfact_scatac/Results/NMF_sparse_bp_K70/subsample/analysis/"
#genename="Nanog"
#genename="Pou5f1"  # all 1kb bins mapped to the same 5kb bin
genes=read.table(paste0(indir,genename,".txt"),stringsAsFactors = F)
genes=genes$V1
dbin=read.table(paste0(indir,genename,"_bin1kbs_bin5kbs.txt"),stringsAsFactors = F)
#dtss=fread("Mus_musculus.GRCm38.74.TSS_trans2gene2bin_morebins.txt",stringsAsFactors = F)
#dtss1=dtss %>% group_by(bin1kb) %>% summarise(gene=paste0(gene,collapse = ";")) %>%as.data.frame()
dpeak=fread(paste0("marker_peaks_ds_",chr,".txt"),stringsAsFactors = F)
dgene=fread(paste0("Mus_musculus.GRCm38.74.TSS_trans2gene2bin_morebins_",chr,".txt"),stringsAsFactors = F)
hicdata=fread(paste0(chr,"_5kb_pairs_duan_scaled_1.txt"))
dhicall=NULL
for(i in 1:nrow(dbin)){
  gene=dbin$V1[i]
  bin5kb=dbin$V2[i]
  ii=which(dgene$bin1kb==gene)
  genetss=dgene$gene[ii]
  gene1kb=unlist(strsplit(as.character(gene),"_"))
  geneloc=(as.numeric(gene1kb[2])+as.numeric(gene1kb[3]))/2
  ## stability selection:
  infile=paste0(outdir,genename,"_",gene,"_consensus_EP_predictions.txt")
  if (!file.exists(infile)){
    next
  }
  pred=read.table(infile)
  names(pred)[1:3]=c("interact_bin1kb","promoter_bin1kb","pred")
  pred$genetss=genetss
  bin1=do.call(rbind,strsplit(as.character(pred$interact_bin1kb),"_"))
  bin1start=as.integer(as.numeric(bin1[,2])/5000)*5000
  mid1kb=(as.numeric(bin1[,2])+as.numeric(bin1[,3]))/2
  pred$distance=mid1kb-geneloc
  pred$pair1kb=paste0(pred$interact_bin1kb,"-",gene)
  id0=which(pred$distance>0)
  pred$pair1kb[id0]=paste0(gene,"-",pred$interact_bin1kb[id0])
  pred$interact_bin5kb=paste0(chr,"_",bin1start,"_",bin1start+5000)
  pred$promoter_bin5kb=bin5kb
  peak=unlist(strsplit(as.character(bin5kb),"_"))
  peakloc=(as.numeric(peak[2])+as.numeric(peak[3]))/2
  pred$pair5kb=paste0(pred$interact_bin5kb,"-",bin5kb)
  id=which(bin1start>peakloc)
  pred$pair5kb[id]=paste0(bin5kb,"-",pred$interact_bin5kb[id])
  pred1=merge(pred,dpeak,by.x="interact_bin1kb",by.y="bin1kb",all.x=T)
  pred2=merge(pred1,hicdata,by.x="pair5kb",by.y="Pair",all.x=T)
  #write.table(pred2,file=paste0(outdir,genename,"_",gene,"_bins_prediction.txt"),quote=F,sep="\t",row.names = F,col.names =T)
  dhicall=rbind(dhicall,pred2)
}
fwrite(dhicall,paste0(outdir,genename,"_hicCount_confidence_predictions_beforemerging.txt"),quote=F,sep="\t",row.names = F,col.names = T)
#dff=dhicall %>% group_by(pair1kb) %>%  summarise(pred=max(pred),distance=max(distance),pair5kb=unique(pair5kb),count=unique(Cnt),BinomQ=unique(BinomQ),gene=unique(gene)) %>% as.data.frame()
#dff=dff[order(dff$pred,decreasing = T),]
#dff$rank=rank(desc(dff$pred))
#cfwrite(dff,paste0(outdir,genename,"_hicCount_confidence_predictions.txt"),quote=F,sep="\t",row.names = F,col.names = T)
