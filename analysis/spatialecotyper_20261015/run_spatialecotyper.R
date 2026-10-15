#!/usr/bin/env Rscript
options(warn=1)
suppressPackageStartupMessages({library(scop);library(Seurat);library(SpatialEcoTyper);library(NMF);library(Matrix);library(data.table)})
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
root<-normalizePath(args[1]);out<-normalizePath(args[2])
for(d in c('objects','tables','figures'))dir.create(file.path(out,d),showWarnings=FALSE)
stopifnot('RunSpatialEcoTyper'%in%getNamespaceExports('scop'))
data.table::setDTthreads(1)
samples<-c('Sham_Day3_rep1','Sham_Day3_rep2',paste0('MI_Day3_rep',1:3),paste0('MI_Day7_rep',1:3),'MI_1HR_rep1','MI_4HR_rep1','NeedlePass_Day3_rep3','NeedlePass_Day3_rep4')
wpath<-system.file('extdata','Bulk_SE_Recovery_W.rds',package='SpatialEcoTyper');W<-readRDS(wpath)
genes<-unique(sub('__.*','',rownames(W)))
stopifnot(length(genes)==499L,ncol(W)==10L)
mapfun<-getFromNamespace('mouse2human_mat','SecAct')
mapping<-read.csv(system.file('extdata','Mouse2Human_filter.csv',package='SecAct'))
fwrite(mapping,file.path(out,'tables/mouse_human_mapping.csv'))
saveRDS(W,file.path(out,'objects/tumor_reference_W.rds'))
cat('START',format(Sys.time(),tz='UTC',usetz=TRUE),'\nSCOP',as.character(packageVersion('scop')),'SpatialEcoTyper',as.character(packageVersion('SpatialEcoTyper')),'\n')
cat('MODE deconvolute; native tumor reference; mouse-human ortholog projection; no cardiac ecotype discovery\n')
status<-list();combined<-list();coverage<-list()
for(k in samples){
 t0<-proc.time()[3];dest<-file.path(out,'objects',paste0(k,'_SpatialEcoTyper.rds'))
 if(file.exists(dest)){
  z<-readRDS(dest);stopifnot(identical(z$sample,k),identical(z$parameters$mode,'deconvolute'))
  cat('REUSE',k,'\n')
 }else{
  cat('RUN',k,format(Sys.time(),tz='UTC',usetz=TRUE),'\n');flush.console()
  o<-readRDS(file.path(root,'objects/phase_09/deconvolution',paste0(k,'_phase09_expanded.rds')))
  counts<-LayerData(o,assay='Spatial',layer='counts')
  stopifnot(all(counts@x>=0),all(abs(counts@x-round(counts@x))<1e-8))
  human<-mapfun(counts);stopifnot(identical(colnames(human),colnames(o)),!anyDuplicated(rownames(human)),all(Matrix::colSums(human)>0))
  present<-intersect(genes,rownames(human));expressed<-present[Matrix::rowSums(human[present,,drop=FALSE])>0]
  gene_coverage<-data.table(sample=k,gene=genes,mapped=genes%in%present,expressed=genes%in%expressed)
  # Full mapped transcriptome determines library sizes before selecting model genes.
  norm<-NormalizeData(human,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
  pred<-scop::RunSpatialEcoTyper(srt=norm,mode='deconvolute',W=W,scale=TRUE,nsample.per.run=500,sum2one=TRUE,ncores=2,verbose=TRUE)
  stopifnot(is.matrix(pred),setequal(rownames(pred),colnames(o)),identical(colnames(pred),colnames(W)),all(is.finite(pred)),all(pred>=0),max(abs(rowSums(pred)-1))<1e-8)
  pred<-pred[colnames(o),,drop=FALSE]
  coords<-GetTissueCoordinates(o);stopifnot(all(colnames(o)%in%rownames(coords)));coords<-coords[colnames(o),,drop=FALSE]
  z<-list(sample=k,abundance=pred,coordinates=coords,gene_coverage=gene_coverage,
    parameters=list(mode='deconvolute',reference='SpatialEcoTyper packaged tumor Bulk_SE_Recovery_W',mapping='SecAct Mouse2Human_filter.csv',normalization='log1p CP10K on full human-mapped counts',scale=TRUE,nsample.per.run=500,sum2one=TRUE,ncores=2,native_seed=39),
    scientific_status='EXPLORATORY_CROSS_SPECIES_TUMOR_REFERENCE_PROJECTION_NOT_CARDIAC_ECOTYPE_DISCOVERY',elapsed_seconds=unname(proc.time()[3]-t0))
  saveRDS(z,dest)
  rm(o,counts,human,norm,pred);gc(verbose=FALSE)
 }
 stopifnot(all(is.finite(z$abundance)),max(abs(rowSums(z$abundance)-1))<1e-8)
 status[[k]]<-data.table(sample=k,status='COMPLETE',spots=nrow(z$abundance),components=ncol(z$abundance),mapped_genes=sum(z$gene_coverage$mapped),expressed_genes=sum(z$gene_coverage$expressed),reference_genes=length(genes),elapsed_seconds=z$elapsed_seconds,mode='deconvolute',interpretation=z$scientific_status)
 combined[[k]]<-data.table(sample=k,barcode=rownames(z$abundance),as.data.frame(z$abundance))
 coverage[[k]]<-z$gene_coverage
 fwrite(rbindlist(status),file.path(out,'tables/sample_status.tsv'),sep='\t')
 cat('DONE',k,'spots',nrow(z$abundance),'mapped',sum(z$gene_coverage$mapped),'/',length(genes),'seconds',round(z$elapsed_seconds,1),'\n');flush.console()
}
a<-rbindlist(combined);stopifnot(nrow(a)==25548L,length(status)==12L,!anyDuplicated(a[,.(sample,barcode)]))
fwrite(a,file.path(out,'tables/spot_abundances.tsv.gz'),sep='\t',compress='gzip')
fwrite(rbindlist(coverage),file.path(out,'tables/model_gene_coverage.tsv'),sep='\t')
l<-melt(a,id.vars=c('sample','barcode'),variable.name='component',value.name='abundance')
sumtab<-l[,.(mean_abundance=mean(abundance),median_abundance=median(abundance),q25=quantile(abundance,.25),q75=quantile(abundance,.75),spots=.N),by=.(sample,component)]
fwrite(sumtab,file.path(out,'tables/tissue_summary.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines('COMPLETE: 12 real Visium tissues, 25548 spots; cross-species tumor-reference projection only',file.path(out,'_SUCCESS'))
cat('ALL_COMPLETE',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
