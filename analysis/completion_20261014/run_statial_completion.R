# Explicit amendment: parent is a set of joint cell-class@context labels,
# never a metadata column name. All 12 fixed directions, r=200 um retained.
# Spot-level contextual L is descriptive; tissue is the replicate unit.
args <- commandArgs(TRUE); stopifnot(length(args)==2L)
root <- normalizePath(args[1]); out <- args[2]
.libPaths(c(file.path(root,'21_scop_expansion/_runtime/statial/R-4.5/site-library'),.libPaths()))
suppressPackageStartupMessages({library(Seurat);library(scop);library(Statial);library(data.table);library(Matrix)})
stopifnot(as.character(packageVersion('Statial'))=='1.12.0',as.character(packageVersion('scop'))=='0.9.1')
for(d in c('tables','objects','logs'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
writeLines(c('Context amendment: cell-class@zone labels with parent equal to all labels in that zone.',
 'All physical-tissue coordinates remain in each call; r=200 micrometres; convex window expanded by relative 1e-8 for floating-point boundary inclusion; border edge correction.',
 '12 originally fixed directions; no outcome-driven pair or radius selection.',
 'Dominant RCTD8 mixture labels are not single-cell identities. Descriptive per-tissue estimates; no spot-level P values.',
 'External BANKSY domains are sample-specific and cannot be matched numerically between tissues.'),file.path(out,'parameters.txt'))
# Local numerical compatibility: expand the convex hull by 1e-8 about its centre.
# This is sub-micrometre padding, preventing boundary points from being dropped.
statial_env <- new.env(parent=asNamespace('Statial'))
statial_env$makeWindow_safe <- function(data,window,window.length) {
 w<-Statial::makeWindow(data,window,window.length)
 ctr<-c(mean(w$xrange),mean(w$yrange))
 w<-spatstat.geom::affine(w,mat=diag(1+1e-8,2),vec=-ctr*1e-8)
 stopifnot(all(spatstat.geom::inside.owin(data$x,data$y,w)))
 w
}
f<-getFromNamespace('Kontextual','Statial')
btxt<-paste(deparse(body(f)),collapse='\n')
stopifnot(length(gregexpr('Statial::makeWindow',btxt,fixed=TRUE)[[1]])==1L)
body(f)<-parse(text=gsub('Statial::makeWindow','makeWindow_safe',btxt,fixed=TRUE))[[1]]
environment(f)<-statial_env
Kontextual_safe<-f
scop_env<-new.env(parent=asNamespace('scop'))
scop_env$get_namespace_fun<-function(pkg,name) {
 if(pkg=='Statial' && name=='Kontextual')return(Kontextual_safe)
 get('get_namespace_fun',envir=asNamespace('scop'),inherits=TRUE)(pkg,name)
}
RunStatial_safe<-getFromNamespace('RunStatialKontextual','scop');environment(RunStatial_safe)<-scop_env
# A one-context encoding must reproduce the ordinary-label result exactly.
g <- expand.grid(x=seq(0,500,100),y=seq(0,500,100));g$imageID<-'test';g$cellType<-rep(c('A','B'),18)
a <- Kontextual_safe(g,r=200,from='A',to='B',parent=c('A','B'),cores=1)
g$cellType<-paste0(g$cellType,'@zone')
b <- Kontextual_safe(g,r=200,from='A@zone',to='B@zone',parent=c('A@zone','B@zone'),cores=1)
stopifnot(isTRUE(all.equal(a$kontextual,b$kontextual)),is.finite(b$kontextual))
cat('Joint-context label equivalence check PASS\n')
pairs <- fread(file.path(root,'21_scop_expansion/P02/tables/preregistered_context_pairs.tsv'))
stopifnot(nrow(pairs)==12L)
contracts <- rbind(fread(file.path(root,'21_scop_expansion/P01/tables/discovery_object_contract.tsv')),
                   fread(file.path(root,'21_scop_expansion/P01/tables/external_object_contract.tsv')))
# Pilot first; then sequential processing, releasing expression objects immediately.
contracts <- contracts[order(biological_sample_id!='MI_D3_rep1')]
all <- list()
for(i in seq_len(nrow(contracts))) {
 c <- contracts[i];id<-c$biological_sample_id
 path<-file.path(out,'objects',paste0(id,'_Statial.rds'))
 if(file.exists(path)){all[[id]]<-readRDS(path)$results;next}
 o<-readRDS(c$object_file);meta<-o[[]];stopifnot(nrow(meta)==c$n_spots)
 props<-strsplit(c$RCTD8_columns,';',fixed=TRUE)[[1]]
 m<-as.matrix(meta[,props]);valid<-apply(is.finite(m)&m>=0,1,all)&rowSums(m,na.rm=TRUE)>0
 cellclass<-rep('Unassigned',nrow(m))
 cellclass[valid]<-sub('RCTD8_prop_','',props)[max.col(m[valid,,drop=FALSE],ties.method='first')]
 if('RCTD8_dominant_type'%in%names(meta))stopifnot(all(cellclass[valid]==as.character(meta$RCTD8_dominant_type[valid])))
 if(!is.na(c$image)) {
  xy<-scop:::spatial_analysis_coords(o,image=c$image,coordinate_space='raw')$data
  stopifnot(identical(rownames(xy),rownames(meta)))
  x<-xy$x;y<-xy$y;context<-as.character(meta[[c$coarse_corridor_region_column]])
  context_source<-c$coarse_corridor_region_column
 } else {
  x<-meta[[sub('metadata:','',c$raw_x_column)]];y<-meta[[sub('metadata:','',c$raw_y_column)]]
  sc<-fread(c$banksy_sidecar_file);stopifnot(!anyDuplicated(sc$barcode))
  at<-match(rownames(meta),sc$barcode);stopifnot(!anyNA(at))
  context<-as.character(sc$banksy_domain[at]);context_source<-'sample_specific_BANKSY_domain'
 }
 stopifnot(length(context)==nrow(meta),!anyNA(context),all(nzchar(context)))
 smallmeta<-data.frame(sample_id=id,x_um=x*c$pixel_to_micron,y_um=y*c$pixel_to_micron,
      cellclass=cellclass,context=context,joint=paste(cellclass,context,sep='@'),row.names=rownames(meta))
 rm(o,meta,m);gc(verbose=FALSE)
 # Metadata-only Seurat avoids retaining expression matrices for a point-pattern analysis.
 dummy<-Matrix(1,nrow=1,ncol=nrow(smallmeta),sparse=TRUE,dimnames=list('metadata_only',rownames(smallmeta)))
 thin<-CreateSeuratObject(dummy,meta.data=smallmeta)
 coords<-scop:::spatial_analysis_coords(thin,coord.cols=c('x_um','y_um'),coordinate_space='raw')$data
 stopifnot(identical(coords$cell_id,rownames(smallmeta)),isTRUE(all.equal(coords$x,smallmeta$x_um)),isTRUE(all.equal(coords$y,smallmeta$y_um)))
 spec<-rbindlist(lapply(sort(unique(context)),function(z){
  p<-copy(pairs[,.(relation_family,direction_id,from,to)])
  p[,`:=`(zone=z,from_joint=paste(from,z,sep='@'),to_joint=paste(to,z,sep='@'))]
  p[,`:=`(n_from=vapply(from_joint,function(v)sum(smallmeta$joint==v),integer(1)),n_to=vapply(to_joint,function(v)sum(smallmeta$joint==v),integer(1)))]
  p
 }))
 parentdf<-data.frame(from=spec$from_joint,to=spec$to_joint)
 parentdf$parent<-I(lapply(spec$zone,function(z)sort(unique(smallmeta$joint[smallmeta$context==z]))))
 ans<-RunStatial_safe(thin,group.by='joint',r=200,parent_df=parentdf,sample.by='sample_id',
  coord.cols=c('x_um','y_um'),coordinate_space='raw',cores=1,inhom=FALSE,edge_correct=TRUE,window='convex',store_input=TRUE)
 raw<-as.data.table(ans@tools$StatialKontextual$table)
 spec[,test:=paste(from_joint,to_joint,sep='__')]
 stopifnot(nrow(raw)==nrow(spec),!anyDuplicated(raw$test),all(raw$r==200),all(raw$imageID==id))
 res<-merge(spec,raw,by='test',sort=FALSE)
 res[,`:=`(sample_id=id,n_unassigned=sum(cellclass=='Unassigned'),cohort=c$cohort,stage=c$timepoint,condition=c$condition,analysis_role=c$analysis_role,
  context_source=context_source,n_context=vapply(zone,function(z)sum(smallmeta$context==z),integer(1)),radius_um=200)]
 res[,status:=ifelse(n_from==0|n_to==0,'NOT_EVALUABLE_ABSENT_CLASS',ifelse(is.finite(kontextual),'EVALUABLE','NOT_EVALUABLE_NONFINITE'))]
 res[,interpretation:='descriptive_mixed_spot_contextual_association']
 saveRDS(list(results=res,backend=ans@tools$StatialKontextual,input_source=c$object_file),path)
 stopifnot(nrow(readRDS(path)$results)==nrow(spec))
 all[[id]]<-res
 fwrite(rbindlist(all),file.path(out,'tables/per_sample_kontextual.tsv'),sep='\t')
 cat('STATIAL_DONE',id,'finite=',sum(is.finite(res$kontextual)),'/',nrow(res),'\n')
 rm(ans,thin);gc(verbose=FALSE)
}
d<-rbindlist(all)
# Cross-tissue summaries only for common anatomical context labels in GSE214611.
sumtab<-d[cohort=='GSE214611',.(n_tissues=.N,n_evaluable=sum(is.finite(kontextual)),
 median_L=if(any(is.finite(kontextual)))median(kontextual[is.finite(kontextual)])else NA_real_,
 n_positive=sum(kontextual>0,na.rm=TRUE),n_negative=sum(kontextual<0,na.rm=TRUE)),by=.(condition,stage,context_source,zone,relation_family,direction_id)]
fwrite(sumtab,file.path(out,'tables/context_pair_stage_summary.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
stopifnot(length(unique(d$sample_id))==22L)
writeLines('Completed 22 physical-tissue runs; non-evaluable relationships retained explicitly.',file.path(out,'_SUCCESS'))
