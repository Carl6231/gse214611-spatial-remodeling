args<-commandArgs(TRUE);stopifnot(length(args)==1L)
root<-normalizePath(args[1]);base<-file.path(root,'22_JTM_completion_20261014')
.libPaths(c(file.path(root,'21_scop_expansion/_runtime/spacet/R-4.5/site-library'),.libPaths()))
suppressPackageStartupMessages(library(SpaCET))
prior<-file.path(base,'P06_SECACT_compat2');out<-file.path(base,'P06_SECACT_parallel')
suppressPackageStartupMessages(library(data.table))
sn<-c('SRR21777596','SRR21777597','SRR21777598','SRR21777581','SRR21777582','SRR21777599')
st<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3));ids<-c(sn,st)
pstat<-fread(file.path(prior,'functions_final_status.tsv'))
old_complete<-pstat[,.(complete=.N==ifelse(startsWith(sample_id[1],'SRR'),2L,15L) && all(status%in%c('PASS','NEGATIVE_COMPLETE','NOT_APPLICABLE'))),by=sample_id][complete==TRUE,sample_id]
source_for<-setNames(vapply(ids,function(id)if(id%in%old_complete)prior else file.path(out,'jobs',id),character(1)),ids)
for(id in setdiff(ids,old_complete))stopifnot(file.exists(file.path(source_for[id],'_SUCCESS')))
status<-rbindlist(lapply(ids,function(id){d<-fread(file.path(source_for[id],'functions_final_status.tsv'))[sample_id==id];d[,source_run:=unname(source_for[id])];d}))
stopifnot(nrow(status)==102L,all(status$status%in%c('PASS','NEGATIVE_COMPLETE','NOT_APPLICABLE')),
 !anyDuplicated(status[,.(function_name=get('function'),sample_id,target,signal_mode)]),
 identical(as.integer(table(factor(status[['function']],levels=c('RunSecAct','RunSecActCCC','RunSecActSignalingPattern','RunSecActPatternGenes','RunSecActVelocity')))),c(12L,6L,6L,6L,72L)))
pr<-fread(file.path(prior,'parameters.tsv'))[,.(parameter,value)]
for(s in unique(source_for))stopifnot(identical(pr,fread(file.path(s,'parameters.tsv'))[,.(parameter,value)]))
write_gz<-function(d,path){con<-gzfile(path,'wt');on.exit(close(con));write.table(d,con,sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
for(f in c('secact_activity.tsv.gz','secact_ccc.tsv','secact_patterns.tsv','secact_pattern_genes.tsv','secact_velocity_sample_effects.tsv')){
 d<-rbindlist(lapply(ids,function(id){z<-fread(file.path(source_for[id],'tables',f));z[sample_id==id]}),fill=TRUE)
 stopifnot(nrow(d)>0)
 if(grepl('gz$',f))write_gz(d,file.path(out,'tables',f))else fwrite(d,file.path(out,'tables',f),sep='\t')
}
a<-fread(file.path(out,'tables/secact_activity.tsv.gz'));v<-fread(file.path(out,'tables/secact_velocity_sample_effects.tsv'))
a_summary<-a[,.(n_biological_samples=uniqueN(sample_id),median_effect=median(zscore),IQR_effect=IQR(zscore),n_positive=sum(zscore>0),n_negative=sum(zscore<0)),by=.(modality,stage,feature=secretedProtein,context=celltype_or_scope)]
a_summary[,`:=`(record_type='activity_zscore',inference='descriptive_sample_level_no_spot_or_nucleus_P_value')]
# The native backend rescales x and y separately. Preserve these as display diagnostics,
# not signal strengths or comparable physical displacements.
v_unique<-unique(v[,.(sample_id,stage,ligand_secact_human,signal_mode,n_vectors,median_vector_length,IQR_vector_length,median_dx,median_dy,status)])
stopifnot(nrow(v_unique)==72L)
coverage<-merge(v_unique,unique(a[modality=='Visium',.(sample_id,n_spots=n_input_units)]),by='sample_id')
coverage[,vector_fraction:=n_vectors/n_spots]
stopifnot(all(coverage$vector_fraction>=0 & coverage$vector_fraction<=1))
coverage[,metric_interpretation:='coverage_is_computational_evaluability;length_dx_dy_are_independently_rescaled_native_display_geometry']
fwrite(coverage,file.path(out,'tables/secact_vector_coverage.tsv'),sep='\t')
v_summary<-coverage[,.(n_biological_samples=uniqueN(sample_id),median_effect=median(vector_fraction),IQR_effect=IQR(vector_fraction),n_positive=sum(vector_fraction>0),n_negative=0L),by=.(stage,feature=ligand_secact_human,context=signal_mode)]
v_summary[,`:=`(record_type='fraction_spots_with_native_vectors',modality='Visium',inference='descriptive_tissue_level_computational_evaluability_not_signal_strength')]
fwrite(rbindlist(list(a_summary,v_summary),use.names=TRUE),file.path(out,'tables/secact_discovery_summary.tsv'),sep='\t')
fwrite(status,file.path(out,'functions_final_status.tsv'),sep='\t');fwrite(pr,file.path(out,'parameters.tsv'),sep='\t')
file.copy(file.path(prior,'tables/ortholog_mapping.tsv'),file.path(out,'tables/ortholog_mapping.tsv'),overwrite=FALSE)
file.copy(file.path(prior,'tables/pilot_output_audit.tsv'),file.path(out,'tables/pilot_output_audit.tsv'),overwrite=FALSE)
for(id in ids)for(folder in c('objects','figures')){
 files<-list.files(file.path(source_for[id],folder),full.names=TRUE)
 files<-files[startsWith(basename(files),paste0(id,'_'))|startsWith(basename(files),paste0('SECACT_',id,'_'))]
 for(f in files){
  target<-file.path(out,folder,basename(f))
  if(!file.exists(target)){
   if(endsWith(f,'_velocity.rds')){
    x<-readRDS(f)
    if(!is.null(x$native_plot)){x$spots<-as.data.frame(x$native_plot$data);x$native_plot<-NULL}
    x$format<-'compact_native_display_source_v1';x$source_rds<-f
    stopifnot(is.data.frame(x$spots),all(c('x','y','value')%in%names(x$spots)))
    saveRDS(x,target);rm(x);gc(verbose=FALSE)
   }else stopifnot(file.copy(f,target))
  }
 }
}
expected_obj<-c(file.path(out,'objects',paste0(sn,'_snRNA_RunSecAct.rds')),file.path(out,'objects',paste0(st,'_SpaCET_SecAct.rds')))
stopifnot(all(file.exists(expected_obj)),length(list.files(file.path(out,'figures'),pattern='[.]pdf$'))==72L)
# Read the model objects once and supply the numerical source data behind the new figures.
pattern_source<-list()
for(id in st){
 z<-readRDS(file.path(out,'objects',paste0(id,'_SpaCET_SecAct.rds')))
 stopifnot(inherits(z,'SpaCET'),all(is.finite(z@results$SecAct_output$SecretedProteinActivity$zscore)),ncol(z@results$SecAct_output$pattern$weight.W)%in%2:5)
 H<-z@results$SecAct_output$pattern$signal.H
 coords<-z@input$spotCoordinates[colnames(H),,drop=FALSE]
 stopifnot(identical(rownames(coords),colnames(H)))
 pattern_source[[id]]<-data.table(sample_id=id,pattern_index=rep(seq_len(nrow(H)),each=ncol(H)),barcode=rep(colnames(H),times=nrow(H)),
  x_um=rep(coords$coordinate_x_um,times=nrow(H)),y_um=rep(coords$coordinate_y_um,times=nrow(H)),loading=as.vector(t(H)))
}
write_gz(rbindlist(pattern_source),file.path(out,'tables/secact_pattern_spot_loadings.tsv.gz'))
vector_source<-list();background_source<-list()
for(id in st)for(g in unique(v$ligand_secact_human))for(mode in c('sending','receiving')){
 x<-readRDS(file.path(out,'objects',paste0(id,'_',g,'_',mode,'_velocity.rds')))
 key<-paste(id,g,mode,sep='|')
 seg<-as.data.table(x$segments);seg[,`:=`(sample_id=id,ligand=g,signal_mode=mode)]
 vector_source[[key]]<-seg
 bg<-as.data.table(x$spots);bg[,`:=`(sample_id=id,ligand=g,signal_mode=mode)]
 background_source[[key]]<-bg
}
write_gz(rbindlist(vector_source,fill=TRUE),file.path(out,'tables/secact_vector_fields_native_display.tsv.gz'))
write_gz(rbindlist(background_source,fill=TRUE),file.path(out,'tables/secact_vector_background_native_display.tsv.gz'))
fwrite(data.table(status='PASS',samples_completed=12L,producer_records=102L,vector_fields=72L,
 biological_interpretation='descriptive_model_outputs_not_causal_validation'),file.path(out,'STATUS.tsv'),sep='\t')
for(s in unique(source_for))file.copy(file.path(s,'sessionInfo.txt'),file.path(out,'logs',paste0('sessionInfo_analysis_',basename(s),'.txt')),overwrite=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
file.create(file.path(out,'_SUCCESS'));cat('SECACT_ALL_COMPLETE 12 samples; 102 producer records; 72 vector fields\n')
