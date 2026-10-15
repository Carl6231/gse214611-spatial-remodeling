# Read completed analysis outputs; never invoke analysis/model functions.
root<-'/home/zhenhanli/myocardial/scop/GSE214611_spatial'
script_dir_collect<-file.path(root,'22_JTM_completion_20261014/figures_reviewed/script')
dest<-file.path(root,'23_figure_structure_20261015/source_data');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
source(file.path(script_dir_collect,'figure_remaining_helpers.R'))
save_orig<-save_spatial_asset_v2
save_spatial_asset_v2<-function(plot,id,width_mm,height_mm){f<-file.path(ASSET_V2,paste0(id,'.png'));if(file.exists(f))return(f);save_orig(plot,id,width_mm,height_mm)}
export_supplement<-function(pages,id,height=225){saveRDS(list(pages=pages,height=height),file.path(dest,paste0(id,'_source_plots.rds')))}
export_remaining<-function(p,id,section='main',height=225){saveRDS(p,file.path(dest,paste0(id,'_source_plot.rds')))}
evaluate<-function(name,cut=NULL){s<-readLines(file.path(script_dir_collect,name));s<-s[!grepl('^source\\(',s)];if(!is.null(cut))s<-s[seq_len(which(grepl(cut,s,fixed=TRUE))[1]-1L)];eval(parse(text=s),envir=.GlobalEnv)}
# Ready-made SCOP panels for image pairs; retain their original pixels.
if(!file.exists(file.path(dest,'Figure_1_source_panels.rds'))) file.copy(file.path(root,'22_JTM_completion_20261014/figures_reviewed/build/Figure_1_refined/panels.rds'),file.path(dest,'Figure_1_source_panels.rds'),overwrite=FALSE)
if(!file.exists(file.path(dest,'Figure_2_values.rds'))){evaluate('render_figure2_v2.R');saveRDS(list(composition=composition,burden=burden,agreement=agreement,markers=markers,uncertain=uncertain),file.path(dest,'Figure_2_values.rds'))}
evaluate('finalize_figure5_6.R','for(n in 5:6)')
saveRDS(list(physical=physical,genes=gene_panel,pathways=path_df,joint=joint,models=model_source,meringue=mer_matrix,axes=axis_plot_data,oracle=oracle,continuum_limit=day7_lim,gene_limits=gene_limits),file.path(dest,'Figure_5_6_values.rds'))
for(f in c(continuum_assets,unname(gene_assets)))file.copy(f,file.path(dest,basename(f)),overwrite=FALSE)
# Extract the saved expression and corrected-niche coordinates once, one tissue at a time.
spatial<-lapply(SAMPLE,function(k){o<-readRDS(file.path(root,'objects/phase_09/niches',paste0(k,'_phase09_niche.rds')))
 xy<-Seurat::GetTissueCoordinates(o,scale=NULL);co<-as.data.table(xy,keep.rownames='barcode')
 print(c(k,nrow(co)));out<-as.data.table(o@meta.data,keep.rownames='barcode');out[,sample_key:=k]
 fields<-c('barcode','sample_key',grep('^P09N_',names(out),value=TRUE));out<-out[,..fields]
 out<-merge(out,co,by='barcode',sort=FALSE)
 if(k%in%d7){for(g in c('Sdf4','Abca1'))out[[g]]<-get_gene_values(k,g)[out$barcode]}
 out})
saveRDS(rbindlist(spatial,fill=TRUE),file.path(dest,'niche_spatial_values.rds'))
# Retain the complete original plot objects for rearranging the contextual and activity figures.
evaluate('refine_statial.R')
saveRDS(list(mi=mi,sel=sel,coverage=cv,maps=maps),file.path(dest,'Figure_7_values.rds'))
# The SecAct script performs plotting and summaries only; all model results have _SUCCESS.
s<-readLines(file.path(script_dir_collect,'render_secact.R'));s<-s[!grepl('^source\\(',s)];s<-s[seq_len(which(grepl("classes<-unique",s,fixed=TRUE))[1]-1L)];eval(parse(text=s),envir=.GlobalEnv)
saveRDS(list(sn=sn,st=st,cc=cc,meta=meta,v=v,sn_ids=sn_ids,st_ids=st_ids,labels_sn=labels_sn,labels_st=labels_st),file.path(dest,'Figure_8_values.rds'))
# Numerical source tables needed for the new overview and paired comparisons.
paths<-c(list.files(file.path(root,'results/figures/Figure_1/source_data'),full.names=TRUE),file.path(root,'metadata',c('phase_10_spot_signature_scores.tsv.gz','phase_09_niche_parameter_summary.tsv','phase_09_niche_grid_slice_summary.tsv','phase_09_niche_spatial_null.tsv','phase_09_niche_continuous_slice_summary.tsv','phase_09_niche_leave_one_slice_out.tsv')))
for(f in paths)if(file.exists(f))file.copy(f,file.path(dest,basename(f)),overwrite=FALSE)
cat('SOURCE_EXPORT_DONE\n')
