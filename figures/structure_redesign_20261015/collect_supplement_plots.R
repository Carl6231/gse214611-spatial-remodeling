root_collect<-'/home/zhenhanli/myocardial/scop/GSE214611_spatial'
script_dir_collect<-file.path(root_collect,'22_JTM_completion_20261014/figures_reviewed/script')
dest_collect<-file.path(root_collect,'23_figure_structure_20261015/source_data/supplement_plots');dir.create(dest_collect,recursive=TRUE,showWarnings=FALSE)
jobs<-list(c('render_supplement_atlases_v2.R','Figure_S1,Figure_S3,Figure_S7'),c('render_supplement_qc_v2.R','Figure_S2'),c('render_supplement_giotto_v2.R','Figure_S4,Figure_S5'),c('render_reference_v2.R','Figure_S6'),c('render_support_atlases.R',''),c('render_supplement_spatial_v2.R','Figure_S9,Figure_S10,Figure_S11'),c('render_supplement_boundary_v2.R','Figure_S12,Figure_S13'),c('render_original_sensitivity.R','Figure_S14'),c('finalize_s15.R','Figure_S15'),c('render_secact.R','Figure_S16'))
for(job in jobs){
 message('COLLECT ',job[1]);local({
  source(file.path(script_dir_collect,'figure_remaining_helpers.R'),local=TRUE)
  requested_figures<-strsplit(job[2],',',fixed=TRUE)[[1]]
  wants<-function(id)!length(requested_figures)||id%in%requested_figures
  original_save<-save_spatial_asset_v2
  save_spatial_asset_v2<-function(plot,id,width_mm,height_mm){f<-file.path(ASSET_V2,paste0(id,'.png'));if(file.exists(f))return(f);original_save(plot,id,width_mm,height_mm)}
  export_supplement<-function(pages,id,height=225){saveRDS(list(pages=pages,height=height),file.path(dest_collect,paste0(id,'.rds')));cat('COLLECTED',id,'pages',length(pages),'\n')}
  export_remaining<-function(...)invisible(NULL)
  s<-readLines(file.path(script_dir_collect,job[1]));s<-s[!grepl('^source\\(',s)]
  eval(parse(text=s),envir=environment())
 });gc()
}
cat('SUPPLEMENT_SOURCE_DONE\n')
