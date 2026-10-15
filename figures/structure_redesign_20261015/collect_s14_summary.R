root<-'/home/zhenhanli/myocardial/scop/GSE214611_spatial';dir<-file.path(root,'22_JTM_completion_20261014/figures_reviewed/script')
source(file.path(dir,'figure_common_v2.R'))
oldsave<-save_spatial_asset_v2
save_spatial_asset_v2<-function(plot,id,width_mm,height_mm){f<-file.path(ASSET_V2,paste0(id,'.png'));if(file.exists(f))return(f);oldsave(plot,id,width_mm,height_mm)}
Sys.setenv(REDESIGN_V2_BUILD_ROOT=file.path(root,'22_JTM_completion_20261014/figures_reviewed/build'))
s<-readLines(file.path(dir,'render_sensitivity_S14.R'));s<-s[!grepl('^source\\(|^if \\(file.exists\\(STAGE_V2\\)\\) stop',s)];s<-s[seq_len(which(grepl('paths <- export_page_v2',s,fixed=TRUE))[1]-1L)];eval(parse(text=s))
saveRDS(list(pages=list(figure4),height=225),file.path(root,'23_figure_structure_20261015/source_data/supplement_plots/Figure_S14_summary.rds'))
