#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
stopifnot(!any(file.exists(file.path(FIG_V2,c('Figure_2.pdf','Figure_2.svg')))))
stopifnot('SpatialDeconvolutionPlot' %in% getNamespaceExports('scop'))
src <- file.path(ROOT,'results','figures','Figure_2','source_data')
composition <- read_table(file.path(src,'Figure_2f_slice_celltype_composition_row_zscore.tsv'))
burden <- read_table(file.path(src,'Figure_2g_slice_macrophage_burden.tsv'))
agreement <- read_table(file.path(src,'Figure_2i_rctd_spatialdwls_agreement.tsv'))
markers <- read_table(file.path(src,'Figure_2j_marker_concordance.tsv'))
uncertain <- read_table(file.path(src,'Figure_2j_uncertain_burden.tsv'))
stopifnot(setequal(uncertain$sample_key,SAMPLE),sum(burden$spots_descriptive)==25548)
colours <- c('#FAFAFA','#EBC7D9','#CC79A7','#6F214F')
maps <- setNames(lapply(FORMAL,function(key) {
  message('Figure 2 spatial: ',key)
  cached <- file.path(ASSET_V2,paste0('F2_',key,'.png'))
  if(file.exists(cached)) return(cached)
  obj <- readRDS(file.path(ROOT,'objects','deconvolution','samples',paste0(key,'_phase04.rds')))
  p <- scop::SpatialDeconvolutionPlot(obj,tool_name='RCTD',cell_types='Macrophage',
    plot_type='point',combine=FALSE,image=Images(obj)[[1]],overlay_image=TRUE,
    image.alpha=.18,crop=TRUE,show_axes=FALSE,pt.size=.43,pt.alpha=.92,stroke=0,
    palcolor=colours,legend.position='none')[[1]]
  p <- spatial_leaf(p)
  p <- suppressMessages(p+scale_colour_gradientn(colours=colours,limits=c(0,1),oob=scales::squish))
  stopifnot(identical(as.numeric(p$scales$get_scales('colour')$limits),c(0,1)))
  out <- save_spatial_asset_v2(p,paste0('F2_',key),42,42)
  rm(obj,p);gc(verbose=FALSE);out
}),FORMAL)
sham <- wrap_plots(lapply(FORMAL[1:2],function(k) heading(img_panel_compact_v2(maps[[k]]),LABEL[k])),ncol=1)
mi3 <- wrap_plots(lapply(FORMAL[3:5],function(k) heading(img_panel_compact_v2(maps[[k]]),LABEL[k])),ncol=3)
mi7 <- wrap_plots(lapply(FORMAL[6:8],function(k) heading(img_panel_compact_v2(maps[[k]]),LABEL[k])),ncol=3)
atlas <- wrap_plots(list(sham,wrap_plots(list(mi3,mi7),ncol=1)),ncol=2,widths=c(1,3.08))
key<-wrap_elements(full=grobTree(
  rasterGrob(matrix(colorRampPalette(colours,space='Lab')(256),nrow=1),x=.5,y=.72,width=.98,height=unit(1.4,'mm')),
  textGrob('0',x=.01,y=.15,just='left',gp=gpar(fontfamily=FONT_V2,fontsize=6.5)),
  textGrob('RCTD macrophage fraction',x=.5,y=.15,gp=gpar(fontfamily=FONT_V2,fontsize=6.5)),
  textGrob('1',x=.99,y=.15,just='right',gp=gpar(fontfamily=FONT_V2,fontsize=6.5))))
atlas <- wrap_plots(list(atlas,key),ncol=1,heights=unit.c(unit(1,'null'),unit(6,'mm')))
a <- block(atlas,'a','Macrophage distribution across formal tissues')

bd <- burden[sample_key %in% FORMAL]
bd[, group:=factor(role(sample_key),levels=c('Sham','MI Day 3','MI Day 7'))]
bd[, x:=as.integer(group)+rep(c(-.09,0,.09),length.out=.N),by=group]
summary_points <- function(value,ylabel) {
  ggplot(bd,aes(x,.data[[value]],colour=group))+geom_point(size=1.65)+
    geom_segment(data=bd[,.(x=as.integer(group[1]),y=median(get(value))),by=group],
      aes(x=x-.2,xend=x+.2,y=y,yend=y),inherit.aes=FALSE,linewidth=.4,colour=COL_AXIS)+
    scale_colour_manual(values=role_cols_v2,guide='none')+
    scale_x_continuous(NULL,breaks=1:3,labels=c('Sham\nn=2','D3\nn=3','D7\nn=3'))+
    scale_y_continuous(ylabel,limits=c(0,NA),expand=expansion(mult=c(0,.12)))+quiet_theme
}
b <- block(summary_points('mean_RCTD_Macrophage','Mean fraction'),'b','Macrophage load')
c <- block(summary_points('macrophage_rich_fraction','Fraction of retained spots'),'c','Macrophage-rich area')
six <- c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage','Pericyte','Schwann_cell')
z <- composition[sample_key %in% FORMAL]
z[,sample:=factor(sample_key,levels=FORMAL)]
z[,cell:=factor(celltype,levels=rev(six))]
d <- ggplot(z,aes(sample,cell,fill=row_z_score))+geom_tile(colour='white',linewidth=.15)+
  scale_fill_gradient2(low=COL_CONTRACTILE,mid='white',high=COL_FIBROTIC,midpoint=0,
    limits=c(-2,2),oob=scales::squish,name='Row z')+
  scale_x_discrete(NULL,labels=short_label[FORMAL])+scale_y_discrete(NULL,labels=CLASS_LABEL)+
  quiet_theme+theme(axis.line=element_blank(),axis.ticks=element_blank(),
    axis.text.x=element_text(angle=45,hjust=1),legend.key.height=unit(2,'mm'))
d <- block(d,'d','Composition context')

others <- z[celltype %in% c('Cardiomyocyte','Fibroblast','Endothelial')]
others[,group:=factor(role(sample_key),levels=c('Sham','MI Day 3','MI Day 7'))]
others[,celltype:=factor(celltype,levels=c('Cardiomyocyte','Fibroblast','Endothelial'))]
others[,x:=as.integer(group)+seq(-.09,.09,length.out=.N),by=.(celltype,group)]
e <- ggplot(others,aes(x,mean_proportion,colour=group))+geom_point(size=1.6)+
  geom_segment(data=others[,.(x=as.integer(group[1]),y=median(mean_proportion)),by=.(celltype,group)],
    aes(x=x-.2,xend=x+.2,y=y,yend=y),inherit.aes=FALSE,colour=COL_AXIS,linewidth=.4)+
  facet_wrap(~celltype,nrow=1,scales='free_y')+scale_colour_manual(values=role_cols_v2,guide='none')+
  scale_x_continuous(NULL,breaks=1:3,labels=c('Sham','Day 3 MI','Day 7 MI'))+
  scale_y_continuous(NULL,breaks=scales::breaks_pretty(n=2),limits=c(0,NA),expand=expansion(mult=c(0,.12)))+quiet_theme
e <- block(e,'e','Other core components · mean RCTD fraction, class-specific scales')

core <- six[1:4]
matrix_plot <- function(dt,value) {
  d <- copy(dt[sample_key %in% FORMAL & celltype %in% core])
  stopifnot(nrow(d)==32L,!anyDuplicated(d[,.(sample_key,celltype)]))
  d[,sample:=factor(sample_key,levels=FORMAL)]
  d[,cell:=factor(celltype,levels=rev(core))]
  ggplot(d,aes(sample,cell,fill=.data[[value]]))+geom_tile(colour='white',linewidth=.2)+
    geom_text(aes(label=sprintf('%.2f',.data[[value]])),family=FONT_V2,size=6.5/GG_PT)+
    scale_fill_gradientn(colours=c('#F3F6FA','#CCDDEE','#80A8CE'),limits=c(-1,1),name='ρ')+
    scale_x_discrete(NULL,labels=short_label[FORMAL])+scale_y_discrete(NULL,labels=CLASS_LABEL)+
    quiet_theme+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(angle=45,hjust=1),legend.position='none')
}
f <- block(matrix_plot(agreement,'spearman_rctd_vs_dwls'),'f','RCTD–SpatialDWLS agreement')
g <- block(matrix_plot(markers,'rctd_marker_spearman'),'g','RCTD–marker concordance')
uncertain[,sample:=factor(sample_key,levels=SAMPLE)]
h <- ggplot(uncertain,aes(sample,uncertain_fraction))+geom_segment(aes(xend=sample,y=0,yend=uncertain_fraction),colour='#C5CBD1',linewidth=.45)+
  geom_point(aes(colour=role(sample_key)),size=1.6)+scale_colour_manual(values=role_cols_v2,guide='none')+
  scale_x_discrete(NULL,labels=short_label[SAMPLE])+scale_y_continuous('Uncertain fraction',limits=c(0,NA),expand=expansion(mult=c(0,.1)))+quiet_theme
h <- block(with_role(h),'h','Uncertainty retained across all physical tissues')
page <- wrap_plots(list(a,
  wrap_plots(list(b,c,d),ncol=3,widths=c(1,1,1.6)),e,
  wrap_plots(list(f,g),ncol=2),h,
  note('One point = one physical tissue. Rich-area cutoff: RCTD macrophage ≥0.20. 1 h/4 h: descriptive; N3/N4: auxiliary.')),
  ncol=1,heights=c(105,34,24,29,24,5))
export_remaining(page,'Figure_2')
