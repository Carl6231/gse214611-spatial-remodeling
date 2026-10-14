#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
stopifnot(all(c('CellDimPlot','FeatureDimPlot')%in%getNamespaceExports('scop')))
src<-file.path(ROOT,'results/figures/Figure_3_corrected/corrected_v1r4/source_data')
ann<-read_table(file.path(src,'panel_a_cluster_correction.tsv'))
counts<-read_table(file.path(src,'phase_09_balanced_reference_class_counts.tsv'))
cal<-read_table(file.path(src,'panel_h_mapping_calibration.tsv'))
ms<-read_table(file.path(src,'panel_i_macrophage_time_support.tsv'))
embedding<-read_table(file.path(V1,'source_data/corrected_broad_embedding.tsv'))
balanced<-readRDS(file.path(ROOT,'objects/phase_09/reference_expanded_immune_balanced.rds'))
stopifnot(nrow(embedding)==18849L,ncol(balanced)==18849L,sum(counts$full_cells)==176526,
          sum(counts$balanced_cells)==18849L,counts[celltype=='Schwann_cell',balanced_cells]==28L,
          !anyDuplicated(embedding$cell),setequal(embedding$cell,colnames(balanced)))
embedding<-embedding[match(colnames(balanced),cell)]
em<-as.matrix(embedding[,.(UMAP_1,UMAP_2)]);rownames(em)<-embedding$cell
balanced[['FrozenBroadUMAP']]<-CreateDimReducObject(em,key='FrozenBroadUMAP_',assay='RNA')
classes<-c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage','Monocyte','Neutrophil','Dendritic','B_cell','T_NK','Pericyte','Schwann_cell')
balanced$P09_expanded_celltype<-factor(balanced$P09_expanded_celltype,levels=classes)
meta<-as.data.table(balanced[[]],keep.rownames='cell')
support<-meta[,.(nuclei=.N,samples=uniqueN(sample_id)),by=.(P09_expanded_celltype,timepoint)]
support[,fraction:=nuclei/sum(nuclei),by=P09_expanded_celltype]
stopifnot(sum(support$nuclei)==18849L)
bytime<-ms[,.(nuclei=sum(Freq),samples=uniqueN(sample_id[Freq>0])),by=.(state,timepoint)]
stopifnot(!anyDuplicated(bytime[,.(state,timepoint)]),sum(bytime$nuclei)==sum(ms$Freq))
native_umap<-function(group,cols){
  p<-scop::CellDimPlot(balanced,group.by=group,reduction='FrozenBroadUMAP',label=FALSE,
    show_stat=FALSE,pt.size=.23,pt.alpha=.95,palcolor=unname(cols),legend.position='none',
    theme_use='theme_void',theme_args=list(base_family=FONT_V2,base_size=7.1),title='',subtitle='',xlab='',ylab='',combine=FALSE,raster=FALSE,verbose=FALSE)[[1]]
  while(inherits(p,'patchwork'))p<-p[[1]]
  aes<-unique(c(intersect(names(p$mapping),c('fill','colour')),unlist(lapply(p$layers,function(l)intersect(names(l$mapping),c('fill','colour'))))))
  stopifnot(length(aes)==1L)
  p<-suppressMessages(p+scale_colour_manual(values=cols,aesthetics=aes,drop=FALSE,guide='none'))
  p+coord_equal(expand=FALSE)+theme_void(base_family=FONT_V2)+theme(plot.title=element_blank(),plot.subtitle=element_blank(),legend.position='none')
}
umap<-native_umap('P09_expanded_celltype',CELL[classes])
balanced_plot<-ggplot(counts,aes(balanced_cells,factor(celltype,levels=rev(classes))))+
  geom_segment(aes(x=0,xend=balanced_cells,yend=factor(celltype,levels=rev(classes))),colour='#C2CBD3')+
  geom_point(aes(colour=celltype),size=1.4)+geom_text(aes(label=scales::comma(balanced_cells)),hjust=-.2,size=6.5/GG_PT,family=FONT_V2)+
  scale_colour_manual(values=CELL,guide='none')+scale_y_discrete(NULL,labels=CLASS_LABEL)+
  scale_x_continuous('Balanced nuclei',expand=expansion(mult=c(.02,.25)),breaks=c(0,1500,3000))+quiet_theme+theme(axis.text.y=element_text(size=6.5))
support_plot<-ggplot(support,aes(factor(timepoint,levels=TIME),factor(P09_expanded_celltype,levels=rev(classes)),size=samples,fill=fraction))+
  geom_point(shape=21,colour='#466377',stroke=.2)+scale_size_area(max_size=3,name='Samples')+
  scale_fill_gradient(low='#F0F4F8',high=COL_SHAM,limits=c(0,1),name='Within-class\nfraction')+
  scale_x_discrete(NULL,labels=TIME_LABEL)+scale_y_discrete(NULL,labels=CLASS_LABEL)+quiet_theme+theme(axis.text.y=element_text(size=6.5))+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right')
cal_long<-melt(cal,id.vars=c('sample_key','analysis_group'),measure.vars=c('median_cosine_similarity','dominant_agreement_rate','uncertain_fraction'),variable.name='metric')
metric_labels<-c(median_cosine_similarity='Cosine similarity',dominant_agreement_rate='Dominant agreement',uncertain_fraction='Uncertain fraction')
cal_plot<-ggplot(cal_long,aes(factor(sample_key,levels=SAMPLE),value))+geom_point(aes(colour=role(sample_key)),size=1.2)+
  facet_wrap(~metric,nrow=1,labeller=as_labeller(metric_labels),scales='free_y')+scale_colour_manual(values=role_cols_v2,guide='none')+
  scale_x_discrete(NULL,labels=short_label)+scale_y_continuous(NULL,breaks=scales::breaks_pretty(n=2))+quiet_theme+
  theme(axis.text.x=element_text(angle=55,hjust=1))
state_support<-ggplot(bytime,aes(factor(timepoint,levels=TIME),factor(state,levels=rev(names(STATE)[1:5])),fill=log10(nuclei+1)))+
  geom_tile(colour='white',linewidth=.15)+geom_text(aes(label=samples),size=6.5/GG_PT,family=FONT_V2)+
  scale_fill_gradient(low='#F0F4F8',high='#7297BD',guide='none')+scale_x_discrete(NULL,labels=TIME_LABEL)+scale_y_discrete(NULL,labels=STATE_LABEL)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank())
markers<-c('Tnnt2','Pecam1','Col1a1','C1qa','Ccr2','S100a8','Flt3','Cd79a','Cd3d','Rgs5','Sox10')
if(wants('Figure_3')){
  capture<-ragg::agg_capture(width=183,height=225,units='mm',res=300)
  compact_key<-function(cols,labels){g<-ggplotGrob(legend_strip(cols,labels,rows=2));
    wrap_elements(full=g$grobs[[which(g$layout$name=='guide-box-bottom')]])}
  old<-readRDS(file.path(ROOT,'results/figures/Figure_3_corrected/corrected_v1/Figure_3_corrected_v1_plot.rds'))
  colors<-attr(old$patches$plots[[2]]$patches$plots[[1]],'grobs')$full$children[['GRID.recordedGrob.595']]$list$x$gp$fill
  colors<-substr(colors[markers,classes,drop=FALSE],1,7)
  image<-png::readPNG(file.path(ROOT,'results/figures/Figure_3_corrected/corrected_v1r4/panels/panel_c_scop_groupheatmap.png'))
  stopifnot(identical(dim(image),c(1196L,2771L,3L)))
  ys<-round(319+(seq_len(11)-.5)*(1159-319)/11);xs<-round(268+(seq_len(11)-.5)*(2069-268)/11)
  displayed<-matrix(vapply(seq_len(121),function(i){row<-(i-1)%%11+1;col<-(i-1)%/%11+1;
    rgb(image[ys[row],xs[col],1],image[ys[row],xs[col],2],image[ys[row],xs[col],3])},character(1)),11,11)
  stopifnot(identical(toupper(unname(colors)),toupper(displayed)))
  heat_data<-as.data.table(as.table(colors));setnames(heat_data,c('gene','celltype','display_color'))
  heat_data[,`:=`(gene=factor(gene,levels=rev(markers)),celltype=factor(celltype,levels=classes))]
  heat<-ggplot(heat_data,aes(celltype,gene,fill=display_color))+geom_tile()+scale_fill_identity()+
    scale_x_discrete(NULL,labels=CLASS_LABEL)+labs(y=NULL)+quiet_theme+theme(axis.text.y=element_text(size=6.5))+
    theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=50,hjust=1,margin=margin(0,0,0,0)))
  heat<-wrap_plots(list(heat,continuous_bar_v2(c(-1,1),c('#2F6FB0','#F7F7F7','#C81D25'),'Frozen row z',midpoint=0,vertical=TRUE)),ncol=2,widths=c(1,.16))
  rm(old,image);gc(verbose=FALSE)
  audit_text<-paste0(length(unique(ann$P6_cluster)),' original myeloid clusters  →  canonical-lineage audit  →  corrected broad classes\n',
    scales::comma(sum(ann$cluster_cells)),' audited nuclei; 176,526 full-reference nuclei; 18,849 in the balanced reference')
  a<-block(note(audit_text,6.5,.5),'a','Frozen reference correction')
  b<-block(wrap_plots(list(heading(umap,'Corrected broad-class UMAP'),heading(heat,'Canonical lineage markers')),ncol=2,widths=c(.72,1.28)),'b','Reference identity and lineage markers')
  keys<-c('Sham_Day3_rep1','MI_Day3_rep2','MI_Day7_rep2')
  maps<-lapply(keys,function(k)heading(img_panel_compact_v2(spatial_category(file.path(ROOT,'objects/phase_09/deconvolution',paste0(k,'_phase09_expanded.rds')),
    'P09E_dominant_type',CELL,paste0('F3_broad_',k),width=54,height=36)),LABEL[k]))
  e<-block(wrap_plots(list(wrap_plots(maps,ncol=3),compact_key(CELL,unname(CLASS_LABEL[names(CELL)]))),ncol=1,
    heights=unit.c(unit(1,'null'),unit(9,'mm'))),'e','Broad-class spatial maps')
  fine_keys<-keys[2:3]
  fine_maps<-lapply(fine_keys,function(k){o<-readRDS(file.path(ROOT,'objects/phase_09/macrophage_mapping',paste0(k,'_phase09_macrophage.rds')));cells<-colnames(o)[o$P09M_candidate%in%TRUE]
    p<-spatial_category(o,'P09M_dominant_state',STATE,paste0('F3_fine_',k),cells=cells,width=54,height=35);rm(o);heading(img_panel_compact_v2(p),LABEL[k])})
  h<-block(wrap_plots(list(wrap_plots(fine_maps,ncol=2),compact_key(STATE,unname(STATE_LABEL))),ncol=1,
    heights=unit.c(unit(1,'null'),unit(9,'mm'))),'h','Fine-state spatial maps')
  support_plot<-support_plot+scale_size_area(max_size=3,name='Samples',breaks=c(1,3,6))+
    scale_fill_gradient(low='#F0F4F8',high=COL_SHAM,limits=c(0,1),breaks=c(0,.5,1),name='Fraction')+
    guides(size=guide_legend(nrow=1),fill=guide_colorbar(barwidth=unit(12,'mm'),barheight=unit(2,'mm')))+
    theme(legend.position='bottom',legend.box='horizontal',legend.title=element_text(size=7.1,margin=margin(r=11)),legend.spacing.x=unit(2,'mm'))
  rows<-list(a,b,wrap_plots(list(block(balanced_plot,'c','Balanced reference counts'),block(support_plot,'d','Reference timepoint support')),ncol=2),
    e,block(cal_plot,'f','Mapping calibration'),
    wrap_plots(list(block(state_support,'g','Five-state support · n = samples'),h),ncol=2,widths=c(.43,.57)),
    note('Nuclei and spots are not biological replicates. 1 h/4 h are descriptive; Needle-pass is auxiliary.'))
  heights<-c(10,50,45,44,29,38,5)
  stopifnot(length(rows)==length(heights),sum(heights)==221)
  grobs<-lapply(seq_along(rows),function(i){p<-wrap_plots(list(rows[[i]])) & theme(plot.margin=margin(0,0,0,0));
    g<-patchworkGrob(p)
    grobTree(g,vp=viewport(layout.pos.row=i,layout.pos.col=1))})
  page<-wrap_elements(full=grobTree(children=do.call(gList,grobs),
    vp=viewport(layout=grid.layout(length(rows),1,heights=unit(heights,'null')))))+theme(plot.margin=margin(0,0,0,0))
  dev.off()
  export_remaining(page,'Figure_3')
}
if(wants('Figure_S6')){
  cluster<-ggplot(ann,aes(cluster_cells,reorder(paste0('C',P6_cluster),cluster_cells)))+geom_point(aes(colour=expanded_celltype),size=1.5)+scale_colour_manual(values=CELL,guide='none')+scale_x_continuous('Nuclei',breaks=c(0,2500,5000),labels=c('0','2.5k','5k'))+labs(y=NULL)+quiet_theme
  calls<-ann[,.(nuclei=sum(cluster_cells)),by=audited_state]
  audit<-ggplot(calls,aes(nuclei,reorder(gsub('_',' ',audited_state),nuclei)))+geom_point(size=1.5,colour=COL_DESC)+labs(x='Nuclei',y=NULL)+quiet_theme
  corrected<-ann[,.(nuclei=sum(cluster_cells)),by=expanded_celltype]
  corr<-ggplot(corrected,aes(nuclei,reorder(expanded_celltype,nuclei),colour=expanded_celltype))+geom_point(size=1.5)+scale_colour_manual(values=CELL,guide='none')+scale_y_discrete(NULL,labels=CLASS_LABEL)+labs(x='Audited-compartment nuclei')+quiet_theme
  balance<-ggplot(counts,aes(full_cells,balanced_cells,colour=celltype))+geom_point(size=2)+scale_x_log10('Full reference nuclei',labels=scales::label_number(scale_cut=scales::cut_short_scale()))+labs(y='Balanced reference nuclei')+scale_colour_manual(values=CELL,guide='none')+quiet_theme
  sample_count<-meta[,.N,by=.(sample_id,timepoint)]
  ps<-ggplot(sample_count,aes(reorder(sample_id,N),N,fill=factor(timepoint,levels=TIME)))+geom_col()+scale_fill_manual(values=c(Day0=COL_SHAM,'1HR'='#78C7DB','4HR'=COL_DESC,Day1='#63BCA3',Day3=COL_D3,Day7=COL_D7),name=NULL)+labs(x=NULL,y='Balanced nuclei')+quiet_theme+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position='bottom')
  pu<-ggplot(cal,aes(factor(sample_key,levels=SAMPLE),uncertain_fraction))+geom_point(aes(colour=role(sample_key)),size=1.4)+scale_colour_manual(values=role_cols_v2,guide='none')+
    scale_x_discrete(NULL,labels=short_label)+labs(y='Uncertain fraction')+quiet_theme+theme(axis.text.x=element_text(angle=45,hjust=1))
  page1<-wrap_plots(list(strip_text_v2('Supplementary Figure S6 | Lineage correction and reference balance',8.2,'bold',0),
    wrap_plots(list(block(cluster,'a','Original clusters'),block(audit,'b','Frozen audit calls'),block(corr,'c','Corrected broad classes')),ncol=3,widths=c(.75,1.5,1)),
    wrap_plots(list(block(balanced_plot,'d','Balanced counts'),block(balance,'e','Full-to-balanced representation')),ncol=2),
    block(ps,'f','Contributing biological samples'),block(cal_plot,'g','Cross-method calibration'),block(pu,'h','Uncertain burden retained'),
    note('Descriptive reference support and mapping sensitivity, not independent experimental validation. Schwann n = 28 remains limited.')),
    ncol=1,heights=c(6,66,46,34,34,30,5))
  raw<-LayerData(balanced,assay='RNA',layer='counts');stopifnot(all(markers%in%rownames(raw)))
  marker_plots<-lapply(seq_along(markers),function(i){g<-markers[i];field<-paste0('raw_view_',g);v<-log1p(as.numeric(raw[g,colnames(balanced)]));balanced[[field]]<-v
    p<-scop::FeatureDimPlot(balanced,features=field,reduction='FrozenBroadUMAP',assay='RNA',layer='counts',
      lower_cutoff=0,upper_cutoff=max(v),lower_quantile=0,upper_quantile=1,bg_cutoff=0,
      palcolor=c('#DCE3EA','#A78CC9','#61468C'),pt.size=.22,pt.alpha=1,combine=FALSE,
      raster=FALSE,show_stat=FALSE,title='',subtitle='',xlab='',ylab='',legend.title='',legend.position='right',
      theme_use='theme_void',theme_args=list(base_family=FONT_V2,base_size=7.1),verbose=FALSE)[[1]]
    while(inherits(p,'patchwork'))p<-p[[1]]
    p<-p+labs(x=NULL,y=NULL)+theme_void(base_family=FONT_V2)+theme(plot.title=element_blank(),plot.subtitle=element_blank(),strip.text=element_blank(),strip.background=element_blank(),legend.text=element_text(size=7.1),legend.key.height=unit(2.4,'mm'))
    heading(p,g)})
  page2<-wrap_plots(list(strip_text_v2('Supplementary Figure S6 | Canonical markers on the frozen UMAP',8.2,'bold',0),
    block(wrap_plots(marker_plots,ncol=3),'i','11 fixed markers · log1p(raw counts), gene-specific ranges'),
    note('Frozen RNA counts displayed as log1p(raw counts). No normalization, feature selection, PCA, clustering or UMAP was rerun.')),
    ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))
  class_plots<-lapply(classes,function(cl){d<-support[P09_expanded_celltype==cl]
    p<-ggplot(d,aes(factor(timepoint,levels=TIME),fraction))+geom_col(fill=CELL[cl],width=.7)+geom_text(aes(label=paste0('n=',samples)),vjust=-.4,size=6.5/GG_PT,family=FONT_V2)+
      scale_x_discrete(NULL,labels=TIME_LABEL,drop=FALSE)+scale_y_continuous('Within-class fraction',limits=c(0,1.08),breaks=c(0,.5,1))+quiet_theme
    heading(p,CLASS_LABEL[cl])})
  page3<-wrap_plots(list(strip_text_v2('Supplementary Figure S6 | Class-by-timepoint support',8.2,'bold',0),
    block(wrap_plots(class_plots,ncol=3),'j','11 broad classes · n denotes contributing biological samples'),
    note('Bars describe nucleus composition; sample counts show support. Nuclei are not independent biological replicates.')),
    ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))
  export_supplement(list(page1,page2,page3),'Figure_S6',245)
}
