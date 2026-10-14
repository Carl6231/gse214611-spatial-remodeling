#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
stopifnot('SpatialSpotPlot' %in% getNamespaceExports('scop'))
if(wants('Figure_S1')){
  blocks<-vector('list',12)
  for(i in seq_along(SAMPLE)){
    k<-SAMPLE[i];message('S1: ',k)
    obj<-readRDS(file.path(ROOT,'objects','phase_09','deconvolution',paste0(k,'_phase09_expanded.rds')))
    obj$.retained<-factor(rep('Retained spots',ncol(obj)))
    he<-scop::SpatialSpotPlot(obj,group.by='.retained',image=Images(obj)[[1]],
      plot_type='point',crop=TRUE,overlay_image=TRUE,image.alpha=1,pt.alpha=0,pt.size=.46,
      stroke=0,palcolor='#4F78A7',legend.position='none',theme_use='theme_blank')
    he<-spatial_leaf(he)
    he_path<-file.path(ASSET_V2,paste0('S1_HE_',k,'.png'))
    if(!file.exists(he_path))he_path<-save_spatial_asset_v2(he,paste0('S1_HE_',k),42,40)
    overlay<-spatial_category(obj,'.retained',c('Retained spots'='#4F78A7'),paste0('S1_OV_',k),width=42,height=40)
    blocks[[i]]<-block(wrap_plots(list(heading(img_panel_compact_v2(he_path),'H&E'),
      heading(img_panel_compact_v2(overlay),'Retained spots')),ncol=2),letters[i],unname(LABEL[k]))
    rm(obj,he);gc(verbose=FALSE)
  }
  pages<-atlas_pages(blocks,'Supplementary Figure S1 | Image and spot alignment',
    'The same tissue coordinates are used in each pair. 1 h/4 h are descriptive; Needle-pass is auxiliary.')
  export_supplement(pages,'Figure_S1',190)
  rm(blocks,pages);gc(verbose=FALSE)
}
if(wants('Figure_S3')){
  dt<-read_table(file.path(ROOT,'results','figures','Figure_S3','source_data','Figure_S3_spot_deconvolution.tsv.gz'))
  stopifnot(nrow(dt)==25548L,!anyDuplicated(dt[,.(sample_key,barcode)]))
  cols<-CELL[c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage','Pericyte','Schwann_cell')]
  cc<-c(High='#1B7837',Moderate='#8B78B6',Uncertain=COL_UNCERTAIN)
  blocks<-vector('list',12)
  for(i in seq_along(SAMPLE)){
    k<-SAMPLE[i];message('S3: ',k)
    obj<-readRDS(file.path(ROOT,'objects','deconvolution','samples',paste0(k,'_phase04.rds')))
    d<-dt[sample_key==k][match(colnames(obj),barcode)]
    stopifnot(identical(d$barcode,colnames(obj)))
    fields<-c('rctd_dominant','dwls_dominant','display_confidence')
    titles<-c('RCTD dominant','SpatialDWLS dominant','Cross-method confidence')
    paths<-character(3)
    for(j in 1:3){
      obj$.display<-d[[fields[j]]]
      paths[j]<-spatial_category(obj,'.display',if(j<3)cols else cc,paste0('S3_',j,'_',k),width=56,height=42)
    }
    blocks[[i]]<-block(wrap_plots(Map(function(p,t)heading(img_panel_compact_v2(p),t),paths,titles),ncol=3),letters[i],unname(LABEL[k]))
    rm(obj,d);gc(verbose=FALSE)
  }
  keys<-legend_strip(c(cols,cc),c(unname(CLASS_LABEL[names(cols)]),'High','Moderate','Uncertain'),rows=2)
  pages<-atlas_pages(blocks,'Supplementary Figure S3 | Deconvolution sensitivity',
    'Dominant components are not pure-cell calls. RCTD is primary; SpatialDWLS is sensitivity evidence. Saved Low is shown as Uncertain.',
    ncol=1,per_page=4,legend=keys)
  export_supplement(pages,'Figure_S3',245)
  rm(blocks,pages,dt);gc(verbose=FALSE)
}
for(figure in c('Figure_S7','Figure_S8'))if(wants(figure)){
  fine<-figure=='Figure_S7'
  blocks<-vector('list',12);confidence<-vector('list',12);sample_summary<-vector('list',12);assignment_counts<-vector('list',12)
  for(i in seq_along(SAMPLE)){
    k<-SAMPLE[i];message(figure,': ',k)
    path<-if(fine)file.path(ROOT,'objects','phase_09','macrophage_mapping',paste0(k,'_phase09_macrophage.rds')) else
      file.path(ROOT,'objects','niches',paste0(k,'_phase07.rds'))
    obj<-readRDS(path);md<-obj[[]]
    selected<-if(fine)rownames(md)[md$P09M_candidate %in% TRUE] else rownames(md)
    group_field<-if(fine)'P09M_dominant_state' else 'P7_niche'
    conf_field<-if(fine)'P09M_confidence' else 'P7_confidence'
    cols<-if(fine)STATE else NICHE
    cc<-c(High='#1B7837',Moderate='#8B78B6',Low='#C81D25',Uncertain=COL_UNCERTAIN,Not_applicable='#D8D8D8')
    p1<-spatial_category(obj,group_field,cols,paste0(figure,'_group_',k),cells=selected,width=42,height=40)
    p2<-spatial_category(obj,conf_field,cc,paste0(figure,'_confidence_',k),cells=selected,width=42,height=40)
    blocks[[i]]<-block(wrap_plots(list(heading(img_panel_compact_v2(p1),if(fine)'Fine state' else 'Niche assignment'),
      heading(img_panel_compact_v2(p2),'Mapping confidence')),ncol=2),letters[i],unname(LABEL[k]))
    centers<-if(fine)selected else rownames(md)[!is.na(md$P7_niche)&md$P7_niche!='Not_macrophage_rich']
    vals<-as.character(md[centers,group_field]);cf<-as.character(md[centers,conf_field])
    confidence[[i]]<-data.table(sample_key=k,level=cf)[,.N,by=.(sample_key,level)]
    sample_summary[[i]]<-data.table(sample_key=k,n_centers=length(centers),uncertain=mean(vals=='Uncertain'))
    assignment_counts[[i]]<-data.table(sample_key=k,niche=names(NICHE)[1:5],fraction=vapply(names(NICHE)[1:5],function(n)mean(vals==n),numeric(1)))
    rm(obj,md);gc(verbose=FALSE)
  }
  cf<-rbindlist(confidence);cf[,`:=`(fraction=N/sum(N),level=factor(level,levels=names(cc))),by=sample_key]
  ss<-rbindlist(sample_summary);ss[,sample:=factor(sample_key,levels=SAMPLE)]
  conf_plot<-ggplot(cf,aes(factor(sample_key,levels=SAMPLE),fraction,fill=level))+geom_col(width=.8)+
    scale_fill_manual(values=cc,drop=FALSE,name=NULL)+scale_x_discrete(NULL,labels=short_label)+
    labs(y=if(fine)'Candidate-spot fraction' else 'Center fraction')+quiet_theme+theme(legend.position='bottom',
      axis.text.x=element_text(angle=45,hjust=1))
  uncertain_plot<-ggplot(ss,aes(sample,uncertain))+geom_point(aes(colour=role(sample_key)),size=1.5)+
    scale_colour_manual(values=role_cols_v2,guide='none')+scale_x_discrete(NULL,labels=short_label)+
    scale_y_continuous('Uncertain fraction',limits=c(0,1))+quiet_theme+
    theme(axis.text.x=element_text(angle=45,hjust=1))
  keys<-wrap_plots(list(legend_strip(cols,if(fine)unname(STATE_LABEL[names(cols)]) else names(cols),rows=2),
                        legend_strip(cc,rows=1)),ncol=1)
  pages<-atlas_pages(blocks,paste0('Supplementary ',sub('Figure_','Figure ',figure),' | ',if(fine)'Fine-state mapping' else 'Niche atlas'),
    if(fine)'Fine-state labels are mixed-spot and method-sensitive. Gray background is not a candidate cell-state assignment.' else
      'Networks and niche assignments are slice-local. Uncertain and low confidence are retained.',legend=keys)
  if(fine){
    src<-file.path(ROOT,'results','figures','Figure_3_corrected','corrected_v1r4','source_data')
    ms<-read_table(file.path(src,'panel_i_macrophage_time_support.tsv'))
    support<-read_table(file.path(src,'phase_09_macrophage_state_support.tsv'))
    ms[,state:=factor(state,levels=rev(names(STATE)[1:5]))]
    bytime<-ms[,.(nuclei=sum(Freq),samples=uniqueN(sample_id[Freq>0])),by=.(state,timepoint)]
    stopifnot(!anyDuplicated(bytime[,.(state,timepoint)]),sum(bytime$nuclei)==sum(ms$Freq))
    p1<-ggplot(ms,aes(sample_id,state,fill=log10(Freq+1)))+geom_tile(colour='white',linewidth=.15)+
      scale_fill_gradient(low='#F3F6FA',high='#4F78A7',name='log10(n+1)')+
      scale_y_discrete(NULL,labels=STATE_LABEL)+labs(x=NULL)+quiet_theme+
      theme(axis.text.x=element_text(angle=90,hjust=1),axis.line=element_blank(),axis.ticks=element_blank())
    p2<-ggplot(bytime,aes(factor(timepoint,levels=TIME),state,fill=log10(nuclei+1)))+
      geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=paste0('n=',samples)),size=6.5/GG_PT,family=FONT_V2)+
      scale_fill_gradient(low='#F3F6FA',high='#80A8CE',name='log10 nuclei')+
      scale_x_discrete(NULL,labels=TIME_LABEL)+scale_y_discrete(NULL,labels=STATE_LABEL)+quiet_theme+
      theme(axis.line=element_blank(),axis.ticks=element_blank())
    p3<-ggplot(support,aes(cells,factor(state,levels=rev(names(STATE)[1:5]))))+geom_point(size=2,colour=COL_D7)+
      geom_text(aes(label=cells),hjust=-.2,size=6.5/GG_PT,family=FONT_V2)+
      scale_x_continuous('Reference nuclei',expand=expansion(mult=c(.03,.2)))+
      scale_y_discrete(NULL,labels=STATE_LABEL)+quiet_theme
    broad<-read_table(file.path(src,'panel_h_mapping_calibration.tsv'))
    p4<-ggplot(broad,aes(factor(sample_key,levels=SAMPLE),uncertain_fraction))+geom_point(size=1.6,colour=COL_UNCERTAIN)+
      scale_x_discrete(NULL,labels=short_label)+labs(y='Uncertain fraction')+quiet_theme+
      theme(axis.text.x=element_text(angle=45,hjust=1))
    summary_page<-wrap_plots(list(strip_text_v2('Supplementary Figure S7 | Support and uncertainty',8.2,'bold',0),
      block(p1,'m','Reference sample support'),block(p2,'n','Reference timepoint support · n = samples'),
      wrap_plots(list(block(p3,'o','Five-state reference'),block(p4,'p','Broad-map uncertainty')),ncol=2),
      wrap_plots(list(block(conf_plot,'q','Fine-map confidence'),block(uncertain_plot,'r','Fine-map uncertainty')),ncol=2),
      note('Support counts are descriptive, not independent nuclei-level evidence. Each spatial summary point represents a physical tissue.')),
      ncol=1,heights=c(6,45,35,34,38,5))
  }else{
    np<-read_table(file.path(ROOT,'metadata','phase_07_niche_summary.tsv'))
    actual<-rbindlist(assignment_counts)
    check<-merge(np,actual,by=c('sample_key','niche'))
    stopifnot(all(abs(check$fraction_of_macrophage_rich_centers-check$fraction)<1e-10))
    stageplots<-lapply(names(NICHE)[1:4],function(n){
      d<-actual[niche==n];d[,sample:=factor(sample_key,levels=SAMPLE)]
      p<-ggplot(d,aes(sample,fraction))+geom_point(aes(colour=role(sample_key)),size=1.5)+
        scale_colour_manual(values=role_cols_v2,guide='none')+scale_x_discrete(NULL,labels=short_label,limits=SAMPLE,drop=FALSE)+
        scale_y_continuous('Center fraction',limits=c(0,NA))+quiet_theme+
        theme(axis.text.x=element_text(angle=45,hjust=1))
      heading(p,paste(strwrap(sub(' niche$','',n),width=24),collapse='\n'))
    })
    co<-read_table(file.path(ROOT,'metadata','phase_07_spatial_coherence.tsv'))
    pc<-ggplot(co,aes(factor(sample_key,levels=SAMPLE),enrichment))+geom_hline(yintercept=0,colour='#A8AFB6',linewidth=.3)+
      geom_point(aes(colour=role(sample_key)),size=1.5)+scale_colour_manual(values=role_cols_v2,guide='none')+
      scale_x_discrete(NULL,labels=short_label)+labs(y='Observed − permuted')+quiet_theme+
      theme(axis.text.x=element_text(angle=45,hjust=1))
    pn<-ggplot(ss,aes(sample,n_centers))+geom_point(aes(colour=role(sample_key)),size=1.5)+
      scale_colour_manual(values=role_cols_v2,guide='none')+scale_x_discrete(NULL,labels=short_label)+labs(y='Macrophage-rich centers')+
      quiet_theme+theme(axis.text.x=element_text(angle=45,hjust=1))
    summary_page<-wrap_plots(list(strip_text_v2('Supplementary Figure S8 | Cross-tissue niche heterogeneity',8.2,'bold',0),
      block(wrap_plots(stageplots,ncol=2),'m','Niche prevalence · one point per tissue'),
      wrap_plots(list(block(pc,'n','Within-slice coherence'),block(uncertain_plot,'o','Uncertain assignments')),ncol=2),
      wrap_plots(list(block(conf_plot,'p','Confidence retained'),block(pn,'q','Center accounting')),ncol=2),
      note('1 h/4 h are descriptive; Needle-pass is auxiliary. No spot-level P values or cross-slice neighborhood calculations.')),
      ncol=1,heights=c(6,78,42,42,5))
  }
  pages<-c(pages,list(summary_page))
  if(fine){
    myeloid<-readRDS(file.path(ROOT,'03_snrna_reference/myeloid_refinement/myeloid_reclustered_annotated.rds'))
    genes<-c('Ccr2','H2-Ab1','Trem2','Spp1')
    expr<-LayerData(myeloid,assay='decontXcounts',layer='data')
    stopifnot(all(genes%in%rownames(expr)),'HarmonyUMAP2D'%in%Reductions(myeloid))
    marker_plots<-lapply(genes,function(g){
      lim<-as.numeric(quantile(as.numeric(expr[g,]),c(.01,.99)))
      p<-scop::FeatureDimPlot(myeloid,features=g,reduction='HarmonyUMAP2D',assay='decontXcounts',layer='data',
        lower_cutoff=lim[1],upper_cutoff=lim[2],lower_quantile=0,upper_quantile=1,
        palcolor=colorRampPalette(c('#EFF5FA','#74B5C1','#234F70'))(9),pt.size=.13,pt.alpha=.85,combine=FALSE,
        raster=FALSE,show_stat=FALSE,title='',subtitle='',xlab='',ylab='',legend.title='',legend.position='right',
        theme_use='theme_void',theme_args=list(base_family=FONT_V2,base_size=6.5),verbose=FALSE)[[1]]
      while(inherits(p,'patchwork'))p<-p[[1]]
      p<-p+labs(x=NULL,y=NULL)+theme_void(base_family=FONT_V2)+theme(plot.title=element_blank(),plot.subtitle=element_blank(),strip.text=element_blank(),strip.background=element_blank(),legend.text=element_text(size=6.5))
      heading(p,g)
    })
    pages<-c(pages,list(wrap_plots(list(strip_text_v2('Supplementary Figure S7 | Frozen myeloid markers',8.2,'bold',0),
      block(wrap_plots(marker_plots,ncol=2),'s','Four fixed markers · existing log1p(CP10K) expression'),
      note('Frozen HarmonyUMAP2D; gene-specific pooled 1–99% display ranges. No normalization or embedding was recalculated.')),
      ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))))
    rm(myeloid,expr,marker_plots);gc(verbose=FALSE)
  }
  export_supplement(pages,figure,210)
  rm(blocks,pages,confidence,sample_summary);gc(verbose=FALSE)
}
