#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
P05<-file.path(ROOT,'21_scop_expansion/P05_SECRETED_AXIS_GATE_AMENDED/tables')
X01<-file.path(ROOT,'21_scop_expansion/X01_OMICVERSE_CELLORACLE')
title<-function(id,s)strip_text_v2(paste0(sub('Figure_','Supplementary Figure ',id),' | ',s),8.2,'bold',0)
matrix_theme<-quiet_theme+theme(axis.line=element_blank(),axis.ticks=element_blank())
if(wants('Figure_S12')){
  axes<-read_table(file.path(P05,'frozen_secreted_axis.tsv'))[order(selection_rank)]
  am<-read_table(file.path(P05,'secreted_axis_candidate_matrix.tsv'))
  cv<-read_table(file.path(P05,'expression_coverage.tsv'))
  nh<-read_table(file.path(P05,'neighborhood_celltype_support.tsv'))
  gr<-read_table(file.path(P05,'gradient_relation_support.tsv'))
  stopifnot(nrow(axes)==10L,identical(axes$selection_rank,1:10),!any(axes$external_used_for_selection))
  tissues<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3))
  tissue_labels<-setNames(c(paste0('MI D3 r',1:3),paste0('MI D7 r',1:3)),tissues)
  alabel<-setNames(paste0(axes$selection_rank,'  ',axes$ligand_secact_human,' → ',axes$receptor),axes$axis_id)
  ay<-function(d){d[,axis:=factor(axis_id,levels=rev(axes$axis_id))];d}
  key<-note(paste(vapply(1:5,function(i)paste0(alabel[i],'     |     ',alabel[i+5]),character(1)),collapse='\n'),6.5,0)
  gates<-data.table(gate=c('Candidate axes','CellChat recurrence','Expression coverage','Extra spatial layer','Unique signature mapping','All frozen gates','Selected axes'),
    n=c(uniqueN(am$axis_id),uniqueN(am[spatialcellchat_support==TRUE,axis_id]),uniqueN(am[expression_coverage_pass==TRUE,axis_id]),
        uniqueN(am[additional_spatial_layer_pass==TRUE,axis_id]),uniqueN(am[secact_mapping_unique==TRUE & ligand_in_secact_signature==TRUE,axis_id]),
        uniqueN(am[stage_pass==TRUE,axis_id]),nrow(axes)))
  funnel<-ggplot(gates,aes(n,factor(gate,levels=rev(gate))))+geom_point(colour=COL_SHAM,size=1.8)+
    geom_text(aes(label=n),hjust=-.4,size=6.5/GG_PT,family=FONT_V2)+scale_x_continuous('Axes meeting each saved criterion',expand=expansion(mult=c(.02,.12)))+labs(y=NULL)+quiet_theme
  aux<-am[axis_id%in%axes$axis_id,.(neighborhood=any(neighborhood_support%in%TRUE),
    gradient=if(any(gradient_support%in%c(TRUE,'TRUE')))'TRUE' else if(all(gradient_support=='NOT_EVALUABLE'|is.na(gradient_support)))'NE' else 'FALSE'),by=axis_id]
  d<-ay(merge(axes,aux,by='axis_id',sort=FALSE))
  d[,size_tissue:=1.5+2*n_discovery_supporting_tissues/6]
  matrix<-ggplot(d,aes(y=axis))+geom_hline(yintercept=1:10,colour='#E5E9EC',linewidth=.2)+
    geom_text(aes(x=1,label=paste0(n_supported_stages,'/2')),size=6.5/GG_PT,family=FONT_V2)+
    geom_point(aes(x=2,size=size_tissue),shape=21,fill='#CBDCEC',colour=COL_AXIS)+
    geom_text(aes(x=2,label=paste0(n_discovery_supporting_tissues,'/6')),size=6.5/GG_PT,family=FONT_V2)+
    geom_tile(aes(x=3,fill=minimum_component_coverage),width=.75,height=.8)+
    geom_text(aes(x=3,label=sprintf('%.0f%%',100*minimum_component_coverage)),size=6.5/GG_PT,family=FONT_V2)+
    geom_point(aes(x=4,shape=neighborhood),size=2,colour='#3A765F')+
    geom_point(data=d[gradient!='NE'],aes(x=5,shape=gradient=='TRUE'),size=2,colour=COL_D7)+
    geom_text(data=d[gradient=='NE'],aes(x=5,label='NE'),size=6.5/GG_PT,family=FONT_V2)+
    geom_text(aes(x=6,label=sprintf('%.2g',spatialcellchat_median_score)),size=6.5/GG_PT,family=FONT_V2)+
    scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+scale_size_identity()+
    scale_fill_gradient(low='#EFF5FA',high='#669AC0',limits=c(0,1),guide='none')+
    scale_x_continuous(NULL,breaks=1:6,labels=c('Stages','Tissues','Min.\ncoverage','Neighbor-\nhood','Gradient','CellChat\nscore'),position='top',limits=c(.5,6.5))+
    scale_y_discrete(NULL,labels=alabel)+matrix_theme+theme(axis.text.x=element_text(angle=0,hjust=.5))
  stageplots<-lapply(c('Day3','Day7'),function(st){z<-ay(gr[stage==st & axis_id%in%axes$axis_id]);z[,x:=match(biological_sample_id,tissues[tissues%like%if(st=='Day3')'D3' else 'D7'])]
    p<-ggplot(z,aes(x,axis))+geom_tile(fill='#F4F6F8',colour='white')+
      geom_text(aes(label=ifelse(is.na(ordered_sign_pair),'NE',ordered_sign_pair)),size=6.5/GG_PT,family=FONT_V2)+
      scale_x_continuous(NULL,breaks=1:3,labels=paste0('r',1:3))+scale_y_discrete(NULL,labels=setNames(axes$selection_rank,axes$axis_id))+matrix_theme
    heading(p,st)})
  recplots<-lapply(c('Day3','Day7'),function(st){z<-ay(am[stage==st & axis_id%in%axes$axis_id]);
    heading(ggplot(z,aes(n_same_ordered_axis,axis))+geom_segment(aes(x=0,xend=n_same_ordered_axis,yend=axis),colour='#CAD3DB')+
      geom_point(colour=if(st=='Day3')COL_D3 else COL_D7,size=1.8)+scale_x_continuous('Supporting tissues / 3',limits=c(0,3),breaks=0:3)+
      scale_y_discrete(NULL,labels=setNames(axes$selection_rank,axes$axis_id))+quiet_theme,st)})
  page1<-wrap_plots(list(title('Figure_S12','Frozen candidate secreted axes'),
    block(funnel,'a','Saved eligibility counts · criteria shown separately, not a cumulative funnel'),
    block(matrix,'b','Ten frozen axes · evidence retains its original units'),
    wrap_plots(list(block(wrap_plots(stageplots,ncol=2),'g','Gradient context · NE = not evaluable'),
      block(wrap_plots(recplots,ncol=2),'h','Ordered-axis recurrence')),ncol=2),
    note('All rows follow frozen selection rank. Filled/open circles: supported/not supported. Candidate evidence, not a mechanism or experimental validation.')),
    ncol=1,heights=c(6,48,88,81,7))
  cc<-read_table(file.path(ROOT,'21_scop_expansion/P02/tables/existing_spatialcellchat_evidence.tsv'))
  cc<-cc[condition=='MI' & timepoint%in%c('Day3','Day7') & analysis_role=='primary_replicated' & annotation=='Secreted Signaling' & passes_filter==TRUE]
  cc[,axis_id:=paste(sender,ligand,receptor,receiver,pathway_name,sep='|')]
  cc<-cc[axis_id%in%axes$axis_id & biological_sample_id%in%tissues,.(value=median(score)),by=.(axis_id,biological_sample_id)]
  stopifnot(nrow(cc)>0,!anyDuplicated(cc[,.(axis_id,biological_sample_id)]))
  covs<-cv[axis_id%in%axes$axis_id,.(value=min(coverage_rate)),by=.(axis_id,biological_sample_id,component_role)]
  neigh<-nh[axis_id%in%axes$axis_id,.(value=as.numeric(edge_count)),by=.(axis_id,biological_sample_id)]
  stopifnot(!anyDuplicated(neigh[,.(axis_id,biological_sample_id)]))
  tissueplots<-function(dt,label,limits=NULL){lapply(tissues,function(t){
    z<-merge(data.table(axis_id=axes$axis_id),dt[biological_sample_id==t],by='axis_id',all.x=TRUE);z<-ay(z)
    p<-ggplot(z,aes(value,axis))+geom_segment(aes(x=0,xend=value,yend=axis),colour='#CBD4DD',na.rm=TRUE)+
      geom_point(colour=if(grepl('D3',t))COL_D3 else COL_D7,size=1.4,na.rm=TRUE)+
      geom_text(data=z[is.na(value)],aes(x=0,label='not recorded'),hjust=0,size=6.5/GG_PT,family=FONT_V2)+
      scale_y_discrete(NULL,labels=setNames(axes$selection_rank,axes$axis_id),drop=FALSE)+
      scale_x_continuous(label,limits=limits,breaks=scales::breaks_pretty(n=2))+quiet_theme
    heading(p,tissue_labels[t])})}
  page2<-wrap_plots(list(title('Figure_S12','Per-tissue CellChat and ligand coverage'),key,
    block(wrap_plots(tissueplots(cc,'CellChat score',c(0,max(cc$value))),ncol=3),'c','Six physical tissues · actual saved CellChat scores'),
    block(wrap_plots(tissueplots(covs[component_role=='ligand'],'Ligand coverage',c(0,1)),ncol=3),'d','Ligand nonzero-count coverage'),
    note('Numbers identify the ten frozen axes. Missing CellChat records are not replaced by stage averages or zero. Spots are not biological replicates.')),
    ncol=1,heights=c(6,25,96,96,7))
  page3<-wrap_plots(list(title('Figure_S12','Per-tissue receptor coverage and neighborhood'),key,
    block(wrap_plots(tissueplots(covs[component_role=='receptor_subunit'],'Minimum subunit coverage',c(0,1)),ncol=3),'e','Receptor coverage · minimum across required subunits'),
    block(wrap_plots(tissueplots(neigh,'Directed edge count',c(0,max(neigh$value))),ncol=3),'f','Saved sender → receiver neighborhood support'),
    note('Every required receptor subunit retains the frozen coverage rule. Neighborhoods are slice-local; edge counts are descriptive, not replicates.')),
    ncol=1,heights=c(6,25,96,96,7))
  export_supplement(list(page1,page2,page3),'Figure_S12',245)
}
if(wants('Figure_S13')){
  targets<-c('Bhlhe40','Runx2','Arid5b','Mef2c','Maf','Spi1')
  readx<-function(f)read_table(file.path(X01,'tables',paste0(f,'.tsv')))
  se<-readx('sample_level_effects');me<-readx('phase13_module_effects');tr<-readx('state_transition_summary')
  de<-readx('target_delta_expression');rs<-readx('run_summary')
  coord<-merge(readx('harmony_umap'),readx('state_coordinate_cells'),by='cell_id',sort=FALSE)
  cs<-readx('state_coordinate_by_sample')
  stopifnot(nrow(coord)==4659L,uniqueN(coord$sample_id)==16L,setequal(rs$target,targets),all(rs$status=='NEGATIVE_COMPLETE'),
    all(targets%in%se$target),all(targets%in%de$target),!anyDuplicated(se[,.(target,sample_id,myeloid_state)]))
  statecols<-c(Monocyte_Ccr2='#4477AA',Macrophage_Spp1='#CC6677')
  timecols<-setNames(viridisLite::viridis(6),TIME)
  umap<-function(field,continuous=FALSE){
    p<-ggplot(coord,aes(HarmonyUMAP2D_1,HarmonyUMAP2D_2,colour=.data[[field]]))+geom_point(size=.18,alpha=.7)+coord_equal()+theme_void(base_family=FONT_V2)+
      theme(legend.position='bottom',legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),legend.key.height=unit(2,'mm'))
    if(continuous)p+scale_colour_viridis_c(name='Direction coordinate') else p+scale_colour_manual(values=if(field=='myeloid_state')statecols else timecols,name=NULL,labels=function(x)gsub('_',' ',x))
  }
  cs[,tx:=match(timepoint,TIME)+(frank(sample_id,ties.method='dense')-(uniqueN(sample_id)+1)/2)*.015,by=.(timepoint,myeloid_state)]
  samplecoord<-ggplot(cs,aes(tx,state_coordinate,colour=myeloid_state))+geom_point(size=1.4)+scale_colour_manual(values=statecols,name=NULL)+
    scale_x_continuous(NULL,breaks=1:6,labels=TIME_LABEL)+labs(y='Sample median coordinate')+quiet_theme+theme(legend.position='none')
  vectorplots<-lapply(seq_along(targets),function(i){
    png<-file.path(ASSET_V2,paste0('S13_vector_',targets[i],'.png'))
    if(!file.exists(png)){
      raw<-file.path(ASSET_V2,paste0('S13_native_',i,'.png'))
      rc<-system2('gs',c('-q','-dBATCH','-dNOPAUSE',paste0('-dFirstPage=',i),paste0('-dLastPage=',i),'-sDEVICE=pngalpha','-r300',
        paste0('-sOutputFile=',raw),file.path(X01,'figures/target_vector_fields.pdf')))
      stopifnot(rc==0L,file.exists(raw))
      im<-magick::image_read(raw);info<-magick::image_info(im)
      # The inspected title ends before 9% of page height; preserve all points and arrows below it.
      crop<-magick::image_crop(im,magick::geometry_area(info$width,floor(info$height*.91),0,ceiling(info$height*.09)),repage=TRUE)
      magick::image_write(crop,png)
    }
    heading(img_panel_compact_v2(png),paste0(targets[i],' · NEGATIVE_COMPLETE'))
  })
  page1<-wrap_plots(list(title('Figure_S13','CellOracle · in silico counterfactual prediction'),
    wrap_plots(list(block(umap('myeloid_state'),'a','Frozen input states'),block(umap('timepoint'),'b','Timepoint context')),ncol=2),
    wrap_plots(list(block(umap('state_coordinate',TRUE),'c','Restricted direction coordinate'),block(samplecoord,'d','Biological-sample summaries')),ncol=2),
    block(wrap_plots(vectorplots,ncol=3),'e','All six frozen computational displacement fields'),
    note('The original fields share the frozen embedding, arrow scale = 8 and cell subsampling. 6/6 NEGATIVE_COMPLETE under the preregistered module-consistency rule.')),
    ncol=1,heights=c(6,44,44,122,8))
  se[,tx:=match(timepoint,TIME)+(frank(sample_id,ties.method='dense')-(uniqueN(sample_id)+1)/2)*.025,by=.(target,timepoint,myeloid_state)]
  displacement<-lapply(targets,function(t)heading(ggplot(se[target==t],aes(tx,median_displacement,colour=myeloid_state))+geom_point(size=1.4)+
    scale_colour_manual(values=statecols,guide='none')+scale_x_continuous(NULL,breaks=1:6,labels=TIME_LABEL)+
    scale_y_continuous('Sample median displacement',limits=c(0,max(se$median_displacement)))+quiet_theme,t))
  transition<-lapply(targets,function(t){d<-melt(tr[target==t],id.vars=c('target','source_state'),measure.vars=names(statecols),variable.name='destination',value.name='mass')
    heading(ggplot(d,aes(destination,source_state,fill=mass))+geom_tile(colour='white')+geom_text(aes(label=sprintf('%.3f',mass)),size=6.5/GG_PT,family=FONT_V2)+
      scale_fill_gradient(low='#F0F5FA',high='#6B99BE',limits=c(0,1),guide='none')+
      scale_x_discrete('Destination',labels=c(Monocyte_Ccr2='Monocyte',Macrophage_Spp1='Spp1 macro'))+
      scale_y_discrete('Source',labels=c(Monocyte_Ccr2='Monocyte',Macrophage_Spp1='Spp1 macro'))+matrix_theme,t)})
  page2<-wrap_plots(list(title('Figure_S13','Displacement and descriptive transitions'),
    block(wrap_plots(displacement,ncol=3),'f','One point per biological sample × state'),legend_strip(statecols,gsub('_',' ',names(statecols)),1),
    block(wrap_plots(transition,ncol=3),'g','Frozen model-level transition summaries · descriptive only'),
    note('All targets: NEGATIVE_COMPLETE. A computed transition mass is not a fate prediction, cell conversion or experimental knockout result.')),
    ncol=1,heights=c(6,109,7,92,8))
  ranges<-range(de$delta,finite=TRUE)
  delta<-lapply(targets,function(t){d<-copy(de[target==t][order(delta,gene)]);d[,rank:=seq_len(.N)]
    ends<-d[c(1,.N)];ends[,hjust:=c(0,1)]
    heading(ggplot(d,aes(rank,delta))+geom_hline(yintercept=0,colour='#A7AFB7',linewidth=.3)+geom_line(colour=COL_D7,linewidth=.45)+
      geom_text(data=ends,aes(label=gene,hjust=hjust),vjust=ifelse(ends$delta>0,1.2,-.4),size=6.5/GG_PT,family=FONT_V2)+
      scale_y_continuous('Mean Δ expression',limits=ranges,expand=expansion(mult=.08))+labs(x='All modeled genes · ordered rank')+quiet_theme,t)})
  me[,module_short:=gsub('Macrophage__','',module_id,fixed=TRUE)]
  m<-me[,.(median_shift=median(module_shift_z)),by=.(target,module_short,myeloid_state)]
  lim<-max(abs(m$median_shift))
  programs<-read_table(file.path(ROOT,'results/phase13/recurrent_program_genes.tsv'))[celltype=='Macrophage' & technical_excluded_from_projection==FALSE]
  features<-readx('frozen_feature_universe')
  eligible<-programs[,.(n_genes=sum(gene%in%features$gene)),by=module_id]
  eligible[,module_short:=gsub('Macrophage__','',module_id,fixed=TRUE)]
  stopifnot(nrow(eligible)==8L,identical(eligible[n_genes<5,module_id],'Macrophage__down__peak_1HR'),
    eligible[n_genes<5,n_genes]==3L,setequal(eligible[n_genes>=5,module_short],unique(m$module_short)))
  m<-merge(CJ(target=targets,module_short=eligible$module_short,myeloid_state=names(statecols)),m,
    by=c('target','module_short','myeloid_state'),all.x=TRUE)
  stopifnot(all(m[is.na(median_shift),module_short]%in%eligible[n_genes<5,module_short]))
  m[,`:=`(target=factor(target,levels=targets),module_short=factor(module_short,levels=eligible$module_short))]
  heats<-lapply(names(statecols),function(st)heading(ggplot(m[myeloid_state==st],aes(target,module_short,fill=median_shift))+geom_tile(colour='white',linewidth=.2)+
    geom_text(data=m[myeloid_state==st & is.na(median_shift)],aes(label='N/A'),size=6.5/GG_PT,family=FONT_V2)+
    scale_fill_gradient2(low='#2F6FB0',mid='white',high='#C81D25',midpoint=0,limits=c(-lim,lim),na.value='#E6E8EA',name='Median Δ z',breaks=c(-lim,0,lim),labels=function(x)sprintf('%.3f',x))+
    guides(fill=guide_colorbar(title.position='top',barwidth=unit(40,'mm'),barheight=unit(2,'mm')))+
    labs(x=NULL,y=NULL)+matrix_theme+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom'),gsub('_',' ',st)))
  page3<-wrap_plots(list(title('Figure_S13','Expression shifts and frozen-module boundary'),
    block(wrap_plots(delta,ncol=3),'h','All modeled genes · fixed targets, shared Δ-expression scale'),
    block(wrap_plots(heats,ncol=2),'i','Frozen modules · biological-sample median effects'),
    note('N/A = NOT_APPLICABLE: down__peak_1HR has 3 genes (<5 required). Missing effects are not zero.\n6/6 NEGATIVE_COMPLETE does not mean absence of every modeled effect; no nucleus-level P values.')),
    ncol=1,heights=c(6,117,92,9))
  export_supplement(list(page1,page2,page3),'Figure_S13',245)
}
