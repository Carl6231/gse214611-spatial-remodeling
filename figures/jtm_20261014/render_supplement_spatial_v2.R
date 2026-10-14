#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
mi6<-SAMPLE[3:8];d7<-tail(mi6,3)
MER<-file.path(V1,'source_data/meringue')
object_path<-function(k)file.path(ROOT,'objects/phase_09/niches',paste0(k,'_phase09_niche.rds'))
purple<-c('#F5F3F7','#B7A1C7','#77549A','#39264F')
page_title<-function(id,t)strip_text_v2(paste0('Supplementary ',sub('_',' ',id),' | ',t),8.2,'bold',0)
matrix_theme<-quiet_theme+theme(axis.line=element_blank(),axis.ticks=element_blank())
gene_rows<-function(genes,keys,prefix,all_keys=keys){
  vals<-setNames(lapply(all_keys,function(k){
    o<-readRDS(object_path(k));a<-DefaultAssay(o)
    mat<-suppressWarnings(tryCatch(LayerData(o,assay=a,layer='data'),error=function(e)NULL))
    raw<-is.null(mat)||!all(genes%in%rownames(mat))
    if(raw)mat<-LayerData(o,assay=a,layer='counts')
    stopifnot(all(genes%in%rownames(mat)))
    x<-as.matrix(mat[genes,colnames(o),drop=FALSE]);if(raw)x<-log1p(x)
    rm(o,mat);gc(verbose=FALSE);x
  }),all_keys)
  rows<-lapply(genes,function(g){
    limits<-as.numeric(quantile(unlist(lapply(vals,function(m)m[g,])),c(.01,.99),na.rm=TRUE))
    if(diff(limits)<=0)limits<-range(unlist(lapply(vals,function(m)m[g,])),finite=TRUE)
    stopifnot(diff(limits)>0)
    images<-lapply(keys,function(k){
      id<-paste(prefix,g,k,sep='_');out<-file.path(ASSET_V2,paste0(id,'.png'))
      if(!file.exists(out))out<-spatial_values(object_path(k),vals[[k]][g,],id,limits,purple,width=50,height=40)
      img_panel_compact_v2(out)
    })
    wrap_plots(c(list(strip_text_v2(g,7.2,'italic')),images,
      list(continuous_bar_v2(limits,purple,'Expression',vertical=TRUE))),ncol=5,widths=c(.20,1,1,1,.31))
  })
  names(rows)<-genes;rows
}
if(wants('Figure_S9')){
  spot<-read_table(file.path(ROOT,'metadata/phase_10_spot_signature_scores.tsv.gz'))
  sig<-read_table(file.path(ROOT,'metadata/phase_10_day7_remodeling_signature.tsv'))
  ext<-read_table(file.path(ROOT,'metadata/phase_10_external_signature_directions.tsv'))
  quality<-read_table(file.path(ROOT,'metadata/phase_10_axis_quality_covariates.tsv'))
  lim<-max(abs(quantile(spot[sample_key%in%d7 & is.finite(signature_score),signature_score],c(.01,.99))))
  continuum<-lapply(mi6,function(k){
    z<-spot[sample_key==k];id<-paste0('S9_continuum_',k);path<-file.path(ASSET_V2,paste0(id,'.png'))
    if(!file.exists(path))path<-spatial_values(object_path(k),setNames(z$signature_score,z$barcode),id,c(-lim,lim),c(COL_CONTRACTILE,COL_ZERO,COL_FIBROTIC),midpoint=0,width=54,height=42)
    heading(img_panel_compact_v2(path),LABEL[k])
  })
  a<-block(wrap_plots(list(wrap_plots(continuum,ncol=3),continuous_bar_v2(c(-lim,lim),c(COL_CONTRACTILE,COL_ZERO,COL_FIBROTIC),'Continuum score',midpoint=0,vertical=FALSE)),ncol=1,heights=c(1,.12)),'a','Frozen MI continuum · shared display range')
  o<-readRDS(object_path('MI_Day7_rep2'));mat<-suppressWarnings(tryCatch(LayerData(o,assay=DefaultAssay(o),layer='data'),error=function(e)NULL))
  raw<-is.null(mat)||any(dim(mat)==0L);if(raw)mat<-LayerData(o,assay=DefaultAssay(o),layer='counts')
  mean_set<-function(gs){x<-mat[intersect(gs,rownames(mat)),,drop=FALSE];if(raw)x@x<-log1p(x@x);v<-Matrix::colMeans(x);names(v)<-colnames(o);v}
  scores<-list(ECM=if('P09N_remodeling_axis'%in%colnames(o[[]]))setNames(o$P09N_remodeling_axis,colnames(o)) else mean_set(c('Serpinh1','Thbs2','Cthrc1')),
    'Contractile / OxPhos'=mean_set(sig[direction=='negative',gene]),'Trem2 / Spp1'=mean_set(c('Trem2','Spp1')))
  rm(o,mat);gc(verbose=FALSE)
  score_plots<-lapply(seq_along(scores),function(i){v<-scores[[i]];limit<-as.numeric(quantile(v,c(.01,.99),na.rm=TRUE));id<-paste0('S9_score_',i);p<-file.path(ASSET_V2,paste0(id,'.png'))
    if(!file.exists(p))p<-spatial_values(object_path('MI_Day7_rep2'),v,id,limit,purple,width=54,height=42)
    heading(wrap_plots(list(img_panel_compact_v2(p),continuous_bar_v2(limit,purple,'Score',vertical=TRUE)),ncol=2,widths=c(1,.22)),names(scores)[i])})
  c<-block(wrap_plots(score_plots,ncol=3),'c','Representative scores · MI D7 r2 · separate scales')
  px<-ggplot(ext,aes(internal_median_rho,external_LAD_D7_vs_D3,colour=signature_direction))+geom_hline(yintercept=0,linewidth=.25,colour='#BAC2C8')+geom_vline(xintercept=0,linewidth=.25,colour='#BAC2C8')+
    geom_point(size=.7,alpha=.7)+scale_colour_manual(values=c(negative=COL_CONTRACTILE,positive=COL_FIBROTIC),guide='none')+labs(x='Internal spatial ρ',y='External LAD D7 − D3')+quiet_theme
  pq<-ggplot(quality[condition=='MI' & timepoint=='Day7'],aes(covariate,spearman_rho,group=sample_key,colour=sample_key))+geom_hline(yintercept=0,linewidth=.25,colour='#BAC2C8')+
    geom_line(linewidth=.3)+geom_point(size=1)+scale_colour_manual(values=setNames(c(COL_D3,COL_DESC,COL_SHAM),d7),labels=LABEL[d7],name=NULL)+
    labs(x=NULL,y='Within-slice ρ')+quiet_theme+theme(axis.text.x=element_text(angle=40,hjust=1),legend.position='bottom')+guides(colour=guide_legend(nrow=1))
  summary<-spot[center==TRUE & is.finite(signature_score),.(value=median(signature_score)),by=sample_key]
  ps<-ggplot(summary,aes(factor(sample_key,levels=SAMPLE),value))+geom_hline(yintercept=0,linewidth=.25,colour='#BAC2C8')+geom_point(aes(colour=role(sample_key)),size=1.3)+
    scale_colour_manual(values=role_cols_v2,guide='none')+scale_x_discrete(NULL,labels=short_label)+labs(y='Median center score')+quiet_theme+theme(axis.text.x=element_text(angle=45,hjust=1))
  d<-block(wrap_plots(list(heading(px,'External challenge'),heading(pq,'Technical challenge'),heading(ps,'Physical-tissue summary')),ncol=3),'d','Saved challenges · no new selection or inference')
  pages<-list(wrap_plots(list(page_title('Figure_S9','Continuum and robustness'),a,c,d,note('One spatial summary point = one physical tissue. 1 h/4 h are descriptive; Needle-pass is auxiliary.')),
    ncol=1,heights=c(6,92,58,60,5)))
  genes<-c('Serpinh1','Thbs2','Cthrc1','Myh6','Myl2','Tnni3')
  rows<-gene_rows(genes,d7,'S9_gene')
  for(i in 1:2){idx<-((i-1)*3+1):(i*3)
    pages[[i+1]]<-wrap_plots(list(page_title('Figure_S9',paste0('Fixed-gene atlas · ',i,'/2')),
      block(wrap_plots(list(wrap_plots(c(list(plot_spacer()),lapply(LABEL[d7],strip_text_v2),list(plot_spacer())),ncol=5,widths=c(.20,1,1,1,.31)),
        wrap_plots(rows[idx],ncol=1)),ncol=1,heights=c(.07,1)),'b',if(i==1)'Fixed genes · Day 7 replicates' else 'Fixed genes · continued'),
      note('Each gene has one shared 1–99% display range across the three Day 7 tissues. Frozen expression layer; no normalization rerun.')),
      ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))
  }
  export_supplement(pages,'Figure_S9',225)
}
if(wants('Figure_S10')||wants('Figure_S11')){
  fu<-read_table(file.path(ROOT,'21_scop_expansion/P03_MERINGUE/tables/feature_universe.tsv'))
  tech<-unique(fu[mitochondrial==TRUE|ribosomal==TRUE|hemoglobin==TRUE,feature])
  technical<-function(x)x%in%tech|grepl('^(mt-|Rpl|Rps|Hb[ab]|Gm42418$|Malat1$)',x)
  normsample<-function(x)sub('^MI_D7_','MI_Day7_',sub('^MI_D3_','MI_Day3_',x))
  auto<-rbindlist(lapply(mi6,function(k)read_table(file.path(MER,paste0(sub('Day','D',k),'_autocorrelation.tsv')))))
  auto[,`:=`(sample_key=normsample(biological_sample_id),technical=technical(feature))]
  auto[,stage:=ifelse(grepl('Day3',sample_key),'Day3','Day7')]
}
if(wants('Figure_S10')){
  genes<-c('Sdf4','Abca1','Tmem106a','Nrep','Ighm');pages<-list()
  for(i in 1:2){keys<-mi6[((i-1)*3+1):(i*3)];rows<-gene_rows(genes,keys,'S10_gene',all_keys=mi6)
    pages[[i]]<-wrap_plots(list(page_title('Figure_S10',paste0('MERINGUE autocorrelation · ',if(i==1)'Day 3' else 'Day 7')),
      block(wrap_plots(list(wrap_plots(c(list(plot_spacer()),lapply(LABEL[keys],strip_text_v2),list(plot_spacer())),ncol=5,widths=c(.20,1,1,1,.31)),
        wrap_plots(rows,ncol=1)),ncol=1,heights=c(.045,1)),'a',if(i==1)'Five frozen genes across MI tissues' else 'Five frozen genes · continued'),
      note('Gene-specific 1–99% ranges are shared across all six tissues. MERINGUE provides descriptive statistic/rank; no P values.')),
      ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))
  }
  recur<-auto[sample_key%in%d7 & !technical & rank<=100,.(n=uniqueN(sample_key),median_rank=as.numeric(median(rank)),worst_rank=max(rank),median_statistic=as.numeric(median(statistic))),by=feature][n==3][order(median_rank,worst_rank,-median_statistic,feature)]
  stopifnot(nrow(recur)==16L,identical(head(recur$feature,2),c('Sdf4','Abca1')))
  dat<-auto[feature%in%recur$feature];dat[,`:=`(sample=factor(sample_key,levels=mi6),gene=factor(feature,levels=rev(recur$feature)))]
  matrix<-ggplot(dat,aes(sample,gene))+geom_tile(aes(fill=statistic),colour='white',linewidth=.15)+geom_point(aes(size=rank),shape=21,fill=NA,colour='#273C4D',stroke=.25)+
    scale_fill_gradientn(colours=c('#EFF5FA','#74B5C1','#234F70'),name='Descriptive\nstatistic')+scale_size_area(max_size=3,name='Rank')+scale_x_discrete(NULL,labels=LABEL)+labs(y=NULL)+matrix_theme+
    theme(axis.text.x=element_text(angle=40,hjust=1))
  comp<-auto[,.N,by=.(sample_key,technical)];comp[,fraction:=N/sum(N),by=sample_key]
  pc<-ggplot(comp,aes(factor(sample_key,levels=mi6),fraction,fill=technical))+geom_col(width=.8)+scale_fill_manual(values=c('FALSE'=COL_DESC,'TRUE'='#A9AFB5'),labels=c('Nontechnical','Technical'),name=NULL)+
    scale_x_discrete(NULL,labels=LABEL)+labs(y='Fraction of saved genes')+quiet_theme+theme(axis.text.x=element_text(angle=40,hjust=1),legend.position='bottom')
  distributions<-unlist(lapply(c('rank','statistic'),function(v)lapply(c('Day3','Day7'),function(s){p<-ggplot(auto[stage==s & !technical],aes(.data[[v]]))+geom_histogram(bins=30,fill=if(s=='Day3')COL_D3 else COL_D7)+
    labs(x=if(v=='rank')'Rank' else 'Descriptive statistic',y='Saved gene records')+quiet_theme;heading(p,paste(s,v))})),recursive=FALSE)
  pages[[3]]<-wrap_plots(list(page_title('Figure_S10','Reproducibility and distributions'),
    block(matrix,'b','Day 7 3/3 top-100 genes · all six tissues'),block(pc,'c','Technical-gene composition'),block(wrap_plots(distributions,ncol=2),'d','Descriptive distributions'),
    note('Exact frozen technical exclusion and ranking rules. These spatial statistics are not biological-replicate significance tests.')),
    ncol=1,heights=c(6,104,45,78,5))
  export_supplement(pages,'Figure_S10',245)
}
if(wants('Figure_S11')){
  cross<-rbindlist(lapply(mi6,function(k)read_table(file.path(MER,paste0(sub('Day','D',k),'_cross_correlation_top1000.tsv')))))
  modules<-rbindlist(lapply(mi6,function(k)read_table(file.path(MER,paste0(sub('Day','D',k),'_modules.tsv')))))
  cross[,`:=`(sample_key=normsample(biological_sample_id),technical=technical(feature1)|technical(feature2))]
  modules[,`:=`(sample_key=normsample(biological_sample_id),technical=technical(feature))]
  cross[,stage:=ifelse(grepl('Day3',sample_key),'Day3','Day7')]
  recurrence<-cross[technical==FALSE,.(tissues=uniqueN(sample_key)),by=.(stage,feature1,feature2)][tissues>=2]
  stopifnot(nrow(recurrence)==0L)
  pairs<-lapply(mi6,function(k){d<-head(cross[sample_key==k & !technical][order(rank)],20)
    if(!nrow(d))return(heading(note('No nontechnical pair in the\nfrozen top-1000 table\n\nNEGATIVE_COMPLETE',7.2,.5),LABEL[k]))
    d[,pair:=factor(paste0(rank,'  ',feature1,' – ',feature2),levels=rev(paste0(rank,'  ',feature1,' – ',feature2)))]
    p<-ggplot(d,aes(correlation,pair))+geom_segment(aes(x=0,xend=correlation,yend=pair),colour='#B5C0CA')+geom_point(colour=sample_colour(k),size=1.3)+
      labs(x='Cross-correlation',y=NULL)+quiet_theme
    heading(p,paste0(LABEL[k],' · ',nrow(d),' available pairs'))})
  overlap<-lapply(c('Day3','Day7'),function(st)heading(note('0 exact nontechnical top-1000 pairs\nrecurred in at least 2/3 tissues',7.2,.5),paste(st,'exact overlap')))
  mags<-lapply(c('Day3','Day7'),function(st)heading(ggplot(cross[stage==st & !technical],aes(correlation))+geom_histogram(bins=35,fill=if(st=='Day3')COL_D3 else COL_D7)+
    labs(x='Cross-correlation',y='Saved pairs')+quiet_theme,paste(st,'pair magnitudes')))
  page1<-wrap_plots(list(page_title('Figure_S11','Cross-correlation boundaries'),block(wrap_plots(pairs,ncol=2),'a','Up to 20 nontechnical pairs · labels retain the original saved rank'),
    note('The frozen top-1000 table is the candidate universe; fewer eligible pairs are not replenished. Descriptive correlations; no P values.')),
    ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))
  counts<-modules[technical==FALSE,.N,by=.(sample_key,module,module_size)]
  reps<-counts[order(-N,as.numeric(as.character(module))),.SD[1],by=sample_key]
  spectra<-lapply(mi6,function(k){d<-counts[sample_key==k];heading(ggplot(d,aes(reorder(factor(module),N),N))+geom_col(fill=sample_colour(k))+
    labs(x='Sample-local module',y='Nontechnical genes')+quiet_theme+theme(axis.text.x=element_blank()),LABEL[k])})
  members<-lapply(mi6,function(k){m<-reps[sample_key==k,module][1];d<-head(modules[sample_key==k & module==m & !technical][order(rank)],25);d[,gene:=factor(feature,levels=rev(feature))]
    heading(ggplot(d,aes(rank,gene))+geom_point(colour=sample_colour(k),size=1.2)+labs(x='Saved MERINGUE rank',y=NULL)+quiet_theme,
      paste0(LABEL[k],' · local module ',m))})
  page2<-wrap_plots(list(page_title('Figure_S11','Overlap and sample-local module sizes'),
    block(wrap_plots(overlap,ncol=2),'b','Exact overlap · negative result retained'),
    block(wrap_plots(mags,ncol=2),'c','Cross-correlation magnitude distributions'),
    block(wrap_plots(spectra,ncol=3),'d','Nontechnical module-size spectra'),
    note('Matching, top-k and technical exclusion are unchanged. Module IDs are local to each tissue; no independent experimental validation.')),
    ncol=1,heights=c(6,32,60,125,6))
  page3<-wrap_plots(list(page_title('Figure_S11','Sample-local module membership'),block(wrap_plots(members,ncol=3),'e','Representative modules · first 25 ranked nontechnical members'),
    note('Representative = most nontechnical members; ties use sample-local module ID. IDs are never renamed or aligned across slices.')),
    ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(6,'mm')))
  sets<-list();maps<-list()
  for(k in mi6){m<-reps[sample_key==k,module][1];gs<-modules[sample_key==k & module==m & !technical,feature]
    o<-readRDS(object_path(k));mat<-suppressWarnings(tryCatch(LayerData(o,assay=DefaultAssay(o),layer='data'),error=function(e)NULL));raw<-is.null(mat)||any(dim(mat)==0L)
    if(raw)mat<-LayerData(o,assay=DefaultAssay(o),layer='counts');use<-intersect(gs,rownames(mat));sets[[k]]<-use;x<-mat[use,,drop=FALSE];if(raw)x@x<-log1p(x@x)
    v<-Matrix::colMeans(x);names(v)<-colnames(o);lim<-as.numeric(quantile(v,c(.01,.99),na.rm=TRUE));rm(o,mat,x);gc(verbose=FALSE)
    id<-paste0('S11_module_',k);p<-file.path(ASSET_V2,paste0(id,'.png'));if(!file.exists(p))p<-spatial_values(object_path(k),v,id,lim,purple,width=52,height=40)
    maps[[k]]<-heading(wrap_plots(list(img_panel_compact_v2(p),continuous_bar_v2(lim,purple,'Score',vertical=TRUE)),ncol=2,widths=c(1,.22)),paste0(LABEL[k],' · ',m))
  }
  jac<-lapply(list(head(mi6,3),tail(mi6,3)),function(keys){d<-as.data.table(expand.grid(a=keys,b=keys,stringsAsFactors=FALSE));d[,j:=mapply(function(a,b)length(intersect(sets[[a]],sets[[b]]))/length(union(sets[[a]],sets[[b]])),a,b)]
    heading(ggplot(d,aes(factor(a,levels=keys),factor(b,levels=rev(keys)),fill=j))+geom_tile(colour='white')+geom_text(aes(label=sprintf('%.2f',j)),size=6.5/GG_PT,family=FONT_V2)+
      scale_fill_viridis_c(limits=c(0,1),name='Jaccard')+scale_x_discrete(NULL,labels=LABEL)+scale_y_discrete(NULL,labels=LABEL)+matrix_theme,
      if(grepl('Day3',keys[1]))'Day 3' else 'Day 7')})
  page4<-wrap_plots(list(page_title('Figure_S11','Module spatial context'),
    block(wrap_plots(maps,ncol=3),'f','Representative-module score maps · slice-specific scales'),block(wrap_plots(jac,ncol=2),'g','Gene-set Jaccard · no shared module numbering'),
    note('Scores are deterministic displays of frozen module members; no MERINGUE refitting or clustering.')),
    ncol=1,heights=c(6,155,70,6))
  export_supplement(list(page1,page2,page3,page4),'Figure_S11',245)
}
