#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
src<-file.path(ROOT,'results','figures','Figure_S2','source_data','v4p2')
qc<-read_table(file.path(src,'Figure_S2_spot_QC.tsv.gz'))
summary<-read_table(file.path(src,'Figure_S2_slice_QC_summary.tsv'))
compat<-read_table(file.path(src,'Figure_S2_reference_compatibility.tsv'))
ref<-read_table(file.path(src,'Figure_S2_reference_celltype_counts.tsv'))
marker<-read_table(file.path(src,'Figure_S2_marker_detection.tsv'))
class_marker<-read_table(file.path(src,'Figure_S2_marker_celltype_summary.tsv'))
stopifnot(nrow(qc)==25548L,setequal(unique(qc$sample_key),SAMPLE),
  all(qc$qc_keep),sum(ref$reference_cells)==176526,all(compat$exact_common_genes==22490))
flag_cols<-c('Retained, no review flag'='#C9D2DB','Review flag retained'='#5487A5')
maps<-vector('list',12)
for(i in seq_along(SAMPLE)){
  k<-SAMPLE[i];message('S2 frozen QC: ',k)
  d<-qc[sample_key==k]
  obj<-readRDS(summary[sample_key==k,qc_object])
  d<-d[match(colnames(obj),barcode)]
  stopifnot(identical(d$barcode,colnames(obj)))
  obj$.review<-ifelse(d$qc_review_flag,'Review flag retained','Retained, no review flag')
  path<-spatial_category(obj,'.review',flag_cols,paste0('S2_',k),width=42,height=42)
  maps[[i]]<-block(img_panel_compact_v2(path),letters[i],unname(LABEL[k]))
  rm(obj,d);gc(verbose=FALSE)
}
sumlong<-melt(summary,id.vars='sample_key',measure.vars=c('spots_after','spots_distribution_flagged_retained'),
  variable.name='measure',value.name='n')
account<-ggplot(sumlong,aes(factor(sample_key,levels=SAMPLE),n,colour=measure))+
  geom_point(position=position_dodge(width=.3),size=1.5)+
  scale_colour_manual(values=c(spots_after='#4F78A7',spots_distribution_flagged_retained='#5487A5'),
    labels=c('Retained spots','Review-flagged retained spots'),name=NULL)+
  scale_x_discrete(NULL,labels=short_label)+labs(y='Spots')+quiet_theme+theme(legend.position='bottom')
page1<-wrap_plots(list(strip_text_v2('Supplementary Figure S2 | Frozen QC and spot retention',8.2,'bold',0),
  wrap_plots(maps,ncol=4),legend_strip(flag_cols,rows=1),block(account,'m','Spot accounting'),
  note('25,548 retained spots. Distributional review flags are retained; no new exclusion rule or QC threshold was applied.')),
  ncol=1,heights=c(6,150,7,37,5))

variables<-c('log10_nCount','log10_nFeature','percent.mt','percent.ribo')
metric_labels<-c('log10(UMI + 1)','log10(features + 1)','Mitochondrial reads (%)','Ribosomal reads (%)')
violins<-Map(function(v,lab){
  p<-ggplot(qc,aes(factor(sample_key,levels=SAMPLE),.data[[v]],fill=role(sample_key)))+
    geom_violin(scale='width',linewidth=.2,colour='#77818A')+
    geom_boxplot(width=.12,outlier.shape=NA,linewidth=.22,fill='white')+
    scale_fill_manual(values=role_cols_v2,guide='none')+scale_x_discrete(NULL,labels=short_label)+
    labs(y=lab)+quiet_theme+theme(axis.text.x=element_text(angle=45,hjust=1))
  heading(p,lab)
},variables,metric_labels)
flags<-c('qc_low_count_outlier','qc_low_feature_outlier','qc_high_mt_outlier','qc_high_ribo_outlier')
flag_names<-c('Low UMI','Low features','High mt','High ribo')
fractions<-rbindlist(lapply(seq_along(flags),function(i)
  qc[,.(fraction=mean(get(flags[i]))),by=sample_key][,flag:=flag_names[i]]))
fraction_plot<-ggplot(fractions,aes(factor(sample_key,levels=SAMPLE),factor(flag,levels=rev(flag_names)),fill=fraction))+
  geom_tile(colour='white',linewidth=.2)+scale_fill_gradient(low='#F4F6F8',high='#4F78A7',name='Fraction')+
  scale_x_discrete(NULL,labels=short_label)+labs(y=NULL)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=45,hjust=1))
mat<-as.matrix(qc[,..flags]);storage.mode(mat)<-'numeric';co<-crossprod(mat)
co_dt<-as.data.table(as.table(co));setnames(co_dt,c('first','second','n'))
co_dt[,first:=factor(first,levels=flags,labels=flag_names)]
co_dt[,second:=factor(second,levels=rev(flags),labels=rev(flag_names))]
co_plot<-ggplot(co_dt,aes(first,second,fill=n))+geom_tile(colour='white',linewidth=.2)+
  geom_text(aes(label=n),size=6.5/GG_PT,family=FONT_V2)+
  scale_fill_gradient(low='#F4F6F8',high='#A7C0D9',guide='none')+labs(x=NULL,y=NULL)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=30,hjust=1))
review<-ggplot(summary,aes(factor(sample_key,levels=SAMPLE),review_flag_fraction))+
  geom_point(aes(colour=role(sample_key)),size=1.7)+scale_colour_manual(values=role_cols_v2,guide='none')+
  scale_x_discrete(NULL,labels=short_label)+labs(y='Review-flag fraction')+quiet_theme
page2<-wrap_plots(list(strip_text_v2('Supplementary Figure S2 | Descriptive QC distributions',8.2,'bold',0),
  block(wrap_plots(violins,ncol=2),'n','Retained-spot distributions'),
  wrap_plots(list(block(fraction_plot,'o','Flag composition by tissue'),block(co_plot,'p','Flag co-occurrence · spots')),ncol=2,widths=c(1.35,1)),
  block(review,'q','Review burden by physical tissue'),
  note('Violin and boxplot distributions are descriptive. Co-occurrence is a count, not a test; no spot-level P values.')),
  ncol=1,heights=c(6,112,53,27,5))

six<-c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage','Pericyte','Schwann_cell')
compat_values<-data.frame(source=c('Spatial genes','Reference genes'),fraction=c(
  unique(compat$fraction_spatial_genes_common),unique(compat$fraction_reference_genes_common)))
stopifnot(nrow(compat_values)==2L)
common<-ggplot(compat_values,aes(fraction,source))+geom_segment(aes(x=0,xend=fraction,yend=source),colour='#C9D2DB')+
  geom_point(size=2,colour=COL_SHAM)+geom_text(aes(label=scales::percent(fraction,accuracy=.1)),hjust=-.2,size=6.5/GG_PT,family=FONT_V2)+
  scale_x_continuous('Shared gene fraction',limits=c(0,1.2),breaks=c(0,.5,1))+labs(y=NULL)+quiet_theme
ref_plot<-ggplot(ref,aes(log10(reference_cells),factor(celltype,levels=rev(six))))+
  geom_point(aes(colour=celltype),size=1.8)+geom_text(aes(label=format(reference_cells,big.mark=',',trim=TRUE)),
    hjust=-.2,size=6.5/GG_PT,family=FONT_V2)+scale_colour_manual(values=CELL,guide='none')+
  scale_x_continuous('log10(reference nuclei)',expand=expansion(mult=c(.03,.22)))+
  scale_y_discrete(NULL,labels=CLASS_LABEL)+quiet_theme
class_marker[,sample:=factor(sample_key,levels=SAMPLE)]
class_marker[,cell:=factor(celltype,levels=rev(six))]
gate<-ggplot(class_marker,aes(sample,cell,fill=markers_detected_in_both/n_prespecified_markers))+
  geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=paste0(markers_detected_in_both,'/',n_prespecified_markers)),
    size=6.5/GG_PT,family=FONT_V2)+scale_fill_gradient(low='#F4F6F8',high='#9DB9D3',limits=c(0,1),guide='none')+
  scale_x_discrete(NULL,labels=short_label)+scale_y_discrete(NULL,labels=CLASS_LABEL)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank())
class_heat<-ggplot(class_marker,aes(sample,cell,fill=median_spatial_positive_fraction))+geom_tile(colour='white',linewidth=.2)+
  scale_fill_gradient(low='#F4F6F8',high='#4F78A7',limits=c(0,1),name='Detection\nfraction')+
  scale_x_discrete(NULL,labels=short_label)+scale_y_discrete(NULL,labels=CLASS_LABEL)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank())
marker[,gene:=paste(CLASS_LABEL[celltype],marker,sep=' · ')]
gene_order<-unique(marker[order(match(celltype,six)),gene])
full_heat<-ggplot(marker,aes(factor(sample_key,levels=SAMPLE),factor(gene,levels=rev(gene_order)),fill=spatial_positive_fraction))+
  geom_tile(colour='white',linewidth=.12)+scale_fill_gradient(low='#F4F6F8',high='#4F78A7',limits=c(0,1),name='Detection\nfraction')+
  scale_x_discrete(NULL,labels=short_label)+labs(y=NULL)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=6.5))
page3<-wrap_plots(list(strip_text_v2('Supplementary Figure S2 | Reference compatibility and marker coverage',8.2,'bold',0),
  wrap_plots(list(block(common,'r','22,490 common genes · 12/12 tissues'),block(ref_plot,'s','Six-class reference')),ncol=2),
  block(gate,'t','Prespecified marker gate · detected / required'),
  block(class_heat,'u','Class-level marker detection'),block(full_heat,'v','All frozen markers'),
  note('Reference nuclei are not biological replicates. Schwann cell n=28: limited reference support. Detection is not a pure-cell assignment.')),
  ncol=1,heights=c(6,38,32,32,89,5))
export_supplement(list(page1,page2,page3),'Figure_S2',225)
