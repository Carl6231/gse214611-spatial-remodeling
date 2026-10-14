source('/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_final/script/figure_remaining_helpers.R')
A<-file.path(ROOT,'22_JTM_completion_20261014/P04_STATIAL_corrected')
stopifnot(file.exists(file.path(A,'_SUCCESS')))
d<-fread(file.path(A,'tables/per_sample_kontextual.tsv'))
ids<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3))
pair_order<-unique(d$direction_id)
shorten<-function(x){for(n in c(Cardiomyocyte='CM',Fibroblast='Fib',Mononuclear_myeloid='Myel',Neutrophil='Neu',Endothelial='Endo',Pericyte='Peri')){}
 x<-gsub('Mononuclear_myeloid','Myel',x,fixed=TRUE);x<-gsub('Cardiomyocyte','CM',x,fixed=TRUE);x<-gsub('Fibroblast','Fib',x,fixed=TRUE);x<-gsub('Neutrophil','Neu',x,fixed=TRUE);x<-gsub('Endothelial','Endo',x,fixed=TRUE);x<-gsub('Pericyte','Peri',x,fixed=TRUE);gsub(' to ',' > ',x,fixed=TRUE)}
d[,pair:=factor(direction_id,levels=rev(pair_order),labels=shorten(rev(pair_order)))]
cc<-c(Cardiomyocyte='#5B9BD5',Fibroblast='#35B89A',Mononuclear_myeloid='#B57AC3',Neutrophil='#F58F8F',Endothelial='#22A7C7',Pericyte='#B5A2D8',B_cell='#60C3A8',T_cell='#9BCFBB',Unassigned='#B8BCC2')
# Complete results; external domain numbers have only within-tissue meaning.
samples<-unique(d$sample_id)
pages<-lapply(seq(1,length(samples),by=4),function(start){
 panels<-lapply(samples[start:min(start+3,length(samples))],function(id){
  x<-d[sample_id==id]
  p<-ggplot(x,aes(factor(zone),pair,fill=kontextual))+geom_tile(colour='white',linewidth=.2)+geom_text(data=x[!is.finite(kontextual)],label='×',size=2,colour='#757D86')+
   scale_fill_gradient2(low=COL_CONTRACTILE,mid='white',high=COL_FIBROTIC,limits=c(-200,max(200,max(d$kontextual,na.rm=TRUE))),na.value='#E9EBEE',breaks=c(-200,0,200,400),guide=guide_colorbar(barwidth=unit(28,'mm'),barheight=unit(2,'mm'),title.position='top',title.hjust=.5),name='L (µm)')+
   scale_x_discrete(labels=function(x)gsub('_',' ',x,fixed=TRUE))+labs(x=NULL,y=NULL)+theme_editorial_v2()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=40,hjust=1),legend.position='bottom')
  heading(p,paste(gsub('_',' ',id),unique(x$context_source),sep='\n'))
 })
 wrap_plots(list(strip_text_v2(paste('Supplementary Figure S15 | All contextual relationships',start,'–',min(start+3,length(samples))),8,'bold',0),wrap_plots(panels,ncol=2,guides='collect') & theme(legend.position='bottom'),
 note('×: not evaluable because a class is absent. Each panel is one tissue. External BANKSY labels are not homologous anatomical regions.')),
 ncol=1,heights=unit.c(unit(7,'mm'),unit(1,'null'),unit(6,'mm')))
})
export_supplement(pages,'Figure_S15',225)
