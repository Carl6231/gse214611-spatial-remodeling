source('/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_redesign/script/figure_remaining_helpers.R')
A<-file.path(ROOT,'22_JTM_completion_20261014/P04_STATIAL_corrected')
stopifnot(file.exists(file.path(A,'_SUCCESS')))
d<-fread(file.path(A,'tables/per_sample_kontextual.tsv'))
ids<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3))
pair_order<-unique(d$direction_id)
shorten<-function(x){
 x<-gsub('Mononuclear_myeloid','Myel',x,fixed=TRUE);x<-gsub('Cardiomyocyte','CM',x,fixed=TRUE);x<-gsub('Fibroblast','Fib',x,fixed=TRUE);x<-gsub('Neutrophil','Neu',x,fixed=TRUE);x<-gsub('Endothelial','Endo',x,fixed=TRUE);x<-gsub('Pericyte','Peri',x,fixed=TRUE);gsub(' to ',' > ',x,fixed=TRUE)}
d[,pair:=factor(direction_id,levels=rev(pair_order),labels=shorten(rev(pair_order)))]
cc<-c(Cardiomyocyte='#5079AD',Fibroblast='#009E73',Mononuclear_myeloid='#B7799F',Neutrophil='#E67682',Endothelial='#0072B2',Pericyte='#B08CBF',B_cell='#1B9E77',T_cell='#61B5A7',Unassigned='#B8BCC2')
maps<-lapply(ids,function(id){x<-readRDS(file.path(A,'objects',paste0(id,'_Statial.rds')))$backend$input
 x$class<-sub('@.*','',x$cellType)
 p<-ggplot(x,aes(x/1000,-y/1000,colour=class))+geom_point(size=.25)+scale_colour_manual(values=cc,drop=FALSE)+coord_equal()+theme_void()+theme(legend.position='none')
 heading(p,gsub('_',' ',sub('MI_','',id)))})
mapblock<-wrap_plots(list(wrap_plots(maps,ncol=3),legend_strip(cc,gsub('_',' ',names(cc)),rows=2)),ncol=1,heights=c(1,.10))
mi<-d[sample_id%in%ids];mi[,sample:=factor(sample_id,levels=ids,labels=c('D3r1','D3r2','D3r3','D7r1','D7r2','D7r3'))]
mi[,zone:=factor(zone,levels=c('RZ_core','Transition_corridor','Injury_core','Off_corridor'))]
grid<-CJ(sample_id=ids,direction_id=pair_order,zone=c('RZ_core','Transition_corridor','Injury_core','Off_corridor'))
mi<-merge(grid,mi,by=c('sample_id','direction_id','zone'),all.x=TRUE)
mi[,pair:=factor(direction_id,levels=rev(pair_order),labels=shorten(rev(pair_order)))]
mi[,sample:=factor(sample_id,levels=ids,labels=c('D3r1','D3r2','D3r3','D7r1','D7r2','D7r3'))]
mi[,zone:=factor(zone,levels=c('RZ_core','Transition_corridor','Injury_core','Off_corridor'))]
mi[,missing_symbol:=ifelse(is.na(status),'–','×')]
lim<-max(abs(mi$kontextual),na.rm=TRUE)
heat<-ggplot(mi,aes(sample,pair,fill=kontextual))+geom_tile(colour='white',linewidth=.2)+
 geom_text(data=mi[!is.finite(kontextual)],aes(label=missing_symbol),colour='#7A828A',size=2)+
facet_wrap(~zone,ncol=2,drop=FALSE)+scale_fill_gradient2(low=COL_CONTRACTILE,mid='white',high=COL_FIBROTIC,limits=c(-lim,lim),na.value='#E9EBEE',breaks=c(-200,0,200),guide=guide_colorbar(barwidth=unit(30,'mm'),barheight=unit(2,'mm'),title.position='top',title.hjust=.5),name='Contextual L (µm)')+
 labs(x=NULL,y=NULL)+theme_editorial_v2()+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='bottom',strip.text=element_text(size=6.8))
sel<-mi[stage=='Day7' & direction_id%in%c('Mononuclear_myeloid to Fibroblast','Fibroblast to Mononuclear_myeloid') & zone%in%c('Transition_corridor','Injury_core')]
sel[,direction:=factor(shorten(direction_id))]
pt<-ggplot(sel,aes(zone,kontextual,colour=sample,shape=direction))+geom_hline(yintercept=0,linewidth=.3,colour='#A8ADB3')+
 geom_point(position=position_dodge(width=.4),size=1.8)+scale_colour_manual(values=c(D7r1='#B69BCD',D7r2='#8460A6',D7r3='#4B306C'))+labs(x=NULL,y='Contextual L (µm)')+
 guides(colour=guide_legend(nrow=1),shape=guide_legend(nrow=1))+
 theme_editorial_v2()+theme(legend.position='bottom',legend.box='vertical',legend.title=element_blank(),axis.text.x=element_text(size=6))
coverage<-d[,.(Evaluable=sum(status=='EVALUABLE'),Unavailable=sum(status!='EVALUABLE')),by=.(sample_id,cohort)]
cv<-melt(coverage,id.vars=c('sample_id','cohort'),variable.name='status',value.name='relationships')
cv[,sample_id:=factor(sample_id,levels=unique(coverage$sample_id))]
coverageplot<-ggplot(cv,aes(sample_id,relationships,fill=status))+geom_col(width=.8)+scale_fill_manual(values=c(Evaluable='#348C90',Unavailable='#D7DDE2'))+
 labs(x=NULL,y='Records')+theme_editorial_v2()+theme(axis.text.x=element_text(angle=70,hjust=1,size=5.8),legend.position='bottom',legend.title=element_blank())
rows<-list(block(mapblock,'a','Dominant RCTD8 classes across six infarct tissues'),
 block(heat,'b','Context-specific relationships at 200 µm','×: absent class; –: absent region. Descriptive values; no significance test'),
 wrap_plots(list(block(pt,'c','Day 7 myeloid–fibroblast relationships'),block(coverageplot,'d','Evaluability across 22 tissues')),ncol=2,widths=c(1,1.15)),
 note('Mixture-dominant spots are not individual cells. CM: cardiomyocyte; Fib: fibroblast; Myel: mononuclear myeloid; Neu: neutrophil.'))
heights<-c(53,99,64,5)
ragg::agg_capture(width=183,height=225,units='mm',res=300)
grobs<-lapply(seq_along(rows),function(i){
 g<-patchworkGrob(wrap_plots(list(rows[[i]])) & theme(plot.margin=margin(0,0,0,0)))
 grobTree(g,vp=viewport(layout.pos.row=i,layout.pos.col=1))})
page<-wrap_elements(full=grobTree(children=do.call(gList,grobs),
 vp=viewport(layout=grid.layout(length(rows),1,heights=unit(heights,'null')))))
dev.off()
export_remaining(page,'Figure_7')
