root<-'/home/zhenhanli/myocardial/scop/GSE214611_spatial'
.libPaths(c(file.path(root,'21_scop_expansion/_runtime/spacet/R-4.5/site-library'),.libPaths()))
source(file.path(root,'22_JTM_completion_20261014/figures_structured/script/figure_remaining_helpers.R'))
A<-file.path(root,'22_JTM_completion_20261014/P06_SECACT_parallel')
stopifnot(file.exists(file.path(A,'_SUCCESS')))
a<-fread(file.path(A,'tables/secact_activity.tsv.gz'))
c<-fread(file.path(A,'tables/secact_ccc.tsv'))
p<-fread(file.path(A,'tables/secact_patterns.tsv'))
v<-unique(fread(file.path(A,'tables/secact_velocity_sample_effects.tsv'))[,.(sample_id,stage,ligand_secact_human,signal_mode,n_vectors,median_vector_length)])
ligands<-sort(unique(v$ligand_secact_human));stopifnot(length(ligands)==6L,nrow(v)==72L)
st_ids<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3))
sn_ids<-c('SRR21777596','SRR21777597','SRR21777598','SRR21777581','SRR21777582','SRR21777599')
labels_st<-setNames(c(paste0('D3r',1:3),paste0('D7r',1:3)),st_ids)
labels_sn<-setNames(c('D3 596','D3 597','D3 598','D7 581','D7 582','D7 599'),sn_ids)
activity_heat<-function(x,ids,labels,title,lim=NULL){
 x<-x[secretedProtein%in%ligands,.(sample_id,secretedProtein,zscore)]
 stopifnot(!anyDuplicated(x[,.(sample_id,secretedProtein)]))
 x<-merge(CJ(sample_id=ids,secretedProtein=ligands),x,by=c('sample_id','secretedProtein'),all.x=TRUE)
 x[,sample:=factor(sample_id,levels=ids)];x[,ligand:=factor(secretedProtein,levels=rev(ligands))]
 if(is.null(lim))lim<-max(abs(x$zscore),na.rm=TRUE)
 h<-ggplot(x,aes(sample,ligand,fill=zscore))+geom_tile(colour='white',linewidth=.3)+
 geom_text(aes(label=ifelse(is.na(zscore),'–',sprintf('%.1f',ifelse(abs(zscore)<.05,0,zscore)))),size=6.7/GG_PT,family=FONT_V2)+
 scale_fill_gradient2(low=COL_CONTRACTILE,mid='white',high=COL_FIBROTIC,limits=c(-lim,lim),na.value='#ECEFF2',name='Activity z',
 guide=guide_colorbar(barwidth=unit(20,'mm'),barheight=unit(2,'mm'),title.position='top',title.hjust=.5))+scale_x_discrete(NULL,labels=labels)+labs(y=NULL)+quiet_theme+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='bottom',axis.text.x=element_text(angle=40,hjust=1))
 heading(h,title)
}
sn<-a[modality=='snRNA' & celltype_or_scope=='Mononuclear_myeloid']
st<-a[modality=='Visium']
activity_pair<-wrap_plots(list(activity_heat(sn,sn_ids,labels_sn,'Myeloid-labelled reference subset'),activity_heat(st,st_ids,labels_st,'Within-tissue median activity')),ncol=2)
cc<-c[record_type=='inferred_edge' & secretedProtein%in%ligands,.N,by=.(sample_id,secretedProtein)]
cc<-merge(CJ(sample_id=sn_ids,secretedProtein=ligands),cc,all.x=TRUE,by=c('sample_id','secretedProtein'));cc[is.na(N),N:=0L]
cc[,sample:=factor(sample_id,levels=sn_ids)]
cccplot<-ggplot(cc,aes(sample,factor(secretedProtein,levels=rev(ligands)),fill=N))+geom_tile(colour='#E3E8ED',linewidth=.3)+
 geom_text(aes(label=N),size=6.7/GG_PT,family=FONT_V2)+scale_fill_gradient(low='white',high='#5798A9',limits=c(0,max(1,cc$N)),guide='none')+
 scale_x_discrete(NULL,labels=labels_sn)+labs(y=NULL)+quiet_theme+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=40,hjust=1))
meta<-p[,.(Patterns=uniqueN(pattern_index),Proteins=uniqueN(secretedProtein)),by=sample_id]
stopifnot(setequal(meta$sample_id,st_ids),all(meta$Patterns%in%2:5))
ml<-melt(meta,id.vars='sample_id',variable.name='metric');ml[,sample:=factor(sample_id,levels=st_ids)]
pattern_summary<-ggplot(ml,aes(sample,value,colour=grepl('D7',sample_id)))+geom_point(size=1.8)+facet_wrap(~metric,ncol=2,scales='free_y')+
 scale_y_continuous(breaks=function(x){b<-pretty(x,n=4);b[b==floor(b)]})+
 scale_colour_manual(values=c('FALSE'=COL_D3,'TRUE'=COL_D7),guide='none')+scale_x_discrete(NULL,labels=labels_st)+labs(y=NULL)+quiet_theme+
 theme(axis.text.x=element_text(angle=45,hjust=1))
v[,`:=`(sample=factor(sample_id,levels=st_ids),ligand=factor(ligand_secact_human,levels=ligands),stage=factor(stage,levels=c('Day3','Day7')))]
v<-merge(v,unique(st[,.(sample_id,n_spots=n_input_units)]),by='sample_id',all.x=TRUE)
stopifnot(all(v$n_vectors<=v$n_spots),all(v$n_spots>0))
v[,vector_fraction:=n_vectors/n_spots]
vp<-ggplot(v,aes(sample,vector_fraction,colour=stage,shape=signal_mode))+geom_point(position=position_dodge(width=.45),size=1.7)+
 facet_wrap(~ligand,ncol=3)+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_colour_manual(values=c(Day3=COL_D3,Day7=COL_D7))+scale_x_discrete(NULL,labels=labels_st)+
 labs(y='Fraction of spots with vectors',shape='Mode',colour='Stage')+quiet_theme+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom')
page<-wrap_plots(list(block(activity_pair,'a','Activity for six frozen ligands','SRR suffixes identify nuclei samples; modalities are not animal-paired'),
 wrap_plots(list(block(cccplot,'b','Candidate-ligand communication','Counts across all source–receiver class pairs'),block(pattern_summary,'c','Tissue-specific patterns')),ncol=2),
 block(vp,'d','Sending and receiving vector coverage','One tissue per point; native arrow lengths are not compared'),
 note('Activity and vector lengths are model quantities, not measured secretion, receptor activation or cellular movement.')),
 ncol=1,heights=c(66,48,101,6))
if(wants('Figure_8'))export_remaining(page,'Figure_8')
classes<-unique(a[modality=='snRNA',celltype_or_scope]);lim<-max(abs(a[modality=='snRNA' & secretedProtein%in%ligands,zscore]))
allclass<-lapply(classes,function(ct)activity_heat(a[modality=='snRNA' & celltype_or_scope==ct],sn_ids,labels_sn,gsub('_',' ',ct),lim))
page1<-wrap_plots(list(strip_text_v2('Supplementary Figure S16 | All reference classes',8,'bold',0),(wrap_plots(allclass,ncol=2,guides='collect') & theme(legend.position='bottom')),
 note('Frozen eight-class reference subset; a common activity scale is used across classes.\nDash: no retained nuclei in that sample and class; no activity estimated.')),ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null'),unit(8,'mm')))
cc_all<-merge(data.table(sample_id=sn_ids),c[record_type=='inferred_edge',.N,by=sample_id],by='sample_id',all.x=TRUE);cc_all[is.na(N),N:=0L];cc_all[,sample:=factor(sample_id,levels=sn_ids)]
edgeplot<-ggplot(cc_all,aes(sample,N,fill=grepl('581|582|599',sample_id)))+geom_col(width=.7)+scale_fill_manual(values=c(COL_D3,COL_D7),guide='none')+
 scale_x_discrete(NULL,labels=labels_sn)+labs(y='All native inferred edges')+quiet_theme
elig<-c[record_type=='celltype_eligibility'];elig[,sample:=factor(sample_id,levels=sn_ids)]
eligplot<-ggplot(elig,aes(sample,sender,fill=n_nuclei))+geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=n_nuclei),size=6.7/GG_PT)+
 scale_fill_gradient(low='white',high='#80A8CE',name='Nuclei')+scale_x_discrete(NULL,labels=labels_sn)+scale_y_discrete(labels=function(x)gsub('_',' ',x,fixed=TRUE))+labs(y=NULL)+quiet_theme+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=40,hjust=1))
page2<-wrap_plots(list(strip_text_v2('Supplementary Figure S16 | Spatial summaries and input support',8,'bold',0),
 activity_heat(st,st_ids,labels_st,'All six spatial tissues'),
 wrap_plots(list(heading(cccplot,'Candidate-ligand edges'),heading(edgeplot,'All native edges')),ncol=2),
 heading(eligplot,'Retained nuclei per sample and class; source-step eligibility requires ≥30'),
 note('Zero edge counts mean no retained native edge. Source-step eligibility and model activity are distinct quantities.')),
 ncol=1,heights=c(6,64,60,82,7))
if(wants('Figure_S16'))export_supplement(list(page1,page2),'Figure_S16',225)
# Show every independently fitted spatial factor, without aligning factor numbers.
pattern_source<-fread(file.path(A,'tables/secact_pattern_spot_loadings.tsv.gz'))
stopifnot(all(is.finite(pattern_source$loading)),setequal(pattern_source$sample_id,st_ids))
pattern_pages<-lapply(c('D3','D7'),function(stage){
 ids<-st_ids[grepl(stage,st_ids)];ncols<-max(meta[sample_id%in%ids,Patterns]);rows<-lapply(ids,function(id){
 h<-pattern_source[sample_id==id];ranks<-sort(unique(h$pattern_index))
 stopifnot(length(ranks)==meta[sample_id==id,Patterns],identical(as.integer(ranks),seq_along(ranks)))
 plots<-lapply(ranks,function(i){u<-h[pattern_index==i]
  stopifnot(!anyDuplicated(u$barcode))
  z<-data.frame(x=u$x_um,y=u$y_um,score=if(max(u$loading)>0)u$loading/max(u$loading)else rep(0,nrow(u)))
  q<-ggplot(z,aes(x,y,colour=score))+geom_point(size=.28)+coord_equal()+scale_colour_gradientn(colours=c('#EDF3F8','#77B7C6','#234F70'),limits=c(0,1),guide='none')+
  theme_void()+theme(plot.margin=margin(0,0,0,0))
  heading(q,paste0(labels_st[id],' / factor ',i))})
 wrap_plots(plots,ncol=length(plots))
 })
 wrap_plots(c(list(strip_text_v2(paste('Supplementary Figure S17 |',stage,'spatial patterns'),8,'bold',0)),rows,
 list(note('Colour: H loading divided by its factor-specific maximum (0–1). Each row shows all factors selected for that tissue; factors are not aligned across rows.'))),
 ncol=1,heights=c(7,65,65,65,8))
})
if(wants('Figure_S17'))export_supplement(pattern_pages,'Figure_S17',225)

# Preserve all 72 native display fields, with no anatomical reorientation beyond
# the backend's own coord_flip and no scientific interpretation of rescaled lengths.
vector_source<-fread(file.path(A,'tables/secact_vector_fields_native_display.tsv.gz'))
background_source<-fread(file.path(A,'tables/secact_vector_background_native_display.tsv.gz'))
stopifnot(uniqueN(background_source[,.(sample_id,ligand,signal_mode)])==72L)
vector_pages<-lapply(ligands,function(g){
 blocks<-lapply(c('sending','receiving'),function(mode){
  panels<-lapply(st_ids,function(id){
   xy<-background_source[sample_id==id & ligand==g & signal_mode==mode]
   sg<-vector_source[sample_id==id & ligand==g & signal_mode==mode]
   stopifnot(nrow(xy)>0L,all(c('x','y','value')%in%names(xy)))
   q<-ggplot(xy,aes(x,y))+geom_point(colour='#DFE7EC',size=.17)+theme_void()+
     coord_flip()+theme(plot.margin=margin(0,0,0,0),aspect.ratio=diff(range(c(xy$x,sg$x_start,sg$x_end)))/diff(range(c(xy$y,sg$y_start,sg$y_end))))
   if(nrow(sg))q<-q+geom_segment(data=sg,aes(x=x_start,y=y_start,xend=x_end,yend=y_end),inherit.aes=FALSE,
    colour='#278DAB',linewidth=.13,arrow=arrow(length=unit(.45,'mm')))
   heading(q,labels_st[id])
  })
  block(wrap_plots(panels,ncol=3),if(mode=='sending')'a'else'b',if(mode=='sending')'Sending fields'else'Receiving fields')
 })
 wrap_plots(c(list(strip_text_v2(paste('Supplementary Figure S18 |',g,'· paired sending and receiving'),8,'bold',0)),blocks,
 list(note('All modeled spots are shown in gray. Native x/y components were independently rescaled; arrow lengths are not activity magnitudes.'))),
 ncol=1,heights=c(7,110,110,8))
})
if(wants('Figure_S18'))export_supplement(vector_pages,'Figure_S18',245)
