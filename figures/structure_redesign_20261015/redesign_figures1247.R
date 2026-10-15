source('11_图件结构重设计_20261015/scripts/design_common.R')
ragg::agg_capture(width=183,height=225,units='mm',res=300)
xy<-readRDS(file.path(SRC,'niche_spatial_values.rds'))
# Figure 1: cohort, tissue identity and the usable reference become the overview.
cohort<-data.table(x=1:6,label=c('Sham\nDay 3\nn = 2','MI\nDay 3\nn = 3','MI\nDay 7\nn = 3','MI\n1 h\nn = 1','MI\n4 h\nn = 1','Needle-pass\nDay 3\nn = 2'),role=c('Sham','MI Day 3','MI Day 7','Descriptive','Descriptive','Needle-pass'))
pa<-ggplot(cohort,aes(x,1))+geom_rect(aes(xmin=x-.44,xmax=x+.44,ymin=.4,ymax=1.6,fill=role),colour=NA,alpha=.12)+geom_segment(aes(x=x-.44,xend=x+.44,y=1.6,yend=1.6,colour=role),linewidth=1.3)+geom_text(aes(label=label),size=2.6,lineheight=1.1,family=FONT)+scale_fill_manual(values=COL,guide='none')+scale_colour_manual(values=COL,guide='none')+annotate('text',x=2,y=.1,label='8 tissues in the formal comparisons',size=2.4,fontface='bold',family=FONT)+annotate('text',x=5,y=.1,label='2 descriptive + 2 auxiliary tissues',size=2.4,fontface='bold',family=FONT)+coord_cartesian(xlim=c(.5,6.5),ylim=c(-.05,1.65),expand=FALSE)+theme_void()
a<-block(pa,'a','Study design and analytical scope','10 GEO captures → 12 physical tissues; 176,526 nuclei in the reference')
z<-spatial_coordinates(xy);counts<-z[,.N,by=sample_key];z[,sample:=factor(sample_key,levels=SAMPLE,labels=paste0(LABEL[SAMPLE],'\n',scales::comma(counts[match(SAMPLE,sample_key),N]),' spots'))]
pb<-ggplot(z,aes(xp,yp,colour=role(sample_key)))+geom_point(size=.40)+facet_wrap(~sample,ncol=4)+coord_equal(expand=FALSE)+scale_colour_manual(values=COL,guide='none')+theme_void(base_family=FONT)+theme(strip.text=element_text(size=7.2,face='bold'),panel.spacing=unit(2,'mm'))
he_cells<-lapply(SAMPLE,function(k)heading(image_file(file.path(SRC,paste0('S1_HE_',k,'.png'))),paste0(LABEL[k],' · ',scales::comma(counts[sample_key==k,N]))))
b<-block(stack_rows(list(side_by_side(he_cells[1:4]),side_by_side(he_cells[5:8]),side_by_side(he_cells[9:12])),c(1,1,1)),'b','The complete spatial cohort','H&E tissue views; labels give retained-spot counts (total 25,548)')
pairs<-lapply(c('GSM6613080','GSM6613081'),function(gsm){side_by_side(list(heading(image_file(file.path(SRC,paste0('F1b_before_',gsm,'.png'))),'Before'),heading(image_file(file.path(SRC,paste0('F1b_resolved_',gsm,'.png'))),'Resolved')))})
c<-block(stack_rows(list(heading(pairs[[1]],'GSM6613080'),heading(pairs[[2]],'GSM6613081'),strip('Blue: Sham; cyan: needle-pass.\nAssignment uses author barcodes.',6.5)),c(1,1,.18)),'c','Mixed-capture resolution','Original → author-resolved')
ref<-fread(file.path(SRC,'Figure_1h_reference_celltypes.tsv'));ref[,type:=factor(celltype,levels=rev(celltype))]
pr<-ggplot(ref,aes(reference_cells,type,colour=celltype))+geom_point(size=1.8)+geom_text(aes(label=scales::comma(reference_cells)),hjust=-.2,size=2.2)+scale_colour_manual(values=CELL,guide='none')+scale_x_log10('Reference nuclei',expand=expansion(mult=c(.03,.22)),labels=scales::label_number(scale_cut=scales::cut_short_scale()))+scale_y_discrete(NULL,labels=CLASS_LABEL)+BASE+theme(axis.text.y=element_text(size=7))
d<-block(pr,'d','Initial reference composition','Six initial labels; Schwann cells: n = 28')
flow<-wrap_elements(full=grobTree(textGrob('25,548 retained spots',x=.05,y=.83,just='left',gp=gpar(fontfamily=FONT,fontsize=10,fontface='bold',col='#2B83D5')),textGrob('5,151 selected myeloid-rich centers',x=.05,y=.58,just='left',gp=gpar(fontfamily=FONT,fontsize=9,fontface='bold',col='#8568C6')),textGrob('3,964 assigned   |   1,187 uncertain',x=.05,y=.33,just='left',gp=gpar(fontfamily=FONT,fontsize=8.5,col='#25343E')),textGrob('22,490 genes shared by reference and spatial inputs',x=.05,y=.10,just='left',gp=gpar(fontfamily=FONT,fontsize=7.5,col='#576975'))))
e<-block(flow,'e','Corrected analysis set','Fixed selection rules; corrected reference in Figure 3')
page<-stack_rows(list(a,side_by_side(list(b,c),c(.72,.28)),side_by_side(list(d,e),c(.48,.52)),strip('Quality distributions and marker compatibility: Supplementary Fig. S2. Physical tissues are the units of inference.',6.8)),c(31,110,45,7))
export(page,'Figure_1',199)
# Figure 2: spatial distributions followed by tissue summaries and calibration.
v<-readRDS(file.path(SRC,'Figure_2_values.rds'));spots<-fread(file.path(SRC,'phase_09_expanded_spot_proportions.tsv.gz'))
sd<-merge(xy[,.(sample_key,barcode,x,y)],spots[,.(sample_key,barcode,old_P4_macrophage)],by=c('sample_key','barcode'))
keys<-SAMPLE[c(1,3,4,5,2,6,7,8)]
mp<-spatial_grid(sd,'old_P4_macrophage',keys,c('#EDF7F5','#90DBCA','#32B8A4','#087D85'),c(0,1),'Fraction\n(sqrt scale)',ncol=4,trans='sqrt')
a<-block(raster_plot(mp,179,77),'a','Macrophage-labelled abundance with the initial reference','Eight formal-comparison tissues; same 0–1 scale; corrected-reference effects in Figure 3')
bd<-copy(v$burden[sample_key%in%SAMPLE[1:8]]);bd[,group:=factor(role(sample_key),levels=c('Sham','MI Day 3','MI Day 7'))]
make_summary<-function(d,field,title){d<-copy(d);d[,x:=as.integer(group)+seq(-.1,.1,length.out=.N),by=group];ggplot(d,aes(x,.data[[field]],colour=group))+geom_point(size=1.65)+geom_segment(data=d[,.(x=as.integer(group[1]),med=median(get(field))),by=group],aes(x=x-.18,xend=x+.18,y=med,yend=med),inherit.aes=FALSE,colour='#25343E',linewidth=.4)+scale_colour_manual(values=COL,guide='none')+scale_x_continuous(NULL,breaks=1:3,labels=c('Sham\nn=2','D3\nn=3','D7\nn=3'))+scale_y_continuous(NULL,limits=c(0,NA),labels=scales::label_percent(accuracy=1),breaks=scales::breaks_pretty(n=3))+labs(title=title)+BASE+theme(plot.title=element_text(size=7.1,hjust=.5),axis.text.x=element_text(size=6.6))}
bs<-list(make_summary(bd,'mean_RCTD_Macrophage','Macrophage load'),make_summary(bd,'macrophage_rich_fraction','Rich-area fraction'))
for(cl in c('Cardiomyocyte','Fibroblast','Endothelial')){dt<-copy(v$composition[celltype==cl & sample_key%in%SAMPLE[1:8]]);dt[,group:=factor(role(sample_key),levels=levels(bd$group))];bs<-c(bs,list(make_summary(dt,'mean_proportion',cl)))}
b<-block(side_by_side(bs),'b','Tissue-level composition summaries','Dots: physical tissues; horizontal bars: medians; rich-area cutoff = 0.20')
co<-copy(v$composition[sample_key%in%SAMPLE[1:8]]);co[,`:=`(sm=factor(sample_key,levels=SAMPLE[1:8]),ct=factor(celltype,levels=rev(c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage','Pericyte','Schwann_cell'))))]
pc<-ggplot(co,aes(sm,ct,fill=row_z_score))+geom_tile(colour='white',linewidth=.25)+scale_fill_gradient2(low='#2B83D5',mid='white',high='#EF637B',limits=c(-2,2),oob=scales::squish,name='Row z\n(capped)')+scale_x_discrete(NULL,labels=SHORT)+scale_y_discrete(NULL,labels=CLASS_LABEL)+MATRIX+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom',legend.key.width=unit(5,'mm'),legend.key.height=unit(2,'mm'))
c<-block(pc,'c','Composition context')
meth<-rbindlist(list(v$agreement[,.(sample_key,celltype,value=spearman_rctd_vs_dwls,metric='RCTD vs SpatialDWLS')],v$markers[,.(sample_key,celltype,value=rctd_marker_spearman,metric='RCTD vs marker expression')]))
meth<-meth[sample_key%in%SAMPLE[1:8]&celltype%in%c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage')];meth[,`:=`(sm=factor(sample_key,levels=SAMPLE[1:8]),ct=factor(celltype,levels=rev(c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage'))))]
pm<-ggplot(meth,aes(sm,ct,fill=value))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=sprintf('%.2f',value)),size=2.2)+facet_wrap(~metric,ncol=1)+scale_fill_gradientn(colours=c('#F3F8F6','#AFDCCC','#44AFAC'),limits=c(-1,1),breaks=c(-1,0,1),name='Spearman ρ')+scale_x_discrete(NULL,labels=SHORT)+scale_y_discrete(NULL,labels=CLASS_LABEL)+MATRIX+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom',legend.key.width=unit(10,'mm'),legend.key.height=unit(2,'mm'),strip.text=element_text(size=7.1))
d<-block(pm,'d','Method and marker concordance')
u<-copy(v$uncertain);u[,sm:=factor(sample_key,levels=SAMPLE)]
pu<-ggplot(u,aes(sm,uncertain_fraction,colour=role(sample_key)))+geom_segment(aes(xend=sm,y=0,yend=uncertain_fraction),colour='#CCD7DE',linewidth=.5)+geom_point(size=1.6)+scale_colour_manual(values=COL,guide='none')+scale_x_discrete(NULL,labels=SHORT)+scale_y_continuous('Uncertain fraction',labels=scales::label_percent(accuracy=1))+BASE
e<-block(pu,'e','Uncertain assignments retained across all tissues')
page<-stack_rows(list(a,b,side_by_side(list(c,d),c(.42,.58)),e),c(78,45,68,30))
export(page,'Figure_2',227)
# Figure 4: show corrected niches first, then prevalence and their evidence limits.
NICH<-c('inflammatory-infarct niche'='#F06D7A','macrophage-fibroblast remodeling niche'='#8568C6','vascular-repair niche'='#19AA98','remote resident niche'='#2B83D5',Uncertain='#8899A5',Not_macrophage_center='#E4E9ED')
NL<-c('Inflammatory','Mac–fibro remodeling','Vascular repair','Remote resident','Uncertain','Not selected');names(NL)<-names(NICH)
z<-spatial_coordinates(xy[sample_key%in%SAMPLE[3:8]]);z[,sample:=factor(sample_key,levels=SAMPLE[3:8],labels=LABEL[SAMPLE[3:8]])]
mp<-ggplot(z,aes(xp,yp,colour=P09N_niche))+geom_point(size=.52,stroke=0)+facet_wrap(~sample,ncol=3)+coord_equal(expand=FALSE)+scale_colour_manual(values=NICH,guide='none',drop=FALSE)+theme_void(base_family=FONT)+theme(strip.text=element_text(size=7.2,face='bold'),panel.spacing=unit(2,'mm'))
a<-block(stack_rows(list(raster_plot(mp,108,76),key(NICH,NL,3)),c(76,12)),'a','Spatially organized niches after reference correction','All six replicated infarct tissues')
prev<-fread(file.path(SRC,'phase_09_niche_grid_slice_summary.tsv'))[setting_id=='S054'];stopifnot(sum(unique(prev[,.(sample_key,centers)])$centers)==5151L)
prev[,`:=`(sample=factor(sample_key,levels=rev(SAMPLE)),niche=factor(niche,levels=names(NICH)[1:5]))]
ct<-unique(prev[,.(sample,centers)])
pp<-ggplot(prev,aes(sample,fraction,fill=niche))+geom_col(width=.82)+coord_flip()+geom_text(data=ct,aes(sample,1.05,label=centers),inherit.aes=FALSE,size=2.15)+scale_fill_manual(values=NICH,guide='none')+scale_x_discrete(NULL,labels=SHORT)+scale_y_continuous('Fraction of selected centers',breaks=c(0,.5,1),labels=scales::label_percent(),limits=c(0,1.13),expand=c(0,0))+BASE+theme(axis.text.x=element_text(angle=45,hjust=1))
b<-block(pp,'b','Niche prevalence','Numbers: selected centers (total 5,151)')
fields<-c('P09N_inflammatory_axis','P09N_remodeling_axis','P09N_vascular_axis','P09N_remote_axis')
ax<-melt(xy[P09N_center==TRUE],id.vars=c('sample_key','P09N_niche'),measure.vars=fields,variable.name='axis',value.name='value')[is.finite(value),.(value=median(value)),by=.(sample_key,P09N_niche,axis)][,.(value=median(value)),by=.(P09N_niche,axis)]
ax[,`:=`(axis=factor(axis,levels=fields,labels=c('Inflam.','Remodel.','Vascular','Remote')),niche=factor(P09N_niche,levels=rev(names(NICH)[1:5])))]
pa<-ggplot(ax,aes(axis,niche,fill=value))+geom_tile(colour='white',linewidth=.3)+scale_fill_gradient2(low='#2B83D5',mid='white',high='#EF637B',limits=rep(max(abs(ax$value)),2)*c(-1,1),breaks=c(-2,0,2),name='Axis score')+scale_y_discrete(NULL,labels=NL)+labs(x=NULL)+MATRIX+theme(axis.text.x=element_text(angle=40,hjust=1),axis.text.y=element_text(size=6.7),legend.position='bottom',legend.key.height=unit(2,'mm'))
c<-block(pa,'c','Axis structure','Median of tissue-level niche medians')
axis<-fread(file.path(SRC,'phase_09_niche_continuous_slice_summary.tsv'));axis[,sample:=factor(sample_key,levels=SAMPLE)]
pd<-ggplot(axis,aes(sample,remodeling_axis,colour=role(sample_key)))+geom_hline(yintercept=0,colour='#CBD6DE',linewidth=.3)+geom_point(size=1.7)+scale_colour_manual(values=COL,guide='none')+scale_x_discrete(NULL,labels=SHORT)+labs(y='Remodeling axis')+BASE+theme(axis.text.x=element_text(angle=60,hjust=1,size=6.4))
d<-block(pd,'d','Remodeling across tissues')
loo<-fread(file.path(SRC,'phase_09_niche_leave_one_slice_out.tsv'))[metric=='remodeling_axis' & contrast%in%c('Day7_MI_vs_Sham','Day7_MI_vs_Day3_MI')]
loo[,contrast:=factor(contrast,levels=c('Day7_MI_vs_Day3_MI','Day7_MI_vs_Sham'),labels=c('D7 vs D3','D7 vs Sham'))]
pe<-ggplot(loo,aes(full_effect,contrast))+geom_linerange(aes(xmin=min_leave_one_out_effect,xmax=max_leave_one_out_effect),colour='#8568C6',linewidth=.8)+geom_point(colour='#8568C6',size=2)+labs(x='Effect and leave-one-out range',y=NULL)+BASE+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.title.x=element_text(size=6.5))
e<-block(pe,'f','Leave-one-tissue-out support')
coh<-fread(file.path(SRC,'phase_09_niche_spatial_null.tsv'));coh[,sample:=factor(sample_key,levels=SAMPLE)]
pf<-ggplot(coh[is.finite(enrichment)],aes(sample,enrichment,colour=role(sample_key)))+geom_hline(yintercept=0,colour='#CBD6DE',linewidth=.3)+geom_segment(aes(xend=sample,y=0,yend=enrichment),colour='#CCD7DE',linewidth=.5)+geom_point(size=1.6)+geom_text(data=coh[!is.finite(enrichment)],aes(sample,0,label='NE'),inherit.aes=FALSE,size=2.2)+scale_colour_manual(values=COL,guide='none')+scale_x_discrete(NULL,labels=SHORT,drop=FALSE)+labs(y='Observed − permuted')+BASE+theme(axis.text.x=element_text(angle=50,hjust=1,size=6.6))
f<-block(pf,'e','Spatial coherence','NE: no center-to-center edges')
page<-stack_rows(list(side_by_side(list(a,b),c(.61,.39)),side_by_side(list(c,d,f),c(1.05,1,1)),e,strip('Selection: myeloid ≥0.20; radius ×1.2; program ≥0.10; margin ≥0.15. Niches are rule-defined spatial assignments.',6.7)),c(103,57,29,7))
export(page,'Figure_4',202)
# Figure 7: spatial contexts, directional heterogeneity, and complete evaluability.
v<-readRDS(file.path(SRC,'Figure_7_values.rds'));d<-copy(v$mi)
# Use the native saved class coordinates, with translation/scale only for page fitting.
mp<-side_by_side(lapply(v$maps,function(p) p))
maps<-lapply(v$maps,function(p){pp<-p[[2]];while(inherits(pp,'patchwork'))pp<-pp[[1]];dd<-as.data.table(pp$data);dd})
cc<-c(Cardiomyocyte='#2B8BC9',Fibroblast='#20AE91',Mononuclear_myeloid='#9270CD',Neutrophil='#F28A7B',Endothelial='#23ACBF',Pericyte='#B591D5',B_cell='#63BF93',T_cell='#A5CC86',Unassigned='#B8C2CA')
mapplots<-lapply(seq_along(maps),function(i){q<-ggplot(maps[[i]],aes(x,-y,colour=class))+geom_point(size=.42)+coord_equal()+scale_colour_manual(values=cc,guide='none')+theme_void();heading(raster_plot(q,33,27),c(paste0('D3 r',1:3),paste0('D7 r',1:3))[i])})
a<-block(stack_rows(list(side_by_side(mapplots[1:3]),side_by_side(mapplots[4:6])),c(1,1)),'a','Mixture-dominant spot classes')
sel<-copy(v$sel);sel[,direction:=factor(direction,levels=levels(direction),labels=c('Fib → Myel','Myel → Fib'))]
sel[,zone:=factor(zone,levels=c('Transition_corridor','Injury_core'),labels=c('Transition','Injury core'))]
pb<-ggplot(sel,aes(zone,kontextual,group=sample,colour=sample))+geom_hline(yintercept=0,colour='#BAC7CF',linewidth=.3)+geom_line(linewidth=.55)+geom_point(size=1.8)+facet_wrap(~direction,ncol=2)+scale_colour_manual(values=c(D7r1='#B29DDA',D7r2='#8568C6',D7r3='#5F439D'),name=NULL)+labs(x=NULL,y='Contextual L (µm)')+BASE+theme(legend.position='bottom',axis.text.x=element_text(size=6.5),strip.text=element_text(size=7))
b<-block(pb,'b','Day 7 myeloid–fibroblast context','One line per tissue; directions shown separately')
lim<-max(abs(d$kontextual),na.rm=TRUE)
ph<-ggplot(d,aes(sample,pair,fill=kontextual))+geom_tile(colour='white',linewidth=.25)+geom_text(data=d[!is.finite(kontextual)],aes(label=missing_symbol),size=2.2,colour='#82929D')+facet_wrap(~zone,ncol=2,labeller=as_labeller(c(RZ_core='Remote-zone core',Transition_corridor='Transition corridor',Injury_core='Injury core',Off_corridor='Off corridor')))+scale_fill_gradient2(low='#2B83D5',mid='white',high='#EF637B',limits=c(-lim,lim),na.value='#E9EEF1',name='Contextual L (µm)')+labs(x=NULL,y=NULL)+MATRIX+theme(legend.position='bottom',legend.key.width=unit(9,'mm'),legend.key.height=unit(2,'mm'),axis.text.y=element_text(size=6.5))
c<-block(ph,'c','Complete contextual relationships at 200 µm','×: absent class; –: absent region; descriptive values without significance tests')
cv<-copy(v$coverage);pco<-ggplot(cv,aes(sample_id,relationships,fill=status))+geom_col(width=.8)+scale_fill_manual(values=c(Evaluable='#19AA98',Unavailable='#D7E1E6'),name=NULL)+labs(x=NULL,y='Relationship records')+scale_x_discrete(labels=function(x){x<-gsub('Day3_Needle_pass','Needle',x,fixed=TRUE);x<-gsub('Day3_Sham','Sham',x,fixed=TRUE);x<-gsub('WT_MI_','WT ',x,fixed=TRUE);x<-gsub('MI_','',x,fixed=TRUE);gsub('_',' ',x)})+BASE+theme(axis.text.x=element_text(angle=55,hjust=1,size=6.5),legend.position='right')
dp<-block(pco,'d','Evaluability across all 22 tissues')
page<-stack_rows(list(side_by_side(list(a,b),c(.55,.45)),key(cc,gsub('_',' ',names(cc)),2),c,dp,strip('CM: cardiomyocyte; Fib: fibroblast; Myel: mononuclear myeloid; Neu: neutrophil. Mixture-dominant spots are not single cells.',6.6)),c(65,9,93,45,5))
export(page,'Figure_7',223)
dev.off()
