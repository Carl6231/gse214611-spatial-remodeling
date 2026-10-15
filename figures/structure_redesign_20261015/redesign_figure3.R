source('11_图件结构重设计_20261015/scripts/design_common.R')
suppressPackageStartupMessages(library(Seurat))
stopifnot('CellDimPlot'%in%getNamespaceExports('scop'))
capture<-ragg::agg_capture(width=183,height=215,units='mm',res=300)
ann<-fread(file.path(SRC,'panel_a_cluster_correction.tsv'))
oldnew<-fread(file.path(SRC,'phase_09_macrophage_reallocation_by_slice.tsv'))
mark<-fread(file.path(SRC,'phase_09_cluster_marker_audit.tsv'))
spots<-fread(file.path(SRC,'phase_09_expanded_spot_proportions.tsv.gz'))
stopifnot(sum(ann$cluster_cells)==18438L,sum(ann[retain_in_refined_macrophage==TRUE,cluster_cells])==8693L,nrow(spots)==25548L)
classes<-c('Macrophage','Fibroblast','Monocyte','Endothelial','Neutrophil','Dendritic','B_cell','T_NK')
flow<-ann[,.(n=sum(cluster_cells)),by=expanded_celltype]
flow[,cl:=factor(expanded_celltype,levels=rev(classes))]
flow[,retained:=expanded_celltype=='Macrophage']
a<-ggplot(flow,aes(n,cl))+geom_col(aes(fill=expanded_celltype),width=.7)+
 geom_text(aes(label=paste0(scales::comma(n),'  (',sprintf('%.1f',100*n/18438),'%)')),hjust=-.12,size=2.25,family=FONT)+
 scale_fill_manual(values=CELL,guide='none')+scale_y_discrete(NULL,labels=CLASS_LABEL)+
 scale_x_continuous(NULL,limits=c(0,13500),expand=c(0,0))+BASE+
 theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.text.y=element_text(size=7.1),plot.margin=margin(0,1,0,0,'mm'))
a<-block(a,'a','Reference reassignment','18,438 macrophage-labelled nuclei; 47.1% retained')

balanced<-readRDS(file.path(SRC,'reference_expanded_immune_balanced.rds'))
em<-fread(file.path(SRC,'corrected_broad_embedding.tsv'))
stopifnot(ncol(balanced)==18849L,setequal(colnames(balanced),em$cell))
em<-em[match(colnames(balanced),cell)];xy<-as.matrix(em[,.(UMAP_1,UMAP_2)]);rownames(xy)<-em$cell
balanced[['FrozenBroadUMAP']]<-CreateDimReducObject(xy,key='FrozenBroadUMAP_',assay='RNA')
meta<-balanced[[]]
# Highlight the old macrophage label and its corrected identities on identical coordinates.
stopifnot('Macrophage'%in%unique(meta$celltype))
balanced$old_view<-factor(ifelse(meta$celltype=='Macrophage','Old label','Other reference'),levels=c('Other reference','Old label'))
balanced$new_view<-factor(ifelse(meta$celltype=='Macrophage',as.character(meta$P09_expanded_celltype),'Other reference'),levels=c('Other reference',classes))
umap<-function(field,colors){
 p<-scop::CellDimPlot(balanced,group.by=field,reduction='FrozenBroadUMAP',label=FALSE,show_stat=FALSE,pt.size=.22,pt.alpha=1,palcolor=unname(colors),legend.position='none',theme_use='theme_void',theme_args=list(base_family=FONT,base_size=7),title='',subtitle='',xlab='',ylab='',combine=FALSE,raster=FALSE,verbose=FALSE)[[1]]
 while(inherits(p,'patchwork'))p<-p[[1]]
 aes<-unique(c(intersect(names(p$mapping),c('fill','colour')),unlist(lapply(p$layers,function(l)intersect(names(l$mapping),c('fill','colour'))))))
 stopifnot(length(aes)==1L)
 # Retain SCOP's exact cell coordinates; remove its empty title/strip layout.
 d<-data.table(x=xy[,1],y=xy[,2],group=balanced[[field]][,1])
 ggplot(d,aes(x,y,colour=group))+geom_point(size=.2,stroke=0)+scale_colour_manual(values=colors,guide='none',drop=FALSE)+coord_equal(expand=FALSE)+theme_void(base_family=FONT)+theme(plot.margin=margin(1,1,1,1,'mm'))
}
b1<-umap('old_view',c('Other reference'='#E6EBEF','Old label'='#687F91'))
b2<-umap('new_view',c('Other reference'='#E6EBEF',CELL[classes]))
b<-block(side_by_side(list(heading(raster_plot(b1,37,37),'Original label'),heading(raster_plot(b2,37,37),'After correction'))),'b','Same nuclei, same embedding','Balanced reference: 18,849 nuclei; other lineages in grey')

ordercl<-c(1,2,3,8,11,0,4,5,6,7,9,10)
genes<-c('Csf1r','Adgre1','C1qa','Col1a1','Pdgfra','Ccr2','Ly6c2','Pecam1','Kdr','Csf3r','S100a8','Flt3','Clec9a','Cd79a','Cd3d','Tnnt2','Myh6','Mki67')
md<-unique(mark[gene%in%genes],by=c('cluster','gene'))
stopifnot(nrow(md)==length(genes)*12L,all(md$cluster_detection>=0&md$cluster_detection<=1))
md[,`:=`(gene=factor(gene,levels=genes),cl=factor(cluster,levels=rev(ordercl)))]
rl<-setNames(paste0('C',ann$P6_cluster,'  ',CLASS_LABEL[ann$expanded_celltype]),as.character(ann$P6_cluster))
c<-ggplot(md,aes(gene,cl))+geom_point(aes(size=cluster_detection,fill=log2_cpm_ratio),shape=21,stroke=.15,colour='#576975')+
 scale_size_area(name='Nuclei detected',max_size=3.5,limits=c(0,1),breaks=c(.1,.5,1),labels=c('10%','50%','100%'))+
 scale_fill_gradient2(name='log2 CPM ratio\n(cluster / rest)',low='#287DB8',mid='#F5F8F9',high='#E76882',midpoint=0,limits=c(-4,4),oob=scales::squish,breaks=c(-4,0,4),labels=c('≤−4','0','≥4'))+
 scale_y_discrete(NULL,labels=rl)+scale_x_discrete(NULL)+BASE+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=50,hjust=1,face='italic',size=7),axis.text.y=element_text(size=6.7),legend.position='right',legend.box='vertical',panel.grid.major.y=element_line(linewidth=.18,colour='#E8EEF2'),legend.key.height=unit(3,'mm'))
c<-block(c,'c','Marker evidence supporting the correction','Detection and CPM contrasts (colour capped at ±4); C0 cardiac transcripts may reflect ambient RNA')

ks<-c('Sham_Day3_rep1','MI_Day3_rep2','MI_Day7_rep2')
coordlist<-lapply(ks,function(k){o<-readRDS(file.path(ROOT,'08_SpatialEcoTyper真实运行_20261015/objects',paste0(k,'_SpatialEcoTyper.rds')));dd<-as.data.table(o$coordinates);setnames(dd,'cell','barcode');dd[,sample_key:=k];dd})
coords<-rbindlist(coordlist)
sp<-merge(coords,spots,suffixes=c('',''),by=c('sample_key','barcode'),sort=FALSE)
stopifnot(nrow(sp)==nrow(coords),!anyDuplicated(sp[,.(sample_key,barcode)]))
limits<-c(0,1)
mapplot<-function(k,field){d<-sp[sample_key==k];ggplot(d,aes(x,y,colour=.data[[field]]))+geom_point(size=.43,shape=16)+scale_y_reverse()+coord_equal(expand=FALSE)+scale_colour_gradientn(colours=c('#EDF7F5','#90DBCA','#32B8A4','#087D85'),limits=limits,trans='sqrt',breaks=c(0,.1,.3,1),name='Macrophage fraction',guide='none')+theme_void()+theme(plot.margin=margin(1,1,1,1,'mm'))}
sp[, span:=max(diff(range(x)),diff(range(y))),by=sample_key]
sp[,`:=`(xp=(x-mean(range(x)))/span,yp=-(y-mean(range(y)))/span),by=sample_key]
spm<-melt(sp,id.vars=c('sample_key','barcode','xp','yp'),measure.vars=c('old_P4_macrophage','RCTD_Macrophage'),variable.name='reference',value.name='fraction')
spm[,`:=`(sample=factor(sample_key,levels=ks,labels=LABEL[ks]),reference=factor(reference,levels=c('old_P4_macrophage','RCTD_Macrophage'),labels=c('Original','Corrected')))]
mp<-ggplot(spm,aes(xp,yp,colour=fraction))+geom_point(size=.47)+facet_grid(reference~sample)+coord_equal(xlim=c(-.53,.53),ylim=c(-.53,.53),expand=FALSE)+
 scale_colour_gradientn(colours=c('#EDF7F5','#90DBCA','#32B8A4','#087D85'),limits=limits,trans='sqrt',breaks=c(0,.1,.3,1),guide='none')+theme_void(base_family=FONT)+theme(strip.text=element_text(size=7,face='bold'),strip.text.y=element_text(angle=0),panel.spacing=unit(1.5,'mm'),plot.margin=margin(1,1,1,1,'mm'))
bar<-colour_key(c('#EDF7F5','#90DBCA','#32B8A4','#087D85'),c('0','0.1','0.3','1'),'Macrophage fraction (sqrt colour scale)',sqrt(c(0,.1,.3,1)))
d<-block(stack_rows(list(raster_plot(mp,110,61),bar),c(61,8)),'d','Spatial effect of reference correction','Identical spots; shared scale, no thresholding')
oldnew[,`:=`(sample=factor(sample_key,levels=rev(SAMPLE)),role=role(sample_key))]
e<-ggplot(oldnew,aes(y=sample))+geom_segment(aes(x=old_P4_macrophage_mean,xend=corrected_macrophage_mean,yend=sample),colour='#B7C4CC',linewidth=.6)+geom_point(aes(x=old_P4_macrophage_mean,shape='Original'),colour='#6C8190',size=1.8,fill='white',stroke=.5)+geom_point(aes(x=corrected_macrophage_mean,colour=role,shape='Corrected'),size=1.8)+
 scale_colour_manual(values=COL,guide='none')+scale_shape_manual(values=c(Original=21,Corrected=16),name=NULL)+scale_x_continuous('Mean macrophage fraction',limits=c(0,.4),breaks=c(0,.2,.4),labels=scales::label_percent(accuracy=1))+
 scale_y_discrete(NULL,labels=LABEL)+BASE+theme(legend.position='bottom',legend.text=element_text(size=6.5),axis.text.y=element_text(size=6.7),legend.key.width=unit(2.5,'mm'),plot.margin=margin(1,1,1,1,'mm'))
e<-block(e,'e','All 12 tissues','One paired estimate per tissue')
page<-stack_rows(list(side_by_side(list(a,b),c(.47,.53)),c,side_by_side(list(d,e),c(.65,.35)),strip('Reference support, method agreement and fine-state maps are retained in Supplementary Figs. S6–S7.',6.8)),c(55,66,83,4))
dev.off()
export(page,'Figure_3',height=214)
