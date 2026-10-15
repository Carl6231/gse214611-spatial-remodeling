#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(scop);library(Seurat);library(data.table);library(ggplot2);library(patchwork)})
a<-commandArgs(TRUE);root<-a[1];out<-a[2];stopifnot(file.exists(file.path(out,'_SUCCESS')))
s<-fread(file.path(out,'tables/sample_status.tsv'));samples<-s$sample
components<-c('NonSE',paste0('SE',1:9));pal<-setNames(c('#AAB5C0','#168FAE','#43B6C5','#47B49C','#8DCFA7','#6B91D4','#9383C6','#C2A4DE','#ED8596','#F5AD8E'),components)
short<-function(x)gsub('_',' ',gsub('NeedlePass','Needle',x))
base<-theme_classic(base_size=8,base_family='sans')+theme(plot.title=element_text(face='bold',size=9),legend.title=element_text(size=8),legend.text=element_text(size=7))
t<-fread(file.path(out,'tables/tissue_summary.tsv'));t[,sample:=factor(sample,levels=rev(samples))];t[,component:=factor(component,levels=components)]
pheat<-ggplot(t,aes(component,sample,fill=mean_abundance))+geom_tile(colour='white',linewidth=.4)+geom_text(aes(label=sprintf('%.2f',mean_abundance)),size=2.2,colour='#234552')+scale_fill_gradientn(colours=c('#F4F7FA','#BDE6E5','#47B5B2','#176C8D'),limits=c(0,1),name='Mean\nfraction')+scale_y_discrete(labels=short)+labs(x=NULL,y=NULL,title='a  Tissue-level mean reference fractions')+base+theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())
plots<-list();detailfile<-file.path(out,'figures/SpatialEcoTyper_AllTissues_Components.pdf')
if(length(a)<3 || a[3]!='overview'){
grDevices::cairo_pdf(detailfile,width=7.48,height=9.84,onefile=TRUE)
for(k in samples){
 z<-readRDS(file.path(out,'objects',paste0(k,'_SpatialEcoTyper.rds')))
 o<-readRDS(file.path(root,'objects/phase_09/deconvolution',paste0(k,'_phase09_expanded.rds')))
 stopifnot(identical(rownames(z$abundance),colnames(o)))
 o$Eco_x<-z$coordinates$x;o$Eco_y<-z$coordinates$y
 dp<-SpatialSpotPlot(o,values=as.data.frame(z$abundance),coord.cols=c('Eco_x','Eco_y'),overlay_image=FALSE,flip.y=TRUE,pt.size=.5,pt.alpha=1,stroke=0,combine=FALSE,theme_use='theme_void',legend.title='Fraction')
 dp<-lapply(seq_along(dp),function(i)suppressMessages(dp[[i]]+scale_colour_gradientn(colours=c('#F4F7FA','#BDE6E5','#47B5B2','#176C8D'),limits=c(0,1),name='Fraction',breaks=c(0,.5,1))+labs(title=components[i])+theme(plot.title=element_text(size=9,face='bold',hjust=.5),legend.title=element_text(size=8),legend.text=element_text(size=7),plot.margin=margin(3,3,3,3))))
 page<-wrap_plots(dp,ncol=2,guides='collect')+plot_annotation(title=short(k),subtitle=paste0(nrow(z$abundance),' Visium spots | native SpatialEcoTyper deconvolution'),caption='Cross-species tumor-reference projection; fractions are not cardiac ecotype identities.',theme=theme(plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=9),plot.caption=element_text(size=8))) & theme(legend.position='right')
 print(page);if(k==samples[1])ggsave(file.path(out,'figures/SpatialEcoTyper_Components_Preview.png'),page,width=190,height=250,units='mm',dpi=180,bg='white')
 rm(o);gc(verbose=FALSE)
}
dev.off()
}
plots<-lapply(samples,function(k){
 z<-readRDS(file.path(out,'objects',paste0(k,'_SpatialEcoTyper.rds')))
 d<-data.frame(x=z$coordinates$x,y=z$coordinates$y,component=factor(components[max.col(z$abundance,ties.method='first')],levels=components))
 ggplot(d,aes(x,y,colour=component))+geom_point(size=.38,alpha=1,stroke=0)+scale_colour_manual(values=pal,limits=components,drop=FALSE)+scale_y_reverse()+coord_equal()+labs(title=gsub("Day","D",gsub("rep","r",short(k))))+theme_void(base_family='sans')+theme(legend.position='none',plot.title=element_text(size=8,face='bold',hjust=.5),plot.margin=margin(2,2,2,2))
})
atlas<-wrap_plots(plots,ncol=4)+plot_annotation(title='b  Strongest reference component per spot',theme=theme(plot.title=element_text(face='bold',size=10)))
legendplot<-ggplot(data.frame(component=factor(components,levels=components)),aes(component,1,colour=component))+geom_point(size=3)+scale_colour_manual(values=pal,drop=FALSE)+guides(colour=guide_legend(nrow=2,title=NULL))+theme_void()+theme(legend.position='bottom',legend.text=element_text(size=7))
legend<-patchwork::wrap_elements(full=cowplot::get_legend(legendplot))
combined<-(pheat / wrap_elements(full=atlas) / legend)+plot_layout(heights=c(1,2.2,.25))+plot_annotation(title='SpatialEcoTyper | 12-tissue reference projection',subtitle='25,548 spots | mouse-to-human ortholog mapping | fixed tumor reference',caption='Exploratory projection only. Each row is one tissue; the strongest component is a relative maximum, not a cell-type call.',theme=theme(plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=9),plot.caption=element_text(size=7)))
ggsave(file.path(out,'figures/SpatialEcoTyper_Overview.pdf'),combined,width=190,height=250,units='mm',device=cairo_pdf,bg='white')
ggsave(file.path(out,'figures/SpatialEcoTyper_Overview.png'),combined,width=190,height=250,units='mm',dpi=180,bg='white')
cat('PLOTS_COMPLETE\n')
