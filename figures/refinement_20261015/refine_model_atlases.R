options(warn=1)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(svglite)})
a<-commandArgs(TRUE);root<-normalizePath(a[1]);out<-normalizePath(a[2]);dest<-file.path(out,'Supplementary_figures');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
font<-'sans';ink<-'#263641';blue<-'#397FC2';teal<-'#1CA89B';coral<-'#F4777F';purple<-'#9A78C6';sequential<-c('#F4F7FA','#BDE6E5','#47B5B2','#176C8D')
headplot<-function(text,size=8.5,bold=TRUE){wrap_elements(full=textGrob(text,x=0,y=.5,just='left',gp=gpar(fontfamily=font,fontsize=size,fontface=if(bold)'bold'else'plain',col=ink)))}
block<-function(plot,letter,title)wrap_plots(list(headplot(paste(letter,' ',title)),plot),ncol=1,heights=unit.c(unit(6,'mm'),unit(1,'null')))
export<-function(pages,id,width=183,height=225){
 cairo_pdf(file.path(dest,paste0(id,'.pdf')),width=width/25.4,height=height/25.4,onefile=TRUE,family=font)
 for(p in pages)print(p);dev.off()
 for(i in seq_along(pages)){
  stem<-file.path(dest,paste0(id,if(length(pages)>1)paste0('_p',i)else''))
  svglite(paste0(stem,'.svg'),width=width/25.4,height=height/25.4);print(pages[[i]]);dev.off()
  ragg::agg_png(paste0(stem,'_render.png'),width=width,height=height,units='mm',res=300,background='white');print(pages[[i]]);dev.off()
 }
 cat('EXPORTED',id,length(pages),'pages\n')
}
# Existing SecAct factor loadings: only the display layout changes.
x<-fread(file.path(root,'07_分析补全与组合图重排_20261014/Additional_files/SecAct_tables/secact_pattern_spot_loadings.tsv.gz'))
ids<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3))
legend<-ggplot(data.frame(x=seq(0,1,length.out=200),y=1),aes(x,y,fill=x))+geom_tile()+scale_fill_gradientn(colours=c('#EDF3F8','#77B7C6','#234F70'),limits=c(0,1),guide='none')+scale_x_continuous('Within-factor relative loading',breaks=c(0,.5,1),expand=c(0,0))+theme_classic(base_size=8,base_family=font)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(0,90,0,90))
pages<-lapply(c('D3','D7'),function(stage){
 rows<-lapply(ids[grepl(stage,ids)],function(id){
  h<-x[sample_id==id];k<-sort(unique(h$pattern_index));stopifnot(identical(as.integer(k),seq_along(k)))
  cells<-lapply(1:5,function(i){
   if(!i%in%k)return(wrap_elements(full=nullGrob()))
   u<-h[pattern_index==i];u[,relative:=if(max(loading)>0)loading/max(loading)else 0]
   ggplot(u,aes(x_um,y_um,colour=relative))+geom_point(size=.34)+coord_equal()+scale_colour_gradientn(colours=c('#EDF3F8','#77B7C6','#234F70'),limits=c(0,1),guide='none')+labs(title=paste('Factor',i))+theme_void(base_family=font)+theme(plot.title=element_text(size=8,hjust=.5),plot.margin=margin(1,1,1,1))
  })
  block(wrap_plots(cells,ncol=5),letters[match(id,ids)],paste(gsub('_',' ',id),paste0(' | ',length(k),' independently fitted factors')))
 })
 wrap_plots(c(list(headplot(paste('Supplementary Figure S17 |',stage,'tissue-local spatial factors'),10)),rows,list(legend),list(headplot('Equal map slots. Blank slots are unfitted factors; factor numbers do not align across tissues.',7.1,FALSE))),ncol=1,heights=c(8,50,50,50,9,7)) & theme(plot.margin=margin(2,2,2,2,'mm'))
})
export(pages,'Figure_S17',height=185)
# Completed SpatialEcoTyper results; no fitted coefficient is changed or recomputed.
eco<-file.path(root,'08_SpatialEcoTyper真实运行_20261015');status<-fread(file.path(eco,'tables/sample_status.tsv'));samples<-status$sample;components<-c('NonSE',paste0('SE',1:9))
z<-lapply(samples,function(k)readRDS(file.path(eco,'objects',paste0(k,'_SpatialEcoTyper.rds'))));names(z)<-samples
stopifnot(length(z)==12,sum(vapply(z,function(v)nrow(v$abundance),integer(1)))==25548)
labels<-setNames(c('Sham r1','Sham r2','MI D3 r1','MI D3 r2','MI D3 r3','MI D7 r1','MI D7 r2','MI D7 r3','MI 1 h','MI 4 h','Needle r3','Needle r4'),samples)
rolecols<-setNames(c(rep(blue,2),rep(coral,3),rep(purple,3),rep(teal,2),rep('#78C7DB',2)),samples)
summary<-fread(file.path(eco,'tables/tissue_summary.tsv'));summary[,sample:=factor(sample,levels=samples)];summary_max<-ceiling(max(summary$q75)*10)/10
map<-function(k,component=NULL){v<-z[[k]];d<-data.table(x=v$coordinates$x,y=v$coordinates$y)
 if(is.null(component)){d[,value:=factor(components[max.col(v$abundance,ties.method='first')],levels=components)]}else d[,value:=v$abundance[,component]]
 p<-ggplot(d,aes(x,y,colour=value))+geom_point(size=.52,stroke=0)+scale_y_reverse()+coord_equal()+labs(title=labels[k])+theme_void(base_family=font)+theme(plot.title=element_text(size=8,hjust=.5,face='bold',colour=rolecols[k]),plot.margin=margin(1,1,1,1),legend.position='none')
 p<-if(is.null(component))p+scale_colour_manual(values=palette,limits=components,drop=FALSE)else p+scale_colour_gradientn(colours=sequential,limits=c(0,1),breaks=c(0,.25,.5,.75,1),name='Coefficient')
 wrap_plots(list(wrap_elements(full=textGrob(labels[k],x=.5,y=.5,gp=gpar(fontfamily=font,fontsize=8,fontface='bold',col=rolecols[k]))),wrap_elements(full=ggplotGrob(p+labs(title=NULL)))),ncol=1,heights=unit.c(unit(5,'mm'),unit(1,'null')))
}
palette<-setNames(c('#AAB5C0','#168FAE','#43B6C5','#47B49C','#8DCFA7','#6B91D4','#9383C6','#C2A4DE','#ED8596','#F5AD8E'),components)
leg<-ggplot(data.frame(x=seq(0,1,length.out=200),y=1),aes(x,y,fill=x))+geom_tile()+scale_fill_gradientn(colours=sequential,limits=c(0,1),guide='none')+scale_x_continuous('Reference coefficient (same 0-1 scale on every page)',breaks=c(0,.25,.5,.75,1),expand=c(0,0))+theme_classic(base_size=8,base_family=font)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(0,100,0,100))
pages<-lapply(components,function(ct){
 formal<-wrap_plots(lapply(samples[1:8],map,component=ct),ncol=4)
 supporting<-wrap_plots(lapply(samples[9:12],map,component=ct),ncol=4)
 ss<-summary[component==ct]
 q<-ggplot(ss,aes(sample,median_abundance,colour=sample))+geom_linerange(aes(ymin=q25,ymax=q75),linewidth=.55)+geom_point(size=1.5)+scale_colour_manual(values=rolecols,guide='none')+scale_x_discrete(labels=labels)+scale_y_continuous(limits=c(0,summary_max),breaks=seq(0,summary_max,.1))+labs(x=NULL,y='Median (IQR)')+theme_classic(base_family=font,base_size=8)+theme(axis.text.x=element_text(angle=35,hjust=1),plot.margin=margin(0,3,0,3))
 wrap_plots(list(headplot(paste0('Supplementary Figure S20 | ',ct,' across all 12 tissues'),10),block(formal,'a','Formal comparisons | Sham, MI Day 3 and MI Day 7'),block(supporting,'b','Descriptive timepoints and auxiliary needle-pass tissues'),leg,block(q,'c','Within-tissue coefficient distributions'),headplot('Exploratory tumor-reference projection. IQR describes spots within a tissue, not biological uncertainty.',7.1,FALSE)),ncol=1,heights=c(8,92,48,12,46,7)) & theme(plot.margin=margin(2,2,2,2,'mm'))
})
export(pages,'Figure_S20')
# Match overview titles, cohort grouping and typography to the rest of the atlas.
t<-copy(summary);t[,sample:=factor(sample,levels=rev(samples))];t[,component:=factor(component,levels=components)]
h<-ggplot(t,aes(component,sample,fill=mean_abundance))+geom_tile(colour='white',linewidth=.25)+geom_text(aes(label=sprintf('%.2f',mean_abundance)),size=2.45,colour=ink)+scale_fill_gradientn(colours=sequential,limits=c(0,1),name='Mean\ncoefficient',breaks=c(0,.5,1))+scale_y_discrete(labels=labels)+labs(x=NULL,y=NULL)+theme_minimal(base_family=font,base_size=8)+theme(panel.grid=element_blank(),legend.key.height=unit(4,'mm'),plot.margin=margin(0,0,0,0))
catlegend<-ggplot(data.frame(ct=factor(components,levels=components)),aes(ct,1,colour=ct))+geom_point(size=2)+scale_colour_manual(values=palette,drop=FALSE,name=NULL)+guides(colour=guide_legend(nrow=2,byrow=TRUE))+theme_void(base_family=font)+theme(legend.position='bottom',legend.text=element_text(size=8))
key<-wrap_elements(full=cowplot::get_legend(catlegend))
overview<-wrap_plots(list(headplot('Supplementary Figure S19 | Exploratory reference projection',10),block(h,'a','Mean reference coefficients by physical tissue'),block(wrap_plots(lapply(samples[1:8],map),ncol=4),'b','Strongest reference component | formal comparison tissues'),block(wrap_plots(lapply(samples[9:12],map),ncol=4),'c','Strongest reference component | descriptive and auxiliary tissues'),key,headplot('25,548 spots. Relative maxima have no confidence threshold; these labels are not cardiac ecotype identities.',7.1,FALSE)),ncol=1,heights=c(8,62,92,46,12,7)) & theme(plot.margin=margin(2,2,2,2,'mm'))
export(list(overview),'Figure_S19',height=235)
