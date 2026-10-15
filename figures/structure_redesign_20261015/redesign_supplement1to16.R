source('11_图件结构重设计_20261015/scripts/design_common.R')
options(warn=1)
args<-commandArgs(TRUE);targets<-if(length(args))as.integer(strsplit(args[1],',')[[1]])else 1:16
ragg::agg_capture(width=183,height=230,units='mm',res=300)
suppsrc<-file.path(SRC,'supplement_plots')
oldcols<-c('#4A90D9','#428BD0','#397FC2','#F4777F','#EE6E79','#ED6A75','#9A78C6','#A48ACB','#36B6A4','#39B7A5','#78C7DB','#AEB7C2','#5B9BD5','#22A7C7','#35B89A','#D97FB2','#B57AC3','#F58F8F','#7B6BC2','#60C3A8','#9BCFBB','#B5A2D8','#81CCEA')
newcols<-c('#2B83D5','#2B83D5','#2B83D5','#F06D7A','#F06D7A','#EF637B','#8568C6','#8568C6','#19AA98','#19AA98','#54BFD4','#B8C2CA',unname(CELL[c('Cardiomyocyte','Endothelial','Fibroblast','Macrophage','Monocyte','Neutrophil','Dendritic','B_cell','T_NK','Pericyte','Schwann_cell')]))
recolour<-function(x){i<-match(toupper(x),toupper(oldcols));x[!is.na(i)]<-newcols[i[!is.na(i)]];x}
clean_text<-function(s){s<-gsub('Frozen','Fixed',s,fixed=TRUE);s<-gsub('frozen','fixed',s,fixed=TRUE);s<-gsub('Saved ','',s,fixed=TRUE);s<-gsub('saved ','',s,fixed=TRUE);s<-gsub('NEGATIVE_COMPLETE','no consistent module response',s,fixed=TRUE);s}
text_content<-function(g){out<-character();if(inherits(g,'text'))out<-as.character(g$label);if(!is.null(g$children))for(z in g$children)out<-c(out,text_content(z));if(!is.null(g$grobs))for(z in g$grobs)out<-c(out,text_content(z));out}
fix_grob<-function(g){
 if(inherits(g,'text')){g$label<-clean_text(g$label);g$gp$fontfamily<-FONT;if(!is.null(g$gp$fontsize))g$gp$fontsize<-max(7,g$gp$fontsize);g$gp$col<-'#25343E'}
 if(!is.null(g$children))for(i in seq_along(g$children))g$children[[i]]<-fix_grob(g$children[[i]])
 if(!is.null(g$grobs))for(i in seq_along(g$grobs))g$grobs[[i]]<-fix_grob(g$grobs[[i]])
 g
}
restyle<-function(p){
 if(inherits(p,'patchwork')){for(i in seq_len(length(p)))p[[i]]<-restyle(p[[i]]);w<-p$patches$layout$widths;if(is.numeric(w)&&length(w)==5L&&isTRUE(all.equal(as.numeric(w),c(.2,1,1,1,.31))))p<-p+plot_layout(widths=c(.42,1,1,1,.20));return(p+plot_annotation(theme=theme(plot.margin=margin(0,0,0,0))))}
 if(inherits(p,'wrapped_patch')){
  g<-attr(p,'grobs')$full
  if(!is.null(g)){
   r<-if(!is.null(g$children))Filter(function(z)inherits(z,'rastergrob'),g$children)else list()
   if(length(r)==1L&&length(dim(r[[1]]$raster))>=2L&&length(text_content(g))==0){im<-r[[1]]$raster;g<-gTree(im=im,ratio=dim(im)[2]/dim(im)[1],cl='fitRaster')}else g<-fix_grob(g)
   return(wrap_elements(full=g))
  }
  return(p)
 }
 if(inherits(p,'ggplot')){
  axis_blank<-inherits(p$theme$axis.text,'element_blank');title_blank<-inherits(p$theme$axis.title,'element_blank')
  p<-p+theme(text=element_text(family=FONT),axis.text=if(axis_blank)element_blank()else element_text(size=7,colour='#25343E'),axis.title=if(title_blank)element_blank()else element_text(size=7.2),legend.title=element_text(size=7.1),legend.text=element_text(size=7),strip.text=element_text(size=7.2),plot.margin=margin(1.2,1.2,1.2,1.2,'mm'))
  for(i in seq_along(p$layers)){
   if(!is.null(p$layers[[i]]$aes_params$family))p$layers[[i]]$aes_params$family<-FONT
   if(inherits(p$layers[[i]]$geom,'GeomText')&&!is.null(p$layers[[i]]$aes_params$size))p$layers[[i]]$aes_params$size<-max(2.35,p$layers[[i]]$aes_params$size)
   # Original raster overlays retain their source palette; matching legends do too.
  }
  if(is.data.frame(p$data)&&all(c('id','x','y')%in%names(p$data))&&length(unique(p$data$y))==1L){g<-ggplotGrob(p);ix<-which(g$layout$name=='guide-box-bottom');if(length(ix))return(wrap_elements(full=g$grobs[[ix]]))}
  return(p)
 }
 p
}
rows_of<-function(p)lapply(seq_len(length(p)),function(i)p[[i]])
grid_blocks<-function(p,ncol){z<-rows_of(p);rr<-split(z,ceiling(seq_along(z)/ncol));stack_rows(lapply(rr,function(row)side_by_side(row)),rep(1,length(rr)))}
# Retain every numbered panel; change grouping, plot areas and independent gutters.
repage<-function(p,id,j){p<-restyle(p);z<-rows_of(p);n<-length(z);h<-NULL
 if(id=='S1'){z[[2]]<-grid_blocks(z[[2]],2);h<-c(8,173,8)}
 if(id=='S2'){
  if(j==1){z[[2]]<-grid_blocks(z[[2]],4);h<-c(8,141,9,36,8)}
  if(j==2){z[[2]][[2]][[1]]<-grid_blocks(z[[2]][[2]][[1]],2);h<-c(8,104,57,31,8)}
  if(j==3){h<-c(8,34,29,29,96,8)}
 }
 if(id=='S3'){z[[2]]<-grid_blocks(z[[2]],1);h<-c(8,200,if(n==4)12 else 8);if(n==4)h<-c(8,196,12,8)}
 if(id%in%c('S4','S5')){if(n==4){z[[2]]<-grid_blocks(z[[2]],if(id=='S4')2 else 2);h<-c(8,190,10,8)}else if(n==3){z[[2]]<-grid_blocks(z[[2]],2);h<-c(8,188,8)}}
 if(id=='S6'){
  if(j==1){
   z[[2]][[2]][[2]][[1]]<-z[[2]][[2]][[2]][[1]]+scale_x_continuous('Nuclei',breaks=c(0,2500,5000),labels=c('0','2.5k','5k'),expand=expansion(mult=c(.03,.1)))
   z[[2]]<-side_by_side(rows_of(z[[2]]),c(.75,1.5,1));z[[3]]<-side_by_side(rows_of(z[[3]]))
   z[[4]][[2]][[1]]<-z[[4]][[2]][[1]]+scale_x_discrete(labels=function(x)sub('^SRR21777','',x))+labs(x='Nucleus sample (SRR21777 suffix)')+theme(axis.text.x=element_text(angle=0,hjust=.5),legend.position='bottom')
   h<-c(8,59,42,42,36,29,8)
  }
  if(j%in%2:3){h<-c(8,201,8)}
  if(j>=4){z[[4]]<-grid_blocks(z[[4]],1);h<-c(8,5,6,190,16,8)}
 }
 if(id=='S7'){h<-if(j<=2)c(8,176,12,8)else if(j==3)c(8,72,43,43,8)else c(8,158,8);if(length(h)!=n)h<-NULL}
 if(id=='S8'){z[[3]]<-grid_blocks(z[[3]],2);h<-c(8,6,170,12,8)}
 if(id=='S9'){
  if(j==1){
   z[[3]][[2]][[1]]<-z[[3]][[2]][[1]]+guides(size=guide_legend(order=2,override.aes=list(fill='#A3B5C2',colour='#576975')))+theme(legend.box='vertical')
   h<-c(8,51,78,40,8)
  }else if(n==4)h<-c(8,183,12,8)else if(n==3)h<-c(8,185,8)
 }
 if(id=='S10'){if(n==3)h<-c(8,185,8)else if(n==4)h<-c(8,170,13,8)}
 if(id=='S11'){if(n==3)h<-c(8,185,8)}
 if(id=='S12'){
  if(j==1){
   gplots<-z[[4]][[1]][[2]][[1]]
   ds<-lapply(seq_len(length(gplots)),function(i){q<-gplots[[i]][[2]];while(inherits(q,'patchwork'))q<-q[[1]];q$data})
   stopifnot(length(ds)==2L,all(vapply(ds,nrow,1L)==30L),all(unlist(lapply(ds,function(x)as.character(x$ordered_sign_pair)))=='+,+'))
   compact<-stack_rows(list(strip('Day 3: 30/30 tissue–axis pairs are +,+',7.5),strip('Day 7: 30/30 tissue–axis pairs are +,+',7.5),strip('All 10 axes in all three tissues per stage.\nOrdered gradient signs are descriptive.',7)),c(1,1,1.3))
   z[[4]][[1]]<-block(compact,'g','Gradient context')
   z[[4]]<-side_by_side(rows_of(z[[4]]),c(.44,.56));h<-c(8,43,72,57,8)
  }else h<-c(8,24,83,83,8)
 }
 if(id=='S13'){h<-switch(as.character(j),'1'=c(8,44,44,102,10),'2'=c(8,96,10,83,9),'3'=c(8,102,83,10))}
 if(id=='S14'){if(j==1){z<-c(list(strip('Supplementary Figure S14 | Original-reference sensitivity',9,'bold')),z);n<-length(z);h<-c(8,40,78,31,35,35)}else if(n==4)h<-c(8,180,12,8)else if(n==3)h<-c(8,185,8)}
 if(id=='S15'){h<-c(8,if(j==6)106 else 174,8);z[[2]]<-(z[[2]]+plot_layout(guides='collect')) & theme(legend.position='bottom')}
 if(id=='S16'){h<-if(j==1)c(8,194,10)else c(8,53,46,82,8)}
 if(is.null(h)||length(h)!=n){orig<-p$patches$layout$heights;h<-tryCatch(as.numeric(orig),error=function(e)NULL);if(is.null(h)||length(h)!=n||any(!is.finite(h)|h<=0))h<-c(8,rep(190/max(1,n-2),max(0,n-2)),8);h<-h/sum(h)*215}
 stopifnot(length(h)==n,all(h>0))
 list(plot=stack_rows(z,h),height=sum(h)+4)
}
for(i in targets){
 id<-paste0('Figure_S',i);f<-file.path(suppsrc,paste0(id,'.rds'));stopifnot(file.exists(f));d<-readRDS(f);pages<-d$pages
 if(i==6)pages<-c(pages,readRDS(file.path(suppsrc,'Figure_S6_extended.rds'))$pages)
 if(i==9)pages<-c(readRDS(file.path(SRC,'Figure_S9_models_source_plots.rds'))$pages,pages)
 if(i==14)pages<-c(readRDS(file.path(suppsrc,'Figure_S14_summary.rds'))$pages,pages)
 outpdf<-character();for(j in seq_along(pages)){
  p<-repage(pages[[j]],paste0('S',i),j);stem<-paste0(id,'_p',j)
  export(p$plot,stem,p$height,section='Supplementary_figures');outpdf<-c(outpdf,file.path(OUT,'Supplementary_figures',paste0(stem,'.pdf')))
 }
 qpdf::pdf_combine(outpdf,output=file.path(OUT,'Supplementary_figures',paste0(id,'.pdf')))
 rm(d,pages,p);gc();cat('COMPLETED',id,'\n')
}
dev.off()
