suppressPackageStartupMessages(suppressWarnings({library(ggplot2);library(patchwork);library(data.table);library(grid);library(scop)}))
ROOT <- normalizePath('.')
OUT <- file.path(ROOT,'11_图件结构重设计_20261015')
SRC <- file.path(OUT,'source_data')
FONT <- 'Arial'
COL <- c(Sham='#2B83D5','MI Day 3'='#F06D7A','MI Day 7'='#8568C6',Descriptive='#19AA98','Needle-pass'='#54BFD4')
CELL <- c(Cardiomyocyte='#2B8BC9',Endothelial='#23ACBF',Fibroblast='#20AE91',Macrophage='#E7688C',Monocyte='#9270CD',Neutrophil='#F28A7B',Dendritic='#596DC0',B_cell='#63BF93',T_NK='#A5CC86',Pericyte='#B591D5',Schwann_cell='#73CBE0',Uncertain='#B8C2CA')
CLASS_LABEL <- c(Cardiomyocyte='Cardiomyocyte',Endothelial='Endothelial',Fibroblast='Fibroblast',Macrophage='Macrophage',Monocyte='Monocyte',Neutrophil='Neutrophil',Dendritic='Dendritic',B_cell='B cell',T_NK='T / NK',Pericyte='Pericyte',Schwann_cell='Schwann cell',Uncertain='Uncertain')
SAMPLE <- c('Sham_Day3_rep1','Sham_Day3_rep2','MI_Day3_rep1','MI_Day3_rep2','MI_Day3_rep3','MI_Day7_rep1','MI_Day7_rep2','MI_Day7_rep3','MI_1HR_rep1','MI_4HR_rep1','NeedlePass_Day3_rep3','NeedlePass_Day3_rep4')
LABEL <- setNames(c('Sham r1','Sham r2','MI D3 r1','MI D3 r2','MI D3 r3','MI D7 r1','MI D7 r2','MI D7 r3','MI 1 h','MI 4 h','Needle r3','Needle r4'),SAMPLE)
role <- function(x) ifelse(grepl('^Sham',x),'Sham',ifelse(grepl('^MI_Day3',x),'MI Day 3',ifelse(grepl('^MI_Day7',x),'MI Day 7',ifelse(grepl('^Needle',x),'Needle-pass','Descriptive'))))
BASE <- theme_classic(base_family=FONT,base_size=7.5)+theme(axis.text=element_text(size=7,colour='#25343E'),axis.title=element_text(size=7.4),axis.line=element_line(linewidth=.3,colour='#3C4B55'),axis.ticks=element_line(linewidth=.25),axis.ticks.length=unit(1,'mm'),plot.title=element_text(size=8,face='bold'),strip.text=element_text(size=7.3,face='bold'),strip.background=element_blank(),legend.title=element_text(size=7),legend.text=element_text(size=6.8),legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'mm'),panel.spacing=unit(2,'mm'),plot.margin=margin(1,1,1,1,'mm'))
strip <- function(text,size=7.5,face='plain',align=0,col='#25343E') wrap_elements(full=textGrob(text,x=unit(align,'npc'),y=.5,just=c(if(align==0)'left' else if(align==1)'right' else 'center','center'),gp=gpar(fontfamily=FONT,fontsize=size,fontface=face,col=col)))
block <- function(p,letter,title,sub=NULL) {
 h<-grobTree(textGrob(letter,x=unit(0,'mm'),y=unit(.6,'npc'),just='left',gp=gpar(fontfamily=FONT,fontsize=11,fontface='bold',col='#172C3A')),textGrob(title,x=unit(4.5,'mm'),y=unit(.6,'npc'),just='left',gp=gpar(fontfamily=FONT,fontsize=8.1,fontface='bold',col='#172C3A')))
 if(!is.null(sub))h<-grobTree(h,textGrob(sub,x=unit(4.5,'mm'),y=unit(.08,'npc'),just=c('left','bottom'),gp=gpar(fontfamily=FONT,fontsize=6.7,col='#576975')))
 wrap_plots(list(wrap_elements(full=h),wrap_plots(list(p))),ncol=1,heights=unit.c(unit(if(is.null(sub))6 else 10,'mm'),unit(1,'null')))
}
heading <- function(p,title) wrap_plots(list(strip(title,7.1,'bold',.5),p),ncol=1,heights=unit.c(unit(4.5,'mm'),unit(1,'null')))
key <- function(cols,labels=names(cols),nrow=1) {
 d<-data.frame(cl=factor(names(cols),levels=names(cols)),x=seq_along(cols),y=1)
 p<-ggplot(d,aes(x,y,colour=cl))+geom_point()+scale_colour_manual(values=cols,labels=labels,name=NULL,drop=FALSE)+guides(colour=guide_legend(nrow=nrow,byrow=TRUE,override.aes=list(size=2)))+BASE+theme(legend.position='bottom',legend.margin=margin(0,0,0,0),legend.spacing.x=unit(1,'mm'),legend.text=element_text(size=6.8))
 g<-ggplotGrob(p);wrap_elements(full=g$grobs[[which(g$layout$name=='guide-box-bottom')]])
}
export <- function(p,id,height=215,section='Figures') {
 dir.create(file.path(OUT,section),recursive=TRUE,showWarnings=FALSE)
 base<-file.path(OUT,section,id)
 p<-p+plot_annotation(theme=theme(plot.margin=margin(2,2,2,2,'mm')))
 grDevices::cairo_pdf(paste0(base,'.pdf'),width=183/25.4,height=height/25.4,family=FONT);print(p);dev.off()
 svglite::svglite(paste0(base,'.svg'),width=183/25.4,height=height/25.4,system_fonts=list(sans=FONT));print(p);dev.off()
 ragg::agg_png(paste0(base,'.png'),width=183,height=height,units='mm',res=300,background='white');print(p);dev.off()
 saveRDS(p,paste0(base,'_plot.rds'))
 stopifnot(all(file.info(paste0(base,c('.pdf','.svg','.png')))$size>1000))
 cat('Rendered',id,'at 183 x',height,'mm\n')
 invisible(p)
}
# Independent viewports prevent unrelated panels from inheriting each other's axis widths.
makeContent.fitRaster <- function(x) {
 w<-convertWidth(unit(1,'npc'),'mm',valueOnly=TRUE);h<-convertHeight(unit(1,'npc'),'mm',valueOnly=TRUE)
 ww<-min(w,h*x$ratio);hh<-ww/x$ratio
 setChildren(x,gList(rasterGrob(x$im,width=unit(ww,'mm'),height=unit(hh,'mm'),interpolate=FALSE)))
}
raster_plot <- function(p,w,h) {
 f<-tempfile(fileext='.png');ragg::agg_png(f,width=w,height=h,units='mm',res=600,background='white');print(p);dev.off();im<-png::readPNG(f);unlink(f)
 wrap_elements(full=gTree(im=im,ratio=w/h,cl='fitRaster'))
}
colour_key <- function(palette,values,label,positions=seq(0,1,length.out=length(values))) {
 wrap_elements(full=grobTree(rasterGrob(matrix(colorRampPalette(palette,space='Lab')(512),nrow=1),x=.5,y=.78,width=unit(.7,'npc'),height=unit(2,'mm')),textGrob(values,x=.15+.7*positions,y=.4,gp=gpar(fontfamily=FONT,fontsize=6.8)),textGrob(label,x=.5,y=.10,just=c('center','bottom'),gp=gpar(fontfamily=FONT,fontsize=6.8))))
}
wrap_elements <- function(...) patchwork::wrap_elements(...)+theme(plot.margin=margin(0,0,0,0))
spatial_coordinates <- function(d) {
 z<-copy(d);z[,span:=max(diff(range(x)),diff(range(y))),by=sample_key]
 stopifnot(all(z$span>0));z[,`:=`(xp=(x-mean(range(x)))/span,yp=-(y-mean(range(y)))/span),by=sample_key];z
}
spatial_grid <- function(d,value,keys,palette,limits,title=NULL,ncol=3,trans='identity',midpoint=NULL) {
 z<-spatial_coordinates(d[sample_key%in%keys]);z[,sample:=factor(sample_key,levels=keys,labels=LABEL[keys])]
 q<-ggplot(z,aes(xp,yp,colour=.data[[value]]))+geom_point(size=.72,stroke=0)+facet_wrap(~sample,ncol=ncol)+coord_equal(xlim=c(-.53,.53),ylim=c(-.53,.53),expand=FALSE)+theme_void(base_family=FONT)+theme(strip.text=element_text(size=7.2,face='bold'),panel.spacing=unit(2,'mm'),plot.margin=margin(0,0,0,0),legend.position='right',legend.title=element_text(size=7),legend.text=element_text(size=7),legend.key.height=unit(3,'mm'))
 if(!is.null(midpoint))q<-q+scale_colour_gradient2(low=palette[1],mid=palette[2],high=palette[3],limits=limits,midpoint=midpoint,oob=scales::squish,name=title) else q<-q+scale_colour_gradientn(colours=palette,limits=limits,oob=scales::squish,trans=trans,name=title)
 q
}
image_file <- function(path) {
 im<-png::readPNG(path);wrap_elements(full=gTree(im=im,ratio=dim(im)[2]/dim(im)[1],cl='fitRaster'))
}
SHORT <- setNames(c('Sh1','Sh2','D3 r1','D3 r2','D3 r3','D7 r1','D7 r2','D7 r3','1 h','4 h','N3','N4'),SAMPLE)
MATRIX <- BASE+theme(axis.line=element_blank(),axis.ticks=element_blank())
panel_grob <- function(p) patchworkGrob(wrap_plots(list(p)) & theme(plot.margin=margin(0,0,0,0)))
stack_rows <- function(rows,heights) {
 stopifnot(length(rows)==length(heights))
 grobs<-lapply(seq_along(rows),function(i)grobTree(panel_grob(rows[[i]]),vp=viewport(layout.pos.row=i)))
 wrap_elements(full=grobTree(children=do.call(gList,grobs),vp=viewport(layout=grid.layout(length(rows),1,heights=unit(heights,'null')))))
}
side_by_side <- function(p,widths=rep(1,length(p))) wrap_plots(lapply(p,function(x)wrap_elements(full=panel_grob(x))),ncol=length(p),widths=widths)
