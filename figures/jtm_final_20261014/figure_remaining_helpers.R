source('/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_final/script/figure_common_v2.R')

STAGE_V2 <- Sys.getenv('REMAINING_STAGE', file.path(V2, 'tmp', 'remaining_build'))
stopifnot(startsWith(normalizePath(STAGE_V2, mustWork=FALSE), paste0(V2, '/tmp/')))
ASSET_V2 <- file.path(STAGE_V2, 'assets')
dir.create(ASSET_V2, recursive=TRUE, showWarnings=FALSE)
requested_figures <- strsplit(Sys.getenv('TARGET_FIGURES', ''), ',', fixed=TRUE)[[1]]
wants <- function(id) !length(requested_figures) || !nzchar(requested_figures[1]) || id %in% requested_figures
SAMPLE <- sample_order_v2
LABEL <- sample_labels_v2
FORMAL <- SAMPLE[1:8]
CELL <- c(Cardiomyocyte='#5B9BD5', Endothelial='#22A7C7', Fibroblast='#35B89A',
          Macrophage='#D97FB2', Monocyte='#B57AC3', Neutrophil='#F58F8F',
          Dendritic='#7B6BC2', B_cell='#60C3A8', T_NK='#9BCFBB',
          Pericyte='#B5A2D8', Schwann_cell='#81CCEA', Uncertain='#AEB7C2')
STATE <- c(Resident_Folr2_macrophage='#428BD0', Trem2_Spp1_macrophage='#EE6E79',
           MHCII_macrophage='#39B7A5', Cycling_macrophage='#B5A2D8',
           Stress_response_macrophage='#9F79BE', Uncertain='#AEB7C2')
NICHE <- c('inflammatory-infarct niche'='#EE6E79',
           'macrophage-fibroblast remodeling niche'='#A48ACB',
           'vascular-repair niche'='#39B7A5', 'remote resident niche'='#428BD0',
           Uncertain='#AEB7C2', Not_macrophage_rich='#D8D8D8')
CLASS_LABEL <- c(Cardiomyocyte='Cardiomyocyte', Endothelial='Endothelial',
  Fibroblast='Fibroblast', Macrophage='Macrophage', Monocyte='Monocyte',
  Neutrophil='Neutrophil', Dendritic='Dendritic', B_cell='B cell', T_NK='T / NK',
  Pericyte='Pericyte', Schwann_cell='Schwann cell', Uncertain='Uncertain')
STATE_LABEL <- c(Resident_Folr2_macrophage='Resident Folr2',
  Trem2_Spp1_macrophage='Trem2 / Spp1', MHCII_macrophage='MHCII',
  Cycling_macrophage='Cycling', Stress_response_macrophage='Stress response',
  Uncertain='Uncertain')
TIME <- c('Day0','1HR','4HR','Day1','Day3','Day7')
TIME_LABEL <- setNames(c('D0','1 h','4 h','D1','D3','D7'), TIME)
role <- function(x) ifelse(grepl('^Sham',x),'Sham',
  ifelse(grepl('^MI_Day3',x),'MI Day 3',ifelse(grepl('^MI_Day7',x),'MI Day 7',
  ifelse(grepl('^Needle',x),'Needle-pass','Descriptive'))))
sample_colour <- function(x) unname(role_cols_v2[role(as.character(x))])
short_label <- setNames(c('Sh1','Sh2','D3 r1','D3 r2','D3 r3','D7 r1','D7 r2','D7 r3',
                         '1 h','4 h','N3','N4'), SAMPLE)
read_table <- read_tsv_v2
quiet_theme <- theme_editorial_v2()
block <- function(p,l,t,sub=NULL) evidence_header_v2(p,l,t,sub,header_mm=if(is.null(sub)) 5 else 9)
note <- function(text, size=7.1, align=0) strip_text_v2(text,size,'plain',align)
heading <- function(p,title) wrap_plots(list(strip_text_v2(title),p),ncol=1,
  heights=unit.c(unit(if(grepl('\n',title,fixed=TRUE))7 else 4,'mm'),unit(1,'null')))
legend_strip <- function(cols, labels=names(cols), rows=2) {
  d <- data.frame(id=factor(names(cols),levels=names(cols)),x=seq_along(cols),y=1)
  p <- ggplot(d,aes(x,y,colour=id))+geom_point(alpha=0)+
    scale_colour_manual(values=cols,labels=labels,drop=FALSE,name=NULL)+
    guides(colour=guide_legend(nrow=rows,byrow=TRUE,override.aes=list(alpha=1,size=2)))+
    theme_void(base_family=FONT_V2)+theme(legend.position='bottom',
      legend.text=element_text(size=7.1),legend.key.width=unit(3,'mm'),
      legend.key.height=unit(2.5,'mm'),legend.margin=margin(0,0,0,0))
  p
}
role_strip <- function(keys=SAMPLE) {
  d <- data.frame(sample=factor(keys,levels=keys),role=role(keys))
  ggplot(d,aes(sample,1,fill=role))+geom_tile()+
    scale_fill_manual(values=role_cols_v2,guide='none')+theme_void()+
    scale_x_discrete(expand=c(0,0))+scale_y_continuous(expand=c(0,0))
}
with_role <- function(p,keys=SAMPLE) wrap_plots(list(role_strip(keys),p),ncol=1,
  heights=unit.c(unit(1.5,'mm'),unit(1,'null')))

export_remaining <- function(p,id,section='main',height=225) {
  dest <- file.path(STAGE_V2,'figures',section)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(STAGE_V2,'qa'),recursive=TRUE,showWarnings=FALSE)
  files <- c(pdf=file.path(dest,paste0(id,'.pdf')),svg=file.path(dest,paste0(id,'.svg')),
             png=file.path(STAGE_V2,'qa',paste0(id,'_render.png')))
  stopifnot(!any(file.exists(files)))
  p <- p+plot_annotation(theme=theme(plot.margin=margin(2,2,2,2,unit='mm')))
  cairo_pdf(files['pdf'],width=183/25.4,height=height/25.4,family='sans',onefile=TRUE)
  print(p); dev.off()
  svglite::svglite(files['svg'],width=183/25.4,height=height/25.4,
    bg='white',system_fonts=list(sans=FONT_V2))
  print(p); dev.off()
  rsvg::rsvg_png(files['svg'],files['png'],width=2161,height=round(2161*height/183))
  stopifnot(all(file.info(files)$size>1000))
  cat('EXPORTED',id,paste(file.info(files)$size,collapse=','),'\n')
  invisible(files)
}
export_supplement <- function(pages,id,height=225) {
  dest<-file.path(STAGE_V2,'figures','supplementary')
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(STAGE_V2,'qa'),recursive=TRUE,showWarnings=FALSE)
  pdf<-file.path(dest,paste0(id,'.pdf'))
  stopifnot(!file.exists(pdf))
  pages<-lapply(pages,function(p)p+plot_annotation(theme=theme(plot.margin=margin(2,2,2,2,unit='mm'))))
  cairo_pdf(pdf,width=183/25.4,height=height/25.4,family='sans',onefile=TRUE)
  for(p in pages)print(p)
  dev.off()
  for(i in seq_along(pages)){
    suffix<-if(length(pages)==1L)'' else paste0('_p',i)
    svg<-file.path(dest,paste0(id,suffix,'.svg'))
    png<-file.path(STAGE_V2,'qa',paste0(id,suffix,'_render.png'))
    stopifnot(!any(file.exists(c(svg,png))))
    svglite::svglite(svg,width=183/25.4,height=height/25.4,bg='white',system_fonts=list(sans=FONT_V2))
    print(pages[[i]]);dev.off()
    rsvg::rsvg_png(svg,png,width=2161,height=round(2161*height/183))
    stopifnot(all(file.info(c(svg,png))$size>1000))
  }
  stopifnot(file.info(pdf)$size>1000)
  cat('EXPORTED',id,'pages=',length(pages),'\n')
}

spatial_leaf <- function(p) {
  if(inherits(p,'patchwork')) p <- p[[1]]
  stopifnot(inherits(p,'ggplot'),all(c('x','y') %in% names(p$data)))
  xr<-range(p$data$x,finite=TRUE);yr<-range(p$data$y,finite=TRUE)
  p<-suppressMessages(p+coord_fixed(xlim=xr+c(-.05,.05)*diff(xr),
    ylim=yr+c(-.05,.05)*diff(yr),expand=FALSE))
  p+labs(title=NULL,subtitle=NULL)+theme_void(base_family=FONT_V2)+
    theme(legend.position='none',plot.title=element_blank(),plot.subtitle=element_blank(),
          strip.text=element_blank(),plot.margin=margin(0,0,0,0))
}
spatial_category <- function(path,field,cols,id,cells=NULL,width=54,height=42) {
  cached <- file.path(ASSET_V2,paste0(id,'.png'))
  if(file.exists(cached)) return(cached)
  obj <- if(inherits(path,'Seurat')) path else readRDS(path)
  stopifnot(field %in% colnames(obj[[]]))
  vals <- as.character(obj[[field]][,1]); names(vals) <- colnames(obj)
  if(is.null(cells)) cells <- names(vals) else {
    vals[!names(vals) %in% cells] <- 'Background'
    cols <- c(cols,Background='#E4E4E4')
    cells <- names(vals)
  }
  actual <- intersect(names(cols),unique(vals[cells]))
  stopifnot(length(actual)>0,all(na.omit(vals[cells]) %in% names(cols)))
  obj[[field]] <- factor(vals,levels=names(cols))
  p <- suppressMessages(scop::SpatialSpotPlot(obj,group.by=field,cells=cells,
    image=Images(obj)[[1]],plot_type='point',crop=TRUE,overlay_image=TRUE,
    image.alpha=.18,pt.alpha=.92,pt.size=.46,stroke=0,palcolor=unname(cols[actual]),
    legend.position='none',theme_use='theme_blank'))
  p <- spatial_leaf(p)
  aesthetic <- unique(unlist(lapply(p$layers,function(l)intersect(names(l$mapping),c('colour','fill')))))
  stopifnot(length(aesthetic)==1L)
  p <- suppressMessages(p+scale_colour_manual(values=cols,aesthetics=aesthetic,drop=FALSE,na.value='#E4E4E4'))
  s <- ggplot_build(p)$plot$scales$get_scales(aesthetic)
  stopifnot(identical(unname(s$map(actual)),unname(cols[actual])))
  out <- save_spatial_asset_v2(p,id,width,height)
  rm(obj,p);gc(verbose=FALSE);out
}
atlas_pages <- function(blocks,title,footer,ncol=2,per_page=6,legend=NULL) {
  n<-ceiling(length(blocks)/per_page)
  lapply(seq_len(n),function(i){
    idx<-seq((i-1)*per_page+1,min(i*per_page,length(blocks)))
    content<-blocks[idx]
    if(length(content)<per_page)content<-c(content,rep(list(plot_spacer()),per_page-length(content)))
    parts<-list(strip_text_v2(paste0(title,'  ·  ',i,'/',n),8.2,'bold',0),wrap_plots(content,ncol=ncol))
    heights<-unit.c(unit(6,'mm'),unit(1,'null'))
    if(!is.null(legend)){parts<-c(parts,list(legend));heights<-unit.c(heights,unit(10,'mm'))}
    parts<-c(parts,list(note(footer)));heights<-unit.c(heights,unit(5,'mm'))
    wrap_plots(parts,ncol=1,heights=heights)
  })
}
spatial_values <- function(path,values,id,limits,palette,width=54,height=42,midpoint=NULL) {
  obj<-readRDS(path);values<-values[colnames(obj)]
  stopifnot(length(values)==ncol(obj),all(is.finite(values)|is.na(values)))
  p<-suppressMessages(scop::SpatialSpotPlot(obj,values=values,image=Images(obj)[[1]],
    plot_type='point',crop=TRUE,overlay_image=TRUE,image.alpha=.18,pt.alpha=.92,
    pt.size=.46,stroke=0,palcolor=palette,legend.position='none',theme_use='theme_blank'))
  p<-spatial_leaf(p)
  aesthetic<-unique(unlist(lapply(p$layers,function(l)intersect(names(l$mapping),c('colour','fill')))))
  stopifnot(length(aesthetic)==1L)
  if(!is.null(midpoint))p<-suppressMessages(p+scale_colour_gradient2(low=palette[1],mid=palette[2],high=palette[3],
    midpoint=midpoint,limits=limits,oob=scales::squish,aesthetics=aesthetic)) else
    p<-suppressMessages(p+scale_colour_gradientn(colours=palette,limits=limits,oob=scales::squish,aesthetics=aesthetic))
  stopifnot(identical(as.numeric(ggplot_build(p)$plot$scales$get_scales(aesthetic)$limits),as.numeric(limits)))
  out<-save_spatial_asset_v2(p,id,width,height)
  rm(obj,p);gc(verbose=FALSE);out
}

# Draw frozen deconvolution proportions against the object's own image coordinates.
# SCOP 0.9.1 rejects older result contracts; no model is refitted for plotting.
frozen_deconvolution_plot <- function(srt,tool_name,cell_types,combine=FALSE,...) {
 mat<-srt@tools[[tool_name]]$proportions
 stopifnot(all(colnames(srt)%in%rownames(mat)),cell_types%in%colnames(mat))
 vals<-mat[colnames(srt),cell_types];names(vals)<-colnames(srt)
 field<-paste0(tool_name,'_prop_',cell_types)
 if(field%in%names(srt@meta.data))stopifnot(isTRUE(all.equal(as.numeric(vals),as.numeric(srt@meta.data[[field]]))))
 p<-scop::SpatialSpotPlot(srt,values=vals,...)
 list(p)
}
