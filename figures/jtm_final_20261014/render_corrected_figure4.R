source('/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_polished/script/figure_common_v2.R')
OUT <- file.path(ROOT,'22_JTM_completion_20261014/figures_polished/corrected4')
dir.create(file.path(OUT,'figures'),recursive=TRUE,showWarnings=FALSE)
ASSET_V2 <- file.path(OUT,'tmp/f4_panels')
dir.create(ASSET_V2,recursive=TRUE,showWarnings=FALSE)
stopifnot('SpatialSpotPlot' %in% getNamespaceExports('scop'))
keys <- sample_order_v2
short <- setNames(c('Sh1','Sh2','D3r1','D3r2','D3r3','D7r1','D7r2','D7r3','1h','4h','N3','N4'),keys)
mi <- keys[3:8]
cols <- c('Not_macrophage_center'='#D0D0D0','inflammatory-infarct niche'='#EE6E79',
          'macrophage-fibroblast remodeling niche'='#A48ACB','vascular-repair niche'='#39B7A5',
          'remote resident niche'='#428BD0','Uncertain'='#6F6F6F')
labs_niche <- c('Inflammatory','Mac-fibro remodeling','Vascular repair','Remote resident','Uncertain')
niches <- names(cols)[-1]
shapes <- setNames(c(1,21,22,24,23,25),names(cols))
grid_summary <- fread(file.path(ROOT,'metadata/phase_09_niche_parameter_summary.tsv'))[setting_id=='S054']
stopifnot(nrow(grid_summary)==1L,grid_summary$centers==5151L,
          grid_summary$center_threshold==.20,grid_summary$radius_multiplier==1.2,
          grid_summary$program_gate==.10,grid_summary$margin_gate==.15)
prev <- fread(file.path(ROOT,'metadata/phase_09_niche_grid_slice_summary.tsv'))[setting_id=='S054']
coh <- fread(file.path(ROOT,'metadata/phase_09_niche_spatial_null.tsv'))
axis <- fread(file.path(ROOT,'metadata/phase_09_niche_continuous_slice_summary.tsv'))
loo <- fread(file.path(ROOT,'metadata/phase_09_niche_leave_one_slice_out.tsv'))
md <- list(); maps <- list()
for (key in keys) {
  o <- readRDS(file.path(ROOT,'objects/phase_09/niches',paste0(key,'_phase09_niche.rds')))
  fields <- c('P09N_center','P09N_niche','P09N_niche_margin','P09N_inflammatory_axis',
              'P09N_remodeling_axis','P09N_vascular_axis','P09N_remote_axis')
  stopifnot(all(fields %in% names(o@meta.data)),length(Images(o))==1L)
  z <- as.data.table(o@meta.data[,fields,drop=FALSE]); z[,sample_key:=key]; md[[key]]<-z
  stopifnot(sum(z$P09N_center)==unique(prev[sample_key==key,centers]))
  if (key %in% mi) {
    p <- suppressMessages(scop::SpatialSpotPlot(o,group.by='P09N_niche',image=Images(o)[[1]],
         assay='Spatial',layer='counts',plot_type='point',crop=TRUE,overlay_image=TRUE,
         image.alpha=.18,pt.alpha=.92,pt.size=.44,stroke=0,legend.position='none',
         theme_use='theme_blank',verbose=FALSE))
    stopifnot(inherits(p,'ggplot'),all(c('x','y','.value') %in% names(p$data)))
    p <- suppressMessages(p+scale_fill_manual(values=cols,limits=names(cols),drop=FALSE))
    stopifnot(identical(unname(ggplot_build(p)$plot$scales$get_scales('fill')$map(names(cols))),unname(cols)))
    pp <- p$data[p$data$.value!='Not_macrophage_center',,drop=FALSE]
    p <- p+geom_point(data=pp,aes(x,y,shape=.value),inherit.aes=FALSE,fill=NA,
                     colour=COL_AXIS,size=.36,stroke=.055,alpha=.82,show.legend=FALSE)+
      scale_shape_manual(values=shapes,limits=names(cols),drop=FALSE)
    p <- suppressMessages(p+coord_fixed())+theme_void(base_family=FONT_V2)+
      theme(legend.position='none',plot.margin=margin(0,0,0,0))
    maps[[key]] <- save_spatial_asset_v2(p,paste0('corrected_',key),53,36)
  }
  rm(o); gc(verbose=FALSE)
}
md <- rbindlist(md); stopifnot(sum(md$P09N_center)==5151L)
actual <- md[P09N_center==TRUE,.N,by=.(sample_key,niche=P09N_niche)]
compare <- merge(actual,prev[,.(sample_key,niche,Freq)],all=TRUE)
stopifnot(all(compare$N==compare$Freq),!anyNA(compare))
cat('SIX_CORRECTED_NATIVE_PANELS_READY\n')

theme_set(theme_editorial_v2())
role <- function(x) ifelse(grepl('^Sham',x),'Sham',ifelse(grepl('Day3',x)&grepl('^MI',x),'MI Day 3',
  ifelse(grepl('Day7',x),'MI Day 7',ifelse(grepl('Needle',x),'Needle-pass','Descriptive'))))
counts <- md[P09N_center==TRUE,.(Assigned=sum(P09N_niche!='Uncertain'),Uncertain=sum(P09N_niche=='Uncertain')),by=sample_key]
counts <- melt(counts,id.vars='sample_key',variable.name='assignment',value.name='centers')
counts[,sample_key:=factor(sample_key,levels=keys)]
p_b <- ggplot(counts,aes(sample_key,centers,fill=assignment))+geom_col(width=.75)+
  scale_fill_manual(values=c(Assigned='#36B6A4',Uncertain='#AEB7C2'),name=NULL)+
  scale_x_discrete(labels=short)+labs(x=NULL,y='Centers')+theme(legend.position='right',axis.text.x=element_text(angle=45,hjust=1))
pm <- merge(CJ(sample_key=keys,niche=niches),prev[,.(sample_key,niche,fraction)],all.x=TRUE)
pm[is.na(fraction),fraction:=0]; pm[,sample_key:=factor(sample_key,levels=keys)]
pm[,niche:=factor(niche,levels=rev(niches),labels=rev(labs_niche))]
p_d <- ggplot(pm,aes(sample_key,niche,fill=fraction))+geom_tile(colour='white',linewidth=.25)+
  scale_fill_gradient(low='#F5F7F8',high='#334E68',limits=c(0,1),name='Fraction')+
  scale_x_discrete(labels=short)+labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank())
axis_fields <- c('P09N_inflammatory_axis','P09N_remodeling_axis','P09N_vascular_axis','P09N_remote_axis')
ax <- melt(md[P09N_center==TRUE],id.vars=c('sample_key','P09N_niche'),measure.vars=axis_fields,variable.name='axis',value.name='value')
ax <- ax[is.finite(value),.(value=median(value)),by=.(sample_key,P09N_niche,axis)][,.(value=median(value)),by=.(P09N_niche,axis)]
ax[,`:=`(axis=factor(axis,levels=axis_fields,labels=c('Inflam.','Remodel.','Vascular','Remote')),
          niche=factor(P09N_niche,levels=rev(niches),labels=rev(labs_niche)))]
lim<-max(abs(ax$value)); p_e<-ggplot(ax,aes(axis,niche,fill=value))+geom_tile(colour='white',linewidth=.25)+
  scale_fill_gradient2(low='#397FC2',mid='#FAFCFD',high='#ED6A75',limits=c(-lim,lim),name='Axis score')+
  labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=35,hjust=1))
coh[,`:=`(sample_key=factor(sample_key,levels=keys),role=role(sample_key))]
p_f<-ggplot(coh[is.finite(enrichment)],aes(sample_key,enrichment,colour=role))+
  geom_hline(yintercept=0,colour='#A2A8AE',linewidth=.3)+geom_segment(aes(xend=sample_key,y=0,yend=enrichment),linewidth=.5)+geom_point(size=1.6)+
  geom_text(data=coh[!is.finite(enrichment)],aes(sample_key,y=0,label='NE'),inherit.aes=FALSE,size=7.1/GG_PT,family=FONT_V2)+
  scale_colour_manual(values=role_cols_v2,guide='none')+scale_x_discrete(labels=short,drop=FALSE)+
  labs(x=NULL,y='Observed - permuted')+theme(axis.text.x=element_text(angle=45,hjust=1))
axis[,`:=`(sample_key=factor(sample_key,levels=keys),role=role(sample_key))]
p_g<-ggplot(axis,aes(sample_key,remodeling_axis,colour=role))+geom_hline(yintercept=0,colour='#A2A8AE',linewidth=.3)+geom_point(size=1.8)+
  scale_colour_manual(values=role_cols_v2,guide='none')+scale_x_discrete(labels=short)+
  labs(x=NULL,y='Remodeling axis')+theme(axis.text.x=element_text(angle=45,hjust=1))
# The stored contrast and leave-one-tissue-out values are displayed, never refitted.
lh <- loo[metric=='remodeling_axis' & contrast %in% c('Day7_MI_vs_Sham','Day7_MI_vs_Day3_MI')]
stopifnot(nrow(lh)==2L)
lh[,contrast:=factor(contrast,levels=c('Day7_MI_vs_Day3_MI','Day7_MI_vs_Sham'),labels=c('D7 vs D3','D7 vs Sham'))]
p_h <- ggplot(lh,aes(full_effect,contrast))+geom_segment(aes(x=min_leave_one_out_effect,xend=max_leave_one_out_effect,yend=contrast),linewidth=.6,colour='#9A78C6')+
  geom_point(size=1.8,colour='#9A78C6')+labs(x='Full effect and leave-one-out range',y=NULL)+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())

draw_plot <- function(p,x,y,w,h,letter,title,subtitle=NULL) {
  pushViewport(viewport(x=unit(x,'mm'),y=unit(225-y,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','top')))
  grid.text(letter,x=unit(0,'mm'),y=unit(1,'npc')-unit(1,'mm'),just=c('left','top'),gp=gpar(fontfamily=FONT_V2,fontsize=10,fontface='bold'))
  grid.text(title,x=unit(4,'mm'),y=unit(1,'npc')-unit(1,'mm'),just=c('left','top'),gp=gpar(fontfamily=FONT_V2,fontsize=7.4,fontface='bold'))
  hh<-if(is.null(subtitle)) 6 else 10
  if(!is.null(subtitle)) grid.text(subtitle,x=unit(4,'mm'),y=unit(1,'npc')-unit(5,'mm'),just=c('left','top'),gp=gpar(fontfamily=FONT_V2,fontsize=7.1,col='#555D66'))
  pushViewport(viewport(x=.5,y=0,width=1,height=unit(h-hh,'mm'),just=c('center','bottom')))
  grid.draw(if(inherits(p,'patchwork')) patchworkGrob(p) else ggplotGrob(p)); upViewport(2)
}
draw_page <- function() {
  grid.newpage()
  p_a<-ggplot()+annotate('text',x=0,y=1,hjust=0,vjust=1,label='11 corrected broad classes\n18,849 balanced reference nuclei\nS054: myeloid >=0.20; radius x1.2\nProgram >=0.10; margin >=0.15',size=7.1/GG_PT,family=FONT_V2,lineheight=1.25)+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+theme_void()
  draw_plot(p_a,2,2,60,30,'a','Reference and fixed settings')
  draw_plot(p_b,66,2,115,30,'b','5,151 corrected centers')
  grid.text('c',x=unit(2,'mm'),y=unit(190,'mm'),just='left',gp=gpar(fontfamily=FONT_V2,fontsize=10,fontface='bold'))
  grid.text('Spatial niches · corrected reference',x=unit(6,'mm'),y=unit(190,'mm'),just='left',gp=gpar(fontfamily=FONT_V2,fontsize=7.4,fontface='bold'))
  for(i in 1:6) {
    x<-17+((i-1)%%3)*54; y<-43+((i-1)%/%3)*36
    if(i<=3) grid.text(paste0('r',i),x=unit(x+26.5,'mm'),y=unit(225-y+2,'mm'),gp=gpar(fontfamily=FONT_V2,fontsize=7.2,fontface='bold'))
    im<-png::readPNG(maps[[mi[i]]]); grid.raster(im,x=unit(x+26.5,'mm'),y=unit(225-y-17,'mm'),width=unit(53,'mm'),height=unit(36,'mm'),interpolate=FALSE)
  }
  grid.text('Day 3',x=unit(9,'mm'),y=unit(165,'mm'),rot=90,gp=gpar(fontfamily=FONT_V2,fontsize=7.2,fontface='bold'))
  grid.text('Day 7',x=unit(9,'mm'),y=unit(129,'mm'),rot=90,gp=gpar(fontfamily=FONT_V2,fontsize=7.2,fontface='bold'))
  for(i in 1:6) {x<-8+((i-1)%%3)*59; y<-118+((i-1)%/%3)*4
    grid.points(x=unit(x,'mm'),y=unit(225-y,'mm'),pch=shapes[i],size=unit(1.7,'mm'),gp=gpar(col=COL_AXIS,fill=cols[i]))
    grid.text(c('Not a selected center',labs_niche)[i],x=unit(x+3,'mm'),y=unit(225-y,'mm'),just='left',gp=gpar(fontfamily=FONT_V2,fontsize=7.1))}
  draw_plot(p_d,2,125,179,29,'d','Niche prevalence')
  draw_plot(p_e,2,156,88,36,'e','Axis structure by niche','Median of tissue-level niche medians')
  draw_plot(p_f,94,156,87,36,'f','Within-slice coherence','NE: no center-to-center edges')
  draw_plot(p_g,2,194,88,29,'g','Remodeling axis')
  draw_plot(p_h,94,194,87,29,'h','Leave-one-tissue-out stability')
}
pdf<-file.path(OUT,'figures/Figure_4_corrected.pdf'); svg<-sub('pdf$','svg',pdf); png<-sub('pdf$','png',pdf)
stopifnot(!any(file.exists(c(pdf,svg,png))))
cairo_pdf(pdf,width=183/25.4,height=225/25.4,family='sans'); draw_page(); dev.off()
svglite::svglite(svg,width=183/25.4,height=225/25.4); draw_page(); dev.off()
rsvg::rsvg_png(svg,png,width=2161,height=2657)
fwrite(actual,file.path(OUT,'tmp/corrected_niche_counts.tsv'),sep='\t')
cat('CORRECTED_FIGURE4_RENDERED\n')
