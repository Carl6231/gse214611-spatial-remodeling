#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
src<-file.path(ROOT,'results/figures/Figure_S4/source_data/v1r5')
registry<-read_table(file.path(src,'Figure_S4_panel_registry.tsv'))[match(SAMPLE,sample_key)]
edges<-read_table(file.path(src,'Figure_S4_delaunay_edges.tsv.gz'))
hmrf<-read_table(file.path(src,'Figure_S4_hmrf_assignments.tsv.gz'))
genes<-read_table(file.path(src,'Figure_S4_spatial_genes.tsv.gz'))
stopifnot(identical(registry$sample_key,SAMPLE),all(!edges$cross_slice_edge),
          all(!hmrf$numeric_ids_compared_across_slices),all(registry$top_svg_display_n==5))
base_palette<-c('#22A7C7','#B5A2D8','#35B89A','#D97FB2','#81CCEA','#5B9BD5')
s4<-s5<-vector('list',12)
for(i in seq_along(SAMPLE)){
  k<-SAMPLE[i];message('Saved Giotto geometry: ',k)
  obj<-readRDS(registry$full_object[i])
  # Frozen plain result components; no GiottoClass instantiation or refit.
  xy<-as.data.table(obj$source$coordinates)
  e<-edges[sample_id==k];h<-hmrf[sample_id==k]
  h<-merge(h,xy,by.x='spot',by.y='cell_ID',all.x=TRUE,sort=FALSE)
  stopifnot(nrow(h)==nrow(xy),all(is.finite(h$sdimx)),all(is.finite(h$sdimy)),
            all(e$from%in%xy$cell_ID),all(e$to%in%xy$cell_ID),
            uniqueN(h$selected_hmrf_column)==1L,
            identical(unique(h$selected_hmrf_column),registry$selected_hmrf_column[i]))
  pal<-base_palette[((seq_along(base_palette)+(i-1L)%%6L-1L)%%6L)+1L]
  if(i%%2L==0L)pal<-rev(pal)
  lev<-sort(unique(h$hmrf_domain_local));stopifnot(length(lev)<=6L)
  h[,domain:=factor(hmrf_domain_local,levels=lev)]
  network<-ggplot(e)+geom_segment(aes(sdimx_begin,sdimy_begin,xend=sdimx_end,yend=sdimy_end),
    colour='#5E7B94',alpha=.52,linewidth=.07)+
    geom_point(data=xy,aes(sdimx,sdimy),size=.08,colour='#5E7B94')+
    coord_fixed(expand=TRUE)+theme_void(base_family=FONT_V2)
  domain<-ggplot(h,aes(sdimx,sdimy,colour=domain))+geom_point(size=.38,alpha=.92)+
    scale_colour_manual(values=setNames(pal[seq_along(lev)],lev),guide='none')+
    coord_fixed(expand=TRUE)+theme_void(base_family=FONT_V2)
  local_key<-legend_strip(setNames(pal[seq_along(lev)],lev),paste0(gsub(' ','',short_label[k]),':',lev),rows=2)
  domain<-wrap_plots(list(domain,local_key),ncol=1,heights=unit.c(unit(1,'null'),unit(7,'mm')))
  s4[[i]]<-block(wrap_plots(list(heading(network,'Delaunay'),heading(domain,'Local HMRF')),ncol=2),letters[i],unname(LABEL[k]))
  top<-genes[sample_id==k][order(saved_rank)][1:5]
  top[,gene:=factor(feats,levels=rev(feats))]
  bar<-ggplot(top,aes(score,gene))+geom_segment(aes(x=0,xend=score,yend=gene),colour='#BBC5CE')+
    geom_point(colour='#536979',size=1.4)+labs(x='Saved SVG score',y=NULL)+quiet_theme+
    scale_x_continuous(breaks=scales::breaks_pretty(n=2),labels=scales::label_number(scale_cut=scales::cut_short_scale()),expand=expansion(mult=c(.05,.15)))+theme(plot.margin=margin(1,3,1,1,unit='mm'),axis.text.x=element_text(size=6.5))
  features<-head(obj$results$spatial_modules$spatial_cor$feat_order,5)
  cor<-as.data.table(obj$results$spatial_modules$spatial_cor$cor_DT)
  cor<-cor[feat_ID%in%features & as.character(variable)%in%features]
  stopifnot(length(features)==5L,nrow(cor)==25L,!anyDuplicated(cor[,.(feat_ID,variable)]))
  cor[,`:=`(a=factor(feat_ID,levels=features),b=factor(variable,levels=rev(features)))]
  heat<-ggplot(cor,aes(a,b,fill=spat_cor))+geom_tile(colour='white',linewidth=.12)+
    scale_fill_gradient2(low=COL_CONTRACTILE,mid='white',high=COL_FIBROTIC,midpoint=0,limits=c(-1,1),guide='none')+
    labs(x=NULL,y=NULL)+quiet_theme+coord_fixed()+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(angle=45,hjust=1))
  s5[[i]]<-block(wrap_plots(list(heading(bar,'Top five SVGs'),heading(heat,'Top-feature correlation')),ncol=2),letters[i],unname(LABEL[k]))
  rm(obj,xy,e,h,cor);gc(verbose=FALSE)
}
pages4<-atlas_pages(s4,'Supplementary Figure S4 | Slice-local spatial organization',
  'Frozen Giotto coordinates and edges; table-derived rendering. Local HMRF labels and colors are not aligned across tissues.',per_page=4)
pages5<-atlas_pages(s5,'Supplementary Figure S5 | Slice-local spatial features',
  'Saved top-five order. SVG scores are slice-specific; correlation uses a shared −1 to +1 color scale. No analysis was rerun.',
  legend=wrap_elements(full=grobTree(rasterGrob(matrix(colorRampPalette(c(COL_CONTRACTILE,'white',COL_FIBROTIC),space='Lab')(256),nrow=1),x=.5,y=.78,width=unit(60,'mm'),height=unit(1.6,'mm')),textGrob(c('-1','-0.5','0','0.5','1'),x=unit(.5,'npc')+unit(c(-30,-15,0,15,30),'mm'),y=.45,gp=gpar(fontfamily=FONT_V2,fontsize=7.1)),textGrob('Spatial correlation',x=.5,y=.08,gp=gpar(fontfamily=FONT_V2,fontsize=7.1)))))
if(wants('Figure_S4'))export_supplement(pages4,'Figure_S4',170)
if(wants('Figure_S5'))export_supplement(pages5,'Figure_S5',190)
