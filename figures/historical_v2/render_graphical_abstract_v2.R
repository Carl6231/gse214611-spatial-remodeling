#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])), 'figure_remaining_helpers.R'))
stopifnot(!any(file.exists(file.path(FIG_V2,c('Graphical_Abstract.pdf','Graphical_Abstract.svg')))))
score <- read_table(file.path(ROOT,'metadata','phase_10_spot_signature_scores.tsv.gz'))
keys <- SAMPLE[6:8]
limit <- max(abs(quantile(score[sample_key %in% keys & is.finite(signature_score),signature_score],c(.01,.99))))
stopifnot(is.finite(limit),limit>0)
assets <- lapply(keys,function(k){
  cached<-file.path(ASSET_V2,paste0('GA_',k,'.png'))
  if(file.exists(cached))return(cached)
  vals<-setNames(score[sample_key==k,signature_score],score[sample_key==k,barcode])
  spatial_values(file.path(ROOT,'objects','phase_09','niches',paste0(k,'_phase09_niche.rds')),
    vals,paste0('GA_',k),limits=c(-limit,limit),palette=c(COL_CONTRACTILE,COL_ZERO,COL_FIBROTIC),
    midpoint=0,width=56,height=38)
})
plate<-spatial_plate_v2(assets,1,3,unname(LABEL[keys]),header_height=.09)
gradient<-ggplot(data.frame(x=seq(-1,1,length.out=300),y=1),aes(x,y,fill=x))+
  geom_raster()+scale_fill_gradient2(low=COL_CONTRACTILE,mid=COL_ZERO,high=COL_FIBROTIC,
    midpoint=0,guide='none')+scale_x_continuous(NULL,breaks=NULL)+
  scale_y_continuous(NULL,breaks=NULL)+coord_cartesian(expand=FALSE)+quiet_theme+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=8.2,face='bold'),
    plot.margin=margin(0,1,0,1,unit='mm'))
gradient<-wrap_plots(list(gradient,wrap_plots(list(strip_text_v2('Contractile-rich',8.2,'bold',0),
  strip_text_v2('Spatial continuum',8.2,'bold',.5),strip_text_v2('Fibrotic-remodeling',8.2,'bold',1)),ncol=3)),
  ncol=1,heights=c(1,.45))
source_block<-function(title,genes,program,colour){
  wrap_elements(full=grobTree(
    textGrob(title,x=.5,y=.78,gp=gpar(fontfamily=FONT_V2,fontsize=8.2,fontface='bold',col=colour)),
    textGrob(genes,x=.5,y=.46,gp=gpar(fontfamily=FONT_V2,fontsize=7.2,col=COL_AXIS)),
    textGrob(program,x=.5,y=.16,gp=gpar(fontfamily=FONT_V2,fontsize=7.2,col=COL_AXIS))))
}
sources<-wrap_plots(list(
  source_block('Cardiomyocyte-associated','Fabp3 · Myl2 · Myh6','Contractile / oxidative programs',COL_CONTRACTILE),
  source_block('Fibroblast-associated','Serpinh1 · Thbs2 · Cthrc1','Extracellular-matrix remodeling',COL_FIBROTIC)),ncol=2)
bracket<-wrap_elements(full=grobTree(
  segmentsGrob(x0=.28,x1=.72,y0=.82,y1=.82,gp=gpar(col='#68717D',lwd=.6,lty=2)),
  segmentsGrob(x0=c(.28,.72),x1=c(.28,.72),y0=.82,y1=.65,gp=gpar(col='#68717D',lwd=.6)),
  textGrob('TREM2 / SPP1 macrophage co-enrichment',x=.5,y=.35,
    gp=gpar(fontfamily=FONT_V2,fontsize=7.5,col=COL_AXIS))))
boundary<-wrap_plots(list(note('Candidate secreted axes  ·  In silico counterfactual prediction',7.2,.5),
  note('CellOracle: 6/6 NEGATIVE_COMPLETE under the preregistered module rule',7.2,.5),
  note('Descriptive and computational evidence; causal testing remains required',6.5,.5)),ncol=1)
page<-wrap_plots(list(strip_text_v2('A reproducible Day 7 MI spatial remodeling continuum',10,'bold'),
  plate,gradient,sources,bracket,boundary),ncol=1,heights=c(7,42,10,17,9,17))
export_remaining(page,'Graphical_Abstract',height=110)
