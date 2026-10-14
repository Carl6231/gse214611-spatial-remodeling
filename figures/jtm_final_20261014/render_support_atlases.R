source('/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_structured/script/figure_remaining_helpers.R')

cohort_labels <- c('Formal comparison · Sham and Day 3',
                   'Formal comparison · Day 3 and Day 7',
                   'Descriptive and auxiliary tissues')
confidence_cols <- c(High='#36B6A4',Moderate='#9A78C6',Low='#F4777F',Unavailable='#D7DDE2')
broad_blocks <- niche_blocks <- vector('list',12)
niche_cols <- NICHE
names(niche_cols)[names(niche_cols)=='Not_macrophage_rich'] <- 'Not_macrophage_center'
niche_names <- c('Inflammatory','Mac–fibro remodeling','Vascular repair','Remote resident','Uncertain','Not selected')
names(niche_names) <- names(niche_cols)
saved <- read_table(file.path(ROOT,'metadata/phase_09_niche_grid_slice_summary.tsv'))[setting_id=='S054']
total_centers <- 0L
for(i in seq_along(SAMPLE)) {
  k <- SAMPLE[i]; message('Rendering ', k)
  o <- readRDS(file.path(ROOT,'objects/phase_09/deconvolution',paste0(k,'_phase09_expanded.rds')))
  stopifnot(all(c('P09E_RCTD_dominant','P09E_DWLS_dominant','P09E_confidence')%in%names(o@meta.data)))
  o$.mapconfidence <- as.character(o$P09E_confidence)
  if('P09E_expanded_inference_available'%in%names(o@meta.data))
    o$.mapconfidence[o$P09E_expanded_inference_available%in%FALSE] <- 'Unavailable'
  fields <- c('P09E_RCTD_dominant','P09E_DWLS_dominant','.mapconfidence')
  paths <- vapply(seq_along(fields),function(j)spatial_category(o,fields[j],
    if(j==3)confidence_cols else CELL,paste0('S6_complete_',j,'_',k),width=52,height=40),character(1))
  broad_blocks[[i]] <- block(wrap_plots(lapply(paths,img_panel_compact_v2),ncol=3),
    letters[i+10],unname(LABEL[k]))
  rm(o);gc(verbose=FALSE)
  o <- readRDS(file.path(ROOT,'objects/phase_09/niches',paste0(k,'_phase09_niche.rds')))
  z <- as.character(o$P09N_niche[o$P09N_center%in%TRUE]);total_centers<-total_centers+length(z)
  tab <- data.table(niche=names(niche_cols)[1:5],n=as.integer(table(factor(z,levels=names(niche_cols)[1:5]))))
  check <- merge(tab,saved[sample_key==k,.(niche,Freq)],by='niche',all=TRUE)
  check[is.na(Freq),Freq:=0L]
  stopifnot(all(check$n==check$Freq),sum(tab$n)==length(z))
  tab[,fraction:=n/sum(n)]
  native <- spatial_category(o,'P09N_niche',niche_cols,paste0('S8_corrected_',k),width=60,height=48)
  proportions <- ggplot(tab,aes(1,fraction,fill=factor(niche,levels=names(niche_cols))))+
    geom_col(width=.5,colour='white',linewidth=.15)+coord_flip()+
    scale_fill_manual(values=niche_cols,guide='none')+scale_y_continuous(limits=c(0,1),expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(0,0,0,0))
  detail <- note(paste0('Assigned ',sum(z!='Uncertain'),' / ',length(z),' centers; uncertain ',sum(z=='Uncertain')),6.6,.5)
  niche_blocks[[i]] <- block(wrap_plots(list(img_panel_compact_v2(native),proportions,detail),ncol=1,
    heights=unit.c(unit(1,'null'),unit(3,'mm'),unit(4,'mm'))),letters[i],unname(LABEL[k]))
  rm(o);gc(verbose=FALSE)
}
stopifnot(total_centers==5151L)
broad_pages <- niche_pages <- vector('list',3)
for(j in 1:3) {
  idx <- ((j-1)*4+1):(j*4)
  broad_pages[[j]] <- wrap_plots(list(
    strip_text_v2(paste('Supplementary Figure S6 | Complete corrected mapping',j,'/ 3'),8.2,'bold',0),
    note(cohort_labels[j],7.1,0),
    wrap_plots(lapply(c('RCTD dominant class','SpatialDWLS dominant class','Mapping confidence'),strip_text_v2),ncol=3),
    wrap_plots(broad_blocks[idx],ncol=1),
    legend_strip(c(CELL[names(CELL)!='Uncertain'],confidence_cols),
      c(unname(CLASS_LABEL[names(CELL)[names(CELL)!='Uncertain']]),names(confidence_cols)),rows=3),
    note('All retained spots remain visible. Broad classes are mixture estimates; Low confidence is not a separate biological class.',6.7)),
    ncol=1,heights=unit.c(unit(c(6,5,5),'mm'),unit(1,'null'),unit(c(15,6),'mm')))
  niche_pages[[j]] <- wrap_plots(list(
    strip_text_v2(paste('Supplementary Figure S8 | Corrected-reference niche atlas',j,'/ 3'),8.2,'bold',0),
    note(cohort_labels[j],7.1,0),wrap_plots(niche_blocks[idx],ncol=2),
    legend_strip(niche_cols,niche_names,rows=2),
    note('Fixed S054 rules; 5,151 centers across 12 tissues. Bars show within-tissue center fractions; gray spots were not selected.',6.7)),
    ncol=1,heights=unit.c(unit(c(6,5),'mm'),unit(1,'null'),unit(c(10,6),'mm')))
}
export_supplement(broad_pages,'Figure_S6_extended',245)
export_supplement(niche_pages,'Figure_S8',225)
