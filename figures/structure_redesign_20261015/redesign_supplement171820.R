source('11_图件结构重设计_20261015/scripts/design_common.R')
ragg::agg_capture(width=183,height=220,units='mm',res=300)
SEC<-file.path(ROOT,'07_分析补全与组合图重排_20261014/Additional_files/SecAct_tables')
export_pages<-function(pages,id,heights){
 for(i in seq_along(pages))export(pages[[i]],paste0(id,if(length(pages)>1)paste0('_p',i)else''),heights[i],section='Supplementary_figures')
 pdfs<-file.path(OUT,'Supplementary_figures',paste0(id,if(length(pages)>1)paste0('_p',seq_along(pages))else'', '.pdf'))
 if(length(pages)>1)qpdf::pdf_combine(pdfs,output=file.path(OUT,'Supplementary_figures',paste0(id,'.pdf')))
}
# S17: each row is its own fitted model, with no empty factor slots.
x<-fread(file.path(SEC,'secact_pattern_spot_loadings.tsv.gz'));ids<-c(paste0('MI_D3_rep',1:3),paste0('MI_D7_rep',1:3))
seqpal<-c('#EEF8F6','#97DCCA','#34B8A4','#0B7F88')
pages<-heights<-list()
for(j in 1:2){rows<-list();rh<-numeric()
 for(k in 1:3){id<-ids[(j-1)*3+k];d<-x[sample_id==id];ranks<-sort(unique(d$pattern_index));stopifnot(identical(as.integer(ranks),seq_along(ranks)))
  maps<-lapply(ranks,function(r){z<-d[pattern_index==r];z[,relative:=if(max(loading)>0)loading/max(loading)else 0];p<-ggplot(z,aes(x_um,y_um,colour=relative))+geom_point(size=.47)+coord_equal()+scale_colour_gradientn(colours=seqpal,limits=c(0,1),guide='none')+theme_void();heading(raster_plot(p,38,38),paste0('r',k,' · F',r))})
  width<-min(44,176/length(ranks));row<-side_by_side(maps)
  rows[[k]]<-block(row,letters[k],paste0('MI Day ',if(j==1)3 else 7,' · tissue ',k,' | ',length(ranks),' independently fitted factors'))
  rh[k]<-width+8
 }
 pages[[j]]<-stack_rows(c(list(strip(paste0('Supplementary Figure S17 | ',if(j==1)'Day 3'else'Day 7',' tissue-local factors'),9,'bold')),rows,list(colour_key(seqpal,c('0','0.5','1'),'Loading / factor-specific maximum'),strip('Factor identities are local to each tissue; columns are not matched biological programs across tissues.',6.8))),c(8,rh,10,7))
 heights[[j]]<-sum(c(8,rh,10,7))+4
}
export_pages(pages,'Figure_S17',unlist(heights))
# S18: sending and receiving are adjacent within every tissue group.
vs<-fread(file.path(SEC,'secact_vector_fields_native_display.tsv.gz'));bg<-fread(file.path(SEC,'secact_vector_background_native_display.tsv.gz'))
lig<-c('C3','CCL8','CXCL12','GAS6','MDK','PTN');stopifnot(uniqueN(bg[,.(sample_id,ligand,signal_mode)])==72L)
pages<-lapply(lig,function(g){blocks<-lapply(seq_along(ids),function(i){id<-ids[i];pair<-lapply(c('sending','receiving'),function(mode){z<-bg[sample_id==id&ligand==g&signal_mode==mode];s<-vs[sample_id==id&ligand==g&signal_mode==mode];stopifnot(nrow(z)>0)
 q<-ggplot(z,aes(x,y))+geom_point(colour='#DDE6EA',size=.18)+geom_segment(data=s,aes(x=x_start,y=y_start,xend=x_end,yend=y_end),inherit.aes=FALSE,colour='#088FA4',linewidth=.15,arrow=arrow(length=unit(.5,'mm')))+coord_flip()+theme_void()+theme(aspect.ratio=diff(range(c(z$x,s$x_start,s$x_end)))/diff(range(c(z$y,s$y_start,s$y_end))))
 heading(raster_plot(q,41,35),if(mode=='sending')'Sending'else'Receiving')})
 block(side_by_side(pair),letters[i],paste0('MI ',gsub('_',' ',sub('MI_','',id))))})
 stack_rows(list(strip(paste0('Supplementary Figure S18 | ',g,' · paired vector fields'),9,'bold'),side_by_side(blocks[1:2]),side_by_side(blocks[3:4]),side_by_side(blocks[5:6]),strip('Native display vectors are retained. Arrow lengths were rescaled by the backend and are not activity or movement magnitudes.',6.6)),c(8,53,53,53,8))})
export_pages(pages,'Figure_S18',rep(181,6))
# S19/S20: all completed tumor-reference projection coefficients remain descriptive.
eco<-file.path(ROOT,'08_SpatialEcoTyper真实运行_20261015');components<-c('NonSE',paste0('SE',1:9))
z<-rbindlist(lapply(SAMPLE,function(k){o<-readRDS(file.path(eco,'objects',paste0(k,'_SpatialEcoTyper.rds')));stopifnot(nrow(o$abundance)==nrow(o$coordinates),max(abs(rowSums(o$abundance)-1))<1e-8)
 d<-as.data.table(o$coordinates);d[,sample_key:=k];d[,largest:=factor(components[max.col(o$abundance,ties.method='first')],levels=components)];d[,maximum:=apply(o$abundance,1,max)];cbind(d,as.data.table(o$abundance))}))
z<-spatial_coordinates(z);stopifnot(nrow(z)==25548L)
summary<-fread(file.path(eco,'tables/tissue_summary.tsv'));summary[,sample:=factor(sample,levels=SAMPLE)]
EPAL<-setNames(c('#ADBCC5','#168FAE','#51BACA','#1CA68C','#8BCC9B','#4D89D4','#8A73C8','#B498DA','#EF7E93','#F4AC9C'),components)
sm<-copy(summary);sm[,sample:=factor(sample,levels=rev(SAMPLE))]
h<-ggplot(sm,aes(component,sample,fill=mean_abundance))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=sprintf('%.2f',mean_abundance)),size=2.2)+scale_fill_gradientn(colours=seqpal,limits=c(0,1),trans='sqrt',breaks=c(0,.1,.3,1),name='Mean coefficient\n(sqrt scale)')+scale_y_discrete(NULL,labels=LABEL)+labs(x=NULL)+MATRIX+theme(axis.text.x=element_text(size=6.5),legend.position='bottom',legend.key.width=unit(8,'mm'),legend.key.height=unit(2,'mm'))
a<-block(h,'a','Mean reference coefficients')
mm<-z[,.(median=median(maximum),q25=quantile(maximum,.25),q75=quantile(maximum,.75)),by=sample_key];mm[,sample:=factor(sample_key,levels=rev(SAMPLE))]
p<-ggplot(mm,aes(median,sample,colour=role(sample_key)))+geom_linerange(aes(xmin=q25,xmax=q75),linewidth=.6)+geom_point(size=1.65)+scale_colour_manual(values=COL,guide='none')+scale_x_continuous('Largest coefficient · median (IQR)',limits=c(0,1),breaks=c(0,.5,1))+scale_y_discrete(NULL,labels=LABEL)+BASE+theme(axis.title.x=element_text(size=6.7),axis.text.y=element_text(size=6.7))
b<-block(p,'b','Projection concentration','Within-tissue spot distribution; not confidence')
z[,sample:=factor(sample_key,levels=SAMPLE,labels=LABEL[SAMPLE])]
p<-ggplot(z,aes(xp,yp,colour=largest))+geom_point(size=.47,stroke=0)+facet_wrap(~sample,ncol=4)+coord_equal(expand=FALSE)+scale_colour_manual(values=EPAL,guide='none',drop=FALSE)+theme_void(base_family=FONT)+theme(strip.text=element_text(size=7.1,face='bold'),panel.spacing=unit(2,'mm'))
c<-block(raster_plot(p,179,118),'c','Largest coefficient at each spot','All 25,548 spots; relative maxima without a confidence threshold')
page<-stack_rows(list(strip('Supplementary Figure S19 | Exploratory tumor-reference projection',9,'bold'),side_by_side(list(a,b),c(.65,.35)),c,key(EPAL,nrow=2),strip('Mouse-to-human tumor-reference projection does not establish cardiac ecotype identities. IQR describes spots within a tissue.',6.7)),c(8,77,119,9,7))
export_pages(list(page),'Figure_S19',224)
# One component per page; twelve tissues share the same scale and identical view geometry.
allpages<-lapply(components,function(ct){
 p<-ggplot(z,aes(xp,yp,colour=.data[[ct]]))+geom_point(size=.51,stroke=0)+facet_wrap(~sample,ncol=4)+coord_equal(expand=FALSE)+scale_colour_gradientn(colours=seqpal,limits=c(0,1),trans='sqrt',guide='none')+theme_void(base_family=FONT)+theme(strip.text=element_text(size=7.1,face='bold'),panel.spacing=unit(2,'mm'))
 ss<-summary[component==ct];q<-ggplot(ss,aes(sample,median_abundance,colour=role(as.character(sample))))+geom_linerange(aes(ymin=q25,ymax=q75),linewidth=.65)+geom_point(size=1.7)+scale_colour_manual(values=COL,guide='none')+scale_x_discrete(NULL,labels=SHORT)+scale_y_continuous('Coefficient · median (IQR)',limits=c(0,ceiling(max(summary$q75)*10)/10),breaks=scales::breaks_pretty(n=3))+BASE
 stack_rows(list(strip(paste0('Supplementary Figure S20 | ',ct,' across all 12 tissues'),9,'bold'),block(raster_plot(p,179,127),'a','Complete component map','Rows 1–2: formal comparison tissues; row 3: descriptive and auxiliary tissues'),colour_key(seqpal,c('0','0.1','0.3','1'),'Reference coefficient (shared sqrt colour scale)',sqrt(c(0,.1,.3,1))),block(q,'b','Within-tissue coefficient distributions'),strip('Coefficients come from the same completed projection. IQR describes spot heterogeneity, not uncertainty between animals.',6.7)),c(8,134,11,40,7))})
export_pages(allpages,'Figure_S20',rep(204,10))
dev.off()
