source('11_图件结构重设计_20261015/scripts/design_common.R')
capture<-ragg::agg_capture(width=183,height=215,units='mm',res=300)
v<-readRDS(file.path(SRC,'Figure_5_6_values.rds'));xy<-readRDS(file.path(SRC,'niche_spatial_values.rds'));stopifnot(nrow(xy)==25548L)
d7<-SAMPLE[6:8]
# Figure 5: expose the saved covariate dependence before the descriptive signature.
j<-copy(v$joint);j[,tissue:=factor(tissue,levels=rev(levels(tissue)))];j[,lab:=sprintf('R² %.3f   |   CV %.3f   |   residual %.1f%%',r2,cross_validated_r2,100*residual_variance_fraction)]
pa<-ggplot(j,aes(y=tissue))+geom_segment(aes(x=0,xend=1,yend=tissue),linewidth=5,colour='#FBE4E8',lineend='butt')+geom_segment(aes(x=0,xend=r2,yend=tissue),linewidth=5,colour='#8F73CB',lineend='butt')+geom_text(aes(x=0,y=as.numeric(tissue)+.29,label=lab),hjust=0,size=2.3,family=FONT)+scale_x_continuous('Fraction of score variance',breaks=c(0,.5,1),labels=scales::label_percent())+labs(y=NULL)+BASE+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=7.2),plot.margin=margin(4,1,1,1,'mm'))
a<-block(pa,'a','Dependence on composition and RNA quality','Joint models: 2.3–4.1% residual variance; 25% criterion')
phys<-copy(v$physical);phys[,sample:=factor(sample_key,levels=rev(SAMPLE))]
pb<-ggplot(phys,aes(median_score,sample,colour=role))+geom_vline(xintercept=0,colour='#CED8DF',linewidth=.3)+geom_point(size=1.65)+scale_colour_manual(values=COL,guide='none')+scale_y_discrete(NULL,labels=LABEL)+labs(x='Median signature score')+BASE+theme(axis.text.y=element_text(size=6.7))
b<-block(pb,'b','Tissue-level expression','Each point represents one physical tissue')
sp<-fread(file.path(SRC,'phase_10_spot_signature_scores.tsv.gz'))
sp<-merge(xy[,.(sample_key,barcode,x,y)],sp[,.(sample_key,barcode,signature_score)],by=c('sample_key','barcode'))
stopifnot(nrow(sp)==25548L)
mp<-spatial_grid(sp,'signature_score',d7,c('#2E82C4','#FCFCFC','#EF637B'),c(-v$continuum_limit,v$continuum_limit),'Score',midpoint=0)
c<-block(raster_plot(mp,179,57),'c','Spatial organization of the Day 7 signature','Same pooled 1–99% symmetric display range for all three tissues')
pg<-ggplot(v$genes,aes(median_rho,gene,colour=pole))+geom_vline(xintercept=0,colour='#CED8DF',linewidth=.3)+geom_segment(aes(x=0,xend=median_rho,yend=gene),linewidth=.6)+geom_point(size=1.8)+scale_colour_manual(values=c(Contractile='#2B83D5',Fibrotic='#EF637B'),guide='none')+scale_x_continuous('Median within-tissue ρ',limits=c(-.57,.57),breaks=c(-.5,0,.5))+labs(y=NULL)+BASE+theme(axis.text.y=element_text(face='italic',size=7.3))
d<-block(pg,'d','Nonconstruction genes','Descriptive correlations')
pw<-copy(v$pathways);pw[,term2:=factor(paste(collection_label,term,sep=' · '),levels=rev(paste(collection_label,term,sep=' · ')))]
pp<-ggplot(pw,aes(NES,term2))+geom_vline(xintercept=0,colour='#CED8DF',linewidth=.3)+geom_segment(aes(x=0,xend=NES,yend=term2),colour='#A7B7C4',linewidth=.5)+geom_point(aes(colour=NES>0),size=1.8)+geom_text(aes(label=sprintf('%.2f',NES)),hjust=ifelse(pw$NES<0,1.2,-.2),size=2.25,family=FONT)+scale_colour_manual(values=c('FALSE'='#2B83D5','TRUE'='#EF637B'),guide='none')+scale_x_continuous('Normalized enrichment score',limits=c(-5,4.1),breaks=c(-4,-2,0,2,4))+scale_y_discrete(NULL,labels=function(x)vapply(x,function(y)paste(strwrap(y,36),collapse='\n'),character(1)))+BASE+theme(axis.text.y=element_text(size=6.7,lineheight=.92))
e<-block(pp,'e','Pathway polarity','Nine prespecified pathways')
page<-stack_rows(list(side_by_side(list(a,b),c(.61,.39)),c,side_by_side(list(d,e),c(.34,.66)),strip('Purple bars: joint-model R²; pale coral: residual. Spatial CV is within tissue, not validation in independent animals.',6.8)),c(67,65,74,5))
export(page,'Figure_5',217)
# Figure 6: recurrence, candidates and perturbations occupy distinct evidence groups.
mer<-copy(v$meringue)
pm<-ggplot(mer,aes(feature,sample_key,fill=statistic))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=rank),size=2.3,family=FONT)+scale_fill_gradientn(colours=c('#EDF8F5','#80D2C1','#158A8D'),name='Spatial statistic')+labs(x=NULL,y=NULL)+BASE+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(face='italic',angle=50,hjust=1,size=7),legend.key.height=unit(3,'mm'))
a<-block(pm,'a','Spatial recurrence across all three Day 7 tissues','Sixteen genes in every top 100; numbers denote within-tissue rank')
map_genes<-lapply(c('Sdf4','Abca1'),function(g){q<-spatial_grid(xy,g,d7,c('#EEF8F5','#6CCDBB','#157B8D'),v$gene_limits[[g]],'log1p');raster_plot(q,120,38)})
b<-block(stack_rows(list(heading(map_genes[[1]],'Sdf4'),heading(map_genes[[2]],'Abca1')),c(1,1)),'b','Representative recurrent genes','Shared display range within each gene; complete gene atlas in Supplementary Fig. S10')
ax<-copy(v$axes);ax[,row:=factor(as.character(label),levels=rev(as.character(label)))]
pax<-ggplot(ax,aes(y=row))+geom_hline(yintercept=1:10,colour='#E8EEF1',linewidth=.2)+geom_text(aes(x=1,label=paste0(n_supported_stages,'/2')),size=2.25)+geom_text(aes(x=2,label=paste0(n_discovery_supporting_tissues,'/6')),size=2.25)+geom_tile(aes(x=3,fill=minimum_component_coverage),width=.8,height=.85)+geom_text(aes(x=3,label=scales::percent(minimum_component_coverage,accuracy=1)),size=2.2)+geom_text(aes(x=4,label=ifelse(neighborhood,'●','○')),size=2.2,colour='#179C88')+geom_text(aes(x=5,label=ifelse(gradient=='TRUE','●',ifelse(gradient=='FALSE','○','NE'))),size=2.2,colour='#8568C6')+geom_text(aes(x=6,label=sprintf('%.2g',spatialcellchat_median_score)),size=2.2)+scale_fill_gradient(low='#EFF8F6',high='#82C9BD',limits=c(0,1),guide='none')+scale_x_continuous(NULL,breaks=1:6,labels=c('Stage','Tissue','Coverage','Neighbor','Gradient','CellChat\nscore'),position='top',limits=c(.5,6.6))+labs(y=NULL)+BASE+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=6.3),axis.text.y=element_text(size=6.6))
c<-block(pax,'d','Candidate ligand–receptor evidence','● supported; ○ unsupported; NE not evaluable')
or<-v$oracle
po<-ggplot(or,aes(median_displacement,target))+geom_segment(aes(x=0,xend=median_displacement,yend=target),colour='#CFD8DF',linewidth=.6)+geom_point(colour='#8568C6',size=1.8)+geom_text(aes(label=sprintf('%.3f',median_displacement)),hjust=-.25,size=2.25)+scale_x_continuous('Median displacement',limits=c(0,.78),breaks=c(0,.3,.6))+labs(y=NULL)+BASE+theme(axis.text.y=element_text(face='italic'))
d<-block(stack_rows(list(po,strip('Transition mass = 0.500 for all six.\nNone met the module-consistency rule.',6.6)),c(1,.18)),'c','In silico perturbations','Six targets; model embedding units')
page<-stack_rows(list(a,side_by_side(list(b,plot_spacer(),d),c(.68,.025,.295)),c,strip('Association, model-based candidate evidence and perturbation predictions remain distinct; no causal validation is implied.',6.8)),c(43,96,64,7))
export(page,'Figure_6',216)
# Figure 8: two paired matrices show all 72 vector coverages without repeated axes.
s<-readRDS(file.path(SRC,'Figure_8_values.rds'));lig<-c('C3','CCL8','CXCL12','GAS6','MDK','PTN')
activity<-function(d,ids,labels,title){d<-d[secretedProtein%in%lig];d[,`:=`(s=factor(sample_id,levels=ids),g=factor(secretedProtein,levels=rev(lig)))];lim<-max(abs(d$zscore))
 q<-ggplot(d,aes(s,g,fill=zscore))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=sprintf('%.1f',ifelse(abs(zscore)<.05,0,zscore))),size=2.35)+scale_fill_gradient2(low='#2B83D5',mid='#FCFCFC',high='#EF637B',limits=c(-lim,lim),name='Activity z')+scale_x_discrete(NULL,labels=labels)+labs(y=NULL)+BASE+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=40,hjust=1),legend.position='right',legend.key.height=unit(3,'mm'))
 heading(q,title)}
a<-block(side_by_side(list(activity(s$sn,s$sn_ids,s$labels_sn,'Myeloid-labelled nuclei'),activity(s$st,s$st_ids,s$labels_st,'Spatial tissue medians'))),'a','Activity estimates for six candidate ligands','Nucleus samples and spatial tissues are not animal-paired; separate activity scales')
coverage<-function(mode){d<-s$v[signal_mode==mode];d[,`:=`(sm=factor(sample_id,levels=s$st_ids),lg=factor(ligand_secact_human,levels=rev(lig)))];q<-ggplot(d,aes(sm,lg,fill=vector_fraction))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=sprintf('%.0f%%',100*vector_fraction)),size=2.3)+scale_fill_gradientn(colours=c('#F0F9F7','#83D4C4','#12888D'),limits=c(0,1),name='Spot\nfraction')+scale_x_discrete(NULL,labels=s$labels_st)+labs(y=NULL)+BASE+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=40,hjust=1),legend.key.height=unit(3,'mm'));heading(q,if(mode=='sending')'Sending vectors' else 'Receiving vectors')}
b<-block(side_by_side(list(coverage('sending'),coverage('receiving'))),'b','Spatial coverage of the inferred vector fields','All six ligands × six tissues × two directions; arrows in Supplementary Fig. S18')
stopifnot(nrow(s$cc)==36L,all(s$cc$N==0L))
zero<-wrap_elements(full=grobTree(rectGrob(x=.5,y=.5,width=.96,height=.83,gp=gpar(fill='#F1F6F8',col=NA)),textGrob('0',x=.14,y=.58,gp=gpar(fontfamily=FONT,fontsize=31,fontface='bold',col='#687F91')),textGrob('retained candidate-ligand edges',x=.32,y=.63,just='left',gp=gpar(fontfamily=FONT,fontsize=8,fontface='bold')),textGrob('All 36 ligand × nucleus-sample combinations\nSix ligands; six samples; all class pairs included',x=.32,y=.32,just='left',gp=gpar(fontfamily=FONT,fontsize=7,lineheight=1.2))))
c<-block(zero,'c','Candidate-edge result')
m<-melt(s$meta,id.vars='sample_id',variable.name='metric');m[,sm:=factor(sample_id,levels=s$st_ids)]
q<-ggplot(m,aes(sm,value,colour=grepl('D7',sample_id)))+geom_point(size=1.7)+facet_wrap(~metric,nrow=1,scales='free_y')+scale_y_continuous(breaks=function(x)unique(round(pretty(x,3))))+scale_colour_manual(values=c('FALSE'=unname(COL['MI Day 3']),'TRUE'=unname(COL['MI Day 7'])),guide='none')+scale_x_discrete(NULL,labels=s$labels_st)+labs(y=NULL)+BASE+theme(axis.text.x=element_text(angle=40,hjust=1))
d<-block(q,'d','Tissue-local spatial factors','All independently fitted patterns in Supplementary Fig. S17')
page<-stack_rows(list(a,b,side_by_side(list(c,d),c(.48,.52)),strip('Activity, edges and vector coverage are model outputs; they do not measure secretion, receptor activation or cell movement.',6.8)),c(58,59,46,5))
export(page,'Figure_8',174)
dev.off()
