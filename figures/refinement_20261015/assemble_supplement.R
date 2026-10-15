library(grid);library(gridExtra);library(jsonlite);library(png)
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
out<-normalizePath(args[1]);root<-dirname(out)
figdir<-file.path(out,'Supplementary_figures')
caps<-fromJSON(file.path(root,'07_分析补全与组合图重排_20261014/recovered_source/supplement_captions.json'))
caps$S17_1 <- 'Supplementary Figure S17. Independently fitted SecAct spatial patterns in the three Day 3 tissues. Each tissue has equal-size map slots; all fitted factors are shown, and blank slots denote factors not fitted at the selected rank. Colors show loading divided by each factor-specific maximum (0-1). Factor numbers do not identify homologous programs across tissues.'
caps$S17_2 <- 'Supplementary Figure S17 (continued). All independently fitted Day 7 spatial factors, using the same map slots and within-factor scaling as the Day 3 page. Each row represents one tissue; factor numbers are tissue-local and are not matched between rows.'
caps$S19_1 <- 'Supplementary Figure S19. Exploratory SpatialEcoTyper tumor-reference projection. a, Tissue-level mean reference coefficients. b, The largest reference component at each spot in the eight formal-comparison tissues. c, The corresponding maps in the two descriptive timepoints and two auxiliary needle-pass tissues. All 25,548 spots are represented. Relative maxima are displayed without a confidence threshold; reference labels are not cardiac cell identities. Coefficients sum to one by model constraint.'
for (j in 1:10) caps[[paste0('S20_',j)]] <- paste0('Supplementary Figure S20 (',j,'/10). ',c('NonSE',paste0('SE',1:9))[j],' reference coefficients across all 12 tissues. a, Eight formal-comparison tissues. b, Descriptive and auxiliary tissues. c, Within-tissue medians and interquartile ranges of the saved spot coefficients; these ranges describe spot heterogeneity, not biological-replicate uncertainty. All spatial panels share a 0-1 coefficient scale. Colors of tissue labels and summary points follow the study groups. Original coordinates and all 120 component-by-tissue maps are retained. This cross-species tumor-reference projection does not establish cardiac ecotype identities.')
titles<-c('Tissue registry and corrected center counts','Original labels and corrected reference composition',
'Distinct denominators in original and corrected niche analyses','Sixteen nontechnical MERINGUE genes in the Day 7 top 100 of all three tissues',
'State-specific CellOracle biological-sample summaries','Evaluability of frozen macrophage modules')
notes<-c(
'The 1-hour and 4-hour tissues are descriptive; needle-pass tissues are auxiliary. GSM6613080 and GSM6613081 each contain two resolved tissues.',
'Initial labels do not establish lineage purity. Both full references contain 176,526 nuclei; the balanced reference contains 18,849. Classes not separately labelled initially were not biologically absent.',
'The initial analysis contains 5,937 centers and corrected S054 contains 5,151. Reference and eligibility schemes differ; this is not a temporal change or tissue loss. Corrected thresholds are total myeloid 0.20, radius multiplier 1.2, program 0.10 and margin 0.15.',
'Order follows median rank, worst rank, median statistic and symbol. Statistics are descriptive; missing p-values and q-values were not reconstructed.',
'Each value is the median of saved within-sample medians. All six targets are NEGATIVE_COMPLETE. Displacement is in model embedding units; transition mass is a model probability, not migration, measured fate or experimental knockout.',
'Modules require at least five evaluable genes. The consistency rule requires at least eight samples, absolute sample-median shift of 0.10 baseline SD or more, and the same sign in at least two-thirds of samples. Evaluable does not mean this criterion was met.')
page_count<-0
page<-function(){grid.newpage();page_count<<-page_count+1L;grid.text(page_count,y=unit(9,'mm'),gp=gpar(fontsize=9,fontfamily='Times'))}
wrapped<-function(x,width=125)paste(strwrap(x,width=width),collapse='\n')
text_at<-function(s,y,size=9,bold=FALSE){grid.text(wrapped(s),x=unit(15,'mm'),y=unit(297-y,'mm'),just=c('left','top'),gp=gpar(fontsize=size,fontfamily='Times',fontface=if(bold)'bold'else'plain',lineheight=1.15))}
# Recreate numerical tables from the supplied CSVs; no screenshot masking or image editing.
table_at<-function(i,y){
 d<-read.csv(file.path(root,'05_JTM投稿定稿_20261014/Additional_files',paste0('Additional_file_',i+1,'_Table_S',i,'.csv')),check.names=FALSE,colClasses='character')
 d[]<-lapply(d,function(x){idx<-grepl('^[0-9]{1,3}(,[0-9]{3})+$',x);x[idx]<-gsub(',','',x[idx],fixed=TRUE);x})
 text_at(paste0('Supplementary Table S',i,'. ',titles[i]),y,10,TRUE)
 widths<-switch(as.character(i),'1'=c(58,48,38,36),'2'=c(48,44,44,44),'3'=c(83,30,35,32),'4'=c(42,23,23,23,23,23,23),'5'=c(30,48,24,36,42),'6'=c(70,45,65))
 headers<-vapply(names(d),function(x)wrapped(x,width=if(i==4)10 else 19),character(1))
 names(d)<-headers
 tg<-tableGrob(d,rows=NULL,theme=ttheme_minimal(base_size=8,base_family='Times',
 core=list(fg_params=list(hjust=0,x=.025),bg_params=list(fill='white',col=NA)),
 colhead=list(fg_params=list(hjust=0,x=.025,fontface='bold'),bg_params=list(fill='white',col=NA)),padding=unit(c(2.3,1.4),'mm')))
 tg$widths<-unit(widths,'mm');hh<-sum(convertHeight(tg$heights,'mm',valueOnly=TRUE));head_h<-convertHeight(tg$heights[1],'mm',valueOnly=TRUE)
 pushViewport(viewport(x=unit(15,'mm'),y=unit(297-y-10,'mm'),width=unit(sum(widths),'mm'),height=unit(hh,'mm'),just=c('left','top')))
 grid.draw(tg)
 for(yy in c(0,hh-head_h,hh))grid.lines(x=unit(c(0,1),'npc'),y=unit(c(yy,yy),'mm'),gp=gpar(lwd=if(yy==hh-head_h).6 else 1))
 upViewport()
 text_at(notes[i],y+12+hh,8)
 y+12+hh+length(strwrap(notes[i],125))*3.5+8
}
# Every requested numbered supplement must be present before final assembly.
files<-lapply(1:20,function(i){
 p<-list.files(figdir,pattern=paste0('^Figure_S',i,'(_p[0-9]+)?_render\\.png$'),full.names=TRUE)
 if(!length(p))stop('Missing supplement S',i)
 ix<-ifelse(grepl('_p',basename(p)),as.integer(sub('.*_p([0-9]+)_render.*','\\1',basename(p))),1)
 p[order(ix)]
})
final<-file.path(out,'JTM_Supplementary_Refined_20261015.pdf')
cairo_pdf(final,width=210/25.4,height=297/25.4,onefile=TRUE)
# A reading guide gives each supplementary figure a distinct role beside the main figures.
page()
text_at('Supplementary figure guide',15,16,TRUE)
text_at('Spatial remodeling and tissue composition after myocardial infarction in mice',27,10)
text_at('Complete atlases, model sensitivity and evidence boundaries',37,10,TRUE)
guide_rows<-list(
 c('S1  |  Fig. 1','Tissue images','All 12 tissues: H&E and retained spots'),
 c('S2  |  Fig. 1','Quality and compatibility','Per-tissue QC and initial reference support'),
 c('S3  |  Fig. 2','Initial deconvolution','Complete two-method maps and confidence'),
 c('S4  |  Fig. 4','Spatial structure','Saved tissue-local networks and domains'),
 c('S5  |  Fig. 4','Spatial features','Frozen feature ranks and correlations'),
 c('S6  |  Fig. 3','Corrected broad reference','Lineage, markers, support; complete 11-class maps'),
 c('S7  |  Fig. 3','Fine-state uncertainty','Macrophage-state maps and method sensitivity'),
 c('S8  |  Fig. 4','Corrected primary niches','All 12 tissues, fractions and center counts'),
 c('S9  |  Fig. 5','Expression and covariates','Models A–D, cellular sources and fixed-gene atlas'),
 c('S10 |  Fig. 6','Gene recurrence','Complete spatial-gene maps and rankings'),
 c('S11 |  Fig. 6','Network and module limits','Cross-correlation and tissue-local modules'),
 c('S12 |  Fig. 6','Candidate-axis evidence','All ten frozen ligand–receptor candidates'),
 c('S13 |  Fig. 6','Perturbation results','All six targets; negative and unevaluable results'),
 c('S14 |  Fig. 4','Original-reference sensitivity','Prior summary and complete original niche atlas'),
 c('S15 |  Fig. 7','Contextual associations','Formal, descriptive/auxiliary and external tissues'),
 c('S16 |  Fig. 8','Secreted-protein activity','Complete ligand and class/tissue summaries'),
 c('S17 |  Fig. 8','Tissue-local patterns','All independently selected spatial factors'),
 c('S18 |  Fig. 8','Paired vector displays','One ligand per page: sending above receiving'),
 c('S19 |  Results','Reference projection','All tissues; formal and supporting groups separated'),
 c('S20 |  Results','Complete component atlas','One component per page; all 12 tissues + summaries'))
figure_starts <- 5L+c(0L,head(cumsum(lengths(files)),-1L))
figure_ends <- figure_starts+lengths(files)-1L
for(k in seq_along(guide_rows)){
 y<-53+(k-1)*9.7
 grid.lines(x=unit(c(15,195),'mm'),y=unit(c(297-y+2,297-y+2),'mm'),gp=gpar(col='#DFE5EA',lwd=.5))
 grid.text(guide_rows[[k]][1],x=unit(15,'mm'),y=unit(297-y,'mm'),just=c('left','top'),gp=gpar(fontsize=8.5,fontfamily='sans',fontface='bold',col='#428BD0'))
 grid.text(paste0('pp. ',figure_starts[k],'–',figure_ends[k]),x=unit(195,'mm'),y=unit(297-y,'mm'),just=c('right','top'),gp=gpar(fontsize=8,fontfamily='sans',col='#505A65'))
 grid.text(guide_rows[[k]][2],x=unit(45,'mm'),y=unit(297-y,'mm'),just=c('left','top'),gp=gpar(fontsize=8.5,fontfamily='sans',fontface='bold'))
 grid.text(guide_rows[[k]][3],x=unit(45,'mm'),y=unit(297-y-4,'mm'),just=c('left','top'),gp=gpar(fontsize=8,fontfamily='sans',col='#505A65'))
}
text_at('Reading order: the six supplementary tables precede Figures S1–S20. S8 is the corrected primary analysis; S14 is the original-reference sensitivity analysis. Each figure retains its own page numbering in the caption.',251,9)
text_at('Upload volumes: Part 1, guide + Tables S1–S6 + Figures S1–S3; Part 2, S4–S6; Part 3, S7–S9; Part 4, S10–S14; Part 5, S15–S18 (Additional file 14); Part 6, S19–S20 (Additional file 15).',269,8)

for(i in seq(1,6,by=2)){page();y<-table_at(i,15);table_at(i+1,y+5)}
for(i in 1:20)for(j in seq_along(files[[i]])){
 page();im<-readPNG(files[[i]][j]);ratio<-dim(im)[2]/dim(im)[1];w<-min(184,238*ratio);h<-w/ratio
 grid.raster(im,x=unit(105,'mm'),y=unit(287-h/2,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),interpolate=FALSE)
 key<-paste0('S',i,'_',j)
 caption<-caps[[key]]
 stopifnot(length(caption)==1L,nzchar(caption))
 text_at(gsub('\n',' ',caption),min(252,16+h),9)
}
dev.off();cat('SUPPLEMENT_PAGES',page_count,'\n')

# Lossless document splitting at complete figure-group boundaries for JTM's 20 MB limit.
groups<-list(1:3,4:6,7:9,10:14,15:18,19:20)
attachment_ids<-c(1,2,3,4,14,15)
page_for<-list();start<-5L
for(i in 1:20){page_for[[i]]<-seq.int(start,length.out=length(files[[i]]));start<-start+length(files[[i]])}
for(i in seq_along(groups)){
 pages<-unlist(page_for[groups[[i]]]);if(i==1)pages<-c(1:4,pages)
 dest<-file.path(out,'Additional_files',paste0('Additional_file_',attachment_ids[i],'_Supplementary_Part_',i,'.pdf'))
 qpdf::pdf_subset(final,pages=pages,output=dest)
 stopifnot(file.info(dest)$size<20e6)
 cat('UPLOAD_PART',i,'pages',length(pages),'bytes',file.info(dest)$size,'\n')
}
