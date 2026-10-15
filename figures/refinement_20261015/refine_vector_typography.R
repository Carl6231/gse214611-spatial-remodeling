suppressPackageStartupMessages({library(xml2);library(rsvg);library(qpdf)})
a<-commandArgs(TRUE);root<-normalizePath(a[1]);out<-normalizePath(a[2]);dir.create(file.path(out,'Figures'),showWarnings=FALSE);dir.create(file.path(out,'Supplementary_figures'),showWarnings=FALSE)
files<-c(list.files(file.path(root,'09_JTM整合EcoTyper_20261015/Figures'),pattern='svg$',full.names=TRUE),list.files(file.path(root,'07_分析补全与组合图重排_20261014/Supplementary_figures'),pattern='svg$',full.names=TRUE))
for(f in files){
 d<-read_xml(f);nodes<-xml_find_all(d,'//*[local-name()="text"]');original<-xml_text(nodes);styles<-xml_attr(nodes,'style');sizes<-as.numeric(regmatches(styles,regexpr('[0-9.]+(?=px)',styles,perl=TRUE)))
 stopifnot(length(sizes)==length(nodes),all(is.finite(sizes)))
 # Keep exact text widths and positions; strengthen the smallest labels without moving scientific marks.
 new<-pmax(sizes,7.1);new[sizes>=7.5 & sizes<9]<-pmax(sizes[sizes>=7.5 & sizes<9],8)
 for(i in seq_along(nodes)){
  st<-sub('font-size: [0-9.]+px',sprintf('font-size: %.2fpx',new[i]),styles[i]);st<-gsub('#555D66','#424B55',st,fixed=TRUE)
  xml_set_attr(nodes[i],'style',st)
 }
 stopifnot(identical(xml_text(nodes),original))
 main<-grepl('^Figure_[1-8]\\.',basename(f));dest<-file.path(out,if(main)'Figures'else'Supplementary_figures',basename(f));write_xml(d,dest)
 pd<-sub('svg$','pdf',dest);rsvg_pdf(dest,pd)
 if(main)rsvg_png(dest,sub('svg$','png',dest),width=2161,height=2657)else rsvg_png(dest,sub('.svg','_render.png',dest,fixed=TRUE),width=2161)
 cat(basename(f),'labels',length(nodes),'raised',sum(new>sizes),'\n')
}
for(i in 1:18){f<-list.files(file.path(out,'Supplementary_figures'),pattern=paste0('^Figure_S',i,'_p[0-9]+\\.pdf$'),full.names=TRUE);if(length(f)){
 f<-f[order(as.integer(sub('.*_p([0-9]+)\\.pdf','\\1',f)))];qpdf::pdf_combine(f,file.path(out,'Supplementary_figures',paste0('Figure_S',i,'.pdf')))
}}
