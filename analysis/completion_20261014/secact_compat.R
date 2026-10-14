# Project-local SCOP 0.9.1 / SecAct 1.1.0 parameter-name compatibility.
# No package namespace or scientific parameter is changed.
stopifnot(as.character(packageVersion("scop")) == "0.9.1",
          as.character(packageVersion("SecAct")) == "1.1.0")
secact_compat_env <- new.env(parent = asNamespace("scop"))
for (name in c("secact_run_scrnaseq", "secact_run_st", "secact_run_matrix")) {
  f <- getFromNamespace(name, "scop")
  changed <- 0L
  repair <- function(x) {
    if (!is.call(x)) return(x)
    if (identical(x[[1L]], as.name("fun")) && "cores" %in% names(x)) {
      names(x)[names(x) == "cores"] <- "ncores"
      changed <<- changed + 1L
    }
    for (i in seq_along(x)[-1L]) if (is.call(x[[i]])) x[[i]] <- repair(x[[i]])
    x
  }
  body(f) <- repair(body(f))
  stopifnot(changed == 1L)
  environment(f) <- secact_compat_env
  assign(name, f, envir = secact_compat_env)
}
RunSecAct_compat <- getFromNamespace("RunSecAct", "scop")
environment(RunSecAct_compat) <- secact_compat_env
message("Compatibility: project-local RunSecAct, backend cores renamed ncores only")

# Preserve a valid zero-edge result in SecAct's native CCC empty-table branches.
ccc_env <- new.env(parent=asNamespace('SecAct'))
ccc <- getFromNamespace('SecAct.CCC.scRNAseq','SecAct')
txt <- paste(deparse(body(ccc)),collapse='\n')
rep_once <- function(old,new) {
 stopifnot(length(gregexpr(old,txt,fixed=TRUE)[[1]])==1L,grepl(old,txt,fixed=TRUE))
 txt <<- sub(old,new,txt,fixed=TRUE)
}
rep_once('smy_deg_comb <- data.frame()', 'smy_deg_comb <- data.frame(cellType=character(),secretedProtein=character(),exp_logFC=numeric(),exp_pv=numeric(),exp_pv.adj=numeric())')
rep_once('smy_act_comb <- data.frame()', 'smy_act_comb <- data.frame(cellType=character(),secretedProtein=character(),act_diff=numeric(),act_pv=numeric(),act_pv.adj=numeric())')
rep_once('ccc <- data.frame()', 'ccc <- data.frame(sender=character(),secretedProtein=character(),receiver=character(),sender_exp_logFC=numeric(),sender_exp_pv=numeric(),sender_exp_pv.adj=numeric(),receiver_act_diff=numeric(),receiver_act_pv=numeric(),receiver_act_pv.adj=numeric(),overall_strength=numeric(),overall_pv=numeric(),overall_pv.adj=numeric())')
rep_once('for (i in 1:nrow(smy_deg_comb))', 'for (i in seq_len(nrow(smy_deg_comb)))')
rep_once('if (ncol(ccc) > 0)', 'if (nrow(ccc) > 0)')
body(ccc)<-parse(text=txt)[[1]];environment(ccc)<-ccc_env
ccc_wrapper_env<-new.env(parent=asNamespace('scop'))
ccc_wrapper_env$get_namespace_fun<-function(package,fun) {
 if(package=='SecAct' && fun=='SecAct.CCC.scRNAseq')return(ccc)
 get('get_namespace_fun',envir=asNamespace('scop'),inherits=TRUE)(package,fun)
}
RunSecActCCC_compat<-getFromNamespace('RunSecActCCC','scop')
environment(RunSecActCCC_compat)<-ccc_wrapper_env
message('Compatibility: native CCC typed empty tables and zero-row guards; filters unchanged')
