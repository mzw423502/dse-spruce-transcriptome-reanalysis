#!/usr/bin/env Rscript
# Portable, read-only replay entry point for the revision's derived outputs.
# Arguments: <frozen_base> <primary_DESeq2_model.rds> <new_output_directory>
# This script does not quantify reads or fit an RNA-seq model. It recreates the
# revision classification/expression exports from frozen result TSVs and an
# explicitly supplied, already fitted primary DESeq2 model.
suppressPackageStartupMessages({library(DESeq2); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
BASE <- normalizePath(args[1], mustWork=TRUE)
RDS <- normalizePath(args[2], mustWork=TRUE)
OUT <- normalizePath(args[3], mustWork=TRUE)
B <- file.path(BASE, 'reproducibility', 'formal_statistics', 'results')
F <- BASE
stopifnot(dir.exists(B), file.exists(file.path(F,'numerical','figure_response.tsv')))
rd <- function(p) read.delim(p, check.names=FALSE, stringsAsFactors=FALSE)
wt <- function(d,n) write.table(d,file.path(OUT,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
bh <- function(p,n=length(p)){v=rep(NA_real_,length(p)); k=is.finite(p); v[k]=p.adjust(p[k],'BH',n=n); v}

P <- rd(file.path(F,'numerical','figure_response.tsv'))
P <- P[P$risk112==0,]
stopifnot(nrow(P)==41139,!anyDuplicated(P$feature_id))
sel <- P[P$response_description=='SUPPORTED_OPPOSITE_DIRECTIONS',]
stopifnot(nrow(sel)==391)
wt(sel,'primary_391.tsv')
rows <- list(); summaries <- list()
for(route in c('allN','competition','edgeR_QL')){
  if(route=='edgeR_QL'){
    z <- rd(file.path(B,'primary','edgeR_QL_sensitivity.tsv'))
    ix <- z$contrast=='interaction'; ss <- !ix; m <- sum(ix)
    stopifnot(m==41221,sum(ss)==4*m,!anyDuplicated(z[ix,'feature_id']))
    z$q_joint4 <- NA_real_; z$q_joint4[ss] <- bh(z$PValue[ss],4*m)
    stopifnot(max(abs(bh(z$PValue[ix],m)-z$padj_contrast_fixed[ix]),na.rm=TRUE)<1e-10)
    w <- data.frame(feature_id=z$feature_id[ix],stringsAsFactors=FALSE)
    for(nm in c('P0','P1','I')){
      cn <- c(P0='pathogen_at_DSE0',P1='pathogen_at_DSE1',I='interaction')[[nm]]
      a <- z[z$contrast==cn,]; a <- a[match(w$feature_id,a$feature_id),]
      w[[paste0('lfc_',nm)]] <- a$logFC; w[[paste0('p_',nm)]] <- a$PValue
      w[[paste0('q_',nm)]] <- if(nm=='I') a$padj_contrast_fixed else a$q_joint4
    }
    finite <- apply(w[,c('lfc_P0','lfc_P1','lfc_I')],1,function(x)all(is.finite(x)))
    qs <- apply(w[,c('q_P0','q_P1','q_I')],1,function(x)all(is.finite(x)&x<.05))
    gate <- apply(w[,c('p_P0','p_P1','p_I')],1,function(x)all(is.finite(x)&x<.05/3))
    w$classification <- 'UNCERTAINTY_OR_CHANGE_UNRESOLVED'; good <- finite&qs&gate
    w$classification[good&w$lfc_P0*w$lfc_P1<0] <- 'SUPPORTED_OPPOSITE_DIRECTIONS'
    w$classification[good&w$lfc_P0*w$lfc_P1>0&w$lfc_I*w$lfc_P0>0] <- 'SUPPORTED_SAME_DIRECTION_LARGER_MAGNITUDE'
    w$classification[good&w$lfc_P0*w$lfc_P1>0&w$lfc_I*w$lfc_P0<0] <- 'SUPPORTED_SAME_DIRECTION_SMALLER_MAGNITUDE'
    w$uncertainty_gate <- gate; w$q_gate <- qs
    w$rule <- 'QL_ANALOGUE: three q<0.05 (m/4m families), three QL p<0.05/3; no reconstructed Wald CI'
  } else {
    z <- rd(file.path(B,route,'pathogen_response_plane.tsv'))
    w <- data.frame(feature_id=z$feature_id,lfc_P0=z$log2FoldChange_P0,lfc_P1=z$log2FoldChange_P1,lfc_I=z$log2FoldChange_I,q_P0=z$padj_simple_four_fixed_P0,q_P1=z$padj_simple_four_fixed_P1,q_I=z$padj_I,classification=z$response_description)
    w$q_gate <- apply(w[,c('q_P0','q_P1','q_I')],1,function(x)all(is.finite(x)&x<.05))
    gate <- rep(TRUE,nrow(z))
    for(nm in c('P0','P1','I')){lo=z[[paste0('joint95_low_',nm)]];hi=z[[paste0('joint95_high_',nm)]];gate=gate&is.finite(lo)&is.finite(hi)&(lo>0|hi<0)}
    w$uncertainty_gate <- gate; w$rule <- 'ORIGINAL_DESeq2_RULE: three q<0.05 and 98.3333% marginal Wald intervals exclude zero'
  }
  stopifnot(!anyDuplicated(w$feature_id))
  idx <- match(sel$feature_id,w$feature_id); a <- w[idx,]; a$feature_id <- sel$feature_id
  a$route <- route; a$present <- !is.na(idx)
  a$finite <- apply(a[,c('lfc_P0','lfc_P1','lfc_I')],1,function(x)all(is.finite(x)))
  a$same_P0 <- is.finite(a$lfc_P0)&sign(a$lfc_P0)==sign(sel$log2FoldChange_P0)
  a$same_P1 <- is.finite(a$lfc_P1)&sign(a$lfc_P1)==sign(sel$log2FoldChange_P1)
  a$same_I <- is.finite(a$lfc_I)&sign(a$lfc_I)==sign(sel$log2FoldChange_I)
  a$original_sign_pattern <- a$same_P0&a$same_P1&a$same_I
  a$opposite_point_estimates <- a$finite&a$lfc_P0*a$lfc_P1<0
  a$retained_supported_opposite <- a$present & !is.na(a$classification) & a$classification=='SUPPORTED_OPPOSITE_DIRECTIONS'
  a$retained_supported_original_pattern <- a$retained_supported_opposite&a$original_sign_pattern
  a$status <- 'POINT_PATTERN_CHANGED'; a$status[a$original_sign_pattern] <- 'ORIGINAL_PATTERN_WITHOUT_FULL_SUPPORT'; a$status[a$retained_supported_original_pattern] <- 'SUPPORTED_ORIGINAL_OPPOSITE'; a$status[!a$finite] <- 'NONFINITE_EFFECT'; a$status[!a$present] <- 'NOT_ELIGIBLE_OR_ABSENT'
  summaries[[route]] <- data.frame(route=route,original_n=nrow(sel),present=sum(a$present),finite=sum(a$finite),same_P0=sum(a$same_P0),same_P1=sum(a$same_P1),same_I=sum(a$same_I),same_original_pattern=sum(a$original_sign_pattern),opposite_points=sum(a$opposite_point_estimates),q_gate=sum(a$q_gate%in%TRUE),uncertainty_gate=sum(a$uncertainty_gate%in%TRUE),supported_opposite=sum(a$retained_supported_opposite),supported_original_pattern=sum(a$retained_supported_original_pattern),missing=sum(!a$present),finite_pattern_changed=sum(a$finite&!a$original_sign_pattern))
  wt(a,paste0('primary391_',route,'_classification.tsv')); rows[[route]] <- a; wt(w,paste0('all_units_',route,'_classification.tsv'))
}
summary <- do.call(rbind,summaries); rownames(summary) <- NULL; wt(summary,'classification_retention_summary.tsv')
dds <- readRDS(RDS)
stopifnot(all(sel$feature_id%in%rownames(dds)),ncol(dds)==12L)
nc <- counts(dds,normalized=TRUE); raw <- counts(dds,normalized=FALSE); meta <- as.data.frame(colData(dds)); stopifnot(all(c('DSE','pathogen')%in%colnames(meta)))
allrows <- lapply(seq_len(ncol(dds)),function(j)data.frame(feature_id=sel$feature_id,sample_id=colnames(dds)[j],DSE=as.character(meta$DSE[j]),pathogen=as.character(meta$pathogen[j]),group=c('Pa','PaPs','PaHp','PaPsHp')[1+as.integer(as.character(meta$DSE[j]))+2*as.integer(as.character(meta$pathogen[j]))],raw_model_count=raw[sel$feature_id,j],normalized_count=nc[sel$feature_id,j],log2_normalized_plus1=log2(nc[sel$feature_id,j]+1)))
wt(do.call(rbind,allrows),'primary391_expression_12pools.tsv')
allcontr <- rd(file.path(B,'primary','DESeq2_all_contrasts.tsv')); wt(allcontr[allcontr$feature_id%in%sel$feature_id,],'primary391_five_contrasts.tsv')
write_json(list(status='COMPLETE',scope='Portable read-only derived classification sensitivity and expression extraction; no refit',arguments=list(frozen_base=BASE,primary_DESeq2_model_rds=RDS,output=OUT),primary_n=391,primary_nonrisk_eligible=41139,summary=summary,edgeR_rule='Method-native sensitivity analogue; not the identical DESeq2 Wald interval classifier. No category-level FDR guarantee.',normalization='Frozen full-data DESeq2 normalization factors',session=capture.output(sessionInfo())),file.path(OUT,'analysis_summary.json'),pretty=TRUE,auto_unbox=TRUE)
cat('COMPLETE\n'); print(summary)
