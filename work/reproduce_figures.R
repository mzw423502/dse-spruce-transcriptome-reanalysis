# Publication graphics rendered on Galaxy/Kylin only.  This script consumes
# already-computed tables and performs no inferential model fit or reclassification.
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript draw_revision_figures.R <analysis_dir> <empty_output_dir>")
root <- normalizePath(args[1],mustWork=TRUE)
out <- args[2]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Output directory must be empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

group_cols <- c(Pa="#7F7F7F", PaPs="#0072B2", PaHp="#E69F00", PaPsHp="#009E73")
group_labs <- c(Pa="Neither fungus", PaPs="DSE only", PaHp="Pathogen only", PaPsHp="Both fungi")

apply_figure_style <- function() {
  par(family="sans", las=1, tcl=-0.25, mgp=c(2.0,0.55,0), xaxs="i", yaxs="i",
      bty="l", cex.axis=0.80, cex.lab=0.90, cex.main=0.90)
}
panel_letter <- function(letter) mtext(letter, side=3, adj=0, line=2.75, font=2, cex=1.08)
open_device <- function(stem, width, height, kind) {
  f <- file.path(out, paste0(stem, ".", kind))
  if (kind == "pdf") pdf(f, width=width, height=height, useDingbats=FALSE, pointsize=11)
  if (kind == "svg") svg(f, width=width, height=height, pointsize=11)
  if (kind == "png") png(f, width=width, height=height, units="in", res=300, pointsize=11, type="cairo")
}

edge <- read.delim(file.path(root, "primary391_edgeR_QL_classification.tsv"), check.names=FALSE)
alln <- read.delim(file.path(root, "primary391_allN_classification.tsv"), check.names=FALSE)
comp <- read.delim(file.path(root, "primary391_competition_classification.tsv"), check.names=FALSE)
expr <- read.delim(file.path(root, "primary391_expression_12pools.tsv"), check.names=FALSE)
cts <- read.delim(file.path(root, "primary391_five_contrasts.tsv"), check.names=FALSE)

draw_fig1 <- function() {
  par(oma=c(0,0,0,0))
  layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), widths=c(1.45,1), heights=c(1,0.92))
  # a: P0/P1 response plane
  apply_figure_style(); par(mar=c(3.4,3.65,4.25,0.5))
  p0 <- cts[cts$contrast == "pathogen_at_DSE0", c("feature_id", "log2FoldChange")]
  p1 <- cts[cts$contrast == "pathogen_at_DSE1", c("feature_id", "log2FoldChange")]
  names(p0)[2] <- "lfc_P0_DESeq2"; names(p1)[2] <- "lfc_P1_DESeq2"
  main_est <- merge(p0, p1, by="feature_id", sort=FALSE)
  plane <- merge(main_est, edge[,c("feature_id","retained_supported_original_pattern")], by="feature_id", sort=FALSE)
  keep <- plane$retained_supported_original_pattern == "TRUE"
  xr <- range(plane$lfc_P0_DESeq2, finite=TRUE); yr <- range(plane$lfc_P1_DESeq2, finite=TRUE)
  xpad <- diff(xr)*0.05; ypad <- diff(yr)*0.05
  plot(plane$lfc_P0_DESeq2, plane$lfc_P1_DESeq2, xlim=xr+c(-xpad,xpad), ylim=yr+c(-ypad,ypad), pch=16, cex=0.33, col=grDevices::adjustcolor("grey70",0.55),
       xlab=expression(P[0]*" pathogen effect without DSE (log"[2]*"FC)"),
       ylab=expression(P[1]*" pathogen effect with DSE (log"[2]*"FC)"),
       main="Primary DESeq2 effects")
  abline(h=0,v=0,lty=3,col="grey55",lwd=0.8)
  points(plane$lfc_P0_DESeq2[keep], plane$lfc_P1_DESeq2[keep], pch=16, cex=0.47, col="#0072B2")
  legend("topright", legend=c("QL support (132)","Other selected units (259)"),
         pch=16, col=c("#0072B2",grDevices::adjustcolor("grey55",0.8)), bty="n", cex=0.75)
  panel_letter("a")
  # b: direction distribution
  apply_figure_style(); par(mar=c(3.4,3.65,4.25,0.5))
  posneg <- sum(edge$lfc_P0 > 0 & edge$lfc_P1 < 0)
  negpos <- sum(edge$lfc_P0 < 0 & edge$lfc_P1 > 0)
  bp <- barplot(c(posneg,negpos), names.arg=c(expression(P[0]*">0, P"[1]*"<0"), expression(P[0]*"<0, P"[1]*">0")),
                col=c("#0072B2", "#E69F00"), border=NA, ylim=c(0,380),
                ylab="Sequence units", main="Response directions")
  text(bp,c(posneg,negpos)+11,c(posneg,negpos),cex=0.8)
  panel_letter("b")
  # c: route retention
  apply_figure_style(); par(mar=c(3.4,3.65,4.25,0.5))
  support <- c(373,391,132); available_not_full <- c(16,0,259); missing <- c(2,0,0)
  y <- barplot(rbind(support,available_not_full,missing), names.arg=c("allN","Competition","edgeR QL\nanalogue"),
               col=c("#0072B2","grey78","white"), border=c(NA,NA,"grey40"), ylim=c(0,525), yaxt="n",
               ylab="Original 391 units", main="Support under each route's rule")
  axis(2, at=c(0,100,200,300,400), las=1)
  legend("top", legend=c("Full support","Available, not full","Not comparable"),
         fill=c("#0072B2","grey78","white"), border=c(NA,NA,"grey40"),
         horiz=TRUE, bty="n", cex=0.78, x.intersp=0.6)
  text(y, support + available_not_full + missing + 12, paste0(support,"/391 full"), cex=0.80)
  text(y[1], 455, "16 not full; 2 not comparable", cex=0.75)
  segments(y[1], 446, y[1], 427, col="grey40", lwd=0.65)
  panel_letter("c")
}

candidates <- data.frame(
  feature_id=c("PaEq_960a1016a6ff6475","PaEq_7aa4e002a1520687","PaEq_534c154a3834c4d4","PaEq_6c314d0bda685b5e","PaEq_f9c6115c5c62134b","PaEq_c045621d823d76f0"),
  family=c("ERF017-like","NHL13-like","CSLD3-like","GST-like","XTH23-like","CHS4-like"),
  maid=c("MA_81979g0010","MA_71687g0010","MA_595526g0010","MA_10426104g0010","MA_10432086g0010","MA_5735g0010"),
  edge_full=c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE), stringsAsFactors=FALSE)

draw_candidate_panel <- function(cid, fam, maid, full, letter) {
  apply_figure_style(); par(mar=c(3.4,3.65,4.25,0.45))
  d <- expr[expr$feature_id==cid,]
  groups <- c("Pa","PaPs","PaHp","PaPsHp")
  vals <- log2(d$normalized_count+1)
  ymax <- max(vals) + 0.9; ymin <- max(0, min(vals)-0.35)
  plot(NA,xlim=c(0.5,4.5),ylim=c(ymin,ymax),xaxt="n",xlab="",ylab=expression(log[2]*"(normalized count + 1)"))
  title(main=fam, line=1.25)
  axis(1,at=1:4,labels=c("Pa","PaPs","PaHp","PaPsHp"),cex.axis=0.80)
  for (i in seq_along(groups)) {
    y <- vals[d$group==groups[i]]
    x <- i+c(-0.10,0,0.10)
    points(x,y,pch=21,bg=group_cols[groups[i]],col="white",cex=0.9,lwd=0.45)
    segments(i-0.17,mean(y),i+0.17,mean(y),lwd=1.35,col="black")
  }
  mtext(maid,side=3,line=0.2,adj=0,cex=0.75,col="grey25")
  panel_letter(letter)
}

draw_fig2 <- function() {
  par(oma=c(0,0,0,0))
  layout(matrix(1:6,nrow=2,byrow=TRUE))
  for (i in seq_len(nrow(candidates))) draw_candidate_panel(candidates$feature_id[i],candidates$family[i],candidates$maid[i],candidates$edge_full[i],letters[i])
}

target_figure <- Sys.getenv("FIGURE_TARGET", unset="all")
if (target_figure %in% c("all","F1")) for (kind in c("pdf","svg","png")) { open_device("Figure1_response_route_sensitivity",6.69,5.15,kind); draw_fig1(); dev.off() }
if (target_figure %in% c("all","F2")) for (kind in c("pdf","svg","png")) { open_device("Figure2_candidate_expression",6.69,5.35,kind); draw_fig2(); dev.off() }

receipt <- data.frame(
  figure=c("Figure1_response_route_sensitivity","Figure2_candidate_expression"),
  width_mm=c(170,170), height_mm=c(130.8,135.9), dpi_png=c(300,300), pointsize_pt=c(11,11),
  source=c("primary391_*_classification.tsv; primary_391.tsv","primary391_expression_12pools.tsv"),
  stringsAsFactors=FALSE)
write.table(receipt,file.path(out,"figure_render_receipt.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
caption1 <- c(
  "Figure 1. Direction-opposite pathogen-response units and sensitivity-route status.",
  "(a) The P0/P1 plane uses the DESeq2 primary-analysis estimates from pathogen_at_DSE0 and pathogen_at_DSE1. Blue points are the 132 units that meet the method-native edgeR QL analogue full-support classifier; grey points are the remaining 259 original-sign units.",
  "(b) Direction distribution among the original 391 units. (c) Stacked route status: blue, full original-pattern support; grey, finite/available but not full support; white, missing. The allN route has 373 full, 16 finite but not full, and 2 missing units. The edgeR route uses its QL analogue rather than reconstructed DESeq2 Wald intervals."
)
writeLines(caption1, file.path(out,"Figure1_caption.txt"), useBytes=TRUE)
caption2 <- c(
  "Figure 2. Normalized RNA-seq counts are displayed as log2(count + 1) for four pooled-library groups (n=3 pooled biological samples per group).",
  "Points are individual pooled libraries; black horizontal segments are group means. Pa: neither fungus; PaPs: DSE only; PaHp: pathogen only; PaPsHp: both fungi.",
  "Panels: (a) PaEq_960a1016a6ff6475 (MA_81979g0010; ERF017-like); (b) PaEq_7aa4e002a1520687 (MA_71687g0010; NHL13-like); (c) PaEq_534c154a3834c4d4 (MA_595526g0010; CSLD3-like); (d) PaEq_6c314d0bda685b5e (MA_10426104g0010; GST-like); (e) PaEq_f9c6115c5c62134b (MA_10432086g0010; XTH23-like); (f) PaEq_c045621d823d76f0 (MA_5735g0010; CHS4-like).",
  "Family labels are source-linked homology annotations rather than host functional verification. Panels a-d meet the method-native edgeR QL full-support classifier; panels e-f preserve directions across routes but do not meet its full-support classifier."
)
writeLines(caption2, file.path(out,"Figure2_caption.txt"), useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"figure_render_sessionInfo.txt"))
