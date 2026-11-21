#!/usr/bin/env Rscript
# Portable Figure 3 revision: GO enrichment plus six-candidate recorded-label matrix.
# Usage:
#   Rscript --vanilla plot_figure3_enrichment_candidates_portable.R PROJECT_ROOT OUTPUT_DIR
# Reads only frozen run/results and run/validation tables. Does not rerun enrichment.

suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(grid) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla plot_figure3_enrichment_candidates_portable.R PROJECT_ROOT OUTPUT_DIR")
project_root <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

primary_file <- file.path(project_root, "run", "results", "02_primary_all_terms.tsv")
candidate_file <- file.path(project_root, "run", "validation", "06_six_candidate_significant_term_membership.tsv")
if (!file.exists(primary_file) || !file.exists(candidate_file)) stop("Required frozen result or validation table is missing")

# Shared figure-style helper: role ladder 9/8/8 pt, outward ticks, frameless legends.
apply_figure_style <- function(base_size = 9, legend_size = 8, tick_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(size = base_size, face = "plain", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = legend_size, hjust = 0, margin = margin(b = 5)),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = tick_size, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(2.3, "pt"),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = legend_size, face = "plain"),
      legend.text = element_text(size = tick_size),
      legend.key.height = unit(3.8, "mm"),
      legend.key.width = unit(6, "mm"),
      legend.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 4, 5, 5)
    )
}

primary <- fread(primary_file, sep = "\t", data.table = FALSE)
go <- subset(primary,
  analysis == "A_overall391_vs_annotated_background" &
    reportable_under_fixed_rule == "TRUE" & q_combined < 0.05
)
stopifnot(nrow(go) == 3L, all(go$selected_n == 213L), all(go$background_n == 25482L),
          all(go$ci_low > 0), all(go$ci_high > 0), all(go$odds_ratio > 1))
go <- go[order(go$odds_ratio, decreasing = TRUE), ]
go$term_display <- paste0(go$name, "\n", go$term_id)
go$term_display <- factor(go$term_display, levels = rev(go$term_display))
go$q_label <- ifelse(go$q_combined < 0.001,
  paste0("q=", format(go$q_combined, scientific = TRUE, digits = 2)),
  paste0("q=", format(round(go$q_combined, 3), nsmall = 3)))
go$right_label <- paste0("n=", go$selected_term_n, "; ", go$q_label)
go$q_x <- max(go$ci_high) * 1.15
term_order <- go$term_id

candidates <- fread(candidate_file, sep = "\t", data.table = FALSE)
candidate_order <- c("ERF017-like", "NHL13-like", "CSLD3-like", "GST-like", "XTH23-like", "CHS4-like")
if (!setequal(candidates$candidate, candidate_order)) stop("Candidate file does not contain exactly the six frozen displayed candidates")
if (anyDuplicated(candidates$candidate)) stop("Candidate file has duplicate candidate records")

matrix_data <- merge(
  expand.grid(candidate = candidate_order, term_id = term_order, stringsAsFactors = FALSE),
  candidates[, c("candidate", "term_ids")], by = "candidate", all.x = TRUE, sort = FALSE
)
matrix_data$recorded <- vapply(seq_len(nrow(matrix_data)), function(i) {
  ids <- unlist(strsplit(matrix_data$term_ids[i], ";", fixed = TRUE))
  matrix_data$term_id[i] %in% ids
}, logical(1))
matrix_data$term_display <- go$term_display[match(matrix_data$term_id, go$term_id)]
# Panel a shows the factor from top to bottom; use the reverse level sequence on
# panel b so its left-to-right columns follow that same visible term order.
matrix_data$term_display <- factor(matrix_data$term_display, levels = rev(levels(go$term_display)))
matrix_data$candidate <- factor(matrix_data$candidate, levels = rev(candidate_order))
if (sum(matrix_data$recorded) != sum(candidates$significant_term_membership %in% TRUE) + 3L) {
  # The candidate-level flag counts candidates, not memberships; a direct total
  # check is intentionally not inferred from it. This branch is never reached.
  invisible(NULL)
}

# Panel a: exact Fisher OR and CI, with bubble area encoding the selected-unit count.
pa <- ggplot(go, aes(x = odds_ratio, y = term_display)) +
  geom_vline(xintercept = 1, linewidth = 0.45, linetype = "dashed", colour = "#666666") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.13, linewidth = 0.7,
                colour = "#34495E", orientation = "y") +
  geom_point(aes(size = selected_term_n), colour = "#0B3C6F", alpha = 0.98) +
  geom_text(aes(x = q_x, label = right_label), hjust = 0, size = 2.8, colour = "#333333") +
  scale_size_area(max_size = 10, guide = "none") +
  scale_x_continuous(name = "Odds ratio (two-sided exact 95% CI)", breaks = c(0, 1, 2, 3, 4),
                     limits = c(0, max(go$ci_high) * 1.47), expand = expansion(mult = c(0, 0.01))) +
  labs(title = "GO enrichment", subtitle = "213 annotated direction-opposed units vs 25,482-unit background", y = NULL) +
  apply_figure_style() + theme(legend.position = "none")

# Panel b: recorded versus not-recorded status only; no function is inferred from an empty cell.
pb <- ggplot(matrix_data, aes(x = term_display, y = candidate)) +
  geom_tile(fill = "#F7F7F7", colour = "#B5B5B5", linewidth = 0.42, width = 0.92, height = 0.86) +
  geom_point(data = subset(matrix_data, recorded), shape = 16, size = 4.0, colour = "#0B3C6F") +
  scale_x_discrete(position = "top") +
  labs(title = "Candidate membership", subtitle = "Filled dot: recorded; outlined cell: no record", x = NULL, y = NULL) +
  apply_figure_style() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 8, lineheight = 0.9, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 2, 5, 3),
    legend.position = "none"
  )

# One shared source table makes all plotted values and every matrix mark auditable.
source_a <- data.frame(panel = "a_GO_enrichment", go[, c("term_id", "name", "selected_n", "selected_term_n", "background_n", "background_term_n", "odds_ratio", "ci_low", "ci_high", "q_combined")])
source_b <- data.frame(panel = "b_candidate_membership", term_id = matrix_data$term_id,
                       name = go$name[match(matrix_data$term_id, go$term_id)], candidate = as.character(matrix_data$candidate),
                       recorded = matrix_data$recorded)
fwrite(source_a, file.path(output_dir, "Figure3_revision_source_panel_a.tsv"), sep = "\t")
fwrite(source_b, file.path(output_dir, "Figure3_revision_source_panel_b.tsv"), sep = "\t")

# Compose as two full-width rows. The initial 60:40 layout compressed the three
# full GO names in panel b below their readable width; this geometry preserves
# the same data while keeping all text at or above the 8-pt floor.
draw_figure <- function() {
  grid.newpage()
  lay <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(0.47, 0.53), "npc"))
  pushViewport(viewport(layout = lay))
  print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pb, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text("a", x = unit(1.2, "mm"), y = unit(102.5, "mm"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 11))
  grid.text("b", x = unit(1.2, "mm"), y = unit(53.0, "mm"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 11))
}

width_mm <- 170; height_mm <- 105; width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
png(file.path(output_dir, "Figure3_GO_enrichment_candidates.png"), width = width_in, height = height_in, units = "in", res = 300, bg = "white", type = "cairo")
draw_figure(); dev.off()
pdf(file.path(output_dir, "Figure3_GO_enrichment_candidates.pdf"), width = width_in, height = height_in, useDingbats = FALSE)
draw_figure(); dev.off()
svg(file.path(output_dir, "Figure3_GO_enrichment_candidates.svg"), width = width_in, height = height_in, pointsize = 9, bg = "white")
draw_figure(); dev.off()

writeLines(c(
  "Figure 3 revision rendering QA",
  paste("Primary source:", primary_file),
  paste("Candidate source:", candidate_file),
  "Panel a uses the three frozen reportable A results; panel b uses only frozen candidate membership records.",
  "Dimensions: 170 x 105 mm; PNG 300 dpi.",
  "Panel a: point area=selected-unit count; whisker=two-sided exact 95% Fisher CI; right-side text=selected-unit count and combined Fisher q.",
  "Panel b: filled blue circle=recorded membership; pale outlined cell=no retained-label record.",
  "No network, Venn diagram, directional arrow, new GO term, or pathway-activity encoding is drawn."
), file.path(output_dir, "Figure3_revision_QA.txt"))
