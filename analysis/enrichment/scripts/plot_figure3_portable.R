#!/usr/bin/env Rscript
# Portable Figure 3 renderer. It preserves the frozen Figure 3 data rule and
# visual specification, but requires an explicit project root and a new output
# directory. It never writes to the original figure directory.
#
# Usage:
#   Rscript --vanilla plot_figure3_portable.R /absolute/project_root /new/output_dir

suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(grid) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla plot_figure3_portable.R /absolute/project_root /new/output_dir")
project_root <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- normalizePath(args[[2]], mustWork = FALSE)
input_file <- file.path(project_root, "run", "results", "02_primary_all_terms.tsv")
if (!file.exists(input_file)) stop("Missing frozen Figure 3 source: ", input_file)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

apply_figure_style <- function(base_size = 9, legend_size = 8, tick_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(size = base_size, face = "plain", hjust = 0, margin = margin(b = 3)),
      plot.subtitle = element_text(size = legend_size, hjust = 0, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = tick_size, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(2.3, "pt"),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      legend.position = "right",
      legend.title = element_text(size = legend_size, face = "plain"),
      legend.text = element_text(size = tick_size),
      legend.key.height = unit(4.2, "mm"),
      legend.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 9, 5.5, 5.5)
    )
}

raw <- fread(input_file, sep = "\t", data.table = FALSE)
need <- subset(raw,
  analysis == "A_overall391_vs_annotated_background" &
    reportable_under_fixed_rule == "TRUE" &
    q_combined < 0.05
)
stopifnot(nrow(need) == 3L, all(need$ci_low > 0), all(need$ci_high > 0), all(need$odds_ratio > 1))
need$minus_log10_q <- -log10(need$q_combined)
need$label <- paste0(need$name, "\n", need$term_id)
need <- need[order(need$odds_ratio, decreasing = FALSE), ]
need$label <- factor(need$label, levels = need$label)
fwrite(need, file.path(figure_dir, "Figure3_source_data.tsv"), sep = "\t")

fig_width_mm <- 170; fig_height_mm <- 90
fig_width_in <- fig_width_mm / 25.4; fig_height_in <- fig_height_mm / 25.4
p <- ggplot(need, aes(x = odds_ratio, y = label)) +
  geom_vline(xintercept = 1, linewidth = 0.45, linetype = "dashed", colour = "#666666") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.13, linewidth = 0.7, colour = "#34495E", orientation = "y") +
  geom_point(aes(size = selected_term_n, colour = minus_log10_q), alpha = 0.98) +
  scale_colour_gradient(low = "#73A8D6", high = "#0B3C6F", name = expression(-log[10](q))) +
  scale_size_area(max_size = 11, breaks = sort(unique(need$selected_term_n)), name = "Selected units\nwith term") +
  scale_x_continuous(name = "Odds ratio (95% exact CI)", breaks = c(0, 1, 2, 3, 4),
                     limits = c(0, max(need$ci_high) * 1.10), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Predicted GO-BP enrichment", y = NULL) +
  scale_y_discrete(expand = expansion(add = c(0.30, 0.68))) +
  annotate("text", x = 1, y = 3.42, label = "OR = 1", vjust = 0, size = 2.8, colour = "#4D4D4D") +
  apply_figure_style()

ggsave(file.path(figure_dir, "Figure3_GO_enrichment.png"), p, width = fig_width_in, height = fig_height_in, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "Figure3_GO_enrichment.pdf"), p, width = fig_width_in, height = fig_height_in, units = "in", device = cairo_pdf, bg = "white")
svg(file.path(figure_dir, "Figure3_GO_enrichment.svg"), width = fig_width_in, height = fig_height_in, pointsize = 9, bg = "white")
print(p)
dev.off()

qa <- c(
  "Figure 3 portable rendering QA",
  paste("Input:", input_file),
  paste("Output:", figure_dir),
  "Eligibility: analysis A; reportable_under_fixed_rule=TRUE; q_combined<0.05.",
  paste("Plotted terms:", nrow(need)),
  paste("Physical dimensions:", fig_width_mm, "x", fig_height_mm, "mm; PNG 300 dpi."),
  "Visual encodings: horizontal position=odds ratio; whisker=two-sided exact 95% Fisher CI; area=selected unit count; colour=-log10(combined Fisher q).",
  "Reference: dashed vertical line marks OR=1; it is not a pathway-activity threshold.",
  "No DAG, gene network, or mechanistic edge is drawn.",
  "Caption must state that term memberships overlap and enrichment does not measure pathway activity."
)
writeLines(qa, file.path(figure_dir, "Figure3_QA.txt"))
