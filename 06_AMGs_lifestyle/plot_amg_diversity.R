#!/usr/bin/env Rscript
# AMG Diversity boxplots: VLP + Bulk, Richness + Shannon (paired Wilcoxon)
# virulent vs temperate phage AMGs

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

script_path <- tryCatch(dirname(rstudioapi::getActiveDocumentContext()$path), error = function(e) NULL)
if (is.null(script_path) || script_path == "") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
}
setwd(script_path)
data_dir <- file.path(script_path, "..", "data")
fig_dir  <- file.path(script_path, "..", "figures")
res_dir  <- file.path(script_path, "..", "results")

make_plot <- function(df, metric, dataset_label, outfile) {
  # Pivot to long
  vir_col <- paste0("virulent_", metric)
  temp_col <- paste0("temperate_", metric)

  plot_data <- df %>%
    select(sample, !!vir_col, !!temp_col) %>%
    rename(virulent = !!vir_col, temperate = !!temp_col) %>%
    pivot_longer(c(virulent, temperate), names_to = "lifestyle", values_to = "value")

  # Paired Wilcoxon
  p <- wilcox.test(df[[vir_col]], df[[temp_col]], paired = TRUE)$p.value
  p_text <- ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, 4)))

  # Line data
  lines <- df %>% select(sample, virulent = !!vir_col, temperate = !!temp_col)

  y_label <- ifelse(metric == "richness", "Richness (unique KOs)", "Shannon diversity")

  ggplot(plot_data, aes(x = lifestyle, y = value)) +
    geom_boxplot(aes(fill = lifestyle), alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_point(aes(color = lifestyle), size = 2.5, alpha = 0.8,
               position = position_jitter(width = 0.08, seed = 42)) +
    scale_fill_manual(values = c("virulent" = "#129880", "temperate" = "#EBBA4C"), guide = "none") +
    scale_color_manual(values = c("virulent" = "#129880", "temperate" = "#EBBA4C"), guide = "none") +
    annotate("text", x = 1.5, y = max(plot_data$value) * 0.95, label = p_text, size = 5, fontface = "bold") +
    labs(x = NULL, y = y_label, title = paste(dataset_label, "-", metric)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 13),
      panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.25)
    )

  ggsave(file.path(fig_dir, outfile), width = 5, height = 5)
  cat(sprintf("[%s] %s: p=%.6f\n", dataset_label, metric, p))
}

# VLP
vlp <- read_tsv(file.path(data_dir, "53vlp_AMG_diversity.tsv"), show_col_types = FALSE)
make_plot(vlp, "richness", "VLP", "AMG_diversity_VLP_richness.pdf")
make_plot(vlp, "shannon", "VLP", "AMG_diversity_VLP_shannon.pdf")

# Bulk
bulk <- read_tsv(file.path(data_dir, "53meta_AMG_diversity.tsv"), show_col_types = FALSE)
make_plot(bulk, "richness", "Bulk", "AMG_diversity_Bulk_richness.pdf")
make_plot(bulk, "shannon", "Bulk", "AMG_diversity_Bulk_shannon.pdf")

cat("[DONE]\n")
