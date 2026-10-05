#!/usr/bin/env Rscript

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
# ============================================================
# Variance-partitioning visualization: predictors only (no residuals), two panels
# Each bar is labeled with the corresponding model's Adj R²
# ============================================================
library(ggplot2); library(dplyr)

build_data <- function(metric_name,
                        vlp_m1_species, vlp_m1_area,
                        vlp_m2_species, vlp_m2_area, vlp_m2_bacteria,
                        bulk_m1_species, bulk_m1_area,
                        bulk_m2_species, bulk_m2_area, bulk_m2_bacteria,
                        vlp_adjr2_m1, vlp_adjr2_m2,
                        bulk_adjr2_m1, bulk_adjr2_m2,
                        vlp_nest_p, bulk_nest_p) {

  per_ds <- c(rep("Species\n+ Area", 2), rep("Species + Area\n+ Bacteria", 3))
  df <- data.frame(
    Dataset   = rep(c("VLPs","Bulk"), each=5),
    Model     = rep(per_ds, 2),
    Predictor = c("Species","Area", "Species","Area","MAGs",
                  "Species","Area", "Species","Area","MAGs"),
    R2 = c(vlp_m1_species, vlp_m1_area,
           vlp_m2_species, vlp_m2_area, vlp_m2_bacteria,
           bulk_m1_species, bulk_m1_area,
           bulk_m2_species, bulk_m2_area, bulk_m2_bacteria),
    adj_r2 = c(rep(vlp_adjr2_m1,2), rep(vlp_adjr2_m2,3),
               rep(bulk_adjr2_m1,2), rep(bulk_adjr2_m2,3)),
    stringsAsFactors = FALSE
  )
  df$Predictor <- factor(df$Predictor, levels=c("Species","Area","MAGs"))
  df$Dataset    <- factor(df$Dataset, levels=c("VLPs","Bulk"))
  df$Model      <- factor(df$Model, levels=c("Species\n+ Area","Species + Area\n+ Bacteria"))
  df$Metric     <- metric_name
  df
}

# Exact Type II R2 values; Adj R2 from lm output
df_pi <- build_data(
  "Nucleotide diversity",
  vlp_m1_species=7.44, vlp_m1_area=36.33,
  vlp_m2_species=5.36, vlp_m2_area=36.98, vlp_m2_bacteria=0.00,
  bulk_m1_species=4.90, bulk_m1_area=20.59,
  bulk_m2_species=1.44, bulk_m2_area=20.87, bulk_m2_bacteria=29.61,
  vlp_adjr2_m1=41.77, vlp_adjr2_m2=40.35,
  bulk_adjr2_m1=18.92, bulk_adjr2_m2=48.60,
  vlp_nest_p="p = 0.96", bulk_nest_p="p < 0.0001"
)

df_sh <- build_data(
  "Shannon diversity",
  vlp_m1_species=1.84, vlp_m1_area=10.29,
  vlp_m2_species=0.41, vlp_m2_area=4.38, vlp_m2_bacteria=5.64,
  bulk_m1_species=0.29, bulk_m1_area=34.10,
  bulk_m2_species=0.64, bulk_m2_area=3.61, bulk_m2_bacteria=22.56,
  vlp_adjr2_m1=6.62, vlp_adjr2_m2=9.95,
  bulk_adjr2_m1=29.93, bulk_adjr2_m2=45.10,
  vlp_nest_p="p = 0.12", bulk_nest_p="p = 0.001"
)

colors <- c("Species"="#F0A000", "Area"="#3070C0", "MAGs"="#208050")

make_plot <- function(df) {
  metric <- df$Metric[1]

  # Per-bar summaries
  bar_info <- df %>%
    group_by(Dataset, Model) %>%
    summarise(bar_h = sum(R2), .groups="drop") %>%
    left_join(df %>% filter(!duplicated(paste(Dataset, Model))) %>%
              select(Dataset, Model, adj_r2), by=c("Dataset","Model"))

  ymax <- max(bar_info$bar_h) * 1.25

  ggplot(df, aes(x=Model, y=R2, fill=Predictor)) +
    geom_col(width=0.55, alpha=0.9) +
    geom_text(data=bar_info,
              aes(x=Model, y=bar_h + ymax*0.08,
                  label=sprintf("Adj R2 = %.1f%%", adj_r2)),
              inherit.aes=FALSE, size=3.0, fontface="italic", color="grey40") +
    scale_fill_manual(values=colors) +
    scale_y_continuous(limits=c(0, ymax), expand=c(0,0)) +
    labs(title=metric, x="", y="Variance explained (%)") +
    facet_wrap(~Dataset) +
    theme_classic(base_size=13) +
    theme(
      legend.position="bottom", legend.title=element_blank(),
      axis.text.x=element_text(size=10, face="bold", color="black"),
      axis.text.y=element_text(size=10, color="black"),
      axis.title.y=element_text(size=12),
      plot.title=element_text(size=14, face="bold", hjust=0.5),
      strip.text=element_text(size=13, face="bold"),
      strip.background=element_blank(),
      panel.grid.major.y=element_line(color="grey90", linewidth=0.2))
}

p_pi <- make_plot(df_pi)
ggsave(file.path(fig_dir, "variance_partitioning_pi.pdf"), p_pi, width=8, height=4.5, dpi=300)
ggsave(file.path(fig_dir, "variance_partitioning_pi.png"), p_pi, width=8, height=4.5, dpi=300)

p_sh <- make_plot(df_sh)
ggsave(file.path(fig_dir, "variance_partitioning_shannon.pdf"), p_sh, width=8, height=4.5, dpi=300)
ggsave(file.path(fig_dir, "variance_partitioning_shannon.png"), p_sh, width=8, height=4.5, dpi=300)

cat("Done.\n")
