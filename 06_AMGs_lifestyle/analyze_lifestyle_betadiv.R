#!/usr/bin/env Rscript
# ============================================================
# Beta-diversity: VLP AMG KO composition split by lifestyle
# Pseudo-sample design: each real sample -> temperate + virulent
# Paired PERMANOVA + PCoA with connecting lines
# ============================================================

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Working directory ──
script_path <- tryCatch(dirname(rstudioapi::getActiveDocumentContext()$path),
                        error = function(e) NULL)
if (is.null(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    script_path <- getwd()
  }
}
setwd(script_path)
data_dir <- file.path(script_path, "..", "data")
fig_dir  <- file.path(script_path, "..", "figures")
res_dir  <- file.path(script_path, "..", "results")

# ── Read data ──
cat("Reading VLP_lifestyle_KO_abundance.tsv ...\n")
df <- read.table(file.path(data_dir, "VLP_lifestyle_KO_abundance.tsv"), header = TRUE,
                 row.names = 1, sep = "\t", check.names = FALSE)

# Exclude Rat201 (zero AMG sample)
if ("Rat201" %in% colnames(df)) df <- df[, colnames(df) != "Rat201"]

n_ko  <- nrow(df)
n_smp <- ncol(df)
cat(sprintf("  %d KO-lifestyle rows x %d samples\n", n_ko, n_smp))

# ── Build pseudo-sample matrix ──
# For each real sample, create 2 columns: sample_temperate, sample_virulent
# Only KOs with BOTH lifestyles are informative for composition comparison

# Extract per-lifestyle KO abundance for each sample
ko_info <- strsplit(rownames(df), "|", fixed = TRUE)
ko_id   <- sapply(ko_info, `[`, 1)
ko_lf   <- sapply(ko_info, `[`, 2)

cat(sprintf("  KO-lifestyle breakdown: temperate=%d, virulent=%d\n",
            sum(ko_lf == "temperate"), sum(ko_lf == "virulent")))

# Create pseudo-samples: rows = KO, cols = sample_lifestyle
# Build pair matrix
temp_rows <- df[ko_lf == "temperate", , drop = FALSE]
vir_rows  <- df[ko_lf == "virulent",  , drop = FALSE]
rownames(temp_rows) <- ko_id[ko_lf == "temperate"]
rownames(vir_rows)  <- ko_id[ko_lf == "virulent"]

# Align to shared KOs
shared_ko <- intersect(rownames(temp_rows), rownames(vir_rows))
cat(sprintf("  KOs with both lifestyles: %d\n", length(shared_ko)))

temp_rows <- temp_rows[shared_ko, , drop = FALSE]
vir_rows  <- vir_rows[shared_ko, , drop = FALSE]

# Build pseudo-sample matrix: columns = sample_temperate, sample_virulent
pseudo_mat <- cbind(temp_rows, vir_rows)
colnames(pseudo_mat) <- paste0(
  rep(colnames(temp_rows), each = 1), "_",
  rep(c("temperate", "virulent"), each = n_smp)
)

# But we want: columns are all pseudo-samples, rows are KOs
# temp_rows and vir_rows already share the same KO order
pseudo_mat <- cbind(temp_rows, vir_rows)
colnames(pseudo_mat) <- c(
  paste0(colnames(temp_rows), "_temperate"),
  paste0(colnames(vir_rows),  "_virulent")
)

cat(sprintf("  Pseudo-sample matrix: %d KOs x %d pseudo-samples\n",
            nrow(pseudo_mat), ncol(pseudo_mat)))

# ── Relative abundance: each pseudo-sample sums to 1 ──
pseudo_rel <- sweep(pseudo_mat, 2, colSums(pseudo_mat), "/")
pseudo_rel[is.na(pseudo_rel)] <- 0  # guard against 0/0

# ── Bray-Curtis distance ──
d <- vegdist(t(pseudo_rel), method = "bray")
cat(sprintf("  Bray-Curtis distance: %d x %d\n", ncol(pseudo_rel), ncol(pseudo_rel)))

# ── Build metadata ──
pseudo_samples <- colnames(pseudo_rel)
meta_pseudo <- data.frame(
  pseudo_id  = pseudo_samples,
  sample     = sub("_(temperate|virulent)$", "", pseudo_samples),
  lifestyle  = ifelse(grepl("_temperate$", pseudo_samples), "temperate", "virulent"),
  row.names  = pseudo_samples
)
meta_pseudo$lifestyle <- factor(meta_pseudo$lifestyle, levels = c("virulent", "temperate"))

cat(sprintf("  Pseudo-samples: temperate=%d, virulent=%d\n",
            sum(meta_pseudo$lifestyle == "temperate"),
            sum(meta_pseudo$lifestyle == "virulent")))

# ============================================================
# PERMANOVA — paired design (strata = sample)
# ============================================================
cat("\n--- PERMANOVA (paired by sample) ---\n")
set.seed(42)
perm_res <- adonis2(d ~ lifestyle, data = meta_pseudo,
                    permutations = 9999, strata = meta_pseudo$sample)
print(perm_res)

# Also unpaired for comparison
cat("\n--- PERMANOVA (unpaired, for comparison) ---\n")
set.seed(42)
perm_unpaired <- adonis2(d ~ lifestyle, data = meta_pseudo, permutations = 9999)
print(perm_unpaired)

# ============================================================
# PCoA
# ============================================================
pcoa <- cmdscale(d, k = 5, eig = TRUE)
eig <- pcoa$eig
pct <- round(100 * eig / sum(eig), 1)

coords <- as.data.frame(pcoa$points[, 1:2])
colnames(coords) <- c("PCoA1", "PCoA2")
coords$pseudo_id <- rownames(coords)
coords$sample     <- meta_pseudo$sample
coords$lifestyle  <- factor(meta_pseudo$lifestyle, levels = c("virulent", "temperate"))

# PCoA variance explained
cat(sprintf("\nPCoA variance: Axis1=%.1f%%, Axis2=%.1f%%, first 2=%.1f%%, first 5=%.1f%%\n",
            pct[1], pct[2], sum(pct[1:2]), sum(pct[1:5])))

# ============================================================
# Plot
# ============================================================
C_TEMP <- "#ecc05d"
C_VIR  <- "#29a28c"

p <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = lifestyle)) +
  # 95% confidence ellipse
  stat_ellipse(aes(fill = lifestyle), level = 0.95, alpha = 0.12,
               geom = "polygon", linewidth = 0.4) +
  # Points
  geom_point(aes(fill = lifestyle), size = 2.8, shape = 21,
             color = "grey30", stroke = 0.3, alpha = 0.85) +
  scale_color_manual(values = c(C_VIR, C_TEMP), breaks = c("virulent", "temperate")) +
  scale_fill_manual(values = c(C_VIR, C_TEMP), breaks = c("virulent", "temperate")) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", pct[1]),
    y = sprintf("PCoA2 (%.1f%%)", pct[2])
  ) +
  annotate("text", x = -Inf, y = -Inf, hjust = 0, vjust = -0.5,
           label = sprintf("PERMANOVA: R² = %.3f, p = %.4f",
                           perm_res$R2[1], perm_res$`Pr(>F)`[1]),
           size = 4, fontface = "italic", color = "grey30") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.88, 0.92),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.background = element_rect(color = "grey85", linewidth = 0.3),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13)
  )

ggsave(file.path(fig_dir, "VLP_lifestyle_KO_betadiv_pcoa.pdf"), p, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "VLP_lifestyle_KO_betadiv_pcoa.png"), p, width = 7.5, height = 5.5, dpi = 300)

cat("\nPlots saved: VLP_lifestyle_KO_betadiv_pcoa.pdf/.png\n")

# ============================================================
# Summary stats
# ============================================================
cat("\n--- Summary ---\n")
all_dist <- as.matrix(d)
# Mean within-pair distance (same sample, temperate vs virulent)
pair_dists <- c()
for (s in unique(meta_pseudo$sample)) {
  idx_t <- which(meta_pseudo$sample == s & meta_pseudo$lifestyle == "temperate")
  idx_v <- which(meta_pseudo$sample == s & meta_pseudo$lifestyle == "virulent")
  if (length(idx_t) == 1 && length(idx_v) == 1) {
    pair_dists <- c(pair_dists, all_dist[idx_t, idx_v])
  }
}
cat(sprintf("Mean within-pair (same sample, temp vs vir): %.4f\n",
            mean(pair_dists)))

# Between-lifestyle distances
t_idx <- which(meta_pseudo$lifestyle == "temperate")
v_idx <- which(meta_pseudo$lifestyle == "virulent")
cross_temp <- all_dist[t_idx, t_idx]; diag(cross_temp) <- NA
cross_vir  <- all_dist[v_idx, v_idx];  diag(cross_vir)  <- NA
cross_tv   <- all_dist[t_idx, v_idx]

cat(sprintf("Mean temperate-temperate: %.4f\n", mean(cross_temp, na.rm = TRUE)))
cat(sprintf("Mean virulent-virulent:   %.4f\n", mean(cross_vir,  na.rm = TRUE)))
cat(sprintf("Mean temperate-virulent:  %.4f\n", mean(cross_tv)))

cat("\nDone!\n")
