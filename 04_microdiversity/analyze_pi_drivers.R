#!/usr/bin/env Rscript
# ============================================================
# Viral π (nucleotide diversity) drivers analysis — Type II nested linear models (consistent with manuscript)
# Base model:  avg_pi ~ Species + Area
# + MAG model: avg_pi ~ Species + Area + MAG_pi
#   (MAG_pi = mean nucleotide diversity of bacterial MAGs, from MAGs_nucl_diversity_matrix.tsv)
# Excludes R. andamanensis (n=46)
# ============================================================

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

library(readr); library(dplyr); library(car)

# Metadata (exclude andamanensis)
meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types = FALSE)
meta <- meta %>% mutate(
  Species = dplyr::recode(Species,
    "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
    "Rattus_norvegicus" = "R. norvegicus", "Rattus_tanezumi" = "R. tanezumi"))
meta <- meta %>% filter(Species != "R. andamanensis")

# Bacterial MAG π: mean π per sample (ignoring NA)
mag <- read_tsv(file.path(data_dir, "MAGs_nucl_diversity_matrix.tsv"), show_col_types = FALSE)
mag_mean <- colMeans(mag[, -1], na.rm = TRUE)
mag_df <- data.frame(sample = names(mag_mean), MAG_pi = as.numeric(mag_mean))

for (lib in c("VLPs", "Bulk")) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("  %s: avg_pi ~ Species + Area + MAG_pi\n", lib))
  cat(sprintf("============================================================\n"))

  csv <- if (lib == "VLPs") file.path(data_dir, "VLPs_global_sample_microdiversity.csv")
         else                  file.path(data_dir, "bulk_global_sample_microdiversity.csv")
  virus <- read_csv(csv, show_col_types = FALSE)

  df <- virus %>%
    inner_join(meta, by = c("sample" = "sample_name")) %>%
    inner_join(mag_df, by = "sample") %>%
    filter(complete.cases(.))
  cat(sprintf("Final n = %d\n", nrow(df)))

  # ---- Model 1: Species + Area ----
  m1 <- lm(avg_pi ~ Species + Area, data = df)
  aov1 <- car::Anova(m1, type = 2)
  ss1 <- sum(aov1$`Sum Sq`)
  r2_1 <- summary(m1)$adj.r.squared

  cat("\n--- Model 1: avg_pi ~ Species + Area ---\n")
  cat(sprintf("Adj R² = %.2f%%\n", r2_1 * 100))
  for (i in 1:nrow(aov1)) {
    nm <- rownames(aov1)[i]
    cat(sprintf("  %-12s  R²=%.2f%%  F=%.3f  p=%.4f\n",
        nm, aov1$`Sum Sq`[i]/ss1*100, aov1$`F value`[i], aov1$`Pr(>F)`[i]))
  }

  # ---- Model 2: Species + Area + MAG_pi ----
  m2 <- lm(avg_pi ~ Species + Area + MAG_pi, data = df)
  aov2 <- car::Anova(m2, type = 2)
  ss2 <- sum(aov2$`Sum Sq`)
  r2_2 <- summary(m2)$adj.r.squared

  cat("\n--- Model 2: avg_pi ~ Species + Area + MAG_pi ---\n")
  cat(sprintf("Adj R² = %.2f%%\n", r2_2 * 100))
  for (i in 1:nrow(aov2)) {
    nm <- rownames(aov2)[i]
    cat(sprintf("  %-12s  R²=%.2f%%  F=%.3f  p=%.4f\n",
        nm, aov2$`Sum Sq`[i]/ss2*100, aov2$`F value`[i], aov2$`Pr(>F)`[i]))
  }

  # ---- Model comparison (nested F-test) ----
  cat("\n--- Model comparison (nested F-test) ---\n")
  comp <- anova(m1, m2)
  print(comp)
  cat(sprintf("ΔAdj R² (MAG_pi added) = %.2f%%\n", (r2_2 - r2_1) * 100))

  # ---- R² change comparison table ----
  cat("\n--- R² comparison ---\n")
  cat(sprintf("%-15s %12s %12s %12s\n", "Factor", "M1 R²", "M2 R²", "Change"))
  cat(strrep("-", 55), "\n")
  cat(sprintf("%-15s %11.2f%% %11.2f%% %11.2f%%\n", "Species",
      aov1["Species","Sum Sq"]/ss1*100, aov2["Species","Sum Sq"]/ss2*100,
      (aov2["Species","Sum Sq"]/ss2 - aov1["Species","Sum Sq"]/ss1)*100))
  cat(sprintf("%-15s %11.2f%% %11.2f%% %11.2f%%\n", "Area",
      aov1["Area","Sum Sq"]/ss1*100, aov2["Area","Sum Sq"]/ss2*100,
      (aov2["Area","Sum Sq"]/ss2 - aov1["Area","Sum Sq"]/ss1)*100))
  cat(sprintf("%-15s %11s %11.2f%% %11s\n", "MAG_pi", "—",
      aov2["MAG_pi","Sum Sq"]/ss2*100, "—"))
  cat(sprintf("%-15s %11.2f%% %11.2f%%\n", "Residuals",
      100 - sum(aov1$`Sum Sq`[1:(nrow(aov1)-1)])/ss1*100,
      100 - sum(aov2$`Sum Sq`[1:(nrow(aov2)-1)])/ss2*100))
}

cat("\nDone.\n")
