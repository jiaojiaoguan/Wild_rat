#!/usr/bin/env Rscript
# ============================================================
# Shannon drivers analysis — Type II ANOVA (consistent with manuscript)
# Base model:  Shannon ~ Species + Area
# + MAG model: Shannon ~ Species + Area + Bac_Shannon
#   (Bac_Shannon = Shannon diversity of bacterial MAGs, from explanatory_variables.csv)
# Excludes R. andamanensis; keeps only samples shared with bacterial data (n=45)
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

meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types = FALSE)
meta <- meta %>% mutate(
  Species = dplyr::recode(Species,
    "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
    "Rattus_norvegicus" = "R. norvegicus", "Rattus_tanezumi" = "R. tanezumi"))
meta <- meta[meta$Species != "R. andamanensis", ]

# Shannon diversity of bacterial MAGs (the "Shannon" column in the explanatory-variable table)
shared_X <- read.csv(file.path(data_dir, "explanatory_variables.csv"), row.names = 1)

for (lib in c("VLPs", "Bulk")) {
  div_file <- if (lib == "VLPs") file.path(data_dir, "VLPs_virus_alpha_diversity.tsv")
              else                  file.path(data_dir, "bulk_virus_alpha_diversity.tsv")
  div <- read_tsv(div_file, show_col_types = FALSE)
  colnames(div) <- c("sample", "Richness", "Shannon")
  df <- inner_join(div, meta, by = c("sample" = "sample_name"))
  df <- cbind(df, Bac_Shannon = shared_X[df$sample, "Shannon"])
  df <- df[complete.cases(df), ]

  cat(sprintf("\n========================================\n"))
  cat(sprintf("  %s (n=%d)\n", lib, nrow(df)))
  cat(sprintf("========================================\n"))

  # ---- Base model: Species + Area (Type II) ----
  m1 <- lm(Shannon ~ Species + Area, data = df)
  aov1 <- car::Anova(m1, type = 2)
  ss1 <- sum(aov1$`Sum Sq`)

  cat("\n--- Base model: Shannon ~ Species + Area (Type II) ---\n")
  cat(sprintf("Adj R² = %.2f%%\n", summary(m1)$adj.r.squared * 100))
  for (i in 1:nrow(aov1)) {
    nm <- rownames(aov1)[i]
    cat(sprintf("  %-15s SS=%.4f  R²=%.2f%%  F=%.3f  p=%.4f\n",
        nm, aov1$`Sum Sq`[i], aov1$`Sum Sq`[i] / ss1 * 100,
        aov1$`F value`[i], aov1$`Pr(>F)`[i]))
  }

  # ---- + MAG Shannon: Species + Area + Bac_Shannon (Type II) ----
  m2 <- lm(Shannon ~ Species + Area + Bac_Shannon, data = df)
  aov2 <- car::Anova(m2, type = 2)
  ss2 <- sum(aov2$`Sum Sq`)

  cat("\n--- + MAG Shannon: Shannon ~ Species + Area + Bac_Shannon (Type II) ---\n")
  cat(sprintf("Adj R² = %.2f%%\n", summary(m2)$adj.r.squared * 100))
  for (i in 1:nrow(aov2)) {
    nm <- rownames(aov2)[i]
    cat(sprintf("  %-15s SS=%.4f  R²=%.2f%%  F=%.3f  p=%.4f\n",
        nm, aov2$`Sum Sq`[i], aov2$`Sum Sq`[i] / ss2 * 100,
        aov2$`F value`[i], aov2$`Pr(>F)`[i]))
  }

  # ---- Nested F-test (base vs +MAG Shannon) ----
  cat("\n--- Nested F-test (base vs +MAG Shannon) ---\n")
  comp <- anova(m1, m2)
  print(comp)
  cat(sprintf("ΔAdj R² = %.2f%%\n",
              (summary(m2)$adj.r.squared - summary(m1)$adj.r.squared) * 100))

  # ---- R² comparison table ----
  cat("\n--- R² comparison ---\n")
  cat(sprintf("%-15s %12s %12s\n", "Factor", "Base R²", "+MAG R²"))
  cat(strrep("-", 42), "\n")
  cat(sprintf("%-15s %11.2f%% %11.2f%%\n", "Species",
      aov1["Species", "Sum Sq"] / ss1 * 100, aov2["Species", "Sum Sq"] / ss2 * 100))
  cat(sprintf("%-15s %11.2f%% %11.2f%%\n", "Area",
      aov1["Area", "Sum Sq"] / ss1 * 100, aov2["Area", "Sum Sq"] / ss2 * 100))
  cat(sprintf("%-15s %11s %11.2f%%\n", "Bac_Shannon", "—",
      aov2["Bac_Shannon", "Sum Sq"] / ss2 * 100))
  cat(sprintf("%-15s %11.2f%% %11.2f%%\n", "Residuals",
      100 - sum(aov1$`Sum Sq`[1:(nrow(aov1) - 1)]) / ss1 * 100,
      100 - sum(aov2$`Sum Sq`[1:(nrow(aov2) - 1)]) / ss2 * 100))
}
cat("\nDone.\n")

