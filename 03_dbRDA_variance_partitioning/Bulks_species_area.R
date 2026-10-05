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

#############################################################################
# Model 0: Partial db-RDA — Bulk virus ~ Species + Area (no bacteria)
# Bray-Curtis distance, varpart() for 2-way decomposition
#############################################################################

library(vegan)
library(readr)
library(dplyr)

set.seed(42)

# ============================================================
# 0. Paths
# ============================================================
out_dir  <- file.path(res_dir, "Bulk_RDA_results_species_area")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sink(file.path(out_dir, "analysis_log.txt"), split = TRUE)
cat("=== Partial RDA: Bulk virus ~ Species + Area ===\n")
cat("Started:", date(), "\n\n")

# ============================================================
# 1. Read data (SAME pipeline as Model 1 / Model 2)
# ============================================================
cat("--- 1. Reading data ---\n")

# Virus table (already normalized version, as in Model 1)
df <- read_tsv(file.path(data_dir,
        "53bulk_virus_map_bulk_sample_RPKM_normalized_data.tsv"),
        show_col_types = FALSE)
feats <- df[[1]]; df[[1]] <- NULL
bulk_mat <- t(as.matrix(df)); colnames(bulk_mat) <- feats
cat("Bulk virus:", nrow(bulk_mat), "samples x", ncol(bulk_mat), "vOTUs\n")

# Metadata
meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"),
                 show_col_types = FALSE)
meta$Species <- dplyr::recode(meta$Species,
  "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
  "Rattus_norvegicus" = "R. norvegicus",
  "Rattus_tanezumi"   = "R. tanezumi")

# Shared sample list (use the same samples as the bacteria-included model
# so all three models are directly comparable)
shared_X <- read.csv(file.path(data_dir,
        "explanatory_variables.csv"),
        row.names = 1)
cat("Shared X:", nrow(shared_X), "samples, cols:",
    paste(colnames(shared_X), collapse = ", "), "\n")

# ============================================================
# 2. Align samples (exclude andamanensis + Rat201)
# ============================================================
cat("\n--- 2. Aligning samples ---\n")
andamanensis_samples <- meta$sample_name[meta$Species == "R. andamanensis"]

common <- Reduce(intersect, list(rownames(bulk_mat),
                                 rownames(shared_X),
                                 meta$sample_name))
common <- setdiff(common, c("Rat201", andamanensis_samples))
cat("Common samples after exclusion:", length(common), "\n")

bulk_mat <- bulk_mat[common, , drop = FALSE]
meta_sub <- meta[match(common, meta$sample_name), ]

# ============================================================
# 3. Build explanatory table X
# ============================================================
X <- data.frame(
  Species = factor(meta_sub$Species),
  Area    = factor(meta_sub$Area),
  row.names = common
)
cat("\nX dimensions:", nrow(X), "x", ncol(X), "\n")
cat("Species table:\n"); print(table(X$Species))
cat("Area table:\n");    print(table(X$Area))
write.csv(X, file.path(out_dir, "explanatory_variables.csv"))

# ============================================================
# 4. Relative abundance + Bray-Curtis (same as Model 1/2)
# ============================================================
cat("\n--- 4. Bray-Curtis distance ---\n")
bulk_rel <- decostand(bulk_mat, method = "total")
d <- vegdist(bulk_rel, method = "bray")
cat("Distance range:", round(range(d), 4), "\n")

# ============================================================
# 5. Full model dbRDA
# ============================================================
cat("\n--- 5. Full dbrda model ---\n")
rda_full <- dbrda(d ~ Species + Area, data = X)
print(rda_full)
adj_full <- RsquareAdj(rda_full)$adj.r.squared
cat(sprintf("\nFull model AdjR2: %.4f (%.2f%%)\n", adj_full, adj_full*100))

# ============================================================
# 6. ANOVA: overall / terms / margin / axis
# ============================================================
cat("\n--- 6a. Overall ANOVA ---\n")
set.seed(42); aov_overall <- anova(rda_full, permutations = 9999)
print(aov_overall)
write.csv(as.data.frame(aov_overall), file.path(out_dir, "anova_overall.csv"))

cat("\n--- 6b. ANOVA by terms (sequential) ---\n")
set.seed(42); aov_terms <- anova(rda_full, by = "terms", permutations = 9999)
print(aov_terms)
write.csv(as.data.frame(aov_terms), file.path(out_dir, "anova_terms.csv"))

cat("\n--- 6c. ANOVA by margin (Type III) ---\n")
set.seed(42); aov_margin <- anova(rda_full, by = "margin", permutations = 9999)
print(aov_margin)
write.csv(as.data.frame(aov_margin), file.path(out_dir, "anova_margin.csv"))

cat("\n--- 6d. ANOVA by axis ---\n")
set.seed(42); aov_axis <- anova(rda_full, by = "axis", permutations = 9999)
print(aov_axis)
write.csv(as.data.frame(aov_axis), file.path(out_dir, "anova_axis.csv"))

# ============================================================
# 7. Partial RDA — Drop-one nested tests (manual marginal)
# ============================================================
cat("\n--- 7. Drop-one nested tests ---\n")

variables <- c("Species", "Area")
partial_results <- data.frame()

for (var in variables) {
  remaining <- setdiff(variables, var)
  fml <- as.formula(paste("d ~", paste(remaining, collapse = " + ")))
  rda_red <- dbrda(fml, data = X)
  adj_red <- RsquareAdj(rda_red)$adj.r.squared
  delta   <- adj_full - adj_red

  set.seed(42)
  nested  <- anova(rda_red, rda_full, permutations = 9999)
  f_val   <- nested$F[2]
  p_val   <- nested$`Pr(>F)`[2]

  partial_results <- rbind(partial_results, data.frame(
    Variable   = var,
    AdjR2_red  = round(adj_red, 6),
    DeltaAdjR2 = round(delta, 6),
    F_stat     = round(f_val, 4),
    p_value    = round(p_val, 4)
  ))
  cat(sprintf("Drop %-10s  AdjR2_red=%.4f  Delta=%.4f  F=%.3f  p=%.4f\n",
              var, adj_red, delta, f_val, p_val))
}
write.csv(partial_results,
          file.path(out_dir, "partial_RDA_marginal_effects.csv"),
          row.names = FALSE)

# ============================================================
# 8. Variance Partitioning — 2 groups via varpart()
#    X1 = Species, X2 = Area
# ============================================================
cat("\n--- 8. Variance Partitioning (varpart) ---\n")

vp <- varpart(d, ~ Species, ~ Area, data = X)
print(vp)

ind <- vp$part$indfract
adj_col <- grep("Adj", colnames(ind), value = TRUE)[1]
cat("\nFractions (rows):", paste(rownames(ind), collapse = " | "), "\n")
cat("Using column:", adj_col, "\n\n")

# For 2 explanatory tables varpart returns 4 rows in fixed order:
#  1 [a]   X1|X2          (pure Species)
#  2 [b]   X2|X1          (pure Area)
#  3 [c]   shared         (shared Species & Area)
#  4 [d]   Residuals
stopifnot(nrow(ind) == 4)
vals <- as.numeric(ind[, adj_col])
names(vals) <- c("[a] Pure Species",
                 "[b] Pure Area",
                 "[c] Shared Sp∩Ar",
                 "[d] Residuals")

explained <- sum(vals[1:3], na.rm = TRUE)

cat("=== Variance Partitioning (Adj R^2) ===\n")
for (nm in names(vals)) {
  cat(sprintf("  %-25s %8.4f  (%6.2f%%)\n", nm, vals[nm], vals[nm]*100))
}
cat(sprintf("  %-25s %8.4f  (%6.2f%%)\n", "TOTAL EXPLAINED", explained, explained*100))

vpa_all <- data.frame(
  Fraction = names(vals),
  AdjR2    = round(vals, 6),
  Pct      = round(vals * 100, 2),
  row.names = NULL
)
vpa_all <- rbind(vpa_all,
                 data.frame(Fraction = "SUM Explained [a..c]",
                            AdjR2 = round(explained, 6),
                            Pct   = round(explained * 100, 2)))
write.csv(vpa_all, file.path(out_dir, "vpa_fractions.csv"), row.names = FALSE)

# ============================================================
# 9. Significance tests for pure (testable) fractions
# ============================================================
cat("\n--- 9. Significance of pure (testable) fractions ---\n")

# Pure Species | Area
set.seed(42); t_Sp <- anova(dbrda(d ~ Species + Condition(Area), data = X),
                            permutations = 9999)
# Pure Area | Species
set.seed(42); t_Ar <- anova(dbrda(d ~ Area + Condition(Species), data = X),
                            permutations = 9999)

cat("\n[a] Pure Species | Area:\n"); print(t_Sp)
cat("\n[b] Pure Area | Species:\n"); print(t_Ar)

pure_tests <- data.frame(
  Fraction = c("[a] Pure Species", "[b] Pure Area"),
  F_stat   = c(t_Sp$F[1], t_Ar$F[1]),
  p_value  = c(t_Sp$`Pr(>F)`[1], t_Ar$`Pr(>F)`[1])
)
write.csv(pure_tests, file.path(out_dir, "vpa_pure_fractions_tests.csv"),
          row.names = FALSE)

# ============================================================
# 10. Biplot variable scores
# ============================================================
cat("\n--- 10. Biplot scores ---\n")
var_scores <- as.data.frame(scores(rda_full, display = "bp"))
print(var_scores)
write.csv(var_scores, file.path(out_dir, "variable_scores.csv"))

# ============================================================
# 11. Done
# ============================================================
cat(sprintf("\n=== Analysis complete: %s ===\n", date()))
cat(sprintf("Output dir: %s\n", out_dir))
sink()

