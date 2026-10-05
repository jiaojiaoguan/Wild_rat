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

options(encoding = "UTF-8")

#############################################################################
# Model 2: Partial db-RDA — VLPs virus ~ Species + Area + Bac_PCoA1-3
# Bray-Curtis distance, varpart() for 3-way decomposition
#############################################################################

library(vegan)
library(readr)
library(dplyr)

set.seed(42)

# ============================================================
# 0. Paths
# ============================================================
out_dir  <- file.path(res_dir, "VLPs_RDA_results_species_area_bacpcoa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sink(file.path(out_dir, "analysis_log.txt"), split = TRUE)
cat("=== Partial RDA: VLPs virus ~ Species + Area + Bac_PCoA1-3 ===\n")
cat("Started:", date(), "\n\n")

# ============================================================
# 1. Read data (SAME pipeline as Model 1)
# ============================================================
cat("--- 1. Reading data ---\n")

df <- read_tsv(file.path(data_dir,
                         "53vlps_virus_map_VLPs_sample_RPKM_normalized_data.tsv"),
               show_col_types = FALSE)
feats <- df[[1]]; df[[1]] <- NULL
vlp_mat <- t(as.matrix(df)); colnames(vlp_mat) <- feats
cat("VLPs virus:", nrow(vlp_mat), "samples x", ncol(vlp_mat), "vOTUs\n")

meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"),
                 show_col_types = FALSE)
meta$Species <- dplyr::recode(meta$Species,
                       "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
                       "Rattus_norvegicus" = "R. norvegicus",
                       "Rattus_tanezumi"   = "R. tanezumi")

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

common <- Reduce(intersect, list(rownames(vlp_mat),
                                 rownames(shared_X),
                                 meta$sample_name))
common <- setdiff(common, c("Rat201", andamanensis_samples))
cat("Common samples after exclusion:", length(common), "\n")

vlp_mat  <- vlp_mat[common, , drop = FALSE]
meta_sub  <- meta[match(common, meta$sample_name), ]
shared_X  <- shared_X[common, , drop = FALSE]

# ============================================================
# 3. Build explanatory table X
# ============================================================
X <- data.frame(
  Species   = factor(meta_sub$Species),
  Area      = factor(meta_sub$Area),
  Bac_PCoA1 = as.numeric(shared_X$Bac_PCoA1),
  Bac_PCoA2 = as.numeric(shared_X$Bac_PCoA2),
  Bac_PCoA3 = as.numeric(shared_X$Bac_PCoA3),
  row.names = common
)
cat("\nX dimensions:", nrow(X), "x", ncol(X), "\n")
cat("Species table:\n"); print(table(X$Species))
cat("Area table:\n");    print(table(X$Area))
write.csv(X, file.path(out_dir, "explanatory_variables.csv"))

# ============================================================
# 4. Relative abundance + Bray-Curtis
# ============================================================
cat("\n--- 4. Bray-Curtis distance ---\n")
vlp_rel <- decostand(vlp_mat, method = "total")
d <- vegdist(vlp_rel, method = "bray")
cat("Distance range:", round(range(d), 4), "\n")

# ============================================================
# 5. Full model dbRDA
# ============================================================
cat("\n--- 5. Full dbrda model ---\n")
rda_full <- dbrda(d ~ Species + Area + Bac_PCoA1 + Bac_PCoA2 + Bac_PCoA3,
                  data = X)
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
# 7. Partial RDA — Drop-one nested tests
# ============================================================
cat("\n--- 7. Drop-one nested tests ---\n")
variables <- c("Species", "Area", "Bac_PCoA1", "Bac_PCoA2", "Bac_PCoA3")
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
# 8. Variance Partitioning — 3 groups via varpart()
# ============================================================
cat("\n--- 8. Variance Partitioning (varpart) ---\n")

vp <- varpart(d,
              ~ Species,
              ~ Area,
              ~ Bac_PCoA1 + Bac_PCoA2 + Bac_PCoA3,
              data = X)
print(vp)

ind <- vp$part$indfract
adj_col <- grep("Adj", colnames(ind), value = TRUE)[1]
cat("\nFractions (rows):", paste(rownames(ind), collapse = " | "), "\n")
cat("Using column:", adj_col, "\n\n")

stopifnot(nrow(ind) == 8)
vals <- as.numeric(ind[, adj_col])
names(vals) <- c("[a] Pure Species", "[b] Pure Area", "[c] Pure Bac_PCoA",
                 "[d] Shared Sp&Ar",
                 "[e] Shared Ar&Bac",
                 "[f] Shared Sp&Bac",
                 "[g] Shared Sp&Ar&Bac",
                 "[h] Residuals")

explained <- sum(vals[1:7], na.rm = TRUE)

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
                 data.frame(Fraction = "SUM Explained [a..g]",
                            AdjR2 = round(explained, 6),
                            Pct   = round(explained * 100, 2)))
write.csv(vpa_all, file.path(out_dir, "vpa_fractions.csv"), row.names = FALSE)

# ============================================================
# 8b. Venn diagram for variance partitioning
# ============================================================
cat("\n--- 8b. Plotting Venn diagram ---\n")

pdf(file.path(fig_dir, "VLPs_vpa_venn_default.pdf"), width = 7, height = 6)
plot(vp,
     Xnames = c("Species", "Area", "Bac_PCoA"),
     bg     = c("#FDB4624D", "#79C36A4D", "#599EC44D"),
     id.size = 1.0, cex = 1.1, digits = 3)
title(main = "Variance partitioning of VLPs virus community\n(Adj R² shown in each fraction)",
      cex.main = 1.0)
mtext(sprintf("Residuals = %.3f   |   Total explained = %.3f",
              vals["[h] Residuals"], explained),
      side = 1, line = 2, cex = 0.9)
dev.off()

png(file.path(fig_dir, "VLPs_vpa_venn_default.png"),
    width = 1800, height = 1500, res = 250)
plot(vp,
     Xnames = c("Species", "Area", "Bac_PCoA"),
     bg     = c("#FDB4624D", "#79C36A4D", "#599EC44D"),
     id.size = 1.0, cex = 1.1, digits = 3)
title(main = "Variance partitioning of VLPs virus community\n(Adj R² shown in each fraction)",
      cex.main = 1.0)
mtext(sprintf("Residuals = %.3f   |   Total explained = %.3f",
              vals["[h] Residuals"], explained),
      side = 1, line = 2, cex = 0.9)
dev.off()

make_venn <- function(){
  par(mar = c(2, 1, 4, 1), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 10), asp = 1)
  
  symbols(x = 3.5, y = 6.5, circles = 2.6, inches = FALSE, add = TRUE,
          fg = "#E07B39", bg = "#FDB46280", lwd = 2)
  symbols(x = 6.5, y = 6.5, circles = 2.6, inches = FALSE, add = TRUE,
          fg = "#4A8B3F", bg = "#79C36A80", lwd = 2)
  symbols(x = 5.0, y = 3.8, circles = 2.6, inches = FALSE, add = TRUE,
          fg = "#2D6B8F", bg = "#599EC480", lwd = 2)
  
  text(1.2, 9.3, "X1: Species", col = "#E07B39", font = 2, cex = 1.2)
  text(8.8, 9.3, "X2: Area",    col = "#4A8B3F", font = 2, cex = 1.2)
  text(5.0, 0.7, "X3: Bac_PCoA",col = "#2D6B8F", font = 2, cex = 1.2)
  
  fmt <- function(x) sprintf("%.4f\n(%.2f%%)", x, x*100)
  
  text(2.0, 6.5, paste0("[a]\n", fmt(vals["[a] Pure Species"])),    cex = 0.85, font = 2)
  text(8.0, 6.5, paste0("[b]\n", fmt(vals["[b] Pure Area"])),       cex = 0.85, font = 2)
  text(5.0, 2.5, paste0("[c]\n", fmt(vals["[c] Pure Bac_PCoA"])),   cex = 0.85, font = 2)
  text(5.0, 7.5, paste0("[d]\n", fmt(vals["[d] Shared Sp&Ar"])),    cex = 0.8)
  text(3.5, 4.7, paste0("[f]\n", fmt(vals["[f] Shared Sp&Bac"])),   cex = 0.8)
  text(6.5, 4.7, paste0("[e]\n", fmt(vals["[e] Shared Ar&Bac"])),   cex = 0.8)
  text(5.0, 5.6, paste0("[g]\n", fmt(vals["[g] Shared Sp&Ar&Bac"])),cex = 0.8, font = 3)
  
  text(5.0, -0.3,
       sprintf("[h] Residuals = %.4f (%.2f%%)   |   Total explained = %.4f (%.2f%%)",
               vals["[h] Residuals"], vals["[h] Residuals"]*100,
               explained, explained*100),
       cex = 0.9, font = 2)
  
  title(main = "Variance partitioning — VLPs virus ~ Species + Area + Bac_PCoA",
        cex.main = 1.05)
}

pdf(file.path(fig_dir, "VLPs_vpa_venn_custom.pdf"), width = 7.5, height = 8); make_venn(); dev.off()
png(file.path(fig_dir, "VLPs_vpa_venn_custom.png"),
    width = 2000, height = 1900, res = 260); make_venn(); dev.off()

cat("Venn diagrams saved.\n")

# ============================================================
# 9. Significance tests for pure (testable) fractions
# ============================================================
cat("\n--- 9. Significance of pure (testable) fractions ---\n")

set.seed(42); t_Sp <- anova(dbrda(d ~ Species +
                                    Condition(Area + Bac_PCoA1 + Bac_PCoA2 + Bac_PCoA3),
                                  data = X), permutations = 9999)
set.seed(42); t_Ar <- anova(dbrda(d ~ Area +
                                    Condition(Species + Bac_PCoA1 + Bac_PCoA2 + Bac_PCoA3),
                                  data = X), permutations = 9999)
set.seed(42); t_Bc <- anova(dbrda(d ~ Bac_PCoA1 + Bac_PCoA2 + Bac_PCoA3 +
                                    Condition(Species + Area),
                                  data = X), permutations = 9999)

cat("\n[a] Pure Species | Area+Bac:\n"); print(t_Sp)
cat("\n[b] Pure Area | Species+Bac:\n"); print(t_Ar)
cat("\n[c] Pure Bac  | Species+Area:\n"); print(t_Bc)

pure_tests <- data.frame(
  Fraction = c("[a] Pure Species", "[b] Pure Area", "[c] Pure Bac_PCoA"),
  F_stat   = c(t_Sp$F[1], t_Ar$F[1], t_Bc$F[1]),
  p_value  = c(t_Sp$`Pr(>F)`[1], t_Ar$`Pr(>F)`[1], t_Bc$`Pr(>F)`[1])
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
# 11. RDA triplot (sites + continuous bp arrows + factor centroids)
#     - color  = Species
#     - shape  = Area  (Livestock=circle, Suburbs=triangle, City regions=diamond)
#     - legend = inside the plot frame
# ============================================================
cat("\n--- 11. RDA triplot ---\n")

draw_rda_triplot <- function(mod, meta_sub, scaling = 2,
                             title = "db-RDA triplot — VLPs virus") {
  # ---- 1) get scores ----
  site_sc <- scores(mod, display = "sites", scaling = scaling, choices = 1:2)
  bp_sc   <- tryCatch(scores(mod, display = "bp",  scaling = scaling, choices = 1:2),
                      error = function(e) NULL)
  cn_sc   <- tryCatch(scores(mod, display = "cn",  scaling = scaling, choices = 1:2),
                      error = function(e) NULL)
  
  # >>> remove factor level rows from bp to avoid overlap with cn <<<
  if (!is.null(bp_sc) && !is.null(cn_sc) && nrow(cn_sc) > 0) {
    bp_sc <- bp_sc[!rownames(bp_sc) %in% rownames(cn_sc), , drop = FALSE]
  }
  
  # ---- 2) axis % explained ----
  eig         <- mod$CCA$eig
  tot_inertia <- mod$tot.chi
  pct         <- 100 * eig / tot_inertia
  xlab <- sprintf("RDA1 (%.1f%%)", pct[1])
  ylab <- sprintf("RDA2 (%.1f%%)", pct[2])
  
  # ---- 3) color / shape mapping ----
  species_cols <- c("R. tanezumi"     = "#E39C3A",
                    "R. norvegicus"   = "#4BA3C3",
                    "R. andamanensis" = "#4E937A")
  area_shapes  <- c("Livestock" = 21,   # circle
                    "Suburbs" = 24,   # triangle
                    "City regions" = 23)   # diamond
  
  sp     <- as.character(meta_sub$Species)
  ar     <- as.character(meta_sub$Area)
  pt_bg  <- species_cols[sp]
  pt_pch <- area_shapes[ar]
  
  # clean centroid labels (remove Species/Area prefix)
  cn_labels <- NULL
  if (!is.null(cn_sc) && nrow(cn_sc) > 0) {
    cn_labels <- rownames(cn_sc)
    cn_labels <- gsub("^Area",    "", cn_labels)
    cn_labels <- gsub("^Species", "", cn_labels)
  }
  
  # ---- 4) axis limits ----
  all_x <- c(site_sc[,1],
             if (!is.null(bp_sc) && nrow(bp_sc) > 0) bp_sc[,1],
             if (!is.null(cn_sc) && nrow(cn_sc) > 0) cn_sc[,1])
  all_y <- c(site_sc[,2],
             if (!is.null(bp_sc) && nrow(bp_sc) > 0) bp_sc[,2],
             if (!is.null(cn_sc) && nrow(cn_sc) > 0) cn_sc[,2])
 xlim <- range(all_x) * 1.30
 ylim <- c(-2.2, max(all_y) * 1.3)
  
  # ---- 5) scale bp arrows proportionally ----
  bp_plot <- bp_sc
  if (!is.null(bp_sc) && nrow(bp_sc) > 0) {
    s <- 0.85 * min(max(abs(xlim)) / max(abs(bp_sc[,1]), 1e-9),
                    max(abs(ylim)) / max(abs(bp_sc[,2]), 1e-9))
    bp_plot <- bp_sc * s
  }
  
  # ---- 6) draw ----
  op <- par(mar = c(4.5, 4.6, 3.2, 1.5), xpd = FALSE,
            family = "sans", cex.axis = 0.95)
  plot(site_sc, type = "n", xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, asp = 1, main = title,
       cex.lab = 1.15, cex.main = 1.25, font.main = 2, las = 1,
       panel.first = {
         abline(h = 0, v = 0, lty = 2, col = "grey55", lwd = 0.8)
       })
  
  # sites
  points(site_sc, pch = pt_pch, bg = pt_bg, col = "grey20",
         cex = 1.8, lwd = 0.7)
  
  # continuous vars (bp) blue arrows
  if (!is.null(bp_plot) && nrow(bp_plot) > 0) {
    arrows(0, 0, bp_plot[,1], bp_plot[,2],
           length = 0.10, angle = 22, col = "#1F3A93", lwd = 1.8)
    text(bp_plot[,1] * 1.12, bp_plot[,2] * 1.12,
         labels = rownames(bp_plot),
         col = "#1F3A93", font = 2, cex = 0.9)
  }
  
  # factor centroid (cn) cross + magenta labels
  if (!is.null(cn_sc) && nrow(cn_sc) > 0) {
    points(cn_sc, pch = 3, col = "#B0306E", cex = 1.6, lwd = 2.2)
    text(cn_sc[,1], cn_sc[,2],
         labels = cn_labels, col = "#B0306E", font = 4,
         cex = 0.95, pos = 3, offset = 0.45)
  }
  
  # ---- 7) legend ----
  sp_levels <- intersect(names(species_cols), unique(sp))
  ar_levels <- intersect(names(area_shapes),  unique(ar))
  
  leg_lab <- c("Species",                  sp_levels,
               "",
               "Area",                     ar_levels)
  leg_pch <- c(NA, rep(21, length(sp_levels)),
               NA,
               NA, area_shapes[ar_levels])
  leg_bg  <- c(NA, species_cols[sp_levels],
               NA,
               NA, rep("grey75", length(ar_levels)))
  leg_col <- c(NA, rep("grey20", length(sp_levels)),
               NA,
               NA, rep("grey20", length(ar_levels)))
  leg_font<- c(2,  rep(3, length(sp_levels)),
               1,
               2,  rep(1, length(ar_levels)))
  
  legend("bottomleft",
         inset      = c(0.015, 0.015),
         legend     = leg_lab,
         pch        = leg_pch,
         pt.bg      = leg_bg,
         col        = leg_col,
         pt.cex     = 1.5,
         text.font  = leg_font,
         bty        = "o",
         bg         = scales::alpha("white", 0.92),
         box.col    = "grey70",
         box.lwd    = 0.8,
         cex        = 0.88,
         y.intersp  = 1.05,
         seg.len    = 1.2)
  
  par(op)
}

# ---- ... PDF + PNG ----
pdf(file.path(fig_dir, "VLPs_rda_triplot.pdf"), width = 7, height = 8.5)
draw_rda_triplot(rda_full, meta_sub, scaling = 2,
                 title = "db-RDA triplot — VLPs virus")
dev.off()

png(file.path(fig_dir, "VLPs_rda_triplot.png"),
    width = 2200, height = 1800, res = 280)
draw_rda_triplot(rda_full, meta_sub, scaling = 2,
                 title = "db-RDA triplot — VLPs virus")
dev.off()

cat("Triplot saved:\n")
cat("  ", file.path(fig_dir, "VLPs_rda_triplot.pdf"), "\n")
cat("  ", file.path(fig_dir, "VLPs_rda_triplot.png"), "\n")

cat("\n=== Analysis complete: ", date(), " ===\n")
sink()