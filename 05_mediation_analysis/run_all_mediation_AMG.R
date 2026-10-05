#!/usr/bin/env Rscript
# ============================================================
# Mediation analysis: Species/Area → Bacteria PC1 → Virus AMG PC1
# Virus AMG PC1 = PCoA axis 1 of AMG KO composition (Bray-Curtis)
# Bacteria PC1 = Bac_PCoA1 from shared explanatory variables
# Bootstrap=10000, FDR-corrected within dataset
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

library(readr); library(dplyr); library(vegan); library(mediation)

# ---- Metadata ----
meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types = FALSE)
meta <- meta %>% mutate(
  Species = dplyr::recode(Species,
    "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
    "Rattus_norvegicus" = "R. norvegicus", "Rattus_tanezumi" = "R. tanezumi"))
meta <- meta %>% filter(Species != "R. andamanensis")

# ---- Shared X (Bacteria PC1) ----
shared_X <- read.csv(file.path(data_dir, "explanatory_variables.csv"), row.names = 1)

results <- data.frame()

for (lib in c("VLP", "Bulk")) {
  # ---- Read AMG KO abundance & compute PCoA PC1 ----
  amg_file <- if (lib == "VLP") file.path(data_dir, "53vlp_AMG_abundance_KO_normalized.tsv")
             else                  file.path(data_dir, "53meta_AMG_abundance_KO_normalized.tsv")
  amg <- read_tsv(amg_file, show_col_types = FALSE)
  feats <- amg[[1]]; amg[[1]] <- NULL
  amg_mat <- t(as.matrix(amg)); colnames(amg_mat) <- feats

  # Relative abundance + Bray-Curtis + PCoA
  amg_rel <- decostand(amg_mat, method = "total")
  amg_dist <- vegdist(amg_rel, method = "bray")
  amg_pcoa <- cmdscale(amg_dist, k = 2, eig = TRUE)
  amg_pc1 <- amg_pcoa$points[, 1]
  amg_pc1_df <- data.frame(sample = names(amg_pc1), AMG_PC1 = as.numeric(amg_pc1))

  # ---- Merge ----
  df <- amg_pc1_df %>%
    inner_join(meta, by = c("sample" = "sample_name")) %>%
    inner_join(data.frame(sample = rownames(shared_X), Bac_PC1 = shared_X$Bac_PCoA1), by = "sample") %>%
    filter(!sample %in% "Rat201", complete.cases(.))
  cat(sprintf("\n=== %s (n=%d) ===\n", lib, nrow(df)))

  # ---- Species pair ----
  cat(sprintf("\n--- %s: Species (R.norvegicus vs R.tanezumi) ---\n", lib))
  m.med <- lm(Bac_PC1 ~ Species + Area, data = df)
  m.out <- lm(AMG_PC1 ~ Species + Bac_PC1 + Area, data = df)
  set.seed(42)
  med <- mediate(m.med, m.out, treat = "Species", mediator = "Bac_PC1", boot = TRUE, sims = 10000)
  s <- summary(med)
  results <- rbind(results, data.frame(
    Dataset = lib, Path = "Species", Pair = "R.norvegicus vs R.tanezumi",
    ACME = s$d0, ACME_ci_low = s$d0.ci[1], ACME_ci_high = s$d0.ci[2], ACME_p = s$d0.p,
    ADE = s$z0, ADE_ci_low = s$z0.ci[1], ADE_ci_high = s$z0.ci[2], ADE_p = s$z0.p,
    Total = s$tau.coef, Total_p = s$tau.p, PropMed = s$n0, PropMed_p = s$n0.p,
    stringsAsFactors = FALSE))

  # ---- Area pairs (3 pairwise) ----
  area_pairs <- list(
    c("Livestock", "City regions"),
    c("Suburbs", "City regions"),
    c("Suburbs", "Livestock")
  )
  for (ap in area_pairs) {
    pair_name <- paste(ap[1], "vs", ap[2])
    cat(sprintf("\n--- %s: Area (%s) ---\n", lib, pair_name))
    df_sub <- df[df$Area %in% ap, ]
    df_sub$Area <- factor(df_sub$Area, levels = ap)

    m.med <- lm(Bac_PC1 ~ Area + Species, data = df_sub)
    m.out <- lm(AMG_PC1 ~ Area + Bac_PC1 + Species, data = df_sub)
    set.seed(42)
    med <- mediate(m.med, m.out, treat = "Area", mediator = "Bac_PC1", boot = TRUE, sims = 10000)
    s <- summary(med)
    results <- rbind(results, data.frame(
      Dataset = lib, Path = "Area", Pair = pair_name,
      ACME = s$d0, ACME_ci_low = s$d0.ci[1], ACME_ci_high = s$d0.ci[2], ACME_p = s$d0.p,
      ADE = s$z0, ADE_ci_low = s$z0.ci[1], ADE_ci_high = s$z0.ci[2], ADE_p = s$z0.p,
      Total = s$tau.coef, Total_p = s$tau.p, PropMed = s$n0, PropMed_p = s$n0.p,
      stringsAsFactors = FALSE))
  }
}

# ---- Save raw ----
write.table(results, file.path(res_dir, "mediation_AMG_all_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---- FDR correction (within each dataset, 4 ACME + 4 ADE) ----
results$ACME_fdr <- NA; results$ADE_fdr <- NA
for (ds in c("VLP", "Bulk")) {
  idx <- results$Dataset == ds
  results$ACME_fdr[idx] <- p.adjust(results$ACME_p[idx], method = "fdr")
  results$ADE_fdr[idx]  <- p.adjust(results$ADE_p[idx], method = "fdr")
}

# ---- Print ----
cat("\n========== ALL RESULTS (Bootstrap=10000, FDR-corrected within dataset) ==========\n")
for (i in 1:nrow(results)) {
  r <- results[i, ]
  sig_acme <- ifelse(r$ACME_fdr < 0.05, "*", " ")
  sig_ade  <- ifelse(r$ADE_fdr < 0.05, "*", " ")
  cat(sprintf("%s | %-10s | %-35s | ACME=%.4f p=%.4f (FDR=%.4f)%s | ADE=%.4f p=%.4f (FDR=%.4f)%s\n",
      r$Dataset, r$Path, r$Pair, r$ACME, r$ACME_p, r$ACME_fdr, sig_acme,
      r$ADE, r$ADE_p, r$ADE_fdr, sig_ade))
}
cat("\nDone!\n")
