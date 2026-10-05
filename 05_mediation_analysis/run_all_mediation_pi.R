#!/usr/bin/env Rscript
# ============================================================
# Mediation analysis for viral π (nucleotide diversity)
# Paths: Species → MAG_pi → Virus_pi (control Area)
#        Area × 3 pairs → MAG_pi → Virus_pi (control Species)
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

library(readr); library(dplyr); library(mediation)

# ---- Metadata ----
meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types = FALSE)
meta <- meta %>% mutate(
  Species = dplyr::recode(Species,
    "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
    "Rattus_norvegicus" = "R. norvegicus", "Rattus_tanezumi" = "R. tanezumi"))
meta <- meta %>% filter(Species != "R. andamanensis")

# ---- MAG pi ----
mag <- read_tsv(file.path(data_dir, "MAGs_nucl_diversity_matrix.tsv"), show_col_types = FALSE)
mag_mean <- colMeans(mag[, -1], na.rm = TRUE)
mag_df <- data.frame(sample = names(mag_mean), MAG_pi = as.numeric(mag_mean))

results <- data.frame()

for (lib in c("VLP", "Bulk")) {
  # Virus pi
  csv <- if (lib == "VLP") file.path(data_dir, "VLPs_global_sample_microdiversity.csv")
         else                  file.path(data_dir, "bulk_global_sample_microdiversity.csv")
  virus <- read_csv(csv, show_col_types = FALSE)

  df <- virus %>%
    inner_join(meta, by = c("sample" = "sample_name")) %>%
    inner_join(mag_df, by = "sample") %>%
    filter(complete.cases(.))
  cat(sprintf("\n=== %s (n=%d) ===\n", lib, nrow(df)))

  # ---- Species pair ----
  cat(sprintf("\n--- %s: Species (R.norvegicus vs R.tanezumi) ---\n", lib))
  m.med <- lm(MAG_pi ~ Species + Area, data = df)
  m.out <- lm(avg_pi ~ Species + MAG_pi + Area, data = df)
  set.seed(42)
  med <- mediate(m.med, m.out, treat = "Species", mediator = "MAG_pi", boot = TRUE, sims = 10000)
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

    m.med <- lm(MAG_pi ~ Area + Species, data = df_sub)
    m.out <- lm(avg_pi ~ Area + MAG_pi + Species, data = df_sub)
    set.seed(42)
    med <- mediate(m.med, m.out, treat = "Area", mediator = "MAG_pi", boot = TRUE, sims = 10000)
    s <- summary(med)
    results <- rbind(results, data.frame(
      Dataset = lib, Path = "Area", Pair = pair_name,
      ACME = s$d0, ACME_ci_low = s$d0.ci[1], ACME_ci_high = s$d0.ci[2], ACME_p = s$d0.p,
      ADE = s$z0, ADE_ci_low = s$z0.ci[1], ADE_ci_high = s$z0.ci[2], ADE_p = s$z0.p,
      Total = s$tau.coef, Total_p = s$tau.p, PropMed = s$n0, PropMed_p = s$n0.p,
      stringsAsFactors = FALSE))
  }
}

# Save raw
write.table(results, file.path(res_dir, "mediation_pi_all_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---- FDR correction (within each dataset, 4 ACME + 4 ADE per dataset) ----
results$ACME_fdr <- NA; results$ADE_fdr <- NA
for (ds in c("VLP", "Bulk")) {
  idx <- results$Dataset == ds
  results$ACME_fdr[idx] <- p.adjust(results$ACME_p[idx], method = "fdr")
  results$ADE_fdr[idx]  <- p.adjust(results$ADE_p[idx], method = "fdr")
}

# Print
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
