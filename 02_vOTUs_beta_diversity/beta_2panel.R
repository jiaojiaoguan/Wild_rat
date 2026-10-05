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

# Beta-diversity PCoA — VLPs / Bulk with PERMANOVA (excluding R. andamanensis)
library(ggplot2); library(dplyr); library(readr); library(tibble); library(patchwork); library(vegan)

out_dir  <- fig_dir

# ---- Read matrices ----
read_sample_mat <- function(path) {
  df <- read_tsv(path, show_col_types=FALSE); feats <- df[[1]]; df[[1]] <- NULL
  mat <- t(as.matrix(df)); colnames(mat) <- feats; mat
}
vlp_mat  <- read_sample_mat(file.path(data_dir, "53vlps_virus_map_VLPs_sample_RPKM_normalized_data.tsv"))
bulk_mat <- read_sample_mat(file.path(data_dir, "53bulk_virus_map_bulk_sample_RPKM_normalized_data.tsv"))

common <- sort(Reduce(intersect, list(rownames(vlp_mat), rownames(bulk_mat))))
common <- setdiff(common, "Rat201")

# ---- Metadata (exclude R. andamanensis) ----
meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types=FALSE)
meta$Species <- dplyr::recode(meta$Species,
  "Rattus_andamanensis_(Rattus_tanezumi_sladeni)" = "R. andamanensis",
  "Rattus_norvegicus" = "R. norvegicus", "Rattus_tanezumi" = "R. tanezumi")
meta <- meta[meta$Species != "R. andamanensis", ]
meta$Species <- factor(meta$Species, levels=c("R. tanezumi", "R. norvegicus"))
meta$Area    <- factor(meta$Area, levels=c("Suburbs", "Livestock", "City regions"))

species_colors <- c("R. tanezumi"="#E39C3A", "R. norvegicus"="#4BA3C3")

# ---- Exclude R. andamanensis ----
andamanensis_samples <- meta$sample_name[as.character(meta$Species) == "R. andamanensis"]
common <- setdiff(common, andamanensis_samples)
meta <- meta[meta$Species != "R. andamanensis", ]
meta$Species <- droplevels(meta$Species)
vlp_mat  <- vlp_mat[common, ]
bulk_mat <- bulk_mat[common, ]

# ---- PCoA ----
compute_pcoa <- function(mat) {
  d <- vegdist(mat, method="bray")
  pcoa <- cmdscale(d, k=2, eig=TRUE)
  scores <- pcoa$points; colnames(scores) <- c("PC1","PC2")
  list(scores=scores, eig=pcoa$eig)
}

# Reference PCoA from VLP
ref <- compute_pcoa(vlp_mat)

# Procrustes rotate Bulk to VLP
procrustes_align <- function(X, Y) {
  Xc <- scale(X, scale=FALSE); Yc <- scale(Y, scale=FALSE)
  svd_res <- svd(t(Yc) %*% Xc)
  Q <- svd_res$u %*% t(svd_res$v)
  scale_factor <- sum(diag(t(Xc) %*% Yc %*% Q)) / sum(diag(t(Yc) %*% Yc))
  Yc %*% Q * scale_factor
}

pc_bulk <- compute_pcoa(bulk_mat)
pc_bulk$scores <- procrustes_align(ref$scores, pc_bulk$scores)
colnames(pc_bulk$scores) <- c("PC1","PC2")

# ---- PERMANOVA ----
run_permanova <- function(mat, lib) {
  df_meta <- meta[match(rownames(mat), meta$sample_name), ]
  df_meta <- df_meta[complete.cases(df_meta), ]
  mat <- mat[df_meta$sample_name, ]
  d <- vegdist(mat, method="bray")
  perm_A <- adonis2(d ~ Species + Area + Species:Area, data=df_meta, permutations=9999)
  perm_B <- adonis2(d ~ Area + Species + Species:Area, data=df_meta, permutations=9999)
  cat(sprintf("\n========== %s: R2 sensitivity (Species-first vs Area-first) ==========\n", lib))
  int_A <- grep(":", rownames(perm_A), value = TRUE)
  int_B <- grep(":", rownames(perm_B), value = TRUE)
  comp <- data.frame(
    Term       = c("Species","Area","Species:Area"),
    R2_SpFirst = round(c(perm_A["Species","R2"], perm_A["Area","R2"], perm_A[int_A,"R2"]), 4),
    R2_ArFirst = round(c(perm_B["Species","R2"], perm_B["Area","R2"], perm_B[int_B,"R2"]), 4),
    Delta      = round(c(perm_A["Species","R2"] - perm_B["Species","R2"],
                         perm_A["Area","R2"]      - perm_B["Area","R2"],
                         perm_A[int_A,"R2"] - perm_B[int_B,"R2"]), 4))
  print(comp, row.names=FALSE)
  list(R2_spec=perm_A["Species","R2"], p_spec=perm_A["Species","Pr(>F)"],
       R2_area=perm_A["Area","R2"],    p_area=perm_A["Area","Pr(>F)"],
       R2_int=perm_A["Species:Area","R2"], p_int=perm_A["Species:Area","Pr(>F)"])
}

pvlp  <- run_permanova(vlp_mat,  "VLPs")
pbulk <- run_permanova(bulk_mat, "Bulk")

# ---- Plot ----
make_plot <- function(scores, lib_label, perm_result, var_exp) {
  scores <- as.data.frame(scores)
  scores$sample <- rownames(ref$scores)
  scores <- merge(scores, meta, by.x="sample", by.y="sample_name")
  if(!"PC1" %in% colnames(scores)) colnames(scores)[2:3] <- c("PC1","PC2")

  anno <- sprintf("Species: R²=%.3f  p=%s\nArea: R²=%.3f  p=%s\nSpecies×Area: R²=%.3f  p=%s",
    perm_result$R2_spec, ifelse(perm_result$p_spec<0.001,"<0.001",sprintf("%.4f",perm_result$p_spec)),
    perm_result$R2_area, ifelse(perm_result$p_area<0.001,"<0.001",sprintf("%.4f",perm_result$p_area)),
    perm_result$R2_int,  ifelse(perm_result$p_int<0.001,"<0.001",sprintf("%.4f",perm_result$p_int)))

  ggplot(scores, aes(x=PC1, y=PC2, color=Species, shape=Area)) +
    geom_point(size=4.0, alpha=0.85) +
    stat_ellipse(aes(group=Species, color=Species, shape=NULL), type="norm", level=0.68, linetype=2, linewidth=0.8) +
    scale_color_manual(values=species_colors) +
    scale_shape_manual(values=c("Suburbs"=17,"Livestock"=16,"City regions"=18)) +
    annotate("text", x=min(scores$PC1)*1.3, y=max(scores$PC2)*1.80,
             label=anno, hjust=0, vjust=1, size=3.5, color="grey30") +
    labs(x=paste0("PC1 (",var_exp[1],"%)"), y=paste0("PC2 (",var_exp[2],"%)")) +
    theme_classic(base_size=16) +
    theme(legend.position="none", axis.title=element_text(face="bold",size=14), axis.text=element_text(size=12))
}

vvlp  <- round(100 * ref$eig[1:2] / sum(ref$eig), 1)
vbulk <- round(100 * pc_bulk$eig[1:2] / sum(pc_bulk$eig), 1)

p1 <- make_plot(ref$scores,     "VLPs", pvlp,  vvlp)
p2 <- make_plot(pc_bulk$scores, "Bulk", pbulk, vbulk)

# Add legend to Bulk panel
p2_legend <- p2 + theme(legend.position="right") +
  guides(color=guide_legend(title="Species"), shape=guide_legend(title="Area"))

all_plot <- p1 | p2_legend
ggsave(file.path(out_dir, "virus_beta_2panel.pdf"), all_plot, width=11, height=6.5)
ggsave(file.path(out_dir, "virus_beta_2panel.png"), all_plot, width=11, height=6.5)
cat("[DONE] virus_beta_2panel.pdf / .png\n")
