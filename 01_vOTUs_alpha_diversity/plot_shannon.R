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
# Shannon 7-boxplot — WITHOUT R. andamanensis
library(ggplot2); library(dplyr); library(readr); library(ggpubr); library(patchwork)

meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types=FALSE)
meta <- meta %>% mutate(
  Species = dplyr::recode(Species, "Rattus_andamanensis_(Rattus_tanezumi_sladeni)"="R. andamanensis",
                   "Rattus_norvegicus"="R. norvegicus","Rattus_tanezumi"="R. tanezumi"),
  Area_label = Area)
meta$Species <- factor(meta$Species, levels=c("R. tanezumi","R. norvegicus"))
meta$Area_label <- factor(meta$Area_label, levels=c("Suburbs","Livestock","City regions"))
meta <- meta[meta$Species != "R. andamanensis", ]
species_colors <- c("R. tanezumi"="#E69F00","R. norvegicus"="#56B4E9")

sig_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

make_panel <- function(tsv, lib_label) {
  div <- read_tsv(tsv, show_col_types=FALSE)
  colnames(div) <- c("sample","Richness","Shannon")
  df <- inner_join(div, meta, by=c("sample"="sample_name"))
  df$group <- interaction(df$Species, df$Area_label, sep=" | ")
  kw <- kruskal.test(Shannon ~ group, data=df)
  pw_spec <- compare_means(Shannon ~ Species, group.by="Area_label", data=df, method="wilcox.test", p.adjust.method="fdr")
  pw_area <- compare_means(Shannon ~ Area_label, group.by="Species", data=df, method="wilcox.test", p.adjust.method="fdr")
  pw_spec$p.signif <- sig_stars(pw_spec$p.adj)
  pw_area$p.signif <- sig_stars(pw_area$p.adj)
  cat(sprintf("%s: KW p=%.4f, Species-sig=%d, Area-sig=%d\n", lib_label, kw$p.value, sum(pw_spec$p.adj<=0.05), sum(pw_area$p.adj<=0.05)))

  sub_lines <- c()
  if(sum(pw_spec$p.adj<=0.05)>0) {
    sigs <- pw_spec %>% filter(p.adj<=0.05)
    for(i in 1:nrow(sigs)) sub_lines <- c(sub_lines, sprintf("%s: %s vs %s padj=%.4f", sigs$Area_label[i], sigs$group1[i], sigs$group2[i], sigs$p.adj[i]))
  }
  if(sum(pw_area$p.adj<=0.05)>0) {
    sigs <- pw_area %>% filter(p.adj<=0.05)
    for(i in 1:nrow(sigs)) sub_lines <- c(sub_lines, sprintf("%s: %s vs %s padj=%.4f", sigs$Species[i], sigs$group1[i], sigs$group2[i], sigs$p.adj[i]))
  }
  sub_text <- if(length(sub_lines)>0) paste(sub_lines, collapse="\n") else ""
  p <- ggplot(df, aes(x=Area_label, y=Shannon, fill=Species)) +
    geom_boxplot(outlier.shape=21, outlier.color="black", outlier.size=1.2, alpha=0.8, linewidth=0.35) +
    scale_fill_manual(values=species_colors) +
    labs(title=lib_label, subtitle=sub_text, x=NULL, y="Shannon diversity") +
    theme_classic(base_size=14) +
    theme(plot.title=element_text(hjust=0.5, size=13, face="bold"),
          plot.subtitle=element_text(size=8, hjust=0.5, color="grey30"),
          axis.text=element_text(size=11, color="black"),
          axis.title.y=element_text(size=13), legend.position="none") +
    annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=2.5,
             label=sprintf("Kruskal-Wallis p = %.4f", kw$p.value), size=3.5, fontface="italic", color="grey25")
  list(plot=p, pw_spec=pw_spec, pw_area=pw_area)
}
r_vlp  <- make_panel(file.path(data_dir,"VLPs_virus_alpha_diversity.tsv"), "VLPs")
r_bulk <- make_panel(file.path(data_dir,"bulk_virus_alpha_diversity.tsv"), "Bulk")

write_tsv(r_vlp$pw_spec,  file.path(res_dir, "VLPs_Shannon_7box_Species_within_Area_without_anda.tsv"))
write_tsv(r_vlp$pw_area,  file.path(res_dir, "VLPs_Shannon_7box_Area_within_Species_without_anda.tsv"))
write_tsv(r_bulk$pw_spec, file.path(res_dir, "Bulk_Shannon_7box_Species_within_Area_without_anda.tsv"))
write_tsv(r_bulk$pw_area, file.path(res_dir, "Bulk_Shannon_7box_Area_within_Species_without_anda.tsv"))

combined <- (r_vlp$plot | r_bulk$plot) + plot_layout(guides="collect") & theme(legend.position="right") & scale_y_continuous(limits=c(3, 7))
ggsave(file.path(fig_dir, "shannon_7box_without_anda.pdf"), combined, width=11, height=6.5, dpi=300)
ggsave(file.path(fig_dir, "shannon_7box_without_anda.png"), combined, width=11, height=6.5, dpi=300)
cat("\nDone!\n")
