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
# 7-boxplot: microdiversity (pi) across Area x Species — VLPs & Bulk
library(ggplot2); library(dplyr); library(readr); library(ggpubr); library(patchwork)

meta <- read_csv(file.path(data_dir, "wild_rat_metadata.csv"), show_col_types=FALSE)
meta <- meta %>% mutate(
  Species = dplyr::recode(Species, "Rattus_andamanensis_(Rattus_tanezumi_sladeni)"="R. andamanensis",
                   "Rattus_norvegicus"="R. norvegicus","Rattus_tanezumi"="R. tanezumi"),
  Area_label = Area)
meta$Species <- factor(meta$Species, levels=c("R. tanezumi","R. norvegicus"))
meta <- meta[meta$Species != "R. andamanensis", ]
meta$Area_label <- factor(meta$Area_label, levels=c("Suburbs","Livestock","City regions"))
species_colors <- c("R. tanezumi"="#E69F00","R. norvegicus"="#56B4E9","R. andamanensis"="#009E73")

make_panel <- function(tsv_path, lib_label) {
  df <- read_csv(tsv_path, show_col_types=FALSE)
  df <- inner_join(df, meta, by=c("sample"="sample_name"))
  df$group <- interaction(df$Species, df$Area_label, sep=" | ")
  # Overall Kruskal-Wallis across all 7+ groups
  kw <- kruskal.test(avg_pi ~ group, data=df)
  # Species within Area
  pw_spec <- compare_means(avg_pi ~ Species, group.by="Area_label", data=df,
                           method="wilcox.test", p.adjust.method="fdr")
  # Area within Species
  pw_area <- compare_means(avg_pi ~ Area_label, group.by="Species", data=df,
                           method="wilcox.test", p.adjust.method="fdr")
  cat(sprintf("\n%s: KW p=%.4f, Species-sig=%d, Area-sig=%d\n",
      lib_label, kw$p.value, sum(pw_spec$p.adj<=0.05), sum(pw_area$p.adj<=0.05)))
  if(sum(pw_spec$p.adj<=0.05)>0) { cat("  Species within Area (sig):\n")
    pw_spec %>% filter(p.adj<=0.05) %>% select(Area_label,group1,group2,p.adj) %>% print() }
  if(sum(pw_area$p.adj<=0.05)>0) { cat("  Area within Species (sig):\n")
    pw_area %>% filter(p.adj<=0.05) %>% select(Species,group1,group2,p.adj) %>% print() }

  # Build subtitle from significant pairwise results
  sub_lines <- c()
  if(sum(pw_spec$p.adj<=0.05)>0) {
    sigs <- pw_spec %>% filter(p.adj<=0.05)
    for(i in 1:nrow(sigs)) {
      sub_lines <- c(sub_lines, sprintf("%s: %s vs %s padj=%.4f",
        sigs$Area_label[i], sigs$group1[i], sigs$group2[i], sigs$p.adj[i]))
    }
  }
  if(sum(pw_area$p.adj<=0.05)>0) {
    sigs <- pw_area %>% filter(p.adj<=0.05)
    for(i in 1:nrow(sigs)) {
      sub_lines <- c(sub_lines, sprintf("%s: %s vs %s padj=%.4f",
        sigs$Species[i], sigs$group1[i], sigs$group2[i], sigs$p.adj[i]))
    }
  }
  sub_text <- if(length(sub_lines)>0) paste(sub_lines, collapse="\n") else ""

  p <- ggplot(df, aes(x=Area_label, y=avg_pi, fill=Species)) +
    geom_boxplot(outlier.shape=21, outlier.color="black", outlier.size=1.2, alpha=0.8, linewidth=0.35) +
    scale_fill_manual(values=species_colors) +
    labs(title=lib_label, subtitle=sub_text, x=NULL, y=expression(Average~pi)) +
    theme_classic(base_size=14) +
    theme(plot.title=element_text(hjust=0.5, size=13, face="bold"),
          plot.subtitle=element_text(size=8, hjust=0.5, color="grey30"),
          axis.text=element_text(size=11, color="black"),
          axis.title.y=element_text(size=13), legend.position="none")

  # Show KW p-value
  p <- p + annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=2.5,
           label=sprintf("Kruskal-Wallis p = %.4f", kw$p.value),
           size=3.5, fontface="italic", color="grey25")
  # If significant, also print significant pairwise to console (already done above)
  list(plot=p, pw_spec=pw_spec, pw_area=pw_area)
}

r_vlp  <- make_panel(file.path(data_dir, "VLPs_global_sample_microdiversity.csv"), "VLPs")
r_bulk <- make_panel(file.path(data_dir, "bulk_global_sample_microdiversity.csv"), "Bulk")

write_tsv(r_vlp$pw_spec,  file.path(res_dir, "VLPs_pi_Species_within_Area.tsv"))
write_tsv(r_vlp$pw_area,  file.path(res_dir, "VLPs_pi_Area_within_Species.tsv"))
write_tsv(r_bulk$pw_spec, file.path(res_dir, "Bulk_pi_Species_within_Area.tsv"))
write_tsv(r_bulk$pw_area, file.path(res_dir, "Bulk_pi_Area_within_Species.tsv"))

combined <- (r_vlp$plot | r_bulk$plot) + plot_layout(guides="collect") &
  theme(legend.position="right") & scale_y_continuous(limits=c(0, 0.001))
ggsave(file.path(fig_dir, "pi_by_species_area_7box.pdf"), combined, width=10, height=6, dpi=300)
ggsave(file.path(fig_dir, "pi_by_species_area_7box.png"), combined, width=10, height=6, dpi=300)
cat("\nDone!\n")
