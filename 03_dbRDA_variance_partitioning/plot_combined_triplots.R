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

# VLPs + Bulk RDA triplots side by side, clean background, matched to beta diversity size
library(vegan); library(readr); library(dplyr)

# ---- Shared triplot function ----
draw_triplot <- function(mod, meta_sub, scaling=2, title="") {
  site_sc <- scores(mod, display="sites", scaling=scaling, choices=1:2)
  bp_sc <- tryCatch(scores(mod, display="bp", scaling=scaling, choices=1:2), error=function(e) NULL)
  cn_sc <- tryCatch(scores(mod, display="cn", scaling=scaling, choices=1:2), error=function(e) NULL)
  if (!is.null(bp_sc) && !is.null(cn_sc) && nrow(cn_sc)>0)
    bp_sc <- bp_sc[!rownames(bp_sc) %in% rownames(cn_sc), , drop=FALSE]

  eig <- mod$CCA$eig; tot <- mod$tot.chi; pct <- 100*eig/tot
  xlab <- sprintf("RDA1 (%.1f%%)", pct[1]); ylab <- sprintf("RDA2 (%.1f%%)", pct[2])

  species_cols <- c("R.tanezumi"="#E39C3A","R.norvegicus"="#4BA3C3")
  area_shapes <- c("Livestock"=21,"Suburbs"=24,"City regions"=23)
  sp <- as.character(meta_sub$Species); ar <- as.character(meta_sub$Area)
  pt_bg <- species_cols[sp]; pt_pch <- area_shapes[ar]

  cn_labels <- NULL
  if (!is.null(cn_sc) && nrow(cn_sc)>0) {
    cn_labels <- rownames(cn_sc)
    cn_labels <- gsub("^Area","",cn_labels); cn_labels <- gsub("^Species","",cn_labels)
  }

  all_x <- c(site_sc[,1], if(!is.null(bp_sc)&&nrow(bp_sc)>0) bp_sc[,1], if(!is.null(cn_sc)&&nrow(cn_sc)>0) cn_sc[,1])
  all_y <- c(site_sc[,2], if(!is.null(bp_sc)&&nrow(bp_sc)>0) bp_sc[,2], if(!is.null(cn_sc)&&nrow(cn_sc)>0) cn_sc[,2])
  xlim <- range(all_x)*1.30; ylim <- c(-2.2, max(all_y)*1.3)

  bp_plot <- bp_sc
  if (!is.null(bp_sc) && nrow(bp_sc)>0) {
    s <- 0.85*min(max(abs(xlim))/max(abs(bp_sc[,1]),1e-9), max(abs(ylim))/max(abs(bp_sc[,2]),1e-9))
    bp_plot <- bp_sc * s
  }

  op <- par(mar=c(4.5,4.6,3.2,1.5), xpd=FALSE, family="sans", cex.axis=0.95)
  plot(site_sc, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, asp=1, main=title,
       cex.lab=1.15, cex.main=1.25, font.main=2, las=1,
       panel.first={ abline(h=0, v=0, lty=2, col="grey55", lwd=0.8) })

  points(site_sc, pch=pt_pch, bg=pt_bg, col="grey20", cex=1.8, lwd=0.7)

  if (!is.null(bp_plot) && nrow(bp_plot)>0) {
    arrows(0,0,bp_plot[,1],bp_plot[,2], length=0.10, angle=22, col="#1F3A93", lwd=1.8)
    text(bp_plot[,1]*1.12, bp_plot[,2]*1.12, labels=rownames(bp_plot), col="#1F3A93", font=2, cex=0.9)
  }

  if (!is.null(cn_sc) && nrow(cn_sc)>0) {
    points(cn_sc, pch=3, col="#B0306E", cex=1.6, lwd=2.2)
    text(cn_sc[,1], cn_sc[,2], labels=cn_labels, col="#B0306E", font=4, cex=0.95, pos=3, offset=0.45)
  }

  par(op)
}

# ---- Prep function for one dataset ----
prep_rda <- function(label) {
  tsv_path <- if(label=="VLPs") {
    file.path(data_dir, "53vlps_virus_map_VLPs_sample_RPKM_normalized_data.tsv")
  } else {
    file.path(data_dir, "53bulk_virus_map_bulk_sample_RPKM_normalized_data.tsv")
  }
  df <- read_tsv(tsv_path, show_col_types=FALSE)
  feats <- df[[1]]; df[[1]] <- NULL
  mat <- t(as.matrix(df)); colnames(mat) <- feats

  meta <- read_csv(file.path(data_dir,"wild_rat_metadata.csv"), show_col_types=FALSE)
  meta$Species <- dplyr::recode(meta$Species, "Rattus_andamanensis_(Rattus_tanezumi_sladeni)"="R.andamanensis","Rattus_norvegicus"="R.norvegicus","Rattus_tanezumi"="R.tanezumi")

  shared_X <- read.csv(file.path(data_dir,"explanatory_variables.csv"), row.names=1)

  anda <- meta$sample_name[meta$Species=="R.andamanensis"]
  common <- Reduce(intersect, list(rownames(mat), rownames(shared_X), meta$sample_name))
  common <- setdiff(common, c("Rat201", anda))
  mat <- mat[common,]; meta_sub <- meta[match(common, meta$sample_name),]; shared_X <- shared_X[common,]
  X <- data.frame(Species=factor(meta_sub$Species), Area=factor(meta_sub$Area),
    Bac_PCoA1=shared_X$Bac_PCoA1, Bac_PCoA2=shared_X$Bac_PCoA2, Bac_PCoA3=shared_X$Bac_PCoA3, row.names=common)
  rel <- decostand(mat, method="total"); d <- vegdist(rel, method="bray")
  mod <- dbrda(d ~ Species + Area + Bac_PCoA1 + Bac_PCoA2 + Bac_PCoA3, data=X)
  list(mod=mod, meta=meta_sub)
}

# ---- Generate combined plot ----
rda_vlp  <- prep_rda("VLPs")
rda_bulk <- prep_rda("Bulk")

draw_combined <- function() {
  layout(matrix(c(1,2,3), nrow=1), widths=c(0.42, 0.42, 0.16))
  draw_triplot(rda_vlp$mod,  rda_vlp$meta,  title="VLPs")
  draw_triplot(rda_bulk$mod, rda_bulk$meta, title="Bulk")
  # Legend panel
  par(mar=c(3,0.5,3,0.5))
  plot.new(); plot.window(xlim=c(0,10), ylim=c(0,10))
  species_cols <- c("R.tanezumi"="#E39C3A","R.norvegicus"="#4BA3C3")
  area_shapes <- c("Livestock"=21,"Suburbs"=24,"City regions"=23)
  legend("center",
    legend=c("Species",names(species_cols),"","Area",names(area_shapes)),
    pch=c(NA,21,21,NA,NA,21,24,23),
    pt.bg=c(NA,species_cols,NA,NA,rep("grey75",3)),
    col=c(NA,rep("grey20",2),NA,NA,rep("grey20",3)),
    text.font=c(2,3,3,1,2,1,1,1),
    pt.cex=1.8, bty="n", cex=1.1, y.intersp=1.2)
}

out_dir <- fig_dir
pdf(file.path(out_dir,"rda_triplots_combined_v2.pdf"), width=12, height=6.5)
draw_combined()
dev.off()
png(file.path(out_dir,"rda_triplots_combined_v2.png"), width=3600, height=1950, res=300)
draw_combined()
dev.off()
cat("Done: rda_triplots_combined_v2.pdf/png\n")
