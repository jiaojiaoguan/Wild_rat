#!/usr/bin/env Rscript
# 3-panel horizontal: AMG mediation for 3 Area pairs
# Area → Bacteria PC1 → Virus AMG PC1, VLP + Bulk
library(ggplot2); library(patchwork)

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

# ---- Values (Bootstrap=10000, raw p-values) ----
vals <- data.frame(
  Pair = c("Livestock vs City", "Suburbs vs City", "Suburbs vs Livestock"),
  VLP_ACME   = c(0.0086,  0.0487,  0.0637),
  VLP_ACME_p = c(0.5598,  0.3152,  0.2578),
  VLP_ADE    = c(0.0282, -0.0776, -0.1347),
  VLP_ADE_p  = c(0.6086,  0.2572,  0.1108),
  Bulk_ACME   = c(-0.0048, -0.0843, -0.0877),
  Bulk_ACME_p = c(0.6168,  0.1366,  0.2280),
  Bulk_ADE    = c(0.0623,  0.1662,  0.1160),
  Bulk_ADE_p  = c(0.1188,  0.0588,  0.1316),
  stringsAsFactors = FALSE
)

fmt_p <- function(p) { if(p<0.001) "<0.001" else sprintf("%.4f", p) }

make_triangle <- function(v) {
  treat_label <- paste0("Area\n(", v$Pair, ")")

  s <- 3.5; h <- s*sqrt(3)/2

  nodes <- data.frame(
    name  = c(treat_label, "Bacterial PC1", "Virus AMGs PC1"),
    color = c("#D46A6A", "#599EC4", "#79C36A"),
    x = c(0, s/2, s), y = c(0, h, 0),
    lx = c(-0.35, s/2, s+0.35),
    ly = c(-0.3, h+0.3, -0.3),
    lhjust = c(0,0.5,1), lvjust = c(1,0,1))

  r <- 0.33
  edges <- data.frame(
    x0=c(nodes$x[1],nodes$x[1],nodes$x[2]),
    y0=c(nodes$y[1],nodes$y[1],nodes$y[2]),
    x1=c(nodes$x[2],nodes$x[3],nodes$x[3]),
    y1=c(nodes$y[2],nodes$y[3],nodes$y[3]),
    color=c("#D46A6A","#4E937A","#7A68A6"))
  edges$dx <- edges$x1-edges$x0; edges$dy <- edges$y1-edges$y0
  edges$len <- sqrt(edges$dx^2+edges$dy^2)
  edges$x0 <- edges$x0 + r*edges$dx/edges$len
  edges$y0 <- edges$y0 + r*edges$dy/edges$len
  edges$x1 <- edges$x1 - r*edges$dx/edges$len
  edges$y1 <- edges$y1 - r*edges$dy/edges$len

  mx <- (edges$x0+edges$x1)/2; my <- (edges$y0+edges$y1)/2
  mx[1] <- mx[1] + 0.85
  my[2] <- my[2] + 0.40

  labels <- c(
    sprintf("VLP: ACME=%.3f, p=%s\nBulk: ACME=%.3f, p=%s",
            abs(v$VLP_ACME), fmt_p(v$VLP_ACME_p),
            abs(v$Bulk_ACME), fmt_p(v$Bulk_ACME_p)),
    sprintf("VLP: ADE=%.3f, p=%s\nBulk: ADE=%.3f, p=%s",
            abs(v$VLP_ADE), fmt_p(v$VLP_ADE_p),
            abs(v$Bulk_ADE), fmt_p(v$Bulk_ADE_p)),
    "")

  ggplot() +
    geom_segment(aes(x=x0,y=y0,xend=x1,yend=y1,color=color), data=edges, linewidth=1.4,
                 arrow=arrow(length=unit(0.22,"cm"),type="closed")) +
    geom_point(aes(x=x,y=y,color=color), data=nodes, size=10, shape=16) +
    geom_text(aes(x=lx,y=ly,label=name,hjust=lhjust,vjust=lvjust), data=nodes,
              size=5.5, fontface="bold", color="grey20", lineheight=0.9) +
    geom_label(data=edges[labels!="",], aes(x=mx[labels!=""], y=my[labels!=""],
              label=labels[labels!=""], color=color[labels!=""]), size=4.5, lineheight=1.15,
              fill="white", label.size=0.4, label.padding=unit(0.22,"cm"),
              fontface="plain", show.legend=FALSE) +
    scale_color_identity() +
    coord_fixed(xlim=c(-1.8,s+1.7), ylim=c(-1.5,h+0.9), expand=FALSE) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(face="bold",hjust=0.5,size=15),
          plot.margin=margin(8,20,8,20),
          panel.background=element_rect(fill="white",color=NA))
}

p1 <- make_triangle(vals[1,])
p2 <- make_triangle(vals[2,])
p3 <- make_triangle(vals[3,])

combined <- p1 | p2 | p3

ggsave(file.path(fig_dir, "mediation_AMG_area_3panel.pdf"), combined, width=22, height=6, dpi=300, device=cairo_pdf)
ggsave(file.path(fig_dir, "mediation_AMG_area_3panel.png"), combined, width=22, height=6, dpi=300, type="cairo")
cat("Done: figures/mediation_AMG_area_3panel.pdf\n")
