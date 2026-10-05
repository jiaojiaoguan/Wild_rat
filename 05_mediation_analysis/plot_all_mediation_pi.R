#!/usr/bin/env Rscript
# 4-panel mediation for viral π: 1 Species + 3 Area pairs, VLP+Bulk per panel
# Mediator: MAG π (bacterial nucleotide diversity)
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

# ---- Values (raw p-values, 10000 bootstrap) ----
vals <- data.frame(
  Path = c("Species","Area","Area","Area"),
  Pair = c("","Livestock vs City","Suburbs vs City","Suburbs vs Livestock"),
  VLP_ACME   = c(1.272e-06, -1.349e-06, -2.286e-07, -3.485e-07),
  VLP_ACME_p = c(0.9176, 0.9566, 0.5968, 0.7996),
  VLP_ADE    = c(-8.342e-05, -1.782e-04, 3.702e-05, 2.032e-04),
  VLP_ADE_p  = c(0.0544, 0.0004, 0.3540, 0.0016),
  Bulk_ACME   = c(2.146e-04, -2.303e-05, 1.171e-05, 2.506e-05),
  Bulk_ACME_p = c(0.0054, 0.6810, 0.9006, 0.6978),
  Bulk_ADE    = c(-8.402e-05, 2.015e-04, 3.267e-04, 1.067e-04),
  Bulk_ADE_p  = c(0.1882, 0.0024, 0.0006, 0.2848),
  stringsAsFactors = FALSE
)

fmt_p <- function(p) { if(p<0.0001) "<0.0001" else sprintf("%.4f", p) }

# Scientific notation: 1.3e-06
fmt_sci <- function(x) {
  if (abs(x) < 1e-15) return("0")
  sprintf("%.1e", x)
}

make_triangle <- function(v) {
  path_name <- v$Path
  pair_text <- if(v$Pair=="") path_name else paste(path_name, v$Pair)
  treat_label <- if(v$Path=="Species") "Host species" else paste0("Area\n(", v$Pair, ")")

  s <- 3.5; h <- s*sqrt(3)/2

  nodes <- data.frame(
    name  = c(treat_label, "Bacterial\nMAG π", "Viral π"),
    color = c(if(v$Path=="Species")"#E8913A" else "#D46A6A", "#599EC4", "#79C36A"),
    x = c(0, s/2, s), y = c(0, h, 0),
    lx = c(-0.5, s/2, s+0.5),
    ly = c(-0.3, h+0.3, -0.3),
    lhjust = c(0.5,0.5,0.5), lvjust = c(1,0,1))

  r <- 0.33
  edges <- data.frame(
    x0=c(nodes$x[1],nodes$x[1],nodes$x[2]),
    y0=c(nodes$y[1],nodes$y[1],nodes$y[2]),
    x1=c(nodes$x[2],nodes$x[3],nodes$x[3]),
    y1=c(nodes$y[2],nodes$y[3],nodes$y[3]),
    color=c(if(v$Path=="Species")"#E8913A" else "#D46A6A","#4E937A","#7A68A6"))
  edges$dx <- edges$x1-edges$x0; edges$dy <- edges$y1-edges$y0
  edges$len <- sqrt(edges$dx^2+edges$dy^2)
  edges$x0 <- edges$x0 + r*edges$dx/edges$len
  edges$y0 <- edges$y0 + r*edges$dy/edges$len
  edges$x1 <- edges$x1 - r*edges$dx/edges$len
  edges$y1 <- edges$y1 - r*edges$dy/edges$len

  mx <- (edges$x0+edges$x1)/2; my <- (edges$y0+edges$y1)/2
  mx[1] <- mx[1] + 0.85  # ACME label → more centered
  my[2] <- my[2] + 0.40  # ADE label pushed up

  labels <- c(
    sprintf("VLP: ACME=%s, p=%s\nBulk: ACME=%s, p=%s",
            fmt_sci(v$VLP_ACME), fmt_p(v$VLP_ACME_p),
            fmt_sci(v$Bulk_ACME), fmt_p(v$Bulk_ACME_p)),
    sprintf("VLP: ADE=%s, p=%s\nBulk: ADE=%s, p=%s",
            fmt_sci(v$VLP_ADE), fmt_p(v$VLP_ADE_p),
            fmt_sci(v$Bulk_ADE), fmt_p(v$Bulk_ADE_p)),
    "")

  ggplot() +
    geom_segment(aes(x=x0,y=y0,xend=x1,yend=y1,color=color), data=edges, linewidth=1.2,
                 arrow=arrow(length=unit(0.2,"cm"),type="closed")) +
    geom_point(aes(x=x,y=y,color=color), data=nodes, size=9, shape=16) +
    geom_text(aes(x=lx,y=ly,label=name,hjust=lhjust,vjust=lvjust), data=nodes,
              size=5.5, fontface="bold", color="grey20", lineheight=0.9) +
    geom_label(data=edges[labels!="",], aes(x=mx[labels!=""], y=my[labels!=""],
              label=labels[labels!=""], color=color[labels!=""]), size=4.5, lineheight=1.15,
              fill="white", label.size=0.4, label.padding=unit(0.22,"cm"),
              fontface="plain", show.legend=FALSE) +
    scale_color_identity() +
    coord_fixed(xlim=c(-2.6,s+1.3), ylim=c(-1.5,h+1.3), expand=FALSE) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(face="bold",hjust=0.5,size=15),
          plot.margin=margin(5,5,5,5),
          panel.background=element_rect(fill="white",color=NA))
}

p1 <- make_triangle(vals[1,])
p2 <- make_triangle(vals[2,])
p3 <- make_triangle(vals[3,])
p4 <- make_triangle(vals[4,])

combined <- (p1 | p2) / (p3 | p4)
ggsave(file.path(fig_dir, "mediation_pi_all.pdf"), combined, width=14, height=12, dpi=300, device=cairo_pdf)
ggsave(file.path(fig_dir, "mediation_pi_all.png"), combined, width=14, height=12, dpi=300, type="cairo")
cat("Done: figures/mediation_pi_all.pdf\n")
