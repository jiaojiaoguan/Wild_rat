#!/usr/bin/env Rscript
# Single-panel mediation triangle: Species → Bacteria PC1 → Virus AMG PC1
# VLP + Bulk, for main text
library(ggplot2)

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
VLP_ACME   <- -0.2442;  VLP_ACME_p <- 0.0990
VLP_ADE    <-  0.2028;  VLP_ADE_p  <- 0.2416
Bulk_ACME  <-  0.3393;  Bulk_ACME_p <- 0.0088
Bulk_ADE   <- -0.1459;  Bulk_ADE_p  <- 0.3544

fmt_p <- function(p) { if(p<0.001) "<0.001" else sprintf("%.4f", p) }

s <- 3.5; h <- s*sqrt(3)/2

nodes <- data.frame(
  name  = c("Host species", "Bacterial PC1", "Virus AMGs PC1"),
  color = c("#E8913A", "#599EC4", "#79C36A"),
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
  color=c("#E8913A","#4E937A","#7A68A6"))
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
  sprintf("VLP: ACME=%.3f, p=%s\nBulk: ACME=%.3f, p=%s",
          abs(VLP_ACME), fmt_p(VLP_ACME_p),
          abs(Bulk_ACME), fmt_p(Bulk_ACME_p)),
  sprintf("VLP: ADE=%.3f, p=%s\nBulk: ADE=%.3f, p=%s",
          abs(VLP_ADE), fmt_p(VLP_ADE_p),
          abs(Bulk_ADE), fmt_p(Bulk_ADE_p)),
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
  coord_fixed(xlim=c(-1.9,s+1.6), ylim=c(-0.9,h+0.9), expand=FALSE) +
  theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold",hjust=0.5,size=15),
        plot.margin=margin(5,5,5,5),
        panel.background=element_rect(fill="white",color=NA))

ggsave(file.path(fig_dir, "mediation_AMG_species.pdf"), width=7, height=5.8, dpi=300, device=cairo_pdf)
ggsave(file.path(fig_dir, "mediation_AMG_species.png"), width=7, height=5.8, dpi=300, type="cairo")
cat("Done: figures/mediation_AMG_species.pdf\n")
