# Wild rat gut virome R code for analysis

R analysis scripts for the wild rat gut virome manuscript (VLP-enriched vs. Bulk
metagenomes). All scripts read their direct input files from `data/` and write
results to `results/`, figures to `figures/`, or the script's own folder.

## Directory layout

```
github_code/
├── data/                           # direct input files (see "Input data" below)
├── 01_alpha_diversity/             # Shannon diversity boxplots + drivers
├── 02_beta_diversity/              # PCoA + PERMANOVA + Procrustes
├── 03_dbRDA_variance_partitioning/ # db-RDA + variance partitioning
├── 04_microdiversity/              # nucleotide diversity (π) + variance partitioning
├── 05_mediation/                   # mediation analyses (viral π, viral AMG PC1)
├── AMGs_lifestyle/                 # AMG lifestyle (temperate/virulent) alpha + beta diversity
├── figures/                        # figure outputs (created at runtime)
├── results/                        # table / intermediate outputs (created at runtime)
└── README.md
```

## Scripts

### 01_alpha_diversity
| Script | What it does | Inputs (in `data/`) | Outputs |
|---|---|---|---|
| `plot_shannon.R` | Shannon boxplot, VLPs + Bulk side by side; Kruskal–Wallis + pairwise Wilcoxon | `wild_rat_metadata.csv`, `VLPs_virus_alpha_diversity.tsv`, `bulk_virus_alpha_diversity.tsv` | `figures/shannon_7box_without_anda.pdf/.png` + 4 pairwise TSVs (`results/`) |
| `analyze_shannon_drivers.R` | Type-II nested linear models for viral Shannon: base `Shannon ~ Species + Area` vs. `+ Bac_Shannon` (MAG Shannon diversity) | `wild_rat_metadata.csv`, `VLPs_virus_alpha_diversity.tsv`, `bulk_virus_alpha_diversity.tsv`, `explanatory_variables.csv` | console |

### 02_beta_diversity
| Script | What it does | Inputs | Outputs |
|---|---|---|---|
| `beta_2panel.R` | VLPs / Bulk PCoA (Bray–Curtis) with PERMANOVA (adonis2) + Procrustes alignment | `53vlps_virus_map_VLPs_sample_RPKM_normalized_data.tsv`, `53bulk_virus_map_bulk_sample_RPKM_normalized_data.tsv`, `wild_rat_metadata.csv` | `figures/virus_beta_2panel.pdf/.png` |

### 03_dbRDA_variance_partitioning
| Script | Model | Inputs | Outputs |
|---|---|---|---|
| `VLPs_species_area.R` | Model 0 — VLPs ~ Species + Area | VLPs virus map, `wild_rat_metadata.csv`, `explanatory_variables.csv` | `results/VLPs_RDA_results_species_area/` |
| `Bulks_species_area.R` | Model 0 — Bulk ~ Species + Area | Bulk virus map, metadata, `explanatory_variables.csv` | `results/Bulk_RDA_results_species_area/` |
| `VLPs_run_partial_RDA_species_area_bacpcoa.R` | VLPs ~ Species + Area + Bac_PCoA1-3 | VLPs virus map, metadata, `explanatory_variables.csv` | `results/VLPs_RDA_results_species_area_bacpcoa/` (tables) + `figures/VLPs_vpa_venn_*.pdf/.png`, `figures/VLPs_rda_triplot.pdf/.png` |
| `Bulk_run_partial_RDA_species_area_bacpcoa.R` | Bulk ~ Species + Area + Bac_PCoA1-3 | Bulk virus map, metadata, `explanatory_variables.csv` | `results/Bulk_RDA_results_species_area_bacpcoa/` (tables) + `figures/Bulk_vpa_venn_*.pdf/.png`, `figures/Bulk_rda_triplot.pdf/.png` |
| `plot_combined_triplots.R` | VLPs + Bulk db-RDA triplots side by side | both virus maps, metadata, `explanatory_variables.csv` | `figures/rda_triplots_combined_v2.pdf/.png` |

All db-RDA scripts use Bray–Curtis distance, `dbrda()` + `varpart()`, and the
Species + Area + Bacteria model reads the shared bacterial PCoA axes
(`Bac_PCoA1-3`) from `explanatory_variables.csv`.

### 04_microdiversity
| Script | What it does | Inputs | Outputs |
|---|---|---|---|
| `plot_pi.R` | π boxplot, VLPs + Bulk | `wild_rat_metadata.csv`, `VLPs_global_sample_microdiversity.csv`, `bulk_global_sample_microdiversity.csv` | `figures/pi_by_species_area_7box.pdf/.png` + 4 pairwise TSVs (`results/`) |
| `analyze_pi_drivers.R` | π drivers — nested type-II LM: `avg_pi ~ Species + Area` vs `+ MAG_pi` (bacterial MAG mean π), without R. andamanensis (n=46) | metadata, microdiversity CSVs, `MAGs_nucl_diversity_matrix.tsv` | console |
| `plot_variance_partitioning.R` | Variance-partitioning barplot for π and Shannon (values hardcoded from the db-RDA / model fits) | none | `figures/variance_partitioning_pi.pdf/.png`, `figures/variance_partitioning_shannon.pdf/.png` |

### 05_mediation
| Script | What it does | Inputs | Outputs |
|---|---|---|---|
| `run_all_mediation_pi.R` | Mediation: Species / Area → bacterial MAG π → viral π (10000 bootstrap, FDR) | `wild_rat_metadata.csv`, `MAGs_nucl_diversity_matrix.tsv`, `VLPs/bulk_global_sample_microdiversity.csv` | `results/mediation_pi_all_results.tsv` + console |
| `plot_all_mediation_pi.R` | 4-panel π mediation triangle (values hardcoded) | none | `figures/mediation_pi_all.pdf/.png` |
| `run_all_mediation_AMG.R` | Mediation: Species / Area → bacterial PC1 → viral AMG PC1 (10000 bootstrap, FDR) | `wild_rat_metadata.csv`, `explanatory_variables.csv`, `53vlp/meta_AMG_abundance_KO_normalized.tsv` | `results/mediation_AMG_all_results.tsv` + console |
| `plot_mediation_AMG_species.R` | Species mediation triangle (values hardcoded) | none | `figures/mediation_AMG_species.pdf/.png` |
| `plot_mediation_AMG_area_3panel.R` | 3 Area-pair mediation triangles (values hardcoded) | none | `figures/mediation_AMG_area_3panel.pdf/.png` |

### AMGs_lifestyle
| Script | What it does | Inputs | Outputs |
|---|---|---|---|
| `plot_amg_diversity.R` | AMG alpha diversity (richness + Shannon) boxplots, temperate vs virulent, VLP + Bulk; paired Wilcoxon | `53vlp_AMG_diversity.tsv`, `53meta_AMG_diversity.tsv` | `figures/AMG_diversity_{VLP,Bulk}_{richness,shannon}.pdf` |
| `analyze_lifestyle_betadiv.R` | VLP AMG KO beta diversity — pseudo-sample (temperate/virulent) paired PERMANOVA + PCoA | `VLP_lifestyle_KO_abundance.tsv` | `figures/VLP_lifestyle_KO_betadiv_pcoa.pdf/.png` |
| `analyze_lifestyle_betadiv_bulk.R` | Bulk AMG KO beta diversity — pseudo-sample paired PERMANOVA + PCoA | `Bulk_lifestyle_KO_abundance.tsv` | `figures/Bulk_lifestyle_KO_betadiv_pcoa.pdf/.png` |

## Input data (`data/`)

| File | Description | Source (upstream step) |
|---|---|---|
| `wild_rat_metadata.csv` | Sample metadata: Species, Area, Age_stage, Sex, etc. | field metadata (`Bacteria_abundance/`) |
| `53vlps_virus_map_VLPs_sample_RPKM_normalized_data.tsv` | VLPs viral abundance matrix (RPKM-normalized) | `virus_Taxa/` |
| `53bulk_virus_map_bulk_sample_RPKM_normalized_data.tsv` | Bulk viral abundance matrix (RPKM-normalized) | `virus_Taxa/` |
| `VLPs_virus_alpha_diversity.tsv` | VLPs sample alpha diversity (richness, Shannon) | `votu_abudance/` |
| `bulk_virus_alpha_diversity.tsv` | Bulk sample alpha diversity (richness, Shannon) | `votu_abudance/` |
| `explanatory_variables.csv` | Explanatory variables: Species, Area, Age_stage, Sex, Shannon, and shared bacterial PCoA axes `Bac_PCoA1-3` | threeway partial RDA (`RDA_analysis/RDA_results_threeway_partial/shared/`) |
| `VLPs_global_sample_microdiversity.csv` | VLPs nucleotide diversity (π) per sample | `Micro_diveristy/` |
| `bulk_global_sample_microdiversity.csv` | Bulk nucleotide diversity (π) per sample | `Micro_diveristy/` |
| `MAGs_nucl_diversity_matrix.tsv` | Bacterial MAG nucleotide diversity (π) matrix (rows = MAGs, cols = samples) | `Micro_diveristy/` |
| `53vlp_AMG_abundance_KO_normalized.tsv` | VLPs AMG KO abundance (normalized) | `AMG_driver/` |
| `53meta_AMG_abundance_KO_normalized.tsv` | Bulk AMG KO abundance (normalized) | `AMG_driver/` |
| `53vlp_AMG_diversity.tsv` | VLPs AMG alpha diversity (virulent/temperate richness + Shannon) | `AMG_0511/` |
| `53meta_AMG_diversity.tsv` | Bulk AMG alpha diversity (virulent/temperate richness + Shannon) | `AMG_0511/` |
| `VLP_lifestyle_KO_abundance.tsv` | VLP AMG KO abundance split by lifestyle (rows = `KO\|lifestyle`) | `AMG_driver/` |
| `Bulk_lifestyle_KO_abundance.tsv` | Bulk AMG KO abundance split by lifestyle (rows = `KO\|lifestyle`) | `AMG_driver/` |

## Required R packages

`ggplot2`, `dplyr`, `tidyr`, `readr`, `ggpubr`, `patchwork`, `vegan`, `car`, `mediation`
(`rstudioapi` is used only to locate the script when run interactively and is
optional — scripts fall back to `commandArgs()`/`getwd()`).

## Running

```r
# from RStudio: open a script and source it, or
Rscript 01_alpha_diversity/plot_shannon.R
```

## Note on model numbering

In the manuscript the bacteria-containing db-RDA is referred to as **Model 1**,
while the script comments label it **"Model 2"** (Species + Area is "Model 0").
The analysis itself is identical; only the comment label differs. Adjust the
comment if you want it to match the manuscript numbering.
