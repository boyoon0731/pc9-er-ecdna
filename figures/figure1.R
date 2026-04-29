# ==============================================================================
# Script:  figure.R
# Purpose: Reproduce figures for LCPC9 ER dataset
#          (genome-wide CNV difference, ecDNA reconstruction,
#           allele-specific CN, expression rank plots,
#           and oncogene CN scatter plot)
#
# Usage:
#   1. Place input files under data/ (see "Input files" below).
#   2. Run the entire script: source("figure.R")
#      or execute section by section in RStudio.
#   3. Figures are saved as PDF/PNG in output/ (created automatically).
#
# Input files (data/):
#   LCPC9_ER_CNV_somatic_purple.tsv  – segment-level copy number 
#   LCPC9_ER_CNV_gene_purple.tsv     – gene-level copy number 
#   LCPC9_ER_SNV_mutect2.rds         – somatic SNVs 
#   LCPC9_ER_RNA_tpm.rds             – gene expression matrix (TPM; rows = genes)
#   LCPC9_ER_SV_bedpe.rds            – structural variants (BEDPE format)
#
# External reference (downloaded automatically at runtime):
#   OncoKB Cancer Gene List  https://www.oncokb.org/api/v1/utils/cancerGeneList.txt
#   (no account required for this endpoint)
#
# Dependencies:
#   install.packages(c("tidyverse", "ggpubr", "ggrepel", "patchwork",
#                      "officer", "rvg"))
#   BiocManager::install(c("GenomicRanges", "GenomeInfoDb"))
# ==============================================================================

# ── 0. Packages ────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(GenomicRanges)
library(GenomeInfoDb)
library(officer)
library(rvg)

# ── 1. Global settings ─────────────────────────────────────────────────────────
SAMPLES       <- c("ER-S", "ER-R")            # parental / resistant
TARGET_GENES  <- c("FANCD2", "PPARG", "RAF1", "VHL")  # ecDNA oncogenes
COLORS_SAMPLE <- c("ER-S" = "#1f77b4", "ER-R" = "#d62728")
COLORS_ECDNA  <- c("TRUE" = "#ff7f0e", "FALSE" = "#1f77b4")

# Create output directory if absent
dir.create("output", showWarnings = FALSE)

# ── 2. Load data ───────────────────────────────────────────────────────────────
cn_somatic <- read_tsv("data/LCPC9_ER_CNV_somatic_purple.tsv", show_col_types = FALSE)
cn_gene    <- read_tsv("data/LCPC9_ER_CNV_gene_purple.tsv",    show_col_types = FALSE)

mut_mutect2_ss <- readRDS("data/LCPC9_ER_SNV_mutect2.rds")
mut_mutect2_ss_filtered <- mut_mutect2_ss %>%
  filter(t_depth >= 30, t_alt_count >= 3, t_alt_count / t_depth >= 0.1)

gene_exp <- readRDS("data/LCPC9_ER_RNA_tpm.rds")
sv_bedpe <- readRDS("data/LCPC9_ER_SV_bedpe.rds")

# OncoKB cancer gene list (downloaded at runtime; no login required)
# Alternatively, place a local copy at data/cancerGeneList.tsv to work offline.
ONCOKB_URL      <- "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt"
ONCOKB_LOCAL    <- "data/cancerGeneList.tsv"

if (!file.exists(ONCOKB_LOCAL)) {
  message("Downloading OncoKB cancer gene list...")
  download.file(ONCOKB_URL, destfile = ONCOKB_LOCAL, quiet = TRUE)
}
oncokb_cancer_genes <- read_tsv(ONCOKB_LOCAL, show_col_types = FALSE)


# ==============================================================================
# Figure a1 – Genome-wide Copy Number Difference (ER-R minus ER-S)
# ==============================================================================

# ── Build disjoint segments with copy-number difference ────────────────────────
make_cn_diff <- function(cn_somatic, sample_parental = "ER-S", sample_resistant = "ER-R") {
  cn_par <- filter(cn_somatic, sample_id == sample_parental)
  cn_res <- filter(cn_somatic, sample_id == sample_resistant)

  gr_par <- GRanges(seqnames = cn_par$chromosome,
                    ranges   = IRanges(cn_par$start, cn_par$end),
                    copyNumber = cn_par$copyNumber)
  gr_res <- GRanges(seqnames = cn_res$chromosome,
                    ranges   = IRanges(cn_res$start, cn_res$end),
                    copyNumber = cn_res$copyNumber)

  gr_comb <- disjoin(c(gr_par, gr_res))
  gr_comb$cn_par <- 0
  gr_comb$cn_res <- 0

  hits_par <- findOverlaps(gr_comb, gr_par, type = "within")
  hits_res <- findOverlaps(gr_comb, gr_res, type = "within")
  gr_comb$cn_par[queryHits(hits_par)] <- gr_par$copyNumber[subjectHits(hits_par)]
  gr_comb$cn_res[queryHits(hits_res)] <- gr_res$copyNumber[subjectHits(hits_res)]

  as.data.frame(gr_comb) %>%
    filter(!seqnames %in% c("X", "Y")) %>%
    mutate(chr            = factor(seqnames, levels = unique(seqnames)),
           copyNumber_diff = cn_res - cn_par)
}

cnv_diff <- make_cn_diff(cn_somatic)

# Minimum visible bar width (4 % of chromosome range)
cnv_diff_vis <- cnv_diff %>%
  group_by(chr) %>%
  mutate(min_width = (max(end) - min(start)) * 0.04,
         plot_end  = pmax(end, start + min_width)) %>%
  ungroup()

# ── Plot ───────────────────────────────────────────────────────────────────────
plot_cn_diff_all_chr <- ggplot(cnv_diff_vis) +
  geom_rect(aes(xmin = start, xmax = plot_end,
                ymin = 0,     ymax = copyNumber_diff,
                fill = copyNumber_diff > 0),
            linewidth = 0) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  facet_grid(. ~ chr, scales = "free_x", switch = "x") +
  scale_fill_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#2471A3"),
    labels = c("TRUE" = "Gain (ER-R vs ER-S)", "FALSE" = "Loss (ER-R vs ER-S)"),
    breaks = c("TRUE", "FALSE")
  ) +
  scale_y_continuous(expand   = c(0, 0),
                     breaks   = seq(-25, 25, 5),
                     labels   = function(x) ifelse(x == 0, "0", x)) +
  coord_cartesian(ylim = c(-25, 25)) +
  labs(y = "Copy Number Change\n(ER-R − ER-S)") +
  theme_classic(base_family = "Arial", base_size = 20) +
  theme(
    strip.background  = element_rect(fill = NA, color = NA),
    strip.placement   = "outside",
    axis.title.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_line(color = "black", linewidth = 0.4),
    axis.line.y       = element_line(color = "black", linewidth = 0.4),
    axis.line.x       = element_blank(),
    panel.spacing.x   = unit(0, "lines"),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
    panel.grid        = element_blank(),
    legend.position   = "top",
    legend.title      = element_blank()
  )


# ==============================================================================
# Figure a2 – ecDNA Reconstruction (SV arcs + copy-number track)
# ==============================================================================

# ── ecDNA interval (chromosome 3) ─────────────────────────────────────────────
amplicon_regions <- data.frame(chr   = "chr3",
                               start = 2970012,
                               end   = 27250043)

# Linear-coordinate helpers
amplicon_regions_linear <- amplicon_regions %>%
  mutate(width        = end - start,
         linear_start = 0,
         linear_end   = width)

x_range <- c(0, max(amplicon_regions_linear$linear_end))

# Map genomic position → linear x coordinate
to_linear <- function(chrom, position) {
  map2_dbl(chrom, position, function(ch, pos) {
    idx <- which(amplicon_regions_linear$chr   == ch &
                 amplicon_regions_linear$start <= pos &
                 amplicon_regions_linear$end   >= pos)
    if (length(idx) == 0) return(NA_real_)
    amplicon_regions_linear$linear_start[idx[1]] +
      (pos - amplicon_regions_linear$start[idx[1]])
  })
}

# Clip copy-number segments to amplicon region
clip_cn_to_regions <- function(cn_df, regions) {
  map_dfr(seq_len(nrow(regions)), function(i) {
    cn_df %>%
      filter(chr == regions$chr[i],
             start <= regions$end[i],
             end   >= regions$start[i]) %>%
      mutate(start = pmax(start, regions$start[i]),
             end   = pmin(end,   regions$end[i]))
  }) %>%
  distinct()
}

# SV display parameters
DIST_THRESHOLD <- 5e5
MAX_CN         <- 25

# ── Panel-building function ────────────────────────────────────────────────────
make_reconplot <- function(sample_id) {
  # Structural variants
  sv_sample <- sv_bedpe %>%
    filter(sample_id == !!sample_id) %>%
    mutate(chr1 = paste0("chr", chr1), chr2 = paste0("chr", chr2))

  sv_region <- sv_sample %>%
    filter(chr1 == amplicon_regions$chr,
           pos1 >= amplicon_regions$start,
           pos1 <= amplicon_regions$end)

  sv_linear <- sv_region %>%
    mutate(l_pos1 = to_linear(chr1, pos1),
           l_pos2 = to_linear(chr2, pos2)) %>%
    filter(!is.na(l_pos1), !is.na(l_pos2)) %>%
    mutate(dist  = abs(l_pos2 - l_pos1),
           type  = if_else(dist >= DIST_THRESHOLD, "curve", "triangle"),
           mid_x = (l_pos1 + l_pos2) / 2,
           height = 0.5)

  # Copy number
  cn_sample <- cn_somatic %>%
    filter(sample_id == !!sample_id) %>%
    mutate(chr = paste0("chr", chromosome))

  cn_linear <- clip_cn_to_regions(cn_sample, amplicon_regions) %>%
    mutate(l_start = to_linear(chr, start),
           l_end   = to_linear(chr, end))

  # SV panel
  p_sv <- ggplot() +
    geom_curve(data    = filter(sv_linear, type == "curve"),
               aes(x = l_pos1, xend = l_pos2, y = 0, yend = 0),
               curvature = -0.2, color = "#4361EE") +
    geom_segment(data = filter(sv_linear, type == "triangle"),
                 aes(x = l_pos1, xend = mid_x,  y = 0,      yend = height),
                 color = "#4361EE") +
    geom_segment(data = filter(sv_linear, type == "triangle"),
                 aes(x = mid_x,  xend = l_pos2, y = height, yend = 0),
                 color = "#4361EE") +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = c(0, 0.6)) +
    labs(y = "SV") +
    theme_classic() +
    theme(axis.title.x     = element_blank(),
          axis.line        = element_blank(),
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          panel.background = element_rect(fill = "grey98"))

  # Copy-number panel
  p_cn <- ggplot() +
    geom_rect(data = cn_linear,
              aes(xmin = l_start, xmax = l_end,
                  ymin = minorAlleleCopyNumber, ymax = copyNumber),
              fill = "#C1121F", alpha = 0.4) +
    geom_rect(data = cn_linear,
              aes(xmin = l_start, xmax = l_end,
                  ymin = 0, ymax = minorAlleleCopyNumber),
              fill = "#C1121F", alpha = 0.6, color = NA) +
    geom_rect(data = cn_linear,
              aes(xmin = l_start, xmax = l_end,
                  ymin = copyNumber, ymax = copyNumber),
              color = "black", linewidth = 1.2) +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = c(0, MAX_CN)) +
    labs(x = "Chromosome 3", y = "Copy Number") +
    theme_classic()

  p_sv / p_cn +
    plot_layout(heights = c(0.7, 2.5)) +
    plot_annotation(title = sample_id)
}

plot_recon_ERS <- make_reconplot("ER-S")
plot_recon_ERR <- make_reconplot("ER-R")
plot_recon_combined <- wrap_plots(plot_recon_ERS, plot_recon_ERR, ncol = 1)


# ==============================================================================
# Figure a4 – ecDNA Oncogene Copy Number
# ==============================================================================

# ── ecDNA Oncogene Copy Number: ER-S vs ER-R Scatter Plot ─────────────────────
# Points coloured by ecDNA status; regression lines per group;
# ecDNA gene labels added with ggrepel.
# P-value annotation (Wilcoxon / t-test result from the original analysis).
cn_oncogene_df <- cn_gene %>%
  filter(gene %in% oncokb_cancer_genes$`Hugo Symbol`) %>%
  select(gene, minCopyNumber, sample_id) %>%
  pivot_wider(names_from  = sample_id,
              values_from = minCopyNumber) %>%
  mutate(in_ecDNA = factor(gene %in% TARGET_GENES,
                           levels = c(TRUE, FALSE)))

plot_CN_oncogene <- ggplot(cn_oncogene_df,
                           aes(x = `ER-S`, y = `ER-R`, color = in_ecDNA)) +
  geom_point(size = 1.4) +
  geom_smooth(aes(group = in_ecDNA),
              method    = "lm",
              se        = FALSE,
              linetype  = "solid",
              linewidth = 0.5,
              fullrange = TRUE) +
  geom_text_repel(
    data          = filter(cn_oncogene_df, in_ecDNA == TRUE),
    aes(label     = gene),
    color         = "#ff7f0e",
    size          = 7,
    box.padding   = 0.5,
    point.padding = 0.5,
    max.overlaps  = 5
  ) +
  annotate("text", x = 15, y = 10,
           label  = expression(italic(P) == 2.8 %*% 10^{-4}),
           size   = 6, hjust = 0, family = "Arial") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_color_manual(values = c("TRUE" = "#ff7f0e", "FALSE" = "#1f77b4"),
                     labels = c("TRUE" = "ecDNA", "FALSE" = "Non-ecDNA")) +
  labs(x = "Copy Number in ER-S",
       y = "Copy Number in ER-R") +
  coord_fixed() +
  theme_classic(base_family = "Arial", base_size = 22) +
  theme(
    plot.title        = element_blank(),
    panel.grid        = element_blank(),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 20),
    legend.position      = c(0.5, 0.1),
    legend.justification = c(0, 0),
    legend.background = element_rect(color = "black", fill = "white",
                                     linewidth = 0.5),
    legend.margin     = margin(5, 8, 5, 8)
  )

# ── Allele-specific Copy Number of ecDNA Oncogenes ────────────────────────────

cn_gene_ecDNA <- cn_gene %>%
  filter(gene %in% TARGET_GENES) %>%
  transmute(sample_id,
            gene,
            totalCN = minCopyNumber,
            minorCN = minMinorAlleleCopyNumber) %>%
  mutate(gene      = factor(gene,      levels = TARGET_GENES),
         sample_id = factor(sample_id, levels = SAMPLES)) %>%
  pivot_longer(cols      = c(totalCN, minorCN),
               names_to  = "AlleleCN",
               values_to = "CN")

plot_alleleCN_ecDNA_oncogene <- ggplot() +
  geom_bar(
    data = filter(cn_gene_ecDNA, AlleleCN != "majorCN"),
    aes(x     = sample_id,
        y     = CN,
        fill  = sample_id,
        alpha = if_else(AlleleCN == "minorCN", 1, 0.5)),
    stat = "identity", position = "identity", width = 0.8
  ) +
  # CN labels – ER-S (right of bar) and ER-R (left of bar)
  geom_text(
    data = filter(cn_gene_ecDNA, AlleleCN == "totalCN", sample_id == "ER-S"),
    aes(x = sample_id, y = CN, label = round(CN, 2)),
    color = "black", size = 5.5, hjust = -0.2, angle = 90
  ) +
  geom_text(
    data = filter(cn_gene_ecDNA, AlleleCN == "totalCN", sample_id == "ER-R"),
    aes(x = sample_id, y = CN, label = round(CN, 2)),
    color = "black", size = 5.5, hjust = 1.2, angle = 90
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_grid(gene ~ ., switch = "y") +
  scale_alpha_identity() +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0),
                     breaks = seq(0, 30, 10)) +
  scale_fill_manual(values = COLORS_SAMPLE) +
  labs(y = "Copy Number (minor allele CN + major allele CN)") +
  theme_classic(base_family = "Arial", base_size = 18) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 20, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    strip.text.y     = element_text(size = 18, face = "italic", color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.placement  = "outside",
    panel.spacing    = unit(1, "lines"),
    panel.grid       = element_blank(),
    legend.position  = "none"
  )


# ==============================================================================
# Figure b1 – Genome-wide Expression Rank + ecDNA Oncogene Expression
# ==============================================================================

# ── Rank plot (all expressed genes) ───────────────────────────────────────────
df_rank_plot <- gene_exp %>%
  rownames_to_column("gene") %>%
  filter(`ER-S` > 0, `ER-R` > 0) %>%
  mutate(is_ecDNA = factor(gene %in% TARGET_GENES, levels = c(TRUE, FALSE)),
         log2FC   = log2((`ER-R` + 0.5) / (`ER-S` + 0.5))) %>%
  arrange(log2FC) %>%
  mutate(gene = factor(gene, levels = gene))

plot_rank_expression <- ggplot(df_rank_plot,
                               aes(x = gene, y = log2FC, fill = is_ecDNA)) +
  geom_col(width = 1) +
  geom_col(data  = filter(df_rank_plot, is_ecDNA == TRUE), width = 40) +
  geom_text(data = filter(df_rank_plot, is_ecDNA == TRUE),
            aes(y = -5.8, label = gene),
            hjust = 0, size = 5, color = "#ff7f0e", family = "Arial") +
  geom_hline(yintercept = 3.5, linetype = "dashed",
             color = "#ff7f0e", linewidth = 0.8) +
  annotate("text", x = 1, y = 3.6, label = "3.5",
           color = "#ff7f0e", size = 5, vjust = 0, family = "Arial") +
  scale_fill_manual(values = COLORS_ECDNA,
                    labels = c("ecDNA", "non-ecDNA")) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) +
  labs(y = "Log2FC (ER-R/ER-S)", x = "Ranked Genes") +
  coord_flip() +
  theme_classic(base_family = "Arial", base_size = 22) +
  theme(
    axis.text.y          = element_blank(),
    axis.ticks.y         = element_blank(),
    legend.title         = element_blank(),
    legend.position      = c(0.02, 0.05),
    legend.justification = c(0, 0),
    legend.background    = element_rect(color = "black", linewidth = 0.5),
    legend.margin        = margin(4, 8, 4, 8)
  )

# ── ecDNA oncogene TPM bar plot ────────────────────────────────────────────────
df_ecdna_exp <- gene_exp %>%
  rownames_to_column("gene_name") %>%
  filter(gene_name %in% TARGET_GENES) %>%
  pivot_longer(cols      = all_of(SAMPLES),
               names_to  = "sample",
               values_to = "expression") %>%
  mutate(sample = factor(sample, levels = SAMPLES))

plot_ecdna_expression <- ggplot(df_ecdna_exp,
                                aes(x = sample, y = expression, fill = sample)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = round(expression, 1),
                hjust = if_else(expression >= 250, 1.1, -0.2)),
            color = "black", size = 5.5, angle = 90) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_grid(gene_name ~ ., switch = "y") +
  scale_y_continuous(limits = c(0, 450), expand = c(0, 0),
                     breaks = seq(0, 400, 100)) +
  scale_fill_manual(values = COLORS_SAMPLE) +
  labs(y = "mRNA Expression (TPM)") +
  theme_classic(base_family = "Arial", base_size = 18) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 20, color = "black",
                                    angle = 90, vjust = 0.5),
    axis.text.y      = element_text(size = 16),
    strip.text.y     = element_text(face = "italic", color = "black", size = 18),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )



# ==============================================================================
# Save all figures
# ==============================================================================
ggsave("output/fig_a1_cnv_diff.pdf",
       plot_cn_diff_all_chr, width = 20, height = 4)

ggsave("output/fig_a2_recon_combined.pdf",
       plot_recon_combined,  width = 12, height = 10)

ggsave("output/fig_a4_oncogene_cn_scatter.pdf",
       plot_CN_oncogene,             width = 6, height = 6)

ggsave("output/fig_a4_allele_cn.pdf",
       plot_alleleCN_ecDNA_oncogene, width = 5, height = 10)

ggsave("output/fig_b1_rank_expression.pdf",
       plot_rank_expression,         width = 6, height = 10)

ggsave("output/fig_b1_ecdna_expression.pdf",
       plot_ecdna_expression,        width = 4, height = 10)



message("Done. Figures saved to output/")
