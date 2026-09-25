library(dplyr)
library(ggplot2)
library(readr)


# ============================================================
# Paths
# ============================================================

WORKDIR <- "/path/to/genome_comparison"

INDIR <- file.path(
  WORKDIR,
  "go_enrichment_fdr"
)

OUTDIR <- file.path(
  WORKDIR,
  "go_enrichment_fdr",
  "plots"
)

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# Read topGO results
# ============================================================

go_all <- read_tsv(
  file.path(
    INDIR,
    "topGO_dnds_gt1_ds001_all_ontologies_results.tsv"
  ),
  show_col_types = FALSE
)


# ============================================================
# Prepare plotting data
# ============================================================

go_all <- go_all %>%
  mutate(
    FoldEnrichment =
      Significant /
      Expected
  )


plot_df <- go_all %>%
  filter(
    weight01_Fisher_numeric <
      0.05
  )


# Restore full GO term names used in the publication figure.
plot_df <- plot_df %>%
  mutate(
    Term = case_when(
      GO.ID == "GO:0140534" ~
        "endoplasmic reticulum protein-containing complex",

      GO.ID == "GO:0042626" ~
        "ATPase-coupled transmembrane transporter activity",

      TRUE ~ Term
    )
  )


# Order terms within each ontology by fold enrichment.
plot_df <- plot_df %>%
  group_by(
    Ontology
  ) %>%
  arrange(
    FoldEnrichment,
    .by_group = TRUE
  ) %>%
  mutate(
    Term = factor(
      Term,
      levels = unique(
        Term
      )
    )
  ) %>%
  ungroup()


# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = FoldEnrichment,
    y = Term,
    size = Significant,
    color = minus_log10_p
  )
) +
  geom_point() +
  facet_grid(
    Ontology ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.05
      )
    )
  ) +
  scale_size_continuous(
    range = c(
      2.5,
      9
    ),
    breaks = sort(
      unique(
        plot_df$Significant
      )
    )
  ) +
  scale_color_gradient(
    low = "#1f4e79",
    high = "#5dade2"
  ) +
  labs(
    x = "Fold enrichment",
    y = NULL,
    color = expression(
      -log[10](p)
    ),
    size = "Gene count"
  ) +
  theme_bw(
    base_size = 14,
    base_family = "sans"
  ) +
  theme(
    text = element_text(
      family = "sans"
    ),

    strip.background = element_rect(
      fill = "grey90"
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    axis.text.y = element_text(
      size = 12
    ),

    axis.text.x = element_text(
      size = 12
    ),

    legend.title = element_text(
      size = 13
    ),

    legend.text = element_text(
      size = 11
    )
  )


# ============================================================
# Save figure
# ============================================================

ggsave(
  file.path(
    OUTDIR,
    "GO_fold_enrichment_dotplot.pdf"
  ),
  p,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "GO_fold_enrichment_dotplot.png"
  ),
  p,
  width = 10,
  height = 7,
  dpi = 300
)