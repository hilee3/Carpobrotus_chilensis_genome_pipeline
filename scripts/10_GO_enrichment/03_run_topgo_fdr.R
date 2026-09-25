library(topGO)


# ============================================================
# Paths
# ============================================================

WORKDIR <- "/path/to/genome_comparison"

INDIR <- file.path(
  WORKDIR,
  "go_enrichment"
)

OUTDIR <- file.path(
  WORKDIR,
  "go_enrichment_fdr"
)

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# Input files
# ============================================================

study_file <- file.path(
  INDIR,
  "dnds_gt1_ds001_genes.txt"
)

background_file <- file.path(
  INDIR,
  "kaks_background_genes.txt"
)

gene2go_file <- file.path(
  INDIR,
  "gene2go.tsv"
)


# ============================================================
# Output files
# ============================================================

bp_out <- file.path(
  OUTDIR,
  "topGO_dnds_gt1_ds001_BP_results.tsv"
)

mf_out <- file.path(
  OUTDIR,
  "topGO_dnds_gt1_ds001_MF_results.tsv"
)

cc_out <- file.path(
  OUTDIR,
  "topGO_dnds_gt1_ds001_CC_results.tsv"
)

combined_out <- file.path(
  OUTDIR,
  "topGO_dnds_gt1_ds001_all_ontologies_results.tsv"
)


# ============================================================
# Read gene lists
# ============================================================

study_genes <- scan(
  study_file,
  what = character(),
  quiet = TRUE
)

background_genes <- scan(
  background_file,
  what = character(),
  quiet = TRUE
)


# ============================================================
# Read gene-to-GO mapping
# ============================================================

gene2go_df <- read.table(
  gene2go_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

colnames(
  gene2go_df
) <- c(
  "gene",
  "GO"
)

gene2go_df <- gene2go_df[
  gene2go_df$gene %in%
    background_genes,
]

gene2go_df <- gene2go_df[
  !is.na(
    gene2go_df$GO
  ) &
    gene2go_df$GO != "",
]

gene2GO <- split(
  gene2go_df$GO,
  gene2go_df$gene
)


# ============================================================
# Define study/background membership
# ============================================================

all_genes <- factor(
  as.integer(
    background_genes %in%
      study_genes
  )
)

names(
  all_genes
) <- background_genes


# ============================================================
# topGO analysis
# ============================================================

run_topgo <- function(
  ontology_name,
  node_size = 10
) {

  GOdata <- new(
    "topGOdata",
    ontology = ontology_name,
    allGenes = all_genes,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO,
    nodeSize = node_size
  )

  result <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )

  n_tested <- length(
    score(
      result
    )
  )

  result_table <- GenTable(
    GOdata,
    weight01_Fisher = result,
    orderBy = "weight01_Fisher",
    topNodes = n_tested,
    numChar = 1000
  )

  result_table$Ontology <- (
    ontology_name
  )

  result_table$weight01_Fisher_numeric <- as.numeric(
    sub(
      "^<\\s*",
      "",
      result_table$weight01_Fisher
    )
  )

  # Benjamini-Hochberg correction within each ontology.
  result_table$FDR_BH <- p.adjust(
    result_table$weight01_Fisher_numeric,
    method = "BH"
  )

  result_table$GeneRatio <- (
    result_table$Significant /
      result_table$Annotated
  )

  result_table$minus_log10_p <- -log10(
    result_table$weight01_Fisher_numeric
  )

  result_table$minus_log10_FDR_BH <- -log10(
    result_table$FDR_BH
  )

  result_table <- result_table[
    order(
      result_table$weight01_Fisher_numeric
    ),
  ]

  rownames(
    result_table
  ) <- NULL

  return(
    result_table
  )
}


# ============================================================
# Run BP, MF, and CC
# ============================================================

table_BP <- run_topgo(
  "BP",
  node_size = 10
)

table_MF <- run_topgo(
  "MF",
  node_size = 10
)

table_CC <- run_topgo(
  "CC",
  node_size = 10
)


# ============================================================
# Save results
# ============================================================

write.table(
  table_BP,
  file = bp_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  table_MF,
  file = mf_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  table_CC,
  file = cc_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


table_all <- rbind(
  table_BP,
  table_MF,
  table_CC
)

write.table(
  table_all,
  file = combined_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)