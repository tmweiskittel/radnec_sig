library(data.table)
library(dplyr)
library(readr)
library(tibble)
library(DSS)
library(bsseq)

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------

dir.create(
  "results/differential_methylation",
  recursive = TRUE,
  showWarnings = FALSE
)

prevalence_threshold <- 0.50

# -------------------------------------------------------------------------
# Load comparison metadata
# -------------------------------------------------------------------------

comparison_1 <- read_tsv(
  "results/comparisons/radnec_vs_recurrent.tsv",
  show_col_types = FALSE
)

comparison_2 <- read_tsv(
  "results/comparisons/necrosis_type.tsv",
  show_col_types = FALSE
)

stopifnot(
  all(comparison_1$surgery == 2),
  all(comparison_2$surgery == 2),
  !anyDuplicated(comparison_1$id),
  !anyDuplicated(comparison_2$id)
)

# -------------------------------------------------------------------------
# Load RRBS data
# -------------------------------------------------------------------------

load("data/klughammer/rrbsCg.RData")

stopifnot(
  exists("ret"),
  is.data.table(ret),
  all(c(
    "id",
    "chr",
    "start",
    "hitCount",
    "readCount"
  ) %in% names(ret))
)

comparison_ids <- unique(c(
  comparison_1$id,
  comparison_2$id
))

rrbs <- ret[
  id %in% comparison_ids,
  .(
    id,
    chr,
    start,
    hitCount,
    readCount
  )
]

rm(ret)
gc()

# Verify assumptions about supplied rrbsCg object
stopifnot(
  all(rrbs$hitCount >= 0),
  all(rrbs$readCount >= 10),
  all(rrbs$hitCount <= rrbs$readCount)
)

cat(
  "\nRRBS analysis object\n",
  "Rows: ", format(nrow(rrbs), big.mark = ","), "\n",
  "Samples: ", uniqueN(rrbs$id), "\n",
  sep = ""
)

# -------------------------------------------------------------------------
# Helper: identify CpGs with >=50% prevalence in every group
# -------------------------------------------------------------------------

get_eligible_cpgs <- function(
  rrbs_data,
  metadata,
  group_col,
  threshold = 0.50
) {

  metadata <- as.data.table(metadata)

  group_sizes <- metadata[
    ,
    .(
      group_n = uniqueN(id)
    ),
    by = group_col
  ]

  dat <- rrbs_data[
    id %in% metadata$id
  ]

  dat[
    metadata,
    group := get(paste0("i.", group_col)),
    on = "id"
  ]

  prevalence <- dat[
    ,
    .(
      n_samples = uniqueN(id)
    ),
    by = .(
      chr,
      start,
      group
    )
  ]

  prevalence <- merge(
    prevalence,
    group_sizes,
    by.x = "group",
    by.y = group_col
  )

  prevalence[
    ,
    sample_prevalence := n_samples / group_n
  ]

  n_groups <- uniqueN(metadata[[group_col]])

  eligible <- prevalence[
    sample_prevalence >= threshold,
    .(
      n_groups_passing = uniqueN(group)
    ),
    by = .(
      chr,
      start
    )
  ][
    n_groups_passing == n_groups,
    .(
      chr,
      start
    )
  ]

  list(
    eligible = eligible,
    prevalence = prevalence,
    group_sizes = group_sizes
  )
}

# -------------------------------------------------------------------------
# Helper: create a BSseq object for two groups
# -------------------------------------------------------------------------

make_bsseq <- function(
  rrbs_data,
  metadata,
  group_col,
  group_1,
  group_2,
  eligible_cpgs
) {

  metadata <- as.data.table(metadata)

  metadata <- metadata[
    get(group_col) %in% c(
      group_1,
      group_2
    )
  ]

  sample_ids <- metadata$id

  dat <- rrbs_data[
    id %in% sample_ids
  ]

  dat <- merge(
    dat,
    eligible_cpgs,
    by = c(
      "chr",
      "start"
    )
  )

  # CpG rows in stable genomic order
  loci <- unique(
    eligible_cpgs[
      order(
        chr,
        start
      )
    ]
  )

  loci[
    ,
    locus_id := paste(
      chr,
      start,
      sep = ":"
    )
  ]

  dat[
    ,
    locus_id := paste(
      chr,
      start,
      sep = ":"
    )
  ]

  # M = methylated reads
  methylated <- dcast(
    dat,
    locus_id ~ id,
    value.var = "hitCount"
  )

  # Cov = total reads
  coverage <- dcast(
    dat,
    locus_id ~ id,
    value.var = "readCount"
  )

  # Ensure all eligible loci are represented
  methylated <- merge(
    loci[
      ,
      .(
        locus_id
      )
    ],
    methylated,
    by = "locus_id",
    all.x = TRUE
  )

  coverage <- merge(
    loci[
      ,
      .(
        locus_id
      )
    ],
    coverage,
    by = "locus_id",
    all.x = TRUE
  )

  # Put columns in metadata order
  methylated <- methylated[
    ,
    c(
      "locus_id",
      sample_ids
    ),
    with = FALSE
  ]

  coverage <- coverage[
    ,
    c(
      "locus_id",
      sample_ids
    ),
    with = FALSE
  ]

  M <- as.matrix(
    methylated[
      ,
      -1,
      with = FALSE
    ]
  )

  Cov <- as.matrix(
    coverage[
      ,
      -1,
      with = FALSE
    ]
  )

  storage.mode(M) <- "integer"
  storage.mode(Cov) <- "integer"

  # DSS/BSseq represents absent RRBS observations as zero coverage.
  M[is.na(M)] <- 0L
  Cov[is.na(Cov)] <- 0L

  stopifnot(
    all(M <= Cov),
    all(M >= 0),
    all(Cov >= 0)
  )

  loci <- loci[
    match(
      methylated$locus_id,
      locus_id
    )
  ]

  bs <- BSseq(
    chr = loci$chr,
    pos = loci$start,
    M = M,
    Cov = Cov,
    sampleNames = sample_ids
  )

  group_vector <- metadata[
    match(
      sample_ids,
      id
    ),
    get(group_col)
  ]

  list(
    bs = bs,
    group = group_vector,
    metadata = metadata
  )
}

# -------------------------------------------------------------------------
# Helper: run a two-group DSS comparison
# -------------------------------------------------------------------------

run_dss_comparison <- function(
  rrbs_data,
  metadata,
  group_col,
  group_1,
  group_2,
  eligible_cpgs,
  output_file
) {

  cat(
    "\n",
    group_1,
    " vs ",
    group_2,
    "\n",
    sep = ""
  )

  bs_data <- make_bsseq(
    rrbs_data = rrbs_data,
    metadata = metadata,
    group_col = group_col,
    group_1 = group_1,
    group_2 = group_2,
    eligible_cpgs = eligible_cpgs
  )

  group_vector <- bs_data$group

  cat(
    "Samples in ",
    group_1,
    ": ",
    sum(group_vector == group_1),
    "\n",
    sep = ""
  )

  cat(
    "Samples in ",
    group_2,
    ": ",
    sum(group_vector == group_2),
    "\n",
    sep = ""
  )

  # DSS requires group labels 1 and 2
  dss_group <- ifelse(
    group_vector == group_1,
    1,
    2
  )

fit <- DMLtest(
  bs_data$bs,
  group1 = which(dss_group == 1),
  group2 = which(dss_group == 2),
  smoothing = FALSE
)

result <- as.data.table(fit)

cat("\nDSS result columns:\n")
print(names(result))

stopifnot(
  all(c(
    "chr",
    "pos",
    "mu1",
    "mu2",
    "diff",
    "pval",
    "fdr"
  ) %in% names(result))
)

result[
  ,
  comparison := paste(
    group_1,
    "vs",
    group_2
  )
]

result[
  ,
  abs_diff := abs(diff)
]

setorder(
  result,
  fdr,
  -abs_diff
)

write_tsv(
  as_tibble(result),
  output_file
)

cat(
  "CpGs tested: ",
  format(nrow(result), big.mark = ","),
  "\n",
  "FDR < 0.05: ",
  sum(result$fdr < 0.05, na.rm = TRUE),
  "\n",
  "FDR < 0.05 and absolute methylation difference >= 0.10: ",
  sum(
    result$fdr < 0.05 &
      result$abs_diff >= 0.10,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

result

# =========================================================================
# Comparison 1
# Radiation necrosis diagnosis vs recurrent tumor
#
# Exploratory because the radiation-necrosis group contains only four
# samples.
# =========================================================================

cat("\nComparison 1 prevalence filtering\n")

prevalence_1 <- get_eligible_cpgs(
  rrbs_data = rrbs,
  metadata = comparison_1,
  group_col = "comparison_1",
  threshold = prevalence_threshold
)

eligible_cpg_1 <- prevalence_1$eligible

cat(
  "CpGs passing >=50% prevalence in both groups: ",
  format(nrow(eligible_cpg_1), big.mark = ","),
  "\n",
  sep = ""
)

stopifnot(
  nrow(eligible_cpg_1) > 0
)

dm_1 <- run_dss_comparison(
  rrbs_data = rrbs,
  metadata = comparison_1,
  group_col = "comparison_1",
  group_1 = "Radiation necrosis",
  group_2 = "Recurrent tumor",
  eligible_cpgs = eligible_cpg_1,
  output_file =
    "results/differential_methylation/radnec_diagnosis_vs_recurrent.tsv"
)

# =========================================================================
# Comparison 2
# Necrosis phenotype
# =========================================================================

cat("\nComparison 2 prevalence filtering\n")

prevalence_2 <- get_eligible_cpgs(
  rrbs_data = rrbs,
  metadata = comparison_2,
  group_col = "comparison_2",
  threshold = prevalence_threshold
)

eligible_cpg_2 <- prevalence_2$eligible

cat(
  "CpGs passing >=50% prevalence in all three groups: ",
  format(nrow(eligible_cpg_2), big.mark = ","),
  "\n",
  sep = ""
)

# This should reproduce the QC result (~42,516 CpGs)
stopifnot(
  nrow(eligible_cpg_2) > 0
)

# -------------------------------------------------------------------------
# Radiation necrosis vs neither
# -------------------------------------------------------------------------

dm_radnec_vs_neither <- run_dss_comparison(
  rrbs_data = rrbs,
  metadata = comparison_2,
  group_col = "comparison_2",
  group_1 = "Radiation necrosis",
  group_2 = "Neither",
  eligible_cpgs = eligible_cpg_2,
  output_file =
    "results/differential_methylation/radnec_vs_neither.tsv"
)

# -------------------------------------------------------------------------
# Radiation necrosis vs pseudopalisading necrosis
# -------------------------------------------------------------------------

dm_radnec_vs_pseudo <- run_dss_comparison(
  rrbs_data = rrbs,
  metadata = comparison_2,
  group_col = "comparison_2",
  group_1 = "Radiation necrosis",
  group_2 = "Pseudopalisading necrosis",
  eligible_cpgs = eligible_cpg_2,
  output_file =
    "results/differential_methylation/radnec_vs_pseudopalisading.tsv"
)

# -------------------------------------------------------------------------
# Pseudopalisading necrosis vs neither
# -------------------------------------------------------------------------

dm_pseudo_vs_neither <- run_dss_comparison(
  rrbs_data = rrbs,
  metadata = comparison_2,
  group_col = "comparison_2",
  group_1 = "Pseudopalisading necrosis",
  group_2 = "Neither",
  eligible_cpgs = eligible_cpg_2,
  output_file =
    "results/differential_methylation/pseudopalisading_vs_neither.tsv"
)

# -------------------------------------------------------------------------
# Save eligible CpGs
# -------------------------------------------------------------------------

write_tsv(
  as_tibble(eligible_cpg_1),
  "results/differential_methylation/comparison_1_eligible_cpgs.tsv"
)

write_tsv(
  as_tibble(eligible_cpg_2),
  "results/differential_methylation/comparison_2_eligible_cpgs.tsv"
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

summary_results <- tibble(
  comparison = c(
    "RN diagnosis vs recurrent tumor",
    "RN vs neither",
    "RN vs pseudopalisading",
    "Pseudopalisading vs neither"
  ),

  n_cpg_tested = c(
    nrow(dm_1),
    nrow(dm_radnec_vs_neither),
    nrow(dm_radnec_vs_pseudo),
    nrow(dm_pseudo_vs_neither)
  ),

  fdr_05 = c(
    sum(
      dm_1$fdr < 0.05,
      na.rm = TRUE
    ),
    sum(
      dm_radnec_vs_neither$fdr < 0.05,
      na.rm = TRUE
    ),
    sum(
      dm_radnec_vs_pseudo$fdr < 0.05,
      na.rm = TRUE
    ),
    sum(
      dm_pseudo_vs_neither$fdr < 0.05,
      na.rm = TRUE
    )
  ),

  fdr_05_delta_10 = c(
    sum(
      dm_1$fdr < 0.05 &
        dm_1$abs_diff >= 0.10,
      na.rm = TRUE
    ),
    sum(
      dm_radnec_vs_neither$fdr < 0.05 &
        dm_radnec_vs_neither$abs_diff >= 0.10,
      na.rm = TRUE
    ),
    sum(
      dm_radnec_vs_pseudo$fdr < 0.05 &
        dm_radnec_vs_pseudo$abs_diff >= 0.10,
      na.rm = TRUE
    ),
    sum(
      dm_pseudo_vs_neither$fdr < 0.05 &
        dm_pseudo_vs_neither$abs_diff >= 0.10,
      na.rm = TRUE
    )
  )
)

cat("\nDifferential methylation summary:\n")
print(summary_results)

write_tsv(
  summary_results,
  "results/differential_methylation/summary.tsv"
)

cat("\nDifferential methylation summary:\n")
print(summary_results)

write_tsv(
  summary_results,
  "results/differential_methylation/summary.tsv"
)

cat("\nDifferential methylation analysis complete.\n")
