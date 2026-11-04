library(data.table)
library(dplyr)
library(readr)
library(tibble)

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------

dir.create(
  "results/methylation_qc",
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------
# Load comparison definitions
# -------------------------------------------------------------------------

comparison_1 <- read_tsv(
  "results/comparisons/radnec_vs_recurrent.tsv",
  show_col_types = FALSE
)

comparison_2 <- read_tsv(
  "results/comparisons/necrosis_type.tsv",
  show_col_types = FALSE
)

# Sample IDs included in either comparison
comparison_ids <- unique(c(
  comparison_1$id,
  comparison_2$id
))

# -------------------------------------------------------------------------
# Load single-CpG RRBS data
# -------------------------------------------------------------------------

load("data/klughammer/rrbsCg.RData")

stopifnot(
  exists("ret"),
  is.data.table(ret)
)

cat(
  "\nFull RRBS object\n",
  "Rows: ", format(nrow(ret), big.mark = ","), "\n",
  "Samples: ", uniqueN(ret$id), "\n",
  sep = ""
)

# -------------------------------------------------------------------------
# Check sample matching
# -------------------------------------------------------------------------

sample_match <- data.table(
  id = comparison_ids,
  in_rrbs = comparison_ids %in% ret$id
)

cat("\nSample matching:\n")

print(
  sample_match[
    ,
    .N,
    by = in_rrbs
  ]
)

if (any(!sample_match$in_rrbs)) {
  cat("\nSamples missing from RRBS data:\n")
  print(
    sample_match[in_rrbs == FALSE]
  )
}

stopifnot(
  all(sample_match$in_rrbs)
)

# -------------------------------------------------------------------------
# Restrict RRBS data to samples used in either comparison
# -------------------------------------------------------------------------

rrbs <- ret[
  id %in% comparison_ids
]

rm(ret)
gc()

cat(
  "\nFiltered RRBS object\n",
  "Rows: ", format(nrow(rrbs), big.mark = ","), "\n",
  "Samples: ", uniqueN(rrbs$id), "\n",
  sep = ""
)

# -------------------------------------------------------------------------
# Confirm methylation calculation
# -------------------------------------------------------------------------

max_methyl_difference <- rrbs[
  ,
  max(
    abs(methyl - hitCount / readCount),
    na.rm = TRUE
  )
]

cat(
  "\nMaximum difference between supplied methylation and hitCount/readCount: ",
  max_methyl_difference,
  "\n",
  sep = ""
)

# -------------------------------------------------------------------------
# Sample-level CpG coverage
# -------------------------------------------------------------------------

sample_coverage <- rrbs[
  ,
  .(
    n_cpg = .N,
    median_read_count = median(
      readCount,
      na.rm = TRUE
    ),
    mean_read_count = mean(
      readCount,
      na.rm = TRUE
    ),
    min_read_count = min(
      readCount,
      na.rm = TRUE
    ),
    max_read_count = max(
      readCount,
      na.rm = TRUE
    )
  ),
  by = id
]

# -------------------------------------------------------------------------
# Sample coverage for comparison 1
# -------------------------------------------------------------------------

comp1_meta <- comparison_1 |>
  select(
    id,
    comparison_1
  ) |>
  distinct()

sample_coverage_1 <- sample_coverage |>
  as_tibble() |>
  inner_join(
    comp1_meta,
    by = "id"
  )

coverage_summary_1 <- sample_coverage_1 |>
  group_by(comparison_1) |>
  summarise(
    n_samples = n(),
    min_cpg = min(n_cpg),
    q1_cpg = quantile(n_cpg, 0.25),
    median_cpg = median(n_cpg),
    q3_cpg = quantile(n_cpg, 0.75),
    max_cpg = max(n_cpg),
    .groups = "drop"
  )

cat("\nComparison 1 sample coverage:\n")
print(coverage_summary_1)

cat("\nLowest-coverage samples in comparison 1:\n")

sample_coverage_1 |>
  arrange(n_cpg) |>
  print(n = 20)

cat("\nRadiation necrosis diagnosis samples:\n")

sample_coverage_1 |>
  filter(
    comparison_1 == "Radiation necrosis"
  ) |>
  arrange(n_cpg) |>
  print()

# -------------------------------------------------------------------------
# Sample coverage for comparison 2
# -------------------------------------------------------------------------

comp2_meta <- comparison_2 |>
  select(
    id,
    comparison_2
  ) |>
  distinct()

sample_coverage_2 <- sample_coverage |>
  as_tibble() |>
  inner_join(
    comp2_meta,
    by = "id"
  )

coverage_summary_2 <- sample_coverage_2 |>
  group_by(comparison_2) |>
  summarise(
    n_samples = n(),
    min_cpg = min(n_cpg),
    q1_cpg = quantile(n_cpg, 0.25),
    median_cpg = median(n_cpg),
    q3_cpg = quantile(n_cpg, 0.75),
    max_cpg = max(n_cpg),
    .groups = "drop"
  )

cat("\nComparison 2 sample coverage:\n")
print(coverage_summary_2)

cat("\nLowest-coverage samples in comparison 2:\n")

sample_coverage_2 |>
  arrange(n_cpg) |>
  print(n = 20)

# -------------------------------------------------------------------------
# Overall CpG overlap
# -------------------------------------------------------------------------

cpg_coverage <- rrbs[
  ,
  .(
    n_samples = uniqueN(id)
  ),
  by = .(
    chr,
    start
  )
]

cat("\nOverall CpG coverage across comparison samples:\n")
print(
  summary(cpg_coverage$n_samples)
)

cat(
  "\nUnique CpGs: ",
  format(nrow(cpg_coverage), big.mark = ","),
  "\n",
  sep = ""
)

# -------------------------------------------------------------------------
# CpG prevalence within comparison 1
# -------------------------------------------------------------------------

rrbs_comp1 <- rrbs[
  id %in% comp1_meta$id
]

rrbs_comp1[
  as.data.table(comp1_meta),
  comparison_1 := i.comparison_1,
  on = "id"
]

cpg_coverage_1 <- rrbs_comp1[
  ,
  .(
    n_samples = uniqueN(id)
  ),
  by = .(
    chr,
    start,
    comparison_1
  )
]

group_sizes_1 <- as.data.table(
  comp1_meta |>
    count(
      comparison_1,
      name = "group_n"
    )
)

cpg_coverage_1 <- merge(
  cpg_coverage_1,
  group_sizes_1,
  by = "comparison_1"
)

cpg_coverage_1[
  ,
  proportion := n_samples / group_n
]

prevalence_summary_1 <- cpg_coverage_1[
  ,
  .(
    cpg_25pct = sum(proportion >= 0.25),
    cpg_50pct = sum(proportion >= 0.50),
    cpg_75pct = sum(proportion >= 0.75),
    cpg_90pct = sum(proportion >= 0.90)
  ),
  by = comparison_1
]

cat("\nComparison 1 CpG prevalence:\n")
print(prevalence_summary_1)

# -------------------------------------------------------------------------
# CpG prevalence within comparison 2
# -------------------------------------------------------------------------

rrbs_comp2 <- rrbs[
  id %in% comp2_meta$id
]

rrbs_comp2[
  as.data.table(comp2_meta),
  comparison_2 := i.comparison_2,
  on = "id"
]

cpg_coverage_2 <- rrbs_comp2[
  ,
  .(
    n_samples = uniqueN(id)
  ),
  by = .(
    chr,
    start,
    comparison_2
  )
]

group_sizes_2 <- as.data.table(
  comp2_meta |>
    count(
      comparison_2,
      name = "group_n"
    )
)

cpg_coverage_2 <- merge(
  cpg_coverage_2,
  group_sizes_2,
  by = "comparison_2"
)

cpg_coverage_2[
  ,
  proportion := n_samples / group_n
]

prevalence_summary_2 <- cpg_coverage_2[
  ,
  .(
    cpg_25pct = sum(proportion >= 0.25),
    cpg_50pct = sum(proportion >= 0.50),
    cpg_75pct = sum(proportion >= 0.75),
    cpg_90pct = sum(proportion >= 0.90)
  ),
  by = comparison_2
]

cat("\nComparison 2 CpG prevalence:\n")
print(prevalence_summary_2)

# -------------------------------------------------------------------------
# Save QC summaries
# -------------------------------------------------------------------------

write_tsv(
  sample_coverage_1,
  "results/methylation_qc/comparison_1_sample_coverage.tsv"
)

write_tsv(
  sample_coverage_2,
  "results/methylation_qc/comparison_2_sample_coverage.tsv"
)

write_tsv(
  as_tibble(prevalence_summary_1),
  "results/methylation_qc/comparison_1_cpg_prevalence.tsv"
)

write_tsv(
  as_tibble(prevalence_summary_2),
  "results/methylation_qc/comparison_2_cpg_prevalence.tsv"
)

cat("\nMethylation QC complete.\n")
