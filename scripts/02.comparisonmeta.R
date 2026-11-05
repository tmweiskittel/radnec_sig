library(dplyr)
library(readr)

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------

dir.create(
  "results/comparisons",
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------
# Read metadata and define first-recurrence cohort
# -------------------------------------------------------------------------

meta <- read_tsv(
  "meta/GBMatch_sampleAnnotation (1).tsv",
  show_col_types = FALSE
)

included_diagnoses <- c(
  "GBM classic",
  "GBM giant cell",
  "Gliosarcoma",
  "Radiation necrosis"
)

# Restrict all subsequent analyses to first-recurrence GBM-related specimens
first_recurrence <- meta |>
  filter(
    surgery == 2,
    category == "GBMatch",
    WHO2016_classification %in% included_diagnoses
  ) |>
  mutate(
    radnec = case_when(
      `Radiation necrosis` == "Present" ~ TRUE,
      `Radiation necrosis` == "Absent" ~ FALSE,
      TRUE ~ NA
    ),
    pseudopalisading = case_when(
      `Pseudopalisading necrosis` %in% c(
        "Present",
        "Abundant"
      ) ~ TRUE,
      `Pseudopalisading necrosis` == "Absent" ~ FALSE,
      TRUE ~ NA
    )
  )

# Confirm cohort restriction
stopifnot(
  all(first_recurrence$surgery == 2)
)

cat(
  "\nFirst-recurrence GBM-related cohort:\n",
  "Samples: ", nrow(first_recurrence), "\n",
  "Patients: ", n_distinct(first_recurrence$patID), "\n",
  sep = ""
)

first_recurrence |>
  count(
    WHO2016_classification,
    name = "n_samples"
  ) |>
  print()

# -------------------------------------------------------------------------
# Comparison 1:
# Radiation necrosis diagnosis vs first-recurrence tumor
# without annotated radiation necrosis
# -------------------------------------------------------------------------

comparison_1 <- first_recurrence |>
  mutate(
    comparison_1 = case_when(

      # Pathologist diagnosis of radiation necrosis
      WHO2016_classification == "Radiation necrosis" ~
        "Radiation necrosis",

      # First-recurrence tumor with explicitly absent RN annotation
      WHO2016_classification != "Radiation necrosis" &
        radnec == FALSE ~
        "Recurrent tumor",

      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(comparison_1)
  ) |>
  mutate(
    comparison_1 = factor(
      comparison_1,
      levels = c(
        "Recurrent tumor",
        "Radiation necrosis"
      )
    )
  )

# Confirm all comparison samples are from surgery 2
stopifnot(
  all(comparison_1$surgery == 2)
)

# Sample counts
comparison_1_counts <- comparison_1 |>
  count(
    comparison_1,
    name = "n_samples"
  )

cat("\nComparison 1:\n")
print(comparison_1_counts)

# Patient counts
comparison_1 |>
  group_by(comparison_1) |>
  summarise(
    n_samples = n(),
    n_patients = n_distinct(patID),
    .groups = "drop"
  ) |>
  print()

# Check pathology annotations among WHO-diagnosed RN samples
cat("\nRadiation necrosis diagnosis pathology annotations:\n")

comparison_1 |>
  filter(
    comparison_1 == "Radiation necrosis"
  ) |>
  count(
    `Radiation necrosis`,
    `Pseudopalisading necrosis`
  ) |>
  print()

# -------------------------------------------------------------------------
# Comparison 2:
# Radiation necrosis vs pseudopalisading necrosis vs neither
# at first recurrence
# -------------------------------------------------------------------------

comparison_2 <- first_recurrence |>
  mutate(
    comparison_2 = case_when(

      # Radiation necrosis without pseudopalisading necrosis
      radnec == TRUE &
        pseudopalisading == FALSE ~
        "Radiation necrosis",

      # Pseudopalisading necrosis without radiation necrosis
      radnec == FALSE &
        pseudopalisading == TRUE ~
        "Pseudopalisading necrosis",

      # Neither form of necrosis annotated
      radnec == FALSE &
        pseudopalisading == FALSE ~
        "Neither",

      # Both-positive and missing annotations are excluded
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(comparison_2)
  ) |>
  mutate(
    comparison_2 = factor(
      comparison_2,
      levels = c(
        "Neither",
        "Pseudopalisading necrosis",
        "Radiation necrosis"
      )
    )
  )

# Confirm all comparison samples are from surgery 2
stopifnot(
  all(comparison_2$surgery == 2)
)

# Sample counts
comparison_2_counts <- comparison_2 |>
  count(
    comparison_2,
    name = "n_samples"
  )

cat("\nComparison 2:\n")
print(comparison_2_counts)

# Patient counts
comparison_2 |>
  group_by(comparison_2) |>
  summarise(
    n_samples = n(),
    n_patients = n_distinct(patID),
    .groups = "drop"
  ) |>
  print()

# -------------------------------------------------------------------------
# Save comparison definitions
# -------------------------------------------------------------------------

comparison_counts <- bind_rows(
  comparison_1_counts |>
    transmute(
      comparison = "Radiation necrosis diagnosis vs first-recurrence tumor",
      group = as.character(comparison_1),
      n_samples
    ),

  comparison_2_counts |>
    transmute(
      comparison = "Necrosis type at first recurrence",
      group = as.character(comparison_2),
      n_samples
    )
)

write_tsv(
  comparison_counts,
  "results/comparisons/sample_counts.tsv"
)

write_tsv(
  comparison_1 |>
    select(
      patID,
      surgery,
      id,
      WHO2016_classification,
      `Radiation necrosis`,
      `Pseudopalisading necrosis`,
      comparison_1
    ),
  "results/comparisons/radnec_vs_recurrent.tsv"
)

write_tsv(
  comparison_2 |>
    select(
      patID,
      surgery,
      id,
      WHO2016_classification,
      `Radiation necrosis`,
      `Pseudopalisading necrosis`,
      comparison_2
    ),
  "results/comparisons/necrosis_type.tsv"
)

# -------------------------------------------------------------------------
# Necrotic tissue content in necrosis-positive samples
# -------------------------------------------------------------------------

necrosis_tissue_counts <- first_recurrence |>
  filter(
    radnec == TRUE |
      pseudopalisading == TRUE
  ) |>
  mutate(
    necrosis_type = case_when(
      radnec == TRUE &
        pseudopalisading == FALSE ~
        "Radiation necrosis",

      radnec == FALSE &
        pseudopalisading == TRUE ~
        "Pseudopalisading necrosis",

      radnec == TRUE &
        pseudopalisading == TRUE ~
        "Both",

      TRUE ~ NA_character_
    )
  ) |>
  group_by(necrosis_type) |>
  summarise(
    n_total = n(),
    n_necrosis_measure = sum(
      !is.na(`Relative share necrosis`)
    ),
    n_ge_25 = sum(
      `Relative share necrosis` >= 0.25,
      na.rm = TRUE
    ),
    n_ge_50 = sum(
      `Relative share necrosis` >= 0.50,
      na.rm = TRUE
    ),
    n_ge_75 = sum(
      `Relative share necrosis` >= 0.75,
      na.rm = TRUE
    ),
    median_necrosis = median(
      `Relative share necrosis`,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(necrosis_tissue_counts)

write_tsv(
  necrosis_tissue_counts,
  "results/comparisons/necrosis_tissue_threshold_counts.tsv"
)

cat("\nComparison metadata complete.\n")
