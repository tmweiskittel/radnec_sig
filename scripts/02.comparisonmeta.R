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
# Read metadata
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

# Restrict immediately to recurrent GBM-related specimens
recurrent <- meta |>
  filter(
    surgery > 1,
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

# -------------------------------------------------------------------------
# Comparison 1:
# Radiation necrosis diagnosis vs recurrent tumor without radiation necrosis
# -------------------------------------------------------------------------

comparison_1 <- recurrent |>
  mutate(
    comparison_1 = case_when(

      # Pathologist diagnosis of radiation necrosis
      WHO2016_classification == "Radiation necrosis" ~
        "Radiation necrosis",

      # Recurrent GBM-related tumor with explicitly absent RN annotation
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

# QA
comparison_1_counts <- comparison_1 |>
  count(
    comparison_1,
    name = "n_samples"
  )

print(comparison_1_counts)

# Check number of unique patients
comparison_1 |>
  summarise(
    n_samples = n(),
    n_patients = n_distinct(patID),
    .by = comparison_1
  ) |>
  print()

# Check the RN pathology flag among WHO-diagnosed RN samples
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
# -------------------------------------------------------------------------

comparison_2 <- recurrent |>
  mutate(
    comparison_2 = case_when(

      radnec == TRUE &
        pseudopalisading == FALSE ~
        "Radiation necrosis",

      radnec == FALSE &
        pseudopalisading == TRUE ~
        "Pseudopalisading necrosis",

      radnec == FALSE &
        pseudopalisading == FALSE ~
        "Neither",

      # Both-positive and missing samples are deliberately excluded
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

# QA
comparison_2_counts <- comparison_2 |>
  count(
    comparison_2,
    name = "n_samples"
  )

print(comparison_2_counts)

comparison_2 |>
  summarise(
    n_samples = n(),
    n_patients = n_distinct(patID),
    .by = comparison_2
  ) |>
  print()

# -------------------------------------------------------------------------
# Save comparison definitions
# -------------------------------------------------------------------------

comparison_counts <- bind_rows(
  comparison_1_counts |>
    transmute(
      comparison = "Radiation necrosis diagnosis vs recurrent tumor",
      group = as.character(comparison_1),
      n_samples
    ),

  comparison_2_counts |>
    transmute(
      comparison = "Necrosis type",
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
