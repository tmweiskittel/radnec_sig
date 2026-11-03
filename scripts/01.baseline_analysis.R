library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(ggplot2)

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------

dir.create(
  "results/baseline",
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------
# Read metadata and define analysis cohort
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

# Restrict all subsequent analyses to recurrent GBM-related samples
recurrent <- meta |>
  filter(
    surgery > 1,
    WHO2016_classification %in% included_diagnoses
  )

# -------------------------------------------------------------------------
# Define necrosis classes
# -------------------------------------------------------------------------

recurrent <- recurrent |>
  mutate(
    radnec = case_when(
      `Radiation necrosis` == "Present" ~ TRUE,
      `Radiation necrosis` == "Absent" ~ FALSE,
      TRUE ~ NA
    ),
    pseudopalisading = case_when(
      `Pseudopalisading necrosis` %in% c("Present", "Abundant") ~ TRUE,
      `Pseudopalisading necrosis` == "Absent" ~ FALSE,
      TRUE ~ NA
    ),
    necrosis_class = case_when(
      radnec & !pseudopalisading ~ "Radiation necrosis only",
      !radnec & pseudopalisading ~ "Pseudopalisading only",
      radnec & pseudopalisading ~ "Both",
      !radnec & !pseudopalisading ~ "Neither",
      TRUE ~ "Missing annotation"
    ),
    necrosis_class = factor(
      necrosis_class,
      levels = c(
        "Neither",
        "Pseudopalisading only",
        "Radiation necrosis only",
        "Both",
        "Missing annotation"
      )
    )
  )

# -------------------------------------------------------------------------
# Sample counts
# -------------------------------------------------------------------------

sample_counts <- bind_rows(
  tibble(
    category = "Total",
    group = "All recurrent GBM-related samples",
    n = nrow(recurrent)
  ),
  recurrent |>
    count(WHO2016_classification, name = "n") |>
    transmute(
      category = "Diagnosis",
      group = WHO2016_classification,
      n
    ),
  recurrent |>
    count(necrosis_class, name = "n") |>
    transmute(
      category = "Necrosis class",
      group = as.character(necrosis_class),
      n
    )
)

print(sample_counts, n = Inf)

write_tsv(
  sample_counts,
  "results/baseline/sample_counts.tsv"
)

# -------------------------------------------------------------------------
# Check radiation necrosis diagnosis against pathology annotation
# -------------------------------------------------------------------------

radnec_diagnosis_check <- recurrent |>
  filter(
    WHO2016_classification == "Radiation necrosis"
  ) |>
  select(
    patID,
    surgery,
    id,
    WHO2016_classification,
    `Radiation necrosis`,
    `Pseudopalisading necrosis`,
    necrosis_class
  )

print(radnec_diagnosis_check, n = Inf)

# -------------------------------------------------------------------------
# Continuous necrosis variables
# -------------------------------------------------------------------------

necrosis_vars <- c(
  "Relative share necrosis",
  "Area necrosis [mm²]",
  "Mean Nuclei # Necrosis",
  "Median Nuclei # Necrosis",
  "Necrotic  (mm3)",
  "Proportion Necrotic (%)"
)

# -------------------------------------------------------------------------
# Correlation of continuous necrosis variables
# -------------------------------------------------------------------------

necrosis_cor <- recurrent |>
  select(all_of(necrosis_vars)) |>
  cor(
    use = "pairwise.complete.obs",
    method = "spearman"
  )

print(round(necrosis_cor, 2))

write_csv(
  as.data.frame(necrosis_cor) |>
    rownames_to_column("variable"),
  "results/baseline/necrosis_correlations.csv"
)

cor_long <- as.data.frame(necrosis_cor) |>
  rownames_to_column("variable_1") |>
  pivot_longer(
    -variable_1,
    names_to = "variable_2",
    values_to = "rho"
  )

p_cor <- ggplot(
  cor_long,
  aes(
    x = variable_1,
    y = variable_2,
    fill = rho
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", rho)),
    size = 3
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0
  ) +
  coord_equal() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Spearman rho"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  "results/baseline/necrosis_correlation_heatmap.png",
  p_cor,
  width = 8,
  height = 7,
  dpi = 300
)

# -------------------------------------------------------------------------
# Continuous necrosis variables by necrosis class
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Distribution of continuous necrosis variables by necrosis class
# -------------------------------------------------------------------------

necrosis_long <- recurrent |>
  filter(
    necrosis_class != "Missing annotation"
  ) |>
  select(
    patID,
    surgery,
    id,
    necrosis_class,
    all_of(necrosis_vars)
  ) |>
  pivot_longer(
    cols = all_of(necrosis_vars),
    names_to = "variable",
    values_to = "value"
  )

# Create one overlaid histogram for each continuous necrosis variable
necrosis_long |>
  filter(!is.na(value)) |>
  group_split(variable) |>
  walk(function(dat) {

    variable_name <- unique(dat$variable)

p <- ggplot(
  dat,
  aes(
    x = value,
    fill = necrosis_class,
    color = necrosis_class
  )
) +
  geom_histogram(
    bins = 20,
    position = "identity",
    alpha = 0.4,
    linewidth = 0.5
  ) +
  labs(
    x = variable_name,
    y = "Number of samples",
    fill = "Necrosis class",
    color = "Necrosis class"
  ) +
  theme_classic()


radnec_variables <- recurrent |>
  filter(
    WHO2016_classification == "Radiation necrosis"
  ) |>
  select(
    patID,
    surgery,
    id,
    `Radiation necrosis`,
    `Pseudopalisading necrosis`,
    all_of(necrosis_vars)
  )

print(
  radnec_variables,
  n = Inf,
  width = Inf
)

