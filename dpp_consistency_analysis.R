# Assess direction and consistency of change for differentially
# phosphorylated peptides (DPP) across experimental conditions.
# Outputs:
#   results/hyperprocessed/dpp_results.csv     — long-form LFC per peptide and condition
#   results/hyperprocessed/dpp_matrix.csv      — wide LFC matrix (Peptide × condition)
#   results/hyperprocessed/dpp_correlation.csv — Spearman correlation between conditions

suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(purrr)
    library(stringr)
    library(tidyr)
    library(dplyr)
})

# ── Load data ────────────────────────────────────────────────────────────────

dpp_files <- list.files("results", "dpp", full.names = TRUE)

combined_data <- dpp_files |>
    set_names(\(x) basename(x)) |>
    map(read_csv, show_col_types = FALSE) |>
    bind_rows(.id = "source_file") |>
    select(Peptide, lfc = totalMeanLFC, source_file) |>
    distinct() |>
    mutate(
        condition = str_extract(source_file, "dpp_(\\w+)-STK\\.csv", group = 1)
    )

# ── Long-form results ─────────────────────────────────────────────────────────

combined_data |>
    select(Peptide, condition, lfc) |>
    arrange(Peptide, condition) |>
    write_csv("results/hyperprocessed/dpp_results.csv")

# ── Build wide comparison matrix ─────────────────────────────────────────────

comparison_matrix <- combined_data |>
    select(-source_file) |>
    pivot_wider(names_from = condition, values_from = lfc, values_fill = 0) |>
    select(
        Peptide,
        starts_with("Neg5m"),  starts_with("Neg15m"), starts_with("Neg30m"),
        starts_with("Neg1h"),  starts_with("Neg4h"),  starts_with("Neg24h"),
        starts_with("Pos5m"),  starts_with("Pos15m"), starts_with("Pos30m"),
        starts_with("Pos1h"),  starts_with("Pos4h"),  starts_with("Pos24h")
    ) |>
    arrange(Peptide)

write_csv(comparison_matrix, "results/hyperprocessed/dpp_matrix.csv")

# ── Spearman correlation between conditions ───────────────────────────────────

cor_matrix <- comparison_matrix |>
    column_to_rownames("Peptide") |>
    as.matrix()

cor(cor_matrix, method = "spearman") |>
    as.data.frame() |>
    rownames_to_column("condition") |>
    write_csv("results/hyperprocessed/dpp_correlation.csv")
