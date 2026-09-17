# Summarise per-kinase LFC results from PamGene mapped-kinase files,
# then assess consistency across experimental conditions.
# Outputs:
#   results/hyperprocessed/kinase_mapped_dpp_results.csv     — long-form per-kinase median LFC
#   results/hyperprocessed/kinase_mapped_dpp_matrix.csv      — wide LFC matrix (Kinase × condition)
#   results/hyperprocessed/kinase_mapped_dpp_correlation.csv — Spearman correlation between conditions

suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(purrr)
    library(stringr)
    library(tidyr)
    library(dplyr)
})

# ── Load data ────────────────────────────────────────────────────────────────

mapped_files <- list.files("results", "mapped_kinases", full.names = TRUE)

mapped_data <- mapped_files |>
    set_names(\(x) basename(x)) |>
    map(read_csv, show_col_types = FALSE) |>
    bind_rows(.id = "source_file") |>
    group_by(source_file, Kinase) |>
    summarise(lfc = median(LFC), .groups = "drop") |>
    mutate(
        condition = str_extract(
            source_file,
            "hpv_(neg|pos)_\\d{2}\\.\\d{2}_(\\w+)-STK",
            group = 2
        ),
        biologically_meaningful = abs(lfc) > 0.2
    ) |>
    select(-source_file)

write_csv(mapped_data, "results/hyperprocessed/kinase_mapped_dpp_results.csv")

# ── Build wide comparison matrix ─────────────────────────────────────────────

comparison_matrix <- mapped_data |>
    mutate(lfc = if_else(biologically_meaningful, lfc, 0)) |>
    select(-biologically_meaningful) |>
    pivot_wider(names_from = condition, values_from = lfc, values_fill = 0) |>
    select(
        Kinase,
        starts_with("Neg5m"),  starts_with("Neg15m"), starts_with("Neg30m"),
        starts_with("Neg1h"),  starts_with("Neg4h"),  starts_with("Neg24h"),
        starts_with("Pos5m"),  starts_with("Pos15m"), starts_with("Pos30m"),
        starts_with("Pos1h"),  starts_with("Pos4h"),  starts_with("Pos24h")
    ) |>
    arrange(Kinase)

write_csv(comparison_matrix, "results/hyperprocessed/kinase_mapped_dpp_matrix.csv")

# ── Spearman correlation between conditions ───────────────────────────────────

cor_matrix <- comparison_matrix |>
    column_to_rownames("Kinase") |>
    as.matrix()

cor(cor_matrix, method = "spearman") |>
    as.data.frame() |>
    rownames_to_column("condition") |>
    write_csv("results/hyperprocessed/kinase_mapped_dpp_correlation.csv")
