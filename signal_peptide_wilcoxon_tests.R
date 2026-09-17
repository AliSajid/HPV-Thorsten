# For each signal file, perform a paired Wilcoxon signed-rank test
# (treatment vs. control) on peptide slopes across matched barcodes.
# Outputs:
#   results/hyperprocessed/signal_peptide_results.csv — per-peptide test results
#                                                        across all comparisons

suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(purrr)
    library(stringr)
    library(tidyr)
    library(dplyr)
})

# ── Helpers ───────────────────────────────────────────────────────────────────

# Run a paired Wilcoxon signed-rank test for one peptide within one file.
# Slopes are paired by barcode so that within-array technical variation is
# accounted for. Returns a one-row tibble of test statistics.
#
# Effect size: matched-pairs rank-biserial correlation,
#   r = W / (n * (n + 1) / 2)
# where n is the number of pairs. Ranges from -1 to 1; sign reflects
# direction of change (positive = treatment > control).
wilcoxon_peptide <- function(df) {
    # Pivot so each barcode is one row with treatment and control slope columns
    paired <- df |>
        select(Barcode, Group, slope) |>
        pivot_wider(names_from = Group, values_from = slope)

    treatment <- paired[[2]]   # first non-Barcode column  = treatment group
    control   <- paired[[3]]   # second non-Barcode column = control group

    n    <- nrow(paired)
    test <- wilcox.test(treatment, control, paired = TRUE, exact = FALSE)

    tibble(
        n              = n,
        mean_diff      = mean(treatment - control),
        statistic      = test$statistic,
        effect_size    = test$statistic / (n * (n + 1) / 2),
        p_value        = test$p.value
    )
}

# ── Load and process signal files ─────────────────────────────────────────────

signal_files <- list.files("results", "signal", full.names = TRUE)

results <- signal_files |>
    set_names(\(x) basename(x)) |>
    map(read_csv, show_col_types = FALSE) |>
    bind_rows(.id = "source_file") |>
    mutate(
        condition = str_extract(
            source_file,
            "signal_(\\w+)-STK\\.csv",
            group = 1
        )
    ) |>
    nest(data = -c(source_file, condition, Peptide)) |>
    mutate(test = map(data, wilcoxon_peptide)) |>
    select(-data) |>
    unnest(test)

# ── Multiple testing correction and output ────────────────────────────────────

signal_peptide_results <- results |>
    group_by(condition) |>
    mutate(p_adjusted = p.adjust(p_value, method = "BH")) |>
    ungroup() |>
    select(condition, Peptide, n, mean_diff, statistic, effect_size, p_value, p_adjusted, source_file) |>
    arrange(condition, p_value)

write_csv(
    signal_peptide_results,
    "results/hyperprocessed/signal_peptide_results.csv"
)
