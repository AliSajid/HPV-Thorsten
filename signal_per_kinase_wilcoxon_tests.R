# For each signal file, join peptide slopes to the kinase mapping, then
# perform a paired Wilcoxon signed-rank test per kinase using the slopes
# of all peptides mapped to that kinase as the observations.
#
# Pairing is by barcode (within-array replicates). Each peptide × barcode
# pair contributes one paired observation, so kinases with more mapped
# peptides have more power.
#
# Outputs:
#   results/hyperprocessed/signal_kinase_results.csv — per-kinase test results
#                                                       across all comparisons

suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(purrr)
    library(stringr)
    library(tidyr)
    library(dplyr)
})

# ── Load kinase → peptide mapping ─────────────────────────────────────────────

# The raw mapped_kinases files carry one row per peptide with its Kinase label
# and the LFC for that file's comparison. We only need the Kinase ↔ Peptide
# mapping here, so we read all files and keep distinct pairs.
kinase_peptide_map <- list.files("results", "mapped_kinases", full.names = TRUE) |>
    map(read_csv, show_col_types = FALSE) |>
    bind_rows() |>
    select(Kinase, Peptide) |>
    distinct()

# ── Helpers ───────────────────────────────────────────────────────────────────

# Run a paired Wilcoxon signed-rank test for one kinase within one file.
# Each peptide mapped to the kinase contributes one paired observation per
# barcode, so the test uses n_peptides × n_barcodes pairs in total.
#
# Effect size: matched-pairs rank-biserial correlation,
#   r = W / (n_pairs * (n_pairs + 1) / 2)
# where n_pairs = n_peptides × n_barcodes. Ranges from -1 to 1; sign
# reflects direction of change (positive = treatment > control).
wilcoxon_kinase <- function(df) {
    paired <- df |>
        select(Barcode, Peptide, Group, slope) |>
        pivot_wider(names_from = Group, values_from = slope)

    treatment  <- paired[[3]]   # first non-Barcode/Peptide column = treatment
    control    <- paired[[4]]   # second                           = control

    n_peptides <- n_distinct(df$Peptide)
    n_pairs    <- nrow(paired)
    test       <- wilcox.test(treatment, control, paired = TRUE, exact = FALSE)

    tibble(
        n_peptides  = n_peptides,
        n_pairs     = n_pairs,
        mean_diff   = mean(treatment - control),
        statistic   = test$statistic,
        effect_size = test$statistic / (n_pairs * (n_pairs + 1) / 2),
        p_value     = test$p.value
    )
}

# ── Load signal files and join kinase mapping ──────────────────────────────────

signal_files <- list.files("results", "signal", full.names = TRUE)

signal_data <- signal_files |>
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
    inner_join(kinase_peptide_map, by = "Peptide")

# ── Per-kinase paired Wilcoxon test ───────────────────────────────────────────

results <- signal_data |>
    nest(data = -c(source_file, condition, Kinase)) |>
    mutate(test = map(data, wilcoxon_kinase)) |>
    select(-data) |>
    unnest(test)

# ── Multiple testing correction and output ────────────────────────────────────

signal_kinase_results <- results |>
    group_by(condition) |>
    mutate(p_adjusted = p.adjust(p_value, method = "BH")) |>
    ungroup() |>
    select(condition, Kinase, n_peptides, n_pairs, mean_diff, statistic, effect_size, p_value, p_adjusted, source_file) |>
    arrange(condition, p_value)

write_csv(
    signal_kinase_results,
    "results/hyperprocessed/signal_kinase_results.csv"
)
