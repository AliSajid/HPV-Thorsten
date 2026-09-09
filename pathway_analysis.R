# Pathway Analysis

suppressPackageStartupMessages({
  library(tidyverse)
  library(enrichR)
  library(writexl)
})

load_reference_data <- function() {
  load_specific_reference <- function(chip) {
    id_file <- file.path(
      "reference_data",
      str_glue("{str_to_lower(chip)}_id_map.Rds")
    )
    hgnc_file <- file.path(
      "reference_data",
      str_glue("{str_to_lower(chip)}_hgnc_map.Rds")
    )

    id <- read_rds(id_file)
    hgnc <- read_rds(hgnc_file)

    peptide_map <- id |>
      inner_join(hgnc)
  }

  stk <- load_specific_reference("STK")
  ptk <- load_specific_reference("PTK")

  stk |>
    bind_rows(ptk)
}

perform_pathway_enrichment <- function(filename, peptide_map,
                                       threshold = 0.10, databases = NULL) {
  dpp_data <- read_csv(file.path("results", filename)) |>
    group_by(Peptide) |>
    filter(abs(LFC) == max(abs(LFC))) |>
    ungroup() |>
    select(Peptide, LFC) |>
    inner_join(peptide_map)

  lower_bound <- quantile(dpp_data[["LFC"]], threshold)
  upper_bound <- quantile(dpp_data[["LFC"]], 1 - threshold)

  selected_genes <- dpp_data |>
    filter(LFC >= upper_bound | LFC <= lower_bound) |>
    pull(Gene)

  if (is.null(databases)) {
    databases <- c(
      "GO_Molecular_Function_2025",
      "GO_Cellular_Component_2025",
      "GO_Biological_Process_2025",
      "KEGG_2021_Human",
      "Reactome_Pathways_2024"
    )
  }

  enriched <- enrichr(selected_genes, databases)
}

peptide_map <- load_reference_data()


dpp_files <- list.files("results", "dpp") |>
  set_names(~ .x |> str_extract("-dpp_(.*)\\.csv", 1L))

dpp_names <- names(dpp_files)

run_prefix <- dpp_files |>
  basename() |>
  str_extract("(.*)-dpp.*", 1L)

outfiles <- file.path(
  "results",
  str_glue_data(
    list(
      run_prefix = run_prefix,
      dpp_names = dpp_names
    ),
    "{run_prefix}-{dpp_names}-pathways.xlsx"
  )
)

results <- dpp_files |>
  map(
    ~ perform_pathway_enrichment(.x, peptide_map)
  ) |>
  map2(
    outfiles,
    ~ write_xlsx(
      .x,
      .y
    )
  )
