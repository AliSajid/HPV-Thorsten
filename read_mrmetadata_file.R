# check_mrmetadata_file.R
library(tidyverse)
read_metadata_file <- function(filepath, pairing_variable_columnname=NULL, run_id_columnname=NULL, median_sigmb_columnname=NULL, signal_saturation_columnname=NULL) {
  tryCatch({
    metadata <- read_csv(filepath)
    message("Metadata file imported successfully: OK")
  }, error = function(e) {
    message(paste("ERROR:", e$message))
  })
  required_vars <- c(run_id_columnname=run_id_columnname, median_sigmb_columnname=median_sigmb_columnname, signal_saturation_columnname=signal_saturation_columnname)
  if (any(is.null(required_vars))) {
    variable_names <- names(required_vars)
    null_variables <- variable_names[which(!is.null(required_vars))]
    stop(paste("ERROR: the variables", paste(null_variables, sep=", "), "are not defined or NULL. Please change these to column names that are present in the metadata file."))
  }
  if (!is.null(pairing_variable_columnname)) {
      metadata <- metadata |> rename(pairing_variable:=!!sym(pairing_variable_columnname), median_sigmb_file:=!!sym(median_sigmb_columnname), signal_saturation_file:=!!sym(signal_saturation_columnname), run_id:=!!sym(run_id_columnname))
  } else {
    metadata <- metadata |> rename(median_sigmb_file:=!!sym(median_sigmb_columnname), signal_saturation_file:=!!sym(signal_saturation_columnname), run_id:=!!sym(run_id_columnname))
  }
  invalid_metadata <- metadata |>
    group_by(run_id) |>
    summarize(
      n_distinct_medsigmb = n_distinct(median_sigmb_file),
      n_distinct_sigsatfile = n_distinct(signal_saturation_file)
    ) |> ungroup() |> filter(unique(n_distinct_medsigmb)!=1 || unique(n_distinct_sigsatfile)!=1)
  if (nrow(invalid_metadata)>0 & any(invalid_metadata$n_distinct_medsigmb!=1)) {
    invalid_medsigmb_runs <- invalid_metadata |> filter(n_distinct_medsigmb!=1) |> pull(run_id) |> unique()
    for (invalid_run in unique(invalid_medsigmb_runs)) {
      warning(paste(invalid_run, "has", invalid_metadata |> filter(run_id==invalid_run) |> pull(n_distinct_medsigmb) |> unique(), "unique medsigmb filenames when it should only have 1."))
    }
  }
  if (nrow(invalid_metadata)>0 & any(invalid_metadata$n_distinct_sigsatfile!=1)) {
    invalid_sigsat_runs <- invalid_metadata |> filter(n_distinct_sigsatfile!=1) |> pull(run_id) |> unique()
    for (invalid_run in unique(invalid_sigsat_runs)) {
      warning(paste(invalid_run, "has", invalid_metadata |> filter(run_id==invalid_run) |> pull(n_distinct_sigsatfile) |> unique(), "unique sigsat filenames when it should only have 1."))
    }
  }
  if (nrow(invalid_metadata)>0) {
    stop("Invalid number of unique files per run. See warnings() for details.")
  }
  return(metadata)
}
