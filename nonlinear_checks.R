mr_krsa_filter_nonLinear <- function(data, threshold, n_threshold_decimal, samples = NULL, groups = NULL) {
  # Filter data based on samples and groups first
  filtered_data <- data %>%
    {
      if (!is.null(samples)) dplyr::filter(., SampleName %in% samples) else .
    } %>%
    {
      if (!is.null(groups)) dplyr::filter(., Group %in% groups) else .
    } %>%
    dplyr::select(SampleName, Peptide, r.seq)

  # Pivot wider to get r.seq values by SampleName
  wide_data <- filtered_data %>%
    tidyr::pivot_wider(names_from = SampleName, values_from = r.seq)

  # Determine the actual number of samples based on the decimal threshold
  # Exclude the 'Peptide' column to count only sample columns
  num_sample_cols <- ncol(wide_data) - 1 # Subtract 1 for the 'Peptide' column
  n_samples_required <- ceiling(n_threshold_decimal * num_sample_cols)

  # Apply the filtering
  p <- wide_data %>%
    dplyr::filter(rowSums(dplyr::select(., -Peptide) >= threshold) >= n_samples_required) %>%
    dplyr::pull(Peptide)

  message(paste("Filtered out", length(data$Peptide %>% unique()) - length(p), "Peptides"))

  p
}

