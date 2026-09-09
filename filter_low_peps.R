#' Filters out peptides with low signals based on the pw data (max exposure)
#'
#' This function takes in the pw data (max exposure), signal threshold, and samples and group names (optional)
#'
#' @param data krsa pw data (max exposure)
#' @param threshold signal threshold
#' @param samples (optional) sample names
#' @param groups (optional) group names
#'
#' @return vector
#'
#' @family QC functions
#'
#'
#' @export
#'
#' @examples
#' TRUE
mr_krsa_filter_lowPeps <- function(data, threshold, n_threshold_decimal, samples = NULL, groups = NULL) {
  # Filter data based on samples and groups first
  filtered_data <- data %>%
    {
      if (!is.null(samples)) dplyr::filter(., SampleName %in% samples) else .
    } %>%
    {
      if (!is.null(groups)) dplyr::filter(., Group %in% groups) else .
    } %>%
    dplyr::select(-Group)

  # Pivot wider to get Signal intensities by SampleName
  wide_data <- filtered_data %>%
    tidyr::pivot_wider(names_from = SampleName, values_from = Signal)

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
