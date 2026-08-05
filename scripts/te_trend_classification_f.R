# Canonical TE senescence trend classification helpers.
# Sourced by TE-Seq generation and result-comparison code so class labels cannot drift.

te_senescence_non_noise_trend_classes <- function() {
  c(
    "increasing",
    "decreasing",
    "pol_peak",
    "esen_peak",
    "lsen_peak",
    "pol_drop",
    "esen_drop",
    "lsen_drop"
  )
}

te_senescence_noise_classes <- function() {
  c("flat_or_noisy", "unclassified")
}

te_senescence_trend_class_colors <- function(classes = NULL, include_noise = FALSE) {
  class_colors <- c(
    increasing = "#8ccf91",
    decreasing = "#e58f8f",
    pol_peak = "#efc96d",
    esen_peak = "#8dc9d8",
    lsen_peak = "#cdb77d",
    pol_drop = "#b8addf",
    esen_drop = "#b58bd1",
    lsen_drop = "#9fb5d9"
  )
  if (isTRUE(include_noise)) {
    class_colors <- c(
      class_colors,
      flat_or_noisy = "#f2f4f7",
      unclassified = "#ffffff"
    )
  }

  if (!is.null(classes)) {
    classes <- sort(unique(as.character(classes)))
    classes <- classes[!is.na(classes) & nzchar(classes)]
    extra_classes <- setdiff(classes, names(class_colors))
    if (length(extra_classes) > 0L) {
      class_colors <- c(
        class_colors,
        stats::setNames(grDevices::hcl.colors(length(extra_classes), palette = "Dark 3"), extra_classes)
      )
    }
  }
  class_colors
}

te_senescence_is_noise_class <- function(x) {
  is.na(x) | !nzchar(as.character(x)) | as.character(x) %in% te_senescence_noise_classes()
}

te_senescence_trend_delta_signs <- function(trend_class) {
  trend_class <- as.character(trend_class)
  switch(
    trend_class,
    increasing = c(ESEN_vs_PRO = 1, LSEN_vs_ESEN = 1),
    decreasing = c(ESEN_vs_PRO = -1, LSEN_vs_ESEN = -1),
    pol_peak = c(ESEN_vs_PRO = -1, LSEN_vs_PRO = -1),
    esen_peak = c(ESEN_vs_PRO = 1, LSEN_vs_ESEN = -1),
    lsen_peak = c(LSEN_vs_PRO = 1, LSEN_vs_ESEN = 1),
    pol_drop = c(ESEN_vs_PRO = 1, LSEN_vs_PRO = 1),
    esen_drop = c(ESEN_vs_PRO = -1, LSEN_vs_ESEN = 1),
    lsen_drop = c(LSEN_vs_PRO = -1, LSEN_vs_ESEN = -1),
    numeric()
  )
}

te_senescence_trend_p_contrasts <- function(trend_class) {
  trend_class <- as.character(trend_class)
  if (trend_class %in% te_senescence_noise_classes()) {
    return(c("ESEN_vs_PRO", "LSEN_vs_PRO", "LSEN_vs_ESEN"))
  }
  names(te_senescence_trend_delta_signs(trend_class))
}

# Classify each TE subfamily by its PRO -> ESEN -> LSEN activity pattern.
classify_te_senescence_trend <- function(mean_PRO,
                                         mean_ESEN,
                                         mean_LSEN,
                                         delta_threshold) {
  if (!is.numeric(delta_threshold) || length(delta_threshold) != 1L || is.na(delta_threshold) || delta_threshold < 0) {
    stop("delta_threshold must be one non-negative number.", call. = FALSE)
  }

  mean_PRO <- as.numeric(mean_PRO)
  mean_ESEN <- as.numeric(mean_ESEN)
  mean_LSEN <- as.numeric(mean_LSEN)
  if (!identical(length(mean_PRO), length(mean_ESEN)) || !identical(length(mean_PRO), length(mean_LSEN))) {
    stop("mean_PRO, mean_ESEN, and mean_LSEN must have the same length.", call. = FALSE)
  }

  deltas <- data.frame(
    ESEN_vs_PRO = mean_ESEN - mean_PRO,
    LSEN_vs_PRO = mean_LSEN - mean_PRO,
    LSEN_vs_ESEN = mean_LSEN - mean_ESEN
  )

  trend_class <- rep("flat_or_noisy", length(mean_PRO))
  complete_rows <- !(is.na(mean_PRO) | is.na(mean_ESEN) | is.na(mean_LSEN))
  assign_class <- function(class_name) {
    signs <- te_senescence_trend_delta_signs(class_name)
    if (length(signs) == 0L) {
      return(invisible(NULL))
    }
    condition <- complete_rows
    for (contrast_name in names(signs)) {
      values <- deltas[[contrast_name]]
      condition <- condition & if (signs[[contrast_name]] > 0) {
        values >= delta_threshold
      } else {
        values <= -delta_threshold
      }
    }
    idx <- trend_class == "flat_or_noisy" & condition
    trend_class[idx] <<- class_name
    invisible(NULL)
  }

  for (class_name in te_senescence_non_noise_trend_classes()) {
    assign_class(class_name)
  }

  trend_class
}
