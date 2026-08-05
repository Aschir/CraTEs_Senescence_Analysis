# TSV loading and validation helpers for result_comparison.ipynb.
# The comparison notebook is intentionally TSV-only and does not recalculate
# craTEs activities, significance models, or target matrices.

comparison_log <- function(message, level = "INFO", quiet = FALSE) {
  msg <- sprintf("[%s] %s", toupper(level), message)
  if (!isTRUE(quiet)) {
    message(msg)
  }
  invisible(msg)
}

comparison_warn <- function(message) {
  warning(sprintf("[comparison] %s", message), call. = FALSE, immediate. = TRUE)
  invisible(message)
}

comparison_required_files <- function() {
  c(
    activity_fit_summary = "activity_fit_summary.tsv",
    activity_scores = "activity_scores.tsv"
  )
}

comparison_optional_files <- function() {
  c(
    te_contrast_summary = "te_contrast_summary.tsv",
    target_summary = "target_summary.tsv"
  )
}

comparison_diagnostic_files <- function() {
  c(
    sample_activity_diagnostics = "sample_activity_diagnostics.tsv",
    activity_lambda_diagnostics = "activity_lambda_diagnostics.tsv",
    signif_model_diagnostics = "signif_model_diagnostics.tsv"
  )
}

comparison_read_tsv <- function(path, required = TRUE, table_name = basename(path), warn_missing = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) {
      stop(
        sprintf("Required TSV is missing for `%s`: %s", table_name, path),
        call. = FALSE
      )
    }
    if (isTRUE(warn_missing)) {
      comparison_warn(sprintf("Optional TSV is missing for `%s`: %s", table_name, path))
    }
    return(NULL)
  }

  out <- tryCatch(
    {
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
      } else {
        utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
      }
    },
    error = function(e) {
      stop(sprintf("Could not read `%s` from %s: %s", table_name, path, conditionMessage(e)), call. = FALSE)
    }
  )

  if (!is.data.frame(out)) {
    stop(sprintf("Loaded `%s` is not a data frame: %s", table_name, path), call. = FALSE)
  }
  out
}

comparison_require_columns <- function(x, required_cols, table_name, result_dir) {
  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "`%s` in %s is missing required column(s): %s",
        table_name,
        result_dir,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

comparison_coerce_numeric_columns <- function(x, cols, table_name, result_dir) {
  cols <- intersect(cols, colnames(x))
  for (col in cols) {
    if (is.numeric(x[[col]])) {
      next
    }
    original <- x[[col]]
    converted <- suppressWarnings(as.numeric(original))
    bad <- !is.na(original) & is.na(converted)
    if (any(bad)) {
      stop(
        sprintf(
          "`%s` in %s has non-numeric values in column `%s`.",
          table_name,
          result_dir,
          col
        ),
        call. = FALSE
      )
    }
    x[[col]] <- converted
  }
  x
}

comparison_validate_activity_fit_summary <- function(x, result_dir) {
  comparison_require_columns(
    x,
    required_cols = c("scope", "condition"),
    table_name = "activity_fit_summary.tsv",
    result_dir = result_dir
  )
  numeric_cols <- intersect(
    c("n_samples", "mean_r2", "mean_r2_ratio", "mean_sparsity_percent", "lambda_best"),
    colnames(x)
  )
  comparison_coerce_numeric_columns(x, numeric_cols, "activity_fit_summary.tsv", result_dir)
}

comparison_validate_sample_activity_diagnostics <- function(x, result_dir) {
  comparison_require_columns(
    x,
    required_cols = c("sample_id", "condition"),
    table_name = "sample_activity_diagnostics.tsv",
    result_dir = result_dir
  )
  numeric_cols <- intersect(
    c(
      "r2",
      "r2_max",
      "r2_ratio",
      "lambda_best",
      "lambda_max",
      "lambda_min_overall",
      "lambda_best_to_global_ratio",
      "zero_fraction"
    ),
    colnames(x)
  )
  comparison_coerce_numeric_columns(x, numeric_cols, "sample_activity_diagnostics.tsv", result_dir)
}

comparison_validate_activity_lambda_diagnostics <- function(x, result_dir) {
  comparison_require_columns(
    x,
    required_cols = "lambda",
    table_name = "activity_lambda_diagnostics.tsv",
    result_dir = result_dir
  )
  numeric_cols <- intersect(
    c(
      "lambda_index",
      "lambda",
      "mean_crossval_mse",
      "sd_crossval_mse",
      "global_lambda_min",
      "global_lambda_max",
      "global_lambda_grid_min"
    ),
    colnames(x)
  )
  x <- comparison_coerce_numeric_columns(x, numeric_cols, "activity_lambda_diagnostics.tsv", result_dir)
  if ("is_global_lambda_min" %in% colnames(x) && !is.logical(x$is_global_lambda_min)) {
    x$is_global_lambda_min <- tolower(as.character(x$is_global_lambda_min)) %in% c("true", "t", "1", "yes")
  }
  x
}

comparison_validate_signif_model_diagnostics <- function(x, result_dir) {
  comparison_require_columns(
    x,
    required_cols = "condition",
    table_name = "signif_model_diagnostics.tsv",
    result_dir = result_dir
  )
  numeric_cols <- intersect(
    c("f_statistic", "f_numdf", "f_dendf", "f_p_value", "r_squared", "rss", "n_aliased"),
    colnames(x)
  )
  comparison_coerce_numeric_columns(x, numeric_cols, "signif_model_diagnostics.tsv", result_dir)
}

comparison_validate_activity_scores <- function(x, result_dir) {
  comparison_require_columns(
    x,
    required_cols = "te_subfamily",
    table_name = "activity_scores.tsv",
    result_dir = result_dir
  )
  sample_cols <- setdiff(colnames(x), "te_subfamily")
  if (length(sample_cols) == 0L) {
    stop(sprintf("`activity_scores.tsv` in %s has no activity-score columns.", result_dir), call. = FALSE)
  }
  comparison_coerce_numeric_columns(x, sample_cols, "activity_scores.tsv", result_dir)
}

comparison_validate_te_contrast_summary <- function(x, result_dir) {
  comparison_require_columns(
    x,
    required_cols = "te_subfamily",
    table_name = "te_contrast_summary.tsv",
    result_dir = result_dir
  )
  numeric_cols <- grep(
    "^(mean_|delta_|padj_|p_adj|p_value|max_abs_delta|n_significant_contrasts|statistic$|t_statistic$|Statistic$)",
    colnames(x),
    value = TRUE
  )
  comparison_coerce_numeric_columns(x, numeric_cols, "te_contrast_summary.tsv", result_dir)
}

comparison_signif_activity_paths <- function(result_dir) {
  paths <- list.files(result_dir, pattern = "_signif_act\\.tsv$", full.names = TRUE)
  if (length(paths) == 0L) {
    return(stats::setNames(character(), character()))
  }
  contrast_ids <- sub("_signif_act\\.tsv$", "", basename(paths))
  stats::setNames(paths, contrast_ids)
}

comparison_validate_signif_activity <- function(x, result_dir, table_name) {
  comparison_require_columns(
    x,
    required_cols = c("te_subfamily", "p_value"),
    table_name = table_name,
    result_dir = result_dir
  )
  if (anyDuplicated(x$te_subfamily)) {
    duplicated_tes <- unique(x$te_subfamily[duplicated(x$te_subfamily)])
    stop(
      sprintf("`%s` in %s contains duplicated TE subfamilies: %s", table_name, result_dir, paste(duplicated_tes, collapse = ", ")),
      call. = FALSE
    )
  }
  comparison_coerce_numeric_columns(x, "p_value", table_name, result_dir)
}

comparison_add_raw_p_values_to_te_contrast <- function(te_contrast_summary, signif_activity_tables) {
  if (is.null(te_contrast_summary) || length(signif_activity_tables) == 0L) {
    return(te_contrast_summary)
  }
  te_ids <- te_contrast_summary$te_subfamily
  for (contrast_id in names(signif_activity_tables)) {
    signif_table <- signif_activity_tables[[contrast_id]]
    if (is.null(signif_table) || !"p_value" %in% colnames(signif_table)) {
      next
    }
    p_col <- paste0("p_value_", contrast_id)
    p_lookup <- stats::setNames(signif_table$p_value, signif_table$te_subfamily)
    p_values <- unname(p_lookup[te_ids])
    if (!p_col %in% colnames(te_contrast_summary)) {
      te_contrast_summary[[p_col]] <- p_values
    } else {
      missing_values <- is.na(te_contrast_summary[[p_col]])
      te_contrast_summary[[p_col]][missing_values] <- p_values[missing_values]
    }
  }
  te_contrast_summary
}

comparison_validate_target_summary <- function(x, result_dir) {
  if (!any(c("row_id", "gene_id", "ensembl_gene_id", "hgnc_symbol") %in% colnames(x))) {
    comparison_warn(sprintf("`target_summary.tsv` in %s has no obvious gene identifier column.", result_dir))
  }
  numeric_cols <- colnames(x)[vapply(x, function(col) !is.character(col) && !is.logical(col), logical(1))]
  comparison_coerce_numeric_columns(x, numeric_cols, "target_summary.tsv", result_dir)
}

comparison_validate_target_gene_matrix <- function(x, result_dir, table_name) {
  gene_col <- intersect(c("row_id", "gene_id", "ensembl_gene_id"), colnames(x))
  if (length(gene_col) == 0L) {
    gene_col <- colnames(x)[1]
    comparison_warn(sprintf(
      "`%s` in %s has no row_id/gene_id/ensembl_gene_id column; using first column `%s` as gene ID.",
      table_name,
      result_dir,
      gene_col
    ))
  } else {
    gene_col <- gene_col[1]
  }

  if (!identical(gene_col, "row_id")) {
    colnames(x)[match(gene_col, colnames(x))] <- "row_id"
  }

  te_cols <- setdiff(colnames(x), "row_id")
  if (length(te_cols) == 0L) {
    stop(sprintf("`%s` in %s has no TE-subfamily columns.", table_name, result_dir), call. = FALSE)
  }
  comparison_coerce_numeric_columns(x, te_cols, table_name, result_dir)
}

comparison_statistic_columns <- function(x) {
  stat_like <- c("statistic", "t_statistic", "Statistic", "T_statistic", "T_Statistic")
  intersect(colnames(x), stat_like)
}

comparison_validate_statistic_column_consistency <- function(comparison_results) {
  stat_cols <- character()
  locations <- character()

  for (result in comparison_results) {
    contrast <- result$tables$te_contrast_summary
    if (is.null(contrast)) {
      next
    }
    cols <- comparison_statistic_columns(contrast)
    if (length(cols) > 0L) {
      stat_cols <- c(stat_cols, cols)
      locations <- c(locations, sprintf("%s: %s", result$result_id, paste(cols, collapse = ", ")))
    }
  }

  variants <- unique(stat_cols)
  if (length(variants) > 1L) {
    stop(
      paste(
        "Selected result directories use inconsistent statistic column names in `te_contrast_summary.tsv`:",
        paste(locations, collapse = "; "),
        "Make the TSV column names consistent before comparison.",
        sep = " "
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

comparison_validate_activity_fit_schema_consistency <- function(comparison_results) {
  has_fit <- vapply(comparison_results, function(result) !is.null(result$tables$activity_fit_summary), logical(1))
  if (sum(has_fit) < 2L) {
    return(invisible(TRUE))
  }

  schemas <- lapply(comparison_results[has_fit], function(result) colnames(result$tables$activity_fit_summary))
  reference <- schemas[[1]]
  mismatched <- vapply(schemas, function(cols) !identical(cols, reference), logical(1))

  if (any(mismatched)) {
    bad <- names(comparison_results[has_fit])[mismatched]
    stop(
      sprintf(
        "`activity_fit_summary.tsv` schemas differ across selected result directories. Mismatched result(s): %s",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

comparison_validate_loaded_set <- function(comparison_results) {
  if (length(comparison_results) == 0L) {
    stop("No comparison results were loaded.", call. = FALSE)
  }
  comparison_validate_activity_fit_schema_consistency(comparison_results)
  comparison_validate_statistic_column_consistency(comparison_results)
  invisible(TRUE)
}

discover_result_dirs <- function(results_parent = "results") {
  if (!dir.exists(results_parent)) {
    stop(sprintf("Results parent directory does not exist: %s", results_parent), call. = FALSE)
  }

  dirs <- list.dirs(results_parent, recursive = FALSE, full.names = TRUE)
  required_files <- comparison_required_files()
  keep <- vapply(
    dirs,
    function(dir) all(file.exists(file.path(dir, required_files))),
    logical(1)
  )
  dirs[keep]
}

load_result_dir <- function(result_dir, required = TRUE, warn_missing_optional = FALSE) {
  if (!dir.exists(result_dir)) {
    stop(sprintf("Result directory does not exist: %s", result_dir), call. = FALSE)
  }

  result_dir <- normalizePath(result_dir, winslash = "/", mustWork = TRUE)
  result_id <- basename(result_dir)
  required_files <- comparison_required_files()
  optional_files <- comparison_optional_files()
  diagnostic_files <- comparison_diagnostic_files()
  loaded_files <- data.frame(
    table = character(),
    file = character(),
    required = logical(),
    loaded = logical(),
    stringsAsFactors = FALSE
  )

  tables <- list()
  missing_required <- character()
  missing_optional <- character()

  for (table_name in names(required_files)) {
    path <- file.path(result_dir, required_files[[table_name]])
    if (!file.exists(path)) {
      missing_required <- c(missing_required, basename(path))
    }
    tables[[table_name]] <- comparison_read_tsv(path, required = required, table_name = basename(path))
    loaded_files <- rbind(
      loaded_files,
      data.frame(table = table_name, file = path, required = TRUE, loaded = !is.null(tables[[table_name]]), stringsAsFactors = FALSE)
    )
  }

  if (!is.null(tables$activity_fit_summary)) {
    tables$activity_fit_summary <- comparison_validate_activity_fit_summary(tables$activity_fit_summary, result_dir)
  }
  if (!is.null(tables$activity_scores)) {
    tables$activity_scores <- comparison_validate_activity_scores(tables$activity_scores, result_dir)
  }

  for (table_name in names(optional_files)) {
    path <- file.path(result_dir, optional_files[[table_name]])
    if (!file.exists(path)) {
      missing_optional <- c(missing_optional, basename(path))
    }
    tables[[table_name]] <- comparison_read_tsv(
      path,
      required = FALSE,
      table_name = basename(path),
      warn_missing = warn_missing_optional
    )
    loaded_files <- rbind(
      loaded_files,
      data.frame(table = table_name, file = path, required = FALSE, loaded = !is.null(tables[[table_name]]), stringsAsFactors = FALSE)
    )
  }

  if (!is.null(tables$te_contrast_summary)) {
    signif_paths <- comparison_signif_activity_paths(result_dir)
    signif_activity_tables <- list()
    for (contrast_id in names(signif_paths)) {
      path <- signif_paths[[contrast_id]]
      table_name <- basename(path)
      signif_activity_tables[[contrast_id]] <- comparison_read_tsv(
        path,
        required = FALSE,
        table_name = table_name,
        warn_missing = FALSE
      )
      signif_activity_tables[[contrast_id]] <- comparison_validate_signif_activity(
        signif_activity_tables[[contrast_id]],
        result_dir = result_dir,
        table_name = table_name
      )
      loaded_files <- rbind(
        loaded_files,
        data.frame(table = paste0(contrast_id, "_signif_act"), file = path, required = FALSE, loaded = TRUE, stringsAsFactors = FALSE)
      )
    }
    tables$te_contrast_summary <- comparison_add_raw_p_values_to_te_contrast(
      tables$te_contrast_summary,
      signif_activity_tables
    )
    tables$te_contrast_summary <- comparison_validate_te_contrast_summary(tables$te_contrast_summary, result_dir)
    tables$signif_activities <- signif_activity_tables
  }
  if (!is.null(tables$target_summary)) {
    tables$target_summary <- comparison_validate_target_summary(tables$target_summary, result_dir)
  }

  for (table_name in names(diagnostic_files)) {
    path <- file.path(result_dir, diagnostic_files[[table_name]])
    tables[[table_name]] <- comparison_read_tsv(
      path,
      required = FALSE,
      table_name = basename(path),
      warn_missing = FALSE
    )
    loaded_files <- rbind(
      loaded_files,
      data.frame(table = table_name, file = path, required = FALSE, loaded = !is.null(tables[[table_name]]), stringsAsFactors = FALSE)
    )
  }

  if (!is.null(tables$sample_activity_diagnostics)) {
    tables$sample_activity_diagnostics <- comparison_validate_sample_activity_diagnostics(tables$sample_activity_diagnostics, result_dir)
  }
  if (!is.null(tables$activity_lambda_diagnostics)) {
    tables$activity_lambda_diagnostics <- comparison_validate_activity_lambda_diagnostics(tables$activity_lambda_diagnostics, result_dir)
  }
  if (!is.null(tables$signif_model_diagnostics)) {
    tables$signif_model_diagnostics <- comparison_validate_signif_model_diagnostics(tables$signif_model_diagnostics, result_dir)
  }

  target_paths <- list.files(result_dir, pattern = "^target_gene_matrix_.*\\.tsv$", full.names = TRUE)
  target_gene_matrices <- list()
  if (length(target_paths) == 0L) {
    missing_optional <- c(missing_optional, "target_gene_matrix_*.tsv")
  } else {
    for (path in target_paths) {
      matrix_name <- sub("\\.tsv$", "", basename(path))
      target_gene_matrices[[matrix_name]] <- comparison_read_tsv(path, required = FALSE, table_name = basename(path))
      target_gene_matrices[[matrix_name]] <- comparison_validate_target_gene_matrix(
        target_gene_matrices[[matrix_name]],
        result_dir = result_dir,
        table_name = basename(path)
      )
      loaded_files <- rbind(
        loaded_files,
        data.frame(table = matrix_name, file = path, required = FALSE, loaded = TRUE, stringsAsFactors = FALSE)
      )
    }
  }

  structure(
    list(
      result_dir = result_dir,
      result_id = result_id,
      tables = tables,
      target_gene_matrices = target_gene_matrices,
      loaded_files = loaded_files,
      missing_required = unique(missing_required),
      missing_optional = unique(missing_optional),
      skipped_files = character()
    ),
    class = "comparison_result"
  )
}

load_comparison_results <- function(result_dirs, required = TRUE, warn_missing_optional = FALSE) {
  if (length(result_dirs) < 2L) {
    stop("Select at least two result directories for comparison.", call. = FALSE)
  }

  comparison_log(sprintf("Loading %d result directories", length(result_dirs)))
  results <- lapply(
    result_dirs,
    load_result_dir,
    required = required,
    warn_missing_optional = warn_missing_optional
  )
  names(results) <- vapply(results, function(result) result$result_id, character(1))
  results <- structure(results, class = "comparison_results")
  comparison_validate_loaded_set(results)
  comparison_log("Loaded comparison result set")
  results
}
