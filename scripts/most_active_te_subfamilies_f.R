# Helpers for ranking inferred craTEs activity against TE expression.

if (!exists("resolve_duplicates", mode = "function") ||
    !exists("select_top_te_subfamilies_by_activity_change", mode = "function")) {
  helper_path <- file.path("scripts", "craTEs_own_function.R")
  if (file.exists(helper_path)) {
    source(helper_path)
  }
}

.most_active_validate_k_values <- function(k_values) {
  if (!is.numeric(k_values) || length(k_values) == 0L ||
      any(is.na(k_values)) || any(k_values <= 0)) {
    stop("k_values must be a non-empty vector of positive numbers.", call. = FALSE)
  }
  unique(as.integer(k_values))
}

.most_active_normalize_duplicate_method <- function(method) {
  if (is.null(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("duplicate_method must be one non-empty character value.", call. = FALSE)
  }

  method_key <- tolower(gsub("[ -]+", "_", as.character(method)))
  method_aliases <- c(
    sum = "sum",
    mean = "mean",
    average = "mean",
    max_sum = "max_sum",
    max_variance = "max_var",
    max_var = "max_var",
    discard = "discard",
    identical = "identical"
  )

  if (!method_key %in% names(method_aliases)) {
    stop(
      sprintf(
        "duplicate_method must be one of: %s",
        paste(unique(method_aliases), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unname(method_aliases[[method_key]])
}

.most_active_numeric_matrix <- function(raw_values,
                                        row_ids,
                                        col_ids,
                                        object_name = "matrix") {
  value_raw <- as.matrix(raw_values)
  value_numeric <- suppressWarnings(as.numeric(value_raw))
  out <- matrix(
    value_numeric,
    nrow = nrow(value_raw),
    dimnames = list(row_ids, col_ids)
  )

  conversion_failed <- is.na(out) & !is.na(value_raw) & value_raw != ""
  if (any(conversion_failed)) {
    failed_idx <- which(conversion_failed, arr.ind = TRUE)
    examples <- utils::head(data.frame(
      feature = rownames(out)[failed_idx[, "row"]],
      sample = colnames(out)[failed_idx[, "col"]],
      value = value_raw[failed_idx],
      stringsAsFactors = FALSE
    ), 10)
    stop(
      paste(
        sprintf("%s contains non-numeric count values.", object_name),
        "Examples:",
        paste(utils::capture.output(print(examples, row.names = FALSE)), collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  storage.mode(out) <- "double"
  out
}

.most_active_numeric_like_column <- function(x) {
  converted <- suppressWarnings(as.numeric(x))
  observed <- !(is.na(x) | x == "")
  any(observed) && all(!observed | !is.na(converted))
}

extract_te_subfamily_from_feature_id <- function(feature_ids,
                                                 reference_te_subfamilies = NULL) {
  feature_ids <- as.character(feature_ids)
  reference_te_subfamilies <- unique(as.character(reference_te_subfamilies))
  reference_te_subfamilies <- reference_te_subfamilies[!is.na(reference_te_subfamilies) & nzchar(reference_te_subfamilies)]

  parsed <- feature_ids
  parsed <- sub("_((chr)?[0-9XYM]+|[0-9]+[pq]).*$", "", parsed, ignore.case = TRUE)

  if (length(reference_te_subfamilies) > 0L) {
    exact_match <- feature_ids %in% reference_te_subfamilies
    parsed[exact_match] <- feature_ids[exact_match]
  }

  parsed
}

read_raw_expression_count_matrix <- function(E_path,
                                             id_col = 1,
                                             sample_cols = NULL,
                                             sample_pattern = "^(GSM|PRO|ESEN|LSEN)") {
  if (is.null(E_path) || length(E_path) != 1L || !nzchar(E_path)) {
    stop("E_path must be one non-empty file path.", call. = FALSE)
  }
  if (!file.exists(E_path)) {
    stop(sprintf("Raw E matrix does not exist: %s", E_path), call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("read_raw_expression_count_matrix() requires package 'data.table'.", call. = FALSE)
  }

  E_table <- data.table::fread(E_path, sep = "\t", header = TRUE, check.names = FALSE)
  if (ncol(E_table) < 2L) {
    stop("Raw E matrix must contain one feature-ID column and at least one count column.", call. = FALSE)
  }

  if (is.character(id_col)) {
    if (!id_col %in% colnames(E_table)) {
      stop(sprintf("id_col '%s' was not found in raw E matrix.", id_col), call. = FALSE)
    }
    id_index <- match(id_col, colnames(E_table))
  } else {
    id_index <- as.integer(id_col)
    if (length(id_index) != 1L || is.na(id_index) ||
        id_index < 1L || id_index > ncol(E_table)) {
      stop("id_col must be a valid column name or one-based column index.", call. = FALSE)
    }
  }

  row_ids <- as.character(E_table[[id_index]])
  if (any(is.na(row_ids) | row_ids == "")) {
    stop("Raw E feature-ID column contains missing or empty values.", call. = FALSE)
  }

  if (is.null(sample_cols)) {
    sample_cols <- grep(sample_pattern, colnames(E_table), value = TRUE)
  }
  if (length(sample_cols) == 0L) {
    non_id_cols <- setdiff(colnames(E_table), colnames(E_table)[id_index])
    numeric_like <- vapply(E_table[, ..non_id_cols], .most_active_numeric_like_column, logical(1))
    sample_cols <- non_id_cols[numeric_like]
  }
  sample_cols <- unique(as.character(sample_cols))
  if (length(sample_cols) == 0L) {
    stop("No sample/count columns could be identified in raw E matrix.", call. = FALSE)
  }
  missing_sample_cols <- setdiff(sample_cols, colnames(E_table))
  if (length(missing_sample_cols) > 0L) {
    stop(
      sprintf("sample_cols not found in raw E matrix: %s", paste(missing_sample_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  count_matrix <- .most_active_numeric_matrix(
    E_table[, ..sample_cols],
    row_ids = row_ids,
    col_ids = sample_cols,
    object_name = "Raw E matrix"
  )

  list(
    counts = count_matrix,
    feature_table = as.data.frame(E_table, check.names = FALSE, stringsAsFactors = FALSE),
    id_col = colnames(E_table)[id_index],
    sample_cols = sample_cols
  )
}

load_activity_scores_matrix <- function(activity_scores_path,
                                        row_id_col = "te_subfamily") {
  if (is.null(activity_scores_path) || length(activity_scores_path) != 1L ||
      !nzchar(activity_scores_path)) {
    stop("activity_scores_path must be one non-empty file path.", call. = FALSE)
  }
  if (!file.exists(activity_scores_path)) {
    stop(sprintf("Activity-score TSV does not exist: %s", activity_scores_path), call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("load_activity_scores_matrix() requires package 'data.table'.", call. = FALSE)
  }

  activity_table <- data.table::fread(
    activity_scores_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    data.table = FALSE
  )
  if (!row_id_col %in% colnames(activity_table)) {
    stop(sprintf("Activity-score TSV is missing row ID column '%s'.", row_id_col), call. = FALSE)
  }

  te_ids <- as.character(activity_table[[row_id_col]])
  if (any(is.na(te_ids) | te_ids == "")) {
    stop("Activity-score TSV contains missing or empty TE subfamily IDs.", call. = FALSE)
  }
  if (anyDuplicated(te_ids)) {
    duplicated_tes <- unique(te_ids[duplicated(te_ids)])
    stop(
      sprintf("Activity-score TSV contains duplicated TE subfamilies: %s", paste(duplicated_tes, collapse = ", ")),
      call. = FALSE
    )
  }

  sample_cols <- setdiff(colnames(activity_table), row_id_col)
  if (length(sample_cols) == 0L) {
    stop("Activity-score TSV has no sample columns.", call. = FALSE)
  }

  activity_scores <- .most_active_numeric_matrix(
    activity_table[, sample_cols, drop = FALSE],
    row_ids = te_ids,
    col_ids = sample_cols,
    object_name = "Activity-score TSV"
  )
  activity_scores
}

infer_senescence_condition_groups <- function(sample_ids,
                                              condition_prefixes = c("PRO", "ESEN", "LSEN")) {
  sample_ids <- as.character(sample_ids)
  sample_ids <- sample_ids[nzchar(sample_ids)]
  if (length(sample_ids) == 0L) {
    stop("sample_ids must contain at least one non-empty sample name.", call. = FALSE)
  }

  condition_groups <- stats::setNames(
    lapply(condition_prefixes, function(prefix) {
      grep(paste0("^", prefix), sample_ids, value = TRUE)
    }),
    condition_prefixes
  )

  missing_conditions <- names(condition_groups)[lengths(condition_groups) == 0L]
  if (length(missing_conditions) > 0L) {
    stop(
      sprintf(
        "Could not infer senescence samples for condition(s): %s",
        paste(missing_conditions, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  condition_groups
}

condition_mean_matrix <- function(values,
                                  condition_groups,
                                  condition_order = c("PRO", "ESEN", "LSEN"),
                                  object_name = "matrix") {
  if (!is.matrix(values) && !is.data.frame(values)) {
    stop(sprintf("%s must be a matrix or data.frame.", object_name), call. = FALSE)
  }
  values <- as.matrix(values)
  if (!is.numeric(values)) {
    stop(sprintf("%s must be numeric.", object_name), call. = FALSE)
  }
  if (is.null(rownames(values)) || any(!nzchar(rownames(values)))) {
    stop(sprintf("%s must have non-empty rownames.", object_name), call. = FALSE)
  }
  if (is.null(colnames(values)) || any(!nzchar(colnames(values)))) {
    stop(sprintf("%s must have non-empty colnames.", object_name), call. = FALSE)
  }
  if (!is.list(condition_groups) || is.null(names(condition_groups)) || any(!nzchar(names(condition_groups)))) {
    stop("condition_groups must be a named list of sample vectors.", call. = FALSE)
  }

  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[condition_order %in% names(condition_groups)]
  if (length(condition_order) == 0L) {
    stop("condition_order must contain at least one condition present in condition_groups.", call. = FALSE)
  }

  condition_groups <- condition_groups[condition_order]
  missing_samples <- setdiff(unique(unlist(condition_groups, use.names = FALSE)), colnames(values))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("Samples not found in %s: %s", object_name, paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- do.call(cbind, lapply(condition_groups, function(samples) {
    rowMeans(values[, samples, drop = FALSE], na.rm = TRUE)
  }))
  colnames(out) <- condition_order
  rownames(out) <- rownames(values)
  storage.mode(out) <- "double"
  out
}

load_reference_te_subfamilies <- function(N_path = NULL,
                                          activity_scores = NULL) {
  reference_tes <- character()

  if (!is.null(activity_scores)) {
    if (!is.matrix(activity_scores) && !is.data.frame(activity_scores)) {
      stop("activity_scores must be a matrix or data.frame when supplied.", call. = FALSE)
    }
    reference_tes <- c(reference_tes, rownames(activity_scores))
  }

  if (!is.null(N_path) && length(N_path) == 1L && nzchar(N_path) && file.exists(N_path)) {
    if (exists("load_count_association_matrix", mode = "function")) {
      N_loaded <- load_count_association_matrix(N_path)
      reference_tes <- c(reference_tes, colnames(N_loaded$N))
    } else {
      N_table <- data.table::fread(N_path, sep = "\t", header = TRUE, check.names = FALSE)
      reference_tes <- c(reference_tes, setdiff(colnames(N_table), colnames(N_table)[1]))
    }
  }

  reference_tes <- unique(as.character(reference_tes))
  reference_tes[!is.na(reference_tes) & nzchar(reference_tes)]
}

load_raw_te_subfamily_counts <- function(E_path,
                                         reference_te_subfamilies = NULL,
                                         duplicate_method = "max_sum",
                                         id_col = 1,
                                         sample_cols = NULL,
                                         sample_pattern = "^(GSM|PRO|ESEN|LSEN)",
                                         te_type_col_candidates = c(
                                           "gene_or_te",
                                           "feature_type",
                                           "feature",
                                           "type",
                                           "Source",
                                           "source"
                                         ),
                                         te_type_values = c(
                                           "te",
                                           "tes",
                                           "repeat",
                                           "repeats",
                                           "transposable_element",
                                           "transposable_elements"
                                         ),
                                         verbose = TRUE) {
  if (!exists("resolve_duplicates", mode = "function") ||
      !exists("duplicate_stats", mode = "function")) {
    stop("load_raw_te_subfamily_counts() requires resolve_duplicates() and duplicate_stats(). Source scripts/craTEs_own_function.R first.", call. = FALSE)
  }

  duplicate_method <- .most_active_normalize_duplicate_method(duplicate_method)
  raw <- read_raw_expression_count_matrix(
    E_path = E_path,
    id_col = id_col,
    sample_cols = sample_cols,
    sample_pattern = sample_pattern
  )

  row_ids <- rownames(raw$counts)
  reference_te_subfamilies <- unique(as.character(reference_te_subfamilies))
  reference_te_subfamilies <- reference_te_subfamilies[!is.na(reference_te_subfamilies) & nzchar(reference_te_subfamilies)]
  parsed_subfamilies <- extract_te_subfamily_from_feature_id(
    row_ids,
    reference_te_subfamilies = reference_te_subfamilies
  )

  reference_keep <- if (length(reference_te_subfamilies) > 0L) {
    row_ids %in% reference_te_subfamilies | parsed_subfamilies %in% reference_te_subfamilies
  } else {
    rep(FALSE, length(row_ids))
  }

  type_col <- intersect(te_type_col_candidates, colnames(raw$feature_table))[1]
  type_keep <- rep(FALSE, length(row_ids))
  if (!is.na(type_col)) {
    type_values <- tolower(gsub("[^a-z0-9]+", "_", as.character(raw$feature_table[[type_col]])))
    type_targets <- tolower(gsub("[^a-z0-9]+", "_", te_type_values))
    type_keep <- type_values %in% type_targets |
      grepl("(^|_)te($|_)|transposable|repeat", type_values)
  }

  keep <- reference_keep | type_keep
  if (!any(keep)) {
    stop(
      paste(
        "No TE-subfamily rows were detected in the raw E matrix.",
        "Provide reference_te_subfamilies matching the row IDs or a TE marker column such as gene_or_te.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  raw_te_counts <- raw$counts[keep, , drop = FALSE]
  rownames(raw_te_counts) <- parsed_subfamilies[keep]
  duplicate_before <- duplicate_stats(rownames(raw_te_counts), raw_te_counts, include_table = FALSE)
  resolved_te_counts <- resolve_duplicates(
    raw_te_counts,
    method = duplicate_method,
    verbose = verbose
  )
  duplicate_after <- duplicate_stats(rownames(resolved_te_counts), resolved_te_counts, include_table = FALSE)

  selection_table <- data.frame(
    feature_id = row_ids[keep],
    te_subfamily = parsed_subfamilies[keep],
    matched_reference = reference_keep[keep],
    matched_te_type_column = type_keep[keep],
    stringsAsFactors = FALSE
  )

  selection_summary <- data.frame(
    raw_rows_total = length(row_ids),
    raw_rows_kept_as_te = sum(keep),
    rows_kept_by_reference = sum(reference_keep),
    rows_kept_by_type_column = sum(type_keep),
    unique_te_subfamilies_before_duplicate_handling = length(unique(parsed_subfamilies[keep])),
    duplicate_groups_before = duplicate_before$duplicate_groups,
    redundant_rows_before = duplicate_before$redundant_rows,
    te_rows_after_duplicate_handling = nrow(resolved_te_counts),
    duplicate_groups_after = duplicate_after$duplicate_groups,
    duplicate_method = duplicate_method,
    stringsAsFactors = FALSE
  )

  list(
    counts = resolved_te_counts,
    counts_before_duplicate_handling = raw_te_counts,
    selection_table = selection_table,
    selection_summary = selection_summary,
    sample_cols = raw$sample_cols,
    id_col = raw$id_col
  )
}

load_te_expression_condition_means <- function(E_path,
                                               reference_te_subfamilies = NULL,
                                               duplicate_method = "max_sum",
                                               id_col = "gene_id",
                                               condition_groups = NULL,
                                               condition_order = c("PRO", "ESEN", "LSEN"),
                                               sample_cols = NULL,
                                               sample_pattern = "^(PRO|ESEN|LSEN)",
                                               te_type_col_candidates = c(
                                                 "gene_or_te",
                                                 "feature_type",
                                                 "feature",
                                                 "type",
                                                 "Source",
                                                 "source"
                                               ),
                                               te_type_values = c(
                                                 "repeat",
                                                 "repeats",
                                                 "te",
                                                 "tes",
                                                 "transposable_element",
                                                 "transposable_elements"
                                               ),
                                               verbose = TRUE) {
  if (!exists("resolve_duplicates", mode = "function") ||
      !exists("duplicate_stats", mode = "function")) {
    stop("load_te_expression_condition_means() requires resolve_duplicates() and duplicate_stats(). Source scripts/craTEs_own_function.R first.", call. = FALSE)
  }

  duplicate_method <- .most_active_normalize_duplicate_method(duplicate_method)
  if (!is.null(condition_groups)) {
    if (!is.list(condition_groups) || is.null(names(condition_groups))) {
      stop("condition_groups must be NULL or a named list of sample vectors.", call. = FALSE)
    }
    if (is.null(sample_cols)) {
      sample_cols <- unique(unlist(condition_groups, use.names = FALSE))
    }
  }

  raw <- read_raw_expression_count_matrix(
    E_path = E_path,
    id_col = id_col,
    sample_cols = sample_cols,
    sample_pattern = sample_pattern
  )

  if (is.null(condition_groups)) {
    condition_groups <- infer_senescence_condition_groups(
      raw$sample_cols,
      condition_prefixes = condition_order
    )
  }

  row_ids <- rownames(raw$counts)
  reference_te_subfamilies <- unique(as.character(reference_te_subfamilies))
  reference_te_subfamilies <- reference_te_subfamilies[!is.na(reference_te_subfamilies) & nzchar(reference_te_subfamilies)]
  parsed_subfamilies <- extract_te_subfamily_from_feature_id(
    row_ids,
    reference_te_subfamilies = reference_te_subfamilies
  )

  reference_keep <- if (length(reference_te_subfamilies) > 0L) {
    row_ids %in% reference_te_subfamilies | parsed_subfamilies %in% reference_te_subfamilies
  } else {
    rep(FALSE, length(row_ids))
  }

  type_col <- intersect(te_type_col_candidates, colnames(raw$feature_table))[1]
  type_keep <- rep(FALSE, length(row_ids))
  if (!is.na(type_col)) {
    type_values <- tolower(gsub("[^a-z0-9]+", "_", as.character(raw$feature_table[[type_col]])))
    type_targets <- tolower(gsub("[^a-z0-9]+", "_", te_type_values))
    type_keep <- type_values %in% type_targets |
      grepl("(^|_)te($|_)|transposable|repeat", type_values)
    keep <- type_keep
  } else {
    keep <- reference_keep
  }

  if (!any(keep)) {
    stop(
      paste(
        "No TE-expression rows were detected.",
        "Expected a TE marker column such as gene_or_te == 'repeat' or row IDs matching reference TE subfamilies.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  te_counts <- raw$counts[keep, , drop = FALSE]
  rownames(te_counts) <- parsed_subfamilies[keep]
  te_condition_means <- condition_mean_matrix(
    te_counts,
    condition_groups = condition_groups,
    condition_order = condition_order,
    object_name = "TE expression matrix"
  )

  duplicate_before <- duplicate_stats(rownames(te_condition_means), te_condition_means, include_table = FALSE)
  resolved_condition_means <- resolve_duplicates(
    te_condition_means,
    method = duplicate_method,
    verbose = verbose
  )
  duplicate_after <- duplicate_stats(rownames(resolved_condition_means), resolved_condition_means, include_table = FALSE)

  selection_table <- data.frame(
    feature_id = row_ids[keep],
    te_subfamily = parsed_subfamilies[keep],
    matched_reference = reference_keep[keep],
    matched_te_type_column = type_keep[keep],
    stringsAsFactors = FALSE
  )

  selection_summary <- data.frame(
    expression_rows_total = length(row_ids),
    expression_rows_kept_as_te = sum(keep),
    rows_kept_by_reference = sum(reference_keep),
    rows_kept_by_type_column = sum(type_keep),
    te_type_column = if (is.na(type_col)) NA_character_ else type_col,
    unique_te_subfamilies_before_duplicate_handling = length(unique(parsed_subfamilies[keep])),
    duplicate_groups_before = duplicate_before$duplicate_groups,
    redundant_rows_before = duplicate_before$redundant_rows,
    te_subfamilies_after_duplicate_handling = nrow(resolved_condition_means),
    duplicate_groups_after = duplicate_after$duplicate_groups,
    duplicate_method = duplicate_method,
    condition_mean_step = "before_duplicate_handling",
    stringsAsFactors = FALSE
  )

  list(
    condition_means = resolved_condition_means,
    condition_means_before_duplicate_handling = te_condition_means,
    counts_before_condition_means = te_counts,
    selection_table = selection_table,
    selection_summary = selection_summary,
    sample_cols = raw$sample_cols,
    id_col = raw$id_col,
    condition_groups = condition_groups
  )
}

activity_condition_mean_table <- function(activity_scores,
                                          condition_groups,
                                          condition_order = c("PRO", "ESEN", "LSEN")) {
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores)) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(activity_scores)) || any(!nzchar(rownames(activity_scores)))) {
    stop("activity_scores must have non-empty TE-subfamily rownames.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }
  if (!is.list(condition_groups) || is.null(names(condition_groups)) || any(!nzchar(names(condition_groups)))) {
    stop("condition_groups must be a named list of sample vectors.", call. = FALSE)
  }

  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[condition_order %in% names(condition_groups)]
  if (length(condition_order) == 0L) {
    stop("condition_order must contain at least one condition present in condition_groups.", call. = FALSE)
  }

  condition_groups <- condition_groups[condition_order]
  missing_samples <- setdiff(unique(unlist(condition_groups, use.names = FALSE)), colnames(activity_scores))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("Samples not found in activity_scores: %s", paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  condition_means <- do.call(cbind, lapply(condition_groups, function(samples) {
    rowMeans(activity_scores[, samples, drop = FALSE], na.rm = TRUE)
  }))
  colnames(condition_means) <- condition_order

  out <- data.frame(
    te_subfamily = rownames(activity_scores),
    stringsAsFactors = FALSE
  )
  for (condition in condition_order) {
    out[[paste0("mean_", condition)]] <- unname(condition_means[, condition])
    out[[paste0("abs_mean_", condition)]] <- abs(unname(condition_means[, condition]))
  }
  out$sum_abs_mean_activity_score <- rowSums(abs(condition_means), na.rm = TRUE)
  out
}

rank_activity_te_subfamilies <- function(activity_scores,
                                         condition_groups,
                                         k_values = c(3, 5, 10, 20, 50, 100),
                                         condition_order = c("PRO", "ESEN", "LSEN"),
                                         ranking_metric = "sum_abs_mean_activity_score",
                                         scope = "global") {
  k_values <- .most_active_validate_k_values(k_values)
  scope <- as.character(scope)
  if (length(scope) != 1L || is.na(scope) || !nzchar(scope)) {
    stop("scope must be one non-empty value.", call. = FALSE)
  }

  if (identical(ranking_metric, "sum_abs_mean_activity_score")) {
    ranked <- activity_condition_mean_table(
      activity_scores = activity_scores,
      condition_groups = condition_groups,
      condition_order = condition_order
    )
    if (identical(scope, "global")) {
      ranked$ranking_metric <- ranking_metric
      ranked$ranking_score <- ranked$sum_abs_mean_activity_score
      ranked$signed_activity_score <- NA_real_
    } else {
      condition_order <- unique(as.character(condition_order))
      if (!scope %in% condition_order) {
        stop(sprintf("scope must be 'global' or one of: %s", paste(condition_order, collapse = ", ")), call. = FALSE)
      }
      score_col <- paste0("abs_mean_", scope)
      signed_col <- paste0("mean_", scope)
      ranked$ranking_metric <- paste0("abs_mean_activity_score_", scope)
      ranked$ranking_score <- ranked[[score_col]]
      ranked$signed_activity_score <- ranked[[signed_col]]
    }
    ranked$scope <- scope
    ranked <- ranked[order(-ranked$ranking_score, ranked$te_subfamily), , drop = FALSE]
    ranked <- utils::head(ranked, min(max(k_values), nrow(ranked)))
    ranked$rank <- seq_len(nrow(ranked))
  } else {
    if (!identical(scope, "global")) {
      stop("Condition-specific activity ranking is supported for ranking_metric = 'sum_abs_mean_activity_score' only.", call. = FALSE)
    }
    if (!exists("select_top_te_subfamilies_by_activity_change", mode = "function")) {
      stop("rank_activity_te_subfamilies() requires select_top_te_subfamilies_by_activity_change(). Source scripts/craTEs_own_function.R first.", call. = FALSE)
    }
    ranked <- select_top_te_subfamilies_by_activity_change(
      activity_scores = activity_scores,
      condition_groups = condition_groups,
      condition_order = condition_order,
      top_k = max(k_values),
      ranking_metric = ranking_metric,
      te_subfamilies = NULL
    )
    ranked$scope <- scope
    ranked$signed_activity_score <- NA_real_
    ranked$rank <- seq_len(nrow(ranked))
  }

  front_cols <- c("scope", "rank", "te_subfamily", "ranking_metric", "ranking_score", "signed_activity_score", "sum_abs_mean_activity_score")
  front_cols <- intersect(front_cols, colnames(ranked))
  ranked <- ranked[, c(front_cols, setdiff(colnames(ranked), front_cols)), drop = FALSE]
  rownames(ranked) <- NULL
  ranked
}

rank_raw_te_expression <- function(raw_te_counts,
                                   k_values = c(3, 5, 10, 20, 50, 100)) {
  k_values <- .most_active_validate_k_values(k_values)
  if (!is.matrix(raw_te_counts) && !is.data.frame(raw_te_counts)) {
    stop("raw_te_counts must be a matrix or data.frame.", call. = FALSE)
  }

  raw_te_counts <- as.matrix(raw_te_counts)
  if (!is.numeric(raw_te_counts)) {
    stop("raw_te_counts must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(raw_te_counts)) || any(!nzchar(rownames(raw_te_counts)))) {
    stop("raw_te_counts must have non-empty TE-subfamily rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(raw_te_counts))) {
    stop("raw_te_counts must be duplicate-resolved before ranking.", call. = FALSE)
  }

  raw_expression_sum <- rowSums(raw_te_counts, na.rm = TRUE)
  raw_expression_mean <- rowMeans(raw_te_counts, na.rm = TRUE)
  ranked <- data.frame(
    te_subfamily = rownames(raw_te_counts),
    raw_expression_sum = unname(raw_expression_sum),
    raw_expression_mean = unname(raw_expression_mean),
    n_samples = ncol(raw_te_counts),
    stringsAsFactors = FALSE
  )
  ranked <- ranked[order(-ranked$raw_expression_sum, ranked$te_subfamily), , drop = FALSE]
  ranked <- utils::head(ranked, min(max(k_values), nrow(ranked)))
  ranked$rank <- seq_len(nrow(ranked))

  front_cols <- c("rank", "te_subfamily", "raw_expression_sum", "raw_expression_mean", "n_samples")
  ranked <- ranked[, front_cols, drop = FALSE]
  rownames(ranked) <- NULL
  ranked
}

rank_te_expression_condition_means <- function(te_expression_condition_means,
                                               k_values = c(3, 5, 10, 20, 50, 100),
                                               condition_order = c("PRO", "ESEN", "LSEN"),
                                               scope = "global") {
  k_values <- .most_active_validate_k_values(k_values)
  scope <- as.character(scope)
  if (length(scope) != 1L || is.na(scope) || !nzchar(scope)) {
    stop("scope must be one non-empty value.", call. = FALSE)
  }
  if (!is.matrix(te_expression_condition_means) && !is.data.frame(te_expression_condition_means)) {
    stop("te_expression_condition_means must be a matrix or data.frame.", call. = FALSE)
  }

  condition_means <- as.matrix(te_expression_condition_means)
  if (!is.numeric(condition_means)) {
    stop("te_expression_condition_means must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(condition_means)) || any(!nzchar(rownames(condition_means)))) {
    stop("te_expression_condition_means must have non-empty TE-subfamily rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(condition_means))) {
    stop("te_expression_condition_means must be duplicate-resolved before ranking.", call. = FALSE)
  }

  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[condition_order %in% colnames(condition_means)]
  if (length(condition_order) == 0L) {
    stop("condition_order must contain at least one column present in te_expression_condition_means.", call. = FALSE)
  }

  ranked <- data.frame(
    te_subfamily = rownames(condition_means),
    stringsAsFactors = FALSE
  )
  for (condition in condition_order) {
    ranked[[paste0("mean_expression_", condition)]] <- unname(condition_means[, condition])
  }
  ranked$sum_condition_mean_expression_score <- rowSums(condition_means[, condition_order, drop = FALSE], na.rm = TRUE)
  ranked$average_condition_mean_expression <- rowMeans(condition_means[, condition_order, drop = FALSE], na.rm = TRUE)

  if (identical(scope, "global")) {
    ranked$ranking_metric <- "sum_condition_mean_expression_score"
    ranked$ranking_score <- ranked$sum_condition_mean_expression_score
  } else {
    if (!scope %in% condition_order) {
      stop(sprintf("scope must be 'global' or one of: %s", paste(condition_order, collapse = ", ")), call. = FALSE)
    }
    score_col <- paste0("mean_expression_", scope)
    ranked$ranking_metric <- paste0("mean_expression_", scope)
    ranked$ranking_score <- ranked[[score_col]]
  }
  ranked$scope <- scope

  ranked <- ranked[order(-ranked$ranking_score, ranked$te_subfamily), , drop = FALSE]
  ranked <- utils::head(ranked, min(max(k_values), nrow(ranked)))
  ranked$rank <- seq_len(nrow(ranked))

  front_cols <- c(
    "scope",
    "rank",
    "te_subfamily",
    "ranking_metric",
    "ranking_score",
    "sum_condition_mean_expression_score",
    "average_condition_mean_expression"
  )
  front_cols <- intersect(front_cols, colnames(ranked))
  ranked <- ranked[, c(front_cols, setdiff(colnames(ranked), front_cols)), drop = FALSE]
  rownames(ranked) <- NULL
  ranked
}

top_k_long_table <- function(ranked_table,
                             k_values = c(3, 5, 10, 20, 50, 100),
                             source = "ranked",
                             score_col = "ranking_score") {
  k_values <- .most_active_validate_k_values(k_values)
  if (!is.data.frame(ranked_table) || !"te_subfamily" %in% colnames(ranked_table)) {
    stop("ranked_table must be a data.frame with column 'te_subfamily'.", call. = FALSE)
  }
  if (!score_col %in% colnames(ranked_table)) {
    stop(sprintf("score_col '%s' was not found in ranked_table.", score_col), call. = FALSE)
  }

  do.call(rbind, lapply(k_values, function(k) {
    selected <- utils::head(ranked_table, min(k, nrow(ranked_table)))
    selected$k <- k
    selected$source <- source
    selected$rank_within_k <- seq_len(nrow(selected))
    front_cols <- c("source", "scope", "k", "rank_within_k", "rank", "te_subfamily", "ranking_metric", score_col)
    front_cols <- intersect(front_cols, colnames(selected))
    selected[, c(front_cols, setdiff(colnames(selected), front_cols)), drop = FALSE]
  }))
}

build_top_k_intersections <- function(activity_ranked,
                                      raw_ranked,
                                      k_values = c(3, 5, 10, 20, 50, 100),
                                      activity_score_col = "ranking_score",
                                      raw_score_col = "raw_expression_sum",
                                      scope = "global") {
  k_values <- .most_active_validate_k_values(k_values)
  scope <- as.character(scope)
  if (length(scope) != 1L || is.na(scope) || !nzchar(scope)) {
    stop("scope must be one non-empty value.", call. = FALSE)
  }
  required_activity_cols <- c("rank", "te_subfamily", activity_score_col)
  required_raw_cols <- c("rank", "te_subfamily", raw_score_col)
  missing_activity_cols <- setdiff(required_activity_cols, colnames(activity_ranked))
  missing_raw_cols <- setdiff(required_raw_cols, colnames(raw_ranked))
  if (length(missing_activity_cols) > 0L) {
    stop(sprintf("activity_ranked is missing columns: %s", paste(missing_activity_cols, collapse = ", ")), call. = FALSE)
  }
  if (length(missing_raw_cols) > 0L) {
    stop(sprintf("raw_ranked is missing columns: %s", paste(missing_raw_cols, collapse = ", ")), call. = FALSE)
  }

  summary_rows <- vector("list", length(k_values))
  member_rows <- list()

  for (i in seq_along(k_values)) {
    k <- k_values[[i]]
    activity_top <- utils::head(activity_ranked$te_subfamily, min(k, nrow(activity_ranked)))
    raw_top <- utils::head(raw_ranked$te_subfamily, min(k, nrow(raw_ranked)))
    shared <- intersect(activity_top, raw_top)
    joined <- union(activity_top, raw_top)

    summary_rows[[i]] <- data.frame(
      scope = scope,
      k = k,
      n_activity_top = length(activity_top),
      n_te_expression_top = length(raw_top),
      n_raw_top = length(raw_top),
      n_intersection = length(shared),
      pct_activity_top = if (length(activity_top) == 0L) NA_real_ else 100 * length(shared) / length(activity_top),
      pct_te_expression_top = if (length(raw_top) == 0L) NA_real_ else 100 * length(shared) / length(raw_top),
      pct_raw_top = if (length(raw_top) == 0L) NA_real_ else 100 * length(shared) / length(raw_top),
      jaccard = if (length(joined) == 0L) NA_real_ else length(shared) / length(joined),
      intersect_te_subfamilies = paste(shared, collapse = ", "),
      stringsAsFactors = FALSE
    )

    if (length(shared) > 0L) {
      activity_idx <- match(shared, activity_ranked$te_subfamily)
      raw_idx <- match(shared, raw_ranked$te_subfamily)
      member_rows[[length(member_rows) + 1L]] <- data.frame(
        scope = scope,
        k = k,
        te_subfamily = shared,
        activity_rank = activity_ranked$rank[activity_idx],
        activity_score = activity_ranked[[activity_score_col]][activity_idx],
        te_expression_rank = raw_ranked$rank[raw_idx],
        te_expression_score = raw_ranked[[raw_score_col]][raw_idx],
        raw_rank = raw_ranked$rank[raw_idx],
        raw_expression_sum = raw_ranked[[raw_score_col]][raw_idx],
        combined_rank_sum = activity_ranked$rank[activity_idx] + raw_ranked$rank[raw_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  members <- if (length(member_rows) > 0L) {
    do.call(rbind, member_rows)
  } else {
    data.frame(
      scope = character(),
      k = integer(),
      te_subfamily = character(),
      activity_rank = integer(),
      activity_score = numeric(),
      te_expression_rank = integer(),
      te_expression_score = numeric(),
      raw_rank = integer(),
      raw_expression_sum = numeric(),
      combined_rank_sum = integer(),
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, summary_rows),
    members = members
  )
}

build_global_intersection_top_n <- function(activity_ranked,
                                            raw_ranked,
                                            top_n = 10L,
                                            activity_score_col = "ranking_score",
                                            raw_score_col = "raw_expression_sum",
                                            scope = "global") {
  if (!is.data.frame(activity_ranked) || !is.data.frame(raw_ranked)) {
    stop("activity_ranked and raw_ranked must be data.frames.", call. = FALSE)
  }
  scope <- as.character(scope)
  if (length(scope) != 1L || is.na(scope) || !nzchar(scope)) {
    stop("scope must be one non-empty value.", call. = FALSE)
  }
  required_activity_cols <- c("rank", "te_subfamily", activity_score_col)
  required_raw_cols <- c("rank", "te_subfamily", raw_score_col)
  missing_activity_cols <- setdiff(required_activity_cols, colnames(activity_ranked))
  missing_raw_cols <- setdiff(required_raw_cols, colnames(raw_ranked))
  if (length(missing_activity_cols) > 0L) {
    stop(sprintf("activity_ranked is missing columns: %s", paste(missing_activity_cols, collapse = ", ")), call. = FALSE)
  }
  if (length(missing_raw_cols) > 0L) {
    stop(sprintf("raw_ranked is missing columns: %s", paste(missing_raw_cols, collapse = ", ")), call. = FALSE)
  }

  shared <- intersect(activity_ranked$te_subfamily, raw_ranked$te_subfamily)
  if (length(shared) == 0L) {
    return(data.frame(
      global_intersection_rank = integer(),
      scope = character(),
      te_subfamily = character(),
      activity_rank = integer(),
      activity_score = numeric(),
      te_expression_rank = integer(),
      te_expression_score = numeric(),
      raw_rank = integer(),
      raw_expression_sum = numeric(),
      combined_rank_sum = integer(),
      stringsAsFactors = FALSE
    ))
  }

  activity_idx <- match(shared, activity_ranked$te_subfamily)
  raw_idx <- match(shared, raw_ranked$te_subfamily)
  out <- data.frame(
    scope = scope,
    te_subfamily = shared,
    activity_rank = activity_ranked$rank[activity_idx],
    activity_score = activity_ranked[[activity_score_col]][activity_idx],
    te_expression_rank = raw_ranked$rank[raw_idx],
    te_expression_score = raw_ranked[[raw_score_col]][raw_idx],
    raw_rank = raw_ranked$rank[raw_idx],
    raw_expression_sum = raw_ranked[[raw_score_col]][raw_idx],
    stringsAsFactors = FALSE
  )
  out$combined_rank_sum <- out$activity_rank + out$raw_rank
  out <- out[order(out$combined_rank_sum, out$activity_rank, out$raw_rank, out$te_subfamily), , drop = FALSE]
  out <- utils::head(out, min(as.integer(top_n), nrow(out)))
  out$global_intersection_rank <- seq_len(nrow(out))
  out <- out[, c("global_intersection_rank", setdiff(colnames(out), "global_intersection_rank")), drop = FALSE]
  rownames(out) <- NULL
  out
}
