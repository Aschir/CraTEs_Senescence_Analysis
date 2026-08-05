# Vinzent Aschir 
# 28.04.2026 

print("sourcing the custom function from crates")

if (!exists("classify_te_senescence_trend", mode = "function")) {
  source("scripts/te_trend_classification_f.R")
}

.restore_integer_counts <- function(mat, tolerance = sqrt(.Machine$double.eps)) {
  if (!is.numeric(mat)) {
    return(mat)
  }

  finite_values <- mat[is.finite(mat)]
  is_integer_like <- length(finite_values) == 0L ||
    all(abs(finite_values - round(finite_values)) <= tolerance)
  fits_integer <- length(finite_values) == 0L ||
    all(finite_values >= -.Machine$integer.max & finite_values <= .Machine$integer.max)

  if (is_integer_like && fits_integer) {
    storage.mode(mat) <- "integer"
  }

  mat
}

check_non_integer_counts <- function(count_table,
                                     tolerance = sqrt(.Machine$double.eps),
                                     max_entries_to_print = 200,
                                     print_entries = TRUE) {
  # Audits a count matrix for fractional values while keeping the full table available.
  if (!is.matrix(count_table) && !is.data.frame(count_table)) {
    stop("count_table must be a matrix or data.frame", call. = FALSE)
  }

  count_matrix <- as.matrix(count_table)
  source_cols <- intersect(c("Source", "source", "mapping_source"), colnames(count_matrix))
  if (length(source_cols) > 0L) {
    count_matrix <- count_matrix[, setdiff(colnames(count_matrix), source_cols), drop = FALSE]
  }

  numeric_matrix <- matrix(
    suppressWarnings(as.numeric(count_matrix)),
    nrow = nrow(count_matrix),
    dimnames = dimnames(count_matrix)
  )

  conversion_failed <- is.na(numeric_matrix) & !is.na(count_matrix) & count_matrix != ""
  if (any(conversion_failed)) {
    stop("count_table contains non-numeric entries after removing source columns", call. = FALSE)
  }

  if (is.null(rownames(numeric_matrix))) {
    rownames(numeric_matrix) <- paste0("row_", seq_len(nrow(numeric_matrix)))
  }
  if (is.null(colnames(numeric_matrix))) {
    colnames(numeric_matrix) <- paste0("col_", seq_len(ncol(numeric_matrix)))
  }

  non_integer_mask <- is.finite(numeric_matrix) &
    abs(numeric_matrix - round(numeric_matrix)) > tolerance
  non_integer_idx <- which(non_integer_mask, arr.ind = TRUE)

  if (nrow(non_integer_idx) == 0L) {
    non_integer_entries <- data.frame(
      gene = character(),
      sample = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    non_integer_entries <- data.frame(
      gene = rownames(numeric_matrix)[non_integer_idx[, "row"]],
      sample = colnames(numeric_matrix)[non_integer_idx[, "col"]],
      value = numeric_matrix[non_integer_idx],
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf(
    "Non-integer entries: %d / %d (%.4f%%)\n",
    nrow(non_integer_entries),
    length(numeric_matrix),
    100 * nrow(non_integer_entries) / length(numeric_matrix)
  ))
  cat(sprintf("Affected genes: %d\n", length(unique(non_integer_entries$gene))))
  cat(sprintf("Affected samples: %d\n", length(unique(non_integer_entries$sample))))

  if (isTRUE(print_entries)) {
    if (nrow(non_integer_entries) == 0L) {
      cat("All count entries are integer-like within tolerance.\n")
    } else if (nrow(non_integer_entries) <= max_entries_to_print) {
      print(non_integer_entries)
    } else {
      cat(sprintf(
        "Showing first %d entries. Full table is returned by check_non_integer_counts().\n",
        max_entries_to_print
      ))
      print(utils::head(non_integer_entries, max_entries_to_print))
    }
  }

  non_integer_entries
}

as_count_association_matrix <- function(N,
                                        id_col = 1,
                                        object_name = "N") {
  if (is.matrix(N) && is.numeric(N)) {
    return(N)
  }

  if (!is.matrix(N) && !is.data.frame(N)) {
    stop(sprintf("%s must be a matrix or data.frame/data.table.", object_name), call. = FALSE)
  }

  if (ncol(N) < 1L) {
    stop(sprintf("%s must contain at least one column.", object_name), call. = FALSE)
  }

  if (is.data.frame(N)) {
    N_table <- N
    if (ncol(N_table) < 2L) {
      stop(
        sprintf("%s table must contain one gene ID column and at least one TE-subfamily column.", object_name),
        call. = FALSE
      )
    }

    if (is.character(id_col)) {
      if (!id_col %in% colnames(N_table)) {
        stop(sprintf("id_col '%s' was not found in %s.", id_col, object_name), call. = FALSE)
      }
      id_index <- match(id_col, colnames(N_table))
    } else {
      id_index <- as.integer(id_col)
      if (length(id_index) != 1L || is.na(id_index) || id_index < 1L || id_index > ncol(N_table)) {
        stop("id_col must be a valid column name or one-based column index.", call. = FALSE)
      }
    }

    genes_N <- as.character(N_table[[id_index]])
    if (any(is.na(genes_N) | genes_N == "")) {
      stop(sprintf("%s gene ID column contains missing or empty identifiers.", object_name), call. = FALSE)
    }

    value_cols <- setdiff(seq_len(ncol(N_table)), id_index)
    if (inherits(N_table, "data.table")) {
      value_col_names <- colnames(N_table)[value_cols]
      value_raw <- as.matrix(N_table[, ..value_col_names])
    } else {
      value_raw <- as.matrix(N_table[, value_cols, drop = FALSE])
    }
  } else {
    value_raw <- N
    genes_N <- rownames(N)

    if (!is.numeric(N) && ncol(N) >= 2L) {
      first_col <- N[, 1]
      first_col_numeric <- suppressWarnings(as.numeric(first_col))
      first_col_has_text <- any(is.na(first_col_numeric) & !is.na(first_col) & first_col != "")

      if (isTRUE(first_col_has_text)) {
        genes_N <- as.character(first_col)
        value_raw <- N[, -1, drop = FALSE]
      }
    }

    if (is.null(genes_N)) {
      genes_N <- paste0("row_", seq_len(nrow(value_raw)))
    }
  }

  value_numeric <- suppressWarnings(as.numeric(value_raw))
  N_numeric <- matrix(
    value_numeric,
    nrow = nrow(value_raw),
    dimnames = list(genes_N, colnames(value_raw))
  )

  conversion_failed <- is.na(N_numeric) & !is.na(value_raw) & value_raw != ""
  if (any(conversion_failed)) {
    failed_idx <- which(conversion_failed, arr.ind = TRUE)
    examples <- utils::head(data.frame(
      gene = rownames(N_numeric)[failed_idx[, "row"]],
      column = colnames(value_raw)[failed_idx[, "col"]],
      value = value_raw[failed_idx],
      stringsAsFactors = FALSE
    ), 10)
    stop(
      paste(
        sprintf("%s contains non-numeric TE-subfamily values after removing the gene ID column.", object_name),
        "Examples:",
        paste(utils::capture.output(print(examples, row.names = FALSE)), collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  N_numeric
}

load_count_association_matrix <- function(N_path,
                                          id_col = 1,
                                          sep = "\t",
                                          header = TRUE) {
  if (is.null(N_path) || !nzchar(N_path)) {
    stop("N_path must be a non-empty file path", call. = FALSE)
  }
  if (!file.exists(N_path)) {
    stop(sprintf("N_path does not exist: %s", N_path), call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("load_count_association_matrix() requires package 'data.table'.", call. = FALSE)
  }

  N_table <- data.table::fread(N_path, header = header, sep = sep, check.names = FALSE)
  if (ncol(N_table) < 2L) {
    stop("N table must contain one gene ID column and at least one TE-subfamily column.", call. = FALSE)
  }

  N <- as_count_association_matrix(N_table, id_col = id_col, object_name = "N table")

  list(
    N = N,
    genes_N = rownames(N)
  )
}

normalize_reference_genome <- function(reference_genome) {
  reference_genome <- as.character(reference_genome)
  if (length(reference_genome) != 1L || is.na(reference_genome) || !nzchar(reference_genome)) {
    stop("reference_genome must be one of 'Hg19' or 'Hg38'.", call. = FALSE)
  }

  reference_key <- tolower(reference_genome)
  if (reference_key %in% c("hg19", "grch37")) {
    return("Hg19")
  }
  if (reference_key %in% c("hg38", "grch38")) {
    return("Hg38")
  }

  stop("reference_genome must be one of 'Hg19' or 'Hg38'.", call. = FALSE)
}

infer_reference_genome_from_path <- function(path) {
  path_key <- tolower(as.character(path))
  if (grepl("hg19|grch37", path_key)) {
    return("Hg19")
  }
  if (grepl("hg38|grch38", path_key)) {
    return("Hg38")
  }

  stop(sprintf("Could not infer reference genome from path: %s", path), call. = FALSE)
}

infer_analysis_mode_from_n_path <- function(N_path) {
  n_base <- basename(as.character(N_path))
  n_key <- tolower(n_base)

  if (grepl("^n_c_|_n_c_|count", n_key)) {
    return("counted")
  }
  if (grepl("^n_w_|_n_w_|weight", n_key)) {
    return("weighted")
  }

  stop(sprintf("Could not infer analysis mode from N matrix path: %s", N_path), call. = FALSE)
}

validate_reference_n_settings <- function(reference_genome,
                                          N_path,
                                          analysis_mode = NULL) {
  selected_reference <- normalize_reference_genome(reference_genome)
  inferred_reference <- infer_reference_genome_from_path(N_path)
  if (!identical(selected_reference, inferred_reference)) {
    stop(
      sprintf(
        "reference_genome is '%s', but N_path appears to be '%s': %s",
        selected_reference,
        inferred_reference,
        N_path
      ),
      call. = FALSE
    )
  }

  inferred_mode <- infer_analysis_mode_from_n_path(N_path)
  if (!is.null(analysis_mode)) {
    requested_mode <- tolower(as.character(analysis_mode))
    mode_alias <- c(
      counted = "counted",
      count = "counted",
      counts = "counted",
      c = "counted",
      weighted = "weighted",
      weight = "weighted",
      weights = "weighted",
      w = "weighted"
    )
    if (!requested_mode %in% names(mode_alias)) {
      stop("analysis_mode must be NULL, 'counted', or 'weighted'.", call. = FALSE)
    }
    requested_mode <- unname(mode_alias[[requested_mode]])
    if (!identical(requested_mode, inferred_mode)) {
      stop(
        sprintf(
          "analysis_mode is '%s', but N_path implies '%s': %s",
          requested_mode,
          inferred_mode,
          N_path
        ),
        call. = FALSE
      )
    }
  }

  list(
    reference_genome = selected_reference,
    org = tolower(selected_reference),
    inferred_reference_genome = inferred_reference,
    analysis_mode = inferred_mode
  )
}

.result_path_stem <- function(path) {
  stem <- basename(path)
  stem <- sub("\\.gz$", "", stem)
  stem <- sub("\\.[^.]+$", "", stem)
  stem <- sub("_gene_names$", "", stem)
  stem <- sub("_genes$", "", stem)
  stem
}

.safe_result_name <- function(name) {
  name <- as.character(name)
  name <- gsub("[^A-Za-z0-9._-]+", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  if (!nzchar(name)) {
    stop("name must contain at least one file-name-safe character", call. = FALSE)
  }
  name
}

.has_meaningful_rownames <- function(x) {
  rn <- rownames(x)
  !is.null(rn) && !identical(rn, as.character(seq_len(nrow(x))))
}

.result_object_to_table <- function(object) {
  if (is.matrix(object)) {
    df <- as.data.frame(object, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(rownames(object))) {
      df <- cbind(row_id = rownames(object), df)
    }
    return(df)
  }

  if (is.data.frame(object)) {
    df <- as.data.frame(object, check.names = FALSE, stringsAsFactors = FALSE)
    if (.has_meaningful_rownames(df)) {
      df <- cbind(row_id = rownames(df), df)
    }
    return(df)
  }

  if (is.atomic(object) && is.null(dim(object))) {
    value_names <- names(object)
    if (is.null(value_names)) {
      value_names <- seq_along(object)
    }
    return(data.frame(
      name = value_names,
      value = as.vector(object),
      stringsAsFactors = FALSE
    ))
  }

  NULL
}

derive_activity_results_dir <- function(E_path,
                                        N_path,
                                        results_parent = NULL,
                                        prefix = "results") {
  if (is.null(E_path) || !nzchar(E_path)) {
    stop("E_path is not defined.", call. = FALSE)
  }
  if (is.null(N_path) || !nzchar(N_path)) {
    stop("N_path is not defined.", call. = FALSE)
  }

  E_path <- normalizePath(path.expand(E_path), mustWork = TRUE)
  N_path <- normalizePath(path.expand(N_path), mustWork = TRUE)

  if (is.null(results_parent)) {
    # Default to the current project root results directory rather than a subfolder
    # under the input data location.
    results_parent <- file.path(getwd(), "results")
  }
  results_parent <- normalizePath(path.expand(results_parent), mustWork = FALSE)

  e_base <- .result_path_stem(E_path)
  n_base <- .result_path_stem(N_path)
  result_id <- paste0(prefix, "_", e_base, "_", n_base)

  file.path(results_parent, result_id)
}

# Write one result object in a VS Code friendly TSV representation by default.
# Set write_rds = TRUE when the exact R object class, nested list structure, or
# attributes should be preserved next to the flat export.
write_result_object <- function(input_dir,
                                name,
                                object,
                                write_rds = FALSE,
                                write_tsv = TRUE,
                                overwrite = TRUE) {
  if (is.null(input_dir) || !nzchar(input_dir)) {
    stop("input_dir must be a non-empty directory path", call. = FALSE)
  }

  input_dir <- path.expand(input_dir)
  if (!dir.exists(input_dir)) {
    dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  }

  safe_name <- .safe_result_name(name)
  written <- data.frame(
    name = character(),
    format = character(),
    path = character(),
    stringsAsFactors = FALSE
  )

  write_one_table <- function(table_object, file_name, result_name) {
    path <- file.path(input_dir, paste0(file_name, ".tsv"))
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop(sprintf("Refusing to overwrite existing file: %s", path), call. = FALSE)
    }
    utils::write.table(
      table_object,
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    data.frame(
      name = result_name,
      format = "tsv",
      path = normalizePath(path, mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(write_rds)) {
    rds_path <- file.path(input_dir, paste0(safe_name, ".rds"))
    if (file.exists(rds_path) && !isTRUE(overwrite)) {
      stop(sprintf("Refusing to overwrite existing file: %s", rds_path), call. = FALSE)
    }
    saveRDS(object, rds_path)
    written <- rbind(
      written,
      data.frame(
        name = safe_name,
        format = "rds",
        path = normalizePath(rds_path, mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }

  if (isTRUE(write_tsv)) {
    table_object <- .result_object_to_table(object)
    if (!is.null(table_object)) {
      written <- rbind(written, write_one_table(table_object, safe_name, safe_name))
    } else if (is.list(object)) {
      for (element_name in names(object)) {
        element <- object[[element_name]]
        element_table <- .result_object_to_table(element)
        if (is.null(element_table)) {
          next
        }
        safe_element_name <- .safe_result_name(paste(safe_name, element_name, sep = "_"))
        written <- rbind(
          written,
          write_one_table(element_table, safe_element_name, safe_element_name)
        )
      }
    }
  }

  written
}

# Save a TE activity score matrix as a TSV while preserving row and column names.
# Row names are written into the first column so the file stays easy to inspect in
# VS Code and can be reconstructed with read.delim(..., check.names = FALSE).
save_activity_score_matrix_tsv <- function(activity_scores,
                                           output_path,
                                           row_id_col = "te_subfamily",
                                           overwrite = TRUE) {
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  if (is.null(output_path) || !nzchar(output_path)) {
    stop("output_path must be a non-empty file path.", call. = FALSE)
  }
  if (is.null(row_id_col) || length(row_id_col) != 1L || !nzchar(row_id_col)) {
    stop("row_id_col must be one non-empty column name.", call. = FALSE)
  }

  row_ids <- rownames(activity_scores)
  if (is.null(row_ids) || any(is.na(row_ids)) || any(row_ids == "")) {
    stop("activity_scores must have non-empty rownames before saving.", call. = FALSE)
  }

  column_ids <- colnames(activity_scores)
  if (is.null(column_ids) || any(is.na(column_ids)) || any(column_ids == "")) {
    stop("activity_scores must have non-empty colnames before saving.", call. = FALSE)
  }
  if (row_id_col %in% column_ids) {
    stop(sprintf("row_id_col '%s' already exists in activity_scores colnames.", row_id_col), call. = FALSE)
  }

  output_path <- path.expand(output_path)
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  if (file.exists(output_path) && !isTRUE(overwrite)) {
    stop(sprintf("Refusing to overwrite existing file: %s", output_path), call. = FALSE)
  }

  activity_table <- as.data.frame(as.matrix(activity_scores), check.names = FALSE, stringsAsFactors = FALSE)
  activity_table <- cbind(
    stats::setNames(data.frame(row_ids, stringsAsFactors = FALSE), row_id_col),
    activity_table
  )

  utils::write.table(
    activity_table,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  normalizePath(output_path, mustWork = TRUE)
}

# Calculate and plot treatment-control delta activity for one condition.
# The result contains the union of top_k TE subfamilies by absolute delta and
# explicitly requested te_subfamilies.
plot_delta_activity <- function(activity_scores,
                                treatment_samples,
                                control_samples,
                                condition_label = "condition",
                                top_k = 10,
                                te_subfamilies = NULL,
                                plot = TRUE,
                                print_table = TRUE,
                                title = "Delta activity",
                                treatment_label = "treatment",
                                control_label = "control") {
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(activity_scores)) || any(!nzchar(rownames(activity_scores)))) {
    stop("activity_scores must have non-empty TE subfamily rownames.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }

  validate_sample_group <- function(samples, arg_name) {
    if (!is.character(samples)) {
      stop(sprintf("%s must be a character vector.", arg_name), call. = FALSE)
    }
    samples <- unique(samples[nzchar(samples)])
    if (length(samples) == 0L) {
      stop(sprintf("%s must contain at least one sample.", arg_name), call. = FALSE)
    }
    samples
  }

  treatment_samples <- validate_sample_group(treatment_samples, "treatment_samples")
  control_samples <- validate_sample_group(control_samples, "control_samples")

  if (is.null(condition_label) || length(condition_label) != 1L || !nzchar(condition_label)) {
    stop("condition_label must be one non-empty label.", call. = FALSE)
  }

  all_group_samples <- unique(c(treatment_samples, control_samples))
  missing_samples <- setdiff(all_group_samples, colnames(activity_scores))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("Samples not found in activity_scores colnames: %s", paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  if (is.null(top_k)) {
    top_k <- 0L
  }
  if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 0L) {
    stop("top_k must be a single non-negative number or NULL.", call. = FALSE)
  }
  top_k <- as.integer(top_k)

  if (!is.null(te_subfamilies)) {
    te_subfamilies <- unique(as.character(te_subfamilies))
    te_subfamilies <- te_subfamilies[nzchar(te_subfamilies)]
  }
  if (top_k == 0L && length(te_subfamilies) == 0L) {
    stop("Set top_k > 0 and/or provide te_subfamilies.", call. = FALSE)
  }

  treatment_mean <- rowMeans(activity_scores[, treatment_samples, drop = FALSE], na.rm = TRUE)
  control_mean <- rowMeans(activity_scores[, control_samples, drop = FALSE], na.rm = TRUE)
  delta <- treatment_mean - control_mean
  delta_all <- data.frame(
    te_subfamily = rownames(activity_scores),
    mean_treatment = unname(treatment_mean),
    mean_control = unname(control_mean),
    delta_activity = unname(delta),
    stringsAsFactors = FALSE
  )

  if (isTRUE(print_table)) {
    cat(sprintf("\n%s\n", condition_label))
    cat(sprintf(
      "Samples: %d treatment, %d control\n\n",
      length(treatment_samples),
      length(control_samples)
    ))
  }

  top_rows <- rep(FALSE, nrow(delta_all))
  if (top_k > 0L) {
    top_rows <- order(abs(delta_all$delta_activity), decreasing = TRUE)[seq_len(min(top_k, nrow(delta_all)))]
  }

  selected_by_top_k <- rep(FALSE, nrow(delta_all))
  selected_by_top_k[top_rows] <- TRUE
  selected_by_request <- delta_all$te_subfamily %in% te_subfamilies

  missing_tes <- setdiff(te_subfamilies, rownames(activity_scores))
  if (length(missing_tes) > 0L) {
    warning(
      sprintf("Requested TE subfamilies not found in activity_scores rownames: %s", paste(missing_tes, collapse = ", ")),
      call. = FALSE
    )
  }

  delta_table <- delta_all[selected_by_top_k | selected_by_request, , drop = FALSE]
  if (nrow(delta_table) == 0L) {
    stop("No TE subfamilies selected. Check top_k and te_subfamilies.", call. = FALSE)
  }
  delta_table <- delta_table[order(-abs(delta_table$delta_activity), delta_table$te_subfamily), ]
  rownames(delta_table) <- NULL

  if (isTRUE(print_table)) {
    print(delta_table, row.names = FALSE)
  }

  if (isTRUE(plot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("plot = TRUE requires package 'ggplot2'.", call. = FALSE)
    }
    plot_table <- delta_table
    te_order <- unique(plot_table$te_subfamily[order(abs(plot_table$delta_activity), decreasing = FALSE)])
    plot_table$te_subfamily <- factor(plot_table$te_subfamily, levels = te_order)
    p <- ggplot2::ggplot(
      plot_table,
      ggplot2::aes(
        x = te_subfamily,
        y = delta_activity
      )
    ) +
      ggplot2::geom_col(width = 0.7, fill = "#2F6F73") +
      ggplot2::geom_hline(yintercept = 0, size = 0.35, colour = "grey35") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = title,
        x = "TE subfamily",
        y = paste0("Mean activity delta (", treatment_label, " - ", control_label, ")")
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "right"
      )
    print(p)
  }

  delta_table
}

# Plot condition means for a TE set selected by pairwise delta activity.
plot_delta_activity_condition_ranges <- function(activity_scores,
                                                 delta_table,
                                                 condition_groups,
                                                 condition_label = "condition contrast",
                                                 condition_order = names(condition_groups),
                                                 condition_cols = NULL,
                                                 title = NULL,
                                                 x_label = "Mean inferred TE activity",
                                                 selection_label = NULL,
                                                 show_subtitle = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_delta_activity_condition_ranges() requires package 'ggplot2'.", call. = FALSE)
  }
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(activity_scores)) || any(!nzchar(rownames(activity_scores)))) {
    stop("activity_scores must have non-empty TE subfamily rownames.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }

  if (!is.data.frame(delta_table) || !all(c("te_subfamily", "delta_activity") %in% colnames(delta_table))) {
    stop("delta_table must be a data.frame with te_subfamily and delta_activity columns.", call. = FALSE)
  }
  delta_table$te_subfamily <- as.character(delta_table$te_subfamily)
  delta_table$delta_activity <- suppressWarnings(as.numeric(delta_table$delta_activity))
  delta_table <- delta_table[nzchar(delta_table$te_subfamily) & !is.na(delta_table$delta_activity), , drop = FALSE]
  delta_table <- delta_table[!duplicated(delta_table$te_subfamily), , drop = FALSE]

  missing_tes <- setdiff(delta_table$te_subfamily, rownames(activity_scores))
  if (length(missing_tes) > 0L) {
    warning(
      sprintf("Dropping TE subfamilies not found in activity_scores rownames: %s", paste(missing_tes, collapse = ", ")),
      call. = FALSE
    )
    delta_table <- delta_table[delta_table$te_subfamily %in% rownames(activity_scores), , drop = FALSE]
  }
  if (nrow(delta_table) == 0L) {
    stop("No TE subfamilies from delta_table are available in activity_scores.", call. = FALSE)
  }

  if (!is.list(condition_groups) || length(condition_groups) == 0L) {
    stop("condition_groups must be a named list of sample vectors.", call. = FALSE)
  }
  if (is.null(names(condition_groups)) || any(!nzchar(names(condition_groups)))) {
    stop("condition_groups must have non-empty names.", call. = FALSE)
  }
  if (anyDuplicated(names(condition_groups))) {
    stop("condition_groups names must be unique.", call. = FALSE)
  }

  if (is.null(condition_order)) {
    condition_order <- names(condition_groups)
  }
  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[condition_order %in% names(condition_groups)]
  if (length(condition_order) == 0L) {
    stop("condition_order must include at least one name from condition_groups.", call. = FALSE)
  }
  condition_groups <- condition_groups[condition_order]

  for (condition in names(condition_groups)) {
    samples <- condition_groups[[condition]]
    if (!is.character(samples)) {
      stop(sprintf("condition_groups[['%s']] must be a character vector.", condition), call. = FALSE)
    }
    samples <- unique(samples[nzchar(samples)])
    if (length(samples) == 0L) {
      stop(sprintf("condition_groups[['%s']] must contain at least one sample.", condition), call. = FALSE)
    }
    condition_groups[[condition]] <- samples
  }

  missing_samples <- setdiff(unique(unlist(condition_groups, use.names = FALSE)), colnames(activity_scores))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("Samples not found in activity_scores colnames: %s", paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  te_order <- unique(delta_table$te_subfamily[order(abs(delta_table$delta_activity), decreasing = TRUE)])
  finite_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else min(x)
  }
  finite_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }
  condition_ranges <- do.call(rbind, lapply(names(condition_groups), function(condition) {
    samples <- condition_groups[[condition]]
    sample_matrix <- activity_scores[te_order, samples, drop = FALSE]
    data.frame(
      te_subfamily = te_order,
      condition = condition,
      min_activity = apply(sample_matrix, 1, finite_min),
      max_activity = apply(sample_matrix, 1, finite_max),
      mean_activity = rowMeans(sample_matrix, na.rm = TRUE),
      n_samples = rowSums(is.finite(sample_matrix)),
      stringsAsFactors = FALSE
    )
  }))
  condition_ranges <- condition_ranges[
    is.finite(condition_ranges$min_activity) &
      is.finite(condition_ranges$max_activity) &
      is.finite(condition_ranges$mean_activity),
    ,
    drop = FALSE
  ]
  if (nrow(condition_ranges) == 0L) {
    stop("No finite condition activity ranges are available for the selected TE subfamilies.", call. = FALSE)
  }

  te_order <- te_order[te_order %in% unique(condition_ranges$te_subfamily)]
  condition_ranges <- condition_ranges[condition_ranges$te_subfamily %in% te_order, , drop = FALSE]
  levels_te <- rev(te_order)
  condition_ranges$condition <- factor(condition_ranges$condition, levels = condition_order)

  base_y <- match(condition_ranges$te_subfamily, levels_te)
  if (length(condition_order) == 1L) {
    condition_offsets <- stats::setNames(0, condition_order)
  } else {
    condition_offsets <- stats::setNames(seq(0.24, -0.24, length.out = length(condition_order)), condition_order)
  }
  condition_ranges$y_position <- base_y + unname(condition_offsets[as.character(condition_ranges$condition)])
  cap_height <- if (length(condition_order) <= 1L) 0.10 else min(0.08, 0.22 / length(condition_order))

  x_max <- max(condition_ranges$max_activity, na.rm = TRUE)
  x_min <- min(condition_ranges$min_activity, na.rm = TRUE)
  x_lower <- min(0, x_min)
  x_upper <- max(0, x_max)
  x_span <- max(x_upper - x_lower, 0.01)
  axis_pad <- x_span * 0.08

  if (is.null(condition_cols)) {
    condition_color_values <- stats::setNames(
      grDevices::hcl.colors(length(condition_order), palette = "Dark 3"),
      condition_order
    )
  } else {
    condition_cols <- as.character(condition_cols)
    if (is.null(names(condition_cols)) || all(!nzchar(names(condition_cols)))) {
      names(condition_cols) <- condition_order[seq_along(condition_cols)]
    }
    condition_color_values <- condition_cols[condition_order]
    missing_cols <- is.na(condition_color_values) | !nzchar(condition_color_values)
    if (any(missing_cols)) {
      condition_color_values[missing_cols] <- grDevices::hcl.colors(sum(missing_cols), palette = "Dark 3")
    }
  }

  if (is.null(title) || length(title) != 1L || !nzchar(title)) {
    title <- condition_label
  }
  if (is.null(selection_label) || length(selection_label) != 1L || !nzchar(selection_label)) {
    selection_label <- sprintf("abs(%s)", condition_label)
  }
  subtitle_text <- if (isTRUE(show_subtitle)) {
    sprintf("Top %s TEs by %s; brackets = replicate min-max, dots = condition means.", length(te_order), selection_label)
  } else {
    NULL
  }

  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey35") +
    ggplot2::geom_segment(
      data = condition_ranges,
      ggplot2::aes(x = min_activity, xend = max_activity, y = y_position, yend = y_position, color = condition),
      linewidth = 1.1,
      alpha = 0.85,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = condition_ranges,
      ggplot2::aes(x = min_activity, xend = min_activity, y = y_position - cap_height, yend = y_position + cap_height, color = condition),
      linewidth = 0.75,
      alpha = 0.85
    ) +
    ggplot2::geom_segment(
      data = condition_ranges,
      ggplot2::aes(x = max_activity, xend = max_activity, y = y_position - cap_height, yend = y_position + cap_height, color = condition),
      linewidth = 0.75,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      data = condition_ranges,
      ggplot2::aes(x = mean_activity, y = y_position, color = condition),
      size = 2.35,
      alpha = 0.95
    ) +
    ggplot2::scale_x_continuous(limits = c(x_lower - axis_pad, x_upper + axis_pad)) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(levels_te),
      labels = levels_te,
      limits = c(0.5, length(levels_te) + 0.5),
      expand = ggplot2::expansion(add = c(0.05, 0.05))
    ) +
    ggplot2::scale_color_manual(values = condition_color_values, breaks = condition_order, drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle_text,
      x = x_label,
      y = "TE subfamily",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = 9),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5.5, 28, 5.5, 5.5)
    )
}

# Select one shared TE list for multi-condition activity-change plots.
select_top_te_subfamilies_by_activity_change <- function(activity_scores,
                                                         condition_groups,
                                                         condition_order = names(condition_groups),
                                                         top_k = 10,
                                                         ranking_metric = c(
                                                           "max_step_abs_delta",
                                                           "sum_step_abs_delta",
                                                           "max_sum_abs_activity_score",
                                                           "abs_delta_ESEN_vs_PRO",
                                                           "abs_delta_LSEN_vs_ESEN",
                                                           "abs_delta_LSEN_vs_PRO",
                                                           "net_abs_delta",
                                                           "mean_activity_range",
                                                           "monotonic_step_abs_delta",
                                                           "PRO_ESEN",
                                                           "ESEN_LSEN",
                                                           "PRO_LSEN"
                                                         ),
                                                         te_subfamilies = NULL) {
  if (missing(ranking_metric)) {
    ranking_metric <- "max_step_abs_delta"
  } else {
    ranking_metric <- as.character(ranking_metric)
  }
  if (length(ranking_metric) != 1L || is.na(ranking_metric) || !nzchar(ranking_metric)) {
    stop("ranking_metric must be a single non-empty character value.", call. = FALSE)
  }
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(activity_scores)) || any(!nzchar(rownames(activity_scores)))) {
    stop("activity_scores must have non-empty TE subfamily rownames.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }
  if (!is.list(condition_groups) || length(condition_groups) == 0L) {
    stop("condition_groups must be a named list of sample vectors.", call. = FALSE)
  }
  if (is.null(names(condition_groups)) || any(!nzchar(names(condition_groups)))) {
    stop("condition_groups must have non-empty names.", call. = FALSE)
  }
  if (anyDuplicated(names(condition_groups))) {
    stop("condition_groups names must be unique.", call. = FALSE)
  }
  if (is.null(condition_order)) {
    condition_order <- names(condition_groups)
  }
  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[condition_order %in% names(condition_groups)]
  if (length(condition_order) < 2L) {
    stop("condition_order must include at least two names from condition_groups.", call. = FALSE)
  }
  condition_groups <- condition_groups[condition_order]
  pair_gap_aliases <- stats::setNames(character(0), character(0))
  if (length(condition_order) >= 2L) {
    for (from_idx in seq_len(length(condition_order) - 1L)) {
      for (to_idx in seq.int(from_idx + 1L, length(condition_order))) {
        from_condition <- condition_order[from_idx]
        to_condition <- condition_order[to_idx]
        canonical_metric <- paste0("abs_delta_", to_condition, "_vs_", from_condition)
        pair_gap_aliases[paste(from_condition, to_condition, sep = "_")] <- canonical_metric
        pair_gap_aliases[paste(to_condition, from_condition, sep = "_")] <- canonical_metric
      }
    }
  }
  if (ranking_metric %in% names(pair_gap_aliases)) {
    ranking_metric <- unname(pair_gap_aliases[ranking_metric])
  }

  for (condition in names(condition_groups)) {
    samples <- condition_groups[[condition]]
    if (!is.character(samples)) {
      stop(sprintf("condition_groups[['%s']] must be a character vector.", condition), call. = FALSE)
    }
    samples <- unique(samples[nzchar(samples)])
    if (length(samples) == 0L) {
      stop(sprintf("condition_groups[['%s']] must contain at least one sample.", condition), call. = FALSE)
    }
    condition_groups[[condition]] <- samples
  }

  missing_samples <- setdiff(unique(unlist(condition_groups, use.names = FALSE)), colnames(activity_scores))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("Samples not found in activity_scores colnames: %s", paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  if (is.null(top_k)) {
    top_k <- 0L
  }
  if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 0L) {
    stop("top_k must be a single non-negative number or NULL.", call. = FALSE)
  }
  top_k <- as.integer(top_k)

  if (!is.null(te_subfamilies)) {
    te_subfamilies <- unique(as.character(te_subfamilies))
    te_subfamilies <- te_subfamilies[nzchar(te_subfamilies)]
  }
  if (top_k == 0L && length(te_subfamilies) == 0L) {
    stop("Set top_k > 0 and/or provide te_subfamilies.", call. = FALSE)
  }

  condition_means <- do.call(cbind, lapply(condition_groups, function(samples) {
    rowMeans(activity_scores[, samples, drop = FALSE], na.rm = TRUE)
  }))
  colnames(condition_means) <- names(condition_groups)

  step_from <- head(condition_order, -1L)
  step_to <- tail(condition_order, -1L)
  step_delta <- do.call(cbind, lapply(seq_along(step_from), function(i) {
    condition_means[, step_to[i]] - condition_means[, step_from[i]]
  }))
  step_delta_cols <- paste0("delta_", step_to, "_vs_", step_from)
  colnames(step_delta) <- step_delta_cols

  finite_row_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }
  finite_row_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else min(x)
  }
  finite_row_sum <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else sum(x)
  }
  same_direction <- function(x) {
    x <- x[is.finite(x) & x != 0]
    length(x) > 0L && (all(x > 0) || all(x < 0))
  }

  abs_step_delta <- abs(step_delta)
  condition_abs_activity_sums <- do.call(cbind, lapply(condition_groups, function(samples) {
    apply(abs(activity_scores[, samples, drop = FALSE]), 1, finite_row_sum)
  }))
  colnames(condition_abs_activity_sums) <- names(condition_groups)

  summary_table <- data.frame(
    te_subfamily = rownames(activity_scores),
    stringsAsFactors = FALSE
  )
  for (condition in condition_order) {
    summary_table[[paste0("mean_", condition)]] <- unname(condition_means[, condition])
    summary_table[[paste0("sum_abs_activity_", condition)]] <- unname(condition_abs_activity_sums[, condition])
  }
  for (delta_col in step_delta_cols) {
    summary_table[[delta_col]] <- unname(step_delta[, delta_col])
  }

  net_delta_col <- paste0("delta_", tail(condition_order, 1), "_vs_", condition_order[1])
  summary_table[[net_delta_col]] <- unname(condition_means[, tail(condition_order, 1)] - condition_means[, condition_order[1]])
  for (delta_col in unique(c(step_delta_cols, net_delta_col))) {
    summary_table[[paste0("abs_", delta_col)]] <- abs(summary_table[[delta_col]])
  }
  summary_table$max_step_abs_delta <- apply(abs_step_delta, 1, finite_row_max)
  summary_table$sum_step_abs_delta <- apply(abs_step_delta, 1, finite_row_sum)
  summary_table$max_sum_abs_activity_score <- apply(condition_abs_activity_sums, 1, finite_row_max)
  summary_table$net_abs_delta <- summary_table[[paste0("abs_", net_delta_col)]]
  summary_table$mean_activity_range <- apply(condition_means, 1, function(x) finite_row_max(x) - finite_row_min(x))
  summary_table$monotonic_step_abs_delta <- ifelse(
    apply(step_delta, 1, same_direction),
    apply(abs_step_delta, 1, finite_row_min),
    0
  )

  available_ranking_metrics <- c(
    "max_step_abs_delta",
    "sum_step_abs_delta",
    "max_sum_abs_activity_score",
    paste0("abs_", unique(c(step_delta_cols, net_delta_col))),
    "net_abs_delta",
    "mean_activity_range",
    "monotonic_step_abs_delta"
  )
  if (!ranking_metric %in% available_ranking_metrics) {
    stop(
      sprintf(
        "ranking_metric must be one of: %s",
        paste(c(available_ranking_metrics, names(pair_gap_aliases)), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  summary_table$ranking_metric <- ranking_metric
  summary_table$ranking_score <- summary_table[[ranking_metric]]

  missing_tes <- setdiff(te_subfamilies, rownames(activity_scores))
  if (length(missing_tes) > 0L) {
    warning(
      sprintf("Requested TE subfamilies not found in activity_scores rownames: %s", paste(missing_tes, collapse = ", ")),
      call. = FALSE
    )
  }

  selected_by_top_k <- rep(FALSE, nrow(summary_table))
  finite_scores <- is.finite(summary_table$ranking_score)
  if (top_k > 0L && any(finite_scores)) {
    finite_rows <- which(finite_scores)
    top_rows <- finite_rows[order(-summary_table$ranking_score[finite_rows], summary_table$te_subfamily[finite_rows])[
      seq_len(min(top_k, length(finite_rows)))
    ]]
    selected_by_top_k[top_rows] <- TRUE
  }
  selected_by_request <- summary_table$te_subfamily %in% te_subfamilies
  summary_table$selected_by_top_k <- selected_by_top_k
  summary_table$selected_by_request <- selected_by_request

  selected_table <- summary_table[selected_by_top_k | selected_by_request, , drop = FALSE]
  if (nrow(selected_table) == 0L) {
    stop("No TE subfamilies selected. Check top_k, ranking_metric, and te_subfamilies.", call. = FALSE)
  }
  selected_table <- selected_table[
    order(-selected_table$ranking_score, selected_table$te_subfamily),
    ,
    drop = FALSE
  ]
  rownames(selected_table) <- NULL
  selected_table
}

# Build one all-TE table that combines mean activities, pairwise deltas,
# significance p-values, and threshold-based senescence trend classes.
build_te_contrast_summary <- function(activity_scores,
                                      sample_groups,
                                      contrasts,
                                      signif_results,
                                      delta_threshold = 0.05,
                                      signif_alpha = 0.05) {
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(activity_scores)) || any(!nzchar(rownames(activity_scores)))) {
    stop("activity_scores must have non-empty TE subfamily rownames.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }
  if (!is.numeric(delta_threshold) || length(delta_threshold) != 1L || is.na(delta_threshold) || delta_threshold < 0) {
    stop("delta_threshold must be one non-negative number.", call. = FALSE)
  }
  if (!is.numeric(signif_alpha) || length(signif_alpha) != 1L || is.na(signif_alpha) || signif_alpha < 0 || signif_alpha > 1) {
    stop("signif_alpha must be one number between 0 and 1.", call. = FALSE)
  }
  if (!is.list(sample_groups) || is.null(names(sample_groups)) || any(!nzchar(names(sample_groups)))) {
    stop("sample_groups must be a named list.", call. = FALSE)
  }
  required_groups <- c("PRO", "ESEN", "LSEN")
  missing_groups <- setdiff(required_groups, names(sample_groups))
  if (length(missing_groups) > 0L) {
    stop(sprintf("sample_groups is missing required groups: %s", paste(missing_groups, collapse = ", ")), call. = FALSE)
  }
  if (!is.list(contrasts) || is.null(names(contrasts)) || any(!nzchar(names(contrasts)))) {
    stop("contrasts must be a named list.", call. = FALSE)
  }
  required_contrasts <- c("ESEN_vs_PRO", "LSEN_vs_PRO", "LSEN_vs_ESEN")
  missing_contrasts <- setdiff(required_contrasts, names(contrasts))
  if (length(missing_contrasts) > 0L) {
    stop(sprintf("contrasts is missing required contrasts: %s", paste(missing_contrasts, collapse = ", ")), call. = FALSE)
  }
  if (!is.list(signif_results) || is.null(names(signif_results))) {
    stop("signif_results must be a named list.", call. = FALSE)
  }
  missing_signif <- setdiff(required_contrasts, names(signif_results))
  if (length(missing_signif) > 0L) {
    stop(sprintf("signif_results is missing required contrasts: %s", paste(missing_signif, collapse = ", ")), call. = FALSE)
  }

  validate_samples <- function(samples, group_name) {
    samples <- unique(as.character(samples))
    samples <- samples[nzchar(samples)]
    if (length(samples) == 0L) {
      stop(sprintf("sample_groups[['%s']] must contain at least one sample.", group_name), call. = FALSE)
    }
    missing_samples <- setdiff(samples, colnames(activity_scores))
    if (length(missing_samples) > 0L) {
      stop(
        sprintf(
          "sample_groups[['%s']] contains samples not found in activity_scores: %s",
          group_name,
          paste(missing_samples, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    samples
  }
  sample_groups <- stats::setNames(
    lapply(names(sample_groups), function(group_name) validate_samples(sample_groups[[group_name]], group_name)),
    names(sample_groups)
  )

  get_contrast_groups <- function(contrast_name) {
    contrast <- contrasts[[contrast_name]]
    if (is.null(contrast$treatment) || is.null(contrast$control)) {
      stop(sprintf("contrasts[['%s']] must contain treatment and control group names.", contrast_name), call. = FALSE)
    }
    treatment <- as.character(contrast$treatment)
    control <- as.character(contrast$control)
    if (!treatment %in% names(sample_groups) || !control %in% names(sample_groups)) {
      stop(sprintf("contrasts[['%s']] references unknown sample group(s).", contrast_name), call. = FALSE)
    }
    list(treatment = treatment, control = control)
  }

  te_subfamilies <- rownames(activity_scores)
  summary_table <- data.frame(
    te_subfamily = te_subfamilies,
    mean_PRO = unname(rowMeans(activity_scores[, sample_groups$PRO, drop = FALSE], na.rm = TRUE)),
    mean_ESEN = unname(rowMeans(activity_scores[, sample_groups$ESEN, drop = FALSE], na.rm = TRUE)),
    mean_LSEN = unname(rowMeans(activity_scores[, sample_groups$LSEN, drop = FALSE], na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  for (contrast_name in required_contrasts) {
    contrast_groups <- get_contrast_groups(contrast_name)
    treatment_mean <- rowMeans(activity_scores[, sample_groups[[contrast_groups$treatment]], drop = FALSE], na.rm = TRUE)
    control_mean <- rowMeans(activity_scores[, sample_groups[[contrast_groups$control]], drop = FALSE], na.rm = TRUE)
    summary_table[[paste0("delta_", contrast_name)]] <- unname(treatment_mean - control_mean)
  }

  for (contrast_name in required_contrasts) {
    signif_table <- signif_activity_table(signif_results[[contrast_name]], include_intercept = FALSE)
    if (anyDuplicated(signif_table$te_subfamily)) {
      duplicated_tes <- unique(signif_table$te_subfamily[duplicated(signif_table$te_subfamily)])
      stop(
        sprintf("signif_results[['%s']] contains duplicated TE subfamilies: %s", contrast_name, paste(duplicated_tes, collapse = ", ")),
        call. = FALSE
      )
    }
    p_value <- stats::setNames(signif_table$p_value, signif_table$te_subfamily)
    padj <- stats::setNames(signif_table$p_adj, signif_table$te_subfamily)
    summary_table[[paste0("p_value_", contrast_name)]] <- unname(p_value[te_subfamilies])
    summary_table[[paste0("padj_", contrast_name)]] <- unname(padj[te_subfamilies])
    summary_table[[paste0("signif_", contrast_name)]] <- !is.na(summary_table[[paste0("padj_", contrast_name)]]) &
      summary_table[[paste0("padj_", contrast_name)]] < signif_alpha
  }

  summary_table$trend_class <- classify_te_senescence_trend(
    mean_PRO = summary_table$mean_PRO,
    mean_ESEN = summary_table$mean_ESEN,
    mean_LSEN = summary_table$mean_LSEN,
    delta_threshold = delta_threshold
  )
  signif_cols <- paste0("signif_", required_contrasts)
  summary_table$n_significant_contrasts <- rowSums(summary_table[, signif_cols, drop = FALSE], na.rm = TRUE)
  delta_cols <- paste0("delta_", required_contrasts)
  summary_table$max_abs_delta <- apply(abs(summary_table[, delta_cols, drop = FALSE]), 1, max, na.rm = TRUE)

  requested_order <- c(
    "te_subfamily",
    "mean_PRO",
    "mean_ESEN",
    "mean_LSEN",
    "delta_ESEN_vs_PRO",
    "delta_LSEN_vs_PRO",
    "delta_LSEN_vs_ESEN",
    "p_value_ESEN_vs_PRO",
    "p_value_LSEN_vs_PRO",
    "p_value_LSEN_vs_ESEN",
    "padj_ESEN_vs_PRO",
    "padj_LSEN_vs_PRO",
    "padj_LSEN_vs_ESEN",
    "signif_ESEN_vs_PRO",
    "signif_LSEN_vs_PRO",
    "signif_LSEN_vs_ESEN",
    "trend_class",
    "n_significant_contrasts",
    "max_abs_delta"
  )
  summary_table <- summary_table[, requested_order, drop = FALSE]
  summary_table <- summary_table[
    order(-summary_table$n_significant_contrasts, -summary_table$max_abs_delta, summary_table$te_subfamily),
    ,
    drop = FALSE
  ]
  rownames(summary_table) <- NULL
  summary_table
}

# Plot the size of each trend class and the delta distribution used to tune the threshold.
plot_te_trend_summary <- function(te_contrast_summary,
                                  delta_threshold = 0.05) {
  required_cols <- c("trend_class", "max_abs_delta")
  missing_cols <- setdiff(required_cols, colnames(te_contrast_summary))
  if (length(missing_cols) > 0L) {
    stop(sprintf("te_contrast_summary is missing columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(delta_threshold) || length(delta_threshold) != 1L || is.na(delta_threshold) || delta_threshold < 0) {
    stop("delta_threshold must be one non-negative number.", call. = FALSE)
  }

  trend_counts <- sort(table(te_contrast_summary$trend_class), decreasing = TRUE)
  n_total <- sum(trend_counts)
  n_noise <- if ("flat_or_noisy" %in% names(trend_counts)) {
    unname(trend_counts[["flat_or_noisy"]])
  } else {
    0L
  }
  n_non_noise <- n_total - n_noise
  pct_non_noise <- if (n_total > 0L) {
    100 * n_non_noise / n_total
  } else {
    NA_real_
  }
  non_noise_label <- sprintf(
    "Non-noise classified: %d/%d TE families (%.1f%%)",
    n_non_noise,
    n_total,
    pct_non_noise
  )

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(8, 4, 4.5, 1))
  trend_ylim <- c(0, max(trend_counts) * 1.18)
  trend_mid <- graphics::barplot(
    trend_counts,
    las = 2,
    col = "#2F6F73",
    ylab = "TE subfamilies",
    main = "TE trend classes",
    ylim = trend_ylim
  )
  graphics::text(
    x = trend_mid,
    y = as.numeric(trend_counts),
    labels = as.integer(trend_counts),
    pos = 3,
    cex = 0.85,
    xpd = TRUE
  )
  graphics::mtext(
    non_noise_label,
    side = 3,
    line = 0.5,
    cex = 0.85
  )
  delta_values <- te_contrast_summary$max_abs_delta
  delta_hist_counts <- graphics::hist(delta_values, breaks = 40, plot = FALSE)$counts
  n_above_threshold <- sum(delta_values >= delta_threshold, na.rm = TRUE)
  n_below_threshold <- sum(delta_values < delta_threshold, na.rm = TRUE)
  threshold_summary <- sprintf(
    "< threshold: %d | >= threshold: %d",
    n_below_threshold,
    n_above_threshold
  )

  delta_hist <- graphics::hist(
    te_contrast_summary$max_abs_delta,
    breaks = 40,
    col = "#9DB4C0",
    border = "white",
    xlab = "Maximum absolute pairwise delta",
    main = "Delta threshold diagnostic",
    ylim = c(0, max(delta_hist_counts) * 1.18)
  )
  graphics::mtext(
    threshold_summary,
    side = 3,
    line = 0.5,
    cex = 0.85
  )
  graphics::abline(v = delta_threshold, col = "#D1495B", lwd = 2, lty = 2)
  graphics::legend(
    "topright",
    legend = sprintf("threshold = %.3f", delta_threshold),
    col = "#D1495B",
    lwd = 2,
    lty = 2,
    bty = "n"
  )

  invisible(data.frame(
    trend_class = names(trend_counts),
    n = as.integer(trend_counts),
    n_total = as.integer(n_total),
    n_non_noise = as.integer(n_non_noise),
    pct_non_noise = as.numeric(pct_non_noise),
    stringsAsFactors = FALSE
  ))
}

# Calculate correlations between activity score vectors for named sample pairs.
# sample_pairs is a named list such as list("Sp VPR rep1 vs rep2" =
# c("LTR5HS-Sp-VPR-rep1", "LTR5HS-Sp-VPR-rep2")). This keeps the helper usable
# for replicate checks and for cross-system checks such as Sp versus Sa.
activity_pair_correlations <- function(activity_scores,
                                       sample_pairs,
                                       method = "pearson",
                                       use = "pairwise.complete.obs",
                                       print_table = TRUE,
                                       heading = "Replicate correlations") {
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }
  if (!is.list(sample_pairs) || length(sample_pairs) == 0L) {
    stop("sample_pairs must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(sample_pairs)) || any(!nzchar(names(sample_pairs)))) {
    stop("sample_pairs must be a named list so each correlation has a label.", call. = FALSE)
  }
  if (anyDuplicated(names(sample_pairs))) {
    duplicated_names <- unique(names(sample_pairs)[duplicated(names(sample_pairs))])
    stop(
      sprintf("sample_pairs contains duplicated labels: %s", paste(duplicated_names, collapse = ", ")),
      call. = FALSE
    )
  }

  pair_table <- do.call(
    rbind,
    lapply(names(sample_pairs), function(label) {
      pair <- as.character(sample_pairs[[label]])
      pair <- pair[nzchar(pair)]
      if (length(pair) != 2L) {
        stop(sprintf("sample_pairs[['%s']] must contain exactly two sample names.", label), call. = FALSE)
      }
      missing_samples <- setdiff(pair, colnames(activity_scores))
      if (length(missing_samples) > 0L) {
        stop(
          sprintf(
            "sample_pairs[['%s']] contains samples not found in activity_scores: %s",
            label,
            paste(missing_samples, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      data.frame(
        comparison = label,
        sample_1 = pair[1],
        sample_2 = pair[2],
        correlation = stats::cor(activity_scores[, pair[1]], activity_scores[, pair[2]], method = method, use = use),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(pair_table) <- NULL

  if (isTRUE(print_table)) {
    cat(sprintf("\n── %s ──────────────────\n", heading))
    for (i in seq_len(nrow(pair_table))) {
      cat(sprintf("%-28s %.3f\n", paste0(pair_table$comparison[i], ":"), pair_table$correlation[i]))
    }
  }

  pair_table
}

# Summarise activity-score correlations within and between condition groups.
activity_condition_correlations <- function(activity_scores,
                                            sample_groups,
                                            method = "pearson",
                                            use = "pairwise.complete.obs",
                                            print_table = TRUE,
                                            heading = "Activity correlation summary") {
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }
  if (!is.list(sample_groups) || length(sample_groups) == 0L) {
    stop("sample_groups must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(sample_groups)) || any(!nzchar(names(sample_groups)))) {
    stop("sample_groups must be a named list.", call. = FALSE)
  }
  if (anyDuplicated(names(sample_groups))) {
    stop("sample_groups names must be unique.", call. = FALSE)
  }

  for (condition in names(sample_groups)) {
    samples <- as.character(sample_groups[[condition]])
    samples <- unique(samples[nzchar(samples)])
    if (length(samples) == 0L) {
      stop(sprintf("sample_groups[['%s']] must contain at least one sample.", condition), call. = FALSE)
    }
    sample_groups[[condition]] <- samples
  }

  missing_samples <- setdiff(unique(unlist(sample_groups, use.names = FALSE)), colnames(activity_scores))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("sample_groups contains samples not found in activity_scores: %s", paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  cor_one_pair <- function(sample_1, sample_2) {
    stats::cor(activity_scores[, sample_1], activity_scores[, sample_2], method = method, use = use)
  }

  make_pair_rows <- function(scope, comparison, sample_pairs) {
    if (nrow(sample_pairs) == 0L) {
      return(NULL)
    }
    data.frame(
      scope = scope,
      comparison = comparison,
      sample_1 = sample_pairs$sample_1,
      sample_2 = sample_pairs$sample_2,
      correlation = mapply(cor_one_pair, sample_pairs$sample_1, sample_pairs$sample_2),
      stringsAsFactors = FALSE
    )
  }

  within_rows <- do.call(rbind, lapply(names(sample_groups), function(condition) {
    samples <- sample_groups[[condition]]
    if (length(samples) < 2L) {
      return(NULL)
    }
    pairs <- utils::combn(samples, 2, simplify = FALSE)
    make_pair_rows(
      scope = "within_condition",
      comparison = condition,
      sample_pairs = data.frame(
        sample_1 = vapply(pairs, `[`, character(1), 1),
        sample_2 = vapply(pairs, `[`, character(1), 2),
        stringsAsFactors = FALSE
      )
    )
  }))

  condition_pairs <- utils::combn(names(sample_groups), 2, simplify = FALSE)
  between_rows <- do.call(rbind, lapply(condition_pairs, function(condition_pair) {
    sample_pairs <- expand.grid(
      sample_1 = sample_groups[[condition_pair[1]]],
      sample_2 = sample_groups[[condition_pair[2]]],
      stringsAsFactors = FALSE
    )
    make_pair_rows(
      scope = "between_conditions",
      comparison = paste(condition_pair, collapse = " vs "),
      sample_pairs = sample_pairs
    )
  }))

  pair_table <- rbind(within_rows, between_rows)
  rownames(pair_table) <- NULL

  summary_table <- do.call(rbind, lapply(split(pair_table, list(pair_table$scope, pair_table$comparison), drop = TRUE), function(x) {
    data.frame(
      scope = unique(x$scope),
      comparison = unique(x$comparison),
      n_pairs = nrow(x),
      mean_correlation = mean(x$correlation, na.rm = TRUE),
      min_correlation = min(x$correlation, na.rm = TRUE),
      max_correlation = max(x$correlation, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_table) <- NULL
  summary_table <- summary_table[order(summary_table$scope, summary_table$comparison), , drop = FALSE]

  sample_summary <- do.call(rbind, lapply(names(sample_groups), function(condition) {
    samples <- sample_groups[[condition]]
    do.call(rbind, lapply(samples, function(sample_id) {
      sample_rows <- within_rows[
        within_rows$comparison == condition &
          (within_rows$sample_1 == sample_id | within_rows$sample_2 == sample_id),
        ,
        drop = FALSE
      ]
      data.frame(
        condition = condition,
        sample_id = sample_id,
        n_pairs = nrow(sample_rows),
        mean_within_correlation = if (nrow(sample_rows) == 0L) NA_real_ else mean(sample_rows$correlation, na.rm = TRUE),
        min_within_correlation = if (nrow(sample_rows) == 0L) NA_real_ else min(sample_rows$correlation, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(sample_summary) <- NULL

  if (isTRUE(print_table)) {
    cat(sprintf("\n── %s ──────────────────\n", heading))
    cat("Within condition means:\n")
    within_summary <- summary_table[summary_table$scope == "within_condition", , drop = FALSE]
    for (i in seq_len(nrow(within_summary))) {
      cat(sprintf(
        "%-14s %.3f  (%d pairs)\n",
        paste0(within_summary$comparison[i], ":"),
        within_summary$mean_correlation[i],
        within_summary$n_pairs[i]
      ))
    }

    cat("\nWithin-condition sample agreement:\n")
    for (condition in names(sample_groups)) {
      condition_samples <- sample_summary[sample_summary$condition == condition, , drop = FALSE]
      cat(sprintf(
        "%s: %s\n",
        condition,
        paste(
          sprintf("%s=%.3f", condition_samples$sample_id, condition_samples$mean_within_correlation),
          collapse = ", "
        )
      ))
    }

    weakest_samples <- do.call(rbind, lapply(split(sample_summary, sample_summary$condition), function(x) {
      x[which.min(x$mean_within_correlation), , drop = FALSE]
    }))
    rownames(weakest_samples) <- NULL
    cat("\nLowest within-condition sample agreement:\n")
    for (i in seq_len(nrow(weakest_samples))) {
      cat(sprintf(
        "%-14s %s %.3f\n",
        paste0(weakest_samples$condition[i], ":"),
        weakest_samples$sample_id[i],
        weakest_samples$mean_within_correlation[i]
      ))
    }

    cat("\nBetween condition means:\n")
    between_summary <- summary_table[summary_table$scope == "between_conditions", , drop = FALSE]
    for (i in seq_len(nrow(between_summary))) {
      cat(sprintf(
        "%-14s %.3f  (%d pairs)\n",
        paste0(between_summary$comparison[i], ":"),
        between_summary$mean_correlation[i],
        between_summary$n_pairs[i]
      ))
    }
  }

  list(summary = summary_table, sample_summary = sample_summary, pairs = pair_table)
}

# Plot selected TE-subfamily activity scores across samples ordered by group.
plot_activity_across_samples <- function(activity_scores,
                                         te_subfamilies,
                                         sample_groups = NULL,
                                         group_cols = c(
                                           "control_KRAB" = "tomato",
                                           "control_VPR" = "steelblue",
                                           "target_KRAB" = "darkred",
                                           "target_VPR" = "darkblue"
                                         ),
                                         control_pattern = "nontarget",
                                         vpr_pattern = "VPR",
                                         krab_pattern = "KRAB",
                                         group_order = names(group_cols),
                                         main_template = "%s activity across samples",
                                         xlab = "Samples (ordered by group)",
                                         ylab = "Activity",
                                         pch = 16,
                                         zero_line = TRUE,
                                         zero_line_col = "grey50",
                                         zero_line_lty = 2,
                                         legend_position = "topright",
                                         legend_cex = 0.8) {
  if (is.null(activity_scores)) {
    stop("activity_scores must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(activity_scores) && !is.data.frame(activity_scores) && !inherits(activity_scores, "Matrix")) {
    stop("activity_scores must be a matrix or data.frame.", call. = FALSE)
  }
  activity_scores <- as.matrix(activity_scores)
  if (!is.numeric(activity_scores)) {
    stop("activity_scores must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(activity_scores)) || any(!nzchar(rownames(activity_scores)))) {
    stop("activity_scores must have non-empty TE subfamily rownames.", call. = FALSE)
  }
  if (is.null(colnames(activity_scores)) || any(!nzchar(colnames(activity_scores)))) {
    stop("activity_scores must have non-empty sample colnames.", call. = FALSE)
  }

  samples <- colnames(activity_scores)
  if (anyDuplicated(samples)) {
    duplicated_samples <- unique(samples[duplicated(samples)])
    stop(
      sprintf("activity_scores contains duplicated sample colnames: %s", paste(duplicated_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  te_subfamilies <- unique(as.character(te_subfamilies))
  te_subfamilies <- te_subfamilies[nzchar(te_subfamilies)]
  if (length(te_subfamilies) == 0L) {
    stop("te_subfamilies must contain at least one TE subfamily.", call. = FALSE)
  }
  missing_tes <- setdiff(te_subfamilies, rownames(activity_scores))
  if (length(missing_tes) > 0L) {
    warning(
      sprintf("Dropping TE subfamilies not found in activity_scores rownames: %s", paste(missing_tes, collapse = ", ")),
      call. = FALSE
    )
  }
  te_subfamilies <- te_subfamilies[te_subfamilies %in% rownames(activity_scores)]
  if (length(te_subfamilies) == 0L) {
    stop("None of the requested TE subfamilies are present in activity_scores.", call. = FALSE)
  }

  if (is.null(sample_groups)) {
    system <- ifelse(
      grepl(vpr_pattern, samples),
      "VPR",
      ifelse(grepl(krab_pattern, samples), "KRAB", NA_character_)
    )
    missing_system <- samples[is.na(system)]
    if (length(missing_system) > 0L) {
      stop(
        sprintf("Could not infer VPR/KRAB system for samples: %s", paste(missing_system, collapse = ", ")),
        call. = FALSE
      )
    }

    perturbation <- ifelse(grepl(control_pattern, samples), "control", "target")
    sample_groups <- stats::setNames(paste(perturbation, system, sep = "_"), samples)
  }

  if (is.null(names(sample_groups))) {
    stop("sample_groups must be a named vector when provided.", call. = FALSE)
  }
  sample_group_names <- names(sample_groups)
  sample_groups <- stats::setNames(as.character(sample_groups), sample_group_names)
  group_samples <- names(sample_groups)
  if (any(!nzchar(group_samples))) {
    stop("Every sample_groups entry must have a non-empty sample name.", call. = FALSE)
  }
  if (anyDuplicated(group_samples)) {
    duplicated_group_samples <- unique(group_samples[duplicated(group_samples)])
    stop(
      sprintf("sample_groups assigns more than one group to: %s", paste(duplicated_group_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  missing_groups <- setdiff(samples, group_samples)
  extra_groups <- setdiff(group_samples, samples)
  if (length(missing_groups) > 0L) {
    stop(
      sprintf("sample_groups is missing samples: %s", paste(missing_groups, collapse = ", ")),
      call. = FALSE
    )
  }
  if (length(extra_groups) > 0L) {
    stop(
      sprintf("sample_groups contains samples not present in activity_scores: %s", paste(extra_groups, collapse = ", ")),
      call. = FALSE
    )
  }

  sample_groups <- sample_groups[samples]
  if (any(is.na(sample_groups) | !nzchar(sample_groups))) {
    stop("Every sample must have a non-empty group.", call. = FALSE)
  }

  if (is.null(names(group_cols)) || any(!nzchar(names(group_cols)))) {
    stop("group_cols must be a named colour vector.", call. = FALSE)
  }
  group_col_names <- names(group_cols)
  group_cols <- stats::setNames(as.character(group_cols), group_col_names)
  missing_group_cols <- setdiff(unique(sample_groups), names(group_cols))
  if (length(missing_group_cols) > 0L) {
    stop(
      sprintf("group_cols is missing colours for groups: %s", paste(missing_group_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  group_order <- unique(c(as.character(group_order), sort(setdiff(unique(sample_groups), as.character(group_order)))))
  group_order <- group_order[group_order %in% unique(sample_groups)]
  metadata <- data.frame(
    sample = samples,
    group = factor(sample_groups, levels = group_order),
    row.names = samples,
    stringsAsFactors = FALSE
  )
  ord <- order(metadata$group)
  ordered_groups <- as.character(metadata$group[ord])
  legend_groups <- group_order[group_order %in% names(group_cols)]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 1))

  for (te in te_subfamilies) {
    graphics::plot(
      seq_along(ord),
      activity_scores[te, ord],
      col = group_cols[ordered_groups],
      pch = pch,
      main = sprintf(main_template, te),
      xlab = xlab,
      ylab = ylab
    )

    if (isTRUE(zero_line)) {
      graphics::abline(h = 0, lty = zero_line_lty, col = zero_line_col)
    }
    graphics::legend(
      legend_position,
      legend_groups,
      col = group_cols[legend_groups],
      pch = pch,
      cex = legend_cex
    )
  }

  invisible(metadata[ord, , drop = FALSE])
}

# Save the main getSignifActivities() result table without expanding diagnostics.
# getSignifActivities() returns a list with model diagnostics (aliased, rss, p,
# N, fstatistic, genes_excluded) plus res$coefs. For routine exports we keep only
# the TE-subfamily coefficient table: estimated activity, standard error, t
# statistic, raw p-value, and adjusted p-value.
signif_activity_table <- function(signif_result,
                                  include_intercept = FALSE) {
  if (is.null(signif_result) || is.null(signif_result$coefs)) {
    stop("signif_result must be a getSignifActivities() result with a coefs table.", call. = FALSE)
  }

  coefs <- as.data.frame(signif_result$coefs, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(rownames(coefs))) {
    stop("signif_result$coefs must have rownames.", call. = FALSE)
  }

  coef_names <- rownames(coefs)
  coef_names <- sub("^`(.*)`$", "\\1", coef_names)
  if (!isTRUE(include_intercept)) {
    keep <- coef_names != "(Intercept)"
    coefs <- coefs[keep, , drop = FALSE]
    coef_names <- coef_names[keep]
  }

  out <- data.frame(
    te_subfamily = coef_names,
    estimate = coefs[["Estimate"]],
    std_error = coefs[["Std. Error"]],
    t_statistic = coefs[["t value"]],
    p_value = coefs[["Pr(>|t|)"]],
    p_adj = coefs[["p_adj"]],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

signif_model_diagnostics <- function(signif_result,
                                     condition = "condition",
                                     print_table = TRUE) {
  # Summarises the model-level diagnostics from getSignifActivities().
  # These diagnostics describe the whole regression model, not individual TEs.
  if (is.null(signif_result)) {
    stop("signif_result must not be NULL.", call. = FALSE)
  }
  if (is.null(signif_result$fstatistic)) {
    stop("signif_result$fstatistic is missing.", call. = FALSE)
  }
  if (is.null(condition) || length(condition) != 1L || !nzchar(condition)) {
    stop("condition must be one non-empty label.", call. = FALSE)
  }

  fstat <- signif_result$fstatistic
  f_value <- unname(fstat[["value"]])
  f_numdf <- unname(fstat[["numdf"]])
  f_dendf <- unname(fstat[["dendf"]])
  f_p_value <- stats::pf(q = f_value, df1 = f_numdf, df2 = f_dendf, lower.tail = FALSE)

  aliased <- signif_result$aliased
  if (is.null(aliased)) {
    aliased_variables <- character()
  } else {
    aliased <- aliased[!is.na(aliased)]
    aliased_variables <- names(aliased)[aliased]
    aliased_variables <- setdiff(aliased_variables, "(Intercept)")
  }

  diagnostics <- data.frame(
    condition = condition,
    f_statistic = as.numeric(f_value),
    f_numdf = as.numeric(f_numdf),
    f_dendf = as.numeric(f_dendf),
    f_p_value = as.numeric(f_p_value),
    r_squared = as.numeric(signif_result$r.squared),
    rss = as.numeric(signif_result$rss),
    n_aliased = length(aliased_variables),
    aliased_variables = if (length(aliased_variables) == 0L) {
      "none"
    } else {
      paste(aliased_variables, collapse = ", ")
    },
    stringsAsFactors = FALSE
  )

  if (isTRUE(print_table)) {
    print(diagnostics, row.names = FALSE)
  }

  diagnostics
}

save_signif_activity_result <- function(input_dir,
                                        condition,
                                        signif_result,
                                        name = NULL,
                                        include_intercept = FALSE,
                                        write_rds = FALSE,
                                        write_tsv = TRUE,
                                        overwrite = TRUE,
                                        cleanup_legacy = TRUE) {
  if (is.null(condition) || length(condition) != 1L || !nzchar(condition)) {
    stop("condition must be one non-empty label.", call. = FALSE)
  }
  if (is.null(name)) {
    name <- paste0(condition, "_signif_act")
  }

  if (isTRUE(cleanup_legacy)) {
    input_dir_expanded <- path.expand(input_dir)
    safe_name <- .safe_result_name(name)
    legacy_suffixes <- c(
      "coefs",
      "aliased",
      "genes_excluded",
      "fstatistic",
      "f_statistic",
      "rss",
      "p",
      "N",
      "n",
      "r",
      "r.squared",
      "ajd.r.squared",
      "adj.r.squared",
      "control_group",
      "treatment_group"
    )
    legacy_paths <- file.path(input_dir_expanded, paste0(safe_name, "_", legacy_suffixes, ".tsv"))
    legacy_paths <- legacy_paths[file.exists(legacy_paths)]
    if (length(legacy_paths) > 0L) {
      unlink(legacy_paths)
    }
  }

  signif_table <- signif_activity_table(signif_result, include_intercept = include_intercept)
  write_result_object(
    input_dir = input_dir,
    name = name,
    object = signif_table,
    write_rds = write_rds,
    write_tsv = write_tsv,
    overwrite = overwrite
  )
}

.write_activity_score_result <- function(input_dir,
                                         name,
                                         activity_scores,
                                         write_rds = FALSE,
                                         write_tsv = TRUE,
                                         overwrite = TRUE,
                                         row_id_col = "te_subfamily") {
  if (is.null(input_dir) || !nzchar(input_dir)) {
    stop("input_dir must be a non-empty directory path", call. = FALSE)
  }

  input_dir <- path.expand(input_dir)
  if (!dir.exists(input_dir)) {
    dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  }

  safe_name <- .safe_result_name(name)
  written <- data.frame(
    name = character(),
    format = character(),
    path = character(),
    stringsAsFactors = FALSE
  )

  if (isTRUE(write_rds)) {
    rds_path <- file.path(input_dir, paste0(safe_name, ".rds"))
    if (file.exists(rds_path) && !isTRUE(overwrite)) {
      stop(sprintf("Refusing to overwrite existing file: %s", rds_path), call. = FALSE)
    }
    saveRDS(activity_scores, rds_path)
    written <- rbind(
      written,
      data.frame(
        name = safe_name,
        format = "rds",
        path = normalizePath(rds_path, mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }

  if (isTRUE(write_tsv)) {
    tsv_path <- file.path(input_dir, paste0(safe_name, ".tsv"))
    written <- rbind(
      written,
      data.frame(
        name = safe_name,
        format = "tsv",
        path = save_activity_score_matrix_tsv(
          activity_scores = activity_scores,
          output_path = tsv_path,
          row_id_col = row_id_col,
          overwrite = overwrite
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  written
}

# Convert an RDS object to TSV when the stored object is tabular.
# Matrices/data.frames/vectors become one TSV. Lists are written as one TSV per
# tabular list element; keep the RDS when exact list reconstruction matters.
rds_to_tsv <- function(rds_path,
                       tsv_path = NULL,
                       overwrite = TRUE) {
  if (!file.exists(rds_path)) {
    stop(sprintf("RDS file not found: %s", rds_path), call. = FALSE)
  }

  object <- readRDS(rds_path)
  object_name <- .result_path_stem(rds_path)
  output_dir <- dirname(rds_path)

  write_tsv_table <- function(table_object, path) {
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop(sprintf("Refusing to overwrite existing file: %s", path), call. = FALSE)
    }
    utils::write.table(
      table_object,
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    normalizePath(path, mustWork = TRUE)
  }

  table_object <- .result_object_to_table(object)
  if (!is.null(table_object)) {
    if (is.null(tsv_path)) {
      tsv_path <- file.path(output_dir, paste0(object_name, ".tsv"))
    }
    return(write_tsv_table(table_object, tsv_path))
  }

  if (!is.list(object)) {
    stop("RDS object cannot be represented as a TSV table.", call. = FALSE)
  }

  if (!is.null(tsv_path) && !dir.exists(tsv_path)) {
    stop("For list objects, tsv_path must be NULL or an existing output directory.", call. = FALSE)
  }
  if (!is.null(tsv_path)) {
    output_dir <- tsv_path
  }

  written <- data.frame(
    name = character(),
    path = character(),
    stringsAsFactors = FALSE
  )

  for (element_name in names(object)) {
    element_table <- .result_object_to_table(object[[element_name]])
    if (is.null(element_table)) {
      next
    }
    file_name <- .safe_result_name(paste(object_name, element_name, sep = "_"))
    path <- file.path(output_dir, paste0(file_name, ".tsv"))
    written <- rbind(
      written,
      data.frame(
        name = file_name,
        path = write_tsv_table(element_table, path),
        stringsAsFactors = FALSE
      )
    )
  }

  if (nrow(written) == 0L) {
    stop("No tabular list elements could be represented as TSV.", call. = FALSE)
  }

  written
}

# Convert a TSV back to an RDS object.
# This reconstructs a data.frame by default. With as_matrix = TRUE, all columns
# except row_id_col are coerced to a matrix and row_id_col becomes rownames.
tsv_to_rds <- function(tsv_path,
                       rds_path = NULL,
                       row_id_col = "row_id",
                       as_matrix = FALSE,
                       overwrite = TRUE) {
  if (!file.exists(tsv_path)) {
    stop(sprintf("TSV file not found: %s", tsv_path), call. = FALSE)
  }
  if (is.null(rds_path)) {
    rds_path <- file.path(dirname(tsv_path), paste0(.result_path_stem(tsv_path), ".rds"))
  }
  if (file.exists(rds_path) && !isTRUE(overwrite)) {
    stop(sprintf("Refusing to overwrite existing file: %s", rds_path), call. = FALSE)
  }

  object <- utils::read.delim(
    tsv_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(as_matrix)) {
    if (!row_id_col %in% colnames(object)) {
      stop(sprintf("row_id_col '%s' is not present in the TSV.", row_id_col), call. = FALSE)
    }
    row_ids <- object[[row_id_col]]
    object[[row_id_col]] <- NULL
    object <- as.matrix(object)
    rownames(object) <- row_ids
  }

  saveRDS(object, rds_path)
  normalizePath(rds_path, mustWork = TRUE)
}

build_activity_lambda_diagnostics <- function(act) {
  # Captures the shared lambda grid used by getActivities(), including the
  # cross-sample optimum used for the final coefficients.
  if (is.null(act$lambdas) && is.null(act$lambda_min_overall)) {
    return(NULL)
  }

  lambda_values <- if (!is.null(act$lambdas)) {
    as.numeric(act$lambdas)
  } else {
    as.numeric(act$lambda_min_overall)[1]
  }
  lambda_values <- lambda_values[!is.na(lambda_values)]
  if (length(lambda_values) == 0L) {
    return(NULL)
  }

  align_lambda_metric <- function(x, metric_name) {
    if (is.null(x)) {
      return(rep(NA_real_, length(lambda_values)))
    }
    x <- as.numeric(x)
    if (length(x) != length(lambda_values)) {
      stop(
        sprintf(
          "act$%s must have length %d to match act$lambdas.",
          metric_name,
          length(lambda_values)
        ),
        call. = FALSE
      )
    }
    x
  }

  lambda_mse <- align_lambda_metric(act$lambda_mse, "lambda_mse")
  lambda_mse_sd <- align_lambda_metric(act$lambda_mse_sd, "lambda_mse_sd")

  lambda_min_overall <- if (!is.null(act$lambda_min_overall)) {
    as.numeric(act$lambda_min_overall)[1]
  } else if (any(!is.na(lambda_mse))) {
    lambda_values[which.min(lambda_mse)]
  } else {
    NA_real_
  }

  global_idx <- NA_integer_
  if (!is.na(lambda_min_overall)) {
    global_idx <- which.min(abs(lambda_values - lambda_min_overall))
  }

  data.frame(
    lambda_index = seq_along(lambda_values),
    lambda = lambda_values,
    mean_crossval_mse = lambda_mse,
    sd_crossval_mse = lambda_mse_sd,
    is_global_lambda_min = seq_along(lambda_values) == global_idx,
    global_lambda_min = lambda_min_overall,
    global_lambda_max = max(lambda_values, na.rm = TRUE),
    global_lambda_grid_min = min(lambda_values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

build_activity_fit_tables <- function(act,
                                      A = NULL,
                                      sample_conditions) {
  # Builds compact sample-level diagnostics and an overall/condition summary.
  if (is.null(A)) {
    if (is.null(act$activities)) {
      stop("A is NULL and act$activities is not available.", call. = FALSE)
    }
    A <- act$activities
  }
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  if (is.null(colnames(A))) {
    stop("Activity matrix A must have sample colnames.", call. = FALSE)
  }

  samples <- colnames(A)
  if (anyDuplicated(samples)) {
    stop("Activity matrix A contains duplicated sample colnames.", call. = FALSE)
  }

  align_metric <- function(x, metric_name, allow_scalar = FALSE) {
    if (is.null(x)) {
      stop(sprintf("act$%s is required.", metric_name), call. = FALSE)
    }
    if (allow_scalar && length(x) == 1L) {
      return(rep(as.numeric(x), length(samples)))
    }

    metric_names <- names(x)
    if (!is.null(metric_names)) {
      if (anyDuplicated(metric_names)) {
        stop(sprintf("act$%s contains duplicated sample names.", metric_name), call. = FALSE)
      }
      missing <- setdiff(samples, metric_names)
      if (length(missing) > 0L) {
        stop(
          sprintf(
            "act$%s is missing samples: %s",
            metric_name,
            paste(missing, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      return(as.numeric(x[samples]))
    }

    if (length(x) != length(samples)) {
      stop(
        sprintf(
          "act$%s must have length %d or names matching A colnames.",
          metric_name,
          length(samples)
        ),
        call. = FALSE
      )
    }

    as.numeric(x)
  }

  if (is.data.frame(sample_conditions)) {
    required_cols <- c("sample_id", "condition")
    missing_cols <- setdiff(required_cols, colnames(sample_conditions))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "sample_conditions data.frame must contain columns: %s",
          paste(required_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    sample_conditions <- stats::setNames(
      as.character(sample_conditions$condition),
      as.character(sample_conditions$sample_id)
    )
  }

  if (is.null(sample_conditions) || is.null(names(sample_conditions))) {
    stop("sample_conditions must be a named vector or a sample_id/condition data.frame.", call. = FALSE)
  }

  condition_samples <- names(sample_conditions)
  if (any(!nzchar(condition_samples))) {
    stop("Every sample_conditions entry must have a non-empty sample name.", call. = FALSE)
  }
  if (anyDuplicated(condition_samples)) {
    duplicated_samples <- unique(condition_samples[duplicated(condition_samples)])
    stop(
      sprintf(
        "sample_conditions assigns more than one condition to: %s",
        paste(duplicated_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_conditions <- setdiff(samples, condition_samples)
  extra_conditions <- setdiff(condition_samples, samples)
  if (length(missing_conditions) > 0L) {
    stop(
      sprintf(
        "sample_conditions is missing samples: %s",
        paste(missing_conditions, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(extra_conditions) > 0L) {
    stop(
      sprintf(
        "sample_conditions contains samples not present in A: %s",
        paste(extra_conditions, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  conditions <- as.character(sample_conditions[samples])
  if (any(is.na(conditions) | !nzchar(conditions))) {
    stop("Every sample must have a non-empty condition.", call. = FALSE)
  }

  sample_activity_diagnostics <- data.frame(
    sample_id = samples,
    condition = conditions,
    r2 = align_metric(act$r2, "r2"),
    r2_max = align_metric(act$r2_max, "r2_max"),
    stringsAsFactors = FALSE
  )
  sample_activity_diagnostics$r2_ratio <- sample_activity_diagnostics$r2 / sample_activity_diagnostics$r2_max

  if (!is.null(act$lambdas_best)) {
    lambda_best <- tryCatch(
      align_metric(act$lambdas_best, "lambdas_best", allow_scalar = TRUE),
      error = function(e) NULL
    )
    if (!is.null(lambda_best)) {
      sample_activity_diagnostics$lambda_best <- lambda_best
    }
  }

  if (!is.null(act$lambda_max_by_sample)) {
    lambda_max_by_sample <- tryCatch(
      align_metric(act$lambda_max_by_sample, "lambda_max_by_sample", allow_scalar = TRUE),
      error = function(e) NULL
    )
    if (!is.null(lambda_max_by_sample)) {
      sample_activity_diagnostics$lambda_max <- lambda_max_by_sample
    }
  }

  if (!is.null(act$lambda_min_overall)) {
    lambda_min_overall <- as.numeric(act$lambda_min_overall)[1]
    sample_activity_diagnostics$lambda_min_overall <- lambda_min_overall
    if ("lambda_best" %in% colnames(sample_activity_diagnostics) && !is.na(lambda_min_overall) && lambda_min_overall != 0) {
      sample_activity_diagnostics$lambda_best_to_global_ratio <- sample_activity_diagnostics$lambda_best / lambda_min_overall
    }
  }

  sample_activity_diagnostics$zero_fraction <- unname(colMeans(A == 0, na.rm = TRUE))

  summary_rows <- function(scope, condition, idx) {
    data.frame(
      scope = scope,
      condition = condition,
      n_samples = length(idx),
      mean_r2 = mean(sample_activity_diagnostics$r2[idx], na.rm = TRUE),
      mean_r2_ratio = mean(sample_activity_diagnostics$r2_ratio[idx], na.rm = TRUE),
      mean_sparsity_percent = 100 * mean(sample_activity_diagnostics$zero_fraction[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  condition_order <- unique(conditions)
  activity_fit_summary <- do.call(
    rbind,
    c(
      list(summary_rows("overall", "all", seq_along(samples))),
      lapply(condition_order, function(condition) {
        summary_rows("condition", condition, which(conditions == condition))
      })
    )
  )
  rownames(activity_fit_summary) <- NULL

  list(
    sample_activity_diagnostics = sample_activity_diagnostics,
    activity_fit_summary = activity_fit_summary,
    activity_lambda_diagnostics = build_activity_lambda_diagnostics(act)
  )
}

# Export the standard activity-analysis objects for the current E/N combination.
# TSV is the default because it is inspectable; RDS is optional for exact R reuse.
save_activity_results <- function(input_dir,
                                  activity_scores,
                                  model_fit,
                                  sparsity,
                                  significance_activities = list(),
                                  target_gene_matrices = list(),
                                  extra = list(),
                                  write_rds = FALSE,
                                  write_tsv = TRUE,
                                  overwrite = TRUE) {
  saved <- list(
    .write_activity_score_result(input_dir, "activity_scores", activity_scores, write_rds = write_rds, write_tsv = write_tsv, overwrite = overwrite),
    write_result_object(input_dir, "model_fit", model_fit, write_rds = write_rds, write_tsv = write_tsv, overwrite = overwrite),
    write_result_object(input_dir, "sparsity", sparsity, write_rds = write_rds, write_tsv = write_tsv, overwrite = overwrite)
  )

  for (name in names(significance_activities)) {
    saved[[length(saved) + 1L]] <- save_signif_activity_result(
      input_dir = input_dir,
      condition = name,
      signif_result = significance_activities[[name]],
      name = paste0(name, "_signif_act"),
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  for (name in names(target_gene_matrices)) {
    saved[[length(saved) + 1L]] <- write_result_object(
      input_dir,
      paste0("target_gene_matrix_", name),
      target_gene_matrices[[name]],
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  for (name in names(extra)) {
    saved[[length(saved) + 1L]] <- write_result_object(
      input_dir,
      name,
      extra[[name]],
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  do.call(rbind, saved)
}

save_crates_activity_results <- function(input_dir,
                                         act,
                                         sample_conditions,
                                         A = NULL,
                                         fit_tables = NULL,
                                         significance_activities = list(),
                                         target_gene_matrices = list(),
                                         extra = list(),
                                         write_rds = FALSE,
                                         write_tsv = TRUE,
                                         overwrite = TRUE) {
  if (is.null(A)) {
    if (is.null(act$activities)) {
      stop("A is NULL and act$activities is not available.", call. = FALSE)
    }
    A <- act$activities
  }
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }

  reserved_extra <- c("sample_activity_diagnostics", "activity_fit_summary", "activity_lambda_diagnostics")
  conflicting_extra <- intersect(reserved_extra, names(extra))
  if (length(conflicting_extra) > 0L) {
    stop(
      sprintf(
        "extra must not contain automatically built result names: %s",
        paste(conflicting_extra, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (is.null(fit_tables)) {
    fit_tables <- build_activity_fit_tables(
      act = act,
      A = A,
      sample_conditions = sample_conditions
    )
  }

  required_fit_tables <- c("sample_activity_diagnostics", "activity_fit_summary")
  missing_fit_tables <- setdiff(required_fit_tables, names(fit_tables))
  if (length(missing_fit_tables) > 0L) {
    stop(
      sprintf(
        "fit_tables is missing required tables: %s",
        paste(missing_fit_tables, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  saved <- list(
    .write_activity_score_result(input_dir, "activity_scores", A, write_rds = write_rds, write_tsv = write_tsv, overwrite = overwrite),
    write_result_object(input_dir, "sample_activity_diagnostics", fit_tables$sample_activity_diagnostics, write_rds = write_rds, write_tsv = write_tsv, overwrite = overwrite),
    write_result_object(input_dir, "activity_fit_summary", fit_tables$activity_fit_summary, write_rds = write_rds, write_tsv = write_tsv, overwrite = overwrite)
  )

  if (!is.null(fit_tables$activity_lambda_diagnostics)) {
    saved[[length(saved) + 1L]] <- write_result_object(
      input_dir,
      "activity_lambda_diagnostics",
      fit_tables$activity_lambda_diagnostics,
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  for (name in names(significance_activities)) {
    saved[[length(saved) + 1L]] <- save_signif_activity_result(
      input_dir = input_dir,
      condition = name,
      signif_result = significance_activities[[name]],
      name = paste0(name, "_signif_act"),
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  for (name in names(target_gene_matrices)) {
    saved[[length(saved) + 1L]] <- write_result_object(
      input_dir,
      paste0("target_gene_matrix_", name),
      target_gene_matrices[[name]],
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  for (name in names(extra)) {
    saved[[length(saved) + 1L]] <- write_result_object(
      input_dir,
      name,
      extra[[name]],
      write_rds = write_rds,
      write_tsv = write_tsv,
      overwrite = overwrite
    )
  }

  do.call(rbind, saved)
}

duplicate_stats <- function(x, data = NULL, sum_thres  = NULL, include_table = TRUE) {
  # print stats by rownames duplicates.
  # Returns a list with summary counts + a per-group data.frame.
  if (!is.null(sum_thres ) && is.null(data)) {
    stop("data must be supplied when sum_thres  is used")
  }
  if (!is.null(data) && nrow(data) != length(x)) {
    stop("data must have the same number of rows as x has identifiers", call. = FALSE)
  }
  if (!is.null(data)) {
    source_cols <- intersect(c("Source", "source", "mapping_source"), colnames(data))
    if (length(source_cols) > 0L) {
      data <- data[, setdiff(colnames(data), source_cols), drop = FALSE]
    }
    if (!is.numeric(data)) {
      numeric_data <- matrix(
        suppressWarnings(as.numeric(data)),
        nrow = nrow(data),
        dimnames = dimnames(data)
      )
      conversion_failed <- is.na(numeric_data) & !is.na(data) & data != ""
      if (any(conversion_failed)) {
        stop("data must contain numeric counts after removing source columns", call. = FALSE)
      }
      data <- numeric_data
    }
  }
  valid <- !is.na(x)
  x_non_na <- x[valid]
  
  tab <- table(x_non_na)
  dup_groups <- tab[tab > 1]
  
  positions <- lapply(split(which(valid), x_non_na), function(idx) { if (length(idx) > 1) idx else NULL })
  positions <- positions[!sapply(positions, is.null)]
  
  if (!is.null(sum_thres )) {
    row_sums <- rowSums(data[valid, , drop = FALSE], na.rm = TRUE)
    positions <- positions[sapply(positions, function(idx) all(row_sums[idx] > sum_thres ))]
  }
  
  dup_groups <- setNames(vapply(positions, length, integer(1)), names(positions))

  summary <- list(
    duplicate_groups = length(dup_groups),
    affected_rows = sum(dup_groups),
    redundant_rows = sum(dup_groups - 1)
  )

  if (!include_table) {
    return(summary)
  }
  
  dup_df <- do.call(rbind, lapply(names(positions), function(id) {
    idx <- positions[[id]]
    
    row_sum <- NA
    row_var <- NA
    
    if (!is.null(data)) {
      sub <- data[idx, , drop = FALSE]
      
      # row-wise stats
      rs <- rowSums(sub, na.rm = TRUE)
      rv <- apply(sub, 1, var, na.rm = TRUE)
      
      row_sum <- paste(rs, collapse = ",")
      row_var <- paste(round(rv, digits= 0), collapse = ",")}
    
    data.frame(
      mapped_id = id,
      count = length(idx),
      positions = paste(idx, collapse = ","),
      row_sum = row_sum,
      row_var = row_var,
      stringsAsFactors = FALSE)}))
  
  c(summary, list(table = dup_df))}

smart_dedup_mapping <- function(map_df,
                                  already_mapped,
                                  N_genes,
                                  biotype_col = NULL,
                                  preferred_biotype = "protein_coding") {
    
    # map_df needs at least:
    # hgnc_symbol, ensembl_gene_id
    # optional: gene_biotype
    
    split_df <- split(map_df, map_df$hgnc_symbol)
    
    result <- sapply(split_df, function(df) {
      
      df <- df[!is.na(df$ensembl_gene_id), , drop = FALSE]
      df <- df[!duplicated(df$ensembl_gene_id), , drop = FALSE]
      
      if (nrow(df) == 0) return(NA_character_)
      if (nrow(df) == 1) return(df$ensembl_gene_id[1])
      
      # Optional: prefer protein_coding
      if (!is.null(biotype_col) && biotype_col %in% colnames(df)) {
        df_pref <- df[df[[biotype_col]] == preferred_biotype, , drop = FALSE]
        if (nrow(df_pref) > 0) df <- df_pref}
      
      ids <- df$ensembl_gene_id
      
      # Priority 1: in N and not yet used
      p1 <- ids[ids %in% N_genes & !ids %in% already_mapped]
      if (length(p1) > 0) return(p1[1])
      
      # Priority 2: in N
      p2 <- ids[ids %in% N_genes]
      if (length(p2) > 0) return(p2[1])
      
      # Priority 3: fallback
      return(ids[1])})
    
  return(result)}  # named vector: symbol -> ensembl_id

smart_pick <- function(ids, already_claimed, N_gene ) {
    #wrapper function for make_smart_pick for smart assignment of Ensembl IDs using N for second mapping
    ids <- unique(ids[!is.na(ids)])
    if (length(ids) == 0) return(NA_character_)
    if (length(ids) == 1) return(ids)
    
    # Priority 1: in N and not yet claimed by biomaRt
    p1 <- ids[ids %in% N_genes & !ids %in% already_claimed]
    if (length(p1) > 0) return(p1[1])
    
    # Priority 2: in N (even if already claimed)
    p2 <- ids[ids %in% N_genes]
    if (length(p2) > 0) return(p2[1])
    
    # Priority 3: fallback
    return(ids[1])}

make_smart_pick <- function(already_claimed, N_genes) {
  # returns a function that can be used in mapIds to smartly pick an Ensembl ID for a given HGNC symbol, prioritizing N genes and avoiding already claimed IDs
  function(ids) {  # mapIds will call this with just ids
    ids <- unique(ids[!is.na(ids)])
    if (length(ids) == 0) return(NA_character_)
    if (length(ids) == 1) return(ids)
    
    # Priority 1: in N and not yet claimed
    p1 <- ids[ids %in% N_genes & !ids %in% already_claimed]
    if (length(p1) > 0) return(p1[1])
    
    # Priority 2: in N (even if already claimed)
    p2 <- ids[ids %in% N_genes]
    if (length(p2) > 0) return(p2[1])
    
    # Priority 3: fallback
    return(ids[1])}}

# ─────────────────────────────────────────────────────────────────
# resolve_duplicates()
#
# method:
#   "sum"     – colSums of all duplicate rows
#   "mean"    – colMeans of all duplicate rows
#   "max_sum" – keep single row with highest rowSum
#   "max_var" – keep single row with highest rowVariance
#   "discard" – remove all duplicate rows entirely, keep singletons only
#   "identical" – collapse duplicate rows only if count profiles are exactly identical
#
# threshold:
#   Numeric. Before aggregating a duplicate group, rows whose
#   rowSum < threshold are dropped. If ALL rows are dropped the
#   gene is removed from the result (set threshold = NULL to keep
#   all rows regardless of expression level).
#   Has no effect when method = "discard".
#
# verbose:
#   Logical. If TRUE, prints a summary of remaining duplicates
#   after threshold filtering, just before aggregation.
#
# source:
#   Optional source reader. If a single column name is supplied and E
#   contains that column (default: "Source"), the column is extracted
#   and removed before the count matrix is processed. A row-wise
#   vector or attr(E, "source") are still accepted as fallbacks.
#   Rows marked "pre_mapped" or "pre-mapped" can be handled specially
#   within duplicate groups.
#
# pre_mapped:
#   "pre-mapped" – prefer pre-mapped rows over translated rows within
#                  duplicate groups
#   "same"       – treat pre-mapped rows like every other mapped row
#   "unmapped"   – discard pre-mapped rows within duplicate groups
#
# Returns a matrix with the same column names as E and one row per
# unique gene.
# ─────────────────────────────────────────────────────────────────
resolve_duplicates <- function(E,
                               method    = c("sum", "mean", "max_sum", "max_var", "discard", "identical"),
                               threshold = NULL,
                               verbose   = TRUE,
                               source    = "Source",
                               pre_mapped = c("pre-mapped", "same", "unmapped", "pre_mapped")) {
  
  method <- match.arg(method)
  pre_mapped <- match.arg(pre_mapped)
  if (pre_mapped == "pre_mapped") pre_mapped <- "pre-mapped"
  stopifnot(is.matrix(E), !is.null(rownames(E)))
  
  source_values <- NULL
  if (!is.null(source) && length(source) == 1L && !is.na(source) && source %in% colnames(E)) {
    source_col <- match(source, colnames(E))
    source_values <- E[, source_col]
    E <- E[, -source_col, drop = FALSE]
  } else if (is.null(source)) {
    source_values <- attr(E, "source", exact = TRUE)
  } else if (length(source) == nrow(E)) {
    source_values <- source
  } else {
    source_values <- attr(E, "source", exact = TRUE)
  }

  if (ncol(E) == 0L) {
    stop("E must contain at least one count column after removing the source column", call. = FALSE)
  }

  if (!is.numeric(E)) {
    count_values <- suppressWarnings(as.numeric(E))
    numeric_E <- matrix(count_values, nrow = nrow(E), dimnames = dimnames(E))
    conversion_failed <- is.na(numeric_E) & !is.na(E) & E != ""
    if (any(conversion_failed)) {
      stop("All non-source columns in E must be numeric counts", call. = FALSE)
    }
    E <- .restore_integer_counts(numeric_E)
  } else {
    E <- .restore_integer_counts(E)
  }
  
  rn <- rownames(E)
  n <- length(rn)

  if (is.null(source_values)) {
    source_values <- rep(NA_character_, n)
  } else {
    source_values <- as.character(source_values)
    if (length(source_values) != n) {
      stop("source must be NULL, a column name in E, or have one entry per row in E", call. = FALSE)
    }
  }
  is_pre_mapped <- source_values %in% c("pre_mapped", "pre-mapped")
  
  # Pre-compute once; reused for threshold filter & max_sum / max_var
  rs <- rowSums(E, na.rm = TRUE)
  rv <- if (method == "max_var") apply(E, 1, var, na.rm = TRUE) else NULL
  
  # Identify duplicated genes, then optionally resolve pre-mapped conflicts
  # before applying the requested aggregation method.
  idx_list <- split(seq_along(rn), rn)
  pre_mapped_adjusted <- 0L
  pre_mapped_removed <- 0L

  if (pre_mapped != "same" && any(is_pre_mapped)) {
    idx_list <- lapply(idx_list, function(idx) {
      if (length(idx) <= 1L) return(idx)

      pre_idx <- idx[is_pre_mapped[idx]]
      if (length(pre_idx) == 0L) return(idx)

      if (pre_mapped == "pre-mapped") {
        pre_mapped_adjusted <<- pre_mapped_adjusted + 1L
        return(pre_idx)
      }

      kept_idx <- idx[!is_pre_mapped[idx]]
      pre_mapped_removed <<- pre_mapped_removed + length(pre_idx)
      if (length(kept_idx) == 0L) NULL else kept_idx
    })
  }

  idx_list <- idx_list[!vapply(idx_list, is.null, logical(1))]
  idx_lengths <- vapply(idx_list, length, integer(1))
  dup_ids  <- names(idx_lengths[idx_lengths > 1L])
  uniq_ids <- names(idx_lengths[idx_lengths == 1L])
  
  # ── Part 1: singletons – copy directly ────────────────────────
  singleton_idx <- unlist(idx_list[uniq_ids], use.names = FALSE)
  singleton_mat <- E[singleton_idx, , drop = FALSE]
  
  # ── Discard: return singletons immediately, skip aggregation ──
  if (method == "discard") {
    if (verbose) {
      if (pre_mapped == "pre-mapped" && pre_mapped_adjusted > 0L) {
        message(sprintf(
          "[resolve_duplicates] pre_mapped = 'pre-mapped': pre-mapped rows preferred in %d duplicate groups.",
          pre_mapped_adjusted
        ))
      }
      if (pre_mapped == "unmapped" && pre_mapped_removed > 0L) {
        message(sprintf(
          "[resolve_duplicates] pre_mapped = 'unmapped': %d pre-mapped duplicate rows discarded before duplicate resolution.",
          pre_mapped_removed
        ))
      }
      message(sprintf(
        "[resolve_duplicates] Method 'discard': %d duplicate groups (%d rows) removed.",
        length(dup_ids), sum(idx_lengths[dup_ids])
      ))
    }
    return(.restore_integer_counts(singleton_mat))
  }

  # ── Build per-gene index list for duplicate groups ───────────────
  idx_list <- idx_list[dup_ids]   # keep only duplicated genes

  # ── Identical: collapse only truly identical duplicate profiles ─
  if (method == "identical") {
    is_identical_group <- vapply(idx_list, function(idx) {
      nrow(unique(E[idx, , drop = FALSE])) == 1L
    }, logical(1))

    identical_ids <- names(is_identical_group)[is_identical_group]
    non_identical_ids <- names(is_identical_group)[!is_identical_group]

    duplicate_mat <- if (length(identical_ids) > 0L) {
      m <- E[vapply(idx_list[identical_ids], `[`, integer(1), 1L), , drop = FALSE]
      rownames(m) <- identical_ids
      m
    } else {
      NULL
    }

    non_identical_mat <- if (length(non_identical_ids) > 0L) {
      E[unlist(idx_list[non_identical_ids], use.names = FALSE), , drop = FALSE]
    } else {
      NULL
    }

    if (verbose) {
      if (pre_mapped == "pre-mapped" && pre_mapped_adjusted > 0L) {
        message(sprintf(
          "[resolve_duplicates] pre_mapped = 'pre-mapped': pre-mapped rows preferred in %d duplicate groups.",
          pre_mapped_adjusted
        ))
      }
      if (pre_mapped == "unmapped" && pre_mapped_removed > 0L) {
        message(sprintf(
          "[resolve_duplicates] pre_mapped = 'unmapped': %d pre-mapped duplicate rows discarded before duplicate resolution.",
          pre_mapped_removed
        ))
      }
      message(sprintf(
        "[resolve_duplicates] Method 'identical': %d duplicate groups collapsed.",
        length(identical_ids)
      ))
      if (length(non_identical_ids) > 0L) {
        message(sprintf(
          "[resolve_duplicates] Method 'identical': %d duplicate groups retained because duplicate rows were not identical.",
          length(non_identical_ids)
        ))
      }
    }

    return(.restore_integer_counts(rbind(singleton_mat, duplicate_mat, non_identical_mat)))
  }
  
  idx_list_filtered <- lapply(idx_list, function(idx) {
    if (!is.null(threshold) && length(idx) > 1L) {
      idx <- idx[rs[idx] >= threshold]
    }
    if (length(idx) == 0L) NULL else idx   # NULL flags gene for removal
  })
  
  # ── Part 3: verbose report after threshold, before aggregation ─
  if (verbose) {
    kept_idx <- unlist(Filter(Negate(is.null), idx_list_filtered))

    if (pre_mapped == "pre-mapped" && pre_mapped_adjusted > 0L) {
      message(sprintf(
        "[resolve_duplicates] pre_mapped = 'pre-mapped': pre-mapped rows preferred in %d duplicate groups.",
        pre_mapped_adjusted
      ))
    }
    if (pre_mapped == "unmapped" && pre_mapped_removed > 0L) {
      message(sprintf(
        "[resolve_duplicates] pre_mapped = 'unmapped': %d pre-mapped duplicate rows discarded before duplicate resolution.",
        pre_mapped_removed
      ))
    }
    
    if (!is.null(threshold)) {
      n_removed <- sum(sapply(idx_list_filtered, is.null))
      message(sprintf(
        "[resolve_duplicates] After threshold (>= %g): %d duplicate groups removed.",
        threshold, n_removed
      ))
    }
    
    if (length(kept_idx) > 0) {
      stats <- duplicate_stats(rn[kept_idx], data = E[kept_idx, , drop = FALSE])
      message(sprintf(
        "[resolve_duplicates] Remaining duplicates: %d groups, %d affected rows, %d redundant rows.",
        stats$duplicate_groups, stats$affected_rows, stats$redundant_rows
      ))
      print(stats$table)
    }
  }
  
  # ── Part 4: duplicates – aggregate ────────────────────────────
  agg_rows <- lapply(dup_ids, function(id) {
    
    idx <- idx_list_filtered[[id]]
    if (is.null(idx)) return(NULL)   # gene was filtered out by threshold
    
    sub <- E[idx, , drop = FALSE]
    
    agg <- switch(method,
                  sum     = colSums(sub,  na.rm = TRUE),
                  mean    = colMeans(sub, na.rm = TRUE),
                  max_sum = sub[which.max(rs[idx]), , drop = TRUE],
                  max_var = sub[which.max(rv[idx]), , drop = TRUE]
    )
    
    matrix(agg, nrow = 1L, dimnames = list(id, colnames(E)))
  })
  
  # Drop NULLs (genes filtered out by threshold)
  agg_rows <- Filter(Negate(is.null), agg_rows)
  
  # ── Combine ────────────────────────────────────────────────────
  if (length(agg_rows) == 0L) return(.restore_integer_counts(singleton_mat))
  
  agg_mat <- do.call(rbind, agg_rows)
  .restore_integer_counts(rbind(singleton_mat, agg_mat))
}

# ─────────────────────────────────────────────────────────────────
# apply_mapping_to_expression()
#
# Applies an HGNC/symbol -> Ensembl mapping to an expression count
# matrix and optionally carries row-level mapping-source metadata
# forward into a Source column.
#
# E:
#   Expression count matrix with one row per original gene identifier.
#
# gene_ids:
#   Original row identifiers used to look up rows in mapping.
#   Defaults to rownames(E).
#
# mapping:
#   Either a named vector (names = original identifiers, values =
#   Ensembl IDs) or a data.frame/data.table loaded from a mapping file.
#   Mapping tables are expected to contain gene_col and id_col.
#
# source_col:
#   Candidate column names for mapping-source metadata. A Source column
#   is appended to the returned matrix only if one of these columns exists
#   and at least one kept row is marked pre-mapped.
#
# Returns:
#   E restricted to mapped rows, with rownames set to Ensembl IDs. If the
#   mapping source metadata contains pre-mapped rows, the returned matrix
#   also contains output_source_col (default: "Source").
# ─────────────────────────────────────────────────────────────────
apply_mapping_to_expression <- function(E,
                                        mapping,
                                        gene_ids = rownames(E),
                                        gene_col = "hgnc_symbol",
                                        id_col = "ensembl_gene_id",
                                        source_col = c("Source", "source", "mapping_source"),
                                        output_source_col = "Source",
                                        pre_mapped_values = c("pre-mapped", "pre_mapped"),
                                        verbose = TRUE) {
  stopifnot(is.matrix(E), !is.null(gene_ids))

  if (length(gene_ids) != nrow(E)) {
    stop("gene_ids must have one entry per row in E", call. = FALSE)
  }

  source_values <- NULL

  if (is.vector(mapping) && is.null(dim(mapping))) {
    if (is.null(names(mapping))) {
      stop("mapping vector must be named by the original gene identifiers", call. = FALSE)
    }
    mapped_ids <- unname(mapping[gene_ids])
  } else if (is.data.frame(mapping)) {
    if (!gene_col %in% colnames(mapping)) {
      stop(sprintf("mapping table is missing gene_col '%s'", gene_col), call. = FALSE)
    }
    if (!id_col %in% colnames(mapping)) {
      stop(sprintf("mapping table is missing id_col '%s'", id_col), call. = FALSE)
    }

    mapping_key <- as.character(mapping[[gene_col]])
    matched_rows <- match(gene_ids, mapping_key)
    mapped_ids <- as.character(mapping[[id_col]])[matched_rows]

    source_name <- intersect(source_col, colnames(mapping))[1]
    if (!is.na(source_name)) {
      source_values <- as.character(mapping[[source_name]])[matched_rows]
    }
  } else {
    stop("mapping must be a named vector or a data.frame/data.table", call. = FALSE)
  }

  keep <- !is.na(mapped_ids) & mapped_ids != ""
  E_mapped <- E[keep, , drop = FALSE]
  rownames(E_mapped) <- mapped_ids[keep]

  if (!is.null(source_values)) {
    kept_source <- source_values[keep]
    has_pre_mapped <- any(kept_source %in% pre_mapped_values, na.rm = TRUE)

    if (has_pre_mapped) {
      E_mapped <- cbind(E_mapped, kept_source)
      colnames(E_mapped)[ncol(E_mapped)] <- output_source_col

      if (verbose) {
        message(sprintf(
          "[apply_mapping_to_expression] Added '%s' column for %d pre-mapped rows.",
          output_source_col,
          sum(kept_source %in% pre_mapped_values, na.rm = TRUE)
        ))
      }
    } else if (verbose) {
      message("[apply_mapping_to_expression] Mapping source column found, but no pre-mapped rows were present.")
    }
  }

  E_mapped
}

# ─────────────────────────────────────────────────────────────────
# apply_cached_mapping_to_expression()
#
# Loads a cached mapping TSV from a file or mapping directory and
# applies it to an expression matrix. This keeps notebooks short when
# switching between multiple E datasets that each have their own cache.
# ─────────────────────────────────────────────────────────────────
load_mapping_cache <- function(mapping_path) {
  if (dir.exists(mapping_path)) {
    candidates <- list.files(
      mapping_path,
      pattern = "^map_.*\\.tsv$",
      full.names = TRUE
    )
    if (length(candidates) != 1L) {
      stop(
        sprintf(
          "Expected exactly one mapping TSV in '%s', found %d.",
          mapping_path,
          length(candidates)
        ),
        call. = FALSE
      )
    }
    mapping_path <- candidates[[1]]
  }

  if (!file.exists(mapping_path)) {
    stop(sprintf("Mapping cache not found: %s", mapping_path), call. = FALSE)
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(mapping_path, sep = "\t")
  } else {
    read.delim(mapping_path, sep = "\t", stringsAsFactors = FALSE)
  }
}

apply_cached_mapping_to_expression <- function(E,
                                               mapping_path,
                                               ...,
                                               verbose = TRUE) {
  mapping_cache <- load_mapping_cache(mapping_path)
  apply_mapping_to_expression(
    E,
    mapping = mapping_cache,
    ...,
    verbose = verbose
  )
}

# ─────────────────────────────────────────────────────────────────
# Duplicate-resolution row-sum comparison helpers
#
# resolve_duplicate_chain() applies several duplicate-resolution
# methods in sequence, e.g. c("identical", "sum").
# compare_rowsum_matrices() compares any named list of matrices in the
# row-sum distribution plots used by the duplicate analysis notebook.
# ─────────────────────────────────────────────────────────────────
strip_source_columns <- function(mat,
                                 source_cols = c("Source", "source", "mapping_source")) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }

  drop_cols <- intersect(source_cols, colnames(mat))
  if (length(drop_cols) > 0L) {
    mat <- mat[, setdiff(colnames(mat), drop_cols), drop = FALSE]
  }

  if (ncol(mat) == 0L) {
    stop("Matrix has no count columns after removing source columns.", call. = FALSE)
  }

  if (!is.numeric(mat)) {
    numeric_mat <- matrix(
      suppressWarnings(as.numeric(mat)),
      nrow = nrow(mat),
      dimnames = dimnames(mat)
    )
    conversion_failed <- is.na(numeric_mat) & !is.na(mat) & mat != ""
    if (any(conversion_failed)) {
      stop("Matrix contains non-numeric values after removing source columns.", call. = FALSE)
    }
    mat <- .restore_integer_counts(numeric_mat)
  } else {
    mat <- .restore_integer_counts(mat)
  }

  mat
}

resolve_duplicate_chain <- function(E,
                                    methods,
                                    ...,
                                    verbose = FALSE) {
  out <- E
  for (method in methods) {
    out <- resolve_duplicates(out, method = method, ..., verbose = verbose)
  }
  out
}

non_identical_duplicate_genes <- function(E,
                                          source_cols = c("Source", "source", "mapping_source")) {
  E_counts <- strip_source_columns(E, source_cols = source_cols)
  idx_list <- split(seq_len(nrow(E_counts)), rownames(E_counts))
  idx_list <- idx_list[vapply(idx_list, length, integer(1)) > 1L]

  names(idx_list)[vapply(idx_list, function(idx) {
    nrow(unique(E_counts[idx, , drop = FALSE])) > 1L
  }, logical(1))]
}

compare_rowsum_matrices <- function(matrix_list,
                                    reference = names(matrix_list)[[1]],
                                    highlight_genes = character(),
                                    prob_grid = seq(0.001, 0.999, length.out = 500),
                                    method_cols = NULL,
                                    highlight_label = "highlighted genes",
                                    other_label = "other genes",
                                    highlight_cols = NULL) {
  if (!is.list(matrix_list) || length(matrix_list) == 0L) {
    stop("matrix_list must be a non-empty named list of matrices.", call. = FALSE)
  }
  if (is.null(names(matrix_list)) || any(names(matrix_list) == "")) {
    stop("matrix_list must be named.", call. = FALSE)
  }
  if (!reference %in% names(matrix_list)) {
    stop(sprintf("Reference '%s' is not present in matrix_list.", reference), call. = FALSE)
  }

  method_names <- names(matrix_list)
  highlight_genes <- unique(as.character(unlist(highlight_genes)))
  if (length(highlight_label) != 1L || is.na(highlight_label) || highlight_label == "") {
    stop("highlight_label must be a non-empty character scalar.", call. = FALSE)
  }
  if (length(other_label) != 1L || is.na(other_label) || other_label == "") {
    stop("other_label must be a non-empty character scalar.", call. = FALSE)
  }

  if (is.null(highlight_cols)) {
    highlight_cols <- c("grey72", "#d7191c")
    names(highlight_cols) <- c(other_label, highlight_label)
  } else if (is.null(names(highlight_cols))) {
    if (length(highlight_cols) != 2L) {
      stop("Unnamed highlight_cols must contain exactly two colors: other, highlighted.", call. = FALSE)
    }
    names(highlight_cols) <- c(other_label, highlight_label)
  } else if (!all(c(other_label, highlight_label) %in% names(highlight_cols))) {
    old_names <- c("other genes", "highlighted genes")
    if (all(old_names %in% names(highlight_cols))) {
      highlight_cols <- stats::setNames(
        unname(highlight_cols[old_names]),
        c(other_label, highlight_label)
      )
    } else {
      stop("highlight_cols must contain colors named by other_label and highlight_label.", call. = FALSE)
    }
  }

  make_rowsum_df <- function(mat, plot_label) {
    mat_numeric <- strip_source_columns(mat)
    row_sum <- rowSums(mat_numeric, na.rm = TRUE)

    rowsum_df <- data.frame(
      row_index = seq_len(nrow(mat_numeric)),
      gene_id = rownames(mat_numeric),
      row_sum = row_sum,
      row_sum_log10 = log10(row_sum + 1),
      plot_label = plot_label,
      stringsAsFactors = FALSE
    )

    rowsum_df$is_highlighted <- rowsum_df$gene_id %in% highlight_genes
    rowsum_df$highlight <- ifelse(rowsum_df$is_highlighted, highlight_label, other_label)
    rowsum_df <- rowsum_df[order(rowsum_df$row_sum, rowsum_df$gene_id, rowsum_df$row_index), ]
    rowsum_df$rank <- seq_len(nrow(rowsum_df))
    rowsum_df$rank_percent <- 100 * (rowsum_df$rank - 0.5) / nrow(rowsum_df)
    rowsum_df
  }

  rowsum_comparison <- do.call(
    rbind,
    Map(make_rowsum_df, matrix_list, method_names)
  )
  rowsum_comparison$plot_label <- factor(rowsum_comparison$plot_label, levels = method_names)

  if (is.null(method_cols)) {
    base_cols <- c("#111111", "#2b8cbe", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e")
    method_cols <- if (length(method_names) <= length(base_cols)) {
      base_cols[seq_along(method_names)]
    } else {
      grDevices::hcl.colors(length(method_names), palette = "Dark 3")
    }
    names(method_cols) <- method_names
    method_cols[reference] <- "#111111"
  }

  ref_values <- rowsum_comparison$row_sum_log10[rowsum_comparison$plot_label == reference]
  ref_quantiles <- quantile(ref_values, probs = prob_grid, na.rm = TRUE, names = FALSE, type = 8)
  comparison_names <- setdiff(method_names, reference)

  quantile_comparison <- if (length(comparison_names) > 0L) {
    do.call(
      rbind,
      lapply(comparison_names, function(method_label) {
        method_values <- rowsum_comparison$row_sum_log10[rowsum_comparison$plot_label == method_label]
        method_quantiles <- quantile(method_values, probs = prob_grid, na.rm = TRUE, names = FALSE, type = 8)

        data.frame(
          method = method_label,
          prob = prob_grid,
          reference_quantile = ref_quantiles,
          method_quantile = method_quantiles,
          diff = method_quantiles - ref_quantiles,
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    data.frame(
      method = character(),
      prob = numeric(),
      reference_quantile = numeric(),
      method_quantile = numeric(),
      diff = numeric(),
      stringsAsFactors = FALSE
    )
  }

  quantile_comparison$abs_diff <- abs(quantile_comparison$diff)

  distance_summary <- if (nrow(quantile_comparison) > 0L) {
    do.call(
      rbind,
      lapply(split(quantile_comparison, quantile_comparison$method), function(df) {
        method_values <- rowsum_comparison$row_sum_log10[rowsum_comparison$plot_label == df$method[1]]

        data.frame(
          method = df$method[1],
          mean_abs_quantile_diff = mean(df$abs_diff, na.rm = TRUE),
          median_abs_quantile_diff = median(df$abs_diff, na.rm = TRUE),
          rmse_quantile_diff = sqrt(mean(df$diff^2, na.rm = TRUE)),
          max_abs_quantile_diff = max(df$abs_diff, na.rm = TRUE),
          ks_statistic = unname(suppressWarnings(stats::ks.test(ref_values, method_values)$statistic)),
          median_signed_diff = median(df$diff, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    data.frame()
  }

  if (nrow(distance_summary) > 0L) {
    distance_summary <- distance_summary[order(distance_summary$mean_abs_quantile_diff), ]
  }

  p_distance <- if (nrow(distance_summary) > 0L) {
    distance_plot_df <- distance_summary
    distance_plot_df$method <- factor(
      distance_plot_df$method,
      levels = rev(distance_plot_df$method)
    )
    distance_plot_df$is_closest <- distance_plot_df$method == distance_summary$method[[1]]

    ggplot2::ggplot(
      distance_plot_df,
      ggplot2::aes(x = method, y = mean_abs_quantile_diff, fill = is_closest)
    ) +
      ggplot2::geom_col(show.legend = FALSE, width = 0.72) +
      ggplot2::geom_text(
        ggplot2::aes(label = round(mean_abs_quantile_diff, 4)),
        hjust = -0.08,
        size = 3.2
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "grey70")) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
      ggplot2::labs(
        title = "Distance to reference distribution",
        subtitle = sprintf(
          "Lower mean absolute quantile difference is closer to '%s'. Closest: %s.",
          reference,
          distance_summary$method[[1]]
        ),
        x = NULL,
        y = "Mean absolute quantile difference"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  } else {
    NULL
  }

  p_density <- ggplot2::ggplot(rowsum_comparison, ggplot2::aes(x = row_sum_log10, colour = plot_label)) +
    ggplot2::geom_density(linewidth = 1) +
    ggplot2::scale_colour_manual(values = method_cols) +
    ggplot2::labs(
      title = "Row-sum distribution after duplicate handling",
      subtitle = sprintf("Curves are compared against '%s'.", reference),
      x = "log10(row sum + 1)",
      y = "Density",
      colour = "Matrix"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")

  p_ecdf <- ggplot2::ggplot(rowsum_comparison, ggplot2::aes(x = row_sum_log10, colour = plot_label)) +
    ggplot2::stat_ecdf(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = method_cols) +
    ggplot2::labs(
      title = "Cumulative row-sum distribution",
      subtitle = sprintf("Closer ECDFs to '%s' indicate a more similar row-sum distribution.", reference),
      x = "log10(row sum + 1)",
      y = "Cumulative fraction of genes",
      colour = "Matrix"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")

  p_qq <- if (nrow(quantile_comparison) > 0L) {
    ggplot2::ggplot(quantile_comparison, ggplot2::aes(x = reference_quantile, y = method_quantile, colour = method)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~ method, ncol = min(3L, length(comparison_names))) +
      ggplot2::coord_equal() +
      ggplot2::scale_colour_manual(values = method_cols[comparison_names]) +
      ggplot2::labs(
        title = sprintf("Q-Q comparison against '%s'", reference),
        subtitle = "Closer to the dashed diagonal means closer to the reference row-sum distribution.",
        x = sprintf("%s quantile: log10(row sum + 1)", reference),
        y = "Method quantile: log10(row sum + 1)",
        colour = "Method"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")
  } else {
    NULL
  }

  position_df <- rowsum_comparison[order(rowsum_comparison$is_highlighted), ]
  p_highlight_position <- ggplot2::ggplot(position_df, ggplot2::aes(x = rank_percent, y = row_sum_log10, colour = highlight)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.4) +
    ggplot2::facet_wrap(~ plot_label, ncol = 1, scales = "free_x") +
    ggplot2::scale_colour_manual(values = highlight_cols) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      title = sprintf("Position of %s within each row-sum distribution", highlight_label),
      subtitle = sprintf("%s are plotted last so they remain visible.", highlight_label),
      x = "Genes ordered by row-sum percentile within matrix",
      y = "log10(row sum + 1)",
      colour = "Gene set"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  list(
    matrices = matrix_list,
    rowsum_comparison = rowsum_comparison,
    quantile_comparison = quantile_comparison,
    distance_summary = distance_summary,
    plots = list(
      distance = p_distance,
      density = p_density,
      ecdf = p_ecdf,
      qq = p_qq,
      highlight_position = p_highlight_position
    )
  )
}

print_rowsum_comparison <- function(comparison) {
  print(comparison$distance_summary)
  if (!is.null(comparison$plots$distance)) {
    print(comparison$plots$distance)
  }
  print(comparison$plots$density)
  print(comparison$plots$ecdf)
  if (!is.null(comparison$plots$qq)) {
    print(comparison$plots$qq)
  }
  print(comparison$plots$highlight_position)
  invisible(comparison)
}






transform_excel_dates <- function(genes) {
  # Regex to find the "Number-Month" pattern for ALL months
  pattern <- "^([0-9]+)-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"
  matches <- grep(pattern, genes, value = TRUE, ignore.case = TRUE)
  
  if(length(matches) > 0) {
    message("Found Excel date-mangled genes: ", paste(matches, collapse = ", "))
  }
  
  # 1. Legacy HGNC Cases (Pre-2020 Nomenclature)
  # "6-Mar" becomes "MARCH6" (NOT MARCHF6)
  genes <- gsub("^([0-9]+)-Mar$", "MARCH\\1", genes, ignore.case = TRUE)
  
  # "9-Sep" becomes "SEPT9" (NOT SEPTIN9)
  genes <- gsub("^([0-9]+)-Sep$", "SEPT\\1", genes, ignore.case = TRUE)
  
  # 2. Standard Reversals for the rest (e.g., "1-Dec" -> "DEC1")
  genes <- gsub("^([0-9]+)-Jan$", "JAN\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Feb$", "FEB\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Apr$", "APR\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-May$", "MAY\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Jun$", "JUN\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Jul$", "JUL\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Aug$", "AUG\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Oct$", "OCT\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Nov$", "NOV\\1", genes, ignore.case = TRUE)
  genes <- gsub("^([0-9]+)-Dec$", "DEC\\1", genes, ignore.case = TRUE)
  
  return(genes)
}

display_intersect = function(
    genes_E,
    genes_N,
    name_E = deparse(substitute(genes_E)),
    name_N = deparse(substitute(genes_N)),
    title = "Gene overlap",
    plot = TRUE) {
  # This function is comparing the interset of two arrays. 
  # we do not considered duplicated named genes here. 
  common_genes   <- intersect(genes_E, genes_N)
  union_genes <- union(genes_E, genes_N)
  only_E         <- setdiff(genes_E, genes_N)
  only_N         <- setdiff(genes_N, genes_E)
  
  # summary: 
  cat(sprintf("%s genes: %d\n", name_E, length(genes_E)))
  cat(sprintf("%s genes: %d\n", name_N, length(genes_N)))
  cat(sprintf("amount of genes we keep in %s: %d (%.2f%% of %s)\n", name_E, length(common_genes), 100 * length(common_genes) / length(genes_E), name_E))
  cat(sprintf("amount of genes we keep in %s: %d (%.2f%% of %s)\n", name_N, length(common_genes), 100 * length(common_genes) / length(genes_N), name_N))
  cat(sprintf("Common genes: %d (%.2f%% of union)\n", length(common_genes), 100 * length(common_genes) / length(union_genes)))
  cat(sprintf("Only in %s: %d (%.2f%% of %s, %.2f%% of union)\n", name_E, length(only_E), 100 * length(only_E) / length(genes_E), name_E, 100 * length(only_E) / length(union_genes)))
  cat(sprintf("Only in %s: %d (%.2f%% of %s, %.2f%% of union)\n", name_N, length(only_N),  100 * length(only_N) / length(genes_N), name_N, 100 * length(only_N) / length(union_genes)))
  
  if (plot) {
    circle_radius <- 0.49
    circle_offset <- 0.35

    circle_E <- data.frame(
      x = circle_radius * cos(seq(0, 2*pi, length.out = 300)) - circle_offset,
      y = circle_radius * sin(seq(0, 2*pi, length.out = 300))
    )
    
    circle_N <- data.frame(
      x = circle_radius * cos(seq(0, 2*pi, length.out = 300)) + circle_offset,
      y = circle_radius * sin(seq(0, 2*pi, length.out = 300))
    )
    
    label_data <- data.frame(
      x = c(-0.62, 0, 0.62),
      y = c(-0.1, -0.1, -0.1),
      label = c(
        sprintf(
          "<b>Only %s</b><br>%d<br>%.2f%% of %s<br>%.2f%% of union",
          name_E,
          length(only_E),
          100 * length(only_E) / length(genes_E),
          name_E,
          100 * length(only_E) / length(union_genes)
        ),
        sprintf(
          "<b>Intersection</b><br>%d<br>%.2f%% of %s<br>%.2f%% of %s<br>%.2f%% of union",
          length(common_genes),
          100 * length(common_genes) / length(genes_E),
          name_E,
          100 * length(common_genes) / length(genes_N),
          name_N,
          100 * length(common_genes) / length(union_genes)
        ),
        sprintf(
          "<b>Only %s</b><br>%d<br>%.2f%% of %s<br>%.2f%% of union",
          name_N,
          length(only_N),
          100 * length(only_N) / length(genes_N),
          name_N,
          100 * length(only_N) / length(union_genes)
        )
      )
    )

    heading_data <- data.frame(
      x = c(-circle_offset, circle_offset),
      y = c(0.31, 0.31),
      label = c(
        sprintf("<b>%s</b><br><span style='font-size:10pt'>%d</span>", name_E, length(genes_E)),
        sprintf("<b>%s</b><br><span style='font-size:10pt'>%d</span>", name_N, length(genes_N))
      )
    )
    
    print(
      ggplot() +
        geom_polygon(data = circle_E, aes(x, y), fill = "steelblue", alpha = 0.35, color = "black") +
        geom_polygon(data = circle_N, aes(x, y), fill = "orange", alpha = 0.35, color = "black") +
        geom_richtext(
          data = label_data,
          aes(x, y, label = label),
          size = 3.3,
          lineheight = 0.9,
          hjust = 0.5,
          label.color = NA,
          fill = NA
        ) +
        geom_richtext(
          data = heading_data,
          aes(x, y, label = label),
          size = 4.3,
          lineheight = 0.95,
          hjust = 0.5,
          label.color = NA,
          fill = "white",
          label.padding = grid::unit(c(2, 4, 2, 4), "pt"),
          label.r = grid::unit(2, "pt")
        ) +
        labs(title = title) +
        coord_equal(xlim = c(-0.84, 0.84), ylim = c(-0.51, 0.51), expand = FALSE, clip = "off") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          plot.margin = margin(3, 3, 3, 3)
        ))}}

check_npm_threshold <- function(N, threshold = 150, plot = TRUE) {
  if (!is.matrix(N) || !is.numeric(N)) {
    N <- as_count_association_matrix(N, object_name = "N")
  }
  
  # column sums = total weighted read count (p_npm) per subfamily
  col_sums <- colSums(N, na.rm = TRUE)
  
  n_total   <- length(col_sums)
  n_below   <- sum(col_sums < threshold)
  n_above   <- sum(col_sums >= threshold)
  pct_below <- 100 * n_below / n_total
  pct_above <- 100 * n_above / n_total
  
  cat(sprintf(
    "\np_npm threshold: %g\n%d / %d subfamilies (%.1f%%) are discarded\n%d / %d subfamilies (%.1f%%) are retained\n",
    threshold,
    n_below, n_total, pct_below,
    n_above, n_total, pct_above
  ))
  
  # names of subfamilies that would be filtered
  filtered <- names(col_sums[col_sums < threshold])
  if (length(filtered) > 0) {
    cat("\nSubfamilies below threshold (filtered):\n")
    print(filtered)
  }
  
  if (plot) {
    df <- data.frame(
      subfamily = names(col_sums),
      p_npm     = as.numeric(col_sums),
      status    = ifelse(col_sums < threshold, "filtered", "kept")
    )
    df <- df[order(df$p_npm, decreasing = TRUE), ]
    df$rank <- seq_len(nrow(df))
    
    p <- ggplot(df, aes(x = rank, y = p_npm, colour = status)) +
      geom_point(alpha = 0.6, size = 0.8) +
      geom_hline(yintercept = threshold, linetype = "dashed", colour = "firebrick", linewidth = 0.6) +
      annotate("text", x = Inf, y = threshold,
               hjust = 1.1, vjust = -0.5,
               label  = sprintf("threshold = %g", threshold),
               colour = "firebrick", size = 3) +
      scale_colour_manual(values = c(filtered = "grey60", kept = "steelblue")) +
      scale_y_continuous(trans = "log2", labels = scales::comma) +
      labs(
        x        = "TE subfamilies ranked by count",
        y        = "Counts per subfamily log 2 scaled",
        colour   = NULL,
        title    = "Filtered TE subfamilies",
        subtitle = sprintf(
          "%d / %d = %.1f%% subfamilies are filtered at threshold (%g)",
          n_below, n_total, pct_below, threshold)
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
    
    print(p)
  }
  
  invisible(list(
    col_sums = col_sums,
    filtered = filtered,
    n_below  = n_below,
    n_above  = n_above,
    n_total  = n_total
  ))
}


#plot the biotypes that were extracted through the mapping: 
plot_mapping_distribution <- function(map_df,
                                      category_col,
                                      highlight = NULL,
                                      heading = "Mapping category distribution",
                                      x_label = NULL,
                                      n_thres = 1) {
  library(ggplot2)
  library(dplyr)
  
  stopifnot(category_col %in% colnames(map_df))
  
  plot_df <- map_df %>%
    filter(!is.na(.data[[category_col]]), .data[[category_col]] != "") %>%
    count(category = .data[[category_col]], name = "count") %>%
    # Filter based on the minimum occurrence threshold
    filter(count >= n_thres) %>% 
    mutate(
      pct    = count / sum(count) * 100,
      label  = paste0(count, "\n(", round(pct, 1), "%)"),
      is_hit = if (is.null(highlight)) FALSE else category == highlight
    ) %>%
    arrange(desc(count))
  
  ordered_levels <- plot_df$category[order(plot_df$count)]
  axis_colours <- if (!is.null(highlight)) {
    ifelse(ordered_levels == highlight, "#185FA5", "grey30")
  } else {
    rep("grey30", length(ordered_levels))
  }
  
  ggplot(plot_df, aes(x = reorder(category, count), y = count, fill = is_hit)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_text(aes(label = label), hjust = -0.05, size = 2.8, lineheight = 0.9) +
    scale_fill_manual(
      values = c("TRUE" = "#185FA5", "FALSE" = "#888780"),
      guide = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    coord_flip() +
    labs(
      title = heading,
      subtitle = paste0("n = ", sum(plot_df$count), " total mappings"),
      x = x_label,
      y = "# genes in mappings"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(colour = axis_colours)
    )
}

# for third mapping: 
# Extraction of Ensembl IDs
extract_ensembl <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_character_)
  if (is.data.frame(x) || is.list(x)) {
    if ("gene" %in% names(x)) return(as.character(x$gene[1]))
    return(NA_character_)
  }
  as.character(x[1])
}


# this function is from the original craTEs paper: 
# We have to extend the matrix just like when the estimation was made. 

# thus e_hat is of dimension n_replicates * n_genes, but actually contains three times the same vector. We can simply rep e_hat by the number of replicates.
# generating the difference in expression vector

getTargetMatrix <- function(res, preprocessed, auto_refit_on_mismatch = TRUE) {
  # preprocessing
  G = nrow(preprocessed$E_centered) # number of genes
  
  M = ncol(preprocessed$N) # number of TE subfams

  if (!is.null(res$p) && !is.na(res$p) && as.integer(res$p) != M) {
    mismatch_message <- sprintf(
      paste(
        "getTargetMatrix() received a getSignifActivities() result fitted with %d TE predictors,",
        "but preprocessed$N has %d TE predictors."
      ),
      as.integer(res$p),
      M
    )

    if (!isTRUE(auto_refit_on_mismatch)) {
      stop(
        paste(
          mismatch_message,
          "This usually means the significance result belongs to an older N matrix.",
          "Rerun getSignifActivities() for the active E_centered/N_filtered pair,",
          "then rebuild preprocessed <- list(E_centered = E_centered, N = N_filtered)."
        ),
        call. = FALSE
      )
    }

    if (!exists("getSignifActivities", mode = "function")) {
      stop(
        paste(
          mismatch_message,
          "Could not auto-refit because getSignifActivities() is not available.",
          "Load craTEs/source the function, rerun getSignifActivities(), then call getTargetMatrix() again."
        ),
        call. = FALSE
      )
    }

    message(paste(
      mismatch_message,
      "Auto-refitting getSignifActivities() with the active preprocessed E/N matrices."
    ))

    res <- getSignifActivities(
      preprocessed$E_centered,
      preprocessed$N,
      treatment_group = res$treatment_group,
      control_group = res$control_group,
      exclude_genes = res$genes_excluded
    )
  }
  
  control = c(sapply(res$control_group, (function(x) preprocessed$E_centered[, x])))
  
  treatment = c(sapply(res$treatment_group, (function(x) preprocessed$E_centered[, x])))
  
  delta = treatment-control # the observations, true Y
  
  S = length(delta) / G # number of treatment/control contrasts after recycling
  stopifnot(S == as.integer(S))
  S = as.integer(S)
  
  intercept = res$coefs[1, 'Estimate']                                       
  
  coef_table <- as.data.frame(res$coefs, check.names = FALSE, stringsAsFactors = FALSE)
  coef_names <- rownames(coef_table)
  coef_names_unquoted <- sub("^`(.*)`$", "\\1", coef_names)
  fitted_coef_names <- coef_names[coef_names != "(Intercept)"]
  fitted_estimates <- coef_table[fitted_coef_names, "Estimate"]
  names(fitted_estimates) <- coef_names_unquoted[coef_names != "(Intercept)"]

  A_hat <- setNames(rep(0, M), colnames(preprocessed$N)) # the fitted activities
  coefficient_was_estimated <- setNames(rep(FALSE, M), colnames(preprocessed$N))

  direct_idx <- match(names(fitted_estimates), names(A_hat))
  direct_hit <- !is.na(direct_idx)
  if (any(direct_hit)) {
    A_hat[direct_idx[direct_hit]] <- fitted_estimates[direct_hit]
    coefficient_was_estimated[direct_idx[direct_hit]] <- TRUE
  }

  remaining <- !direct_hit
  if (any(remaining)) {
    syntactic_n_names <- make.names(colnames(preprocessed$N), unique = TRUE)
    syntactic_idx <- match(names(fitted_estimates)[remaining], syntactic_n_names)
    syntactic_hit <- !is.na(syntactic_idx)
    if (any(syntactic_hit)) {
      remaining_estimates <- fitted_estimates[remaining]
      A_hat[syntactic_idx[syntactic_hit]] <- remaining_estimates[syntactic_hit]
      coefficient_was_estimated[syntactic_idx[syntactic_hit]] <- TRUE
    }

    unmatched <- names(fitted_estimates)[remaining][!syntactic_hit]
    if (length(unmatched) > 0) {
      stop(
        sprintf(
          paste(
            "Could not align getSignifActivities() coefficients to preprocessed$N columns:",
            "%s"
          ),
          paste(head(unmatched, 10), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  unestimated_predictors <- names(coefficient_was_estimated)[!coefficient_was_estimated]
  if (length(unestimated_predictors) > 0) {
    message(sprintf(
      "getTargetMatrix(): %d TE predictors had no independent lm coefficient; their target columns will be NA.",
      length(unestimated_predictors)
    ))
  }
  
  delta_hat = c(rep(x = intercept + (preprocessed$N %*% A_hat), times = S)) # estimated Y according to the fitted model, repeated the number of contrasts
  
  err = sapply(delta-delta_hat, function(x) return(x*x)) # absolute error of the fitted model
  
  # err_mat is a Genes x Replicates matrix that contains the squared errors across all paired deltas
  err_mat = matrix(data = err, nrow = G, ncol = S, byrow=F)
  err_row = rowSums(err_mat)
  err_col = colSums(err_mat)
  
  dev_mat = matrix(0, nrow = G, ncol = S)
  for (i in 1:G) {
    for (j in 1:S) {
      dev_mat[i, j] = err_row[i] + err_col[j]
    }
  }
  dev_mat = dev_mat/(G*S) # dev_mat contains the averaged errors for each gene x sample pair
  
  N_prime = preprocessed$N # the matrix N that will contain the "mutations", i.e. rows set to zero
  
  target_genes = matrix(0, nrow = G, ncol = M) # res matrix for the impact of mutations on predictions
  rownames(target_genes) = rownames(preprocessed$E_centered)
  colnames(target_genes) = colnames(preprocessed$N)
  
  # start of loop
  pb <- txtProgressBar(min = 0, max = ncol(N_prime))
  i = 0
  for (s in colnames(N_prime)) { # s is the TE subfamily of interest
    if (s %in% unestimated_predictors) {
      target_genes[, s] = NA_real_
      i = i+1
      setTxtProgressBar(pb, i)
      next
    }
    
    N_prime[, s] = rep(0, times = nrow(N_prime)) # mutating N_prime: "removing" TEs
    
    delta_hat_prime = c(rep(x = intercept + (N_prime %*% A_hat), times = S)) # estimated Y according to the "mutated" model
    
    err_prime = sapply(delta-delta_hat_prime, function(x) return(x*x)) # absolute error according to the "mutated" model
    
    delta_err = (err_prime - err) # absolute difference in prediction between the fitted and mutated models. still too highly dimensional
    
    delta_err_rel = matrix(data = delta_err, nrow = nrow(target_genes), ncol = S)/dev_mat
    
    delta_err_rel_avg = apply(delta_err_rel, MARGIN = 1, mean) # relative difference between errors, averaged over replicates
    
    target_genes[, s] = delta_err_rel_avg
    
    N_prime[, s] = preprocessed$N[, s] # resetting N_prime
    i = i+1
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  return(target_genes)
}

#reverse mapping 
annotate_targets <- function(x, mapping) {
  data.frame(
    ensembl_gene_id = names(x),
    hgnc_symbol = mapping[names(x)],
    score = unname(x),
    stringsAsFactors = FALSE
  )
}

# what kind of genes are still unmapped+
classify_unmapped <- function(genes) {
  list(
    LOC      = genes[grepl("^LOC\\d+", genes)],
    SNORD    = genes[grepl("^SNORD|^SNORA", genes)],
    FLJ      = genes[grepl("^FLJ\\d+", genes)],
    LINC     = genes[grepl("^LINC\\d+", genes)],
    KIAA     = genes[grepl("^KIAA\\d+", genes)],
    MIR      = genes[grepl("^MIR\\d+", genes)],
    other    = genes[!grepl("^LOC\\d+|^SNORD|^SNORA|^FLJ|^LINC|^KIAA|^MIR", genes)]
  )
}

summarise_unmapped_classes <- function(classes) {
  if (!is.list(classes) || is.null(names(classes))) {
    classes <- classify_unmapped(classes)
  }

  counts <- lengths(classes)

  data.frame(
    class = names(counts),
    n = as.integer(counts),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# my gene mapping in smaller chunks
query_mygene_chunked <- function(ids,
                                 chunk_size = 100,
                                 scopes = c("symbol", "entrezgene", "alias"),
                                 fields = c("symbol", "name", "ensembl.gene", "type_of_gene"),
                                 species = "human") {
  
  chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))
  pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)
  
  res_list <- lapply(seq_along(chunks), function(i) {
    
    res <- NULL
    
    invisible(capture.output({
      res <- suppressMessages(
        queryMany(
          chunks[[i]],
          scopes = scopes,
          fields = fields,
          species = species,
          returnall = TRUE,
          as_dataframe = FALSE
        )
      )
    }))
    
    setTxtProgressBar(pb, i)
    
    as.data.frame(res$response)
  })
  
  close(pb)
  
  as.data.frame(data.table::rbindlist(res_list, fill = TRUE))
}

# ── 11. Full target matrix plot ───────────────────────────────────────────────
#
# Visualises which genes are most dependent on each TE subfamily as a heatmap.
# Rows = top genes (by mean absolute dependency across all TEs).
# Cols = top TE subfamilies (by cross-gene variance); highlight_TEs are always
#        included regardless of their variance rank.
# Color encodes row z-score: relative preference of a gene for a TE subfamily
#        compared to that gene's average across all shown TEs.
#
# Args:
#   target_mat_raw    – output of getTargetMatrix() (genes x TE subfamilies)
#   ensembl_to_symbol – named vector: Ensembl ID -> HGNC symbol
#   condition_label   – string appended to title and filename ("VPR" / "KRAB")
#   top_n_genes       – how many top genes to show (ranked by mean |dependency|)
#   top_n_TE          – how many top TE subfamilies to show (ranked by variance)
#   highlight_TEs     – TE subfamilies always included and colour-annotated
#   out_dir           – output directory for the PDF

plot_target_matrix_full <- function(
    target_mat_raw,
    ensembl_to_symbol,
    condition_label  = "VPR",
    top_n_genes      = 60,
    top_n_TE         = 30,
    highlight_TEs    = c("LTR5_Hs", "LTR5A", "LTR5B"),
    selected_genes   = NULL,
    selected_TEs     = NULL,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    treeheight_row   = 18,
    treeheight_col   = 50,
    legend           = TRUE,
    annotation_legend = TRUE,
    show_col_annotation = TRUE,
    show_row_annotation = TRUE,
    color_breaks     = NULL,
    legend_breaks    = NULL,
    legend_labels    = NULL,
    main             = NULL,
    gene_axis_label  = "Top Protein Coding Genes (ranked by mean |target dependency|)",
    te_axis_label    = "Top TE-subfamilies (ranked by cross-gene variance)",
    out_dir          = getOption("crates.figures_dir", "figures"),
    show             = FALSE
) {
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  mat <- as.matrix(target_mat_raw)
  if (is.null(rownames(mat))) {
    stop("target_mat_raw must have gene IDs as row names.", call. = FALSE)
  }
  
  # ── 1. Filter: drop rows/cols that are entirely NA; fill remaining NA -> 0 ──
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    stop("target_mat_raw has no non-empty gene/TE space to plot.", call. = FALSE)
  }
  mat[is.na(mat)] <- 0

  if (!is.null(selected_genes)) {
    selected_genes <- selected_genes[selected_genes %in% rownames(mat)]
    if (length(selected_genes) == 0L) {
      stop("selected_genes do not overlap target_mat_raw row names.", call. = FALSE)
    }
  }

  if (!is.null(selected_TEs)) {
    selected_TEs <- selected_TEs[selected_TEs %in% colnames(mat)]
    if (length(selected_TEs) == 0L) {
      stop("selected_TEs do not overlap target_mat_raw column names.", call. = FALSE)
    }
  }

  # Total dimensions after filtering (used in plot title)
  n_genes_total <- nrow(mat)
  n_TEs_total   <- ncol(mat)
  
  # ── 3. Select top N genes by mean absolute dependency ─────────────────────
  gene_score <- rowMeans(abs(mat), na.rm = TRUE)
  if (!is.null(selected_genes)) {
    top_genes <- selected_genes
  } else {
    top_genes  <- names(sort(gene_score, decreasing = TRUE))[
      seq_len(min(top_n_genes, nrow(mat)))
    ]
  }
  
  # ── 4. Select top N TEs by cross-gene variance ────────────────────────────
  #       highlight_TEs are forced in even if below the variance threshold
  te_var  <- apply(mat, 2, var, na.rm = TRUE)
  if (!is.null(selected_TEs)) {
    top_TEs <- selected_TEs
  } else {
    top_TEs <- names(sort(te_var, decreasing = TRUE))[
      seq_len(min(top_n_TE, ncol(mat)))
    ]
  }
  forced  <- highlight_TEs[highlight_TEs %in% colnames(mat) &
                             !highlight_TEs %in% top_TEs]
  top_TEs <- unique(c(forced, top_TEs))
  
  plot_mat <- mat[top_genes, top_TEs, drop = FALSE]

  # ── 5. Gene labels: replace Ensembl IDs with HGNC symbols where available ──
  if (is.null(ensembl_to_symbol)) {
    ensembl_to_symbol <- stats::setNames(rownames(plot_mat), rownames(plot_mat))
  }
  gene_labels <- ensembl_to_symbol[rownames(plot_mat)]
  gene_labels[is.na(gene_labels) | gene_labels == ""] <-
    rownames(plot_mat)[is.na(gene_labels) | gene_labels == ""]
  rownames(plot_mat) <- make.unique(gene_labels)   # suffix .1/.2 for true duplicates
  
  # ── 6. Column annotation: colour-code highlight_TEs, rest = "other" ───────
  highlight_present <- intersect(unique(highlight_TEs), colnames(plot_mat))
  fixed_hl_colors <- c(LTR5_Hs = "#E63946", LTR5A = "#F4A261", LTR5B = "#E9C46A")
  missing_hl_colors <- setdiff(highlight_present, names(fixed_hl_colors))
  dynamic_hl_colors <- if (length(missing_hl_colors) > 0L) {
    stats::setNames(grDevices::hcl.colors(length(missing_hl_colors), "Dark 3"), missing_hl_colors)
  } else {
    character(0)
  }
  hl_colors <- c(fixed_hl_colors[names(fixed_hl_colors) %in% highlight_present],
                 dynamic_hl_colors,
                 other = "#ADB5BD")
  hl_vec <- setNames(rep("other", ncol(plot_mat)), colnames(plot_mat))
  for (te in highlight_present) hl_vec[colnames(plot_mat) == te] <- te
  
  highlight_levels <- c(highlight_present, "other")
  col_ann    <- data.frame(
    Highlight = factor(hl_vec, levels = highlight_levels),
    row.names = colnames(plot_mat)
  )
  ann_colors <- list(
    Highlight  = hl_colors[highlight_levels],
    Rank_group = c("Top 10" = "#2D6A4F", "Top 25" = "#74C69D", "Rest" = "#D8F3DC")
  )
  
  # ── 7. Row annotation: rank-group band (Top 10 / Top 25 / Rest) ───────────
  rank_vec <- rank(-gene_score[top_genes], ties.method = "min")
  row_ann  <- data.frame(
    Rank_group = cut(rank_vec, breaks = c(0, 10, 25, Inf),
                     labels = c("Top 10", "Top 25", "Rest")),
    row.names  = rownames(plot_mat)
  )
  
  # ── 8. Sparse label logic: above THRESHOLD only the most relevant are named ─
  #       Row labels: only top 20 genes labelled when top_n_genes >= 100
  #       Col labels: only highlight_TEs labelled when top_n_TE   >= 100
  THRESHOLD <- 100
  
  if (top_n_genes >= THRESHOLD) {
    top20 <- names(sort(gene_score[top_genes], decreasing = TRUE))[
      seq_len(min(20L, length(top_genes)))
    ]
    labels_row <- ifelse(top_genes %in% top20, rownames(plot_mat), "")
  } else {
    labels_row <- rownames(plot_mat)
  }
  
  if (top_n_TE >= THRESHOLD) {
    labels_col <- ifelse(colnames(plot_mat) %in% highlight_TEs,
                         colnames(plot_mat), "")
  } else {
    labels_col <- colnames(plot_mat)
  }

  stack_target_matrix_legends <- function(gt, spacer_height = grid::unit(2, "mm")) {
    if (!requireNamespace("gtable", quietly = TRUE)) {
      return(gt)
    }

    legend_idx <- which(gt$layout$name == "legend")
    annotation_legend_idx <- which(gt$layout$name == "annotation_legend")
    if (length(legend_idx) == 0L || length(annotation_legend_idx) == 0L) {
      return(gt)
    }

    legend_idx <- legend_idx[[1]]
    annotation_legend_idx <- annotation_legend_idx[[1]]
    legend_grob <- gt$grobs[[legend_idx]]
    annotation_legend_grob <- gt$grobs[[annotation_legend_idx]]

    stacked_legend <- gtable::gtable(
      widths = grid::unit.pmax(
        grid::grobWidth(legend_grob),
        grid::grobWidth(annotation_legend_grob)
      ),
      heights = grid::unit.c(
        grid::grobHeight(legend_grob),
        spacer_height,
        grid::grobHeight(annotation_legend_grob)
      )
    )
    stacked_legend <- gtable::gtable_add_grob(
      stacked_legend,
      legend_grob,
      t = 1,
      l = 1,
      name = "legend"
    )
    stacked_legend <- gtable::gtable_add_grob(
      stacked_legend,
      annotation_legend_grob,
      t = 3,
      l = 1,
      name = "annotation_legend"
    )

    gt$grobs[[legend_idx]] <- stacked_legend
    gt$grobs[[annotation_legend_idx]] <- grid::nullGrob()
    gt$layout$clip[legend_idx] <- "off"

    annotation_cols <- seq(
      gt$layout$l[annotation_legend_idx],
      gt$layout$r[annotation_legend_idx]
    )
    legend_cols <- seq(gt$layout$l[legend_idx], gt$layout$r[legend_idx])
    collapse_cols <- setdiff(annotation_cols, legend_cols)
    if (length(collapse_cols) > 0L) {
      gt$widths[collapse_cols] <- grid::unit(0, "pt")
    }

    gt
  }

  draw_target_matrix_axis_labels <- function() {
    has_label <- function(label) {
      length(label) == 1L && !is.na(label) && nzchar(label)
    }

    draw_axis_label <- function(label, x, y, rot = 0) {
      label <- gsub("\\*\\*", "", label)
      split_match <- regexec("^([^()]*)\\s*(\\(.*\\))$", label)
      parts <- regmatches(label, split_match)[[1]]

      if (length(parts) == 3L) {
        bold_label <- trimws(parts[[2]])
        plain_label <- paste0(" ", parts[[3]])
        label_expr <- as.expression(bquote(bold(.(bold_label)) * .(plain_label)))
      } else {
        label_expr <- label
      }

      grid::grid.text(
        label_expr,
        x = x,
        y = y,
        rot = rot,
        gp = grid::gpar(fontsize = 8, col = "grey35")
      )
    }

    if (has_label(te_axis_label)) {
      draw_axis_label(
        te_axis_label,
        x = grid::unit(0.50, "npc"),
        y = grid::unit(0.018, "npc")
      )
    }

    if (has_label(gene_axis_label)) {
      draw_axis_label(
        gene_axis_label,
        x = grid::unit(0.982, "npc"),
        y = grid::unit(0.50, "npc"),
        rot = 270
      )
    }
  }
  
  # ── 9. Color palette: divergent blue-white-red (row z-score) ──────────────
  pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  
  plot_main <- if (!is.null(main)) {
    main
  } else {
    paste0(
      "Gene x TE-subfamily target dependency - ", condition_label,
      "\n(top ", length(top_genes), " of ", n_genes_total, " genes  |  ",
      "top ", length(top_TEs), " of ", n_TEs_total, " TE subfamilies  |  row z-score)"
    )
  }

  effective_cluster_rows <- isTRUE(cluster_rows) && nrow(plot_mat) > 1L
  effective_cluster_cols <- isTRUE(cluster_cols) && ncol(plot_mat) > 1L

  # ── 10. Build the pheatmap object silently and optionally save/show it ─────
  out_file   <- if (!is.null(out_dir)) {
    file.path(out_dir, paste0("target_matrix_full_", condition_label, ".pdf"))
  } else {
    NULL
  }
  height_val <- 0.18 * top_n_genes + 4   # figure height scales with gene count
  
  p <- pheatmap::pheatmap(
    plot_mat,
    scale             = "row",               # row z-score (see function header)
    color             = pal,
    breaks            = color_breaks,
    cluster_rows      = effective_cluster_rows,
    cluster_cols      = effective_cluster_cols,
    treeheight_row    = if (effective_cluster_rows) treeheight_row else 0,
    treeheight_col    = if (effective_cluster_cols) treeheight_col else 0,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    labels_row        = labels_row,
    labels_col        = labels_col,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    angle_col         = "45",
    fontsize_row      = ifelse(top_n_genes >= THRESHOLD, 6, 7),
    fontsize_col      = ifelse(top_n_TE    >= THRESHOLD, 6, 7),
    annotation_col    = if (isTRUE(show_col_annotation)) col_ann else NULL,
    annotation_row    = if (isTRUE(show_row_annotation)) row_ann else NULL,
    annotation_colors = if (isTRUE(show_col_annotation) || isTRUE(show_row_annotation)) ann_colors else NA,
    border_color      = NA,
    legend            = legend,
    legend_breaks     = legend_breaks,
    legend_labels     = legend_labels,
    annotation_legend = annotation_legend,
    main              = plot_main,
    silent = TRUE
  )
  p$gtable <- stack_target_matrix_legends(p$gtable)
  
  if (!is.null(out_file)) {
    Cairo::CairoPDF(out_file, width = 14, height = height_val)
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
    draw_target_matrix_axis_labels()
    dev.off()
    message("Saved: ", out_file)
  }
  
  if (show) {
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
    draw_target_matrix_axis_labels()
  }

  invisible(p)
}


classify_unmapped <- function(genes) {
  list(
    LOC      = genes[grepl("^LOC\\d+", genes)],
    SNORD    = genes[grepl("^SNORD|^SNORA", genes)],
    FLJ      = genes[grepl("^FLJ\\d+", genes)],
    LINC     = genes[grepl("^LINC\\d+", genes)],
    KIAA     = genes[grepl("^KIAA\\d+", genes)],
    MIR      = genes[grepl("^MIR\\d+", genes)],
    other    = genes[!grepl("^LOC\\d+|^SNORD|^SNORA|^FLJ|^LINC|^KIAA|^MIR", genes)]
  )
}
