# Comparison and plotting helpers for result_comparison.ipynb.
# The workflow classifies the selected result set globally before running
# comparison blocks. Pairwise calculations are only used internally where a
# similarity matrix requires them.

if (!exists("load_comparison_results", mode = "function")) {
  source("scripts/comparison_csv_f.R")
}
if (!exists("classify_te_senescence_trend", mode = "function")) {
  source("scripts/te_trend_classification_f.R")
}

comparison_need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for comparison plots.", call. = FALSE)
  }
  invisible(TRUE)
}

comparison_theme <- function() {
  comparison_need_ggplot()
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey35"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

comparison_zero_reference_line <- function() {
  comparison_need_ggplot()
  ggplot2::geom_vline(
    xintercept = 0,
    color = "#98a2b3",
    linewidth = 0.5,
    alpha = 0.85,
    show.legend = FALSE
  )
}

comparison_result_ids <- function(comparison_results) {
  ids <- names(comparison_results)
  if (is.null(ids) || any(!nzchar(ids))) {
    ids <- vapply(comparison_results, function(result) result$result_id, character(1))
  }
  ids
}

comparison_clean_matrix_id <- function(x, prefix) {
  x <- sub(sprintf("^%s_", prefix), "", x)
  x <- sub("_craTEs$", "", x)
  x
}

comparison_pretty_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

comparison_window_label <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) {
    return(x)
  }
  if (value >= 1000) {
    kb <- value / 1000
    return(sprintf("%gkb", kb))
  }
  sprintf("%gkb", value)
}

comparison_n_components <- function(n_id) {
  if (is.na(n_id) || !nzchar(n_id)) {
    return(data.frame(
      n_window = NA_character_,
      n_mode_code = NA_character_,
      n_genome = NA_character_,
      n_origin = NA_character_,
      n_source = NA_character_,
      n_base_label = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  clean <- sub("^N_", "", n_id)
  parts <- strsplit(clean, "_", fixed = TRUE)[[1]]
  mode_code <- parts[1]
  mode_label <- if (identical(mode_code, "c")) {
    "count"
  } else if (identical(mode_code, "w")) {
    "weighted"
  } else {
    comparison_pretty_label(mode_code)
  }
  window_label <- if (length(parts) >= 2L) comparison_window_label(parts[2]) else NA_character_
  genome_idx <- grep("^Hg[0-9]+$", parts)
  genome <- parts[genome_idx]
  genome <- if (length(genome) == 0L) NA_character_ else genome[1]
  origin_idx <- setdiff(seq_along(parts), c(1L, 2L, genome_idx))
  origin <- if (length(origin_idx) > 0L) paste(parts[origin_idx], collapse = "_") else NA_character_
  source <- if (!is.na(origin) && nzchar(origin)) {
    comparison_pretty_label(origin)
  } else {
    NA_character_
  }

  data.frame(
    n_window = window_label,
    n_mode_code = mode_label,
    n_genome = genome,
    n_origin = origin,
    n_source = source,
    n_base_label = paste(na.omit(c(window_label, mode_label)), collapse = " "),
    stringsAsFactors = FALSE
  )
}

comparison_e_source_label <- function(e_id) {
  clean <- comparison_clean_matrix_id(e_id, "E")
  if (is.na(clean) || !nzchar(clean)) {
    return(NA_character_)
  }
  comparison_pretty_label(clean)
}

comparison_shortest_unique_token_label <- function(ids, prefix) {
  clean_ids <- vapply(ids, comparison_clean_matrix_id, character(1), prefix = prefix)
  clean_ids[is.na(clean_ids)] <- NA_character_
  labels <- clean_ids
  unique_values <- unique(clean_ids[!is.na(clean_ids)])

  for (value in unique_values) {
    idx <- which(clean_ids == value)
    tokens <- strsplit(value, "_", fixed = TRUE)[[1]]
    candidates <- unique(c(tokens, vapply(seq_along(tokens), function(i) paste(tokens[seq_len(i)], collapse = "_"), character(1)), value))
    candidates <- candidates[nzchar(candidates)]
    other_values <- setdiff(unique_values, value)
    is_unique <- vapply(candidates, function(candidate) !any(grepl(candidate, other_values, fixed = TRUE)), logical(1))
    labels[idx] <- if (any(is_unique)) candidates[which(is_unique)[1]] else value
  }
  comparison_pretty_label(labels)
}

comparison_make_unique_labels <- function(base_labels, fallback_labels) {
  out <- base_labels
  duplicated_labels <- unique(out[duplicated(out) | duplicated(out, fromLast = TRUE)])
  for (label in duplicated_labels) {
    idx <- which(out == label)
    out[idx] <- fallback_labels[idx]
  }
  out
}

comparison_add_matrix_plot_ids <- function(meta) {
  n_components <- do.call(rbind, lapply(meta$n_id, comparison_n_components))
  meta <- cbind(meta, n_components)

  n_label_genome <- ifelse(
    is.na(meta$n_genome) | !nzchar(meta$n_genome),
    meta$n_base_label,
    paste(meta$n_base_label, meta$n_genome, sep = " ")
  )
  meta$n_plot_id <- comparison_make_unique_labels(meta$n_base_label, n_label_genome)
  meta$n_plot_id <- comparison_make_unique_labels(meta$n_plot_id, comparison_pretty_label(comparison_clean_matrix_id(meta$n_id, "N")))
  meta$n_plot_id <- comparison_pretty_label(meta$n_plot_id)

  meta$e_plot_id <- comparison_shortest_unique_token_label(meta$e_id, "E")
  meta$e_source <- vapply(meta$e_id, comparison_e_source_label, character(1))
  meta$matrix_plot_id <- paste(meta$e_plot_id, meta$n_plot_id, sep = " | ")
  meta$matrix_plot_id <- comparison_pretty_label(meta$matrix_plot_id)
  meta
}

parse_result_metadata <- function(result_dir) {
  result_id <- basename(normalizePath(result_dir, winslash = "/", mustWork = FALSE))
  core_id <- sub("^results_", "", result_id)

  e_id <- NA_character_
  n_id <- NA_character_
  if (grepl("_N_", core_id, fixed = TRUE)) {
    parts <- strsplit(core_id, "_N_", fixed = TRUE)[[1]]
    e_id <- parts[1]
    n_id <- paste0("N_", paste(parts[-1], collapse = "_N_"))
  }

  genome_match <- regmatches(core_id, regexpr("Hg[0-9]+", core_id, ignore.case = FALSE))
  genome <- if (length(genome_match) > 0L && nzchar(genome_match)) genome_match else NA_character_
  n_mode <- if (grepl("^N_c_", n_id)) {
    "count"
  } else if (grepl("^N_w_", n_id)) {
    "weighted"
  } else {
    NA_character_
  }

  data.frame(
    result_id = result_id,
    result_dir = result_dir,
    e_id = e_id,
    n_id = n_id,
    genome = genome,
    n_mode = n_mode,
    stringsAsFactors = FALSE
  )
}

comparison_metadata_table <- function(comparison_results) {
  meta <- do.call(
    rbind,
    lapply(comparison_results, function(result) parse_result_metadata(result$result_dir))
  )
  rownames(meta) <- NULL
  comparison_add_matrix_plot_ids(meta)
}

store_matrix_plot_ids <- function(comparison_results) {
  meta <- comparison_metadata_table(comparison_results)
  for (i in seq_along(comparison_results)) {
    comparison_results[[i]]$metadata <- meta[i, , drop = FALSE]
  }
  attr(comparison_results, "matrix_metadata") <- meta
  comparison_results
}

get_matrix_plot_metadata <- function(comparison_results) {
  meta <- attr(comparison_results, "matrix_metadata", exact = TRUE)
  if (is.null(meta)) {
    meta <- comparison_metadata_table(comparison_results)
  }
  meta
}

classify_comparison_set <- function(comparison_results) {
  meta <- get_matrix_plot_metadata(comparison_results)
  e_ids <- unique(meta$e_id)
  n_ids <- unique(meta$n_id)
  same_E <- length(e_ids[!is.na(e_ids)]) == 1L && !any(is.na(meta$e_id))
  same_N <- length(n_ids[!is.na(n_ids)]) == 1L && !any(is.na(meta$n_id))

  class_label <- if (same_E && same_N) {
    "same E, same N"
  } else if (same_E && !same_N) {
    "same E, different N"
  } else if (!same_E && same_N) {
    "different E, same N"
  } else {
    "different E, different N"
  }

  structure(
    list(
      class_label = class_label,
      same_E = same_E,
      same_N = same_N,
      e_ids = e_ids,
      n_ids = n_ids,
      metadata = meta,
      n_results = nrow(meta)
    ),
    class = "comparison_class"
  )
}

comparison_notebook_context <- function() {
  isTRUE(getOption("jupyter.in_kernel")) || nzchar(Sys.getenv("JPY_PARENT_PID"))
}

comparison_class_table <- function(comparison_class) {
  meta <- comparison_class$metadata
  meta$result <- comparison_short_result_label(meta$result_id)
  meta[, c("result", "result_id", "result_dir", "e_source", "n_genome", "n_mode", "n_source"), drop = FALSE]
}

comparison_class_html <- function(comparison_class) {
  meta <- comparison_class_table(comparison_class)
  badge_class <- gsub("[^a-z0-9]+", "-", tolower(comparison_class$class_label))

  flag_cell <- function(value) {
    class_name <- if (isTRUE(value)) "present" else "missing"
    sprintf(
      "<span class=\"%s\">%s %s</span>",
      class_name,
      comparison_presence_mark(value),
      if (isTRUE(value)) "TRUE" else "FALSE"
    )
  }

  rows <- vapply(seq_len(nrow(meta)), function(i) {
    x <- meta[i, , drop = FALSE]
    paste0(
      "<tr>",
      "<td class=\"result-name\">", comparison_result_link(x$result_dir, x$result_id), "</td>",
      "<td>", comparison_escape_html(x$e_source), "</td>",
      "<td>", comparison_escape_html(x$n_genome), "</td>",
      "<td>", comparison_escape_html(x$n_mode), "</td>",
      "<td>", comparison_escape_html(x$n_source), "</td>",
      "</tr>"
    )
  }, character(1))

  paste0(
    "<style>",
    ".comparison-class-wrap{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:100%;}",
    ".comparison-class-summary{display:flex;gap:6px;align-items:center;flex-wrap:wrap;margin:0 0 6px 0;}",
    ".comparison-class-badge{background:#eef4ff;color:#1849a9;border:1px solid #b2ccff;border-radius:4px;padding:3px 6px;font-weight:700;}",
    ".comparison-class-badge.", badge_class, "{}",
    ".comparison-class-flag{background:#f8fafc;border:1px solid #e5e7eb;border-radius:4px;padding:3px 6px;color:#344054;}",
    ".comparison-class-flag .present{color:#177245;font-weight:700;}",
    ".comparison-class-flag .missing{color:#b42318;font-weight:700;}",
    ".comparison-class{border-collapse:collapse;width:auto;max-width:100%;font-size:12px;line-height:1.14;}",
    ".comparison-class th{background:#f4f6f8;color:#1f2933;text-align:left;border-bottom:1px solid #d8dee4;padding:4px 5px;white-space:nowrap;}",
    ".comparison-class td{border-bottom:1px solid #e5e7eb;padding:3px 5px;vertical-align:top;text-align:left;}",
    ".comparison-class th:first-child{text-align:left;}",
    ".comparison-class .result-name{font-weight:600;min-width:160px;max-width:360px;text-align:left;}",
    ".comparison-class .result-name a{color:#0969da;text-decoration:none;}",
    ".comparison-class .id-cell{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:#344054;max-width:260px;overflow-wrap:anywhere;}",
    "</style>",
    "<div class=\"comparison-class-wrap\">",
    "<div class=\"comparison-class-summary\">",
    "<div class=\"comparison-class-badge ", badge_class, "\">Dataset similarities: ",
    comparison_escape_html(comparison_class$class_label), "</div>",
    "<div class=\"comparison-class-flag\">same E: ", flag_cell(comparison_class$same_E), "</div>",
    "<div class=\"comparison-class-flag\">same N: ", flag_cell(comparison_class$same_N), "</div>",
    "</div>",
    "<table class=\"comparison-class\">",
    "<thead><tr><th>result</th><th>E source</th><th>ref. genome</th><th>N mode</th><th>N source</th></tr></thead>",
    "<tbody>", paste(rows, collapse = ""), "</tbody>",
    "</table>",
    "</div>"
  )
}

print_comparison_class <- function(comparison_class, html = TRUE) {
  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_class_html(comparison_class))
    return(invisible(comparison_class))
  }

  cat("Dataset similarities:", comparison_class$class_label, "\n")
  cat("All selected result directories share E:", comparison_class$same_E, "\n")
  cat("All selected result directories share N:", comparison_class$same_N, "\n\n")
  meta <- comparison_class_table(comparison_class)
  meta <- meta[, c("result", "e_source", "n_genome", "n_mode", "n_source"), drop = FALSE]
  print(meta, row.names = FALSE)
  invisible(comparison_class)
}

comparison_short_result_label <- function(result_id) {
  sub("^results_", "", result_id)
}

comparison_escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

comparison_result_link <- function(result_dir, result_id = basename(result_dir)) {
  label <- comparison_short_result_label(result_id)
  normalized <- normalizePath(result_dir, winslash = "/", mustWork = FALSE)
  href <- paste0(
    "https://file+.vscode-resource.vscode-cdn.net",
    utils::URLencode(normalized, reserved = FALSE)
  )
  sprintf(
    "<a href=\"%s\" title=\"%s\">%s</a>",
    comparison_escape_html(href),
    comparison_escape_html(normalized),
    comparison_escape_html(label)
  )
}

comparison_presence_mark <- function(present) {
  ifelse(isTRUE(present), "&#10003;", "&#10007;")
}

comparison_presence_text <- function(present) {
  ifelse(isTRUE(present), "[x]", "[ ]")
}

result_manifest_table <- function(comparison_results) {
  rows <- lapply(seq_along(comparison_results), function(i) {
    result <- comparison_results[[i]]
    fit <- result$tables$activity_fit_summary
    sample_diag <- result$tables$sample_activity_diagnostics
    lambda_diag <- result$tables$activity_lambda_diagnostics
    activity <- result$tables$activity_scores
    contrast <- result$tables$te_contrast_summary
    contrast_diag <- result$tables$signif_model_diagnostics
    target_summary <- result$tables$target_summary
    n_target_matrices <- length(result$target_gene_matrices)

    data.frame(
      result_label = comparison_short_result_label(result$result_id),
      result_id = result$result_id,
      result_dir = result$result_dir,
      has_activity_fit_summary = !is.null(fit),
      has_sample_activity_diagnostics = !is.null(sample_diag),
      has_activity_lambda_diagnostics = !is.null(lambda_diag),
      has_activity_scores = !is.null(activity),
      has_te_contrast_summary = !is.null(contrast),
      has_signif_model_diagnostics = !is.null(contrast_diag),
      has_target_gene_matrix = n_target_matrices > 0L,
      has_target_summary = !is.null(target_summary),
      n_activity_te_subfamilies = if (is.null(activity)) NA_integer_ else nrow(activity),
      n_activity_samples = if (is.null(activity)) NA_integer_ else max(0L, ncol(activity) - 1L),
      n_activity_fit_rows = if (is.null(fit)) NA_integer_ else nrow(fit),
      n_sample_diagnostic_rows = if (is.null(sample_diag)) NA_integer_ else nrow(sample_diag),
      n_lambda_diagnostic_rows = if (is.null(lambda_diag)) NA_integer_ else nrow(lambda_diag),
      n_te_contrast_subfamilies = if (is.null(contrast)) NA_integer_ else nrow(contrast),
      n_signif_model_diagnostic_rows = if (is.null(contrast_diag)) NA_integer_ else nrow(contrast_diag),
      n_target_gene_matrices = n_target_matrices,
      missing_optional_tsv_files = paste(result$missing_optional, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

comparison_manifest_html <- function(manifest) {
  status_cell <- function(present, detail = NULL) {
    class_name <- if (isTRUE(present)) "present" else "missing"
    mark <- comparison_presence_mark(present)
    detail <- if (is.null(detail) || is.na(detail) || !nzchar(as.character(detail))) "" else sprintf("<span>%s</span>", comparison_escape_html(detail))
    sprintf("<td class=\"%s\">%s%s</td>", class_name, mark, detail)
  }

  detail_count <- function(value) {
    value <- suppressWarnings(as.integer(value))
    if (length(value) != 1L || is.na(value) || value <= 0L) {
      return(NULL)
    }
    paste0("n=", value)
  }

  count_cell <- function(value) {
    value <- ifelse(is.na(value), "&mdash;", as.character(value))
    sprintf("<td class=\"count\">%s</td>", value)
  }

  rows <- vapply(seq_len(nrow(manifest)), function(i) {
    x <- manifest[i, , drop = FALSE]
    paste0(
      "<tr>",
      "<td class=\"result-name\">", comparison_result_link(x$result_dir, x$result_id), "</td>",
      status_cell(x$has_activity_fit_summary),
      status_cell(x$has_sample_activity_diagnostics, detail_count(x$n_sample_diagnostic_rows)),
      status_cell(x$has_activity_lambda_diagnostics, detail_count(x$n_lambda_diagnostic_rows)),
      status_cell(x$has_activity_scores),
      status_cell(x$has_te_contrast_summary),
      status_cell(x$has_signif_model_diagnostics, detail_count(x$n_signif_model_diagnostic_rows)),
      status_cell(x$has_target_gene_matrix, detail_count(x$n_target_gene_matrices)),
      status_cell(x$has_target_summary),
      count_cell(x$n_activity_te_subfamilies),
      count_cell(x$n_activity_samples),
      count_cell(x$n_activity_fit_rows),
      count_cell(x$n_te_contrast_subfamilies),
      count_cell(x$n_target_gene_matrices),
      "</tr>"
    )
  }, character(1))

  paste0(
    "<style>",
    ".comparison-manifest-wrap{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:100%;}",
    ".comparison-manifest{border-collapse:collapse;width:auto;max-width:100%;font-size:12px;line-height:1.14;}",
    ".comparison-manifest th{background:#f4f6f8;color:#1f2933;text-align:left;border-bottom:1px solid #d8dee4;padding:4px 5px;white-space:nowrap;}",
    ".comparison-manifest td{border-bottom:1px solid #e5e7eb;padding:3px 5px;vertical-align:top;text-align:left;}",
    ".comparison-manifest th:first-child{text-align:left;}",
    ".comparison-manifest .result-name{font-weight:600;min-width:160px;max-width:360px;text-align:left;}",
    ".comparison-manifest .result-name a{color:#0969da;text-decoration:none;}",
    ".comparison-manifest .id-cell{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:#344054;max-width:210px;overflow-wrap:anywhere;}",
    ".comparison-manifest .present{color:#177245;font-weight:700;text-align:center;}",
    ".comparison-manifest .missing{color:#b42318;font-weight:700;text-align:center;}",
    ".comparison-manifest .present span,.comparison-manifest .missing span{display:block;color:#667085;font-weight:500;font-size:9px;margin-top:1px;}",
    ".comparison-manifest .count{text-align:right;font-variant-numeric:tabular-nums;color:#344054;}",
    ".comparison-manifest-note{margin-top:5px;color:#667085;font-size:10px;}",
    "</style>",
    "<div class=\"comparison-manifest-wrap\">",
    "<table class=\"comparison-manifest\">",
    "<thead><tr>",
    "<th>result</th>",
    "<th>model fit<br><code>activity_fit_summary</code></th>",
    "<th>sample fit<br><code>sample_activity_diagnostics</code></th>",
    "<th>lambda grid<br><code>activity_lambda_diagnostics</code></th>",
    "<th>activity scores<br><code>activity_scores</code></th>",
    "<th>TE contrasts<br><code>te_contrast_summary</code></th>",
    "<th>contrast fit<br><code>signif_model_diagnostics</code></th>",
    "<th>target matrix<br><code>target_gene_matrix_*</code></th>",
    "<th>target summary<br><code>target_summary</code></th>",
    "<th>activity TEs</th><th>samples</th><th>fit rows</th><th>contrast TEs</th><th>target matrices</th>",
    "</tr></thead><tbody>",
    paste(rows, collapse = ""),
    "</tbody></table>",
    "<div class=\"comparison-manifest-note\">&#10003; present, &#10007; missing. ",
    "<code>*_signif_act.tsv</code> is loaded only to recover raw p-values for combined evidence labels.</div>",
    "</div>"
  )
}

print_result_manifest <- function(comparison_results, html = TRUE) {
  manifest <- result_manifest_table(comparison_results)

  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_manifest_html(manifest))
  } else {
    console_manifest <- data.frame(
      result = manifest$result_label,
      model_fit = vapply(manifest$has_activity_fit_summary, comparison_presence_text, character(1)),
      sample_fit = vapply(manifest$has_sample_activity_diagnostics, comparison_presence_text, character(1)),
      lambda_grid = vapply(manifest$has_activity_lambda_diagnostics, comparison_presence_text, character(1)),
      activity_scores = vapply(manifest$has_activity_scores, comparison_presence_text, character(1)),
      TE_contrasts = vapply(manifest$has_te_contrast_summary, comparison_presence_text, character(1)),
      contrast_fit = vapply(manifest$has_signif_model_diagnostics, comparison_presence_text, character(1)),
      target_matrix = ifelse(
        manifest$has_target_gene_matrix,
        paste0("[x] n=", manifest$n_target_gene_matrices),
        "[ ]"
      ),
      target_summary = vapply(manifest$has_target_summary, comparison_presence_text, character(1)),
      activity_TEs = manifest$n_activity_te_subfamilies,
      samples = manifest$n_activity_samples,
      fit_rows = manifest$n_activity_fit_rows,
      contrast_TEs = manifest$n_te_contrast_subfamilies,
      stringsAsFactors = FALSE
    )
    print(console_manifest, row.names = FALSE)
    cat("*_signif_act.tsv is loaded only for raw p-values\n")
  }
  invisible(manifest)
}

comparison_activity_te_ids <- function(result) {
  unique(result$tables$activity_scores$te_subfamily)
}

comparison_activity_sample_ids <- function(result) {
  setdiff(colnames(result$tables$activity_scores), "te_subfamily")
}

get_shared_te_subfamilies <- function(comparison_results) {
  te_sets <- lapply(comparison_results, comparison_activity_te_ids)
  Reduce(intersect, te_sets)
}

get_shared_samples <- function(comparison_results) {
  sample_sets <- lapply(comparison_results, comparison_activity_sample_ids)
  Reduce(intersect, sample_sets)
}

get_shared_conditions <- function(comparison_results) {
  condition_sets <- lapply(comparison_results, function(result) {
    fit <- result$tables$activity_fit_summary
    if (!"condition" %in% colnames(fit)) {
      return(character())
    }
    setdiff(unique(as.character(fit$condition)), c("all", NA_character_))
  })
  Reduce(intersect, condition_sets)
}

get_common_target_matrix_names <- function(comparison_results) {
  matrix_sets <- lapply(comparison_results, function(result) names(result$target_gene_matrices))
  if (any(vapply(matrix_sets, length, integer(1)) == 0L)) {
    return(character())
  }
  Reduce(intersect, matrix_sets)
}

get_shared_target_genes <- function(comparison_results) {
  matrix_names <- get_common_target_matrix_names(comparison_results)
  if (length(matrix_names) == 0L) {
    return(character())
  }

  gene_sets <- lapply(comparison_results, function(result) {
    unique(unlist(lapply(result$target_gene_matrices[matrix_names], function(x) x$row_id), use.names = FALSE))
  })
  Reduce(intersect, gene_sets)
}

comparison_dataset_labels <- function(comparison_results, comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  meta <- comparison_class$metadata
  fallback <- comparison_short_result_label(meta$result_id)
  labels <- if (isTRUE(comparison_class$same_E) && !isTRUE(comparison_class$same_N)) {
    meta$n_plot_id
  } else if (!isTRUE(comparison_class$same_E) && isTRUE(comparison_class$same_N)) {
    meta$e_plot_id
  } else if (!isTRUE(comparison_class$same_E) && !isTRUE(comparison_class$same_N)) {
    meta$matrix_plot_id
  } else {
    fallback
  }

  labels[is.na(labels) | !nzchar(labels)] <- fallback[is.na(labels) | !nzchar(labels)]
  labels <- comparison_make_unique_labels(labels, fallback)
  names(labels) <- meta$result_id
  labels
}

comparison_named_feature_sets <- function(comparison_results, labels, extractor) {
  ids <- comparison_result_ids(comparison_results)
  sets <- lapply(comparison_results, extractor)
  names(sets) <- unname(labels[ids])
  sets
}

comparison_clean_feature_values <- function(x) {
  x <- unique(as.character(x))
  x[!is.na(x) & nzchar(x)]
}

comparison_target_summary_gene_ids <- function(result) {
  target_summary <- result$tables$target_summary
  if (is.null(target_summary)) {
    return(character())
  }
  gene_col <- intersect(c("row_id", "gene_id", "ensembl_gene_id"), colnames(target_summary))
  if (length(gene_col) == 0L) {
    return(character())
  }
  comparison_clean_feature_values(target_summary[[gene_col[1]]])
}

comparison_target_matrix_gene_ids <- function(result, matrix_names = NULL) {
  if (length(result$target_gene_matrices) == 0L) {
    return(character())
  }

  use_names <- names(result$target_gene_matrices)
  if (!is.null(matrix_names)) {
    use_names <- intersect(use_names, matrix_names)
  }
  if (length(use_names) == 0L) {
    return(character())
  }

  comparison_clean_feature_values(unlist(
    lapply(result$target_gene_matrices[use_names], function(x) x$row_id),
    use.names = FALSE
  ))
}

comparison_gene_id_sets <- function(comparison_results, labels) {
  ids <- comparison_result_ids(comparison_results)
  common_matrix_names <- get_common_target_matrix_names(comparison_results)

  if (length(common_matrix_names) > 0L) {
    sets <- lapply(
      comparison_results,
      comparison_target_matrix_gene_ids,
      matrix_names = common_matrix_names
    )
    source <- sprintf(
      "target_gene_matrix_*.tsv (%s)",
      paste(sub("^target_gene_matrix_", "", common_matrix_names), collapse = ", ")
    )
  } else {
    summary_sets <- lapply(comparison_results, comparison_target_summary_gene_ids)
    if (all(vapply(summary_sets, length, integer(1)) > 0L)) {
      sets <- summary_sets
      source <- "target_summary.tsv fallback"
    } else {
      sets <- lapply(comparison_results, comparison_target_matrix_gene_ids)
      missing_matrix <- vapply(sets, length, integer(1)) == 0L
      if (any(missing_matrix)) {
        sets[missing_matrix] <- summary_sets[missing_matrix]
      }
      source <- "target_gene_matrix_*.tsv / target_summary.tsv fallback"
    }
  }

  names(sets) <- unname(labels[ids])
  list(sets = sets, source = source)
}

comparison_shared_feature_set_status <- function(sets) {
  if (is.null(sets)) {
    return("skipped")
  }
  totals <- vapply(sets, length, integer(1))
  if (length(totals) < 2L) {
    return("skipped")
  }
  if (all(totals == 0L)) {
    return("unavailable")
  }
  if (any(totals == 0L)) {
    return("partial")
  }
  "ok"
}

comparison_feature_class_note <- function(category_key, comparison_class) {
  switch(
    category_key,
    samples = if (isTRUE(comparison_class$same_E)) {
      "Same E: sample overlap is mainly a sanity check."
    } else {
      "Different E: sample overlap is descriptive only."
    },
    gene_ids = if (isTRUE(comparison_class$same_E) && !isTRUE(comparison_class$same_N)) {
      "Same E, different N: useful for target-support stability."
    } else if (isTRUE(comparison_class$same_N)) {
      "Same N: target support can be compared on shared genes."
    } else {
      "Different E/N: descriptive target-gene compatibility."
    },
    te_subfamilies = if (isTRUE(comparison_class$same_N)) {
      "Same N: TE-space overlap should be high."
    } else {
      "Different N: shared TEs define the comparable activity rows."
    },
    significant_genes = "Skipped: no dedicated significant-gene source is defined.",
    "Descriptive overlap."
  )
}

comparison_common_matrix_attributes <- function(comparison_class, matrix = c("all", "E", "N")) {
  matrix <- match.arg(matrix)
  meta <- comparison_class$metadata
  attrs <- switch(
    matrix,
    E = list(source = meta$e_source),
    N = list(`ref. genome` = meta$n_genome, mode = meta$n_mode, source = meta$n_source),
    all = list(
      `E source` = meta$e_source,
      `N ref. genome` = meta$n_genome,
      `N mode` = meta$n_mode,
      `N source` = meta$n_source
    )
  )

  values <- vapply(attrs, function(x) {
    x <- unique(comparison_pretty_label(x[!is.na(x) & nzchar(x)]))
    if (length(x) == 1L) {
      x
    } else {
      NA_character_
    }
  }, character(1))

  values[!is.na(values) & nzchar(values)]
}

comparison_common_matrix_attribute_line <- function(comparison_class, matrix = c("E", "N")) {
  matrix <- match.arg(matrix)
  values <- comparison_common_matrix_attributes(comparison_class, matrix = matrix)
  if (length(values) == 0L) {
    return(sprintf("Common %s: none detected.", matrix))
  }
  paste0(
    "Common ", matrix, ": ",
    paste(sprintf("%s = %s", names(values), unname(values)), collapse = "; "),
    "."
  )
}

comparison_common_matrix_attribute_text <- function(comparison_class) {
  paste(
    comparison_common_matrix_attribute_line(comparison_class, "E"),
    comparison_common_matrix_attribute_line(comparison_class, "N"),
    sep = "\n"
  )
}

comparison_plot_subtitle <- function(main, comparison_class, extra = NULL) {
  paste(
    c(main, comparison_common_matrix_attribute_text(comparison_class), extra),
    collapse = "\n"
  )
}

comparison_feature_attribute <- function(category_key) {
  switch(
    category_key,
    samples = "activity-score columns",
    gene_ids = "target-gene row IDs",
    te_subfamilies = "activity-score TE rows",
    significant_genes = "not analysed",
    category_key
  )
}

comparison_feature_display_name <- function(category_key) {
  switch(
    category_key,
    samples = "Samples",
    gene_ids = "Gene IDs",
    te_subfamilies = "TE subfamilies",
    significant_genes = "Significant genes",
    category_key
  )
}

comparison_feature_key <- function(category) {
  if (category %in% c("samples", "gene_ids", "te_subfamilies", "significant_genes")) {
    return(category)
  }
  key <- tolower(gsub("[^a-z0-9]+", "_", category))
  key <- sub("_$", "", key)
  aliases <- c(
    sample = "samples",
    samples = "samples",
    gene = "gene_ids",
    genes = "gene_ids",
    gene_id = "gene_ids",
    gene_ids = "gene_ids",
    te = "te_subfamilies",
    tes = "te_subfamilies",
    te_subfamily = "te_subfamilies",
    te_subfamilies = "te_subfamilies",
    significant_gene = "significant_genes",
    significant_genes = "significant_genes"
  )
  if (!key %in% names(aliases)) {
    stop(sprintf("Unknown shared feature category: %s", category), call. = FALSE)
  }
  aliases[[key]]
}

comparison_feature_union <- function(sets) {
  comparison_clean_feature_values(unlist(sets, use.names = FALSE))
}

comparison_feature_intersection <- function(sets) {
  if (length(sets) == 0L) {
    return(character())
  }
  Reduce(intersect, sets)
}

comparison_feature_pct <- function(n, denom) {
  if (is.na(denom) || denom == 0L) {
    return(NA_real_)
  }
  100 * n / denom
}

comparison_feature_absolute_table <- function(category_key, sets) {
  if (is.null(sets)) {
    return(data.frame())
  }

  union_values <- comparison_feature_union(sets)
  shared_values <- comparison_feature_intersection(sets)
  rows <- lapply(names(sets), function(dataset) {
    other_names <- setdiff(names(sets), dataset)
    other_values <- if (length(other_names) == 0L) character() else comparison_feature_union(sets[other_names])
    only_values <- setdiff(sets[[dataset]], other_values)
    data.frame(
      category_key = category_key,
      category = comparison_feature_display_name(category_key),
      dataset = dataset,
      total_n = length(sets[[dataset]]),
      union_n = length(union_values),
      shared_all_n = length(shared_values),
      only_dataset_n = length(only_values),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

comparison_feature_relative_table <- function(category_key, sets) {
  absolute <- comparison_feature_absolute_table(category_key, sets)
  if (nrow(absolute) == 0L) {
    return(absolute)
  }

  absolute$shared_pct_union <- mapply(comparison_feature_pct, absolute$shared_all_n, absolute$union_n)
  absolute$shared_pct_dataset <- mapply(comparison_feature_pct, absolute$shared_all_n, absolute$total_n)
  absolute$only_pct_union <- mapply(comparison_feature_pct, absolute$only_dataset_n, absolute$union_n)
  absolute$only_pct_dataset <- mapply(comparison_feature_pct, absolute$only_dataset_n, absolute$total_n)
  absolute
}

collect_shared_feature_sets <- function(comparison_results, comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  labels <- comparison_dataset_labels(comparison_results, comparison_class)
  gene_sets <- comparison_gene_id_sets(comparison_results, labels)

  sets <- list(
    samples = comparison_named_feature_sets(comparison_results, labels, comparison_activity_sample_ids),
    gene_ids = gene_sets$sets,
    te_subfamilies = comparison_named_feature_sets(comparison_results, labels, comparison_activity_te_ids),
    significant_genes = NULL
  )

  sources <- c(
    samples = "activity_scores.tsv",
    gene_ids = gene_sets$source,
    te_subfamilies = "activity_scores.tsv",
    significant_genes = "skipped"
  )

  list(
    comparison_class = comparison_class,
    labels = labels,
    sets = sets,
    sources = sources
  )
}

build_shared_feature_space <- function(comparison_results, comparison_class = NULL) {
  feature_sets <- collect_shared_feature_sets(comparison_results, comparison_class)
  comparison_class <- feature_sets$comparison_class
  category_keys <- c("samples", "gene_ids", "te_subfamilies", "significant_genes")

  overview <- do.call(rbind, lapply(category_keys, function(category_key) {
    sets <- feature_sets$sets[[category_key]]
    status <- comparison_shared_feature_set_status(sets)
    union_n <- if (is.null(sets)) NA_integer_ else length(comparison_feature_union(sets))
    shared_all_n <- if (is.null(sets)) NA_integer_ else length(comparison_feature_intersection(sets))
    totals <- if (is.null(sets)) {
      "not analysed"
    } else {
      paste(sprintf("%s n=%s", names(sets), vapply(sets, length, integer(1))), collapse = "; ")
    }

    data.frame(
      category_key = category_key,
      category = comparison_feature_display_name(category_key),
      attribute = comparison_feature_attribute(category_key),
      source = unname(feature_sets$sources[[category_key]]),
      status = status,
      union_n = union_n,
      shared_all_n = shared_all_n,
      dataset_totals = totals,
      interpretation = comparison_feature_class_note(category_key, comparison_class),
      stringsAsFactors = FALSE
    )
  }))

  active_keys <- category_keys[!vapply(feature_sets$sets[category_keys], is.null, logical(1))]
  absolute_counts <- do.call(rbind, lapply(active_keys, function(key) {
    comparison_feature_absolute_table(key, feature_sets$sets[[key]])
  }))
  relative_counts <- do.call(rbind, lapply(active_keys, function(key) {
    comparison_feature_relative_table(key, feature_sets$sets[[key]])
  }))

  if (is.null(absolute_counts)) {
    absolute_counts <- data.frame()
  }
  if (is.null(relative_counts)) {
    relative_counts <- data.frame()
  }

  list(
    comparison_class = comparison_class,
    labels = feature_sets$labels,
    sets = feature_sets$sets,
    sources = feature_sets$sources,
    overview = overview,
    absolute_counts = absolute_counts,
    relative_counts = relative_counts
  )
}

comparison_format_percent <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.1f%%", x))
}

comparison_shared_feature_table_html <- function(x, note = NULL) {
  x[] <- lapply(x, as.character)

  rows <- vapply(seq_len(nrow(x)), function(i) {
    cells <- vapply(x[i, , drop = FALSE], function(value) {
      sprintf("<td>%s</td>", comparison_escape_html(value))
    }, character(1))
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  }, character(1))

  headers <- paste(sprintf("<th>%s</th>", comparison_escape_html(colnames(x))), collapse = "")
  note_html <- if (is.null(note) || !nzchar(note)) {
    ""
  } else {
    sprintf("<div class=\"shared-feature-note\">%s</div>", comparison_escape_html(note))
  }

  paste0(
    "<style>",
    ".shared-feature-wrap{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:100%;margin:2px 0 8px 0;}",
    ".shared-feature{border-collapse:collapse;width:auto;max-width:100%;font-size:12px;line-height:1.14;}",
    ".shared-feature th{background:#f4f6f8;color:#1f2933;text-align:left;border-bottom:1px solid #d8dee4;padding:4px 5px;white-space:nowrap;}",
    ".shared-feature td{border-bottom:1px solid #e5e7eb;padding:3px 5px;vertical-align:top;text-align:left;max-width:360px;overflow-wrap:anywhere;}",
    ".shared-feature-note{margin-top:4px;color:#667085;font-size:10px;}",
    "</style>",
    "<div class=\"shared-feature-wrap\">",
    "<table class=\"shared-feature\"><thead><tr>", headers, "</tr></thead><tbody>",
    paste(rows, collapse = ""),
    "</tbody></table>",
    note_html,
    "</div>"
  )
}

display_shared_feature_overview <- function(shared_feature_space, html = TRUE) {
  overview <- shared_feature_space$overview
  display <- overview[, c(
    "category", "attribute", "source", "status",
    "union_n", "shared_all_n", "dataset_totals", "interpretation"
  ), drop = FALSE]
  colnames(display) <- c(
    "category", "attribute", "source", "status",
    "union", "shared by all", "dataset totals", "interpretation"
  )

  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_shared_feature_table_html(
      display,
      note = sprintf("Dataset similarity class: %s", shared_feature_space$comparison_class$class_label)
    ))
  } else {
    print(display, row.names = FALSE)
  }
  invisible(overview)
}

display_shared_feature_relative_table <- function(shared_feature_space, category, html = TRUE) {
  category_key <- comparison_feature_key(category)
  rows <- shared_feature_space$relative_counts[
    shared_feature_space$relative_counts$category_key == category_key,
    ,
    drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    message(sprintf("%s: no relative overlap table is available.", comparison_feature_display_name(category_key)))
    return(invisible(rows))
  }

  display <- data.frame(
    dataset = rows$dataset,
    `shared of dataset` = comparison_format_percent(rows$shared_pct_dataset),
    `shared of union` = comparison_format_percent(rows$shared_pct_union),
    `only here of dataset` = comparison_format_percent(rows$only_pct_dataset),
    `only here of union` = comparison_format_percent(rows$only_pct_union),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  note <- "Relative values compare the Venn counts to each dataset total and to the union."

  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_shared_feature_table_html(display, note = note))
  } else {
    print(display, row.names = FALSE)
    cat(note, "\n")
  }
  invisible(rows)
}

display_shared_feature_skip_note <- function(shared_feature_space, category, html = TRUE) {
  category_key <- comparison_feature_key(category)
  row <- shared_feature_space$overview[shared_feature_space$overview$category_key == category_key, , drop = FALSE]
  if (nrow(row) == 0L) {
    return(invisible(NULL))
  }
  display <- row[, c("category", "status", "source", "interpretation"), drop = FALSE]

  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_shared_feature_table_html(display))
  } else {
    print(display, row.names = FALSE)
  }
  invisible(display)
}

comparison_venn_label_text <- function(label, total) {
  label <- gsub("\\s*\\|\\s*", "\n", label)
  sprintf("%s\nn = %s", label, total)
}

comparison_venn_set_label_positions <- function(plot_sets) {
  plot_data <- ggVennDiagram::process_data(ggVennDiagram::Venn(plot_sets))
  set_edge <- plot_data$setEdge
  ids <- unique(set_edge$id)

  rows <- lapply(seq_along(ids), function(i) {
    edge <- set_edge[set_edge$id == ids[i], , drop = FALSE]
    y_span <- diff(range(edge$Y, na.rm = TRUE))
    data.frame(
      id = ids[i],
      X = mean(edge$X, na.rm = TRUE),
      Y = mean(edge$Y, na.rm = TRUE),
      y_span = y_span,
      label = names(plot_sets)[i],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  middle_y <- stats::median(out$Y, na.rm = TRUE)
  out$Y <- ifelse(
    out$Y >= middle_y,
    out$Y + 0.28 * out$y_span,
    out$Y - 0.28 * out$y_span
  )
  out
}

plot_shared_feature_venn <- function(shared_feature_space, category, max_sets = 7L) {
  category_key <- comparison_feature_key(category)
  sets <- shared_feature_space$sets[[category_key]]
  status <- comparison_shared_feature_set_status(sets)

  if (identical(status, "skipped") || identical(status, "unavailable")) {
    comparison_warn(sprintf(
      "%s Venn skipped: category status is `%s`.",
      comparison_feature_display_name(category_key),
      status
    ))
    return(invisible(NULL))
  }
  if (length(sets) < 2L) {
    comparison_warn(sprintf("%s Venn requires at least two result directories.", comparison_feature_display_name(category_key)))
    return(invisible(NULL))
  }
  if (length(sets) > max_sets) {
    comparison_warn(sprintf(
      "%s Venn skipped: %d result directories selected; fixed-layout Venn plotting is limited to %d sets. Use the stored overlap tables instead.",
      comparison_feature_display_name(category_key),
      length(sets),
      max_sets
    ))
    return(invisible(NULL))
  }
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    comparison_warn("Package `ggVennDiagram` is required for shared feature Venn plots.")
    return(invisible(NULL))
  }

  totals <- vapply(sets, length, integer(1))
  plot_sets <- sets
  names(plot_sets) <- mapply(comparison_venn_label_text, names(sets), totals, USE.NAMES = FALSE)
  label_positions <- comparison_venn_set_label_positions(plot_sets)
  union_n <- length(comparison_feature_union(sets))
  shared_n <- length(comparison_feature_intersection(sets))

  ggVennDiagram::ggVennDiagram(
    plot_sets,
    category.names = rep("", length(plot_sets)),
    label = "count",
    label_geom = "text",
    label_alpha = 1,
    label_size = 4,
    label_color = "#1f2933"
  ) +
    ggplot2::geom_label(
      data = label_positions,
      ggplot2::aes(x = X, y = Y, label = label),
      inherit.aes = FALSE,
      linewidth = 0,
      fill = "white",
      alpha = 0.86,
      color = "#1f2933",
      size = 3.2,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_gradient(low = "#f8fbff", high = "#4c78a8") +
    ggplot2::labs(
      title = sprintf("%s overlap", comparison_feature_display_name(category_key)),
      subtitle = comparison_common_matrix_attribute_text(shared_feature_space$comparison_class)
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey35", hjust = 0),
      legend.position = "none"
    )
}

comparison_activity_matrix <- function(result, te_subfamilies = NULL, samples = NULL) {
  x <- result$tables$activity_scores
  if (anyDuplicated(x$te_subfamily)) {
    stop(sprintf("`activity_scores.tsv` in %s contains duplicated te_subfamily rows.", result$result_id), call. = FALSE)
  }

  rownames(x) <- x$te_subfamily
  x$te_subfamily <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "double"

  if (!is.null(te_subfamilies)) {
    mat <- mat[intersect(te_subfamilies, rownames(mat)), , drop = FALSE]
  }
  if (!is.null(samples)) {
    mat <- mat[, intersect(samples, colnames(mat)), drop = FALSE]
  }
  mat
}

comparison_condition_activity_matrix <- function(result, te_subfamilies = NULL, conditions = NULL) {
  mat <- comparison_activity_matrix(result, te_subfamilies = te_subfamilies)
  if (is.null(conditions) || length(conditions) == 0L) {
    return(mat[, FALSE, drop = FALSE])
  }

  condition_cols <- lapply(conditions, function(condition) {
    grep(condition, colnames(mat), value = TRUE, fixed = TRUE)
  })
  names(condition_cols) <- conditions
  condition_cols <- condition_cols[vapply(condition_cols, length, integer(1)) > 0L]

  if (length(condition_cols) == 0L) {
    return(mat[, FALSE, drop = FALSE])
  }

  out <- do.call(cbind, lapply(condition_cols, function(cols) {
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  }))
  colnames(out) <- names(condition_cols)
  rownames(out) <- rownames(mat)
  out
}

comparison_melt_matrix <- function(mat, row_name, col_name, value_name) {
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(out) <- c(row_name, col_name, value_name)
  out
}

comparison_plot_intersection_counts <- function(sets, title, item_label) {
  comparison_need_ggplot()
  shared <- Reduce(intersect, sets)
  df <- data.frame(
    result_id = names(sets),
    total = vapply(sets, length, integer(1)),
    shared = length(shared),
    stringsAsFactors = FALSE
  )
  df_long <- rbind(
    data.frame(result_id = df$result_id, count_type = "total", count = df$total, stringsAsFactors = FALSE),
    data.frame(result_id = df$result_id, count_type = "shared across all", count = df$shared, stringsAsFactors = FALSE)
  )
  ggplot2::ggplot(df_long, ggplot2::aes(x = result_id, y = count, fill = count_type)) +
    ggplot2::geom_col(position = "dodge", width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = sprintf("Number of %s", item_label), fill = NULL) +
    comparison_theme()
}

plot_shared_feature_space <- function(comparison_results) {
  shared_feature_space <- build_shared_feature_space(comparison_results)
  display_shared_feature_overview(shared_feature_space)

  for (category in c("samples", "gene_ids", "te_subfamilies")) {
    plot <- plot_shared_feature_venn(shared_feature_space, category)
    if (!is.null(plot)) {
      print(plot)
    }
    display_shared_feature_relative_table(shared_feature_space, category)
  }

  display_shared_feature_skip_note(shared_feature_space, "significant_genes")
  invisible(shared_feature_space)
}

summarise_activity_fit <- function(comparison_results) {
  out <- do.call(rbind, lapply(seq_along(comparison_results), function(i) {
    fit <- comparison_results[[i]]$tables$activity_fit_summary
    cbind(
      result_id = comparison_results[[i]]$result_id,
      fit,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarise_activity_scores <- function(comparison_results, target_tes = NULL) {
  shared_tes <- get_shared_te_subfamilies(comparison_results)
  shared_samples <- get_shared_samples(comparison_results)
  shared_conditions <- get_shared_conditions(comparison_results)

  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }

  if (length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared TE subfamilies are available across all selected result directories.",
      shared_te_subfamilies = shared_tes,
      shared_samples = shared_samples
    ))
  }

  comparison_space <- "shared samples"
  axis_label <- "Sample"
  if (length(shared_samples) > 0L) {
    matrices <- lapply(comparison_results, comparison_activity_matrix, te_subfamilies = shared_tes, samples = shared_samples)
  } else if (length(shared_conditions) > 0L) {
    matrices <- lapply(comparison_results, comparison_condition_activity_matrix, te_subfamilies = shared_tes, conditions = shared_conditions)
    shared_samples <- Reduce(intersect, lapply(matrices, colnames))
    matrices <- lapply(matrices, function(mat) mat[, shared_samples, drop = FALSE])
    comparison_space <- "shared condition means"
    axis_label <- "Condition"
  } else {
    return(list(
      status = "skipped",
      reason = "No shared samples or shared conditions are available across all selected result directories.",
      shared_te_subfamilies = shared_tes,
      shared_samples = character(),
      shared_conditions = shared_conditions
    ))
  }

  if (length(shared_samples) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared sample or condition columns remain after activity-score alignment.",
      shared_te_subfamilies = shared_tes,
      shared_samples = character(),
      shared_conditions = shared_conditions
    ))
  }

  names(matrices) <- comparison_result_ids(comparison_results)

  vectorized <- lapply(matrices, as.vector)
  activity_correlation <- stats::cor(do.call(cbind, vectorized), use = "pairwise.complete.obs")

  total_abs_activity <- do.call(rbind, lapply(names(matrices), function(result_id) {
    data.frame(
      result_id = result_id,
      te_subfamily = rownames(matrices[[result_id]]),
      total_abs_activity = rowSums(abs(matrices[[result_id]]), na.rm = TRUE),
      zero_fraction = rowMeans(abs(matrices[[result_id]]) <= .Machine$double.eps, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  heatmap_data <- do.call(rbind, lapply(names(matrices), function(result_id) {
    x <- comparison_melt_matrix(matrices[[result_id]], "te_subfamily", "sample_id", "activity_score")
    x$result_id <- result_id
    x
  }))

  list(
    status = "ok",
    comparison_space = comparison_space,
    axis_label = axis_label,
    shared_te_subfamilies = shared_tes,
    shared_samples = shared_samples,
    shared_conditions = shared_conditions,
    activity_correlation = activity_correlation,
    total_abs_activity = total_abs_activity,
    heatmap_data = heatmap_data
  )
}

comparison_te_majority_classification <- function(comparison_results, te_subfamilies, class_col = "trend_class") {
  votes <- do.call(rbind, lapply(comparison_results, function(result) {
    contrast <- result$tables$te_contrast_summary
    if (is.null(contrast) || !(class_col %in% colnames(contrast)) || !("te_subfamily" %in% colnames(contrast))) {
      return(NULL)
    }
    x <- contrast[contrast$te_subfamily %in% te_subfamilies, c("te_subfamily", class_col), drop = FALSE]
    colnames(x)[2] <- "te_class"
    x$te_class <- as.character(x$te_class)
    x <- x[!is.na(x$te_class) & nzchar(x$te_class), , drop = FALSE]
    if (nrow(x) == 0L) {
      return(NULL)
    }
    x$result_id <- result$result_id
    x
  }))

  if (is.null(votes) || nrow(votes) == 0L) {
    return(data.frame(
      te_subfamily = te_subfamilies,
      majority_class = NA_character_,
      n_class_votes = 0L,
      stringsAsFactors = FALSE
    ))
  }

  majority <- do.call(rbind, lapply(te_subfamilies, function(te) {
    te_votes <- votes$te_class[votes$te_subfamily == te]
    if (length(te_votes) == 0L) {
      return(data.frame(
        te_subfamily = te,
        majority_class = NA_character_,
        n_class_votes = 0L,
        stringsAsFactors = FALSE
      ))
    }
    counts <- sort(table(te_votes), decreasing = TRUE)
    top_count <- max(as.integer(counts))
    top_classes <- sort(names(counts)[as.integer(counts) == top_count])
    data.frame(
      te_subfamily = te,
      majority_class = top_classes[1],
      n_class_votes = top_count,
      stringsAsFactors = FALSE
    )
  }))
  rownames(majority) <- NULL
  majority
}

comparison_order_activity_heatmap_tes <- function(plot_data, comparison_results, target_tes = NULL, max_tes = 30L) {
  totals <- aggregate(abs(activity_score) ~ te_subfamily, data = plot_data, FUN = sum)
  colnames(totals)[2] <- "total_abs_activity"

  if (is.null(target_tes)) {
    keep <- head(totals$te_subfamily[order(totals$total_abs_activity, decreasing = TRUE)], max_tes)
    return(list(
      te_order = keep,
      class_data = NULL,
      has_classification = FALSE
    ))
  }

  class_data <- comparison_te_majority_classification(
    comparison_results = comparison_results,
    te_subfamilies = unique(plot_data$te_subfamily)
  )
  class_data <- merge(class_data, totals, by = "te_subfamily", all.x = TRUE, sort = FALSE)
  class_data$total_abs_activity[is.na(class_data$total_abs_activity)] <- 0
  class_data$majority_class_display <- class_data$majority_class
  class_data$majority_class_display[is.na(class_data$majority_class_display) | !nzchar(class_data$majority_class_display)] <- "unclassified"

  class_totals <- aggregate(total_abs_activity ~ majority_class_display, data = class_data, FUN = sum)
  class_totals$class_rank <- order(order(-class_totals$total_abs_activity, class_totals$majority_class_display))
  class_data <- merge(
    class_data,
    class_totals[, c("majority_class_display", "class_rank"), drop = FALSE],
    by = "majority_class_display",
    all.x = TRUE,
    sort = FALSE
  )
  class_data <- class_data[order(class_data$class_rank, -class_data$total_abs_activity, class_data$te_subfamily), , drop = FALSE]

  list(
    te_order = class_data$te_subfamily,
    class_data = class_data[, c("te_subfamily", "majority_class_display"), drop = FALSE],
    has_classification = any(class_data$majority_class_display != "unclassified")
  )
}

comparison_activity_matrices_shared_tes <- function(comparison_results, target_tes = NULL) {
  shared_tes <- get_shared_te_subfamilies(comparison_results)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared TE subfamilies are available across all selected result directories.",
      shared_te_subfamilies = shared_tes
    ))
  }

  matrices <- lapply(comparison_results, comparison_activity_matrix, te_subfamilies = shared_tes)
  names(matrices) <- comparison_result_ids(comparison_results)
  matrices <- lapply(matrices, function(mat) mat[shared_tes, , drop = FALSE])

  list(
    status = "ok",
    shared_te_subfamilies = shared_tes,
    matrices = matrices,
    shared_samples = Reduce(intersect, lapply(matrices, colnames))
  )
}

summarise_sparsity <- function(comparison_results) {
  fit <- summarise_activity_fit(comparison_results)
  fit_sparsity <- if ("mean_sparsity_percent" %in% colnames(fit)) {
    fit[, c("result_id", "scope", "condition", "mean_sparsity_percent"), drop = FALSE]
  } else {
    data.frame()
  }

  activity_sparsity <- do.call(rbind, lapply(comparison_results, function(result) {
    mat <- comparison_activity_matrix(result)
    data.frame(
      result_id = result$result_id,
      n_values = length(mat),
      n_near_zero = sum(abs(mat) <= .Machine$double.eps, na.rm = TRUE),
      n_active = sum(abs(mat) > .Machine$double.eps, na.rm = TRUE),
      zero_fraction_from_activity_scores = mean(abs(mat) <= .Machine$double.eps, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(
    fit_sparsity = fit_sparsity,
    activity_sparsity = activity_sparsity
  )
}

summarise_lambda_best <- function(comparison_results) {
  rows <- lapply(comparison_results, function(result) {
    sample_diag <- result$tables$sample_activity_diagnostics
    lambda_diag <- result$tables$activity_lambda_diagnostics
    fit <- result$tables$activity_fit_summary

    source_table <- NA_character_
    lambda_column <- NA_character_
    values <- numeric()

    if (!is.null(sample_diag) && "lambda_best" %in% colnames(sample_diag)) {
      values <- as.numeric(sample_diag$lambda_best)
      source_table <- "sample_activity_diagnostics.tsv"
      lambda_column <- "lambda_best"
    } else {
      lambda_cols <- if (is.null(fit)) character() else grep("lambda", colnames(fit), value = TRUE, ignore.case = TRUE)
      if (length(lambda_cols) > 0L) {
        values <- unlist(fit[lambda_cols], use.names = FALSE)
        source_table <- "activity_fit_summary.tsv"
        lambda_column <- paste(lambda_cols, collapse = ", ")
      }
    }
    values <- values[!is.na(values)]

    global_lambda <- NA_real_
    global_source <- NA_character_
    if (!is.null(lambda_diag)) {
      if ("global_lambda_min" %in% colnames(lambda_diag)) {
        global_values <- as.numeric(lambda_diag$global_lambda_min)
        global_values <- global_values[!is.na(global_values)]
        if (length(global_values) > 0L) {
          global_lambda <- global_values[1]
          global_source <- "activity_lambda_diagnostics.tsv"
        }
      }
      if (is.na(global_lambda) && all(c("lambda", "is_global_lambda_min") %in% colnames(lambda_diag))) {
        marked <- lambda_diag[isTRUE(lambda_diag$is_global_lambda_min) | lambda_diag$is_global_lambda_min %in% TRUE, , drop = FALSE]
        if (nrow(marked) > 0L) {
          global_lambda <- as.numeric(marked$lambda[1])
          global_source <- "activity_lambda_diagnostics.tsv"
        }
      }
      if (is.na(global_lambda) && all(c("lambda", "mean_crossval_mse") %in% colnames(lambda_diag))) {
        mse <- as.numeric(lambda_diag$mean_crossval_mse)
        if (any(!is.na(mse))) {
          global_lambda <- as.numeric(lambda_diag$lambda[which.min(mse)])
          global_source <- "activity_lambda_diagnostics.tsv"
        }
      }
    }
    if (is.na(global_lambda) && !is.null(sample_diag) && "lambda_min_overall" %in% colnames(sample_diag)) {
      global_values <- as.numeric(sample_diag$lambda_min_overall)
      global_values <- global_values[!is.na(global_values)]
      if (length(global_values) > 0L) {
        global_lambda <- global_values[1]
        global_source <- "sample_activity_diagnostics.tsv"
      }
    }

    data.frame(
      result_id = result$result_id,
      lambda_best_available = length(values) > 0L,
      global_lambda_available = !is.na(global_lambda),
      source_table = source_table,
      lambda_column = lambda_column,
      global_lambda_source = global_source,
      n_values = length(values),
      mean_lambda_best = if (length(values) > 0L) mean(values) else NA_real_,
      median_lambda_best = if (length(values) > 0L) stats::median(values) else NA_real_,
      min_lambda_best = if (length(values) > 0L) min(values) else NA_real_,
      max_lambda_best = if (length(values) > 0L) max(values) else NA_real_,
      global_lambda_min = global_lambda,
      note = if (length(values) > 0L) NA_character_ else "No sample-specific lambda_best column is available in loaded comparison TSVs.",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!any(out$lambda_best_available)) {
    comparison_warn("No sample-specific `lambda_best` column is available in the loaded comparison TSVs.")
  }
  out
}

summarise_signif_model_diagnostics <- function(comparison_results) {
  diagnostics <- do.call(rbind, lapply(comparison_results, function(result) {
    x <- result$tables$signif_model_diagnostics
    if (is.null(x) || nrow(x) == 0L) {
      return(NULL)
    }
    x$result_id <- result$result_id
    x
  }))

  if (is.null(diagnostics) || nrow(diagnostics) == 0L) {
    return(list(
      status = "skipped",
      reason = "`signif_model_diagnostics.tsv` is not available in any selected result directory.",
      diagnostics = data.frame()
    ))
  }

  result_ids <- comparison_result_ids(comparison_results)
  contrast_sets <- lapply(comparison_results, function(result) {
    x <- result$tables$signif_model_diagnostics
    if (is.null(x) || !"condition" %in% colnames(x)) {
      return(character())
    }
    unique(as.character(x$condition))
  })
  contrast_sets <- contrast_sets[vapply(contrast_sets, length, integer(1)) > 0L]
  shared_contrasts <- if (length(contrast_sets) == 0L) character() else Reduce(intersect, contrast_sets)

  if (length(shared_contrasts) > 0L) {
    diagnostics <- diagnostics[as.character(diagnostics$condition) %in% shared_contrasts, , drop = FALSE]
  }

  metric_cols <- intersect(c("r_squared", "rss", "f_statistic", "n_aliased"), colnames(diagnostics))
  summary_table <- if (length(metric_cols) == 0L) {
    data.frame()
  } else {
    aggregate(
      diagnostics[, metric_cols, drop = FALSE],
      by = list(result_id = diagnostics$result_id),
      FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
    )
  }

  list(
    status = "ok",
    diagnostics = diagnostics,
    shared_contrasts = shared_contrasts,
    summary_table = summary_table,
    result_ids = result_ids
  )
}

plot_signif_model_diagnostics_comparison <- function(comparison_results,
                                                     metrics = c("r_squared", "rss", "f_statistic", "n_aliased"),
                                                     comparison_class = NULL) {
  comparison_need_ggplot()
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  diagnostic_summary <- summarise_signif_model_diagnostics(comparison_results)
  if (!identical(diagnostic_summary$status, "ok")) {
    comparison_warn(diagnostic_summary$reason)
    return(invisible(NULL))
  }

  diagnostics <- diagnostic_summary$diagnostics
  metrics <- intersect(metrics, colnames(diagnostics))
  if (length(metrics) == 0L) {
    comparison_warn("No supported numeric columns are available in `signif_model_diagnostics.tsv`.")
    return(invisible(NULL))
  }

  result_ids <- comparison_result_ids(comparison_results)
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  diagnostics$result_label <- unname(label_map[diagnostics$result_id])
  diagnostics$result_label[is.na(diagnostics$result_label) | !nzchar(diagnostics$result_label)] <- comparison_short_result_label(diagnostics$result_id[is.na(diagnostics$result_label) | !nzchar(diagnostics$result_label)])
  result_labels <- unname(label_map[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- comparison_short_result_label(result_ids[is.na(result_labels) | !nzchar(result_labels)])
  diagnostics$result_label <- factor(diagnostics$result_label, levels = result_labels)

  contrast_levels <- diagnostic_summary$shared_contrasts
  if (length(contrast_levels) == 0L) {
    contrast_levels <- unique(as.character(diagnostics$condition))
  }
  diagnostics$condition <- factor(as.character(diagnostics$condition), levels = contrast_levels)

  metric_labels <- c(
    r_squared = "R2",
    rss = "RSS",
    f_statistic = "F statistic",
    n_aliased = "Aliased predictors"
  )
  plot_data <- do.call(rbind, lapply(metrics, function(metric) {
    data.frame(
      result_id = diagnostics$result_id,
      result_label = diagnostics$result_label,
      condition = diagnostics$condition,
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      value = as.numeric(diagnostics[[metric]]),
      stringsAsFactors = FALSE
    )
  }))
  plot_data$metric_label[is.na(plot_data$metric_label) | !nzchar(plot_data$metric_label)] <- comparison_pretty_label(plot_data$metric[is.na(plot_data$metric_label) | !nzchar(plot_data$metric_label)])

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = condition, y = value, color = result_label, group = result_label)
  ) +
    ggplot2::geom_line(linewidth = 0.65, alpha = 0.85) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = "Contrast model-fit diagnostics",
      subtitle = comparison_plot_subtitle(
        "Loaded from signif_model_diagnostics.tsv for each pairwise contrast",
        comparison_class
      ),
      x = "Contrast",
      y = NULL,
      color = "Result"
    ) +
    comparison_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}

comparison_contrast_padj_columns <- function(common_cols) {
  grep("^padj_", common_cols, value = TRUE)
}

comparison_neg_log10_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x) | x < 0 | x > 1] <- NA_real_
  -log10(pmax(x, .Machine$double.xmin, na.rm = FALSE))
}

comparison_row_mean <- function(mat) {
  apply(mat, 1, function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  })
}

comparison_row_min <- function(mat) {
  apply(mat, 1, function(x) {
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  })
}

comparison_row_max <- function(mat) {
  apply(mat, 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
}

comparison_contrast_raw_p_columns <- function(common_cols) {
  grep("^p_value_", common_cols, value = TRUE)
}

comparison_contrast_delta_columns <- function(common_cols) {
  grep("^delta_", common_cols, value = TRUE)
}

comparison_non_noise_trend_classes <- function() {
  te_senescence_non_noise_trend_classes()
}

comparison_noise_classes <- function() {
  te_senescence_noise_classes()
}

comparison_trend_class_colors <- function(classes = NULL, include_noise = FALSE) {
  te_senescence_trend_class_colors(classes = classes, include_noise = include_noise)
}

comparison_is_noise_class <- function(x) {
  te_senescence_is_noise_class(x)
}

comparison_trend_class_delta_signs <- function(trend_class) {
  te_senescence_trend_delta_signs(trend_class)
}

comparison_trend_class_p_contrasts <- function(trend_class) {
  te_senescence_trend_p_contrasts(trend_class)
}

comparison_trend_class_raw_p <- function(x, trend_class, raw_p_cols) {
  contrasts <- comparison_trend_class_p_contrasts(trend_class)
  if (length(contrasts) == 0L || length(raw_p_cols) == 0L) {
    return(numeric())
  }
  cols <- intersect(paste0("p_value_", contrasts), raw_p_cols)
  suppressWarnings(as.numeric(unlist(x[cols], use.names = FALSE)))
}

comparison_trend_class_fisher_p <- function(x, trend_class, raw_p_cols) {
  comparison_fisher_p(comparison_trend_class_raw_p(x, trend_class, raw_p_cols))
}

comparison_class_direction_ok <- function(x, trend_class, delta_cols) {
  signs <- comparison_trend_class_delta_signs(trend_class)
  if (length(signs) == 0L || length(delta_cols) == 0L) {
    return(FALSE)
  }
  cols <- paste0("delta_", names(signs))
  if (!all(cols %in% delta_cols) || !all(cols %in% colnames(x))) {
    return(FALSE)
  }
  deltas <- suppressWarnings(as.numeric(unlist(x[cols], use.names = FALSE)))
  all(!is.na(deltas) & sign(deltas) == unname(signs))
}

comparison_adjust_finite_p <- function(x, method = "BH") {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  keep <- !is.na(x) & x >= 0 & x <= 1
  if (any(keep)) {
    out[keep] <- stats::p.adjust(x[keep], method = method)
  }
  out
}

comparison_contrast_aggregated_candidate_table <- function(result_summary, raw_p_cols, delta_cols) {
  if (length(raw_p_cols) == 0L || length(delta_cols) == 0L || nrow(result_summary) == 0L) {
    return(data.frame())
  }
  do.call(rbind, lapply(seq_len(nrow(result_summary)), function(i) {
    row <- result_summary[i, , drop = FALSE]
    do.call(rbind, lapply(comparison_non_noise_trend_classes(), function(class_name) {
      direction_ok <- comparison_class_direction_ok(row, class_name, delta_cols)
      data.frame(
        result_id = row$result_id,
        result_label = row$result_label,
        te_subfamily = row$te_subfamily,
        contrast_aggregated_class_candidate = class_name,
        contrast_aggregated_direction_ok = direction_ok,
        contrast_aggregated_raw_p_fisher = if (isTRUE(direction_ok)) {
          comparison_trend_class_fisher_p(row, class_name, raw_p_cols)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }))
  }))
}

comparison_best_contrast_aggregated_class <- function(candidate_table, result_summary, alpha) {
  if (nrow(result_summary) == 0L) {
    return(data.frame())
  }
  empty_best <- result_summary[, c("result_id", "result_label", "te_subfamily"), drop = FALSE]
  empty_best$contrast_aggregated_class <- "flat_or_noisy"
  empty_best$contrast_aggregated_raw_p_fisher <- NA_real_
  empty_best$contrast_aggregated_raw_p_fisher_bh_evidence <- NA_real_

  if (is.null(candidate_table) || nrow(candidate_table) == 0L) {
    return(empty_best)
  }
  candidate_table$contrast_aggregated_raw_p_fisher_bh_evidence <- comparison_adjust_finite_p(
    candidate_table$contrast_aggregated_raw_p_fisher
  )
  key <- paste(candidate_table$result_id, candidate_table$te_subfamily, sep = "\r")
  best <- do.call(rbind, lapply(split(candidate_table, key), function(x) {
    valid <- x[
      !is.na(x$contrast_aggregated_raw_p_fisher_bh_evidence) &
        x$contrast_aggregated_raw_p_fisher_bh_evidence <= alpha,
      ,
      drop = FALSE
    ]
    if (nrow(valid) == 0L) {
      return(data.frame(
        result_id = x$result_id[1],
        result_label = x$result_label[1],
        te_subfamily = x$te_subfamily[1],
        contrast_aggregated_class = "flat_or_noisy",
        contrast_aggregated_raw_p_fisher = NA_real_,
        contrast_aggregated_raw_p_fisher_bh_evidence = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    valid <- valid[order(valid$contrast_aggregated_raw_p_fisher_bh_evidence, valid$contrast_aggregated_class_candidate), , drop = FALSE]
    data.frame(
      result_id = valid$result_id[1],
      result_label = valid$result_label[1],
      te_subfamily = valid$te_subfamily[1],
      contrast_aggregated_class = valid$contrast_aggregated_class_candidate[1],
      contrast_aggregated_raw_p_fisher = valid$contrast_aggregated_raw_p_fisher[1],
      contrast_aggregated_raw_p_fisher_bh_evidence = valid$contrast_aggregated_raw_p_fisher_bh_evidence[1],
      stringsAsFactors = FALSE
    )
  }))
  rownames(best) <- NULL

  out_key <- paste(empty_best$result_id, empty_best$te_subfamily, sep = "\r")
  best_key <- paste(best$result_id, best$te_subfamily, sep = "\r")
  idx <- match(out_key, best_key)
  has_best <- !is.na(idx)
  empty_best[has_best, ] <- best[idx[has_best], colnames(empty_best), drop = FALSE]
  empty_best
}

comparison_global_class_assignment <- function(x, raw_p_cols) {
  classes <- as.character(x$trend_class)
  classes[is.na(classes) | !nzchar(classes)] <- "unclassified"
  non_noise <- classes[!comparison_is_noise_class(classes)]

  if (length(non_noise) == 0L) {
    class_rows <- x[classes == "flat_or_noisy", , drop = FALSE]
    if (nrow(class_rows) == 0L) {
      class_rows <- x
    }
    global_raw_p <- unlist(lapply(seq_len(nrow(class_rows)), function(i) {
      comparison_trend_class_raw_p(class_rows[i, , drop = FALSE], "flat_or_noisy", raw_p_cols)
    }), use.names = FALSE)
    return(data.frame(
      majority_class = "flat_or_noisy",
      global_assigned_class = "flat_or_noisy",
      global_class_votes = sum(classes == "flat_or_noisy"),
      global_non_noise_votes = 0L,
      global_class_tie = FALSE,
      global_assignment_rule = "all result dirs are flat_or_noisy or unclassified",
      global_combined_raw_p_fisher = comparison_fisher_p(global_raw_p),
      stringsAsFactors = FALSE
    ))
  }

  counts <- sort(table(non_noise), decreasing = TRUE)
  top_count <- max(as.integer(counts))
  top_classes <- sort(names(counts)[as.integer(counts) == top_count])
  class_p <- stats::setNames(rep(NA_real_, length(top_classes)), top_classes)
  for (class_name in top_classes) {
    class_rows <- x[classes == class_name, , drop = FALSE]
    class_raw_p <- unlist(lapply(seq_len(nrow(class_rows)), function(i) {
      comparison_trend_class_raw_p(class_rows[i, , drop = FALSE], class_name, raw_p_cols)
    }), use.names = FALSE)
    class_p[[class_name]] <- comparison_fisher_p(class_raw_p)
  }
  if (all(is.na(class_p))) {
    assigned_class <- top_classes[1]
  } else {
    assigned_class <- names(class_p)[order(class_p, names(class_p), na.last = TRUE)][1]
  }
  data.frame(
    majority_class = assigned_class,
    global_assigned_class = assigned_class,
    global_class_votes = top_count,
    global_non_noise_votes = length(non_noise),
    global_class_tie = length(top_classes) > 1L,
    global_assignment_rule = if (length(top_classes) > 1L) {
      "non-noise majority tie; selected lower same-class Fisher p-value"
    } else {
      "non-noise majority class"
    },
    global_combined_raw_p_fisher = unname(class_p[[assigned_class]]),
    stringsAsFactors = FALSE
  )
}

summarise_te_padj_evidence <- function(comparison_results,
                                       target_tes = NULL,
                                       p_score_mode = c("mean_padj", "mean_neg_log10_padj"),
                                       alpha = 0.05,
                                       comparison_class = NULL) {
  p_score_mode <- match.arg(p_score_mode)
  score_label <- if (identical(p_score_mode, "mean_neg_log10_padj")) {
    "mean -log10(adjusted p-value)"
  } else {
    "mean adjusted p-value"
  }
  lower_is_stronger <- identical(p_score_mode, "mean_padj")
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  has_contrast <- vapply(comparison_results, function(result) !is.null(result$tables$te_contrast_summary), logical(1))
  if (!all(has_contrast)) {
    missing_ids <- comparison_result_ids(comparison_results)[!has_contrast]
    return(list(
      status = "skipped",
      reason = sprintf("Missing te_contrast_summary.tsv in: %s", paste(missing_ids, collapse = ", "))
    ))
  }

  contrast_sets <- lapply(comparison_results, function(result) unique(result$tables$te_contrast_summary$te_subfamily))
  shared_tes <- Reduce(intersect, contrast_sets)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(status = "skipped", reason = "No shared TE subfamilies are available in te_contrast_summary.tsv."))
  }

  common_cols <- Reduce(intersect, lapply(comparison_results, function(result) colnames(result$tables$te_contrast_summary)))
  padj_cols <- comparison_contrast_padj_columns(common_cols)
  if (length(padj_cols) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared pairwise adjusted p-value columns matching `padj_*` are available in te_contrast_summary.tsv."
    ))
  }
  raw_p_cols <- comparison_contrast_raw_p_columns(common_cols)
  delta_cols <- comparison_contrast_delta_columns(common_cols)
  class_col <- if ("trend_class" %in% common_cols) "trend_class" else NA_character_

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  result_summary <- do.call(rbind, lapply(comparison_results, function(result) {
    x <- result$tables$te_contrast_summary
    x <- x[x$te_subfamily %in% shared_tes, , drop = FALSE]
    x <- x[match(shared_tes, x$te_subfamily), , drop = FALSE]
    p_mat <- as.matrix(x[, padj_cols, drop = FALSE])
    storage.mode(p_mat) <- "double"
    score_mat <- if (identical(p_score_mode, "mean_neg_log10_padj")) comparison_neg_log10_p(p_mat) else p_mat
    dim(score_mat) <- dim(p_mat)
    dimnames(score_mat) <- dimnames(p_mat)

    out <- data.frame(
      result_id = result$result_id,
      te_subfamily = x$te_subfamily,
      p_score = comparison_row_mean(score_mat),
      min_p_score = comparison_row_min(score_mat),
      max_p_score = comparison_row_max(score_mat),
      n_p_values = rowSums(!is.na(score_mat)),
      n_significant_contrasts = rowSums(p_mat <= alpha, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    out$any_significant <- out$n_significant_contrasts > 0L
    out$trend_class <- if (!is.na(class_col)) as.character(x[[class_col]]) else NA_character_
    if (length(raw_p_cols) > 0L) {
      out[raw_p_cols] <- x[, raw_p_cols, drop = FALSE]
      out$combined_raw_p_fisher <- vapply(
        seq_len(nrow(x)),
        function(i) comparison_trend_class_fisher_p(x[i, , drop = FALSE], out$trend_class[i], raw_p_cols),
        numeric(1)
      )
    } else {
      out$combined_raw_p_fisher <- NA_real_
    }
    if (length(delta_cols) > 0L) {
      out[delta_cols] <- x[, delta_cols, drop = FALSE]
    }
    out$result_label <- unname(label_map[out$result_id])
    out$result_label[is.na(out$result_label) | !nzchar(out$result_label)] <- out$result_id[is.na(out$result_label) | !nzchar(out$result_label)]
    out
  }))
  rownames(result_summary) <- NULL
  result_summary$combined_raw_p_fisher_bh_evidence <- comparison_adjust_finite_p(result_summary$combined_raw_p_fisher)

  contrast_candidate_summary <- comparison_contrast_aggregated_candidate_table(result_summary, raw_p_cols, delta_cols)
  contrast_class_summary <- comparison_best_contrast_aggregated_class(contrast_candidate_summary, result_summary, alpha)
  result_key <- paste(result_summary$result_id, result_summary$te_subfamily, sep = "\r")
  contrast_key <- paste(contrast_class_summary$result_id, contrast_class_summary$te_subfamily, sep = "\r")
  contrast_idx <- match(result_key, contrast_key)
  result_summary$contrast_aggregated_class <- contrast_class_summary$contrast_aggregated_class[contrast_idx]
  result_summary$contrast_aggregated_raw_p_fisher <- contrast_class_summary$contrast_aggregated_raw_p_fisher[contrast_idx]
  result_summary$contrast_aggregated_raw_p_fisher_bh_evidence <- contrast_class_summary$contrast_aggregated_raw_p_fisher_bh_evidence[contrast_idx]
  result_summary$activity_score_class <- result_summary$trend_class
  result_summary$activity_score_class[is.na(result_summary$activity_score_class) | !nzchar(result_summary$activity_score_class)] <- "unclassified"
  result_summary$activity_vs_contrast_class_match <- result_summary$activity_score_class == result_summary$contrast_aggregated_class

  te_summary <- do.call(rbind, lapply(split(result_summary, result_summary$te_subfamily), function(x) {
    global_assignment <- comparison_global_class_assignment(x, raw_p_cols)
    class_values <- x$trend_class[!is.na(x$trend_class) & nzchar(x$trend_class)]
    class_values_no_noise <- class_values[!comparison_is_noise_class(class_values)]
    classification_agreement <- if (length(class_values_no_noise) > 0L) {
      global_assignment$global_class_votes[1] / length(class_values_no_noise)
    } else if (length(class_values) > 0L) {
      sum(class_values == "flat_or_noisy") / length(class_values)
    } else {
      NA_real_
    }
    data.frame(
      te_subfamily = x$te_subfamily[1],
      mean_p_score = if (all(is.na(x$p_score))) NA_real_ else mean(x$p_score, na.rm = TRUE),
      min_result_p_score = if (all(is.na(x$p_score))) NA_real_ else min(x$p_score, na.rm = TRUE),
      max_result_p_score = if (all(is.na(x$p_score))) NA_real_ else max(x$p_score, na.rm = TRUE),
      p_score_range = if (all(is.na(x$p_score))) NA_real_ else max(x$p_score, na.rm = TRUE) - min(x$p_score, na.rm = TRUE),
      n_results = length(unique(x$result_id)),
      n_results_significant = sum(x$any_significant, na.rm = TRUE),
      any_result_significant = any(x$any_significant, na.rm = TRUE),
      all_results_significant = all(x$any_significant),
      majority_class = global_assignment$majority_class[1],
      global_assigned_class = global_assignment$global_assigned_class[1],
      global_class_votes = global_assignment$global_class_votes[1],
      global_non_noise_votes = global_assignment$global_non_noise_votes[1],
      global_class_tie = global_assignment$global_class_tie[1],
      global_assignment_rule = global_assignment$global_assignment_rule[1],
      classification_agreement = classification_agreement,
      global_combined_raw_p_fisher = global_assignment$global_combined_raw_p_fisher[1],
      stringsAsFactors = FALSE
    )
  }))
  rownames(te_summary) <- NULL
  te_summary$global_combined_raw_p_fisher_bh_evidence <- comparison_adjust_finite_p(te_summary$global_combined_raw_p_fisher)
  te_summary <- if (isTRUE(lower_is_stronger)) {
    te_summary[order(te_summary$mean_p_score, -te_summary$p_score_range, te_summary$te_subfamily), , drop = FALSE]
  } else {
    te_summary[order(-te_summary$mean_p_score, -te_summary$p_score_range, te_summary$te_subfamily), , drop = FALSE]
  }

  list(
    status = "ok",
    shared_te_subfamilies = shared_tes,
    p_score_mode = p_score_mode,
    alpha = alpha,
    p_columns = padj_cols,
    raw_p_columns = raw_p_cols,
    delta_columns = delta_cols,
    score_label = score_label,
    lower_is_stronger = lower_is_stronger,
    result_summary = result_summary,
    te_summary = te_summary,
    contrast_candidate_summary = contrast_candidate_summary,
    contrast_class_summary = contrast_class_summary
  )
}

summarise_activity_vs_contrast_classification <- function(te_padj_evidence) {
  if (!is.list(te_padj_evidence) || !identical(te_padj_evidence$status, "ok")) {
    return(data.frame())
  }
  result_summary <- te_padj_evidence$result_summary
  required <- c("result_label", "te_subfamily", "trend_class", "contrast_aggregated_class")
  if (is.null(result_summary) || nrow(result_summary) == 0L || !all(required %in% colnames(result_summary))) {
    return(data.frame())
  }

  activity_class <- as.character(result_summary$trend_class)
  activity_class[is.na(activity_class) | !nzchar(activity_class)] <- "unclassified"
  contrast_class <- as.character(result_summary$contrast_aggregated_class)
  contrast_class[is.na(contrast_class) | !nzchar(contrast_class)] <- "flat_or_noisy"
  result_summary$activity_score_class <- activity_class
  result_summary$contrast_aggregated_class <- contrast_class
  result_summary$activity_vs_contrast_class_match <- activity_class == contrast_class

  labels <- unique(as.character(result_summary$result_label))
  rows <- lapply(labels, function(label) {
    x <- result_summary[as.character(result_summary$result_label) == label, , drop = FALSE]
    activity_non_noise <- !comparison_is_noise_class(x$activity_score_class)
    contrast_non_noise <- !comparison_is_noise_class(x$contrast_aggregated_class)
    both_non_noise <- activity_non_noise & contrast_non_noise
    data.frame(
      result_label = label,
      n_te_subfamilies = nrow(x),
      n_matching_classes = sum(x$activity_vs_contrast_class_match, na.rm = TRUE),
      fraction_matching = if (nrow(x) == 0L) NA_real_ else mean(x$activity_vs_contrast_class_match, na.rm = TRUE),
      n_activity_score_non_noise = sum(activity_non_noise, na.rm = TRUE),
      n_contrast_aggregated_non_noise = sum(contrast_non_noise, na.rm = TRUE),
      n_both_non_noise_matching = sum(x$activity_vs_contrast_class_match & both_non_noise, na.rm = TRUE),
      n_both_non_noise_discordant = sum(!x$activity_vs_contrast_class_match & both_non_noise, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

activity_vs_contrast_classification_examples <- function(te_padj_evidence, max_rows = 20L) {
  if (!is.list(te_padj_evidence) || !identical(te_padj_evidence$status, "ok")) {
    return(data.frame())
  }
  result_summary <- te_padj_evidence$result_summary
  required <- c("result_label", "te_subfamily", "trend_class", "contrast_aggregated_class")
  if (is.null(result_summary) || nrow(result_summary) == 0L || !all(required %in% colnames(result_summary))) {
    return(data.frame())
  }

  activity_class <- as.character(result_summary$trend_class)
  activity_class[is.na(activity_class) | !nzchar(activity_class)] <- "unclassified"
  contrast_class <- as.character(result_summary$contrast_aggregated_class)
  contrast_class[is.na(contrast_class) | !nzchar(contrast_class)] <- "flat_or_noisy"
  result_summary$activity_score_class <- activity_class
  result_summary$contrast_aggregated_class <- contrast_class
  result_summary$activity_vs_contrast_class_match <- activity_class == contrast_class

  out <- result_summary[!result_summary$activity_vs_contrast_class_match, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(data.frame())
  }
  activity_p <- if ("combined_raw_p_fisher_bh_evidence" %in% colnames(out)) out$combined_raw_p_fisher_bh_evidence else rep(NA_real_, nrow(out))
  contrast_p <- if ("contrast_aggregated_raw_p_fisher_bh_evidence" %in% colnames(out)) out$contrast_aggregated_raw_p_fisher_bh_evidence else rep(NA_real_, nrow(out))
  best_p <- rep(Inf, nrow(out))
  for (i in seq_len(nrow(out))) {
    vals <- suppressWarnings(as.numeric(c(activity_p[i], contrast_p[i])))
    vals <- vals[!is.na(vals) & vals >= 0 & vals <= 1]
    if (length(vals) > 0L) {
      best_p[i] <- min(vals)
    }
  }
  out$best_adjusted_evidence <- best_p
  out <- out[order(out$best_adjusted_evidence, out$result_label, out$te_subfamily), , drop = FALSE]
  out <- head(out, as.integer(max_rows))
  keep_cols <- intersect(
    c(
      "result_label",
      "te_subfamily",
      "activity_score_class",
      "contrast_aggregated_class",
      "combined_raw_p_fisher_bh_evidence",
      "contrast_aggregated_raw_p_fisher_bh_evidence",
      "n_significant_contrasts"
    ),
    colnames(out)
  )
  rownames(out) <- NULL
  out[, keep_cols, drop = FALSE]
}

comparison_activity_ranking_metric_choices <- function() {
  c(
    "sum_step_abs_delta",
    "max_step_abs_delta",
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
  )
}

comparison_activity_ranking_metric_labels <- function() {
  c(
    sum_step_abs_delta = "sum of adjacent mean activity changes",
    max_step_abs_delta = "maximum adjacent mean activity change",
    max_sum_abs_activity_score = "max condition sum of absolute activity scores",
    abs_delta_ESEN_vs_PRO = "absolute ESEN-vs-PRO mean activity gap",
    abs_delta_LSEN_vs_ESEN = "absolute LSEN-vs-ESEN mean activity gap",
    abs_delta_LSEN_vs_PRO = "absolute LSEN-vs-PRO mean activity gap",
    net_abs_delta = "net PRO-to-LSEN mean activity change",
    mean_activity_range = "condition mean activity range",
    monotonic_step_abs_delta = "monotonic adjacent mean activity change",
    PRO_ESEN = "absolute PRO-to-ESEN mean activity gap",
    ESEN_LSEN = "absolute ESEN-to-LSEN mean activity gap",
    PRO_LSEN = "absolute PRO-to-LSEN mean activity gap"
  )
}

comparison_activity_ranking_metric_label <- function(ranking_metric) {
  labels <- comparison_activity_ranking_metric_labels()
  label <- unname(labels[ranking_metric])
  if (length(label) == 0L || is.na(label) || !nzchar(label)) ranking_metric else label
}

summarise_comparison_te_activity_ranking <- function(comparison_results,
                                                     target_tes = NULL,
                                                     ranking_metric = "sum_step_abs_delta",
                                                     top_k = 30L,
                                                     aggregation = c("mean", "max"),
                                                     condition_order = NULL,
                                                     comparison_class = NULL) {
  ranking_metric <- match.arg(ranking_metric, comparison_activity_ranking_metric_choices())
  aggregation <- match.arg(aggregation)
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  if (!exists("select_top_te_subfamilies_by_activity_change", mode = "function") &&
      file.exists(file.path("scripts", "craTEs_own_function.R"))) {
    try(source("scripts/craTEs_own_function.R"), silent = TRUE)
  }
  if (!exists("select_top_te_subfamilies_by_activity_change", mode = "function")) {
    return(list(
      status = "skipped",
      reason = "select_top_te_subfamilies_by_activity_change() is unavailable. Source scripts/craTEs_own_function.R first."
    ))
  }

  shared_tes <- get_shared_te_subfamilies(comparison_results)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared TE subfamilies are available for activity-change ranking.",
      shared_te_subfamilies = shared_tes
    ))
  }

  shared_conditions <- get_shared_conditions(comparison_results)
  requested_condition_order <- condition_order
  if (is.null(condition_order)) {
    condition_order <- shared_conditions
  }
  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[nzchar(condition_order)]
  condition_order <- condition_order[condition_order %in% shared_conditions]
  if (length(condition_order) < 2L) {
    return(list(
      status = "skipped",
      reason = sprintf(
        "Activity-change ranking needs at least two shared conditions. Requested: %s. Shared: %s.",
        ifelse(is.null(requested_condition_order), "all shared conditions", paste(requested_condition_order, collapse = ", ")),
        ifelse(length(shared_conditions) == 0L, "none", paste(shared_conditions, collapse = ", "))
      ),
      shared_te_subfamilies = shared_tes,
      shared_conditions = shared_conditions
    ))
  }

  if (is.null(top_k)) {
    top_k <- length(shared_tes)
  }
  if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 1L) {
    stop("top_k must be one positive number or NULL.", call. = FALSE)
  }
  top_k <- as.integer(top_k)

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  condition_groups <- stats::setNames(as.list(condition_order), condition_order)
  result_tables <- lapply(comparison_results, function(result) {
    mat <- comparison_condition_activity_matrix(
      result,
      te_subfamilies = shared_tes,
      conditions = condition_order
    )
    mat <- mat[intersect(shared_tes, rownames(mat)), intersect(condition_order, colnames(mat)), drop = FALSE]
    if (nrow(mat) == 0L || ncol(mat) < 2L) {
      return(NULL)
    }
    ranked <- select_top_te_subfamilies_by_activity_change(
      activity_scores = mat,
      condition_groups = condition_groups[colnames(mat)],
      condition_order = colnames(mat),
      top_k = nrow(mat),
      ranking_metric = ranking_metric,
      te_subfamilies = NULL
    )
    ranked$result_id <- result$result_id
    ranked$result_label <- unname(label_map[result$result_id])
    ranked$result_label[is.na(ranked$result_label) | !nzchar(ranked$result_label)] <- result$result_id
    ranked
  })
  result_tables <- result_tables[!vapply(result_tables, is.null, logical(1))]
  if (length(result_tables) == 0L) {
    return(list(
      status = "skipped",
      reason = "No finite activity-change rankings could be computed for the selected result directories.",
      shared_te_subfamilies = shared_tes,
      shared_conditions = condition_order
    ))
  }

  result_summary <- do.call(rbind, result_tables)
  rownames(result_summary) <- NULL
  result_summary <- result_summary[is.finite(result_summary$ranking_score), , drop = FALSE]
  if (nrow(result_summary) == 0L) {
    return(list(
      status = "skipped",
      reason = "No finite ranking scores are available after activity-change ranking.",
      shared_te_subfamilies = shared_tes,
      shared_conditions = condition_order
    ))
  }

  te_summary <- do.call(rbind, lapply(split(result_summary, result_summary$te_subfamily), function(x) {
    scores <- suppressWarnings(as.numeric(x$ranking_score))
    scores <- scores[is.finite(scores)]
    global_score <- if (identical(aggregation, "max")) max(scores, na.rm = TRUE) else mean(scores, na.rm = TRUE)
    data.frame(
      te_subfamily = x$te_subfamily[1],
      global_ranking_score = global_score,
      mean_result_ranking_score = mean(scores, na.rm = TRUE),
      max_result_ranking_score = max(scores, na.rm = TRUE),
      min_result_ranking_score = min(scores, na.rm = TRUE),
      result_score_range = max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE),
      n_results = length(unique(x$result_id)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(te_summary) <- NULL
  te_summary <- te_summary[order(-te_summary$global_ranking_score, te_summary$te_subfamily), , drop = FALSE]
  te_summary$rank <- seq_len(nrow(te_summary))
  selected_te_summary <- te_summary[seq_len(min(top_k, nrow(te_summary))), , drop = FALSE]
  selected_te_summary$selected_by_global_top_k <- TRUE

  list(
    status = "ok",
    ranking_metric = ranking_metric,
    ranking_label = comparison_activity_ranking_metric_label(ranking_metric),
    ranking_aggregation = aggregation,
    top_k = top_k,
    condition_order = condition_order,
    shared_te_subfamilies = shared_tes,
    shared_conditions = condition_order,
    result_summary = result_summary,
    te_summary = te_summary,
    selected_te_summary = selected_te_summary,
    selected_te_subfamilies = selected_te_summary$te_subfamily,
    comparison_class = comparison_class
  )
}

comparison_load_absolute_activity_plot_helpers <- function() {
  required <- c("select_top_te_subfamilies_by_activity_change", "plot_delta_activity_condition_ranges")
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing) > 0L && file.exists(file.path("scripts", "craTEs_own_function.R"))) {
    try(source("scripts/craTEs_own_function.R"), silent = TRUE)
  }
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  list(ok = length(missing) == 0L, missing = missing)
}

summarise_cross_result_absolute_activity_signal <- function(comparison_results,
                                                           target_tes = NULL,
                                                           top_k = 30L,
                                                           condition_order = NULL,
                                                           comparison_class = NULL) {
  helpers <- comparison_load_absolute_activity_plot_helpers()
  if (!isTRUE(helpers$ok)) {
    return(list(
      status = "skipped",
      reason = sprintf("Missing helper function(s): %s. Source scripts/craTEs_own_function.R first.", paste(helpers$missing, collapse = ", "))
    ))
  }

  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  shared_tes <- get_shared_te_subfamilies(comparison_results)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared TE subfamilies are available for cross-result absolute activity ranking.",
      shared_te_subfamilies = shared_tes
    ))
  }

  shared_conditions <- get_shared_conditions(comparison_results)
  requested_condition_order <- condition_order
  if (is.null(condition_order)) {
    preferred <- c("PRO", "ESEN", "LSEN")
    condition_order <- preferred[preferred %in% shared_conditions]
    if (length(condition_order) == 0L) {
      condition_order <- shared_conditions
    }
  }
  condition_order <- unique(as.character(condition_order))
  condition_order <- condition_order[nzchar(condition_order)]
  condition_order <- condition_order[condition_order %in% shared_conditions]
  if (length(condition_order) == 0L) {
    return(list(
      status = "skipped",
      reason = sprintf(
        "No requested conditions are shared across result directories. Requested: %s. Shared: %s.",
        ifelse(is.null(requested_condition_order), "default/shared conditions", paste(requested_condition_order, collapse = ", ")),
        ifelse(length(shared_conditions) == 0L, "none", paste(shared_conditions, collapse = ", "))
      ),
      shared_te_subfamilies = shared_tes,
      shared_conditions = shared_conditions
    ))
  }

  if (is.null(top_k)) {
    top_k <- length(shared_tes)
  }
  if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 1L) {
    stop("top_k must be one positive number or NULL.", call. = FALSE)
  }
  top_k <- as.integer(top_k)

  result_ids <- comparison_result_ids(comparison_results)
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  result_labels <- unname(label_map[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- comparison_short_result_label(result_ids[is.na(result_labels) | !nzchar(result_labels)])
  result_labels <- make.unique(result_labels)
  names(result_labels) <- result_ids

  condition_matrices <- lapply(comparison_results, function(result) {
    mat <- comparison_condition_activity_matrix(
      result,
      te_subfamilies = shared_tes,
      conditions = condition_order
    )
    missing_conditions <- setdiff(condition_order, colnames(mat))
    if (length(missing_conditions) > 0L) {
      stop(
        sprintf("Condition mean(s) missing for %s: %s", result$result_id, paste(missing_conditions, collapse = ", ")),
        call. = FALSE
      )
    }
    mat[shared_tes, condition_order, drop = FALSE]
  })
  names(condition_matrices) <- result_ids

  shared_tes <- Reduce(intersect, lapply(condition_matrices, rownames))
  if (length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No TE subfamilies remain after cross-result condition-mean alignment.",
      shared_conditions = condition_order
    ))
  }
  condition_matrices <- lapply(condition_matrices, function(mat) mat[shared_tes, condition_order, drop = FALSE])

  aggregate_activity_scores <- do.call(cbind, lapply(result_ids, function(result_id) {
    mat <- condition_matrices[[result_id]]
    colnames(mat) <- paste(result_labels[[result_id]], condition_order, sep = "__")
    mat
  }))
  storage.mode(aggregate_activity_scores) <- "double"

  aggregate_condition_groups <- stats::setNames(
    lapply(condition_order, function(condition) paste(result_labels, condition, sep = "__")),
    condition_order
  )

  selected <- select_top_te_subfamilies_by_activity_change(
    activity_scores = aggregate_activity_scores,
    condition_groups = aggregate_condition_groups,
    condition_order = condition_order,
    top_k = top_k,
    ranking_metric = "max_sum_abs_activity_score",
    te_subfamilies = NULL
  )
  selected$rank <- seq_len(nrow(selected))
  ranked_delta_table <- data.frame(
    te_subfamily = selected$te_subfamily,
    delta_activity = selected$ranking_score,
    stringsAsFactors = FALSE
  )

  list(
    status = "ok",
    ranking_metric = "max_sum_abs_activity_score",
    ranking_label = comparison_activity_ranking_metric_label("max_sum_abs_activity_score"),
    top_k = top_k,
    condition_order = condition_order,
    shared_te_subfamilies = shared_tes,
    aggregate_activity_scores = aggregate_activity_scores,
    aggregate_condition_groups = aggregate_condition_groups,
    selected_te_summary = selected,
    ranked_delta_table = ranked_delta_table,
    result_ids = result_ids,
    result_labels = result_labels,
    comparison_class = comparison_class
  )
}

plot_cross_result_absolute_activity_signal_summary <- function(activity_signal_summary,
                                                               condition_cols = NULL) {
  comparison_need_ggplot()
  helpers <- comparison_load_absolute_activity_plot_helpers()
  if (!isTRUE(helpers$ok)) {
    skipped <- list(
      status = "skipped",
      reason = sprintf("Missing helper function(s): %s. Source scripts/craTEs_own_function.R first.", paste(helpers$missing, collapse = ", "))
    )
    comparison_warn(skipped$reason)
    return(invisible(skipped))
  }

  if (!is.list(activity_signal_summary) || is.null(activity_signal_summary$status)) {
    stop("activity_signal_summary must be a summary returned by summarise_cross_result_absolute_activity_signal().", call. = FALSE)
  }
  if (!identical(activity_signal_summary$status, "ok")) {
    comparison_warn(activity_signal_summary$reason)
    return(invisible(activity_signal_summary))
  }

  if (is.null(condition_cols)) {
    default_condition_cols <- c(PRO = "grey35", ESEN = "steelblue", LSEN = "tomato")
    condition_cols <- default_condition_cols[activity_signal_summary$condition_order]
    condition_cols <- condition_cols[!is.na(condition_cols)]
    if (length(condition_cols) == 0L) {
      condition_cols <- NULL
    }
  }

  plot <- plot_delta_activity_condition_ranges(
    activity_scores = activity_signal_summary$aggregate_activity_scores,
    delta_table = activity_signal_summary$ranked_delta_table,
    condition_groups = activity_signal_summary$aggregate_condition_groups,
    condition_order = activity_signal_summary$condition_order,
    condition_cols = condition_cols,
    condition_label = activity_signal_summary$ranking_label,
    title = sprintf("Top %s TE subfamilies by cross-result absolute activity signal", nrow(activity_signal_summary$selected_te_summary)),
    x_label = "Condition mean inferred TE activity score",
    selection_label = activity_signal_summary$ranking_label,
    show_subtitle = FALSE
  ) +
    ggplot2::labs(
      subtitle = comparison_plot_subtitle(
        sprintf(
          "Brackets show result-directory min-max ranges; dots show cross-result condition means across %s result directories.",
          length(activity_signal_summary$result_ids)
        ),
        activity_signal_summary$comparison_class
      )
    )

  attr(plot, "summary") <- activity_signal_summary
  plot
}

plot_cross_result_absolute_activity_signal <- function(comparison_results,
                                                       target_tes = NULL,
                                                       top_k = 30L,
                                                       condition_order = NULL,
                                                       condition_cols = NULL,
                                                       comparison_class = NULL) {
  activity_signal_summary <- summarise_cross_result_absolute_activity_signal(
    comparison_results = comparison_results,
    target_tes = target_tes,
    top_k = top_k,
    condition_order = condition_order,
    comparison_class = comparison_class
  )
  plot_cross_result_absolute_activity_signal_summary(
    activity_signal_summary,
    condition_cols = condition_cols
  )
}

build_te_selection_space <- function(comparison_results,
                                     target_tes = NULL,
                                     filter_mode = c("all", "significant"),
                                     significance_scope = c("any", "all"),
                                     p_score_mode = c("mean_padj", "mean_neg_log10_padj"),
                                     alpha = 0.05,
                                     comparison_class = NULL,
                                     ranking_mode = c("none", "activity_change"),
                                     ranking_metric = "sum_step_abs_delta",
                                     ranking_top_k = NULL,
                                     ranking_aggregation = c("mean", "max"),
                                     ranking_condition_order = NULL) {
  filter_mode <- match.arg(filter_mode)
  significance_scope <- match.arg(significance_scope)
  p_score_mode <- match.arg(p_score_mode)
  ranking_mode <- match.arg(ranking_mode)
  ranking_aggregation <- match.arg(ranking_aggregation)

  base_tes <- get_shared_te_subfamilies(comparison_results)
  if (!is.null(target_tes)) {
    base_tes <- target_tes[target_tes %in% base_tes]
  }
  if (length(base_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared TE subfamilies are available after applying target_tes.",
      te_subfamilies = character(),
      filter_mode = filter_mode,
      significance_scope = significance_scope
    ))
  }

  evidence <- summarise_te_padj_evidence(
    comparison_results = comparison_results,
    target_tes = base_tes,
    p_score_mode = p_score_mode,
    alpha = alpha,
    comparison_class = comparison_class
  )

  selected_tes <- base_tes
  reason <- "All shared TE subfamilies are retained."
  activity_ranking <- NULL
  if (identical(filter_mode, "significant")) {
    if (!identical(evidence$status, "ok")) {
      return(list(
        status = "skipped",
        reason = sprintf("Significant TE filtering was requested, but contrast evidence is unavailable: %s", evidence$reason),
        te_subfamilies = character(),
        filter_mode = filter_mode,
        significance_scope = significance_scope,
        evidence = evidence
      ))
    }
    keep_col <- if (identical(significance_scope, "all")) "all_results_significant" else "any_result_significant"
    selected_tes <- evidence$te_summary$te_subfamily[evidence$te_summary[[keep_col]]]
    reason <- sprintf(
      "TE subfamilies are retained when they are significant in %s selected result director%s.",
      ifelse(identical(significance_scope, "all"), "all", "at least one"),
      ifelse(identical(significance_scope, "all"), "ies", "y")
    )
  }

  selected_tes <- selected_tes[selected_tes %in% base_tes]
  if (identical(ranking_mode, "activity_change")) {
    activity_ranking <- summarise_comparison_te_activity_ranking(
      comparison_results = comparison_results,
      target_tes = selected_tes,
      ranking_metric = ranking_metric,
      top_k = ranking_top_k,
      aggregation = ranking_aggregation,
      condition_order = ranking_condition_order,
      comparison_class = comparison_class
    )
    if (!identical(activity_ranking$status, "ok")) {
      return(list(
        status = "skipped",
        reason = sprintf("Activity-change TE ranking was requested, but ranking is unavailable: %s", activity_ranking$reason),
        te_subfamilies = character(),
        filter_mode = filter_mode,
        significance_scope = significance_scope,
        p_score_mode = p_score_mode,
        alpha = alpha,
        ranking_mode = ranking_mode,
        activity_ranking = activity_ranking,
        evidence = evidence
      ))
    }
    selected_tes <- activity_ranking$selected_te_subfamilies
    reason <- paste(
      reason,
      sprintf(
        "Then globally ranked by %s (%s across result directories) and top %s retained.",
        activity_ranking$ranking_label,
        activity_ranking$ranking_aggregation,
        length(selected_tes)
      )
    )
  }
  selection_summary <- data.frame(
    filter_mode = filter_mode,
    significance_scope = ifelse(identical(filter_mode, "significant"), significance_scope, "not used"),
    p_score_mode = p_score_mode,
    alpha = alpha,
    ranking_mode = ranking_mode,
    ranking_metric = ifelse(identical(ranking_mode, "activity_change"), ranking_metric, "not used"),
    ranking_aggregation = ifelse(identical(ranking_mode, "activity_change"), ranking_aggregation, "not used"),
    ranking_top_k = ifelse(identical(ranking_mode, "activity_change"), length(selected_tes), NA_integer_),
    n_shared_activity_tes = length(base_tes),
    n_selected_tes = length(selected_tes),
    note = reason,
    stringsAsFactors = FALSE
  )

  list(
    status = "ok",
    te_subfamilies = selected_tes,
    filter_mode = filter_mode,
    significance_scope = significance_scope,
    p_score_mode = p_score_mode,
    alpha = alpha,
    ranking_mode = ranking_mode,
    selection_summary = selection_summary,
    activity_ranking = activity_ranking,
    evidence = evidence
  )
}

summarise_te_contrasts <- function(comparison_results, target_tes = NULL, alpha = 0.05) {
  has_contrast <- vapply(comparison_results, function(result) !is.null(result$tables$te_contrast_summary), logical(1))
  if (!all(has_contrast)) {
    missing_ids <- comparison_result_ids(comparison_results)[!has_contrast]
    return(list(
      status = "skipped",
      reason = sprintf("Missing te_contrast_summary.tsv in: %s", paste(missing_ids, collapse = ", "))
    ))
  }

  contrast_sets <- lapply(comparison_results, function(result) unique(result$tables$te_contrast_summary$te_subfamily))
  shared_tes <- Reduce(intersect, contrast_sets)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(status = "skipped", reason = "No shared TE subfamilies are available in te_contrast_summary.tsv."))
  }

  common_cols <- Reduce(intersect, lapply(comparison_results, function(result) colnames(result$tables$te_contrast_summary)))
  effect_cols <- grep("^(delta_|estimate)", common_cols, value = TRUE)
  p_cols <- grep("^(padj_|p_adj|p_value)", common_cols, value = TRUE)
  signif_cols <- grep("^(signif_|significant)", common_cols, value = TRUE)
  class_cols <- intersect("trend_class", common_cols)
  stat_cols <- comparison_statistic_columns(data.frame(setNames(rep(list(numeric()), length(common_cols)), common_cols), check.names = FALSE))

  combined <- do.call(rbind, lapply(comparison_results, function(result) {
    x <- result$tables$te_contrast_summary
    x <- x[x$te_subfamily %in% shared_tes, common_cols, drop = FALSE]
    x$result_id <- result$result_id
    x
  }))
  rownames(combined) <- NULL

  significant_sets <- lapply(comparison_results, function(result) {
    x <- result$tables$te_contrast_summary
    if (length(signif_cols) > 0L) {
      is_sig <- rowSums(x[signif_cols] == TRUE, na.rm = TRUE) > 0L
    } else if (length(p_cols) > 0L) {
      is_sig <- apply(x[p_cols], 1, function(row) any(as.numeric(row) <= alpha, na.rm = TRUE))
    } else {
      is_sig <- rep(FALSE, nrow(x))
    }
    x$te_subfamily[is_sig]
  })
  names(significant_sets) <- comparison_result_ids(comparison_results)

  significant_summary <- data.frame(
    result_id = names(significant_sets),
    n_significant_te_subfamilies = vapply(significant_sets, length, integer(1)),
    stringsAsFactors = FALSE
  )
  shared_significant <- if (length(significant_sets) > 0L) Reduce(intersect, significant_sets) else character()

  direction_consistency <- data.frame()
  if (length(effect_cols) > 0L) {
    direction_consistency <- do.call(rbind, lapply(effect_cols, function(effect_col) {
      by_te <- split(combined[, c("te_subfamily", "result_id", effect_col)], combined$te_subfamily)
      consistent <- vapply(by_te, function(x) {
        signs <- sign(as.numeric(x[[effect_col]]))
        signs <- signs[!is.na(signs) & signs != 0]
        length(unique(signs)) <= 1L
      }, logical(1))
      data.frame(
        effect_column = effect_col,
        n_te_subfamilies = length(consistent),
        n_direction_consistent = sum(consistent),
        fraction_direction_consistent = if (length(consistent) == 0L) NA_real_ else mean(consistent),
        stringsAsFactors = FALSE
      )
    }))
  }

  list(
    status = "ok",
    shared_te_subfamilies = shared_tes,
    common_columns = common_cols,
    effect_columns = effect_cols,
    statistic_columns = stat_cols,
    p_columns = p_cols,
    significance_columns = signif_cols,
    class_columns = class_cols,
    combined = combined,
    significant_summary = significant_summary,
    shared_significant_te_subfamilies = shared_significant,
    direction_consistency = direction_consistency
  )
}

summarise_target_gene_matrices <- function(comparison_results, target_tes = NULL, target_genes = NULL) {
  matrix_names <- get_common_target_matrix_names(comparison_results)
  if (length(matrix_names) == 0L) {
    return(list(
      status = "skipped",
      reason = "No target_gene_matrix_*.tsv file is shared by all selected result directories."
    ))
  }

  all_heatmap <- list()
  all_support <- list()

  for (matrix_name in matrix_names) {
    matrices <- lapply(comparison_results, function(result) result$target_gene_matrices[[matrix_name]])
    shared_genes <- Reduce(intersect, lapply(matrices, function(x) unique(x$row_id)))
    shared_tes <- Reduce(intersect, lapply(matrices, function(x) setdiff(colnames(x), "row_id")))

    if (!is.null(target_genes)) {
      shared_genes <- intersect(shared_genes, target_genes)
    }
    if (!is.null(target_tes)) {
      shared_tes <- target_tes[target_tes %in% shared_tes]
    }

    if (length(shared_genes) == 0L || length(shared_tes) == 0L) {
      next
    }

    matrix_by_result <- lapply(matrices, function(x) {
      rownames(x) <- x$row_id
      x <- x[shared_genes, c("row_id", shared_tes), drop = FALSE]
      mat <- as.matrix(x[, shared_tes, drop = FALSE])
      storage.mode(mat) <- "double"
      mat
    })
    names(matrix_by_result) <- comparison_result_ids(comparison_results)

    all_heatmap[[matrix_name]] <- do.call(rbind, lapply(names(matrix_by_result), function(result_id) {
      x <- comparison_melt_matrix(matrix_by_result[[result_id]], "gene_id", "te_subfamily", "target_score")
      x$result_id <- result_id
      x$target_matrix <- matrix_name
      x
    }))

    all_support[[matrix_name]] <- do.call(rbind, lapply(shared_tes, function(te) {
      supported_by_result <- lapply(matrix_by_result, function(mat) rownames(mat)[abs(mat[, te]) > .Machine$double.eps])
      data.frame(
        target_matrix = matrix_name,
        te_subfamily = te,
        n_shared_genes = length(shared_genes),
        n_supported_in_all_results = length(Reduce(intersect, supported_by_result)),
        stringsAsFactors = FALSE
      )
    }))
  }

  if (length(all_heatmap) == 0L) {
    return(list(status = "skipped", reason = "No shared target gene matrix space remains after filtering."))
  }

  list(
    status = "ok",
    common_target_matrices = matrix_names,
    heatmap_data = do.call(rbind, all_heatmap),
    support_summary = do.call(rbind, all_support)
  )
}

comparison_fit_metric_label <- function(metric) {
  switch(
    metric,
    mean_r2 = "Mean R2",
    mean_r2_ratio = "Mean R2 ratio",
    mean_sparsity_percent = "Mean sparsity (%)",
    metric
  )
}

comparison_fit_condition_mode <- function(fit, comparison_class) {
  condition_rows <- fit[
    tolower(as.character(fit$scope)) == "condition" &
      !tolower(as.character(fit$condition)) %in% c("all", "overall", ""),
    ,
    drop = FALSE
  ]

  if (nrow(condition_rows) == 0L) {
    return(list(mode = "overall", conditions = character(), note = "No condition-wise rows are available."))
  }

  if (isTRUE(comparison_class$same_E)) {
    conditions <- unique(as.character(condition_rows$condition))
    return(list(
      mode = "condition",
      conditions = conditions,
      note = "Same E: condition-wise fit rows are shown inside each overall bar."
    ))
  }

  by_result <- split(as.character(condition_rows$condition), condition_rows$result_id)
  if (length(by_result) == 0L) {
    return(list(mode = "overall", conditions = character(), note = "No condition-wise rows are available."))
  }
  shared_conditions <- Reduce(intersect, lapply(by_result, unique))
  if (length(shared_conditions) == 0L) {
    return(list(
      mode = "overall",
      conditions = character(),
      note = sprintf(
        "%s: condition names are not shared across all result directories, so only overall fit is shown.",
        comparison_class$class_label
      )
    ))
  }

  list(
    mode = "condition",
    conditions = shared_conditions,
    note = sprintf(
      "%s: condition-wise fit is restricted to shared conditions.",
      comparison_class$class_label
    )
  )
}

comparison_fit_palette <- function(groups) {
  base <- c(
    overall = "#4c78a8",
    "#b2182b", "#1b7837", "#762a83", "#b35806",
    "#0072b2", "#008c8c", "#5e3c99", "#8c510a",
    "#01665e", "#c51b7d", "#666666"
  )
  if (length(groups) <= length(base)) {
    out <- base[seq_along(groups)]
  } else {
    out <- grDevices::hcl.colors(length(groups), palette = "Dark 3")
    out[1] <- base[["overall"]]
  }
  names(out) <- groups
  out
}

comparison_prepare_fit_metric_plot <- function(comparison_results, fit, comparison_class, metric) {
  labels <- comparison_dataset_labels(comparison_results, comparison_class)
  result_ids <- comparison_result_ids(comparison_results)
  result_labels <- unname(labels[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- comparison_short_result_label(result_ids[is.na(result_labels) | !nzchar(result_labels)])

  y_lookup <- data.frame(
    result_id = result_ids,
    result_label = result_labels,
    y = rev(seq_along(result_ids)),
    stringsAsFactors = FALSE
  )

  overall <- fit[
    tolower(as.character(fit$scope)) == "overall" |
      tolower(as.character(fit$condition)) %in% c("all", "overall"),
    ,
    drop = FALSE
  ]
  if (nrow(overall) == 0L) {
    overall <- fit[!duplicated(fit$result_id), , drop = FALSE]
  }
  overall <- overall[match(result_ids, overall$result_id), , drop = FALSE]
  overall <- merge(overall, y_lookup, by = "result_id", all.x = TRUE, sort = FALSE)
  overall$value <- as.numeric(overall[[metric]])
  overall$fit_group <- "overall"

  condition_mode <- comparison_fit_condition_mode(fit, comparison_class)
  condition <- data.frame()
  if (identical(condition_mode$mode, "condition")) {
    condition <- fit[
      tolower(as.character(fit$scope)) == "condition" &
        as.character(fit$condition) %in% condition_mode$conditions,
      ,
      drop = FALSE
    ]
    condition <- merge(condition, y_lookup, by = "result_id", all.x = TRUE, sort = FALSE)
    condition$value <- as.numeric(condition[[metric]])
    condition$fit_group <- as.character(condition$condition)

    condition <- do.call(rbind, lapply(split(condition, condition$result_id), function(x) {
      x <- x[match(condition_mode$conditions, x$condition), , drop = FALSE]
      x <- x[!is.na(x$result_id), , drop = FALSE]
      n <- nrow(x)
      if (n == 0L) {
        return(x)
      }
      offsets <- if (n == 1L) 0 else seq(-0.24, 0.24, length.out = n)
      x$y_condition <- x$y + offsets
      x$condition_height <- min(0.14, 0.48 / n)
      x
    }))
    rownames(condition) <- NULL
  }

  overall_height <- if (nrow(condition) == 0L) 0.42 else 0.72
  overall$ymin <- overall$y - overall_height / 2
  overall$ymax <- overall$y + overall_height / 2
  overall$xmin <- 0
  overall$xmax <- overall$value

  if (nrow(condition) > 0L) {
    condition$ymin <- condition$y_condition - condition$condition_height / 2
    condition$ymax <- condition$y_condition + condition$condition_height / 2
    condition$xmin <- 0
    condition$xmax <- condition$value
  }

  values <- c(overall$value, condition$value)
  values <- values[!is.na(values)]
  x_max <- if (length(values) == 0L) 1 else max(values, 0, na.rm = TRUE)
  x_min <- if (length(values) == 0L) 0 else min(values, 0, na.rm = TRUE)
  pad <- max((x_max - x_min) * 0.08, 0.01)

  condition_levels <- if (nrow(condition) == 0L) character() else unique(as.character(condition$condition))
  groups <- c("overall", condition_levels)

  list(
    overall = overall,
    condition = condition,
    y_lookup = y_lookup,
    condition_mode = condition_mode,
    groups = groups,
    colors = comparison_fit_palette(groups),
    x_limits = c(min(0, x_min) - pad * 2.8, x_max + pad * 3.2),
    label_pad = pad
  )
}

plot_activity_fit_metric <- function(comparison_results, metric, comparison_class = NULL) {
  comparison_need_ggplot()
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  fit <- summarise_activity_fit(comparison_results)
  if (!metric %in% colnames(fit)) {
    comparison_warn(sprintf("Model-fit metric `%s` is not available for plotting.", metric))
    return(invisible(NULL))
  }

  plot_data <- comparison_prepare_fit_metric_plot(comparison_results, fit, comparison_class, metric)
  overall <- plot_data$overall
  condition <- plot_data$condition
  show_conditions <- nrow(condition) > 0L

  overall_label <- data.frame(
    x = -plot_data$label_pad * 0.35,
    y = overall$y,
    label = sprintf("%.4f", overall$value),
    fit_group = "overall",
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = overall,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fit_group),
      alpha = if (show_conditions) 0.52 else 0.72,
      color = NA
    )

  if (show_conditions) {
    condition_label <- data.frame(
      x = condition$value + plot_data$label_pad * 0.2,
      y = condition$y_condition,
      label = sprintf("%.4f", condition$value),
      fit_group = condition$fit_group,
      stringsAsFactors = FALSE
    )

    p <- p +
      ggplot2::geom_rect(
        data = condition,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fit_group),
        alpha = 0.48,
        color = NA
      ) +
      ggplot2::geom_text(
        data = condition_label,
        ggplot2::aes(x = x, y = y, label = label, color = fit_group),
        hjust = 0,
        size = 3.2,
        fontface = "bold"
      )
  }

  p +
    ggplot2::geom_text(
      data = overall_label,
      ggplot2::aes(x = x, y = y, label = label, color = fit_group),
      hjust = 1,
      size = 3.4,
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      breaks = plot_data$y_lookup$y,
      labels = plot_data$y_lookup$result_label,
      expand = ggplot2::expansion(mult = if (show_conditions) 0.08 else 0.03)
    ) +
    ggplot2::scale_x_continuous(limits = plot_data$x_limits) +
    ggplot2::scale_fill_manual(values = plot_data$colors, breaks = plot_data$groups, drop = FALSE) +
    ggplot2::scale_color_manual(values = plot_data$colors, breaks = plot_data$groups, drop = FALSE, guide = "none") +
    ggplot2::labs(
      title = sprintf("%s model fit", comparison_fit_metric_label(metric)),
      subtitle = comparison_plot_subtitle(
        plot_data$condition_mode$note,
        comparison_class
      ),
      x = comparison_fit_metric_label(metric),
      y = NULL,
      fill = NULL
    ) +
    comparison_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = if (show_conditions) "bottom" else "none"
    )
}

plot_activity_fit_comparison <- function(comparison_results, comparison_class = NULL) {
  fit <- summarise_activity_fit(comparison_results)
  metric_cols <- intersect(c("mean_r2", "mean_r2_ratio"), colnames(fit))
  if (length(metric_cols) == 0L) {
    comparison_warn("No model-fit metric columns are available for plotting.")
    return(invisible(list()))
  }

  plots <- lapply(metric_cols, function(metric) {
    plot_activity_fit_metric(comparison_results, metric = metric, comparison_class = comparison_class)
  })
  names(plots) <- metric_cols
  plots
}

comparison_activity_result_label_map <- function(comparison_results, comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  labels <- comparison_dataset_labels(comparison_results, comparison_class)
  ids <- comparison_result_ids(comparison_results)
  out <- unname(labels[ids])
  out[is.na(out) | !nzchar(out)] <- comparison_short_result_label(ids[is.na(out) | !nzchar(out)])
  names(out) <- ids
  out
}

plot_activity_score_heatmap <- function(comparison_results, target_tes = NULL, max_tes = 30L, comparison_class = NULL) {
  comparison_need_ggplot()
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  summary <- summarise_activity_scores(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }

  plot_data <- summary$heatmap_data
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  plot_data$result_label <- unname(label_map[plot_data$result_id])
  plot_data$result_label[is.na(plot_data$result_label) | !nzchar(plot_data$result_label)] <- plot_data$result_id[is.na(plot_data$result_label) | !nzchar(plot_data$result_label)]

  te_order <- unique(plot_data$te_subfamily)
  heatmap_order <- comparison_order_activity_heatmap_tes(
    plot_data = plot_data,
    comparison_results = comparison_results,
    target_tes = target_tes,
    max_tes = max_tes
  )
  if (is.null(target_tes) && length(unique(plot_data$te_subfamily)) > max_tes) {
    plot_data <- plot_data[plot_data$te_subfamily %in% heatmap_order$te_order, , drop = FALSE]
  }
  if (!is.null(heatmap_order$class_data)) {
    plot_data <- merge(plot_data, heatmap_order$class_data, by = "te_subfamily", all.x = TRUE, sort = FALSE)
    plot_data$majority_class_display[is.na(plot_data$majority_class_display) | !nzchar(plot_data$majority_class_display)] <- "unclassified"
  }
  te_order <- heatmap_order$te_order[heatmap_order$te_order %in% unique(plot_data$te_subfamily)]
  if (length(te_order) > 0L) {
    plot_data$te_subfamily <- factor(plot_data$te_subfamily, levels = rev(te_order))
  }

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_id, y = te_subfamily, fill = activity_score)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    ggplot2::labs(
      title = "Activity-score heatmap",
      subtitle = comparison_plot_subtitle(
        sprintf(
          "Shared TE subfamilies and %s only; top k TE = %s%s.",
          summary$comparison_space,
          ifelse(is.null(target_tes), max_tes, length(unique(plot_data$te_subfamily))),
          ifelse(!is.null(target_tes) && isTRUE(heatmap_order$has_classification), "; target TEs grouped by majority trend class", "")
        ),
        comparison_class
      ),
      x = summary$axis_label,
      y = "TE subfamily",
      fill = "Activity"
    ) +
    comparison_theme()

  if (!is.null(target_tes) && isTRUE(heatmap_order$has_classification)) {
    plot + ggplot2::facet_grid(result_label + majority_class_display ~ ., scales = "free_y", space = "free_y")
  } else {
    plot + ggplot2::facet_wrap(~ result_label, ncol = 1)
  }
}

plot_activity_score_similarity <- function(comparison_results, target_tes = NULL, comparison_class = NULL) {
  comparison_need_ggplot()
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  summary <- summarise_activity_scores(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  cor <- summary$activity_correlation
  rownames(cor) <- unname(label_map[rownames(cor)])
  colnames(cor) <- unname(label_map[colnames(cor)])
  cor_df <- comparison_melt_matrix(cor, "result_a", "result_b", "correlation")
  label_levels <- colnames(cor)
  cor_df$result_a <- factor(cor_df$result_a, levels = label_levels)
  cor_df$result_b <- factor(cor_df$result_b, levels = rev(label_levels))

  ggplot2::ggplot(cor_df, ggplot2::aes(x = result_a, y = result_b, fill = correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", correlation)), size = 3.4) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Activity-score similarity",
      subtitle = comparison_plot_subtitle(
        "Correlation on the shared activity-score space.",
        comparison_class
      ),
      x = NULL,
      y = NULL,
      fill = "Pearson r"
    ) +
    comparison_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )
}

comparison_one_line_label <- function(x) {
  trimws(gsub("\\s+", " ", as.character(x)))
}

summarise_activity_sample_similarity <- function(comparison_results,
                                                target_tes = NULL,
                                                method = "pearson",
                                                comparison_class = NULL) {
  method <- match.arg(method, choices = c("pearson", "spearman", "kendall"))
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  summary <- comparison_activity_matrices_shared_tes(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    return(summary)
  }
  if (length(summary$shared_samples) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared sample columns are available across all selected result directories.",
      shared_te_subfamilies = summary$shared_te_subfamilies,
      shared_samples = character()
    ))
  }

  result_ids <- names(summary$matrices)
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  result_labels <- unname(label_map[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- result_ids[is.na(result_labels) | !nzchar(result_labels)]
  names(result_labels) <- result_ids

  pair_idx <- utils::combn(result_ids, 2, simplify = FALSE)
  sample_rows <- do.call(rbind, lapply(pair_idx, function(pair) {
    mat_a <- summary$matrices[[pair[1]]][, summary$shared_samples, drop = FALSE]
    mat_b <- summary$matrices[[pair[2]]][, summary$shared_samples, drop = FALSE]
    result_pair_label <- comparison_one_line_label(paste(
      result_labels[[pair[1]]],
      result_labels[[pair[2]]],
      sep = " vs "
    ))

    do.call(rbind, lapply(summary$shared_samples, function(sample_id) {
      correlation <- suppressWarnings(stats::cor(
        mat_a[, sample_id],
        mat_b[, sample_id],
        use = "pairwise.complete.obs",
        method = method
      ))
      data.frame(
        result_a = pair[1],
        result_b = pair[2],
        result_a_label = result_labels[[pair[1]]],
        result_b_label = result_labels[[pair[2]]],
        result_pair = result_pair_label,
        result_pair_label = result_pair_label,
        sample_id = sample_id,
        correlation = correlation,
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(sample_rows) <- NULL

  list(
    status = "ok",
    method = method,
    shared_te_subfamilies = summary$shared_te_subfamilies,
    shared_samples = summary$shared_samples,
    sample_similarity = sample_rows
  )
}

plot_activity_sample_similarity <- function(activity_sample_similarity) {
  comparison_need_ggplot()
  if (!identical(activity_sample_similarity$status, "ok")) {
    comparison_warn(activity_sample_similarity$reason)
    return(invisible(activity_sample_similarity))
  }

  plot_data <- activity_sample_similarity$sample_similarity
  if (!"result_pair_label" %in% colnames(plot_data)) {
    plot_data$result_pair_label <- comparison_one_line_label(plot_data$result_pair)
  }
  pair_levels <- unique(plot_data$result_pair_label)
  plot_data$result_pair_label <- factor(plot_data$result_pair_label, levels = rev(pair_levels))
  plot_data$sample_id <- factor(plot_data$sample_id, levels = activity_sample_similarity$shared_samples)
  method_label <- switch(
    activity_sample_similarity$method,
    pearson = "Pearson r",
    spearman = "Spearman rho",
    kendall = "Kendall tau",
    "Correlation"
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_id, y = result_pair_label, fill = correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(correlation), "NA", sprintf("%.2f", correlation))), size = 3.3) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), na.value = "#f2f4f7") +
    ggplot2::scale_y_discrete(labels = comparison_one_line_label) +
    ggplot2::labs(
      title = "Sample-wise activity-score similarity",
      subtitle = "Correlation per shared sample across shared TE subfamilies.",
      x = "Sample",
      y = NULL,
      fill = method_label
    ) +
    comparison_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )
}

comparison_activity_sample_conditions <- function(result) {
  sample_ids <- comparison_activity_sample_ids(result)
  fit <- result$tables$activity_fit_summary
  if (is.null(fit) || !"condition" %in% colnames(fit)) {
    return(data.frame(sample_id = sample_ids, condition = NA_character_, stringsAsFactors = FALSE))
  }

  condition_rows <- fit[
    tolower(as.character(fit$scope)) == "condition" &
      !tolower(as.character(fit$condition)) %in% c("all", "overall", ""),
    ,
    drop = FALSE
  ]
  conditions <- unique(as.character(condition_rows$condition))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  conditions <- conditions[order(nchar(conditions), decreasing = TRUE)]

  sample_condition <- vapply(sample_ids, function(sample_id) {
    matched <- conditions[vapply(conditions, function(condition) {
      grepl(condition, sample_id, fixed = TRUE)
    }, logical(1))]
    if (length(matched) == 0L) NA_character_ else matched[[1]]
  }, character(1))

  data.frame(sample_id = sample_ids, condition = sample_condition, stringsAsFactors = FALSE)
}

summarise_activity_cross_sample_similarity <- function(comparison_results,
                                                      target_tes = NULL,
                                                      method = "pearson",
                                                      comparison_class = NULL) {
  method <- match.arg(method, choices = c("pearson", "spearman", "kendall"))
  if (length(comparison_results) != 2L) {
    return(list(
      status = "skipped",
      reason = sprintf("Cross-sample activity similarity requires exactly two result directories; got %d.", length(comparison_results))
    ))
  }
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  summary <- comparison_activity_matrices_shared_tes(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    return(summary)
  }

  result_ids <- names(summary$matrices)
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  result_labels <- unname(label_map[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- comparison_short_result_label(result_ids[is.na(result_labels) | !nzchar(result_labels)])
  names(result_labels) <- result_ids

  sample_meta <- do.call(rbind, lapply(seq_along(comparison_results), function(i) {
    result <- comparison_results[[i]]
    result_id <- result_ids[[i]]
    x <- comparison_activity_sample_conditions(result)
    x$result_id <- result_id
    x$result_label <- result_labels[[result_id]]
    x$profile_id <- paste(result_id, x$sample_id, sep = "::")
    x$profile_label <- comparison_one_line_label(paste(x$result_label, x$sample_id, sep = " | "))
    x[, c("profile_id", "result_id", "result_label", "sample_id", "condition", "profile_label"), drop = FALSE]
  }))
  rownames(sample_meta) <- NULL

  sample_a <- sample_meta[sample_meta$result_id == result_ids[[1]], , drop = FALSE]
  sample_b <- sample_meta[sample_meta$result_id == result_ids[[2]], , drop = FALSE]
  if (nrow(sample_a) == 0L || nrow(sample_b) == 0L) {
    return(list(
      status = "skipped",
      reason = "Both selected result directories must contain activity-score sample columns.",
      shared_te_subfamilies = summary$shared_te_subfamilies
    ))
  }

  shared_samples <- sample_a$sample_id[sample_a$sample_id %in% sample_b$sample_id]
  if (length(shared_samples) == 0L) {
    return(list(
      status = "skipped",
      reason = "Cross-sample half-matrix requires shared sample IDs across the two selected result directories.",
      shared_te_subfamilies = summary$shared_te_subfamilies,
      shared_samples = character()
    ))
  }
  if (!setequal(sample_a$sample_id, sample_b$sample_id)) {
    return(list(
      status = "skipped",
      reason = "Cross-sample half-matrix requires the two selected result directories to have the same sample IDs.",
      shared_te_subfamilies = summary$shared_te_subfamilies,
      shared_samples = shared_samples
    ))
  }

  sample_a <- sample_a[match(shared_samples, sample_a$sample_id), , drop = FALSE]
  sample_b <- sample_b[match(shared_samples, sample_b$sample_id), , drop = FALSE]
  mat_a <- summary$matrices[[result_ids[[1]]]][, sample_a$sample_id, drop = FALSE]
  mat_b <- summary$matrices[[result_ids[[2]]]][, sample_b$sample_id, drop = FALSE]
  pair_grid <- expand.grid(
    sample_a_index = seq_len(nrow(sample_a)),
    sample_b_index = seq_len(nrow(sample_b)),
    KEEP.OUT.ATTRS = FALSE
  )

  pairwise <- do.call(rbind, lapply(seq_len(nrow(pair_grid)), function(i) {
    row_a <- sample_a[pair_grid$sample_a_index[[i]], , drop = FALSE]
    row_b <- sample_b[pair_grid$sample_b_index[[i]], , drop = FALSE]
    correlation <- suppressWarnings(stats::cor(
      mat_a[, row_a$sample_id],
      mat_b[, row_b$sample_id],
      use = "pairwise.complete.obs",
      method = method
    ))
    same_sample <- identical(row_a$sample_id[[1]], row_b$sample_id[[1]])
    same_condition <- !is.na(row_a$condition[[1]]) &&
      !is.na(row_b$condition[[1]]) &&
      nzchar(row_a$condition[[1]]) &&
      nzchar(row_b$condition[[1]]) &&
      identical(row_a$condition[[1]], row_b$condition[[1]])

    data.frame(
      profile_a = row_a$profile_id,
      profile_b = row_b$profile_id,
      result_a = row_a$result_id,
      result_b = row_b$result_id,
      result_a_label = row_a$result_label,
      result_b_label = row_b$result_label,
      sample_a = row_a$sample_id,
      sample_b = row_b$sample_id,
      condition_a = row_a$condition,
      condition_b = row_b$condition,
      same_condition = same_condition,
      same_sample = same_sample,
      correlation = correlation,
      stringsAsFactors = FALSE
    )
  }))
  rownames(pairwise) <- NULL

  list(
    status = "ok",
    method = method,
    result_ids = result_ids,
    result_labels = result_labels,
    profiles = sample_meta,
    pairwise_similarity = pairwise,
    shared_te_subfamilies = summary$shared_te_subfamilies,
    shared_samples = shared_samples,
    comparison_class = comparison_class
  )
}

plot_activity_cross_sample_similarity <- function(activity_cross_sample_similarity,
                                                  triangle = c("upper", "full"),
                                                  show_values = TRUE,
                                                  text_size = 2.7) {
  comparison_need_ggplot()
  triangle <- match.arg(triangle)
  if (!identical(activity_cross_sample_similarity$status, "ok")) {
    comparison_warn(activity_cross_sample_similarity$reason)
    return(invisible(activity_cross_sample_similarity))
  }

  sample_levels <- activity_cross_sample_similarity$shared_samples
  if (is.null(sample_levels) || length(sample_levels) == 0L) {
    comparison_warn("Cross-sample activity similarity has no shared sample axis for half-matrix plotting.")
    return(invisible(activity_cross_sample_similarity))
  }
  pairwise <- activity_cross_sample_similarity$pairwise_similarity
  if (identical(triangle, "full")) {
    mirrored <- pairwise
    swap_cols <- list(
      profile_a = "profile_b",
      profile_b = "profile_a",
      result_a = "result_b",
      result_b = "result_a",
      result_a_label = "result_b_label",
      result_b_label = "result_a_label",
      sample_a = "sample_b",
      sample_b = "sample_a",
      condition_a = "condition_b",
      condition_b = "condition_a"
    )
    mirrored_original <- mirrored
    for (to_col in names(swap_cols)) {
      mirrored[[to_col]] <- mirrored_original[[swap_cols[[to_col]]]]
    }
    pairwise <- rbind(pairwise, mirrored)
  } else {
    pairwise$x_index <- match(pairwise$sample_b, sample_levels)
    pairwise$y_index <- match(pairwise$sample_a, sample_levels)
    pairwise <- pairwise[
      !is.na(pairwise$x_index) &
        !is.na(pairwise$y_index) &
        pairwise$x_index >= pairwise$y_index,
      ,
      drop = FALSE
    ]
  }

  pairwise$x_profile <- if (identical(triangle, "upper")) pairwise$sample_b else pairwise$sample_a
  pairwise$y_profile <- if (identical(triangle, "upper")) pairwise$sample_a else pairwise$sample_b
  frame_colors <- c(
    "same condition" = "#00A676",
    "exact sample" = "#FF8C00"
  )
  pairwise$plot_value <- pairwise$correlation
  pairwise$tile_label <- ifelse(is.na(pairwise$plot_value), "NA", sprintf("%.2f", pairwise$plot_value))
  pairwise$x_profile <- factor(pairwise$x_profile, levels = sample_levels)
  pairwise$y_profile <- factor(pairwise$y_profile, levels = rev(sample_levels))
  pairwise$x_num <- as.numeric(pairwise$x_profile)
  pairwise$y_num <- as.numeric(pairwise$y_profile)

  empty_border_segments <- data.frame(
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric(),
    orientation = character(),
    fixed = numeric(),
    start = numeric(),
    end = numeric(),
    normal_x = numeric(),
    normal_y = numeric(),
    frame_label = factor(character(), levels = names(frame_colors)),
    stringsAsFactors = FALSE
  )

  merge_collinear_segments <- function(segments) {
    if (nrow(segments) == 0L) {
      return(segments)
    }
    has_shifted <- "shifted" %in% colnames(segments)
    if (!isTRUE(has_shifted)) {
      segments$shifted <- FALSE
    }

    segments$fixed <- ifelse(segments$orientation == "h", segments$y, segments$x)
    segments$start <- ifelse(segments$orientation == "h", pmin(segments$x, segments$xend), pmin(segments$y, segments$yend))
    segments$end <- ifelse(segments$orientation == "h", pmax(segments$x, segments$xend), pmax(segments$y, segments$yend))
    split_key <- paste(
      segments$frame_label,
      segments$orientation,
      sprintf("%.6f", segments$fixed),
      sprintf("%.6f", segments$normal_x),
      sprintf("%.6f", segments$normal_y),
      sep = "::"
    )

    merged <- do.call(rbind, lapply(split(segments, split_key), function(x) {
      x <- x[order(x$start, x$end), , drop = FALSE]
      out <- list()
      current_start <- x$start[[1]]
      current_end <- x$end[[1]]

      if (nrow(x) > 1L) {
        for (i in 2:nrow(x)) {
          if (x$start[[i]] <= current_end + 1e-8) {
            current_end <- max(current_end, x$end[[i]])
          } else {
            out[[length(out) + 1L]] <- data.frame(
              frame_label = x$frame_label[[1]],
              orientation = x$orientation[[1]],
              fixed = x$fixed[[1]],
              start = current_start,
              end = current_end,
              normal_x = x$normal_x[[1]],
              normal_y = x$normal_y[[1]],
              shifted = x$shifted[[1]],
              stringsAsFactors = FALSE
            )
            current_start <- x$start[[i]]
            current_end <- x$end[[i]]
          }
        }
      }

      out[[length(out) + 1L]] <- data.frame(
        frame_label = x$frame_label[[1]],
        orientation = x$orientation[[1]],
        fixed = x$fixed[[1]],
        start = current_start,
        end = current_end,
        normal_x = x$normal_x[[1]],
        normal_y = x$normal_y[[1]],
        shifted = x$shifted[[1]],
        stringsAsFactors = FALSE
      )
      do.call(rbind, out)
    }))

    out <- data.frame(
      x = ifelse(merged$orientation == "h", merged$start, merged$fixed),
      y = ifelse(merged$orientation == "h", merged$fixed, merged$start),
      xend = ifelse(merged$orientation == "h", merged$end, merged$fixed),
      yend = ifelse(merged$orientation == "h", merged$fixed, merged$end),
      orientation = merged$orientation,
      fixed = merged$fixed,
      start = merged$start,
      end = merged$end,
      normal_x = merged$normal_x,
      normal_y = merged$normal_y,
      shifted = merged$shifted,
      frame_label = factor(merged$frame_label, levels = names(frame_colors)),
      stringsAsFactors = FALSE
    )
    if (!isTRUE(has_shifted)) {
      out$shifted <- NULL
    }
    out
  }

  make_border_segments <- function(cells, frame_label) {
    if (nrow(cells) == 0L) {
      return(empty_border_segments)
    }

    xmin <- cells$x_num - 0.5
    xmax <- cells$x_num + 0.5
    ymin <- cells$y_num - 0.5
    ymax <- cells$y_num + 0.5
    segments <- rbind(
      data.frame(orientation = "h", x = xmin, y = ymin, xend = xmax, yend = ymin, normal_x = 0, normal_y = 1, stringsAsFactors = FALSE),
      data.frame(orientation = "h", x = xmin, y = ymax, xend = xmax, yend = ymax, normal_x = 0, normal_y = -1, stringsAsFactors = FALSE),
      data.frame(orientation = "v", x = xmin, y = ymin, xend = xmin, yend = ymax, normal_x = 1, normal_y = 0, stringsAsFactors = FALSE),
      data.frame(orientation = "v", x = xmax, y = ymin, xend = xmax, yend = ymax, normal_x = -1, normal_y = 0, stringsAsFactors = FALSE)
    )
    segments$key <- paste(
      segments$orientation,
      sprintf("%.6f", segments$x),
      sprintf("%.6f", segments$y),
      sprintf("%.6f", segments$xend),
      sprintf("%.6f", segments$yend),
      sep = "::"
    )
    key_counts <- table(segments$key)
    segments <- segments[key_counts[segments$key] == 1L, , drop = FALSE]
    if (nrow(segments) == 0L) {
      return(empty_border_segments)
    }
    segments$frame_label <- frame_label
    merge_collinear_segments(segments)
  }

  offset_overlapping_edges <- function(segments, base_segments, frame_label, offset = 0.04) {
    use_segments <- segments[as.character(segments$frame_label) == frame_label, , drop = FALSE]
    other_segments <- segments[as.character(segments$frame_label) != frame_label, , drop = FALSE]
    if (nrow(use_segments) == 0L || nrow(base_segments) == 0L) {
      return(segments)
    }

    split_segments <- do.call(rbind, lapply(seq_len(nrow(use_segments)), function(i) {
      segment <- use_segments[i, , drop = FALSE]
      overlaps <- base_segments[
        base_segments$orientation == segment$orientation &
          abs(base_segments$fixed - segment$fixed) < 1e-8 &
          pmax(base_segments$start, segment$start) < pmin(base_segments$end, segment$end) - 1e-8,
        ,
        drop = FALSE
      ]
      if (nrow(overlaps) == 0L) {
        segment$overlaps_base <- FALSE
        return(segment)
      }

      breaks <- sort(unique(c(
        segment$start,
        segment$end,
        pmax(overlaps$start, segment$start),
        pmin(overlaps$end, segment$end)
      )))
      pieces <- lapply(seq_len(length(breaks) - 1L), function(j) {
        piece_start <- breaks[[j]]
        piece_end <- breaks[[j + 1L]]
        if (piece_end <= piece_start + 1e-8) {
          return(NULL)
        }
        piece <- segment
        piece$start <- piece_start
        piece$end <- piece_end
        piece$overlaps_base <- any(
          pmax(overlaps$start, piece_start) < pmin(overlaps$end, piece_end) - 1e-8
        )
        piece
      })
      do.call(rbind, Filter(Negate(is.null), pieces))
    }))

    split_segments$fixed <- ifelse(
      split_segments$overlaps_base,
      split_segments$fixed + ifelse(split_segments$orientation == "h", split_segments$normal_y, split_segments$normal_x) * offset,
      split_segments$fixed
    )
    split_segments$shifted <- split_segments$overlaps_base
    split_segments$x <- ifelse(split_segments$orientation == "h", split_segments$start, split_segments$fixed)
    split_segments$y <- ifelse(split_segments$orientation == "h", split_segments$fixed, split_segments$start)
    split_segments$xend <- ifelse(split_segments$orientation == "h", split_segments$end, split_segments$fixed)
    split_segments$yend <- ifelse(split_segments$orientation == "h", split_segments$fixed, split_segments$end)
    split_segments$overlaps_base <- NULL

    shifted_segments <- merge_collinear_segments(split_segments)
    shifted_idx <- isTRUE(shifted_segments$shifted) | shifted_segments$shifted
    if (any(shifted_idx)) {
      segment_length <- shifted_segments$end[shifted_idx] - shifted_segments$start[shifted_idx]
      shrink <- pmin(offset * 1.45, pmax(segment_length / 2 - 1e-6, 0))
      shifted_segments$start[shifted_idx] <- shifted_segments$start[shifted_idx] + shrink
      shifted_segments$end[shifted_idx] <- shifted_segments$end[shifted_idx] - shrink
      shifted_segments$x[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$start[shifted_idx],
        shifted_segments$fixed[shifted_idx]
      )
      shifted_segments$y[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$fixed[shifted_idx],
        shifted_segments$start[shifted_idx]
      )
      shifted_segments$xend[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$end[shifted_idx],
        shifted_segments$fixed[shifted_idx]
      )
      shifted_segments$yend[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$fixed[shifted_idx],
        shifted_segments$end[shifted_idx]
      )
    }
    shifted_segments$shifted <- NULL

    rbind(other_segments, shifted_segments)
  }

  same_condition_cells <- pairwise[pairwise$same_condition, , drop = FALSE]
  exact_sample_cells <- pairwise[pairwise$same_sample, , drop = FALSE]
  same_condition_segments <- make_border_segments(same_condition_cells, "same condition")
  exact_sample_segments <- make_border_segments(exact_sample_cells, "exact sample")
  border_segments <- rbind(
    same_condition_segments,
    exact_sample_segments
  )
  border_segments <- offset_overlapping_edges(border_segments, same_condition_segments, "exact sample")
  border_segments$frame_label <- factor(border_segments$frame_label, levels = names(frame_colors))
  frame_count <- length(frame_colors)
  matrix_legend_height_mm <- max(52, min(76, length(sample_levels) * 14))
  legend_block_gap_mm <- max(1.5, min(4, matrix_legend_height_mm * 0.045))
  frame_gap_mm <- max(1, min(3, matrix_legend_height_mm * 0.025))
  fill_barheight_mm <- max(14, min(20, matrix_legend_height_mm * 0.25))
  frame_keyheight_mm <- max(
    18,
    min(
      28,
      (matrix_legend_height_mm * 0.72 - frame_gap_mm * max(frame_count - 1, 0)) /
        max(frame_count, 1)
    )
  )
  frame_legend_data <- data.frame(
    x = rep(pairwise$x_num[[1]], length(frame_colors)),
    y = rep(pairwise$y_num[[1]], length(frame_colors)),
    frame_label = factor(names(frame_colors), levels = names(frame_colors)),
    stringsAsFactors = FALSE
  )
  method_label <- switch(
    activity_cross_sample_similarity$method,
    pearson = "Pearson r",
    spearman = "Spearman rho",
    kendall = "Kendall tau",
    "Correlation"
  )
  result_labels <- activity_cross_sample_similarity$result_labels

  subtitle <- sprintf(
    "%s correlation between %s samples and %s samples (%d shared TE subfamilies).",
    comparison_pretty_label(activity_cross_sample_similarity$method),
    result_labels[[1]],
    result_labels[[2]],
    length(activity_cross_sample_similarity$shared_te_subfamilies)
  )
  if (!is.null(activity_cross_sample_similarity$comparison_class)) {
    subtitle <- comparison_plot_subtitle(subtitle, activity_cross_sample_similarity$comparison_class)
  }

  plot <- ggplot2::ggplot(pairwise, ggplot2::aes(x = x_profile, y = y_profile, fill = plot_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::geom_segment(
      data = border_segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = frame_label),
      inherit.aes = FALSE,
      linewidth = 1.1,
      lineend = "butt",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = frame_legend_data,
      ggplot2::aes(x = x, y = y, color = frame_label),
      inherit.aes = FALSE,
      shape = 22,
      fill = "white",
      size = 3.2,
      stroke = 1.25,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), na.value = "#f2f4f7") +
    ggplot2::scale_color_manual(
      values = frame_colors,
      labels = c("Same condition", "Exact same sample"),
      drop = FALSE
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Cross-sample activity-score similarity",
      subtitle = subtitle,
      x = result_labels[[2]],
      y = result_labels[[1]],
      fill = method_label,
      color = "Frame"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        order = 1,
        title.position = "top",
        barheight = grid::unit(fill_barheight_mm, "mm"),
        barwidth = grid::unit(4, "mm")
      ),
      color = ggplot2::guide_legend(
        order = 2,
        title.position = "top",
        direction = "vertical",
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(frame_keyheight_mm, "mm"),
        keywidth = grid::unit(4.5, "mm"),
        override.aes = list(fill = "white", alpha = 1, shape = 22, size = 3.2, stroke = 1.25),
        theme = ggplot2::theme(
          legend.text = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8.5),
          legend.key.width = grid::unit(4.5, "mm"),
          legend.key.height = grid::unit(frame_keyheight_mm, "mm"),
          legend.spacing.y = grid::unit(frame_gap_mm, "mm")
        )
      )
    ) +
    comparison_theme() +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 8.5, lineheight = 0.95),
      legend.key.width = grid::unit(4.5, "mm"),
      legend.key.height = grid::unit(frame_keyheight_mm, "mm"),
      legend.spacing.y = grid::unit(frame_gap_mm, "mm"),
      legend.box.spacing = grid::unit(legend_block_gap_mm, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 4),
      axis.text.y = ggplot2::element_text(angle = 35, hjust = 1, vjust = 0.5),
      plot.margin = ggplot2::margin(5.5, 4, 5.5, 5.5)
    )

  if (isTRUE(show_values)) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes(label = tile_label), size = text_size, lineheight = 0.9)
  }
  plot
}

plot_activity_top_change <- function(comparison_results, target_tes = NULL, top_k = 25L, use_abs = FALSE, comparison_class = NULL) {
  comparison_need_ggplot()
  if (isTRUE(use_abs)) {
    comparison_warn("`use_abs = TRUE` is deprecated for `plot_activity_top_change()` and ignored; signed mean activity is plotted.")
  }
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  summary <- comparison_activity_matrices_shared_tes(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }

  totals <- do.call(rbind, lapply(names(summary$matrices), function(result_id) {
    mat <- summary$matrices[[result_id]]
    selection_value <- rowMeans(abs(mat), na.rm = TRUE)
    activity_value <- rowMeans(mat, na.rm = TRUE)
    data.frame(
      result_id = result_id,
      te_subfamily = rownames(mat),
      selection_magnitude = selection_value,
      activity_magnitude = activity_value,
      n_activity_columns = ncol(mat),
      stringsAsFactors = FALSE
    )
  }))

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  totals$result_label <- unname(label_map[totals$result_id])
  totals$result_label[is.na(totals$result_label) | !nzchar(totals$result_label)] <- totals$result_id[is.na(totals$result_label) | !nzchar(totals$result_label)]

  selection_stats <- aggregate(selection_magnitude ~ te_subfamily, data = totals, FUN = max)
  colnames(selection_stats)[2] <- "max_selection_magnitude"

  te_stats <- aggregate(activity_magnitude ~ te_subfamily, data = totals, FUN = function(x) {
    c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
  })
  te_stats <- data.frame(
    te_subfamily = te_stats$te_subfamily,
    min_activity = te_stats$activity_magnitude[, "min"],
    max_activity = te_stats$activity_magnitude[, "max"],
    mean_activity = te_stats$activity_magnitude[, "mean"],
    stringsAsFactors = FALSE
  )
  te_stats <- merge(te_stats, selection_stats, by = "te_subfamily", all.x = TRUE, sort = FALSE)
  te_stats$activity_range <- te_stats$max_activity - te_stats$min_activity
  te_stats <- te_stats[order(te_stats$max_selection_magnitude, te_stats$activity_range, decreasing = TRUE), , drop = FALSE]
  keep <- head(te_stats$te_subfamily, top_k)
  te_stats <- te_stats[te_stats$te_subfamily %in% keep, , drop = FALSE]
  totals <- totals[totals$te_subfamily %in% keep, , drop = FALSE]

  levels_te <- rev(te_stats$te_subfamily[order(te_stats$max_selection_magnitude, decreasing = TRUE)])
  totals$te_subfamily <- factor(totals$te_subfamily, levels = levels_te)
  te_stats$te_subfamily <- factor(te_stats$te_subfamily, levels = levels_te)
  totals$result_label <- factor(totals$result_label, levels = unname(label_map))

  x_max <- max(c(totals$activity_magnitude, te_stats$max_activity), na.rm = TRUE)
  x_min <- min(c(totals$activity_magnitude, te_stats$min_activity), na.rm = TRUE)
  x_lower <- min(0, x_min)
  x_upper <- max(0, x_max)
  x_span <- max(x_upper - x_lower, 0.01)
  axis_pad <- x_span * 0.08
  label_pad <- x_span * 0.05
  delta_x <- x_upper + label_pad
  te_stats$delta_x <- delta_x
  te_stats$range_label <- sprintf("%.3g", te_stats$activity_range)
  delta_header <- data.frame(
    delta_x = delta_x,
    te_subfamily = factor(levels_te[length(levels_te)], levels = levels_te),
    label = "max. d.",
    stringsAsFactors = FALSE
  )
  result_labels <- levels(totals$result_label)
  color_values <- c(
    `delta max` = "#98a2b3",
    stats::setNames(grDevices::hcl.colors(length(result_labels), palette = "Dark 3"), result_labels)
  )
  activity_label <- "mean signed activity"

  ggplot2::ggplot() +
    comparison_zero_reference_line() +
    ggplot2::geom_segment(
      data = te_stats,
      ggplot2::aes(x = min_activity, xend = max_activity, y = te_subfamily, yend = te_subfamily, color = "delta max"),
      linewidth = 1.15,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = totals,
      ggplot2::aes(x = activity_magnitude, y = te_subfamily, color = result_label),
      size = 2.4,
      alpha = 0.9
    ) +
    ggplot2::geom_text(
      data = te_stats,
      ggplot2::aes(x = delta_x, y = te_subfamily, label = range_label),
      hjust = 0,
      color = "#667085",
      size = 2.8
    ) +
    ggplot2::geom_text(
      data = delta_header,
      ggplot2::aes(x = delta_x, y = te_subfamily, label = label),
      hjust = 0,
      vjust = 0,
      nudge_y = 0.72,
      color = "#475467",
      fontface = "bold",
      size = 2.8
    ) +
    ggplot2::scale_x_continuous(limits = c(x_lower - axis_pad, x_upper + axis_pad)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5, 1.0))) +
    ggplot2::scale_color_manual(values = color_values, breaks = names(color_values), drop = FALSE) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Top activity magnitude and change",
      subtitle = comparison_plot_subtitle(
        sprintf("Top %s shared TEs; ranked by max mean abs(activity).", length(keep)),
        comparison_class
      ),
      x = paste0(comparison_pretty_label(activity_label), " across available samples or conditions"),
      y = "TE subfamily",
      color = NULL
    ) +
    comparison_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5.5, 28, 5.5, 5.5)
    )
}

comparison_available_te_contrasts <- function(comparison_results) {
  has_contrast <- vapply(comparison_results, function(result) !is.null(result$tables$te_contrast_summary), logical(1))
  if (!all(has_contrast)) {
    return(character())
  }
  common_cols <- Reduce(intersect, lapply(comparison_results, function(result) colnames(result$tables$te_contrast_summary)))
  sub("^padj_", "", grep("^padj_", common_cols, value = TRUE))
}

comparison_clean_contrast_id <- function(contrast_id) {
  contrast_id <- as.character(contrast_id)
  contrast_id <- sub("^padj_", "", contrast_id)
  contrast_id <- sub("^delta_", "", contrast_id)
  contrast_id <- sub("^signif_", "", contrast_id)
  contrast_id
}

comparison_resolve_te_contrast <- function(comparison_results, contrast_id = NULL) {
  available <- comparison_available_te_contrasts(comparison_results)
  if (length(available) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared adjusted p-value contrast columns matching `padj_*` are available in te_contrast_summary.tsv.",
      available_contrasts = available
    ))
  }

  if (is.null(contrast_id) || length(contrast_id) == 0L || is.na(contrast_id[1]) || !nzchar(as.character(contrast_id[1]))) {
    contrast_id <- available[1]
  }
  contrast_id <- comparison_clean_contrast_id(contrast_id[1])
  if (!contrast_id %in% available) {
    return(list(
      status = "skipped",
      reason = sprintf(
        "Selected contrast `%s` is not shared by all result directories. Available shared contrasts: %s",
        contrast_id,
        paste(available, collapse = ", ")
      ),
      available_contrasts = available
    ))
  }

  parts <- strsplit(contrast_id, "_vs_", fixed = TRUE)[[1]]
  condition_a <- if (length(parts) >= 1L) parts[1] else NA_character_
  condition_b <- if (length(parts) >= 2L) paste(parts[-1], collapse = "_vs_") else NA_character_
  list(
    status = "ok",
    contrast_id = contrast_id,
    available_contrasts = available,
    condition_a = condition_a,
    condition_b = condition_b,
    padj_col = paste0("padj_", contrast_id),
    delta_col = paste0("delta_", contrast_id),
    signif_col = paste0("signif_", contrast_id),
    mean_a_col = if (!is.na(condition_a) && nzchar(condition_a)) paste0("mean_", condition_a) else NA_character_,
    mean_b_col = if (!is.na(condition_b) && nzchar(condition_b)) paste0("mean_", condition_b) else NA_character_
  )
}

comparison_fisher_p <- function(p_values) {
  p_values <- suppressWarnings(as.numeric(p_values))
  p_values <- p_values[!is.na(p_values) & p_values >= 0 & p_values <= 1]
  if (length(p_values) == 0L) {
    return(NA_real_)
  }
  p_values <- pmax(p_values, .Machine$double.xmin)
  stats::pchisq(-2 * sum(log(p_values)), df = 2 * length(p_values), lower.tail = FALSE)
}

comparison_as_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

comparison_fixed_contrast_class <- function(delta, padj, alpha, condition_a = NA_character_, condition_b = NA_character_) {
  delta <- suppressWarnings(as.numeric(delta))
  padj <- suppressWarnings(as.numeric(padj))
  condition_a <- if (!is.na(condition_a) && nzchar(condition_a)) condition_a else "positive delta"
  condition_b <- if (!is.na(condition_b) && nzchar(condition_b)) condition_b else "negative delta"
  out <- rep("flat_or_noisy", length(delta))
  is_signal <- !is.na(delta) & !is.na(padj) & padj <= alpha
  out[is_signal & delta > .Machine$double.eps] <- paste(condition_a, "higher")
  out[is_signal & delta < -.Machine$double.eps] <- paste(condition_b, "higher")
  out
}

comparison_format_p_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(signif(x, 3), format = "fg", digits = 3))
  )
}

comparison_format_compact_number <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "NA", formatC(signif(x, 3), format = "fg", digits = 3))
}

comparison_available_te_delta_contrasts <- function(comparison_results) {
  has_contrast <- vapply(comparison_results, function(result) !is.null(result$tables$te_contrast_summary), logical(1))
  if (!all(has_contrast)) {
    return(character())
  }
  common_cols <- Reduce(intersect, lapply(comparison_results, function(result) colnames(result$tables$te_contrast_summary)))
  padj_contrasts <- sub("^padj_", "", grep("^padj_", common_cols, value = TRUE))
  delta_contrasts <- sub("^delta_", "", grep("^delta_", common_cols, value = TRUE))
  intersect(padj_contrasts, delta_contrasts)
}

summarise_activity_contrast_change <- function(comparison_results,
                                               contrasts = NULL,
                                               target_tes = NULL,
                                               top_k = 25L,
                                               comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  has_contrast <- vapply(comparison_results, function(result) !is.null(result$tables$te_contrast_summary), logical(1))
  if (!all(has_contrast)) {
    missing_ids <- comparison_result_ids(comparison_results)[!has_contrast]
    return(list(
      status = "skipped",
      reason = sprintf("Missing te_contrast_summary.tsv in: %s", paste(missing_ids, collapse = ", "))
    ))
  }

  available <- comparison_available_te_delta_contrasts(comparison_results)
  if (length(available) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared contrast columns with both `delta_*` and `padj_*` are available in te_contrast_summary.tsv."
    ))
  }

  if (is.null(contrasts) || length(contrasts) == 0L) {
    contrasts <- available
  } else {
    contrasts <- unique(vapply(contrasts, comparison_clean_contrast_id, character(1)))
    missing_contrasts <- setdiff(contrasts, available)
    contrasts <- intersect(contrasts, available)
    if (length(contrasts) == 0L) {
      return(list(
        status = "skipped",
        reason = sprintf(
          "None of the requested contrasts are shared with both `delta_*` and `padj_*` columns. Available contrasts: %s",
          paste(available, collapse = ", ")
        ),
        available_contrasts = available,
        missing_contrasts = missing_contrasts
      ))
    }
  }

  contrast_sets <- lapply(comparison_results, function(result) unique(result$tables$te_contrast_summary$te_subfamily))
  shared_tes <- Reduce(intersect, contrast_sets)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(status = "skipped", reason = "No shared TE subfamilies are available for contrast-specific change plots."))
  }

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  result_ids <- comparison_result_ids(comparison_results)
  result_levels <- unname(label_map[result_ids])
  result_levels[is.na(result_levels) | !nzchar(result_levels)] <- comparison_short_result_label(result_ids[is.na(result_levels) | !nzchar(result_levels)])

  contrast_summaries <- lapply(contrasts, function(contrast_id) {
    delta_col <- paste0("delta_", contrast_id)
    totals <- do.call(rbind, lapply(comparison_results, function(result) {
      x <- result$tables$te_contrast_summary
      x <- x[x$te_subfamily %in% shared_tes, , drop = FALSE]
      x <- x[match(shared_tes, x$te_subfamily), , drop = FALSE]
      delta <- suppressWarnings(as.numeric(x[[delta_col]]))
      result_label <- unname(label_map[result$result_id])
      if (is.na(result_label) || !nzchar(result_label)) {
        result_label <- comparison_short_result_label(result$result_id)
      }
      data.frame(
        result_id = result$result_id,
        result_label = result_label,
        contrast_id = contrast_id,
        te_subfamily = x$te_subfamily,
        contrast_delta = delta,
        selection_magnitude = abs(delta),
        stringsAsFactors = FALSE
      )
    }))
    totals <- totals[!is.na(totals$contrast_delta), , drop = FALSE]
    if (nrow(totals) == 0L) {
      return(NULL)
    }

    selection_stats <- aggregate(selection_magnitude ~ te_subfamily, data = totals, FUN = max)
    colnames(selection_stats)[2] <- "max_selection_magnitude"
    te_stats <- aggregate(contrast_delta ~ te_subfamily, data = totals, FUN = function(x) {
      c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
    })
    te_stats <- data.frame(
      te_subfamily = te_stats$te_subfamily,
      min_delta = te_stats$contrast_delta[, "min"],
      max_delta = te_stats$contrast_delta[, "max"],
      mean_delta = te_stats$contrast_delta[, "mean"],
      stringsAsFactors = FALSE
    )
    te_stats <- merge(te_stats, selection_stats, by = "te_subfamily", all.x = TRUE, sort = FALSE)
    te_stats$delta_range <- te_stats$max_delta - te_stats$min_delta
    te_stats <- te_stats[order(te_stats$max_selection_magnitude, te_stats$delta_range, decreasing = TRUE), , drop = FALSE]
    keep <- head(te_stats$te_subfamily, top_k)
    te_stats <- te_stats[te_stats$te_subfamily %in% keep, , drop = FALSE]
    totals <- totals[totals$te_subfamily %in% keep, , drop = FALSE]

    levels_te <- rev(te_stats$te_subfamily[order(te_stats$max_selection_magnitude, decreasing = TRUE)])
    totals$te_subfamily <- factor(totals$te_subfamily, levels = levels_te)
    te_stats$te_subfamily <- factor(te_stats$te_subfamily, levels = levels_te)
    totals$result_label <- factor(totals$result_label, levels = result_levels)

    list(
      contrast_id = contrast_id,
      totals = totals,
      te_stats = te_stats,
      levels_te = levels_te,
      n_selected = length(keep)
    )
  })
  names(contrast_summaries) <- contrasts
  contrast_summaries <- contrast_summaries[!vapply(contrast_summaries, is.null, logical(1))]

  if (length(contrast_summaries) == 0L) {
    return(list(status = "skipped", reason = "No finite contrast delta values are available after filtering."))
  }

  list(
    status = "ok",
    contrasts = names(contrast_summaries),
    available_contrasts = available,
    shared_te_subfamilies = shared_tes,
    top_k = top_k,
    result_levels = result_levels,
    contrast_summaries = contrast_summaries,
    comparison_class = comparison_class
  )
}

plot_activity_contrast_change <- function(contrast_change_summary, contrast_id) {
  comparison_need_ggplot()
  contrast_data <- contrast_change_summary$contrast_summaries[[contrast_id]]
  if (is.null(contrast_data)) {
    comparison_warn(sprintf("No contrast-change summary is available for `%s`.", contrast_id))
    return(invisible(NULL))
  }

  totals <- contrast_data$totals
  te_stats <- contrast_data$te_stats
  levels_te <- contrast_data$levels_te

  x_max <- max(c(totals$contrast_delta, te_stats$max_delta), na.rm = TRUE)
  x_min <- min(c(totals$contrast_delta, te_stats$min_delta), na.rm = TRUE)
  x_lower <- min(0, x_min)
  x_upper <- max(0, x_max)
  x_span <- max(x_upper - x_lower, 0.01)
  axis_pad <- x_span * 0.08
  label_pad <- x_span * 0.05
  delta_x <- x_upper + label_pad
  te_stats$delta_x <- delta_x
  te_stats$range_label <- sprintf("%.3g", te_stats$delta_range)
  delta_header <- data.frame(
    delta_x = delta_x,
    te_subfamily = factor(levels_te[length(levels_te)], levels = levels_te),
    label = "max. d.",
    stringsAsFactors = FALSE
  )

  color_values <- c(
    `delta max` = "#98a2b3",
    stats::setNames(grDevices::hcl.colors(length(contrast_change_summary$result_levels), palette = "Dark 3"), contrast_change_summary$result_levels)
  )

  ggplot2::ggplot() +
    comparison_zero_reference_line() +
    ggplot2::geom_segment(
      data = te_stats,
      ggplot2::aes(x = min_delta, xend = max_delta, y = te_subfamily, yend = te_subfamily, color = "delta max"),
      linewidth = 1.15,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = totals,
      ggplot2::aes(x = contrast_delta, y = te_subfamily, color = result_label),
      size = 2.25,
      alpha = 0.9
    ) +
    ggplot2::geom_text(
      data = te_stats,
      ggplot2::aes(x = delta_x, y = te_subfamily, label = range_label),
      hjust = 0,
      color = "#667085",
      size = 2.55
    ) +
    ggplot2::geom_text(
      data = delta_header,
      ggplot2::aes(x = delta_x, y = te_subfamily, label = label),
      hjust = 0,
      vjust = 0,
      nudge_y = 0.72,
      color = "#475467",
      fontface = "bold",
      size = 2.55
    ) +
    ggplot2::scale_x_continuous(limits = c(x_lower - axis_pad, x_upper + axis_pad)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5, 1.0))) +
    ggplot2::scale_color_manual(values = color_values, breaks = names(color_values), drop = FALSE) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = contrast_id,
      subtitle = sprintf("Top %s TEs by max abs(delta).", contrast_data$n_selected),
      x = "signed contrast delta",
      y = "TE subfamily",
      color = NULL
    ) +
    comparison_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = 9),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5.5, 28, 5.5, 5.5)
    )
}

print_activity_contrast_change_plots <- function(contrast_change_summary) {
  comparison_need_ggplot()
  if (!identical(contrast_change_summary$status, "ok")) {
    comparison_warn(contrast_change_summary$reason)
    return(invisible(contrast_change_summary))
  }

  plots <- lapply(contrast_change_summary$contrasts, function(contrast_id) {
    plot_activity_contrast_change(contrast_change_summary, contrast_id)
  })
  names(plots) <- contrast_change_summary$contrasts
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (length(plots) == 0L) {
    comparison_warn("No contrast-change plots could be built.")
    return(invisible(contrast_change_summary))
  }

  for (plot in plots) {
    print(plot)
  }
  invisible(plots)
}

summarise_fixed_contrast_te_evidence <- function(comparison_results,
                                                 contrast_id = NULL,
                                                 selected_te_space = NULL,
                                                 target_tes = NULL,
                                                 top_k = 30L,
                                                 alpha = 0.05,
                                                 comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  if (!is.null(selected_te_space) && identical(selected_te_space$status, "ok")) {
    target_tes <- selected_te_space$te_subfamilies
  }

  resolved <- comparison_resolve_te_contrast(comparison_results, contrast_id = contrast_id)
  if (!identical(resolved$status, "ok")) {
    return(resolved)
  }

  has_contrast <- vapply(comparison_results, function(result) !is.null(result$tables$te_contrast_summary), logical(1))
  if (!all(has_contrast)) {
    missing_ids <- comparison_result_ids(comparison_results)[!has_contrast]
    return(list(
      status = "skipped",
      reason = sprintf("Missing te_contrast_summary.tsv in: %s", paste(missing_ids, collapse = ", "))
    ))
  }

  contrast_sets <- lapply(comparison_results, function(result) unique(result$tables$te_contrast_summary$te_subfamily))
  shared_tes <- Reduce(intersect, contrast_sets)
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }
  if (length(shared_tes) == 0L) {
    return(list(status = "skipped", reason = "No shared TE subfamilies are available for the selected contrast."))
  }

  common_cols <- Reduce(intersect, lapply(comparison_results, function(result) colnames(result$tables$te_contrast_summary)))
  required_cols <- c("te_subfamily", resolved$padj_col, resolved$delta_col)
  missing_cols <- setdiff(required_cols, common_cols)
  if (length(missing_cols) > 0L) {
    return(list(
      status = "skipped",
      reason = sprintf("Selected contrast `%s` is missing required shared column(s): %s", resolved$contrast_id, paste(missing_cols, collapse = ", "))
    ))
  }

  has_activity_means <- all(c(resolved$mean_a_col, resolved$mean_b_col) %in% common_cols)
  has_signif <- resolved$signif_col %in% common_cols
  has_trend <- "trend_class" %in% common_cols
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)

  result_summary <- do.call(rbind, lapply(comparison_results, function(result) {
    x <- result$tables$te_contrast_summary
    x <- x[x$te_subfamily %in% shared_tes, , drop = FALSE]
    x <- x[match(shared_tes, x$te_subfamily), , drop = FALSE]
    padj <- suppressWarnings(as.numeric(x[[resolved$padj_col]]))
    delta <- suppressWarnings(as.numeric(x[[resolved$delta_col]]))
    is_significant <- if (isTRUE(has_signif)) comparison_as_logical(x[[resolved$signif_col]]) else padj <= alpha

    out <- data.frame(
      result_id = result$result_id,
      result_label = unname(label_map[result$result_id]),
      te_subfamily = x$te_subfamily,
      contrast_id = resolved$contrast_id,
      mean_a = if (isTRUE(has_activity_means)) suppressWarnings(as.numeric(x[[resolved$mean_a_col]])) else NA_real_,
      mean_b = if (isTRUE(has_activity_means)) suppressWarnings(as.numeric(x[[resolved$mean_b_col]])) else NA_real_,
      delta = delta,
      padj = padj,
      is_significant = is_significant,
      trend_class = if (isTRUE(has_trend)) as.character(x$trend_class) else NA_character_,
      stringsAsFactors = FALSE
    )
    out$result_label[is.na(out$result_label) | !nzchar(out$result_label)] <- result$result_id
    out$contrast_class <- comparison_fixed_contrast_class(
      delta = out$delta,
      padj = out$padj,
      alpha = alpha,
      condition_a = resolved$condition_a,
      condition_b = resolved$condition_b
    )
    out
  }))
  rownames(result_summary) <- NULL

  te_summary <- do.call(rbind, lapply(split(result_summary, result_summary$te_subfamily), function(x) {
    classes <- x$contrast_class[!is.na(x$contrast_class) & nzchar(x$contrast_class)]
    class_counts <- if (length(classes) > 0L) sort(table(classes), decreasing = TRUE) else integer()
    top_count <- if (length(class_counts) > 0L) max(as.integer(class_counts)) else 0L
    top_classes <- if (length(class_counts) > 0L) sort(names(class_counts)[as.integer(class_counts) == top_count]) else character()
    majority_class <- if (length(top_classes) == 0L) {
      "unclassified"
    } else if (top_count > length(classes) / 2) {
      top_classes[1]
    } else {
      "mixed"
    }
    noise_votes <- sum(classes == "flat_or_noisy")
    non_noise_votes <- sum(classes != "flat_or_noisy")
    marker_status <- if (length(classes) == 0L) {
      "no class votes"
    } else if (noise_votes > length(classes) / 2) {
      "majority noise"
    } else {
      "not majority noise"
    }
    activity_values <- c(x$mean_a, x$mean_b)
    activity_values <- activity_values[!is.na(activity_values)]
    data.frame(
      te_subfamily = x$te_subfamily[1],
      aggregated_padj = comparison_fisher_p(x$padj),
      mean_padj = if (all(is.na(x$padj))) NA_real_ else mean(x$padj, na.rm = TRUE),
      min_padj = if (all(is.na(x$padj))) NA_real_ else min(x$padj, na.rm = TRUE),
      n_result_padj = sum(!is.na(x$padj)),
      mean_delta = if (all(is.na(x$delta))) NA_real_ else mean(x$delta, na.rm = TRUE),
      min_delta = if (all(is.na(x$delta))) NA_real_ else min(x$delta, na.rm = TRUE),
      max_delta = if (all(is.na(x$delta))) NA_real_ else max(x$delta, na.rm = TRUE),
      n_results = length(unique(x$result_id)),
      n_significant_results = sum(x$is_significant, na.rm = TRUE),
      majority_contrast_class = majority_class,
      majority_class_votes = top_count,
      n_noise_votes = noise_votes,
      n_non_noise_votes = non_noise_votes,
      te_marker_status = marker_status,
      min_activity = if (length(activity_values) == 0L) NA_real_ else min(activity_values, na.rm = TRUE),
      max_activity = if (length(activity_values) == 0L) NA_real_ else max(activity_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(te_summary) <- NULL

  global_evidence <- summarise_te_padj_evidence(
    comparison_results = comparison_results,
    target_tes = shared_tes,
    p_score_mode = "mean_padj",
    alpha = alpha,
    comparison_class = comparison_class
  )
  if (identical(global_evidence$status, "ok")) {
    global_cols <- intersect(
      c(
        "te_subfamily",
        "global_assigned_class",
        "global_class_votes",
        "global_non_noise_votes",
        "global_class_tie",
        "classification_agreement",
        "global_combined_raw_p_fisher_bh_evidence"
      ),
      colnames(global_evidence$te_summary)
    )
    global_summary <- global_evidence$te_summary[, global_cols, drop = FALSE]
    te_summary <- merge(te_summary, global_summary, by = "te_subfamily", all.x = TRUE, sort = FALSE)
  } else {
    te_summary$global_assigned_class <- NA_character_
    te_summary$global_class_votes <- NA_integer_
    te_summary$global_non_noise_votes <- NA_integer_
    te_summary$global_class_tie <- NA
    te_summary$classification_agreement <- NA_real_
    te_summary$global_combined_raw_p_fisher_bh_evidence <- NA_real_
  }
  if (!"global_assigned_class" %in% colnames(te_summary)) {
    te_summary$global_assigned_class <- NA_character_
  }
  te_summary$global_assigned_class[is.na(te_summary$global_assigned_class) | !nzchar(te_summary$global_assigned_class)] <- "unclassified"

  te_summary <- te_summary[order(te_summary$aggregated_padj, te_summary$mean_padj, -abs(te_summary$mean_delta), te_summary$te_subfamily, na.last = TRUE), , drop = FALSE]

  keep <- head(te_summary$te_subfamily, top_k)
  selected_te_summary <- te_summary[te_summary$te_subfamily %in% keep, , drop = FALSE]
  result_summary_top <- result_summary[result_summary$te_subfamily %in% keep, , drop = FALSE]

  list(
    status = "ok",
    contrast_id = resolved$contrast_id,
    condition_a = resolved$condition_a,
    condition_b = resolved$condition_b,
    alpha = alpha,
    top_k = top_k,
    available_contrasts = resolved$available_contrasts,
    has_activity_means = has_activity_means,
    result_summary = result_summary,
    te_summary = te_summary,
    selected_te_summary = selected_te_summary,
    result_summary_top = result_summary_top,
    comparison_class = comparison_class
  )
}

comparison_fixed_contrast_plot_theme <- function() {
  comparison_need_ggplot()
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = 0),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.height = grid::unit(3.5, "mm"),
      plot.margin = ggplot2::margin(3, 5, 3, 5)
    )
}

plot_fixed_contrast_te_evidence <- function(fixed_contrast_summary) {
  comparison_need_ggplot()
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    comparison_warn("Package `patchwork` is required for the fixed-contrast TE evidence panel.")
    return(invisible(fixed_contrast_summary))
  }
  if (!identical(fixed_contrast_summary$status, "ok")) {
    comparison_warn(fixed_contrast_summary$reason)
    return(invisible(fixed_contrast_summary))
  }

  te_summary <- fixed_contrast_summary$selected_te_summary
  result_summary <- fixed_contrast_summary$result_summary_top
  if (nrow(te_summary) == 0L || nrow(result_summary) == 0L) {
    comparison_warn("No TE subfamilies remain for fixed-contrast plotting.")
    return(invisible(fixed_contrast_summary))
  }

  te_summary$y <- rev(seq_len(nrow(te_summary)))
  y_lookup <- te_summary[, c("te_subfamily", "y"), drop = FALSE]
  result_order <- unique(result_summary$result_id)
  result_levels <- vapply(result_order, function(result_id) {
    label <- result_summary$result_label[result_summary$result_id == result_id]
    label <- label[!is.na(label) & nzchar(label)]
    if (length(label) == 0L) result_id else label[1]
  }, character(1))
  result_summary <- merge(result_summary, y_lookup, by = "te_subfamily", all.x = TRUE, sort = FALSE)
  result_offsets <- if (length(result_levels) <= 1L) 0 else seq(-0.16, 0.16, length.out = length(result_levels))
  names(result_offsets) <- result_levels
  result_summary$y_result <- result_summary$y + unname(result_offsets[result_summary$result_label])

  y_limits <- c(0.5, nrow(te_summary) + 0.5)
  result_colors <- stats::setNames(grDevices::hcl.colors(length(result_levels), palette = "Dark 3"), result_levels)
  if (!"global_assigned_class" %in% colnames(te_summary)) {
    te_summary$global_assigned_class <- "unclassified"
  }
  te_summary$global_assigned_class[is.na(te_summary$global_assigned_class) | !nzchar(te_summary$global_assigned_class)] <- "unclassified"
  te_summary$label_fill_class <- as.character(te_summary$global_assigned_class)
  te_summary$label_fill_class[comparison_is_noise_class(te_summary$label_fill_class)] <- NA_character_
  global_class_colors <- comparison_trend_class_colors(te_summary$label_fill_class, include_noise = FALSE)
  global_class_breaks <- names(global_class_colors)[names(global_class_colors) %in% unique(na.omit(te_summary$label_fill_class))]
  global_class_legend_rows <- if (length(global_class_breaks) > 4L) 2L else 1L
  te_summary$label_fill_class <- factor(te_summary$label_fill_class, levels = names(global_class_colors))
  selected_contrast_significant <- if ("n_significant_results" %in% colnames(te_summary)) {
    !is.na(te_summary$n_significant_results) & te_summary$n_significant_results > 0L
  } else if ("aggregated_padj" %in% colnames(te_summary)) {
    !is.na(te_summary$aggregated_padj) & te_summary$aggregated_padj <= fixed_contrast_summary$alpha
  } else {
    rep(FALSE, nrow(te_summary))
  }
  te_summary$selected_contrast_significance <- factor(
    ifelse(selected_contrast_significant, "significant", "not significant"),
    levels = c("significant", "not significant")
  )
  selected_contrast_significance_breaks <- if (any(selected_contrast_significant)) "significant" else character()

  label_plot <- ggplot2::ggplot(te_summary, ggplot2::aes(x = 0, y = y)) +
    ggplot2::geom_label(
      ggplot2::aes(label = te_subfamily, fill = label_fill_class, linetype = selected_contrast_significance),
      hjust = 0,
      linewidth = 0.5,
      label.r = grid::unit(2, "pt"),
      label.padding = grid::unit(1.4, "pt"),
      border.colour = "#111827",
      text.colour = "#1f2933",
      size = 2.7
    ) +
    ggplot2::scale_fill_manual(
      values = global_class_colors,
      breaks = global_class_breaks,
      name = "Global trend class",
      drop = FALSE,
      na.value = "white",
      na.translate = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(significant = "solid", `not significant` = "blank"),
      breaks = selected_contrast_significance_breaks,
      labels = "significant",
      name = "Selected contrast",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = te_summary$y) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "TE subfamily", x = NULL, y = NULL) +
    comparison_fixed_contrast_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 1, 3, 5)
    )

  if (isTRUE(fixed_contrast_summary$has_activity_means)) {
    activity_long <- rbind(
      data.frame(
        result_summary[, c("te_subfamily", "result_id", "result_label", "y_result"), drop = FALSE],
        condition = fixed_contrast_summary$condition_a,
        activity_score = result_summary$mean_a,
        stringsAsFactors = FALSE
      ),
      data.frame(
        result_summary[, c("te_subfamily", "result_id", "result_label", "y_result"), drop = FALSE],
        condition = fixed_contrast_summary$condition_b,
        activity_score = result_summary$mean_b,
        stringsAsFactors = FALSE
      )
    )
    activity_plot <- ggplot2::ggplot() +
      comparison_zero_reference_line() +
      ggplot2::geom_segment(
        data = result_summary,
        ggplot2::aes(x = mean_b, xend = mean_a, y = y_result, yend = y_result, color = result_label),
        linewidth = 0.55,
        alpha = 0.45
      ) +
      ggplot2::geom_point(
        data = activity_long,
        ggplot2::aes(x = activity_score, y = y_result, color = result_label, shape = condition),
        size = 1.75,
        alpha = 0.92
      ) +
      ggplot2::scale_color_manual(values = result_colors, breaks = names(result_colors), name = "Result dir", drop = FALSE, guide = "none") +
      ggplot2::scale_shape_manual(values = c(16, 17), name = "Mean activity", drop = FALSE, guide = "none") +
      ggplot2::scale_y_continuous(limits = y_limits, breaks = te_summary$y) +
      ggplot2::labs(title = "Mean activity", x = "Activity score", y = NULL) +
      comparison_fixed_contrast_plot_theme() +
      ggplot2::theme(plot.margin = ggplot2::margin(3, 18, 3, 2))
  } else {
    activity_plot <- ggplot2::ggplot(te_summary, ggplot2::aes(x = 0, y = y)) +
      ggplot2::geom_text(label = "mean_* columns missing", size = 2.7, color = "#667085") +
      ggplot2::scale_y_continuous(limits = y_limits, breaks = te_summary$y) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::labs(title = "Mean activity", x = NULL, y = NULL) +
      comparison_fixed_contrast_plot_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(3, 18, 3, 2)
      )
  }

  delta_plot <- ggplot2::ggplot() +
    comparison_zero_reference_line() +
    ggplot2::geom_segment(
      data = te_summary,
      ggplot2::aes(x = min_delta, xend = max_delta, y = y, yend = y),
      color = "#98a2b3",
      linewidth = 1.05,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = result_summary,
      ggplot2::aes(x = delta, y = y_result, color = result_label),
      size = 1.9,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = result_colors, breaks = names(result_colors), name = "Result dir", drop = FALSE, guide = "none") +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = te_summary$y) +
    ggplot2::labs(
      title = "Signed distance",
      x = sprintf("%s - %s", fixed_contrast_summary$condition_a, fixed_contrast_summary$condition_b),
      y = NULL
    ) +
    comparison_fixed_contrast_plot_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 18, 3, 12))

  te_summary$neg_log10_aggregated_padj <- comparison_neg_log10_p(te_summary$aggregated_padj)
  p_span <- max(diff(range(te_summary$neg_log10_aggregated_padj, na.rm = TRUE)), 1)
  p_pad <- max(p_span * 0.025, 0.1)
  te_summary$p_label_x <- pmax(te_summary$neg_log10_aggregated_padj - p_pad * 0.35, 0)
  p_plot <- ggplot2::ggplot(te_summary, ggplot2::aes(x = neg_log10_aggregated_padj, y = y)) +
    ggplot2::geom_point(color = "#344054", fill = "#f9fafb", shape = 21, size = 2.1, stroke = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(x = p_label_x, label = comparison_format_p_value(aggregated_padj)),
      hjust = 1.08,
      color = "#475467",
      size = 2.55
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = te_summary$y) +
    ggplot2::scale_x_continuous(limits = c(0, max(te_summary$neg_log10_aggregated_padj, na.rm = TRUE) + p_pad * 0.8), expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "Aggregated p", x = "Fisher -log10(padj)", y = NULL) +
    comparison_fixed_contrast_plot_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 14, 3, 12))

  legend_groups <- list()
  if (length(global_class_breaks) > 0L) {
    legend_groups[[length(legend_groups) + 1L]] <- data.frame(
      group_title = "Global trend class",
      type = "fill",
      label = global_class_breaks,
      fill = unname(global_class_colors[global_class_breaks]),
      color = NA_character_,
      shape = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (length(selected_contrast_significance_breaks) > 0L) {
    legend_groups[[length(legend_groups) + 1L]] <- data.frame(
      group_title = "Selected contrast",
      type = "outline",
      label = "significant",
      fill = NA_character_,
      color = NA_character_,
      shape = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  legend_groups[[length(legend_groups) + 1L]] <- data.frame(
    group_title = "Result dir",
    type = "result",
    label = names(result_colors),
    fill = NA_character_,
    color = unname(result_colors),
    shape = NA_real_,
    stringsAsFactors = FALSE
  )
  condition_shapes <- stats::setNames(c(16, 17), c(fixed_contrast_summary$condition_a, fixed_contrast_summary$condition_b))
  legend_groups[[length(legend_groups) + 1L]] <- data.frame(
    group_title = "Mean activity",
    type = "shape",
    label = names(condition_shapes),
    fill = NA_character_,
    color = "#111827",
    shape = unname(condition_shapes),
    stringsAsFactors = FALSE
  )

  legend_title_width <- 2.15
  legend_item_width <- 1.85
  legend_items_per_row <- 3L
  legend_line_gap <- 0.34
  legend_group_gap <- 0.08
  legend_items_by_group <- vector("list", length(legend_groups))
  legend_titles_by_group <- vector("list", length(legend_groups))
  current_y <- 0
  for (i in seq_along(legend_groups)) {
    group <- legend_groups[[i]]
    item_idx <- seq_len(nrow(group))
    item_line <- (item_idx - 1L) %/% legend_items_per_row
    item_col <- (item_idx - 1L) %% legend_items_per_row
    group$x_symbol <- legend_title_width + item_col * legend_item_width
    group$x_text <- group$x_symbol + 0.18
    group$y <- current_y - item_line * legend_line_gap
    legend_items_by_group[[i]] <- group
    legend_titles_by_group[[i]] <- data.frame(
      title = group$group_title[1],
      x = 0,
      y = current_y,
      stringsAsFactors = FALSE
    )
    current_y <- current_y - (max(item_line) + 1L) * legend_line_gap - legend_group_gap
  }
  legend_items <- do.call(rbind, legend_items_by_group)
  legend_titles <- do.call(rbind, legend_titles_by_group)
  legend_fill_items <- legend_items[legend_items$type == "fill", , drop = FALSE]
  legend_outline_items <- legend_items[legend_items$type == "outline", , drop = FALSE]
  legend_result_items <- legend_items[legend_items$type == "result", , drop = FALSE]
  legend_shape_items <- legend_items[legend_items$type == "shape", , drop = FALSE]
  legend_plot <- ggplot2::ggplot() +
    ggplot2::geom_text(
      data = legend_titles,
      ggplot2::aes(x = x, y = y, label = title),
      hjust = 0,
      color = "#111827",
      fontface = "bold",
      size = 2.55
    ) +
    ggplot2::geom_point(
      data = legend_fill_items,
      ggplot2::aes(x = x_symbol, y = y, fill = fill),
      shape = 22,
      size = 2.4,
      stroke = 0.45,
      color = "#98a2b3"
    ) +
    ggplot2::geom_point(
      data = legend_outline_items,
      ggplot2::aes(x = x_symbol, y = y),
      shape = 22,
      size = 2.4,
      stroke = 0.7,
      fill = "white",
      color = "#111827"
    ) +
    ggplot2::geom_segment(
      data = legend_result_items,
      ggplot2::aes(x = x_symbol - 0.08, xend = x_symbol + 0.08, y = y, yend = y, color = color),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      data = legend_result_items,
      ggplot2::aes(x = x_symbol, y = y, color = color),
      size = 1.8
    ) +
    ggplot2::geom_point(
      data = legend_shape_items,
      ggplot2::aes(x = x_symbol, y = y, shape = shape),
      color = "#111827",
      size = 1.9
    ) +
    ggplot2::geom_text(
      data = legend_items,
      ggplot2::aes(x = x_text, y = y, label = label),
      hjust = 0,
      color = "#111827",
      size = 2.55
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_shape_identity() +
    ggplot2::coord_cartesian(
      xlim = c(0, legend_title_width + legend_items_per_row * legend_item_width + 0.45),
      ylim = c(current_y + legend_group_gap - 0.06, 0.16),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 5, 0, 5))

  combined_panels <- label_plot + activity_plot + delta_plot + p_plot +
    patchwork::plot_layout(widths = c(1.25, 2.45, 1.9, 2.05), guides = "keep") +
    ggplot2::theme(legend.position = "none")

  combined <- combined_panels / legend_plot +
    patchwork::plot_layout(heights = c(1, 0.16)) +
    patchwork::plot_annotation(
      title = sprintf("Fixed contrast TE evidence: %s", fixed_contrast_summary$contrast_id),
      subtitle = comparison_plot_subtitle(
        sprintf(
          "Top %s TEs by Fisher-combined adjusted p-values.",
          nrow(te_summary)
        ),
        fixed_contrast_summary$comparison_class
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = ggplot2::element_text(color = "grey35", size = 10)
      )
    )

  combined
}

summarise_fixed_contrast_te_evidence_all <- function(comparison_results,
                                                     contrast_ids = NULL,
                                                     selected_te_space = NULL,
                                                     target_tes = NULL,
                                                     top_k = 30L,
                                                     alpha = 0.05,
                                                     comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  available <- comparison_available_te_delta_contrasts(comparison_results)
  if (length(available) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared contrast columns with both `delta_*` and `padj_*` are available in te_contrast_summary.tsv.",
      available_contrasts = available
    ))
  }

  if (is.null(contrast_ids) || length(contrast_ids) == 0L) {
    contrast_ids <- available
  } else {
    contrast_ids <- unique(vapply(contrast_ids, comparison_clean_contrast_id, character(1)))
    contrast_ids <- contrast_ids[!is.na(contrast_ids) & nzchar(contrast_ids)]
    missing_contrasts <- setdiff(contrast_ids, available)
    contrast_ids <- intersect(contrast_ids, available)
    if (length(contrast_ids) == 0L) {
      return(list(
        status = "skipped",
        reason = sprintf(
          "None of the requested fixed contrasts are shared with both `delta_*` and `padj_*` columns. Available contrasts: %s",
          paste(available, collapse = ", ")
        ),
        available_contrasts = available,
        missing_contrasts = missing_contrasts
      ))
    }
  }

  summaries <- lapply(contrast_ids, function(contrast_id) {
    summarise_fixed_contrast_te_evidence(
      comparison_results = comparison_results,
      contrast_id = contrast_id,
      selected_te_space = selected_te_space,
      target_tes = target_tes,
      top_k = top_k,
      alpha = alpha,
      comparison_class = comparison_class
    )
  })
  names(summaries) <- contrast_ids

  ok <- vapply(summaries, function(x) identical(x$status, "ok"), logical(1))
  if (!any(ok)) {
    reasons <- vapply(summaries, function(x) {
      reason <- x$reason
      if (is.null(reason) || !nzchar(reason)) "unknown reason" else reason
    }, character(1))
    return(list(
      status = "skipped",
      reason = paste(sprintf("%s: %s", names(reasons), reasons), collapse = "; "),
      available_contrasts = available,
      summaries = summaries
    ))
  }

  list(
    status = "ok",
    contrasts = names(summaries)[ok],
    available_contrasts = available,
    summaries = summaries[ok],
    skipped = summaries[!ok],
    top_k = top_k,
    alpha = alpha,
    comparison_class = comparison_class
  )
}

print_fixed_contrast_te_evidence_plots <- function(fixed_contrast_summaries) {
  comparison_need_ggplot()
  if (!identical(fixed_contrast_summaries$status, "ok")) {
    comparison_warn(fixed_contrast_summaries$reason)
    return(invisible(fixed_contrast_summaries))
  }

  plots <- lapply(fixed_contrast_summaries$summaries, plot_fixed_contrast_te_evidence)
  plots <- plots[vapply(plots, function(plot) inherits(plot, c("ggplot", "patchwork")), logical(1))]
  if (length(plots) == 0L) {
    comparison_warn("No fixed-contrast TE evidence plots could be built.")
    return(invisible(fixed_contrast_summaries))
  }

  for (plot in plots) {
    print(plot)
  }
  invisible(plots)
}

comparison_wide_result_value <- function(result_summary, te_order, value_col) {
  result_levels <- unique(result_summary$result_label[order(match(result_summary$result_id, unique(result_summary$result_id)))])
  out <- data.frame(te_subfamily = te_order, stringsAsFactors = FALSE)
  for (result_label in result_levels) {
    x <- result_summary[result_summary$result_label == result_label, , drop = FALSE]
    key <- paste(x$te_subfamily, x$result_label, sep = "\r")
    lookup <- x[[value_col]]
    names(lookup) <- key
    out[[result_label]] <- unname(lookup[paste(te_order, result_label, sep = "\r")])
  }
  out
}

fixed_contrast_class_table <- function(fixed_contrast_summary) {
  if (!identical(fixed_contrast_summary$status, "ok")) {
    return(data.frame())
  }
  te_order <- fixed_contrast_summary$selected_te_summary$te_subfamily
  class_wide <- comparison_wide_result_value(fixed_contrast_summary$result_summary_top, te_order, "contrast_class")
  te_summary <- fixed_contrast_summary$selected_te_summary
  out <- merge(
    te_summary[, c("te_subfamily", "te_marker_status", "majority_contrast_class", "majority_class_votes", "n_non_noise_votes", "n_results"), drop = FALSE],
    class_wide,
    by = "te_subfamily",
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[match(te_order, out$te_subfamily), , drop = FALSE]
  vote_label <- paste0(out$n_non_noise_votes, "/", out$n_results)
  out <- cbind(
    data.frame(
      `TE subfamily` = out$te_subfamily,
      `TE marker` = out$te_marker_status,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    out[, setdiff(colnames(class_wide), "te_subfamily"), drop = FALSE],
    data.frame(
      `majority class` = out$majority_contrast_class,
      `class votes` = out$majority_class_votes,
      `non-noise dirs` = vote_label,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

fixed_contrast_delta_table <- function(fixed_contrast_summary) {
  if (!identical(fixed_contrast_summary$status, "ok")) {
    return(data.frame())
  }
  te_order <- fixed_contrast_summary$selected_te_summary$te_subfamily
  delta_wide <- comparison_wide_result_value(fixed_contrast_summary$result_summary_top, te_order, "delta")
  result_cols <- setdiff(colnames(delta_wide), "te_subfamily")
  delta_wide[result_cols] <- lapply(delta_wide[result_cols], comparison_format_compact_number)
  te_summary <- fixed_contrast_summary$selected_te_summary
  out <- merge(
    te_summary[, c("te_subfamily", "mean_delta", "min_delta", "max_delta", "n_significant_results", "n_results"), drop = FALSE],
    delta_wide,
    by = "te_subfamily",
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[match(te_order, out$te_subfamily), , drop = FALSE]
  out <- cbind(
    data.frame(`TE subfamily` = out$te_subfamily, check.names = FALSE, stringsAsFactors = FALSE),
    out[, result_cols, drop = FALSE],
    data.frame(
      `mean signed distance` = comparison_format_compact_number(out$mean_delta),
      `min signed distance` = comparison_format_compact_number(out$min_delta),
      `max signed distance` = comparison_format_compact_number(out$max_delta),
      `significant dirs` = paste0(out$n_significant_results, "/", out$n_results),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

fixed_contrast_pvalue_table <- function(fixed_contrast_summary) {
  if (!identical(fixed_contrast_summary$status, "ok")) {
    return(data.frame())
  }
  te_order <- fixed_contrast_summary$selected_te_summary$te_subfamily
  p_wide <- comparison_wide_result_value(fixed_contrast_summary$result_summary_top, te_order, "padj")
  result_cols <- setdiff(colnames(p_wide), "te_subfamily")
  p_wide[result_cols] <- lapply(p_wide[result_cols], comparison_format_p_value)
  te_summary <- fixed_contrast_summary$selected_te_summary
  out <- merge(
    te_summary[, c("te_subfamily", "aggregated_padj", "mean_padj", "min_padj", "n_result_padj"), drop = FALSE],
    p_wide,
    by = "te_subfamily",
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[match(te_order, out$te_subfamily), , drop = FALSE]
  out <- cbind(
    data.frame(`TE subfamily` = out$te_subfamily, check.names = FALSE, stringsAsFactors = FALSE),
    out[, result_cols, drop = FALSE],
    data.frame(
      `Fisher aggregated padj` = comparison_format_p_value(out$aggregated_padj),
      `mean padj` = comparison_format_p_value(out$mean_padj),
      `min padj` = comparison_format_p_value(out$min_padj),
      `p-value dirs` = out$n_result_padj,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

comparison_fixed_contrast_table_html <- function(x, title, subtitle = NULL) {
  if (nrow(x) == 0L) {
    return("")
  }
  sanitize_class <- function(value) {
    value <- tolower(gsub("[^a-z0-9]+", "-", as.character(value)))
    paste0("cell-", gsub("(^-|-$)", "", value))
  }
  class_like <- grepl("class|marker", colnames(x), ignore.case = TRUE)
  rows <- vapply(seq_len(nrow(x)), function(i) {
    cells <- vapply(seq_along(x), function(j) {
      value <- as.character(x[i, j])
      extra_class <- if (isTRUE(class_like[j])) sanitize_class(value) else ""
      sprintf("<td class=\"%s\">%s</td>", extra_class, comparison_escape_html(value))
    }, character(1))
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  }, character(1))
  headers <- paste(sprintf("<th>%s</th>", comparison_escape_html(colnames(x))), collapse = "")
  subtitle_html <- if (is.null(subtitle) || !nzchar(subtitle)) "" else sprintf("<div class=\"fixed-contrast-subtitle\">%s</div>", comparison_escape_html(subtitle))

  paste0(
    "<style>",
    ".fixed-contrast-wrap{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:100%;margin:8px 0 12px 0;}",
    ".fixed-contrast-title{font-weight:700;color:#1f2933;margin:0 0 2px 0;font-size:13px;}",
    ".fixed-contrast-subtitle{color:#667085;font-size:10px;margin:0 0 5px 0;}",
    ".fixed-contrast-table{border-collapse:collapse;width:auto;max-width:100%;font-size:11px;line-height:1.12;}",
    ".fixed-contrast-table th{background:#f4f6f8;color:#1f2933;text-align:left;border-bottom:1px solid #d8dee4;padding:4px 5px;white-space:nowrap;}",
    ".fixed-contrast-table td{border-bottom:1px solid #e5e7eb;padding:3px 5px;vertical-align:top;text-align:left;max-width:210px;overflow-wrap:anywhere;}",
    ".fixed-contrast-table td:first-child{font-weight:600;color:#1f2933;}",
    ".fixed-contrast-table .cell-not-majority-noise{background:#fff1e8;color:#9a3412;font-weight:700;}",
    ".fixed-contrast-table .cell-majority-noise{background:#f2f4f7;color:#475467;font-weight:700;}",
    ".fixed-contrast-table .cell-flat-or-noisy{color:#667085;}",
    ".fixed-contrast-table .cell-mixed{background:#eef4ff;color:#1849a9;font-weight:700;}",
    "</style>",
    "<div class=\"fixed-contrast-wrap\">",
    "<div class=\"fixed-contrast-title\">", comparison_escape_html(title), "</div>",
    subtitle_html,
    "<table class=\"fixed-contrast-table\"><thead><tr>", headers, "</tr></thead><tbody>",
    paste(rows, collapse = ""),
    "</tbody></table></div>"
  )
}

print_fixed_contrast_class_table <- function(fixed_contrast_summary, html = TRUE) {
  out <- fixed_contrast_class_table(fixed_contrast_summary)
  subtitle <- sprintf(
    "Selected contrast: %s; class is derived from selected delta sign and adjusted p <= %s.",
    fixed_contrast_summary$contrast_id,
    fixed_contrast_summary$alpha
  )
  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_fixed_contrast_table_html(out, "Fixed-contrast assigned classes", subtitle))
  } else {
    print(out, row.names = FALSE)
  }
  invisible(out)
}

print_fixed_contrast_delta_table <- function(fixed_contrast_summary, html = TRUE) {
  out <- fixed_contrast_delta_table(fixed_contrast_summary)
  subtitle <- sprintf(
    "Selected contrast: %s; signed distance is mean activity %s - %s, so positive values point to %s and negative values to %s.",
    fixed_contrast_summary$contrast_id,
    fixed_contrast_summary$condition_a,
    fixed_contrast_summary$condition_b,
    fixed_contrast_summary$condition_a,
    fixed_contrast_summary$condition_b
  )
  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_fixed_contrast_table_html(out, "Fixed-contrast signed distance", subtitle))
  } else {
    print(out, row.names = FALSE)
  }
  invisible(out)
}

print_fixed_contrast_pvalue_table <- function(fixed_contrast_summary, html = TRUE) {
  out <- fixed_contrast_pvalue_table(fixed_contrast_summary)
  subtitle <- sprintf(
    "Selected contrast: %s; aggregated p-value uses Fisher's method over result-directory adjusted p-values.",
    fixed_contrast_summary$contrast_id
  )
  if (isTRUE(html) && comparison_notebook_context() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_html(comparison_fixed_contrast_table_html(out, "Fixed-contrast adjusted p-values", subtitle))
  } else {
    print(out, row.names = FALSE)
  }
  invisible(out)
}

summarise_activity_jaccard <- function(comparison_results, target_tes = NULL, threshold = 0, level = c("auto", "cell", "te_subfamily"), comparison_class = NULL) {
  level <- match.arg(level)
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  summary <- comparison_activity_matrices_shared_tes(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    return(summary)
  }

  actual_level <- level
  if (identical(level, "auto")) {
    actual_level <- if (length(summary$shared_samples) > 0L) "cell" else "te_subfamily"
  } else if (identical(level, "cell") && length(summary$shared_samples) == 0L) {
    comparison_warn("No shared activity-score columns are available for TE-sample Jaccard; using TE-subfamily Jaccard instead.")
    actual_level <- "te_subfamily"
  }

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  active_sets <- lapply(summary$matrices, function(mat) {
    if (actual_level == "cell") {
      mat <- mat[, summary$shared_samples, drop = FALSE]
      active_idx <- which(abs(mat) > threshold, arr.ind = TRUE)
      if (nrow(active_idx) == 0L) {
        return(character())
      }
      unique(paste(rownames(mat)[active_idx[, "row"]], colnames(mat)[active_idx[, "col"]], sep = "||"))
    } else {
      rownames(mat)[rowSums(abs(mat) > threshold, na.rm = TRUE) > 0L]
    }
  })
  active_sets <- active_sets[comparison_result_ids(comparison_results)]
  names(active_sets) <- unname(label_map[names(active_sets)])

  labels <- names(active_sets)
  out <- do.call(rbind, lapply(labels, function(a) {
    do.call(rbind, lapply(labels, function(b) {
      intersection_n <- length(intersect(active_sets[[a]], active_sets[[b]]))
      union_n <- length(union(active_sets[[a]], active_sets[[b]]))
      data.frame(
        result_a = a,
        result_b = b,
        intersection_n = intersection_n,
        union_n = union_n,
        jaccard = if (union_n == 0L) NA_real_ else intersection_n / union_n,
        active_a_n = length(active_sets[[a]]),
        active_b_n = length(active_sets[[b]]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(out) <- NULL

  list(
    status = "ok",
    requested_level = level,
    level = actual_level,
    threshold = threshold,
    jaccard = out,
    active_sets = active_sets,
    summary = summary
  )
}

plot_activity_jaccard <- function(comparison_results, target_tes = NULL, threshold = 0, level = c("auto", "cell", "te_subfamily"), comparison_class = NULL) {
  comparison_need_ggplot()
  level <- match.arg(level)
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  summary <- summarise_activity_jaccard(
    comparison_results,
    target_tes = target_tes,
    threshold = threshold,
    level = level,
    comparison_class = comparison_class
  )
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }

  plot_data <- summary$jaccard
  label_levels <- names(summary$active_sets)
  plot_data$result_a <- factor(plot_data$result_a, levels = label_levels)
  plot_data$result_b <- factor(plot_data$result_b, levels = rev(label_levels))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = result_a, y = result_b, fill = jaccard)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(jaccard), "NA", sprintf("%.2f", jaccard))), size = 3.3) +
    ggplot2::scale_fill_gradient(low = "#f8fbff", high = "#2166ac", limits = c(0, 1), na.value = "#f2f4f7") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Active activity-score Jaccard overlap",
      subtitle = comparison_plot_subtitle(
        sprintf("Jaccard overlap of active %s sets; active means abs(activity score) > %s.", ifelse(summary$level == "cell", "TE-sample", "TE-subfamily"), threshold),
        comparison_class
      ),
      x = NULL,
      y = NULL,
      fill = "Jaccard"
    ) +
    comparison_theme()
}

plot_sparsity_comparison <- function(comparison_results, comparison_class = NULL) {
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  fit <- summarise_activity_fit(comparison_results)
  if (!"mean_sparsity_percent" %in% colnames(fit)) {
    comparison_warn("Model-fit sparsity column `mean_sparsity_percent` is not available for plotting.")
    return(invisible(NULL))
  }

  plot_activity_fit_metric(
    comparison_results,
    metric = "mean_sparsity_percent",
    comparison_class = comparison_class
  ) +
    ggplot2::labs(
      title = "Mean sparsity model fit",
      x = "Mean sparsity (%)"
    )
}

plot_lambda_best_comparison <- function(comparison_results, comparison_class = NULL) {
  comparison_need_ggplot()
  summary <- summarise_lambda_best(comparison_results)
  if (!any(summary$lambda_best_available)) {
    comparison_warn("No lambda_best column is available in the loaded comparison TSVs; lambda plot skipped.")
    return(invisible(NULL))
  }

  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  result_ids <- comparison_result_ids(comparison_results)
  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)

  sample_data <- do.call(rbind, lapply(comparison_results, function(result) {
    sample_diag <- result$tables$sample_activity_diagnostics
    if (is.null(sample_diag) || !"lambda_best" %in% colnames(sample_diag)) {
      return(NULL)
    }
    data.frame(
      result_id = result$result_id,
      sample_id = as.character(sample_diag$sample_id),
      condition = if ("condition" %in% colnames(sample_diag)) as.character(sample_diag$condition) else NA_character_,
      lambda_best = as.numeric(sample_diag$lambda_best),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(sample_data) || nrow(sample_data) == 0L) {
    comparison_warn("No sample-specific lambda_best values are available for plotting.")
    return(invisible(NULL))
  }
  sample_data <- sample_data[!is.na(sample_data$lambda_best), , drop = FALSE]
  if (nrow(sample_data) == 0L) {
    comparison_warn("All loaded sample-specific lambda_best values are NA; lambda plot skipped.")
    return(invisible(NULL))
  }
  sample_data$result_label <- unname(label_map[sample_data$result_id])
  sample_data$result_label[is.na(sample_data$result_label) | !nzchar(sample_data$result_label)] <- comparison_short_result_label(sample_data$result_id[is.na(sample_data$result_label) | !nzchar(sample_data$result_label)])
  result_labels <- unname(label_map[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- comparison_short_result_label(result_ids[is.na(result_labels) | !nzchar(result_labels)])
  sample_data$result_label <- factor(sample_data$result_label, levels = rev(result_labels))

  global_data <- summary[summary$global_lambda_available, c("result_id", "global_lambda_min"), drop = FALSE]
  if (nrow(global_data) > 0L) {
    global_data$result_label <- unname(label_map[global_data$result_id])
    global_data$result_label[is.na(global_data$result_label) | !nzchar(global_data$result_label)] <- comparison_short_result_label(global_data$result_id[is.na(global_data$result_label) | !nzchar(global_data$result_label)])
    global_data$result_label <- factor(global_data$result_label, levels = rev(result_labels))
  }

  all_lambdas <- c(sample_data$lambda_best, global_data$global_lambda_min)
  use_log_scale <- all(all_lambdas[!is.na(all_lambdas)] > 0)

  p <- ggplot2::ggplot(sample_data, ggplot2::aes(x = result_label, y = lambda_best)) +
    ggplot2::geom_boxplot(width = 0.48, fill = "#d6e4f0", color = "#344054", outlier.shape = NA, alpha = 0.95) +
    ggplot2::geom_point(
      ggplot2::aes(color = condition),
      position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1),
      size = 2.2,
      alpha = 0.78
    )

  if (nrow(global_data) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data = global_data,
        ggplot2::aes(x = result_label, y = global_lambda_min),
        inherit.aes = FALSE,
        shape = 23,
        size = 3.4,
        stroke = 0.9,
        fill = "#e45756",
        color = "#7a1f1f"
      )
  }

  if (use_log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Lambda comparison",
      subtitle = if (nrow(global_data) > 0L) {
        "Boxes show per-sample lambda_best; red diamonds show the global lambda selected across samples"
      } else {
        "Boxes show per-sample lambda_best; global lambda diagnostics are not available in these result dirs"
      },
      x = NULL,
      y = if (use_log_scale) "Lambda (log10 scale)" else "Lambda",
      color = "Condition"
    ) +
    comparison_theme() +
    ggplot2::theme(legend.position = "bottom")
}

plot_te_padj_evidence_change <- function(comparison_results,
                                         selected_te_space = NULL,
                                         target_tes = NULL,
                                         top_k = 30L,
                                         p_score_mode = c("mean_padj", "mean_neg_log10_padj"),
                                         alpha = 0.05,
                                         comparison_class = NULL) {
  comparison_need_ggplot()
  p_score_mode <- match.arg(p_score_mode)
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  if (!is.null(selected_te_space) && identical(selected_te_space$status, "ok")) {
    target_tes <- selected_te_space$te_subfamilies
    p_score_mode <- selected_te_space$p_score_mode
    alpha <- selected_te_space$alpha
  }

  evidence <- summarise_te_padj_evidence(
    comparison_results = comparison_results,
    target_tes = target_tes,
    p_score_mode = p_score_mode,
    alpha = alpha,
    comparison_class = comparison_class
  )
  if (!identical(evidence$status, "ok")) {
    comparison_warn(evidence$reason)
    return(invisible(evidence))
  }

  te_stats <- evidence$te_summary
  te_stats <- te_stats[!is.na(te_stats$mean_p_score), , drop = FALSE]
  if (nrow(te_stats) == 0L) {
    comparison_warn("No finite adjusted p-value evidence scores are available after filtering.")
    return(invisible(evidence))
  }
  keep <- head(te_stats$te_subfamily, top_k)
  te_stats <- te_stats[te_stats$te_subfamily %in% keep, , drop = FALSE]
  plot_data <- evidence$result_summary[evidence$result_summary$te_subfamily %in% keep, , drop = FALSE]

  levels_te <- rev(te_stats$te_subfamily)
  te_stats$te_subfamily <- factor(te_stats$te_subfamily, levels = levels_te)
  plot_data$te_subfamily <- factor(plot_data$te_subfamily, levels = levels_te)
  result_levels <- unique(unname(comparison_activity_result_label_map(comparison_results, comparison_class)))
  plot_data$result_label <- factor(plot_data$result_label, levels = result_levels)

  x_max <- max(c(plot_data$p_score, te_stats$max_result_p_score), na.rm = TRUE)
  x_min <- min(c(plot_data$p_score, te_stats$min_result_p_score), na.rm = TRUE)
  x_span <- max(x_max - min(0, x_min), 0.01)
  pad <- max(x_span * 0.04, 0.05)
  delta_x <- x_max + pad * 3.2
  te_stats$delta_x <- delta_x
  te_stats$range_label <- sprintf("%.3g", te_stats$p_score_range)
  color_values <- c(
    `delta max` = "#98a2b3",
    stats::setNames(grDevices::hcl.colors(length(result_levels), palette = "Dark 3"), result_levels)
  )

  ggplot2::ggplot() +
    comparison_zero_reference_line() +
    ggplot2::geom_segment(
      data = te_stats,
      ggplot2::aes(x = min_result_p_score, xend = max_result_p_score, y = te_subfamily, yend = te_subfamily, color = "delta max"),
      linewidth = 1.15,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = p_score, y = te_subfamily, color = result_label),
      size = 2.4,
      alpha = 0.9
    ) +
    ggplot2::geom_text(
      data = te_stats,
      ggplot2::aes(x = delta_x, y = te_subfamily, label = range_label),
      hjust = 0,
      color = "#667085",
      size = 2.8
    ) +
    ggplot2::scale_x_continuous(limits = c(min(0, x_min) - pad, delta_x + pad * 3)) +
    ggplot2::scale_color_manual(values = color_values, breaks = names(color_values), drop = FALSE) +
    ggplot2::labs(
      title = "TE adjusted p-value evidence and change",
      subtitle = comparison_plot_subtitle(
        sprintf(
          "Top %s TE subfamilies ranked by %s aggregated across pairwise contrast columns; %s values indicate stronger evidence.",
          length(keep),
          evidence$score_label,
          ifelse(isTRUE(evidence$lower_is_stronger), "smaller", "larger")
        ),
        comparison_class
      ),
      x = paste0(comparison_pretty_label(evidence$score_label), " across pairwise contrasts"),
      y = "TE subfamily",
      color = NULL
    ) +
    comparison_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom"
    )
}

plot_te_classification_similarity <- function(comparison_results,
                                              selected_te_space = NULL,
                                              target_tes = NULL,
                                              top_k = 30L,
                                              p_score_mode = c("mean_padj", "mean_neg_log10_padj"),
                                              alpha = 0.05,
                                              comparison_class = NULL) {
  comparison_need_ggplot()
  p_score_mode <- match.arg(p_score_mode)
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }
  if (!is.null(selected_te_space) && identical(selected_te_space$status, "ok")) {
    target_tes <- selected_te_space$te_subfamilies
    p_score_mode <- selected_te_space$p_score_mode
    alpha <- selected_te_space$alpha
  }

  evidence <- summarise_te_padj_evidence(
    comparison_results = comparison_results,
    target_tes = target_tes,
    p_score_mode = p_score_mode,
    alpha = alpha,
    comparison_class = comparison_class
  )
  if (!identical(evidence$status, "ok")) {
    comparison_warn(evidence$reason)
    return(invisible(evidence))
  }
  if (!"trend_class" %in% colnames(evidence$result_summary) || all(is.na(evidence$result_summary$trend_class))) {
    comparison_warn("No shared `trend_class` column is available for classification similarity.")
    return(invisible(evidence))
  }

  te_stats <- evidence$te_summary
  keep <- head(te_stats$te_subfamily, top_k)
  plot_data <- evidence$result_summary[evidence$result_summary$te_subfamily %in% keep, , drop = FALSE]
  plot_data$trend_class[is.na(plot_data$trend_class) | !nzchar(plot_data$trend_class)] <- "unclassified"

  te_stats <- te_stats[te_stats$te_subfamily %in% keep, , drop = FALSE]
  levels_te <- rev(te_stats$te_subfamily)
  result_levels <- unique(unname(comparison_activity_result_label_map(comparison_results, comparison_class)))
  global_label <- "Global\n(all dirs)"
  result_levels <- c(result_levels, global_label)

  empty_class <- function(x) is.na(x) | !nzchar(x) | x %in% c("flat_or_noisy", "unclassified")
  format_class_evidence <- function(p, trend_class) {
    p <- suppressWarnings(as.numeric(p))
    out <- rep("", length(p))
    missing_supported_class <- is.na(p) & !empty_class(trend_class)
    out[missing_supported_class] <- "NA"
    show <- !is.na(p) & p <= alpha & !empty_class(trend_class)
    out[show] <- comparison_format_p_value(p[show])
    out
  }

  plot_data$cell_type <- "result_dir"
  plot_data$cell_label <- format_class_evidence(plot_data$combined_raw_p_fisher_bh_evidence, plot_data$trend_class)
  plot_data <- plot_data[, c("result_label", "te_subfamily", "trend_class", "combined_raw_p_fisher_bh_evidence", "cell_type", "cell_label"), drop = FALSE]

  global_data <- data.frame(
    result_label = global_label,
    te_subfamily = te_stats$te_subfamily,
    trend_class = if ("global_assigned_class" %in% colnames(te_stats)) te_stats$global_assigned_class else te_stats$majority_class,
    combined_raw_p_fisher_bh_evidence = te_stats$global_combined_raw_p_fisher_bh_evidence,
    cell_type = "global",
    stringsAsFactors = FALSE
  )
  global_data$trend_class[is.na(global_data$trend_class) | !nzchar(global_data$trend_class)] <- "unclassified"
  global_data$cell_label <- format_class_evidence(
    global_data$combined_raw_p_fisher_bh_evidence,
    global_data$trend_class
  )
  plot_data <- rbind(plot_data, global_data[, colnames(plot_data), drop = FALSE])

  plot_data$fill_class <- plot_data$trend_class
  plot_data$fill_class[empty_class(plot_data$fill_class)] <- NA_character_
  plot_data$te_subfamily <- factor(plot_data$te_subfamily, levels = levels_te)
  plot_data$result_label <- factor(plot_data$result_label, levels = result_levels)
  result_dir_text <- plot_data[plot_data$cell_type == "result_dir", , drop = FALSE]
  global_text <- plot_data[plot_data$cell_type == "global", , drop = FALSE]

  class_colors <- comparison_trend_class_colors(unique(na.omit(plot_data$fill_class)))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = result_label, y = te_subfamily, fill = fill_class)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = length(result_levels) - 0.5, color = "#111827", linewidth = 0.75, linetype = "dashed") +
    ggplot2::geom_tile(
      data = global_text,
      ggplot2::aes(x = result_label, y = te_subfamily),
      fill = NA,
      color = "#111827",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = result_dir_text,
      ggplot2::aes(label = cell_label),
      color = "#1f2933",
      size = 2.35,
      lineheight = 0.92
    ) +
    ggplot2::geom_text(
      data = global_text,
      ggplot2::aes(label = cell_label),
      color = "#111827",
      size = 2.45,
      lineheight = 0.92,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = class_colors, na.value = "white", drop = FALSE, na.translate = FALSE) +
    ggplot2::labs(
      title = "TE trend-class similarity",
      subtitle = comparison_plot_subtitle(
        sprintf(
          "Top %s by %s; text = Fisher+BH p; Global = same-class.",
          length(keep),
          evidence$score_label
        ),
        comparison_class
      ),
      x = NULL,
      y = "TE subfamily",
      fill = "Trend class"
    ) +
    comparison_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

plot_te_contrast_comparison <- function(comparison_results, target_tes = NULL, max_tes = 25L) {
  comparison_need_ggplot()
  summary <- summarise_te_contrasts(comparison_results, target_tes = target_tes)
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }
  if (length(summary$effect_columns) == 0L) {
    comparison_warn("No shared effect-size columns are available in te_contrast_summary.tsv.")
    return(invisible(summary))
  }

  plot_cols <- summary$effect_columns
  plot_data <- do.call(rbind, lapply(plot_cols, function(col) {
    data.frame(
      result_id = summary$combined$result_id,
      te_subfamily = summary$combined$te_subfamily,
      contrast = col,
      effect_size = summary$combined[[col]],
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(target_tes) && length(unique(plot_data$te_subfamily)) > max_tes) {
    rank_data <- aggregate(abs(effect_size) ~ te_subfamily, data = plot_data, FUN = max)
    colnames(rank_data)[2] <- "max_abs_effect"
    keep <- head(rank_data$te_subfamily[order(rank_data$max_abs_effect, decreasing = TRUE)], max_tes)
    plot_data <- plot_data[plot_data$te_subfamily %in% keep, , drop = FALSE]
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = te_subfamily, y = effect_size, fill = result_id)) +
    ggplot2::geom_col(position = "dodge", width = 0.72) +
    ggplot2::facet_wrap(~ contrast, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "TE contrast comparison",
      subtitle = "Loaded only from te_contrast_summary.tsv",
      x = "TE subfamily",
      y = "Effect size",
      fill = "Result"
    ) +
    comparison_theme()
}

plot_target_gene_matrix_heatmap <- function(comparison_results,
                                            target_tes = NULL,
                                            target_genes = NULL,
                                            max_tes = 20L,
                                            max_genes = 30L,
                                            target_summary = NULL) {
  comparison_need_ggplot()
  summary <- if (is.null(target_summary)) {
    summarise_target_gene_matrices(comparison_results, target_tes = target_tes, target_genes = target_genes)
  } else {
    target_summary
  }
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }

  plot_data <- summary$heatmap_data
  if (is.null(target_tes) && length(unique(plot_data$te_subfamily)) > max_tes) {
    te_rank <- aggregate(abs(target_score) ~ te_subfamily, data = plot_data, FUN = sum)
    colnames(te_rank)[2] <- "total_abs_score"
    keep_te <- head(te_rank$te_subfamily[order(te_rank$total_abs_score, decreasing = TRUE)], max_tes)
    plot_data <- plot_data[plot_data$te_subfamily %in% keep_te, , drop = FALSE]
  }
  if (is.null(target_genes) && length(unique(plot_data$gene_id)) > max_genes) {
    gene_rank <- aggregate(abs(target_score) ~ gene_id, data = plot_data, FUN = sum)
    colnames(gene_rank)[2] <- "total_abs_score"
    keep_genes <- head(gene_rank$gene_id[order(gene_rank$total_abs_score, decreasing = TRUE)], max_genes)
    plot_data <- plot_data[plot_data$gene_id %in% keep_genes, , drop = FALSE]
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = te_subfamily, y = gene_id, fill = target_score)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(target_matrix ~ result_id, scales = "free", space = "free") +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    ggplot2::labs(
      title = "Target gene matrix heatmap",
      subtitle = "Shared TE subfamilies and shared gene IDs only",
      x = "TE subfamily",
      y = "Gene ID",
      fill = "Target score"
    ) +
    comparison_theme()
}

plot_target_gene_support_overlap <- function(comparison_results,
                                             target_tes = NULL,
                                             target_genes = NULL,
                                             target_summary = NULL) {
  comparison_need_ggplot()
  summary <- if (is.null(target_summary)) {
    summarise_target_gene_matrices(comparison_results, target_tes = target_tes, target_genes = target_genes)
  } else {
    target_summary
  }
  if (!identical(summary$status, "ok")) {
    comparison_warn(summary$reason)
    return(invisible(summary))
  }

  ggplot2::ggplot(summary$support_summary, ggplot2::aes(x = te_subfamily, y = n_supported_in_all_results, fill = target_matrix)) +
    ggplot2::geom_col(position = "dodge", width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Shared target gene support",
      subtitle = "Genes with non-zero target support in all selected result directories",
      x = "TE subfamily",
      y = "Supported genes shared by all results",
      fill = "Target matrix"
    ) +
    comparison_theme()
}

comparison_target_matrix_to_matrix <- function(x, genes, tes) {
  row_idx <- match(genes, x$row_id)
  if (anyNA(row_idx)) {
    stop("Some requested genes are missing from a target gene matrix.", call. = FALSE)
  }

  missing_tes <- setdiff(tes, colnames(x))
  if (length(missing_tes) > 0L) {
    stop(
      sprintf("Some requested TE subfamilies are missing from a target gene matrix: %s", paste(missing_tes, collapse = ", ")),
      call. = FALSE
    )
  }

  mat <- as.matrix(x[row_idx, tes, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat
}

comparison_matrix_col_var <- function(mat) {
  n_non_na <- colSums(!is.na(mat))
  means <- colMeans(mat, na.rm = TRUE)
  centered <- sweep(mat, 2, means, FUN = "-")
  centered[is.na(centered)] <- 0
  out <- colSums(centered * centered) / pmax(n_non_na - 1L, 1L)
  out[n_non_na < 2L] <- NA_real_
  out
}

comparison_target_matrix_contrast_label <- function(matrix_name) {
  comparison_pretty_label(sub("^target_gene_matrix_", "", matrix_name))
}

summarise_target_matrix_similarity <- function(comparison_results,
                                               target_tes = NULL,
                                               target_genes = NULL,
                                               top_k = NULL,
                                               highlight_TEs = character(),
                                               method = c("spearman", "pearson"),
                                               comparison_class = NULL) {
  method <- match.arg(method)
  if (!is.null(top_k)) {
    if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 1L) {
      stop("top_k must be NULL or one positive number.", call. = FALSE)
    }
    top_k <- as.integer(top_k)
  }
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  matrix_names <- get_common_target_matrix_names(comparison_results)
  if (length(matrix_names) == 0L) {
    return(list(
      status = "skipped",
      reason = "No target_gene_matrix_*.tsv file is shared by all selected result directories."
    ))
  }

  result_ids <- comparison_result_ids(comparison_results)
  matrices <- list()
  gene_sets <- list()
  te_sets <- list()

  for (i in seq_along(comparison_results)) {
    result <- comparison_results[[i]]
    result_id <- result_ids[[i]]
    matrices[[result_id]] <- list()

    for (matrix_name in matrix_names) {
      x <- result$target_gene_matrices[[matrix_name]]
      matrices[[result_id]][[matrix_name]] <- x
      profile_key <- paste(result_id, matrix_name, sep = "::")
      gene_sets[[profile_key]] <- unique(x$row_id)
      te_sets[[profile_key]] <- setdiff(colnames(x), "row_id")
    }
  }

  shared_genes <- Reduce(intersect, gene_sets)
  shared_tes <- Reduce(intersect, te_sets)

  if (!is.null(target_genes)) {
    shared_genes <- target_genes[target_genes %in% shared_genes]
  }
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }

  if (length(shared_genes) == 0L || length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared target gene matrix space remains after filtering.",
      common_target_matrices = matrix_names,
      shared_genes = shared_genes,
      shared_te_subfamilies = shared_tes
    ))
  }

  if (!is.null(top_k)) {
    gene_score <- stats::setNames(rep(0, length(shared_genes)), shared_genes)
    te_score <- stats::setNames(rep(0, length(shared_tes)), shared_tes)
    n_panels <- length(result_ids) * length(matrix_names)

    for (result_id in result_ids) {
      for (matrix_name in matrix_names) {
        mat <- comparison_target_matrix_to_matrix(matrices[[result_id]][[matrix_name]], shared_genes, shared_tes)
        panel_gene_score <- rowMeans(abs(mat), na.rm = TRUE)
        panel_te_score <- comparison_matrix_col_var(mat)
        panel_gene_score[!is.finite(panel_gene_score)] <- 0
        panel_te_score[!is.finite(panel_te_score)] <- 0
        gene_score <- gene_score + panel_gene_score
        te_score <- te_score + panel_te_score
      }
    }

    gene_score <- gene_score / n_panels
    te_score <- te_score / n_panels
    shared_genes <- names(sort(gene_score, decreasing = TRUE))[
      seq_len(min(top_k, length(gene_score)))
    ]
    selected_tes <- names(sort(te_score, decreasing = TRUE))[
      seq_len(min(top_k, length(te_score)))
    ]
    forced_tes <- highlight_TEs[highlight_TEs %in% names(te_score) & !highlight_TEs %in% selected_tes]
    shared_tes <- unique(c(forced_tes, selected_tes))
  }

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  profiles <- do.call(rbind, lapply(result_ids, function(result_id) {
    do.call(rbind, lapply(matrix_names, function(matrix_name) {
      result_label <- unname(label_map[result_id])
      if (is.na(result_label) || !nzchar(result_label)) {
        result_label <- comparison_short_result_label(result_id)
      }
      contrast_label <- comparison_target_matrix_contrast_label(matrix_name)
      data.frame(
        result_id = result_id,
        result_label = result_label,
        target_matrix = matrix_name,
        contrast_label = contrast_label,
        profile_id = paste(result_label, contrast_label, sep = " | "),
        stringsAsFactors = FALSE
      )
    }))
  }))
  profiles$profile_index <- seq_len(nrow(profiles))
  profile_ids <- profiles$profile_id

  target_correlation <- matrix(
    NA_real_,
    nrow = length(profile_ids),
    ncol = length(profile_ids),
    dimnames = list(profile_ids, profile_ids)
  )

  profile_vector <- function(profile_row) {
    as.vector(comparison_target_matrix_to_matrix(
      matrices[[profile_row$result_id]][[profile_row$target_matrix]],
      shared_genes,
      shared_tes
    ))
  }

  comparison_type <- function(same_result, same_contrast) {
    if (isTRUE(same_result) && isTRUE(same_contrast)) {
      "same result, same contrast"
    } else if (isTRUE(same_result)) {
      "same result, different contrast"
    } else if (isTRUE(same_contrast)) {
      "different result, same contrast"
    } else {
      "different result, different contrast"
    }
  }

  pair_rows <- vector("list", length(profile_ids) * (length(profile_ids) + 1L) / 2L)
  pair_idx <- 0L
  profile_vectors <- lapply(seq_len(nrow(profiles)), function(i) profile_vector(profiles[i, , drop = FALSE]))

  for (i in seq_along(profile_ids)) {
    vi <- profile_vectors[[i]]

    for (j in seq(from = i, to = length(profile_ids))) {
      vj <- if (i == j) vi else profile_vectors[[j]]
      correlation <- if (i == j) {
        1
      } else {
        suppressWarnings(stats::cor(vi, vj, use = "pairwise.complete.obs", method = method))
      }
      target_correlation[i, j] <- correlation
      target_correlation[j, i] <- correlation

      same_result <- identical(profiles$result_id[[i]], profiles$result_id[[j]])
      same_contrast <- identical(profiles$target_matrix[[i]], profiles$target_matrix[[j]])
      pair_idx <- pair_idx + 1L
      pair_rows[[pair_idx]] <- data.frame(
        profile_a = profiles$profile_id[[i]],
        result_a = profiles$result_id[[i]],
        result_a_label = profiles$result_label[[i]],
        target_matrix_a = profiles$target_matrix[[i]],
        contrast_a = profiles$contrast_label[[i]],
        profile_b = profiles$profile_id[[j]],
        result_b = profiles$result_id[[j]],
        result_b_label = profiles$result_label[[j]],
        target_matrix_b = profiles$target_matrix[[j]],
        contrast_b = profiles$contrast_label[[j]],
        comparison_type = comparison_type(same_result, same_contrast),
        correlation = correlation,
        method = method,
        n_genes = length(shared_genes),
        n_tes = length(shared_tes),
        n_values = length(shared_genes) * length(shared_tes),
        stringsAsFactors = FALSE
      )
    }
  }

  pairwise_similarity <- do.call(rbind, pair_rows)
  same_contrast_cross_result <- pairwise_similarity[
    pairwise_similarity$comparison_type == "different result, same contrast",
    ,
    drop = FALSE
  ]
  if (nrow(same_contrast_cross_result) > 0L) {
    same_contrast_cross_result <- same_contrast_cross_result[
      order(
        same_contrast_cross_result$contrast_a,
        same_contrast_cross_result$result_a_label,
        same_contrast_cross_result$result_b_label
      ),
      ,
      drop = FALSE
    ]
  }

  list(
    status = "ok",
    metric = "correlation",
    method = method,
    top_k = top_k,
    comparison_class = comparison_class,
    common_target_matrices = matrix_names,
    shared_genes = shared_genes,
    shared_te_subfamilies = shared_tes,
    profiles = profiles,
    target_correlation = target_correlation,
    pairwise_similarity = pairwise_similarity,
    same_contrast_cross_result = same_contrast_cross_result
  )
}

summarise_target_matrix_jaccard <- function(comparison_results,
                                            target_tes = NULL,
                                            target_genes = NULL,
                                            top_k = NULL,
                                            highlight_TEs = character(),
                                            threshold = 0,
                                            sign_sensitive = TRUE,
                                            include_zero_state = FALSE,
                                            level = c("gene_te_link", "te_subfamily", "gene"),
                                            comparison_class = NULL) {
  level <- match.arg(level)
  if (!is.null(top_k)) {
    if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 1L) {
      stop("top_k must be NULL or one positive number.", call. = FALSE)
    }
    top_k <- as.integer(top_k)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold < 0) {
    stop("threshold must be one non-negative number.", call. = FALSE)
  }
  if (!is.logical(sign_sensitive) || length(sign_sensitive) != 1L || is.na(sign_sensitive)) {
    stop("sign_sensitive must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_zero_state) || length(include_zero_state) != 1L || is.na(include_zero_state)) {
    stop("include_zero_state must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  matrix_names <- get_common_target_matrix_names(comparison_results)
  if (length(matrix_names) == 0L) {
    return(list(
      status = "skipped",
      reason = "No target_gene_matrix_*.tsv file is shared by all selected result directories."
    ))
  }

  result_ids <- comparison_result_ids(comparison_results)
  matrices <- list()
  gene_sets <- list()
  te_sets <- list()

  for (i in seq_along(comparison_results)) {
    result <- comparison_results[[i]]
    result_id <- result_ids[[i]]
    matrices[[result_id]] <- list()

    for (matrix_name in matrix_names) {
      x <- result$target_gene_matrices[[matrix_name]]
      matrices[[result_id]][[matrix_name]] <- x
      profile_key <- paste(result_id, matrix_name, sep = "::")
      gene_sets[[profile_key]] <- unique(x$row_id)
      te_sets[[profile_key]] <- setdiff(colnames(x), "row_id")
    }
  }

  shared_genes <- Reduce(intersect, gene_sets)
  shared_tes <- Reduce(intersect, te_sets)

  if (!is.null(target_genes)) {
    shared_genes <- target_genes[target_genes %in% shared_genes]
  }
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }

  if (length(shared_genes) == 0L || length(shared_tes) == 0L) {
    return(list(
      status = "skipped",
      reason = "No shared target gene matrix space remains after filtering.",
      common_target_matrices = matrix_names,
      shared_genes = shared_genes,
      shared_te_subfamilies = shared_tes
    ))
  }

  if (!is.null(top_k)) {
    gene_score <- stats::setNames(rep(0, length(shared_genes)), shared_genes)
    te_score <- stats::setNames(rep(0, length(shared_tes)), shared_tes)
    n_panels <- length(result_ids) * length(matrix_names)

    for (result_id in result_ids) {
      for (matrix_name in matrix_names) {
        mat <- comparison_target_matrix_to_matrix(matrices[[result_id]][[matrix_name]], shared_genes, shared_tes)
        panel_gene_score <- rowMeans(abs(mat), na.rm = TRUE)
        panel_te_score <- comparison_matrix_col_var(mat)
        panel_gene_score[!is.finite(panel_gene_score)] <- 0
        panel_te_score[!is.finite(panel_te_score)] <- 0
        gene_score <- gene_score + panel_gene_score
        te_score <- te_score + panel_te_score
      }
    }

    gene_score <- gene_score / n_panels
    te_score <- te_score / n_panels
    shared_genes <- names(sort(gene_score, decreasing = TRUE))[
      seq_len(min(top_k, length(gene_score)))
    ]
    selected_tes <- names(sort(te_score, decreasing = TRUE))[
      seq_len(min(top_k, length(te_score)))
    ]
    forced_tes <- highlight_TEs[highlight_TEs %in% names(te_score) & !highlight_TEs %in% selected_tes]
    shared_tes <- unique(c(forced_tes, selected_tes))
  }

  label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  profiles <- do.call(rbind, lapply(result_ids, function(result_id) {
    do.call(rbind, lapply(matrix_names, function(matrix_name) {
      result_label <- unname(label_map[result_id])
      if (is.na(result_label) || !nzchar(result_label)) {
        result_label <- comparison_short_result_label(result_id)
      }
      contrast_label <- comparison_target_matrix_contrast_label(matrix_name)
      data.frame(
        result_id = result_id,
        result_label = result_label,
        target_matrix = matrix_name,
        contrast_label = contrast_label,
        profile_id = paste(result_label, contrast_label, sep = " | "),
        stringsAsFactors = FALSE
      )
    }))
  }))
  profiles$profile_index <- seq_len(nrow(profiles))
  profile_ids <- profiles$profile_id

  feature_n <- switch(
    level,
    gene_te_link = length(shared_genes) * length(shared_tes),
    te_subfamily = length(shared_tes),
    gene = length(shared_genes)
  )
  feature_states <- if (isTRUE(include_zero_state) && isTRUE(sign_sensitive)) {
    c("positive", "zero", "negative")
  } else if (isTRUE(include_zero_state)) {
    c("active", "zero")
  } else if (isTRUE(sign_sensitive)) {
    c("positive", "negative")
  } else {
    "active"
  }
  jaccard_feature_n <- feature_n * length(feature_states)

  profile_active_set <- function(profile_row) {
    mat <- comparison_target_matrix_to_matrix(
      matrices[[profile_row$result_id]][[profile_row$target_matrix]],
      shared_genes,
      shared_tes
    )
    valid_values <- is.finite(mat)
    state_mat <- matrix(NA_character_, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    if (isTRUE(sign_sensitive)) {
      state_mat[valid_values & mat > threshold] <- "positive"
      state_mat[valid_values & mat < -threshold] <- "negative"
    } else {
      state_mat[valid_values & abs(mat) > threshold] <- "active"
    }
    if (isTRUE(include_zero_state)) {
      state_mat[valid_values & abs(mat) <= threshold] <- "zero"
    }

    feature_idx <- which(!is.na(state_mat), arr.ind = TRUE)
    if (nrow(feature_idx) == 0L) {
      return(character())
    }
    link_ids <- switch(
      level,
      gene_te_link = paste(rownames(mat)[feature_idx[, "row"]], colnames(mat)[feature_idx[, "col"]], sep = "||"),
      te_subfamily = colnames(mat)[feature_idx[, "col"]],
      gene = rownames(mat)[feature_idx[, "row"]]
    )
    if (isTRUE(sign_sensitive) || isTRUE(include_zero_state)) {
      link_ids <- paste(link_ids, state_mat[feature_idx], sep = "||")
    }
    unique(link_ids)
  }

  comparison_type <- function(same_result, same_contrast) {
    if (isTRUE(same_result) && isTRUE(same_contrast)) {
      "same result, same contrast"
    } else if (isTRUE(same_result)) {
      "same result, different contrast"
    } else if (isTRUE(same_contrast)) {
      "different result, same contrast"
    } else {
      "different result, different contrast"
    }
  }

  target_jaccard <- matrix(
    NA_real_,
    nrow = length(profile_ids),
    ncol = length(profile_ids),
    dimnames = list(profile_ids, profile_ids)
  )

  pair_rows <- vector("list", length(profile_ids) * (length(profile_ids) + 1L) / 2L)
  pair_idx <- 0L
  profile_sets <- lapply(seq_len(nrow(profiles)), function(i) profile_active_set(profiles[i, , drop = FALSE]))

  for (i in seq_along(profile_ids)) {
    set_i <- profile_sets[[i]]

    for (j in seq(from = i, to = length(profile_ids))) {
      set_j <- if (i == j) set_i else profile_sets[[j]]
      intersection_n <- length(intersect(set_i, set_j))
      union_n <- length(union(set_i, set_j))
      jaccard <- if (union_n == 0L) NA_real_ else intersection_n / union_n
      target_jaccard[i, j] <- jaccard
      target_jaccard[j, i] <- jaccard

      same_result <- identical(profiles$result_id[[i]], profiles$result_id[[j]])
      same_contrast <- identical(profiles$target_matrix[[i]], profiles$target_matrix[[j]])
      pair_idx <- pair_idx + 1L
      pair_rows[[pair_idx]] <- data.frame(
        profile_a = profiles$profile_id[[i]],
        result_a = profiles$result_id[[i]],
        result_a_label = profiles$result_label[[i]],
        target_matrix_a = profiles$target_matrix[[i]],
        contrast_a = profiles$contrast_label[[i]],
        profile_b = profiles$profile_id[[j]],
        result_b = profiles$result_id[[j]],
        result_b_label = profiles$result_label[[j]],
        target_matrix_b = profiles$target_matrix[[j]],
        contrast_b = profiles$contrast_label[[j]],
        comparison_type = comparison_type(same_result, same_contrast),
        intersection_n = intersection_n,
        union_n = union_n,
        jaccard = jaccard,
        active_a_n = length(set_i),
        active_b_n = length(set_j),
        threshold = threshold,
        sign_sensitive = sign_sensitive,
        include_zero_state = include_zero_state,
        jaccard_level = level,
        jaccard_feature_states = paste(feature_states, collapse = ","),
        n_genes = length(shared_genes),
        n_tes = length(shared_tes),
        n_values = feature_n,
        n_jaccard_features = jaccard_feature_n,
        stringsAsFactors = FALSE
      )
    }
  }

  pairwise_similarity <- do.call(rbind, pair_rows)
  same_contrast_cross_result <- pairwise_similarity[
    pairwise_similarity$comparison_type == "different result, same contrast",
    ,
    drop = FALSE
  ]
  if (nrow(same_contrast_cross_result) > 0L) {
    same_contrast_cross_result <- same_contrast_cross_result[
      order(
        same_contrast_cross_result$contrast_a,
        same_contrast_cross_result$result_a_label,
        same_contrast_cross_result$result_b_label
      ),
      ,
      drop = FALSE
    ]
  }

  list(
    status = "ok",
    metric = "jaccard",
    threshold = threshold,
    sign_sensitive = sign_sensitive,
    include_zero_state = include_zero_state,
    jaccard_feature_states = feature_states,
    jaccard_level = level,
    top_k = top_k,
    comparison_class = comparison_class,
    common_target_matrices = matrix_names,
    shared_genes = shared_genes,
    shared_te_subfamilies = shared_tes,
    profiles = profiles,
    target_jaccard = target_jaccard,
    pairwise_similarity = pairwise_similarity,
    same_contrast_cross_result = same_contrast_cross_result
  )
}

format_target_matrix_similarity_table <- function(target_similarity,
                                                  type = c("same_contrast_cross_result", "all_nonself_pairs", "all_pairs"),
                                                  digits = 3L) {
  type <- match.arg(type)
  if (!identical(target_similarity$status, "ok")) {
    return(data.frame())
  }

  x <- switch(
    type,
    same_contrast_cross_result = target_similarity$same_contrast_cross_result,
    all_nonself_pairs = target_similarity$pairwise_similarity[
      target_similarity$pairwise_similarity$profile_a != target_similarity$pairwise_similarity$profile_b,
      ,
      drop = FALSE
    ],
    all_pairs = target_similarity$pairwise_similarity
  )

  if (nrow(x) == 0L) {
    return(data.frame())
  }

  metric <- target_similarity$metric
  if (is.null(metric) || !nzchar(metric)) {
    metric <- "correlation"
  }
  value_col <- if (identical(metric, "jaccard")) "jaccard" else "correlation"
  if (!value_col %in% colnames(x)) {
    stop(sprintf("Pairwise target matrix summary has no `%s` column.", value_col), call. = FALSE)
  }

  if (identical(type, "same_contrast_cross_result")) {
    out <- data.frame(
      contrast = x$contrast_a,
      result_a = x$result_a_label,
      result_b = x$result_b_label,
      value = round(x[[value_col]], digits),
      n_genes = x$n_genes,
      n_tes = x$n_tes,
      n_values = x$n_values,
      stringsAsFactors = FALSE
    )
    names(out)[names(out) == "value"] <- value_col
    if (identical(metric, "jaccard")) {
      out$threshold <- x$threshold
      out$sign_sensitive <- x$sign_sensitive
      out$include_zero_state <- if ("include_zero_state" %in% colnames(x)) x$include_zero_state else FALSE
      out$jaccard_level <- if ("jaccard_level" %in% colnames(x)) x$jaccard_level else NA_character_
      out$n_jaccard_features <- if ("n_jaccard_features" %in% colnames(x)) x$n_jaccard_features else NA_integer_
      out <- out[, c("contrast", "result_a", "result_b", "jaccard", "threshold", "sign_sensitive", "include_zero_state", "jaccard_level", "n_genes", "n_tes", "n_values", "n_jaccard_features"), drop = FALSE]
    } else {
      out$method <- x$method
      out <- out[, c("contrast", "result_a", "result_b", "correlation", "method", "n_genes", "n_tes", "n_values"), drop = FALSE]
    }
  } else {
    out <- data.frame(
      profile_a = x$profile_a,
      profile_b = x$profile_b,
      comparison_type = x$comparison_type,
      value = round(x[[value_col]], digits),
      n_genes = x$n_genes,
      n_tes = x$n_tes,
      n_values = x$n_values,
      stringsAsFactors = FALSE
    )
    names(out)[names(out) == "value"] <- value_col
    if (identical(metric, "jaccard")) {
      out$threshold <- x$threshold
      out$sign_sensitive <- x$sign_sensitive
      out$include_zero_state <- if ("include_zero_state" %in% colnames(x)) x$include_zero_state else FALSE
      out$jaccard_level <- if ("jaccard_level" %in% colnames(x)) x$jaccard_level else NA_character_
      out$n_jaccard_features <- if ("n_jaccard_features" %in% colnames(x)) x$n_jaccard_features else NA_integer_
      out <- out[, c("profile_a", "profile_b", "comparison_type", "jaccard", "threshold", "sign_sensitive", "include_zero_state", "jaccard_level", "n_genes", "n_tes", "n_values", "n_jaccard_features"), drop = FALSE]
    } else {
      out$method <- x$method
      out <- out[, c("profile_a", "profile_b", "comparison_type", "correlation", "method", "n_genes", "n_tes", "n_values"), drop = FALSE]
    }
  }

  rownames(out) <- NULL
  out
}

plot_target_matrix_similarity <- function(target_similarity,
                                          triangle = c("upper", "full"),
                                          show_values = TRUE,
                                          show_comparison_type = TRUE,
                                          text_size = 3.0) {
  comparison_need_ggplot()
  triangle <- match.arg(triangle)
  if (!identical(target_similarity$status, "ok")) {
    comparison_warn(target_similarity$reason)
    return(invisible(target_similarity))
  }

  profiles <- target_similarity$profiles
  profile_levels <- profiles$profile_id
  profile_axis_labels <- stats::setNames(
    paste(profiles$result_label, profiles$contrast_label, sep = "\n"),
    profiles$profile_id
  )
  metric <- target_similarity$metric
  if (is.null(metric) || !nzchar(metric)) {
    metric <- "correlation"
  }
  value_col <- if (identical(metric, "jaccard")) "jaccard" else "correlation"
  pairwise <- target_similarity$pairwise_similarity
  if (!value_col %in% colnames(pairwise)) {
    stop(sprintf("Pairwise target matrix summary has no `%s` column.", value_col), call. = FALSE)
  }
  if (identical(triangle, "full")) {
    mirrored <- pairwise[pairwise$profile_a != pairwise$profile_b, , drop = FALSE]
    if (nrow(mirrored) > 0L) {
      swap_cols <- list(
        profile_a = "profile_b",
        result_a = "result_b",
        result_a_label = "result_b_label",
        target_matrix_a = "target_matrix_b",
        contrast_a = "contrast_b",
        profile_b = "profile_a",
        result_b = "result_a",
        result_b_label = "result_a_label",
        target_matrix_b = "target_matrix_a",
        contrast_b = "contrast_a"
      )
      mirrored_original <- mirrored
      for (to_col in names(swap_cols)) {
        mirrored[[to_col]] <- mirrored_original[[swap_cols[[to_col]]]]
      }
      pairwise <- rbind(pairwise, mirrored)
    }
  }
  pairwise$x_profile <- if (identical(triangle, "upper")) pairwise$profile_b else pairwise$profile_a
  pairwise$y_profile <- if (identical(triangle, "upper")) pairwise$profile_a else pairwise$profile_b

  relation_colors <- c(
    "same dir" = "#00A676",
    "same contrast" = "#FFB000"
  )
  pairwise$plot_value <- pairwise[[value_col]]
  pairwise$tile_label <- ifelse(is.na(pairwise$plot_value), "NA", sprintf("%.2f", pairwise$plot_value))
  pairwise$x_profile <- factor(pairwise$x_profile, levels = profile_levels)
  pairwise$y_profile <- factor(pairwise$y_profile, levels = rev(profile_levels))
  pairwise$x_num <- as.numeric(pairwise$x_profile)
  pairwise$y_num <- as.numeric(pairwise$y_profile)

  empty_border_segments <- data.frame(
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric(),
    orientation = character(),
    fixed = numeric(),
    start = numeric(),
    end = numeric(),
    normal_x = numeric(),
    normal_y = numeric(),
    relation_label = factor(character(), levels = names(relation_colors)),
    stringsAsFactors = FALSE
  )

  merge_collinear_segments <- function(segments) {
    if (nrow(segments) == 0L) {
      return(segments)
    }
    has_shifted <- "shifted" %in% colnames(segments)
    if (!isTRUE(has_shifted)) {
      segments$shifted <- FALSE
    }

    segments$fixed <- ifelse(segments$orientation == "h", segments$y, segments$x)
    segments$start <- ifelse(segments$orientation == "h", pmin(segments$x, segments$xend), pmin(segments$y, segments$yend))
    segments$end <- ifelse(segments$orientation == "h", pmax(segments$x, segments$xend), pmax(segments$y, segments$yend))
    split_key <- paste(
      segments$relation_label,
      segments$orientation,
      sprintf("%.6f", segments$fixed),
      sprintf("%.6f", segments$normal_x),
      sprintf("%.6f", segments$normal_y),
      sep = "::"
    )

    merged <- do.call(rbind, lapply(split(segments, split_key), function(x) {
      x <- x[order(x$start, x$end), , drop = FALSE]
      out <- list()
      current_start <- x$start[[1]]
      current_end <- x$end[[1]]

      if (nrow(x) > 1L) {
        for (i in 2:nrow(x)) {
          if (x$start[[i]] <= current_end + 1e-8) {
            current_end <- max(current_end, x$end[[i]])
          } else {
            out[[length(out) + 1L]] <- data.frame(
              relation_label = x$relation_label[[1]],
              orientation = x$orientation[[1]],
              fixed = x$fixed[[1]],
              start = current_start,
              end = current_end,
              normal_x = x$normal_x[[1]],
              normal_y = x$normal_y[[1]],
              shifted = x$shifted[[1]],
              stringsAsFactors = FALSE
            )
            current_start <- x$start[[i]]
            current_end <- x$end[[i]]
          }
        }
      }

      out[[length(out) + 1L]] <- data.frame(
        relation_label = x$relation_label[[1]],
        orientation = x$orientation[[1]],
        fixed = x$fixed[[1]],
        start = current_start,
        end = current_end,
        normal_x = x$normal_x[[1]],
        normal_y = x$normal_y[[1]],
        shifted = x$shifted[[1]],
        stringsAsFactors = FALSE
      )
      do.call(rbind, out)
    }))

    out <- data.frame(
      x = ifelse(merged$orientation == "h", merged$start, merged$fixed),
      y = ifelse(merged$orientation == "h", merged$fixed, merged$start),
      xend = ifelse(merged$orientation == "h", merged$end, merged$fixed),
      yend = ifelse(merged$orientation == "h", merged$fixed, merged$end),
      orientation = merged$orientation,
      fixed = merged$fixed,
      start = merged$start,
      end = merged$end,
      normal_x = merged$normal_x,
      normal_y = merged$normal_y,
      shifted = merged$shifted,
      relation_label = factor(merged$relation_label, levels = names(relation_colors)),
      stringsAsFactors = FALSE
    )
    if (!isTRUE(has_shifted)) {
      out$shifted <- NULL
    }
    out
  }

  make_border_segments <- function(cells, relation_label) {
    if (nrow(cells) == 0L) {
      return(empty_border_segments)
    }

    xmin <- cells$x_num - 0.5
    xmax <- cells$x_num + 0.5
    ymin <- cells$y_num - 0.5
    ymax <- cells$y_num + 0.5
    segments <- rbind(
      data.frame(orientation = "h", x = xmin, y = ymin, xend = xmax, yend = ymin, normal_x = 0, normal_y = 1, stringsAsFactors = FALSE),
      data.frame(orientation = "h", x = xmin, y = ymax, xend = xmax, yend = ymax, normal_x = 0, normal_y = -1, stringsAsFactors = FALSE),
      data.frame(orientation = "v", x = xmin, y = ymin, xend = xmin, yend = ymax, normal_x = 1, normal_y = 0, stringsAsFactors = FALSE),
      data.frame(orientation = "v", x = xmax, y = ymin, xend = xmax, yend = ymax, normal_x = -1, normal_y = 0, stringsAsFactors = FALSE)
    )
    segments$key <- paste(
      segments$orientation,
      sprintf("%.6f", segments$x),
      sprintf("%.6f", segments$y),
      sprintf("%.6f", segments$xend),
      sprintf("%.6f", segments$yend),
      sep = "::"
    )
    key_counts <- table(segments$key)
    segments <- segments[key_counts[segments$key] == 1L, , drop = FALSE]
    if (nrow(segments) == 0L) {
      return(empty_border_segments)
    }
    segments$relation_label <- relation_label
    merge_collinear_segments(segments)
  }

  offset_overlapping_edges <- function(segments, base_segments, relation_label, offset = 0.04) {
    use_segments <- segments[as.character(segments$relation_label) == relation_label, , drop = FALSE]
    other_segments <- segments[as.character(segments$relation_label) != relation_label, , drop = FALSE]
    if (nrow(use_segments) == 0L || nrow(base_segments) == 0L) {
      return(segments)
    }

    split_segments <- do.call(rbind, lapply(seq_len(nrow(use_segments)), function(i) {
      segment <- use_segments[i, , drop = FALSE]
      overlaps <- base_segments[
        base_segments$orientation == segment$orientation &
          abs(base_segments$fixed - segment$fixed) < 1e-8 &
          pmax(base_segments$start, segment$start) < pmin(base_segments$end, segment$end) - 1e-8,
        ,
        drop = FALSE
      ]
      if (nrow(overlaps) == 0L) {
        segment$overlaps_base <- FALSE
        return(segment)
      }

      breaks <- sort(unique(c(
        segment$start,
        segment$end,
        pmax(overlaps$start, segment$start),
        pmin(overlaps$end, segment$end)
      )))
      pieces <- lapply(seq_len(length(breaks) - 1L), function(j) {
        piece_start <- breaks[[j]]
        piece_end <- breaks[[j + 1L]]
        if (piece_end <= piece_start + 1e-8) {
          return(NULL)
        }
        piece <- segment
        piece$start <- piece_start
        piece$end <- piece_end
        piece$overlaps_base <- any(
          pmax(overlaps$start, piece_start) < pmin(overlaps$end, piece_end) - 1e-8
        )
        piece
      })
      do.call(rbind, Filter(Negate(is.null), pieces))
    }))

    split_segments$fixed <- ifelse(
      split_segments$overlaps_base,
      split_segments$fixed + ifelse(split_segments$orientation == "h", split_segments$normal_y, split_segments$normal_x) * offset,
      split_segments$fixed
    )
    split_segments$shifted <- split_segments$overlaps_base
    split_segments$x <- ifelse(split_segments$orientation == "h", split_segments$start, split_segments$fixed)
    split_segments$y <- ifelse(split_segments$orientation == "h", split_segments$fixed, split_segments$start)
    split_segments$xend <- ifelse(split_segments$orientation == "h", split_segments$end, split_segments$fixed)
    split_segments$yend <- ifelse(split_segments$orientation == "h", split_segments$fixed, split_segments$end)
    split_segments$overlaps_base <- NULL

    shifted_segments <- merge_collinear_segments(split_segments)
    shifted_idx <- isTRUE(shifted_segments$shifted) | shifted_segments$shifted
    if (any(shifted_idx)) {
      segment_length <- shifted_segments$end[shifted_idx] - shifted_segments$start[shifted_idx]
      shrink <- pmin(offset * 1.45, pmax(segment_length / 2 - 1e-6, 0))
      shifted_segments$start[shifted_idx] <- shifted_segments$start[shifted_idx] + shrink
      shifted_segments$end[shifted_idx] <- shifted_segments$end[shifted_idx] - shrink
      shifted_segments$x[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$start[shifted_idx],
        shifted_segments$fixed[shifted_idx]
      )
      shifted_segments$y[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$fixed[shifted_idx],
        shifted_segments$start[shifted_idx]
      )
      shifted_segments$xend[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$end[shifted_idx],
        shifted_segments$fixed[shifted_idx]
      )
      shifted_segments$yend[shifted_idx] <- ifelse(
        shifted_segments$orientation[shifted_idx] == "h",
        shifted_segments$fixed[shifted_idx],
        shifted_segments$end[shifted_idx]
      )
    }
    shifted_segments$shifted <- NULL

    rbind(other_segments, shifted_segments)
  }

  same_dir_cells <- pairwise[pairwise$result_a == pairwise$result_b, , drop = FALSE]
  same_contrast_cells <- pairwise[pairwise$target_matrix_a == pairwise$target_matrix_b, , drop = FALSE]
  same_dir_segments <- make_border_segments(same_dir_cells, "same dir")
  same_contrast_segments <- make_border_segments(same_contrast_cells, "same contrast")
  border_segments <- rbind(
    same_dir_segments,
    same_contrast_segments
  )
  border_segments <- offset_overlapping_edges(border_segments, same_dir_segments, "same contrast")
  border_segments$relation_label <- factor(border_segments$relation_label, levels = names(relation_colors))
  relation_count <- length(relation_colors)
  matrix_legend_height_mm <- max(52, min(76, length(profile_levels) * 14))
  legend_block_gap_mm <- max(1.5, min(4, matrix_legend_height_mm * 0.045))
  relation_gap_mm <- max(1, min(3, matrix_legend_height_mm * 0.025))
  fill_barheight_mm <- max(14, min(20, matrix_legend_height_mm * 0.25))
  relation_keyheight_mm <- max(
    18,
    min(
      28,
      (matrix_legend_height_mm * 0.72 - relation_gap_mm * max(relation_count - 1, 0)) /
        max(relation_count, 1)
    )
  )
  relation_legend_data <- data.frame(
    x = rep(pairwise$x_num[[1]], length(relation_colors)),
    y = rep(pairwise$y_num[[1]], length(relation_colors)),
    relation_label = factor(names(relation_colors), levels = names(relation_colors)),
    stringsAsFactors = FALSE
  )

  scope_label <- if (is.null(target_similarity$top_k)) "shared" else sprintf("top-%d", target_similarity$top_k)
  if (identical(metric, "jaccard")) {
    threshold <- if (is.null(target_similarity$threshold)) 0 else target_similarity$threshold
    threshold_label <- format(threshold, trim = TRUE, scientific = FALSE)
    jaccard_level <- target_similarity$jaccard_level
    if (is.null(jaccard_level) || !nzchar(jaccard_level)) {
      jaccard_level <- "gene_te_link"
    }
    level_title <- switch(
      jaccard_level,
      te_subfamily = "TE subfamilies",
      gene = "genes",
      gene_te_link = "Gene-TE links",
      "Gene-TE links"
    )
    include_zero_state <- isTRUE(target_similarity$include_zero_state)
    level_description <- if (isTRUE(include_zero_state)) {
      switch(
        jaccard_level,
        te_subfamily = sprintf(
          "Jaccard overlap of TE-subfamily state sets over %s matrices (%d TE subfamilies); states are derived from target-gene links with threshold %s",
          scope_label,
          length(target_similarity$shared_te_subfamilies),
          threshold_label
        ),
        gene = sprintf(
          "Jaccard overlap of gene state sets over %s matrices (%d genes); states are derived from TE-subfamily links with threshold %s",
          scope_label,
          length(target_similarity$shared_genes),
          threshold_label
        ),
        gene_te_link = sprintf(
          "Jaccard overlap of Gene-TE state links over %s matrices (%d genes x %d TE subfamilies); each value is assigned to a target-score state with threshold %s",
          scope_label,
          length(target_similarity$shared_genes),
          length(target_similarity$shared_te_subfamilies),
          threshold_label
        ),
        sprintf(
          "Jaccard overlap over %s matrices; each value is assigned to a target-score state with threshold %s",
          scope_label,
          threshold_label
        )
      )
    } else {
      switch(
        jaccard_level,
        te_subfamily = sprintf(
          "Jaccard overlap of active TE-subfamily sets over %s matrices (%d TE subfamilies); a TE is active if at least one target-gene link has abs(target score) > %s",
          scope_label,
          length(target_similarity$shared_te_subfamilies),
          threshold_label
        ),
        gene = sprintf(
          "Jaccard overlap of active gene sets over %s matrices (%d genes); a gene is active if at least one TE-subfamily link has abs(target score) > %s",
          scope_label,
          length(target_similarity$shared_genes),
          threshold_label
        ),
        gene_te_link = sprintf(
          "Jaccard overlap of active Gene-TE links over %s matrices (%d genes x %d TE subfamilies); active means abs(target score) > %s",
          scope_label,
          length(target_similarity$shared_genes),
          length(target_similarity$shared_te_subfamilies),
          threshold_label
        ),
        sprintf(
          "Jaccard overlap over %s matrices; active means abs(target score) > %s",
          scope_label,
          threshold_label
        )
      )
    }
    sign_label <- if (isTRUE(include_zero_state) && isTRUE(target_similarity$sign_sensitive)) {
      "positive, zero, and negative contributions are treated as distinct state features"
    } else if (isTRUE(include_zero_state)) {
      "zero and non-zero contributions are treated as distinct state features"
    } else if (isTRUE(target_similarity$sign_sensitive)) {
      "positive and negative contributions are treated as different signed features"
    } else {
      "positive and negative contributions are both counted as active features"
    }
    title <- sprintf("Target-matrix Jaccard overlap: %s", level_title)
    fill_label <- "Jaccard"
    subtitle <- sprintf("%s; %s.", level_description, sign_label)
    fill_scale <- ggplot2::scale_fill_gradient(
      low = "#f8fbff",
      high = "#2166ac",
      limits = c(0, 1),
      na.value = "#f2f4f7"
    )
  } else {
    method_title <- if (identical(target_similarity$method, "pearson")) "Pearson" else "Spearman"
    title <- "Target-matrix similarity"
    fill_label <- if (identical(target_similarity$method, "pearson")) "Pearson r" else "Spearman rho"
    subtitle <- sprintf(
      "%s correlation over %s target-score vectors (%d genes x %d TE subfamilies).",
      method_title,
      scope_label,
      length(target_similarity$shared_genes),
      length(target_similarity$shared_te_subfamilies)
    )
    fill_scale <- ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-1, 1)
    )
  }
  if (!is.null(target_similarity$comparison_class)) {
    subtitle <- comparison_plot_subtitle(subtitle, target_similarity$comparison_class)
  }

  plot <- ggplot2::ggplot(pairwise, ggplot2::aes(x = x_profile, y = y_profile, fill = plot_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::geom_segment(
      data = border_segments,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        color = relation_label
      ),
      inherit.aes = FALSE,
      linewidth = 1.15,
      lineend = "butt",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = relation_legend_data,
      ggplot2::aes(x = x, y = y, color = relation_label),
      inherit.aes = FALSE,
      shape = 22,
      fill = "white",
      size = 3.2,
      stroke = 1.25,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_x_discrete(labels = profile_axis_labels) +
    ggplot2::scale_y_discrete(labels = profile_axis_labels) +
    fill_scale +
    ggplot2::scale_color_manual(
      values = relation_colors,
      drop = FALSE
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      fill = fill_label,
      color = "Relation"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        order = 1,
        title.position = "top",
        barheight = grid::unit(fill_barheight_mm, "mm"),
        barwidth = grid::unit(4, "mm")
      ),
      color = ggplot2::guide_legend(
        order = 2,
        title.position = "top",
        direction = "vertical",
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(relation_keyheight_mm, "mm"),
        keywidth = grid::unit(4.5, "mm"),
        override.aes = list(fill = "white", alpha = 1, shape = 22, size = 3.2, stroke = 1.25),
        theme = ggplot2::theme(
          legend.text = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8.5),
          legend.key.width = grid::unit(4.5, "mm"),
          legend.key.height = grid::unit(relation_keyheight_mm, "mm"),
          legend.spacing.y = grid::unit(relation_gap_mm, "mm")
        )
      )
    ) +
    comparison_theme() +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 8.5, lineheight = 0.95),
      legend.key.width = grid::unit(4.5, "mm"),
      legend.key.height = grid::unit(relation_keyheight_mm, "mm"),
      legend.spacing.y = grid::unit(relation_gap_mm, "mm"),
      legend.box.spacing = grid::unit(legend_block_gap_mm, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 4),
      plot.margin = ggplot2::margin(5.5, 4, 5.5, 5.5)
    )

  if (isTRUE(show_values)) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes(label = tile_label), size = text_size, lineheight = 0.9)
  }

  plot
}

plot_target_matrix_jaccard <- function(target_jaccard,
                                       triangle = c("upper", "full"),
                                       show_values = TRUE,
                                       show_comparison_type = TRUE,
                                       text_size = 3.0) {
  plot_target_matrix_similarity(
    target_jaccard,
    triangle = triangle,
    show_values = show_values,
    show_comparison_type = show_comparison_type,
    text_size = text_size
  )
}

comparison_target_gene_label_map <- function(comparison_results) {
  mapping_paths <- unique(unlist(lapply(comparison_results, function(result) {
    mapping_dir <- file.path(result$result_dir, "mapping")
    if (!dir.exists(mapping_dir)) {
      return(character())
    }
    list.files(mapping_dir, pattern = "^map_.*\\.tsv$", full.names = TRUE)
  }), use.names = FALSE))

  if (length(mapping_paths) == 0L) {
    return(character())
  }

  maps <- lapply(mapping_paths, function(path) {
    x <- comparison_read_tsv(path, required = FALSE, table_name = basename(path), warn_missing = FALSE)
    if (is.null(x) || !all(c("ensembl_gene_id", "hgnc_symbol") %in% colnames(x))) {
      return(NULL)
    }
    x <- x[!is.na(x$ensembl_gene_id) & nzchar(x$ensembl_gene_id) &
             !is.na(x$hgnc_symbol) & nzchar(x$hgnc_symbol),
           c("ensembl_gene_id", "hgnc_symbol"), drop = FALSE]
    x[!duplicated(x$ensembl_gene_id), , drop = FALSE]
  })
  maps <- Filter(Negate(is.null), maps)
  if (length(maps) == 0L) {
    return(character())
  }

  merged <- do.call(rbind, maps)
  merged <- merged[!duplicated(merged$ensembl_gene_id), , drop = FALSE]
  stats::setNames(merged$hgnc_symbol, merged$ensembl_gene_id)
}

comparison_row_z_values <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  as.numeric(z)
}

comparison_draw_target_color_legend <- function(limit,
                                                pal = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)) {
  if (!is.finite(limit) || limit <= 0) {
    limit <- 2
  }

  x0 <- grid::unit(0.36, "npc")
  y0 <- grid::unit(0.038, "npc")
  legend_width <- grid::unit(0.28, "npc")
  legend_height <- grid::unit(0.015, "npc")
  n <- length(pal)

  for (i in seq_along(pal)) {
    grid::grid.rect(
      x = x0 + legend_width * ((i - 0.5) / n),
      y = y0,
      width = legend_width / n,
      height = legend_height,
      gp = grid::gpar(fill = pal[i], col = NA)
    )
  }

  tick_x <- x0 + legend_width * c(0, 0.5, 1)
  tick_labels <- c(sprintf("%.1f", -limit), "0", sprintf("%.1f", limit))
  grid::grid.segments(
    x0 = tick_x,
    x1 = tick_x,
    y0 = y0 - grid::unit(0.011, "npc"),
    y1 = y0 + grid::unit(0.011, "npc"),
    gp = grid::gpar(col = "grey35", linewidth = 0.35)
  )
  grid::grid.text(
    tick_labels,
    x = tick_x,
    y = y0 - grid::unit(0.025, "npc"),
    gp = grid::gpar(fontsize = 7, col = "grey35")
  )
  grid::grid.text(
    "row z-score",
    x = x0 - grid::unit(0.02, "npc"),
    y = y0,
    just = "right",
    gp = grid::gpar(fontsize = 8, col = "grey35", fontface = "bold")
  )
}

comparison_draw_target_annotation_legend <- function(selected_tes,
                                                     highlight_TEs = c("LTR5_Hs", "LTR5A", "LTR5B")) {
  fixed_hl_colors <- c(LTR5_Hs = "#E63946", LTR5A = "#F4A261", LTR5B = "#E9C46A")
  highlight_present <- intersect(unique(highlight_TEs), selected_tes)
  missing_hl_colors <- setdiff(highlight_present, names(fixed_hl_colors))
  dynamic_hl_colors <- if (length(missing_hl_colors) > 0L) {
    stats::setNames(grDevices::hcl.colors(length(missing_hl_colors), "Dark 3"), missing_hl_colors)
  } else {
    character(0)
  }
  highlight_cols <- c(
    fixed_hl_colors[names(fixed_hl_colors) %in% highlight_present],
    dynamic_hl_colors,
    other = "#ADB5BD"
  )
  rank_cols <- c("Top 10" = "#2D6A4F", "Top 25" = "#74C69D", "Rest" = "#D8F3DC")

  draw_swatches <- function(title, cols, x, y) {
    grid::grid.text(
      title,
      x = x,
      y = y,
      just = "right",
      gp = grid::gpar(fontsize = 8, col = "grey35", fontface = "bold")
    )
    x_cursor <- x + grid::unit(0.01, "npc")
    for (nm in names(cols)) {
      grid::grid.rect(
        x = x_cursor,
        y = y,
        width = grid::unit(0.012, "npc"),
        height = grid::unit(0.012, "npc"),
        just = "left",
        gp = grid::gpar(fill = cols[[nm]], col = "grey55", linewidth = 0.25)
      )
      grid::grid.text(
        nm,
        x = x_cursor + grid::unit(0.016, "npc"),
        y = y,
        just = "left",
        gp = grid::gpar(fontsize = 7, col = "grey35")
      )
      label_width <- max(0.045, 0.006 * nchar(nm))
      x_cursor <- x_cursor + grid::unit(0.026 + label_width, "npc")
    }
  }

  draw_swatches("Highlight", highlight_cols, grid::unit(0.18, "npc"), grid::unit(0.078, "npc"))
  draw_swatches("Global rank group", rank_cols, grid::unit(0.18, "npc"), grid::unit(0.058, "npc"))
}

plot_target_matrix_full_frame <- function(comparison_results,
                                          top_k = 30L,
                                          target_tes = NULL,
                                          target_genes = NULL,
                                          highlight_TEs = c("LTR5_Hs", "LTR5A", "LTR5B"),
                                          ensembl_to_symbol = NULL,
                                          gene_label = c("symbol", "row_id"),
                                          comparison_class = NULL,
                                          color_limit = NULL,
                                          show_panel_annotations = FALSE,
                                          out_file = NULL,
                                          show = TRUE) {
  if (!exists("plot_target_matrix_full", mode = "function")) {
    stop("plot_target_matrix_full() is not available. Source scripts/craTEs_own_function.R first.", call. = FALSE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package `grid` is required for target matrix frame plots.", call. = FALSE)
  }
  if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 1L) {
    stop("top_k must be one positive number.", call. = FALSE)
  }
  top_k <- as.integer(top_k)
  gene_label <- match.arg(gene_label)
  if (is.null(comparison_class)) {
    comparison_class <- classify_comparison_set(comparison_results)
  }

  matrix_names <- get_common_target_matrix_names(comparison_results)
  if (length(matrix_names) == 0L) {
    reason <- "No target_gene_matrix_*.tsv file is shared by all selected result directories."
    comparison_warn(reason)
    return(invisible(list(status = "skipped", reason = reason)))
  }

  result_ids <- comparison_result_ids(comparison_results)
  matrices <- list()
  gene_sets <- list()
  te_sets <- list()

  for (result_id in result_ids) {
    result <- comparison_results[[result_id]]
    matrices[[result_id]] <- list()
    for (matrix_name in matrix_names) {
      x <- result$target_gene_matrices[[matrix_name]]
      matrices[[result_id]][[matrix_name]] <- x
      gene_sets[[paste(result_id, matrix_name, sep = "::")]] <- unique(x$row_id)
      te_sets[[paste(result_id, matrix_name, sep = "::")]] <- setdiff(colnames(x), "row_id")
    }
  }

  shared_genes <- Reduce(intersect, gene_sets)
  shared_tes <- Reduce(intersect, te_sets)

  if (!is.null(target_genes)) {
    shared_genes <- intersect(shared_genes, target_genes)
  }
  if (!is.null(target_tes)) {
    shared_tes <- target_tes[target_tes %in% shared_tes]
  }

  if (length(shared_genes) == 0L || length(shared_tes) == 0L) {
    reason <- "No shared target gene matrix space remains after filtering."
    comparison_warn(reason)
    return(invisible(list(status = "skipped", reason = reason)))
  }

  gene_score <- stats::setNames(rep(0, length(shared_genes)), shared_genes)
  te_score <- stats::setNames(rep(0, length(shared_tes)), shared_tes)
  n_panels <- length(result_ids) * length(matrix_names)

  for (result_id in result_ids) {
    for (matrix_name in matrix_names) {
      mat <- comparison_target_matrix_to_matrix(matrices[[result_id]][[matrix_name]], shared_genes, shared_tes)
      panel_gene_score <- rowMeans(abs(mat), na.rm = TRUE)
      panel_te_score <- comparison_matrix_col_var(mat)
      panel_gene_score[!is.finite(panel_gene_score)] <- 0
      panel_te_score[!is.finite(panel_te_score)] <- 0
      gene_score <- gene_score + panel_gene_score
      te_score <- te_score + panel_te_score
    }
  }

  gene_score <- gene_score / n_panels
  te_score <- te_score / n_panels
  selected_genes <- names(sort(gene_score, decreasing = TRUE))[
    seq_len(min(top_k, length(gene_score)))
  ]
  selected_tes <- names(sort(te_score, decreasing = TRUE))[
    seq_len(min(top_k, length(te_score)))
  ]
  forced_tes <- highlight_TEs[highlight_TEs %in% shared_tes & !highlight_TEs %in% selected_tes]
  selected_tes <- unique(c(forced_tes, selected_tes))

  if (identical(gene_label, "row_id")) {
    gene_label_map <- stats::setNames(selected_genes, selected_genes)
  } else if (!is.null(ensembl_to_symbol) && length(ensembl_to_symbol) > 0L) {
    gene_label_map <- ensembl_to_symbol
  } else {
    gene_label_map <- comparison_target_gene_label_map(comparison_results)
  }
  if (length(gene_label_map) == 0L) {
    gene_label_map <- stats::setNames(selected_genes, selected_genes)
  }

  selected_panel_mats <- list()
  z_values <- numeric()
  for (result_id in result_ids) {
    selected_panel_mats[[result_id]] <- list()
    for (matrix_name in matrix_names) {
      plot_mat <- comparison_target_matrix_to_matrix(
        matrices[[result_id]][[matrix_name]],
        selected_genes,
        selected_tes
      )
      selected_panel_mats[[result_id]][[matrix_name]] <- plot_mat
      z_values <- c(z_values, comparison_row_z_values(plot_mat))
    }
  }
  z_values <- z_values[is.finite(z_values)]
  global_color_limit <- if (!is.null(color_limit)) {
    color_limit
  } else if (length(z_values) > 0L) {
    max(2, ceiling(max(abs(z_values), na.rm = TRUE) * 10) / 10)
  } else {
    2
  }
  if (!is.numeric(global_color_limit) || length(global_color_limit) != 1L || is.na(global_color_limit) || global_color_limit <= 0) {
    stop("color_limit must be one positive number.", call. = FALSE)
  }
  color_breaks <- seq(-global_color_limit, global_color_limit, length.out = 102)
  legend_breaks <- c(-global_color_limit, 0, global_color_limit)
  legend_labels <- c(sprintf("%.1f", -global_color_limit), "0", sprintf("%.1f", global_color_limit))

  plot_gtables <- vector("list", length(result_ids))
  names(plot_gtables) <- result_ids
  for (result_id in result_ids) {
    plot_gtables[[result_id]] <- vector("list", length(matrix_names))
    names(plot_gtables[[result_id]]) <- matrix_names
    for (matrix_name in matrix_names) {
      p <- plot_target_matrix_full(
        target_mat_raw = selected_panel_mats[[result_id]][[matrix_name]],
        ensembl_to_symbol = gene_label_map,
        condition_label = sub("^target_gene_matrix_", "", matrix_name),
        top_n_genes = top_k,
        top_n_TE = top_k,
        highlight_TEs = highlight_TEs,
        selected_genes = selected_genes,
        selected_TEs = selected_tes,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        legend = FALSE,
        annotation_legend = FALSE,
        show_col_annotation = show_panel_annotations,
        show_row_annotation = show_panel_annotations,
        color_breaks = color_breaks,
        legend_breaks = legend_breaks,
        legend_labels = legend_labels,
        main = "",
        out_dir = NULL,
        show = FALSE
      )
      plot_gtables[[result_id]][[matrix_name]] <- p$gtable
    }
  }

  result_label_map <- comparison_activity_result_label_map(comparison_results, comparison_class)
  result_labels <- unname(result_label_map[result_ids])
  result_labels[is.na(result_labels) | !nzchar(result_labels)] <- comparison_pretty_label(result_ids[is.na(result_labels) | !nzchar(result_labels)])
  matrix_labels <- comparison_pretty_label(sub("^target_gene_matrix_", "", matrix_names))
  width <- max(7, 3.8 * length(matrix_names) + 1.4)
  height <- max(6, 3.2 * length(result_ids) + 0.9)

  draw_frame <- function() {
    grid::grid.newpage()
    grid::grid.text(
      sprintf("Comparable target gene matrix heatmaps: global top %d", top_k),
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.985, "npc"),
      gp = grid::gpar(fontsize = 13, fontface = "bold")
    )
    grid::grid.text(
      sprintf("Rows: result directories | columns: target matrices/contrasts | fixed %d genes x %d TE subfamilies | global gene/TE order | row z-score within each panel",
              length(selected_genes), length(selected_tes)),
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.955, "npc"),
      gp = grid::gpar(fontsize = 9, col = "grey35")
    )

    grid::pushViewport(grid::viewport(
      y = grid::unit(0.505, "npc"),
      height = grid::unit(0.84, "npc"),
      layout = grid::grid.layout(
        nrow = length(result_ids) + 1L,
        ncol = length(matrix_names) + 1L,
        widths = grid::unit.c(grid::unit(1.15, "in"), rep(grid::unit(1, "null"), length(matrix_names))),
        heights = grid::unit.c(grid::unit(0.35, "in"), rep(grid::unit(1, "null"), length(result_ids)))
      )
    ))

    for (j in seq_along(matrix_names)) {
      grid::grid.text(
        matrix_labels[j],
        vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = j + 1L),
        gp = grid::gpar(fontsize = 9, fontface = "bold")
      )
    }

    for (i in seq_along(result_ids)) {
      grid::grid.text(
        result_labels[i],
        rot = 90,
        vp = grid::viewport(layout.pos.row = i + 1L, layout.pos.col = 1L),
        gp = grid::gpar(fontsize = 8, fontface = "bold")
      )

      for (j in seq_along(matrix_names)) {
        grid::pushViewport(grid::viewport(layout.pos.row = i + 1L, layout.pos.col = j + 1L, clip = "on"))
        grid::grid.draw(plot_gtables[[result_ids[i]]][[matrix_names[j]]])
        grid::upViewport()
      }
    }

    grid::upViewport()
    if (!isTRUE(show_panel_annotations)) {
      comparison_draw_target_annotation_legend(selected_tes, highlight_TEs = highlight_TEs)
    }
    comparison_draw_target_color_legend(global_color_limit)
  }

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    if (requireNamespace("Cairo", quietly = TRUE)) {
      Cairo::CairoPDF(out_file, width = width, height = height)
    } else {
      grDevices::pdf(out_file, width = width, height = height)
    }
    draw_frame()
    grDevices::dev.off()
    message("Saved: ", out_file)
  }

  if (isTRUE(show)) {
    draw_frame()
  }

  invisible(list(
    status = "ok",
    common_target_matrices = matrix_names,
    selected_genes = selected_genes,
    selected_tes = selected_tes,
    gene_score = gene_score,
    te_score = te_score,
    color_limit = global_color_limit,
    plot_gtables = plot_gtables
  ))
}

summarise_comparison_final <- function(comparison_results) {
  class_info <- classify_comparison_set(comparison_results)
  fit <- summarise_activity_fit(comparison_results)
  sparsity <- summarise_sparsity(comparison_results)
  lambda <- summarise_lambda_best(comparison_results)
  contrasts <- summarise_te_contrasts(comparison_results)
  targets <- summarise_target_gene_matrices(comparison_results)

  overall_fit <- fit[fit$scope == "overall", , drop = FALSE]
  if (nrow(overall_fit) == 0L) {
    overall_fit <- fit
  }

  rows <- lapply(comparison_results, function(result) {
    fit_row <- overall_fit[overall_fit$result_id == result$result_id, , drop = FALSE]
    lambda_row <- lambda[lambda$result_id == result$result_id, , drop = FALSE]
    activity_sparse <- sparsity$activity_sparsity[sparsity$activity_sparsity$result_id == result$result_id, , drop = FALSE]

    n_sig <- NA_integer_
    if (identical(contrasts$status, "ok")) {
      sig_row <- contrasts$significant_summary[contrasts$significant_summary$result_id == result$result_id, , drop = FALSE]
      if (nrow(sig_row) > 0L) {
        n_sig <- sig_row$n_significant_te_subfamilies[1]
      }
    }

    data.frame(
      result_id = result$result_id,
      comparison_class = class_info$class_label,
      mean_r2 = if ("mean_r2" %in% colnames(fit_row) && nrow(fit_row) > 0L) fit_row$mean_r2[1] else NA_real_,
      mean_r2_ratio = if ("mean_r2_ratio" %in% colnames(fit_row) && nrow(fit_row) > 0L) fit_row$mean_r2_ratio[1] else NA_real_,
      mean_sparsity_percent = if ("mean_sparsity_percent" %in% colnames(fit_row) && nrow(fit_row) > 0L) fit_row$mean_sparsity_percent[1] else NA_real_,
      activity_zero_fraction = if (nrow(activity_sparse) > 0L) activity_sparse$zero_fraction_from_activity_scores[1] else NA_real_,
      lambda_best_summary = if (nrow(lambda_row) > 0L && isTRUE(lambda_row$lambda_best_available[1])) {
        if (isTRUE(lambda_row$global_lambda_available[1])) {
          sprintf("sample median=%g; global=%g", lambda_row$median_lambda_best[1], lambda_row$global_lambda_min[1])
        } else {
          sprintf("sample median=%g", lambda_row$median_lambda_best[1])
        }
      } else {
        "not available"
      },
      n_significant_te_subfamilies = n_sig,
      target_gene_support_summary = if (identical(targets$status, "ok")) "available" else targets$reason,
      warnings_or_skipped_blocks = paste(
        c(result$missing_optional, if (!identical(contrasts$status, "ok")) "te_contrast_summary block skipped" else NULL),
        collapse = ", "
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
