# Helper functions for scripts/gene_mapping.Rmd.
# Kept separate so the Rmd report can focus on the mapping narrative.

parse_cli_args <- function(args) {
  if ("--args" %in% args) {
    args <- args[(match("--args", args) + 1):length(args)]
  }
  if (length(args) == 0) {
    return(list())
  }
  if (length(args) %% 2 != 0) {
    stop("Command-line arguments must be provided as --flag value pairs.", call. = FALSE)
  }

  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (identical(key, args[[i]])) {
      stop(sprintf("Invalid argument '%s'; expected a --flag.", args[[i]]), call. = FALSE)
    }
    parsed[[gsub("-", "_", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  parsed
}

arg_or_null <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || identical(value, "")) NULL else value
}

require_arg <- function(args, name, flag) {
  value <- arg_or_null(args, name)
  if (is.null(value)) {
    stop(sprintf("Missing required argument %s.", flag), call. = FALSE)
  }
  value
}

path_stem <- function(path) {
  stem <- basename(path)
  stem <- sub("\\.gz$", "", stem)
  stem <- sub("\\.[^.]+$", "", stem)
  stem <- sub("_gene_names$", "", stem)
  stem <- sub("_genes$", "", stem)
  stem
}

absolute_path <- function(path, must_work = FALSE) {
  path <- path.expand(path)
  if (file.exists(path) || must_work) {
    return(normalizePath(path, mustWork = must_work))
  }

  file.path(
    normalizePath(dirname(path), mustWork = FALSE),
    basename(path)
  )
}

derive_results_section_dir <- function(e_path,
                                       n_path,
                                       results_parent = NULL,
                                       prefix = "results") {
  e_base <- path_stem(e_path)
  n_base <- path_stem(n_path)

  if (is.null(results_parent)) {
    results_parent <- file.path(getwd(), "results")
  }
  results_parent <- normalizePath(path.expand(results_parent), mustWork = FALSE)

  file.path(results_parent, paste0(prefix, "_", e_base, "_", n_base))
}

derive_mapping_dir <- function(e_path,
                               n_path,
                               results_parent = NULL,
                               prefix = "results") {
  file.path(
    derive_results_section_dir(
      e_path = e_path,
      n_path = n_path,
      results_parent = results_parent,
      prefix = prefix
    ),
    "mapping"
  )
}

derive_output_path <- function(e_path,
                               n_path,
                               results_parent = NULL,
                               prefix = "results") {
  e_base <- path_stem(e_path)
  n_base <- path_stem(n_path)
  mapping_id <- paste0("map_", e_base, "_", n_base)

  file.path(
    derive_mapping_dir(
      e_path = e_path,
      n_path = n_path,
      results_parent = results_parent,
      prefix = prefix
    ),
    paste0(mapping_id, ".tsv")
  )
}

derive_legacy_output_path <- function(e_path, n_path) {
  e_base <- path_stem(e_path)
  n_base <- path_stem(n_path)
  mapping_id <- paste0("map_", e_base, "_", n_base)
  e_dir <- dirname(e_path)

  if (basename(e_dir) == mapping_id) {
    return(file.path(e_dir, paste0(mapping_id, ".tsv")))
  }

  file.path(e_dir, "mappings", mapping_id, paste0(mapping_id, ".tsv"))
}

copy_legacy_mapping_cache <- function(mapping_cache_path, legacy_mapping_cache_path) {
  if (is.null(legacy_mapping_cache_path) ||
      identical(mapping_cache_path, legacy_mapping_cache_path) ||
      file.exists(mapping_cache_path) ||
      !file.exists(legacy_mapping_cache_path)) {
    return(invisible(FALSE))
  }

  dir.create(dirname(mapping_cache_path), showWarnings = FALSE, recursive = TRUE)

  src_paths <- c(
    cache = legacy_mapping_cache_path,
    report = derive_report_path(legacy_mapping_cache_path),
    pdf_report = derive_pdf_report_path(legacy_mapping_cache_path)
  )
  dest_paths <- c(
    cache = mapping_cache_path,
    report = derive_report_path(mapping_cache_path),
    pdf_report = derive_pdf_report_path(mapping_cache_path)
  )

  for (name in names(src_paths)) {
    src <- src_paths[[name]]
    dest <- dest_paths[[name]]
    if (file.exists(src) && !file.exists(dest)) {
      file.copy(src, dest, overwrite = FALSE)
    }
  }

  legacy_diagnostics_dir <- file.path(dirname(legacy_mapping_cache_path), "diagnostics")
  diagnostics_dir <- file.path(dirname(mapping_cache_path), "diagnostics")
  if (dir.exists(legacy_diagnostics_dir) && !dir.exists(diagnostics_dir)) {
    dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(
      list.files(legacy_diagnostics_dir, full.names = TRUE, all.files = FALSE, no.. = TRUE),
      diagnostics_dir,
      overwrite = FALSE,
      recursive = FALSE
    )
  }

  message("Copied legacy mapping cache into results mapping dir: ", mapping_cache_path)
  invisible(TRUE)
}

derive_report_path <- function(mapping_cache_path) {
  cache_stem <- path_stem(mapping_cache_path)
  mapping_pair_id <- sub("^map_", "", cache_stem)

  file.path(
    dirname(mapping_cache_path),
    paste0("report_", mapping_pair_id, ".html")
  )
}

derive_pdf_report_path <- function(mapping_cache_path) {
  cache_stem <- path_stem(mapping_cache_path)
  mapping_pair_id <- sub("^map_", "", cache_stem)

  file.path(
    dirname(mapping_cache_path),
    paste0("report_", mapping_pair_id, ".pdf")
  )
}

derive_status_path <- function(mapping_cache_path) {
  file.path(dirname(mapping_cache_path), "mapping_status.tsv")
}

latex_engine_available <- function() {
  engines <- Sys.which(c("pdflatex", "xelatex", "lualatex"))
  any(nzchar(unname(engines)))
}

configure_quarto_pandoc <- function() {
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    return(invisible(TRUE))
  }

  quarto_bin <- Sys.which("quarto")
  if (identical(unname(quarto_bin), "")) {
    return(invisible(FALSE))
  }

  quarto_paths <- tryCatch(
    system2(quarto_bin, "--paths", stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(quarto_paths) == 0) {
    return(invisible(FALSE))
  }

  quarto_tools <- file.path(quarto_paths[[1]], "tools")
  pandoc_dirs <- unique(c(
    file.path(quarto_tools, Sys.info()[["machine"]]),
    file.path(quarto_tools, R.version$arch),
    file.path(quarto_tools, "aarch64"),
    file.path(quarto_tools, "x86_64")
  ))
  pandoc_dirs <- pandoc_dirs[file.exists(file.path(pandoc_dirs, "pandoc"))]

  if (length(pandoc_dirs) == 0) {
    return(invisible(FALSE))
  }

  Sys.setenv(RSTUDIO_PANDOC = pandoc_dirs[[1]])
  invisible(TRUE)
}

mapping_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
}

write_mapping_status <- function(mapping_cache_path,
                                 status,
                                 stage,
                                 message = "",
                                 report_path = derive_report_path(mapping_cache_path),
                                 started_at = NULL,
                                 started_at_epoch = NULL) {
  status_path <- derive_status_path(mapping_cache_path)
  dir.create(dirname(status_path), showWarnings = FALSE, recursive = TRUE)

  if (is.null(started_at) || is.null(started_at_epoch)) {
    old_status <- tryCatch(
      data.table::fread(status_path, sep = "\t"),
      error = function(e) NULL
    )

    if (!is.null(old_status) && nrow(old_status) > 0) {
      if (is.null(started_at) && "started_at" %in% colnames(old_status)) {
        started_at <- as.character(old_status$started_at[[1]])
      }
      if (is.null(started_at_epoch) && "started_at_epoch" %in% colnames(old_status)) {
        started_at_epoch <- suppressWarnings(as.numeric(old_status$started_at_epoch[[1]]))
      }
    }
  }

  now <- Sys.time()
  if (is.null(started_at)) {
    started_at <- format(now, "%Y-%m-%d %H:%M:%S %Z")
  }
  if (is.null(started_at_epoch) || is.na(started_at_epoch)) {
    started_at_epoch <- as.numeric(now)
  }

  status_table <- data.table::data.table(
    status = status,
    stage = stage,
    started_at = started_at,
    updated_at = format(now, "%Y-%m-%d %H:%M:%S %Z"),
    started_at_epoch = started_at_epoch,
    updated_at_epoch = as.numeric(now),
    cache_path = mapping_cache_path,
    report_path = report_path,
    message = message
  )

  data.table::fwrite(status_table, status_path, sep = "\t")
  invisible(status_path)
}

read_mapping_status <- function(mapping_cache_path) {
  status_path <- derive_status_path(mapping_cache_path)
  if (!file.exists(status_path)) {
    return(NULL)
  }

  status_table <- tryCatch(
    data.table::fread(status_path, sep = "\t"),
    error = function(e) NULL
  )

  if (is.null(status_table) || nrow(status_table) == 0) {
    return(NULL)
  }

  status_table[1]
}

cleanup_mapping_transient_files <- function(mapping_cache_path,
                                            e_gene_names_path,
                                            n_gene_names_path) {
  transient_paths <- as.character(c(
    derive_status_path(mapping_cache_path),
    e_gene_names_path,
    n_gene_names_path
  ))
  transient_paths <- unique(transient_paths[!is.na(transient_paths) & nzchar(transient_paths)])

  # Remove only rebuild scratch files after the final cache has been written.
  transient_paths <- transient_paths[file.exists(transient_paths) & !dir.exists(transient_paths)]
  unlink(transient_paths)

  invisible(transient_paths)
}

read_gene_tsv <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s gene-name path does not exist: %s", label, path), call. = FALSE)
  }
  if (file.access(path, mode = 4) != 0) {
    stop(sprintf("%s gene-name path is not readable: %s", label, path), call. = FALSE)
  }

  dt <- tryCatch(
    data.table::fread(path, sep = "\t", header = TRUE),
    error = function(e) {
      stop(sprintf("%s gene-name TSV could not be read: %s", label, conditionMessage(e)), call. = FALSE)
    }
  )
  if (ncol(dt) < 1) {
    stop(sprintf("%s gene-name TSV has no columns: %s", label, path), call. = FALSE)
  }

  genes <- as.character(dt[[1]])
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) {
    stop(sprintf("%s gene-name TSV contains no usable gene identifiers: %s", label, path), call. = FALSE)
  }
  genes
}

read_default_e_gene_names <- function(path, handling = "discard", duplicate_threshold = 10) {
  if (!file.exists(path)) {
    stop(sprintf("Default E matrix path does not exist: %s", path), call. = FALSE)
  }

  e_dt <- data.table::fread(path)
  gene_col <- if ("Geneid" %in% colnames(e_dt)) "Geneid" else colnames(e_dt)[1]
  gene_ids <- as.character(e_dt[[gene_col]])

  drop_idx <- match(c(gene_col, "Length"), colnames(e_dt), nomatch = 0)
  drop_idx <- unique(drop_idx[drop_idx > 0])
  e_mat <- as.matrix(e_dt[, -drop_idx, with = FALSE])
  rownames(e_mat) <- gene_ids

  e_mat <- resolve_duplicates(e_mat, method = handling, threshold = duplicate_threshold)
  transform_excel_dates(rownames(e_mat))
}

read_default_n_gene_names <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Default N matrix path does not exist: %s", path), call. = FALSE)
  }

  n_dt <- data.table::fread(path, header = TRUE, sep = "\t")
  genes <- as.character(n_dt[[1]])
  genes[!is.na(genes) & genes != ""]
}

write_gene_name_inputs <- function(mapping_cache_path,
                                   e_source_path,
                                   n_source_path,
                                   gene_ids,
                                   genes_N) {
  mapping_dir <- dirname(mapping_cache_path)
  dir.create(mapping_dir, showWarnings = FALSE, recursive = TRUE)

  e_gene_names_path <- file.path(mapping_dir, paste0(path_stem(e_source_path), "_gene_names.tsv"))
  n_gene_names_path <- file.path(mapping_dir, paste0(path_stem(n_source_path), "_gene_names.tsv"))

  data.table::fwrite(
    data.table::data.table(hgnc_symbol = gene_ids),
    e_gene_names_path,
    sep = "\t"
  )
  data.table::fwrite(
    data.table::data.table(ensembl_gene_id = genes_N),
    n_gene_names_path,
    sep = "\t"
  )

  list(
    e_gene_names_path = e_gene_names_path,
    n_gene_names_path = n_gene_names_path
  )
}

mapping_cache_missing_message <- function(mapping_cache_path,
                                          mapping_report_path = derive_report_path(mapping_cache_path),
                                          stale_after_hours = 12) {
  status <- read_mapping_status(mapping_cache_path)
  mart_hint <- if (grepl("hg19|grch37", tolower(mapping_cache_path))) {
    "For Hg19/GRCh37 caches, rebuild with --primary-mart grch37."
  } else {
    "For Hg38/current caches, rebuild with --primary-mart current."
  }
  rebuild_hint <- paste(
    "Uncomment and run the rebuild command in section 3.2 once with server access,",
    "then rerun section 3.",
    mart_hint
  )

  if (is.null(status)) {
    return(paste(
      "Mapping cache has not been built yet:", mapping_cache_path,
      "\n", rebuild_hint
    ))
  }

  state <- if ("status" %in% names(status)) {
    unname(as.character(status[["status"]][1]))
  } else {
    "unknown"
  }
  stage <- if ("stage" %in% names(status)) {
    unname(as.character(status[["stage"]][1]))
  } else {
    "unknown stage"
  }
  updated <- if ("updated_at" %in% names(status)) {
    unname(as.character(status[["updated_at"]][1]))
  } else {
    "unknown time"
  }
  status_message <- if ("message" %in% names(status)) {
    unname(as.character(status[["message"]][1]))
  } else {
    ""
  }
  report_line <- if (file.exists(mapping_report_path)) {
    paste0("\nReport: ", mapping_report_path)
  } else {
    paste0("\nExpected report: ", mapping_report_path)
  }

  if (isTRUE(state == "running")) {
    age_hours <- NA_real_
    if ("updated_at_epoch" %in% names(status)) {
      updated_epoch <- suppressWarnings(as.numeric(status$updated_at_epoch[1]))
      if (!is.na(updated_epoch)) {
        age_hours <- as.numeric(difftime(
          Sys.time(),
          as.POSIXct(updated_epoch, origin = "1970-01-01"),
          units = "hours"
        ))
      }
    }

    if (!is.na(age_hours) && age_hours > stale_after_hours) {
      return(paste(
        "Mapping rebuild may have been interrupted before writing the cache:", mapping_cache_path,
        "\nLast status:", stage, "at", updated,
        "\n", status_message,
        "\nRerun the rebuild cell in section 3.2.",
        report_line
      ))
    }

    return(paste(
      "Mapping rebuild appears to still be in progress; wait until it finishes, then rerun this cell.",
      "\nCache:", mapping_cache_path,
      "\nLast status:", stage, "at", updated,
      "\n", status_message,
      report_line
    ))
  }

  if (isTRUE(state == "failed")) {
    return(paste(
      "Last mapping rebuild failed before writing the cache:", mapping_cache_path,
      "\nLast failed stage:", stage, "at", updated,
      "\n", status_message,
      "\nCheck the mapping report/console output, rerun the rebuild cell, then rerun section 3.",
      report_line
    ))
  }

  if (isTRUE(state == "complete")) {
    return(paste(
      "Mapping status says complete, but the cache TSV is missing:", mapping_cache_path,
      "\nThe file may have been moved or deleted. Rerun the rebuild cell in section 3.2.",
      report_line
    ))
  }

  paste(
    "Mapping cache not found:", mapping_cache_path,
    "\nLast rebuild status:", state, "-", stage, "at", updated,
    "\n", rebuild_hint,
    report_line
  )
}

prepare_mapping_cache_inputs <- function(E_path,
                                         N_path,
                                         gene_ids,
                                         genes_N) {
  if (is.null(E_path) || !nzchar(E_path)) {
    stop("E_path is not defined. Set E_path in the active E loading cell before section 3.", call. = FALSE)
  }
  if (is.null(N_path) || !nzchar(N_path)) {
    stop("N_path is not defined. Set N_path in the active N loading cell before section 3.", call. = FALSE)
  }

  E_path <- absolute_path(E_path, must_work = TRUE)
  N_path <- absolute_path(N_path, must_work = TRUE)

  mapping_cache_path <- derive_output_path(E_path, N_path)
  mapping_cache_path <- absolute_path(mapping_cache_path, must_work = FALSE)
  copy_legacy_mapping_cache(
    mapping_cache_path = mapping_cache_path,
    legacy_mapping_cache_path = absolute_path(derive_legacy_output_path(E_path, N_path), must_work = FALSE)
  )
  mapping_report_path <- derive_report_path(mapping_cache_path)
  mapping_status_path <- derive_status_path(mapping_cache_path)

  gene_name_paths <- write_gene_name_inputs(
    mapping_cache_path = mapping_cache_path,
    e_source_path = E_path,
    n_source_path = N_path,
    gene_ids = gene_ids,
    genes_N = genes_N
  )

  list(
    mapping_cache_path = mapping_cache_path,
    mapping_report_path = mapping_report_path,
    mapping_status_path = mapping_status_path,
    e_gene_names_path = gene_name_paths$e_gene_names_path,
    n_gene_names_path = gene_name_paths$n_gene_names_path
  )
}

ensure_mapping_cache_inputs <- function(E_path,
                                        N_path,
                                        gene_ids,
                                        genes_N,
                                        mapping_inputs = NULL) {
  missing_inputs <- is.null(mapping_inputs) ||
    is.null(mapping_inputs$e_gene_names_path) ||
    is.null(mapping_inputs$n_gene_names_path) ||
    !file.exists(mapping_inputs$e_gene_names_path) ||
    !file.exists(mapping_inputs$n_gene_names_path)

  if (isTRUE(missing_inputs)) {
    mapping_inputs <- prepare_mapping_cache_inputs(
      E_path = E_path,
      N_path = N_path,
      gene_ids = gene_ids,
      genes_N = genes_N
    )
  }

  handoff_paths <- c(
    E = mapping_inputs$e_gene_names_path,
    N = mapping_inputs$n_gene_names_path
  )
  missing_paths <- handoff_paths[!file.exists(handoff_paths)]

  if (length(missing_paths) > 0) {
    stop(
      paste(
        "Mapping gene-name handoff TSVs are missing and could not be written:",
        paste(sprintf("%s=%s", names(missing_paths), missing_paths), collapse = ", "),
        "\nRerun section 3.1 after loading the active E and N matrices."
      ),
      call. = FALSE
    )
  }

  mapping_inputs
}

primary_mart_for_reference_genome <- function(reference_genome) {
  reference_key <- tolower(as.character(reference_genome))
  if (length(reference_key) != 1L || is.na(reference_key) || !nzchar(reference_key)) {
    stop("reference_genome must be one of 'Hg19' or 'Hg38'.", call. = FALSE)
  }

  if (reference_key %in% c("hg19", "grch37")) {
    return("grch37")
  }
  if (reference_key %in% c("hg38", "grch38")) {
    return("current")
  }

  stop("reference_genome must be one of 'Hg19' or 'Hg38'.", call. = FALSE)
}

mapping_reference_genome_from_primary_mart <- function(primary_mart) {
  primary_key <- tolower(as.character(primary_mart))
  if (length(primary_key) != 1L || is.na(primary_key) || !nzchar(primary_key)) {
    stop("primary_mart must be a non-empty string.", call. = FALSE)
  }

  if (primary_key %in% c("grch37", "grch_37", "grch-37", "hg19", "75", "v75")) {
    return("Hg19")
  }

  "Hg38"
}

prepare_mapping_inputs <- function(cli_args,
                                   default_e_path,
                                   default_n_path,
                                   default_primary_mart = "current",
                                   default_mirror = "www",
                                   default_handling = "discard",
                                   default_duplicate_threshold = 10) {
  args <- parse_cli_args(cli_args)
  use_defaults <- length(args) == 0

  primary_mart <- arg_or_null(args, "primary_mart")
  mirror <- arg_or_null(args, "mirror")
  if (is.null(primary_mart)) primary_mart <- default_primary_mart
  if (is.null(mirror)) mirror <- default_mirror
  mapping_reference_genome <- mapping_reference_genome_from_primary_mart(primary_mart)

  if (use_defaults) {
    default_e_path <- absolute_path(default_e_path, must_work = TRUE)
    default_n_path <- absolute_path(default_n_path, must_work = TRUE)

    gene_ids <- read_default_e_gene_names(
      default_e_path,
      handling = default_handling,
      duplicate_threshold = default_duplicate_threshold
    )
    genes_N <- unique(read_default_n_gene_names(default_n_path))
    mapping_cache_path <- absolute_path(
      derive_output_path(default_e_path, default_n_path),
      must_work = FALSE
    )
    copy_legacy_mapping_cache(
      mapping_cache_path = mapping_cache_path,
      legacy_mapping_cache_path = absolute_path(derive_legacy_output_path(default_e_path, default_n_path), must_work = FALSE)
    )

    gene_name_paths <- write_gene_name_inputs(
      mapping_cache_path = mapping_cache_path,
      e_source_path = default_e_path,
      n_source_path = default_n_path,
      gene_ids = gene_ids,
      genes_N = genes_N
    )

    return(c(
      list(
        default_inputs = TRUE,
        e_source_path = default_e_path,
        n_source_path = default_n_path,
        mapping_cache_path = mapping_cache_path,
        primary_mart = primary_mart,
        mapping_reference_genome = mapping_reference_genome,
        mirror = mirror,
        gene_ids = gene_ids,
        genes_N = genes_N
      ),
      gene_name_paths
    ))
  }

  e_gene_names_path <- absolute_path(
    require_arg(args, "e_gene_names", "--e-gene-names"),
    must_work = TRUE
  )
  n_gene_names_path <- absolute_path(
    require_arg(args, "n_gene_names", "--n-gene-names"),
    must_work = TRUE
  )
  mapping_cache_path <- arg_or_null(args, "out")
  if (is.null(mapping_cache_path)) {
    mapping_cache_path <- derive_output_path(e_gene_names_path, n_gene_names_path)
  }
  mapping_cache_path <- absolute_path(mapping_cache_path, must_work = FALSE)

  gene_ids <- read_gene_tsv(e_gene_names_path, "E")
  genes_N <- unique(read_gene_tsv(n_gene_names_path, "N"))

  list(
    default_inputs = FALSE,
    e_source_path = e_gene_names_path,
    n_source_path = n_gene_names_path,
    e_gene_names_path = e_gene_names_path,
    n_gene_names_path = n_gene_names_path,
    mapping_cache_path = mapping_cache_path,
    primary_mart = primary_mart,
    mapping_reference_genome = mapping_reference_genome,
    mirror = mirror,
    gene_ids = gene_ids,
    genes_N = genes_N
  )
}

use_ensembl_with_mirror_fallback <- function(dataset = "hsapiens_gene_ensembl",
                                             mirror = "www",
                                             version = NULL,
                                             grch37 = FALSE) {
  if (isTRUE(grch37)) {
    mart <- tryCatch(
      biomaRt::useEnsembl(
        biomart = "genes",
        dataset = dataset,
        GRCh = 37
      ),
      error = function(e) {
        stop(
          sprintf(
            "Unable to query Ensembl GRCh37 BioMart for dataset '%s'. Last error: %s",
            dataset,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
    return(mart)
  }

  mirrors <- unique(c(mirror, "www", "useast", "asia"))
  mart_label <- if (is.null(version)) {
    "current Ensembl BioMart"
  } else {
    sprintf("Ensembl BioMart version %s", version)
  }
  last_error <- NULL

  for (candidate_mirror in mirrors) {
    mart <- tryCatch(
      {
        if (is.null(version)) {
          biomaRt::useEnsembl("genes", dataset, mirror = candidate_mirror)
        } else {
          biomaRt::useEnsembl("genes", dataset, version = version, mirror = candidate_mirror)
        }
      },
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(mart)) {
      if (!identical(candidate_mirror, mirror)) {
        message(sprintf(
          "Ensembl mirror '%s' was unavailable; using mirror '%s' instead.",
          mirror,
          candidate_mirror
        ))
      }
      return(mart)
    }
  }

  stop(
    sprintf(
      paste(
        "Unable to query %s for dataset '%s' using mirrors: %s.",
        "If a specific Ensembl version was requested, that BioMart archive may be temporarily unavailable",
        "or no longer served by Ensembl.",
        "Last error: %s"
      ),
      mart_label,
      dataset,
      paste(mirrors, collapse = ", "),
      conditionMessage(last_error)
    ),
    call. = FALSE
  )
}

make_ensembl_mart <- function(primary_mart = "current", mirror = "www") {
  primary_mart <- tolower(primary_mart)
  if (primary_mart %in% c("current", "latest")) {
    return(use_ensembl_with_mirror_fallback(mirror = mirror))
  }
  if (primary_mart %in% c("grch37", "grch_37", "grch-37", "hg19")) {
    return(use_ensembl_with_mirror_fallback(mirror = mirror, grch37 = TRUE))
  }

  version <- sub("^v", "", primary_mart)
  if (!grepl("^[0-9]+$", version)) {
    stop("--primary-mart must be 'current', 'latest', 'grch37', 'hg19', 'v80', or a numeric Ensembl version.", call. = FALSE)
  }

  use_ensembl_with_mirror_fallback(mirror = mirror, version = as.integer(version))
}

empty_biomart_biotype_table <- function(include_symbol = FALSE) {
  out <- data.frame(
    ensembl_gene_id = character(),
    gene_biotype = character(),
    stringsAsFactors = FALSE
  )
  if (include_symbol) {
    out$hgnc_symbol <- character()
    out <- out[, c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), drop = FALSE]
  }
  out
}

empty_biomart_symbol_table <- function() {
  data.frame(
    hgnc_symbol = character(),
    ensembl_gene_id = character(),
    gene_biotype = character(),
    stringsAsFactors = FALSE
  )
}

query_biomart_symbols <- function(symbols,
                                  mart,
                                  chunk_size = 2000) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0L || is.null(mart)) {
    return(empty_biomart_symbol_table())
  }

  attributes <- c("hgnc_symbol", "ensembl_gene_id", "gene_biotype")
  symbol_chunks <- split(symbols, ceiling(seq_along(symbols) / chunk_size))
  result_chunks <- lapply(symbol_chunks, function(chunk_symbols) {
    tryCatch(
      biomaRt::getBM(
        attributes = attributes,
        filters = "hgnc_symbol",
        values = chunk_symbols,
        mart = mart
      ),
      error = function(e) {
        warning(
          sprintf(
            "Skipping a BioMart symbol chunk after query failure: %s",
            conditionMessage(e)
          ),
          call. = FALSE
        )
        empty_biomart_symbol_table()
      }
    )
  })

  result <- do.call(rbind, result_chunks)
  if (is.null(result) || nrow(result) == 0L) {
    return(empty_biomart_symbol_table())
  }

  result <- as.data.frame(result, stringsAsFactors = FALSE)
  result <- result[!is.na(result$hgnc_symbol) & result$hgnc_symbol != "", , drop = FALSE]
  unique(result)
}

query_biomart_biotypes <- function(ensembl_ids,
                                   mart,
                                   include_symbol = FALSE,
                                   chunk_size = 2000) {
  ensembl_ids <- unique(as.character(ensembl_ids))
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids) & ensembl_ids != ""]
  if (length(ensembl_ids) == 0L || is.null(mart)) {
    return(empty_biomart_biotype_table(include_symbol = include_symbol))
  }

  attributes <- c("ensembl_gene_id", "gene_biotype")
  if (include_symbol) {
    attributes <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
  }

  id_chunks <- split(ensembl_ids, ceiling(seq_along(ensembl_ids) / chunk_size))
  result_chunks <- lapply(id_chunks, function(ids) {
    tryCatch(
      biomaRt::getBM(
        attributes = attributes,
        filters = "ensembl_gene_id",
        values = ids,
        mart = mart
      ),
      error = function(e) {
        warning(
          sprintf(
            "Skipping a BioMart biotype chunk after query failure: %s",
            conditionMessage(e)
          ),
          call. = FALSE
        )
        empty_biomart_biotype_table(include_symbol = include_symbol)
      }
    )
  })

  result <- do.call(rbind, result_chunks)

  if (is.null(result) || nrow(result) == 0L) {
    return(empty_biomart_biotype_table(include_symbol = include_symbol))
  }

  result <- as.data.frame(result, stringsAsFactors = FALSE)
  result <- result[!is.na(result$ensembl_gene_id) & result$ensembl_gene_id != "", , drop = FALSE]
  unique(result)
}

recover_biomart_biotypes <- function(ensembl_ids,
                                     mart_first,
                                     mart_fallback = NULL,
                                     include_symbol = FALSE,
                                     chunk_size = 2000) {
  ensembl_ids <- unique(as.character(ensembl_ids))
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids) & ensembl_ids != ""]
  if (length(ensembl_ids) == 0L) {
    return(empty_biomart_biotype_table(include_symbol = include_symbol))
  }

  first <- query_biomart_biotypes(
    ensembl_ids,
    mart = mart_first,
    include_symbol = include_symbol,
    chunk_size = chunk_size
  )

  missing_from_first <- ensembl_ids[!ensembl_ids %in% first$ensembl_gene_id]
  fallback <- if (!is.null(mart_fallback) && length(missing_from_first) > 0L) {
    query_biomart_biotypes(
      missing_from_first,
      mart = mart_fallback,
      include_symbol = include_symbol,
      chunk_size = chunk_size
    )
  } else {
    empty_biomart_biotype_table(include_symbol = include_symbol)
  }

  combined <- rbind(first, fallback)
  if (nrow(combined) == 0L) {
    return(empty_biomart_biotype_table(include_symbol = include_symbol))
  }

  combined <- combined[!is.na(combined$gene_biotype) & combined$gene_biotype != "", , drop = FALSE]
  combined[!duplicated(combined$ensembl_gene_id), , drop = FALSE]
}

write_cache_atomically <- function(mapping_table, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  tmp_path <- tempfile(
    pattern = paste0(basename(out_path), "."),
    tmpdir = dirname(out_path),
    fileext = ".tmp"
  )
  backup_path <- tempfile(
    pattern = paste0(basename(out_path), ".backup."),
    tmpdir = dirname(out_path)
  )
  on.exit(unlink(c(tmp_path, backup_path)), add = TRUE)

  data.table::fwrite(mapping_table, tmp_path, sep = "\t", na = "")

  had_existing <- file.exists(out_path)
  if (had_existing && !file.rename(out_path, backup_path)) {
    stop(sprintf("Could not prepare existing cache for replacement: %s", out_path), call. = FALSE)
  }

  if (!file.rename(tmp_path, out_path)) {
    if (had_existing && file.exists(backup_path)) {
      file.rename(backup_path, out_path)
    }
    stop(sprintf("Could not move temporary cache into final path: %s", out_path), call. = FALSE)
  }

  if (had_existing) {
    unlink(backup_path)
  }
}
