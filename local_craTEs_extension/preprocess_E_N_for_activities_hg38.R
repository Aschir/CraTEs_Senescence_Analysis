# hg38-compatible local preprocessing for craTEs workflows.
# Source this file after library(craTEs), then call preprocess_E_N_for_activities()
# without the craTEs:: namespace.

.get_pseudocounts_compat <- function(count_table, q = 0.05) {
  stopifnot(is.matrix(count_table))

  if (typeof(count_table) != "integer") {
    warning("The count table provided to get_pseudocounts() is not an integer matrix, use at your own risks")
  }

  pseudocounts <- apply(count_table, 2, function(x) {
    non_zero <- x[x != 0]
    if (length(non_zero) == 0) {
      return(0)
    }
    as.numeric(stats::quantile(non_zero, q, names = FALSE))
  })

  t(t(count_table) + pseudocounts)
}

.row_center_compat <- function(X) {
  stopifnot(is.matrix(X))
  t(scale(t(X), center = TRUE, scale = FALSE))
}

.make_txdb_from_ucsc_compat <- function(genome, tablename = "ensGene") {
  if (!requireNamespace("RMariaDB", quietly = TRUE)) {
    stop(
      "UCSC exon annotation requires package 'RMariaDB'. ",
      "Either install RMariaDB or use hg38_source = 'EnsDb.Hsapiens.v86' instead.",
      call. = FALSE
    )
  }

  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    return(txdbmaker::makeTxDbFromUCSC(genome = genome, tablename = tablename))
  }

  GenomicFeatures::makeTxDbFromUCSC(genome = genome, tablename = tablename)
}

get_exonic_lengths_for_org <- function(org = c("hg19", "hg38"),
                                       hg38_source = c("EnsDb.Hsapiens.v86", "UCSC")) {
  org <- match.arg(org)
  hg38_source <- match.arg(hg38_source)

  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required to compute exon lengths.", call. = FALSE)
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required to compute exon lengths.", call. = FALSE)
  }

  if (org == "hg38" && hg38_source == "UCSC") {
    stop(
      "UCSC ensGene exon annotation is not supported for hg38 in this workflow. ",
      "Use hg38_source = 'EnsDb.Hsapiens.v86' instead.",
      call. = FALSE
    )
  }

  if (org == "hg19" || hg38_source == "UCSC") {
    txdb <- .make_txdb_from_ucsc_compat(genome = org, tablename = "ensGene")
    exonic <- GenomicFeatures::exonsBy(txdb, by = "gene")
  } else {
    if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      stop("Package 'EnsDb.Hsapiens.v86' is required for hg38_source = 'EnsDb.Hsapiens.v86'.", call. = FALSE)
    }
    exonic <- GenomicFeatures::exonsBy(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, by = "gene")
  }

  reduced_exonic <- GenomicRanges::reduce(exonic)
  exon_lengths <- sum(GenomicRanges::width(reduced_exonic))
  names(exon_lengths) <- sub("\\..*$", "", names(exon_lengths))
  exon_lengths[!duplicated(names(exon_lengths))]
}

get_tpm_local_for_org <- function(countTable,
                                  org = c("hg19", "hg38"),
                                  exonic_lengths_path = NULL,
                                  redownload_exons = FALSE,
                                  exons_dir = "../data/temp/",
                                  hg38_source = c("EnsDb.Hsapiens.v86", "UCSC")) {
  stopifnot(is.matrix(countTable))

  org <- match.arg(org)
  hg38_source <- match.arg(hg38_source)

  if (is.null(exonic_lengths_path)) {
    source_tag <- if (org == "hg38") paste0("_", hg38_source) else ""
    exonic_lengths_path <- file.path(exons_dir, paste0(org, source_tag, "_exon_lengths.RDS"))
  }

  if (!dir.exists(dirname(exonic_lengths_path))) {
    dir.create(dirname(exonic_lengths_path), showWarnings = FALSE, recursive = TRUE)
  }

  if (!file.exists(exonic_lengths_path) || isTRUE(redownload_exons)) {
    if (file.exists(exonic_lengths_path)) {
      unlink(exonic_lengths_path)
    }
    exonic_lengths <- get_exonic_lengths_for_org(org = org, hg38_source = hg38_source)
    saveRDS(exonic_lengths, exonic_lengths_path)
  } else {
    exonic_lengths <- readRDS(exonic_lengths_path)
  }

  with_exons <- rownames(countTable)[rownames(countTable) %in% names(exonic_lengths)]
  if (length(with_exons) == 0) {
    stop("No count_table rownames were found in the exon-length annotation for org = '", org, "'.", call. = FALSE)
  }

  exonic_lengths_filtered <- exonic_lengths[with_exons]
  length_norm <- sweep(countTable[with_exons, , drop = FALSE], 1, exonic_lengths_filtered, "/")
  library_size <- colSums(length_norm)

  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("Cannot compute TPM: at least one sample has non-positive or non-finite length-normalized counts.", call. = FALSE)
  }

  sweep(length_norm, 2, library_size, "/") * 1e6
}

preprocess_E_N_for_activities <- function(count_table,
                                          N,
                                          count_threshold = c(10, 1),
                                          N_threshold = 150,
                                          pseudocount_quantile = 0.05,
                                          log_tpm_plot_path,
                                          preds_protected = NULL,
                                          redownload_exons = FALSE,
                                          exons_dir = "../data/temp/",
                                          org = c("hg19", "hg38"),
                                          hg38_source = c("EnsDb.Hsapiens.v86", "UCSC")) {
  org <- match.arg(org)
  hg38_source <- match.arg(hg38_source)

  stopifnot(is.matrix(count_table) && all(count_table >= 0))
  stopifnot(is.matrix(N) && all(N >= 0))

  if (typeof(count_table) != "integer") {
    warning("The count table contains non-integer values, proceed with caution")
  }

  genes_kept <- colSums(apply(count_table, 1, function(x) x >= count_threshold[1])) >= count_threshold[2]
  count_table <- count_table[genes_kept, , drop = FALSE]

  in_N <- rownames(count_table) %in% rownames(N)
  count_table <- count_table[in_N, , drop = FALSE]
  N <- N[rownames(count_table), , drop = FALSE]

  protected <- colnames(N) %in% preds_protected
  N_sums <- colSums(N)
  lost_cols <- N_sums[N_sums < N_threshold & !protected]
  N <- N[, N_sums >= N_threshold | protected, drop = FALSE]

  if (ncol(N) == 0) {
    stop("No TE predictors remained after applying N_threshold.", call. = FALSE)
  }

  to_keep <- rowSums(N) > 0
  N <- N[to_keep, , drop = FALSE]
  count_table <- count_table[rownames(N), , drop = FALSE]

  count_table <- .get_pseudocounts_compat(count_table, pseudocount_quantile)
  log_tpm <- log2(get_tpm_local_for_org(
    countTable = count_table,
    org = org,
    redownload_exons = redownload_exons,
    exons_dir = exons_dir,
    hg38_source = hg38_source
  ))

  N <- N[rownames(log_tpm), , drop = FALSE]
  gc()

  plot_dir <- dirname(log_tpm_plot_path)
  if (!dir.exists(plot_dir) && plot_dir != ".") {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  }
  grDevices::pdf(log_tpm_plot_path)
  graphics::hist(log_tpm)
  grDevices::dev.off()

  E <- scale(log_tpm, center = TRUE, scale = FALSE)
  gc()

  E <- .row_center_compat(E)

  stopifnot(all(rownames(E) == rownames(N)))

  list(
    E_centered = E,
    N = N,
    lost_cols = lost_cols
  )
}

preprocess_E_N_for_activities_hg38 <- function(..., org = "hg38") {
  preprocess_E_N_for_activities(..., org = org)
}
