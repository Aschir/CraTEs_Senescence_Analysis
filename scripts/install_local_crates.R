# Install the vendored upstream craTEs R sources as a local package namespace.
# Run from the repository root with:
#   Rscript scripts/install_local_crates.R

find_repo_dir <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (
      file.exists(file.path(current, "R", "getActivities.R")) &&
        file.exists(file.path(current, "LICENSES", "craTEs-MIT.md"))
    ) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the CraTEs_Senescence_Analysis repository root.", call. = FALSE)
    }
    current <- parent
  }
}

install_local_crates <- function(repo_dir = find_repo_dir()) {
  pkg_dir <- file.path(tempdir(), "craTEs_local_package")
  if (dir.exists(pkg_dir)) {
    unlink(pkg_dir, recursive = TRUE, force = TRUE)
  }

  dir.create(file.path(pkg_dir, "R"), recursive = TRUE, showWarnings = FALSE)

  r_files <- list.files(file.path(repo_dir, "R"), pattern = "[.]R$", full.names = TRUE)
  if (length(r_files) == 0L) {
    stop("No vendored craTEs R source files found in R/.", call. = FALSE)
  }
  file.copy(r_files, file.path(pkg_dir, "R"), overwrite = TRUE)

  writeLines(
    c(
      "Package: craTEs",
      "Title: Vendored craTEs Sources for CraTEs Senescence Analysis",
      "Version: 0.0.0.9000",
      "Authors@R: c(",
      "    person(\"Didier\", \"Trono\", role = c(\"fnd\")),",
      "    person(\"Cyril\", \"Pulver\", email = \"pulver.cyril@hotmail.com\", role = c(\"aut\", \"cre\")),",
      "    person(\"Vinzent\", \"Aschir\", role = c(\"ctb\"), comment = \"Analysis repository packaging only\")",
      "    )",
      "Description: Local installation wrapper for selected upstream craTEs source files bundled with the analysis repository.",
      "License: MIT + file LICENSE",
      "Encoding: UTF-8",
      "Imports:",
      "    data.table,",
      "    Biobase,",
      "    doParallel,",
      "    GenomicFeatures,",
      "    GenomicRanges,",
      "    glmnet,",
      "    purrr,",
      "    RColorBrewer,",
      "    ensembldb,",
      "    EnsDb.Hsapiens.v86"
    ),
    file.path(pkg_dir, "DESCRIPTION")
  )

  writeLines(
    c(
      "export(N_from_tsv)",
      "export(estimateRMSE)",
      "export(estimate_rmse_crossval)",
      "export(getActivities)",
      "export(getSignifActivities)",
      "export(get_pseudocounts)",
      "export(plotActivities)",
      "export(plot_signif_subfams)",
      "export(plot_signif_subfams_functional)",
      "export(preprocess_E_N_for_activities)",
      "export(pval_from_fstat)",
      "export(row_center)",
      "export(split_N_to_functional)"
    ),
    file.path(pkg_dir, "NAMESPACE")
  )

  file.copy(
    file.path(repo_dir, "LICENSES", "craTEs-MIT.md"),
    file.path(pkg_dir, "LICENSE"),
    overwrite = TRUE
  )

  utils::install.packages(pkg_dir, repos = NULL, type = "source")
}

if (sys.nframe() == 0L) {
  install_local_crates()
}
