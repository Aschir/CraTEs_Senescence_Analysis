#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("scripts/render_gene_mapping.R", mustWork = TRUE)
}

repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(repo_dir)

local_r_cache <- file.path(repo_dir, ".cache", "R")
dir.create(local_r_cache, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(R_USER_CACHE_DIR = normalizePath(local_r_cache, mustWork = TRUE))

source("scripts/own_mapping_f.R")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render scripts/gene_mapping.Rmd.", call. = FALSE)
}

if (!configure_quarto_pandoc() || !rmarkdown::pandoc_available()) {
  stop("Pandoc was not found. Install pandoc or Quarto, then rerun the mapping rebuild.", call. = FALSE)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
out_path <- arg_or_null(args, "out")
if (is.null(out_path)) {
  e_gene_names_path <- require_arg(args, "e_gene_names", "--e-gene-names")
  n_gene_names_path <- require_arg(args, "n_gene_names", "--n-gene-names")
  out_path <- derive_output_path(e_gene_names_path, n_gene_names_path)
}
out_path <- absolute_path(out_path, must_work = FALSE)

report_path <- derive_report_path(out_path)
dir.create(dirname(report_path), showWarnings = FALSE, recursive = TRUE)

rmarkdown::render(
  input = "scripts/gene_mapping.Rmd",
  output_format = "html_document",
  output_file = basename(report_path),
  output_dir = dirname(report_path),
  quiet = TRUE
)

pdf_report_path <- derive_pdf_report_path(out_path)
if (latex_engine_available()) {
  rmarkdown::render(
    input = "scripts/gene_mapping.Rmd",
    output_format = "pdf_document",
    output_file = basename(pdf_report_path),
    output_dir = dirname(pdf_report_path),
    quiet = TRUE
  )
} else {
  message(
    "Skipping PDF mapping report because no LaTeX engine was found. ",
    "Install TinyTeX/MacTeX or put pdflatex, xelatex, or lualatex on PATH to create ",
    pdf_report_path,
    "."
  )
}
