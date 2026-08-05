<div style="border:1px solid #d8dee9; border-radius:8px; padding:16px 18px; background:#fbfcff; margin-bottom:14px;">
  <h1 style="margin:0 0 6px 0;">CraTEs Senescence Analysis</h1>
  <h5 style="margin:0 0 6px 0;">M1 Internship Project · 07.04.–07.08.2026 · IRCAN / Université Côte d'Azur</h5>
  <p style="margin:0; color:#586069;"><strong>Author:</strong> Vinzent Aimad Vivek Aschir &nbsp; | &nbsp; <strong>Main interface:</strong> R notebooks</p>
</div>

<details>
<summary><strong>Project metadata</strong></summary>

- Report title: `Inferring the Cis-Regulatory Influence of Transposable Elements on Gene Expression during Cellular Aging in Human`
- Internship report date: `31 July 2026`
- Supervisor: `Gaël Cristofari, IRCAN`
- GitHub repository: <https://github.com/Aschir/CraTEs_Senescence_Analysis>
- Keywords: `Senescence`, `craTEs`, `cis-regulatory activity`, `transposable elements`

</details>

## Pipeline overview

This repository is the code surface for my M1 internship project. The report explains the biological question and the interpretation in more depth; this README is meant to help you land in the code quickly, understand what each folder is for, and rerun the project without guessing which file belongs where.

The project reconstructs, generalises, and extends a craTEs-based workflow for estimating transposable-element (TE) subfamily activity from gene expression. In short:

> **E ≈ N × A**

`E` is the experiment-specific gene-expression matrix, `N` is the gene × TE-subfamily association matrix from the reference genome, and `A` is the inferred TE-subfamily activity matrix. The analysis uses this framework for two connected purposes:

- A CRISPRa/CRISPRi LTR5-Hs perturbation dataset is used as a reference system to check whether the reconstructed workflow recovers the expected VPR activation and KRAB repression signal as well as in the original craTEs analysis.
- A TE-Seq senescence dataset in human lung fibroblasts is used to analyse TE-associated cis-regulatory activity and its changes between proliferating cells (`PRO`), early senescence (`ESEN`) and late senescence (`LSEN`).

---

## 0 · Quick start

### ── Fresh checkout ──

```bash
git clone https://github.com/Aschir/CraTEs_Senescence_Analysis.git
cd CraTEs_Senescence_Analysis
```

### ── Restore the R environment ──

The R package snapshot is pinned in `renv.lock`. The lockfile records R `4.5.3` and Bioconductor `3.21`.

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::activate(); renv::restore(lockfile = "renv.lock", exclude = "craTEs")'
```

The clean GitHub repository carries selected upstream craTEs source files in `R/`. Install them as a local `library(craTEs)` namespace after restoring the other R packages:

```bash
Rscript -e 'renv::load("."); source("scripts/install_local_crates.R"); install_local_crates()'
```

**Package-loading note.** R packages are restored from `renv.lock`; Python/Jupyter packages are restored from `requirements.txt`. The local `craTEs` install must happen after `renv::restore()` so `library(craTEs)` in the notebooks resolves to the vendored package namespace built from this repository.

### ── Restore Jupyter/Python ──

The notebooks use an R kernel, but Jupyter itself is managed with Python. The tested Python/Jupyter environment is pinned in `requirements.txt` and was tested with Python `3.13.1`.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
Rscript -e 'renv::load("."); IRkernel::installspec(name = "ir-crates-senescence", displayname = "R (craTEs senescence)", user = TRUE)'
jupyter lab
```

### ── Add the input data ──

Raw matrices are intentionally not committed. Place them under `data/Hg19/` and `data/Hg38/` exactly as described in [data/README.md](data/README.md).

### ── Run the notebooks ──

Recommended order:

1. `notebooks/CraTEs_LTR5_CRISPR.ipynb`
2. `notebooks/TE_Seq_craTEs.ipynb`
3. `notebooks/result_comparison.ipynb`

Run `notebooks/paper_mapping_craTEs.R` only when you want the older, more linear mapping/paper-reconstruction script for reference.

---

## 1 · Repository structure

### ── What lives where ──

| Path | Role |
|---|---|
| `notebooks/` | Main executable analysis surface. The notebooks are structured as numbered sections with setup, data loading, mapping, preprocessing, activity estimation, contrasts, target matrices and comparison. |
| `scripts/` | Reusable project helpers for mapping, plotting, comparison, trend classification and exports. |
| `local_craTEs_extension/` | Local drop-in overrides sourced after `library(craTEs)`. These keep the craTEs function names but add hg38 support, stricter checks and clearer diagnostics. |
| `R/` | Selected upstream craTEs R source files bundled for reproducibility and attribution. |
| `data/` | Data access notes and local data slots. Raw input files are ignored by Git. |
| `results/` | Compact result tables, mapping provenance, diagnostics and figures that can be inspected without rerunning everything. Large target-gene matrices are ignored. |
| `requirements.txt` | Python/Jupyter package snapshot. |
| `renv.lock` | R/Bioconductor package snapshot. |
| `LICENSES/` | Third-party license text for bundled upstream craTEs code. |

### ── Main notebooks ──

| Notebook | When to open it |
|---|---|
| `CraTEs_LTR5_CRISPR.ipynb` | Validation workflow. It uses the LTR5-Hs CRISPRa/CRISPRi dataset and checks whether VPR shifts activity upward and KRAB shifts activity downward. |
| `TE_Seq_craTEs.ipynb` | Main senescence workflow. It estimates activity scores and contrasts for `PRO`, `ESEN` and `LSEN`. |
| `result_comparison.ipynb` | Cross-configuration comparison. It compares result directories, N-matrix choices, activity-score similarity, trend classes, p-value evidence and target matrices. |

---

## 2 · Data access

The canonical data documentation is [data/README.md](data/README.md).
It lists every required local file, the source database, GEO/Zenodo accession or DOI, the extraction route used in the notebooks, and the local filename expected by the analysis. Raw matrices stay local and are ignored by Git, while the empty `data/Hg19/` and `data/Hg38/` folders are kept as slots with `.gitkeep` files.

---

## 3 · Local package versions

### ── R snapshot ──

The full R environment is locked in `renv.lock`.

| Layer | Version / source |
|---|---|
| R | `4.5.3` |
| Bioconductor | `3.21` |
| `craTEs` | `0.0.0.9000`, local source namespace |
| Upstream craTEs source commit | `d2923eb7df79a5eb1d343582739781ac9acbb70d` |
| Key modelling/data packages | `glmnet 5.0`, `data.table 1.18.4`, `dplyr 1.2.1`, `ggplot2 4.0.3`, `biomaRt 2.64.0`, `EnsDb.Hsapiens.v86 2.99.0`, `GEOquery 2.76.0` |
| Notebook/report packages | `IRkernel 1.3.2`, `rmarkdown 2.31`, `ggtext 0.1.2`, `ggVennDiagram 1.5.7` |

### ── Python/Jupyter snapshot ──

The full Python environment is locked in `requirements.txt`.

| Layer | Version |
|---|---|
| Python | tested with `3.13.1` |
| Jupyter | `jupyter 1.1.1` |
| JupyterLab | `jupyterlab 4.5.7` |
| Notebook | `notebook 7.5.6` |
| Kernel bridge | `ipykernel 7.2.0` |

**Environment note.** If `renv::restore()` complains about the original local craTEs path from the lockfile, use the `exclude = "craTEs"` restore command above and then run the local install command from the quick start. That command rebuilds the local namespace from the bundled `R/` files inside the active `renv` library.

---

## 4 · Mapping and preprocessing

### ── Why mapping exists ──

The expression matrices and the `N` matrices do not always speak the same
identifier language. `E` often starts with HGNC symbols, while `N` is indexed by
Ensembl gene IDs. The report separates preprocessing into two layers: first the
E-matrix-specific identifier cleanup, then the joint `E`/`N` preprocessing used
for activity estimation.

### ── Report-level preprocessing order ──

The full preprocessing order from the report is:

1. Load the raw expression count matrix and keep the sample/condition labels
   attached to the active analysis.
2. Correct Excel-converted HGNC symbols, for example date-like `MARCH` or
   `SEPT` gene names, before database mapping.
3. Map HGNC symbols and already-present Ensembl IDs to standard Ensembl gene
   IDs and gene biotypes. The final mapping route uses the genome-matched
   BioMart release first, then `org.Hs.eg.db`, then MyGene for unresolved cases.
4. Resolve duplicated gene identifiers that were already present in the raw
   matrix or produced by the Ensembl mapping. The final analyses use
   `handling <- "discard"`, the conservative craTEs-style strategy.
5. Keep protein-coding genes using the mapped biotype annotations.
6. Exclude directly manipulated genes where the experimental design requires
   it, so the perturbation itself does not dominate the inferred TE activity.
7. Restrict `E` and `N` to the shared Ensembl gene IDs and arrange their rows in
   the same order.
8. Remove genes whose expression is below `10` in every sample.
9. Remove TE-subfamily predictor columns in `N` with total association count
   below `150`.
10. Remove genes without any remaining TE-subfamily association.
11. Add a sample-wise pseudocount computed from the `q = 0.05` non-zero count
    quantile.
12. Convert expression counts to TPM using reference-specific exonic gene
    lengths.
13. Apply the `log2` transformation.
14. Center the expression matrix column-wise and then row-wise.

The duplicate-handling strategies (`discard`, `max_sum`, `max_variance`, `sum`,
`average`, `identical`) were compared as diagnostics; the reproducible analysis
uses `discard`.

### ── Where this happens in code ──

The identifier cleanup, mapping, duplicate handling and protein-coding filter
are run in the notebooks with helpers from `scripts/`. The joint preprocessing
is implemented in `local_craTEs_extension/preprocess_E_N_for_activities_hg38.R`,
which generalises the original craTEs preprocessing to hg38 while preserving the
hg19 route.

Inside that function, the E-only low-expression filter is applied immediately
before the `E`/`N` intersection, then the rows are aligned into the common gene
space. After TPM conversion, `N` is aligned once more to the genes retained by
the exon-length annotation. The returned object contains `E_centered`, filtered
`N`, and removed TE predictors (`lost_cols`) for diagnostics. Exon-length caches
are stored under `data/<Hg19|Hg38>/temp/` and are regenerated when absent.

### ── Mapping provenance files ──

The GitHub export intentionally keeps the compact mapping files needed to make
the notebooks reproducible without rebuilding BioMart/MyGene mapping on every
fresh checkout:

```text
results/<result_dir>/mapping/E_*_gene_names.tsv
results/<result_dir>/mapping/N_*_gene_names.tsv
results/<result_dir>/mapping/map_<E_name>_<N_name>.tsv
results/<result_dir>/mapping/report_<E_name>_<N_name>.html
results/<result_dir>/mapping/diagnostics/*.tsv
```

Normally you can run from the cached `map_*.tsv` file. Rebuild the mapping only
when the raw `E` matrix, selected `N` matrix, genome build or mapping strategy
changes.

For command-line rebuilds, the same logic is available through:

```bash
Rscript --vanilla scripts/render_gene_mapping.R \
  --e-gene-names results/<result_dir>/mapping/<E_name>_gene_names.tsv \
  --n-gene-names results/<result_dir>/mapping/<N_name>_gene_names.tsv \
  --out results/<result_dir>/mapping/map_<E_name>_<N_name>.tsv \
  --primary-mart current \
  --mirror www
```

Use `--primary-mart grch37` for hg19/GRCh37-style mapping and
`--primary-mart current` for hg38/current analyses. If the default Ensembl
mirror is slow, try `useast`, `uswest` or `asia`.

---

## 5 · Analysis workflow

### ── CRISPR validation workflow ──

Open `notebooks/CraTEs_LTR5_CRISPR.ipynb`.

Main idea: LTR5-Hs was targeted with two opposite dCas9-based systems. VPR is expected to increase regulatory activity; KRAB is expected to repress it. Recovering this bidirectional pattern is a qualitative validation that the reconstructed workflow still behaves like the original craTEs logic.

Typical output directory:

```text
results/results_E_VPR-KRAB_CraTEs_N_c_100_Hg19_craTEs/
```

### ── TE-Seq senescence workflow ──

Open `notebooks/TE_Seq_craTEs.ipynb`.

The user-editable settings are near the top:

```r
handling <- "discard"
reference_genome <- "Hg38"       # "Hg19" also supported
n_matrix_key <- "counted_100kb"  # or "weighted_250kb"
analysis_mode <- NULL
```

The biological comparison is built around:

- `PRO`: proliferating cells
- `ESEN`: early senescence
- `LSEN`: late senescence

The main pairwise contrasts are:

- `ESEN_vs_PRO`
- `LSEN_vs_PRO`
- `LSEN_vs_ESEN`

### ── Result comparison workflow ──

Open `notebooks/result_comparison.ipynb`.

This notebook does not rerun craTEs. It reads stored result directories and compares:

- shared samples, genes and TE subfamilies,
- model-fit diagnostics,
- activity-score correlations,
- active-signal Jaccard overlap,
- trend-class agreement,
- contrast evidence,
- target-gene matrix similarity.

The main editable block is:

```r
result_dirs <- c(
  "results/results_E_Senescene_TE_Seq_raw_N_c_100_Hg38_craTEs",
  "results/results_E_Senescene_TE_Seq_raw_N_w_2.5e5_Hg38_craTEs"
)

te_global_top_k <- 10L
te_ranking_mode <- "activity_change"
te_ranking_metric <- "sum_step_abs_delta"
```

### ── Reusing the pipeline for new data ──

For a new expression dataset, the main practical template should be `notebooks/TE_Seq_craTEs.ipynb`. It is the most complete notebook for a matrix-pair-specific analysis because it keeps dataset setup, condition definition, identifier mapping, preprocessing diagnostics, activity estimation, contrast testing, target-gene summaries and exports in one coherent workflow.

Do not reuse the senescence conditions mechanically. Define the conditions, control/treatment groups and pairwise contrasts again near the top of the notebook according to the new experimental design. The raw expression import can be adapted to the new file format, but the downstream structure should stay close to the TE-Seq workflow so each `E`/`N` matrix pair remains interpretable and comparable.

The preprocessing diagnostics should also be read again in a dataset-dependent way. In particular, inspect the Excel-symbol repair, mapping recovery, duplicate handling, protein-coding filter, `E`/`N` gene overlap, low-expression gene filtering, TE-subfamily association threshold, final matrix dimensions, model-fit summaries and activity diagnostics. These checks decide whether a specific matrix pair is technically usable and biologically meaningful.

If several genome builds or `N` matrix settings are plausible, run each reasonable matrix pair separately and then compare them with `notebooks/result_comparison.ipynb`. Stable directions, trend classes and target summaries across matrix pairs are more convincing than a single isolated result. The CRISPR notebook is mainly a validation/reference workflow; for new biological datasets, the TE-Seq notebook should be the primary guide.

---

## 6 · Results

### ── What each result directory contains ──

| File | Meaning |
|---|---|
| `activity_scores.tsv` | TE-subfamily × sample activity matrix `A`. |
| `sample_activity_diagnostics.tsv` | Sample-level model fit, sparsity and lambda diagnostics. |
| `activity_fit_summary.tsv` | Overall and condition-wise model-fit summary. |
| `activity_lambda_diagnostics.tsv` | Lambda grid and selected regularisation information when available. |
| `*_signif_act.tsv` | Pairwise significant activity table with coefficients and p-values. |
| `signif_model_diagnostics.tsv` | Full unpenalised contrast-model diagnostics. |
| `te_contrast_summary.tsv` | All-TE contrast summary with mean activities, deltas, adjusted p-values and trend classes. |
| `target_summary.tsv` | Compact target-gene summary for selected/high-ranking TE subfamilies. |
| `target_gene_matrix_*.tsv` | Full gene × TE target-dependency matrices. These are large and ignored in the clean GitHub export. |
| `figures/` | Exported heatmaps and selected activity-change figures. |
| `mapping/E_*_gene_names.tsv` and `mapping/N_*_gene_names.tsv` | Gene-name snapshots used to build each mapping cache. |
| `mapping/map_*.tsv` | Cached identifier mapping used by the notebooks. |
| `mapping/report_*.html` | Rendered provenance report for the mapping run. |
| `mapping/diagnostics/` | Mapping source counts, unmapped genes and biotype summaries. |

### ── Biological reading guide ──

- Higher activity in `ESEN` than `PRO` suggests a TE subfamily is more strongly associated with relative gene expression during senescence entry.
- Similar `ESEN` and `LSEN` profiles suggest an activity pattern established early and maintained later.
- A strong target-gene dependency score means the TE subfamily contributed to the model prediction for that gene; it does not prove direct regulation by a specific TE locus.
- Counted and distance-weighted `N` matrices are sensitivity checks. Stable direction, class and ranking are more interpretable across configurations than the absolute coefficient size.

---

## 7 · Attribution and licenses

The files in `R/` are selected source files from the original craTEs package by Cyril Pulver

The analysis can be inspected and reproduced without presenting upstream package code as original analysis code. The craTEs MIT license text is stored in `LICENSES/craTEs-MIT.md`.

Project-specific analysis helpers live in `scripts/` and local compatibility overrides live in `local_craTEs_extension/`.

---

## 8 · Troubleshooting

### ── `library(craTEs)` fails ──

Run:

```bash
Rscript -e 'renv::load("."); source("scripts/install_local_crates.R"); install_local_crates()'
```

Then restart the R kernel.

### ── Mapping cache is missing ──

This should not happen in the clean export, because the compact `map_*.tsv`
files are intentionally tracked. If you deleted them locally or changed the
input matrices, run the notebook until the mapping input paths are printed, then
execute the optional mapping rebuild cell. The rebuild needs internet access to
BioMart, `org.Hs.eg.db` and MyGene.

### ── BioMart mirror fails ──

Retry the mapping rebuild with another mirror:

```text
www, useast, uswest, asia
```

For hg19 mappings, use `--primary-mart grch37`. For hg38 mappings, use `--primary-mart current`.

---

## 9 · Questions / credentials

For questions about the analysis, please use the GitHub repository or contact the project author:

- Author: Vinzent Aimad Vivek Aschir
- GitHub: <https://github.com/Aschir/CraTEs_Senescence_Analysis>
- Project context: M1 Data Science and Artificial Intelligence internship at IRCAN / Université Côte d'Azur
- Supervisor: Gaël Cristofari, IRCAN
