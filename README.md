# CraTEs Senescence Analysis

Analysis notebooks and helper scripts for estimating transposable-element
subfamily activity in senescence and CRISPR perturbation datasets with craTEs.

## Repository Contents

- `notebooks/`: the three main executable analysis notebooks plus the mapping
  R notebook/script.
- `scripts/`: analysis, mapping, comparison, plotting, and export helpers.
- `local_craTEs_extension/`: local backward-compatible overrides used after
  loading craTEs.
- `R/`: selected upstream craTEs R source files included for reproducibility and
  attribution. These files are not authored by Vinzent Aschir.
- `data/`: data access instructions only. Raw input matrices are not committed.
- `results/`: compact output tables, diagnostics, and figures that are small
  enough for normal GitHub storage.
- `LICENSES/`: license text for third-party code included in this repository.

## External Code Attribution

The files in `R/` are selected source files from the original craTEs package by
Cyril Pulver, included from the local upstream checkout at commit
`d2923eb7df79a5eb1d343582739781ac9acbb70d`. They are included so the analysis
can be inspected and reproduced without misrepresenting package code as original
analysis code. The craTEs license text is stored in
`LICENSES/craTEs-MIT.md`.

Local analysis-specific changes and compatibility helpers live in
`local_craTEs_extension/` and `scripts/`.

## Data

Raw data files are intentionally not tracked in Git because of size and data
access constraints. Add the required files locally under `data/Hg19/` and
`data/Hg38/` before rerunning the notebooks. See `data/data_acces.md`.

## Environment

Restore the R environment from the lockfile:

```r
install.packages("renv")
renv::restore(lockfile = "renv.lock")
```

Install the Python/Jupyter environment when needed:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The notebooks currently call `library(craTEs)` and then source the local
overrides. If craTEs is not already installed in the restored R environment,
install the bundled upstream sources as a local package namespace:

```bash
Rscript scripts/install_local_crates.R
```

## Expected Workflow

1. Add the raw data files described in `data/data_acces.md`.
2. Restore R and Jupyter dependencies.
3. Run `notebooks/TE_Seq_craTEs.ipynb`.
4. Run `notebooks/CraTEs_LTR5_CRISPR.ipynb`.
5. Run `notebooks/result_comparison.ipynb`.
6. Use `notebooks/paper_mapping_craTEs.R` or `scripts/render_gene_mapping.R`
   when regenerating mapping diagnostics and reports.
