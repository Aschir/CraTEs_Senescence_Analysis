# Data Access and Local Layout

Raw input matrices are intentionally not committed to this repository. After a
fresh checkout, create the two data folders below and place the required files
under the exact local names used by the notebooks.

```text
data/
|-- Hg19/
|   |-- E_VPR-KRAB_CraTEs.txt.gz
|   |-- N_c_100_Hg19_craTEs.tsv
|   |-- N_w_2.5e5_Hg19_craTEs.tsv
|   `-- 2110_hg19_annotation.csv
`-- Hg38/
    |-- E_Senescene_TE_Seq_raw.tsv
    |-- GSE280365_sample_metadata.tsv
    |-- GSE280365_resultsdf.tsv
    |-- E_TE-Seq_prep_results.tsv
    |-- N_c_100_Hg38_craTEs.tsv
    `-- N_w_2.5e5_Hg38_craTEs.tsv
```

## Hg19: LTR5-Hs CRISPRa/CRISPRi Validation Data

These files are used mainly in `notebooks/CraTEs_LTR5_CRISPR.ipynb` and, for
the older linear reconstruction workflow, `notebooks/paper_mapping_craTEs.R`.

| Local file | Source and extraction | IDs / accessions | Notebook role |
|---|---|---|---|
| `data/Hg19/E_VPR-KRAB_CraTEs.txt.gz` | Download the processed GEO supplementary count table `GSE111332_gene_counts.txt.gz` and save it under this shorter analysis-specific name. | GEO RNA-seq subseries `GSE111332`; umbrella/superseries `GSE111337`; 12 samples `GSM3028661`-`GSM3028672`. | Expression matrix `E` for NCCIT LTR5-Hs CRISPR activation/repression. The notebook reads the table, uses `Geneid` as gene identifiers and drops `Geneid` plus `Length` before modelling. |
| `data/Hg19/N_c_100_Hg19_craTEs.tsv` | Original craTEs hg19 promoter TE-subfamily count matrix. In the upstream craTEs workflow this is also referred to as `N_prom_x_te_subfam.tsv`. | Zenodo DOI `10.5281/zenodo.6707955`. | Gene x TE-subfamily association matrix `N`, counted in 100 kb windows around promoters. This is the default hg19 validation matrix. |
| `data/Hg19/N_w_2.5e5_Hg19_craTEs.tsv` | Original craTEs hg19 weighted association matrix from the public craTEs data archive. | Zenodo DOI `10.5281/zenodo.8117257`. | Alternative distance-weighted hg19 `N` matrix for sensitivity checks. |
| `data/Hg19/2110_hg19_annotation.csv` | Gene annotation table from the original local craTEs paper/data bundle. The notebooks do not list a separate public accession for this helper table. | Local upstream annotation table; combined downstream with Ensembl BioMart GRCh37/Ensembl 80, `org.Hs.eg.db` and MyGene lookups. | Authoritative local HGNC-to-Ensembl helper used by `paper_mapping_craTEs.R` and mapping diagnostics. Expected columns include `ensembl_gene_id`, `entrezgene_id`, `hgnc_symbol` and `description`. |

Legacy names used by older scripts map to the cleaned repository names as
follows:

| Older/local name | Clean repository name |
|---|---|
| `data/GSE111332_gene_counts.txt.gz` | `data/Hg19/E_VPR-KRAB_CraTEs.txt.gz` |
| `data/N_prom_x_te_subfam.tsv` | `data/Hg19/N_c_100_Hg19_craTEs.tsv` |
| `data/2110_hg19_annotation.csv` | `data/Hg19/2110_hg19_annotation.csv` |

## Hg38: LF1 TE-Seq Senescence Data

These files are used mainly in `notebooks/TE_Seq_craTEs.ipynb` and compared in
`notebooks/result_comparison.ipynb`.

| Local file | Source and extraction | IDs / accessions | Notebook role |
|---|---|---|---|
| `data/Hg38/E_Senescene_TE_Seq_raw.tsv` | Download raw gene counts from ARCHS4 for the TE-Seq/GEO sample set. The notebook cites ARCHS4 and loads this file as the main raw expression matrix. | GEO series `GSE280365`; ARCHS4/GEO sample accessions `GSM8595549`-`GSM8595560`. | Main senescence expression matrix `E`. The notebook keeps columns matching `^GSM`, then renames samples with the GEO metadata cache. |
| `data/Hg38/GSE280365_sample_metadata.tsv` | Generated locally by `GEOquery::getGEO("GSE280365", GSEMatrix = TRUE, getGPL = FALSE)` and saved as a small cache with `geo_accession` and `title`. | GEO series `GSE280365`; same samples `GSM8595549`-`GSM8595560`. | Maps raw GSM columns to readable condition/replicate labels (`PRO`, `ESEN`, `LSEN`). This file can be regenerated from GEO. |
| `data/Hg38/GSE280365_resultsdf.tsv` | Download and decompress the GEO supplementary table `GSE280365_resultsdf.tsv.gz`. | GEO series `GSE280365`. | Paper/TE-Seq processed result table used for reference and raw-vs-normalized exploration. It is not the main raw-count input for craTEs. |
| `data/Hg38/E_TE-Seq_prep_results.tsv` | Same table content as `GSE280365_resultsdf.tsv`, kept under an older local analysis name. | GEO series `GSE280365`. | Optional legacy/local name referenced by exploratory or commented notebook paths. |
| `data/Hg38/N_c_100_Hg38_craTEs.tsv` | hg38-compatible craTEs promoter TE-subfamily count matrix. | Zenodo DOI `10.5281/zenodo.10390247`. | Gene x TE-subfamily association matrix `N`, counted in 100 kb promoter windows. This is the default hg38 senescence matrix. |
| `data/Hg38/N_w_2.5e5_Hg38_craTEs.tsv` | hg38-compatible craTEs distance-weighted association matrix. | Zenodo DOI `10.5281/zenodo.10390262`. | Alternative weighted hg38 `N` matrix used for sensitivity comparison against the counted matrix. |

The senescence sample metadata used in the notebook resolves to:

| GEO sample | Local condition label |
|---|---|
| `GSM8595549` | `PRO1`, LF1 cells, proliferating, replicate 1 |
| `GSM8595550` | `PRO2`, LF1 cells, proliferating, replicate 2 |
| `GSM8595551` | `PRO3`, LF1 cells, proliferating, replicate 3 |
| `GSM8595552` | `PRO4`, LF1 cells, proliferating, replicate 4 |
| `GSM8595553` | `ESEN1`, LF1 cells, early senescence, replicate 1 |
| `GSM8595554` | `ESEN2`, LF1 cells, early senescence, replicate 2 |
| `GSM8595555` | `ESEN3`, LF1 cells, early senescence, replicate 3 |
| `GSM8595556` | `ESEN4`, LF1 cells, early senescence, replicate 4 |
| `GSM8595557` | `LSEN1`, LF1 cells, late senescence, replicate 1 |
| `GSM8595558` | `LSEN2`, LF1 cells, late senescence, replicate 2 |
| `GSM8595559` | `LSEN3`, LF1 cells, late senescence, replicate 3 |
| `GSM8595560` | `LSEN4`, LF1 cells, late senescence, replicate 4 |

## Source Links

- NCBI GEO `GSE111332`: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111332>
- NCBI GEO `GSE111337`: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111337>
- NCBI GEO `GSE280365`: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280365>
- ARCHS4 data portal: <https://archs4.org/data>
- Zenodo `10.5281/zenodo.6707955`: <https://doi.org/10.5281/zenodo.6707955>
- Zenodo `10.5281/zenodo.8117257`: <https://doi.org/10.5281/zenodo.8117257>
- Zenodo `10.5281/zenodo.10390247`: <https://doi.org/10.5281/zenodo.10390247>
- Zenodo `10.5281/zenodo.10390262`: <https://doi.org/10.5281/zenodo.10390262>

## Notes For Reuse

- Keep raw data local unless redistribution is explicitly allowed by the
  original source.
- The large weighted `N_w_2.5e5_*` files and full
  `target_gene_matrix_*.tsv` outputs are not suitable for normal GitHub storage.
- `requirements.txt` belongs in the repository root, not inside `data/`,
  because it describes the Python/Jupyter environment rather than a dataset.
- Empty `data/Hg19/` and `data/Hg38/` folders are kept in Git with `.gitkeep`;
  all real data files under those folders are ignored by `.gitignore`.
