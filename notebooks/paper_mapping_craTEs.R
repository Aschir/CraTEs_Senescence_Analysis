# version paper mapping
repo_dir <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "scripts", "craTEs_own_function.R"))) {
  stop("Run this script from the craTEs repository root.", call. = FALSE)
}

## libraries
library(Matrix)
library(data.table)
library(craTEs)
ls("package:craTEs")
library(parallel)
source("scripts/craTEs_own_function.R") 
library(ggplot2)
library(dplyr)




handling = "discard"

## E - Expression matrix for NCCIT cells with CRISPRa/i against LTR5-Hs, LTR5A, LTR5B.
E <- fread("./data/GSE111332_gene_counts.txt.gz")
gene_ids <- E$Geneid
E <- as.matrix(E[, -c(1,2), with = FALSE])  # drop Geneid and Length
rownames(E) <- gene_ids


duplicate_stats(rownames(E),E)
### critical: how to deal with duplicate gene names
E <-  resolve_duplicates(E, method = handling, threshold = 10)  

#handle excel gene naming problem 
gene_ids <- transform_excel_dates(rownames(E)) 
rownames(E) <- gene_ids
duplicate_stats(rownames(E),E)



## N - gene / TE subfamily - Count Association matrix
N <- fread("./data/N_prom_x_te_subfam.tsv", header = TRUE, sep = "\t")
N <- as.matrix(N)
rownames(N) <- N[,1]  # first col = promoters
genes_N  <- rownames(N)
N <- N[,-1] 
storage.mode(N) <- "double" 

#duplicated genes: 
N_dub_stat = duplicate_stats(rownames(N))
sprintf("We have %d redundant duplicated Gene-rows in N in duplicated groups of %d ", N_dub_stat$redundant_rows, N_dub_stat$duplicate_groups)

# how dense is our N matrix? 
sprintf("the amount of 0 vals in the N matrix is %.2f%%", 100 * mean(N == 0))

# what which te fams would be discarded in prerpocessing: 
check_npm_threshold(N, threshold = 150)


########## Mapping from paper CSV + biotype annotation #################

library(biomaRt)
library(mygene)

# ── 1. Load paper mapping CSV ─────────────────────────────
mapping_df <- read.csv(
  "./data/2110_hg19_annotation.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

mapping_df <- mapping_df %>%
  mutate(
    hgnc_symbol = trimws(hgnc_symbol),
    ensembl_gene_id = trimws(ensembl_gene_id)
  ) %>%
  filter(
    hgnc_symbol != "",
    ensembl_gene_id != "",
    !is.na(hgnc_symbol),
    !is.na(ensembl_gene_id)
  ) %>%
  distinct(hgnc_symbol, ensembl_gene_id, .keep_all = TRUE)

# ── 2. Use CSV as authoritative gene-symbol → Ensembl mapping ────────
map_csv <- smart_dedup_mapping(
  map_df = mapping_df %>%
    select(hgnc_symbol, ensembl_gene_id),
  already_mapped = character(0),
  N_genes = genes_N
)

mapped_ids <- map_csv[gene_ids]
unmapped_final <- gene_ids[is.na(mapped_ids)]

sprintf(
  "%d of %d genes (%.2f%%) could not be mapped by CSV",
  length(unmapped_final),
  length(gene_ids),
  100 * length(unmapped_final) / length(gene_ids)
)

display_intersect(unname(map_csv[!is.na(map_csv)]), genes_N)
duplicate_stats(unname(map_csv[!is.na(map_csv)]))

# ── 3. Get biotypes for the CSV Ensembl IDs ──────────────────────────
ensembl_ids_csv <- unique(unname(map_csv[!is.na(map_csv)]))

mart_1 <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 80,
  mirror = "www"
)

biotype_biomart <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_csv,
  mart = mart_1
) %>%
  filter(!is.na(gene_biotype), gene_biotype != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ── 4. Fallback biotypes with mygene for missing Ensembl IDs ─────────
missing_biotype_ids <- setdiff(
  ensembl_ids_csv,
  biotype_biomart$ensembl_gene_id
)

mygene_results <- queryMany(
  missing_biotype_ids,
  scopes = "ensembl.gene",
  fields = c("ensembl.gene", "type_of_gene"),
  species = "human"
)

mygene_df <- as.data.frame(mygene_results)

biotype_mygene <- mygene_df %>%
  transmute(
    ensembl_gene_id = query,
    gene_biotype = gsub("-", "_", as.character(type_of_gene))
  ) %>%
  filter(
    !is.na(ensembl_gene_id),
    !is.na(gene_biotype),
    gene_biotype != ""
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ── 5. Combined biotype table ─────────────────────────────
biotype_combined <- bind_rows(
  biotype_biomart,
  biotype_mygene %>%
    filter(!ensembl_gene_id %in% biotype_biomart$ensembl_gene_id)
)

sprintf(
  "%d of %d mapped Ensembl IDs still have no biotype",
  length(setdiff(ensembl_ids_csv, biotype_combined$ensembl_gene_id)),
  length(ensembl_ids_csv)
)

plot_mapping_distribution(
  biotype_combined,
  category_col = "gene_biotype",
  highlight = "protein_coding",
  heading = "Biotype distribution — CSV mapping with database annotation"
)

# ── 6. Reverse mapping for later target-gene annotation ──────────────
ensembl_to_symbol <- setNames(names(map_csv), unname(map_csv))

# ── 7. Apply CSV mapping to E ─────────────────────────────
keep <- !is.na(mapped_ids)

E_mapped <- E[keep, , drop = FALSE]
rownames(E_mapped) <- unname(mapped_ids[keep])

duplicate_stats(rownames(E_mapped), E_mapped)

E_resolved <- resolve_duplicates(
  E_mapped,
  method = handling,
  threshold = 10
)

duplicate_stats(rownames(E_resolved), E_resolved)

# ── 8. Optional: keep only protein-coding genes ──────────────────────
protein_coding_ensembl <- biotype_combined %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(ensembl_gene_id)

E_filtered <- E_resolved[
  rownames(E_resolved) %in% protein_coding_ensembl,
  ,
  drop = FALSE
]

sprintf("E before biotype filter: %d genes", nrow(E_resolved))
sprintf("E after  biotype filter: %d genes", nrow(E_filtered))

plot_mapping_distribution(
  biotype_combined %>%
    filter(ensembl_gene_id %in% rownames(E_resolved)),
  category_col = "gene_biotype",
  highlight = "protein_coding",
  heading = "Biotype distribution — E before protein_coding filter"
)


res <- craTEs::preprocess_E_N_for_activities(
  count_table = E_filtered,
  N = N,
  log_tpm_plot_path = "log_tpm.pdf"
)


E_centered <- res$E_centered
N_filtered  <- res$N



# ── Common Gene Display ─────────────────────────────────────────
# change here for E_centered  N_filtered to see after the preprocessing function


display_intersect(rownames(E_centered), rownames(N_filtered))

duplicate_stats(rownames(E_centered), E_centered) # E_mapped)


# ── Unmapped Genes ─────────────────────────────────────────
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

classes <- classify_unmapped(unmapped_final)
cat("\nUnmapped gene classes:\n")
sapply(classes, length)

# ── CraTEs Activity ─────────────────────────────────────────

# craTEs::getActivities
act <- getActivities(E_centered , N_filtered)
A <- act$activities # TE activity matrix
A <- A[rownames(A) != "(Intercept)", , drop = FALSE]
colnames(A) <- colnames(E_centered)

# Elastic net - penalised linear regression. (Activity coefficients close to 0) 
boxplot(as.matrix(A), main = "Activity distributions per sample")


# ── 2. Model fit diagnostics (paper Suppl.) ───────────────
# how much varinace in gene expression (E) is explained by TE activity (N). 
cat("\n── Model fit ───────────────────────────────\n")
cat(sprintf("Mean R²:              %.3f\n", mean(act$r2)))
cat(sprintf("Mean R²/R²_max:       %.3f  (expected: 0.6-0.9)\n", mean(act$r2 / act$r2_max)))

# R² per condition
cat(sprintf("Mean R² VPR samples:  %.3f\n",
            mean(act$r2[grepl("VPR",  names(act$r2))])))
cat(sprintf("Mean R² KRAB samples: %.3f\n",
            mean(act$r2[grepl("KRAB", names(act$r2))])))

cat(sprintf("Sparsity of A:        %.2f%% zeros\n", 100 * mean(act$activities == 0)))

# act$r2$lambdas_best is the best regularisation coefficients 
# ── 3. Define groups ──────────────────────────────────────
treatment_VPR  <- grep("LTR5HS.*VPR",      colnames(E_centered), value = TRUE)
control_VPR    <- grep("nontarget.*VPR",    colnames(E_centered), value = TRUE)
treatment_KRAB <- grep("LTR5HS.*KRAB",     colnames(E_centered), value = TRUE)
control_KRAB   <- grep("nontarget.*KRAB",  colnames(E_centered), value = TRUE)

# ── 4. Delta activity for targeted TEs (paper Fig 1C) ────
# difference in activity when we activate or repress compared to normal. 
cat("\n── Delta activity (treatment - control) ────\n")
tes <- c("LTR5_Hs", "LTR5A", "LTR5B")
for (te in tes) {
  delta_vpr  <- mean(A[te, treatment_VPR])  - mean(A[te, control_VPR])
  delta_krab <- mean(A[te, treatment_KRAB]) - mean(A[te, control_KRAB])
  cat(sprintf("%-10s  VPR (activation):  %+.4f\n", te, delta_vpr))
  cat(sprintf("%-10s  KRAB (repression): %+.4f\n", te, delta_krab))
}
# Expected: LTR5_Hs positive for VPR, negative for KRAB

# ── 5. Quality Control - Replicate correlation ─────────────
# Do the experiments in same condition conclude the same activity? 
cat("\n── Replicate correlations ──────────────────\n")
# sp + standart cas9 and sa targetting cas9. 
cat(sprintf("Sp vs Sa KRAB rep1-rep2: %.3f\n", cor(A[, "LTR5HS-Sp-KRAB-rep1"], A[, "LTR5HS-Sa-KRAB-rep1"]))) # different experiments dont conclude same activity. 
cat(sprintf("Sp KRAB rep1 vs rep2:    %.3f\n", cor(A[, "LTR5HS-Sp-KRAB-rep1"], A[, "LTR5HS-Sp-KRAB-rep2"]))) # different replicates same experiment Krab high correlation. 
cat(sprintf("Sp VPR rep1 vs rep2:     %.3f\n", cor(A[, "LTR5HS-Sp-VPR-rep1"],  A[, "LTR5HS-Sp-VPR-rep2"]))) # different replicates same experiment VPR high correlation.

# ── 6. getSignifActivities ────────────────────────────────

# craTEs::getSignifActivities
#shift in gene expression when we activate or repress compared to normal.
signif_VPR  <- getSignifActivities(E_centered, N_filtered, treatment_group = treatment_VPR, control_group   = control_VPR)

signif_KRAB <- getSignifActivities(E_centered, N_filtered, treatment_group = treatment_KRAB, control_group   = control_KRAB)

preprocessed <- list( E_centered = E_centered, N = N_filtered)

# ── 7. PRIMARY PLOT: Fig 1C equivalent ───────────────────

# black dots: significant TE subfamilies (p-val < alpha = 0.05)
# line standard error of these activity 
# X-Axis: Effect-Size (mean difference in activity between treatment and control) as log2 fold change single TE insertion contributes to a gene.
# Y-Axis: -log10(adj. p-value)
par(mfrow = c(1, 1))

plot_signif_subfams(
  res     = signif_VPR,
  alpha   = 0.05,
  main    = "CRISPRa (VPR) — LTR5-Hs activation",
  n_label = 4,
  n_subplots = 2, idx = 1)

plot_signif_subfams(
  res     = signif_KRAB,
  alpha   = 0.05,
  main    = "CRISPRi (KRAB) — LTR5-Hs repression",
  n_label = 3,
  n_subplots = 2, idx = 2)

par(mfrow = c(1, 1))  

# ── 8. Activity boxplots per targeted TE ─────────────────
# Plots: specific TE sequences are active or repressed or normal.
# X-Axsis: all samples grouped by sample group (treatment vs. control)
# one point in boxplot: activity coefficient of all genes close to this specific TE per sample-group for the specific TE
# Y-Axsis: regulatory coefficients of target TE for samples within one sample-group (up or down regulated,tar vs ctrl)
par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
for (te in tes) {
  # CRISPRa data
  data_vpr  <- list(A[te, control_VPR], A[te, treatment_VPR])
  names_vpr <- c(sprintf("control\n(%d samples)", length(control_VPR)),
                 sprintf("treatment\n(%d samples)", length(treatment_VPR)))
  
  boxplot(data_vpr,
          names = names_vpr, # Hier werden die neuen Namen gesetzt
          main  = sprintf("%s — CRISPRa", te),
          ylab  = "Activity",
          col   = c("grey80", "steelblue"))
  
  abline(h = 0, lty = 2, col = "red")
  
  # 2. CRISPRi Daten vorbereiten
  data_krab  <- list(A[te, control_KRAB], A[te, treatment_KRAB])
  names_krab <- c(sprintf("control\n(%d samples)", length(control_KRAB)),
                  sprintf("treatment\n(%d samples)", length(treatment_KRAB)))
  
  boxplot(data_krab,
          names = names_krab,
          main  = sprintf("%s — CRISPRi", te),
          ylab  = "Activity",
          col   = c("grey80", "tomato"))
  
  abline(h = 0, lty = 2, col = "red")}

mtext("Inferred cis-regulatory activity of targeted TE subfamilies (LTR5_Hs, LTR5A, LTR5B)\nacross all proximal genes (≤100kb) — CRISPRa (VPR) vs CRISPRi (KRAB)", 
      outer = TRUE, cex = 1.5, font = 2, line = 1)
par(mfrow = c(1,1), oma = c(0,0,0,0))  # reset

# ── 9. Activity across all samples ordered by group ───────
# for each TE subfamily, plot the activity across all samples
# X-Axs
metadata <- data.frame(
  sample       = colnames(A),
  perturbation = ifelse(grepl("nontarget", colnames(A)), "control", "target"),
  system       = ifelse(grepl("VPR", colnames(A)), "VPR", "KRAB"),
  row.names    = colnames(A),
  stringsAsFactors = FALSE)

metadata$group <- factor(paste(metadata$perturbation, metadata$system, sep = "_"))

par(mfrow = c(1, 1))
ord <- order(metadata$group)
cols <- c("control_KRAB" = "tomato",
          "control_VPR"  = "steelblue",
          "target_KRAB"  = "darkred",
          "target_VPR"   = "darkblue")

for (te in tes) {
  plot(
    seq_along(ord),
    A[te, ord],
    col  = cols[as.character(metadata$group[ord])],
    pch  = 16,
    main = sprintf("%s activity across samples", te),
    xlab = "Samples (ordered by group)",
    ylab = "Activity")
  
  abline(h = 0, lty = 2, col = "grey50")
  legend("topright", names(cols), col = cols, pch = 16, cex = 0.8)}

# ── 10. Top TE subfamilies by mean absolute activity ──────
te_importance <- sort(rowMeans(abs(A), na.rm = TRUE), decreasing = TRUE)
top_te <- data.frame(
  TE_family  = names(te_importance),
  importance = te_importance
)
cat("\n── Top 10 TE subfamilies by mean |activity| ─\n")
print(head(top_te, 10) |> dplyr::mutate(across(where(is.numeric), \(x) round(x, 4))))

# Are LTR5_Hs, LTR5A, LTR5B in the top?
cat(sprintf("\nRank of LTR5_Hs: %d of %d\n",
            which(top_te$TE_family == "LTR5_Hs"), nrow(top_te)))
cat(sprintf("Rank of LTR5A:   %d of %d\n",
            which(top_te$TE_family == "LTR5A"),   nrow(top_te)))
cat(sprintf("Rank of LTR5B:   %d of %d\n",
            which(top_te$TE_family == "LTR5B"),   nrow(top_te)))


# ── 11. Target gene matrix ────────────────────────────────
# Which genes are most dependent on each TE-subfamily signal?
# Run separately for VPR (activation) and KRAB (repression)

# Computing target gene matrix (VPR) 
target_genes_VPR  <- getTargetMatrix(signif_VPR,  preprocessed)

# Computing target gene matrix (KRAB)
target_genes_KRAB <- getTargetMatrix(signif_KRAB, preprocessed)

# Top genes most impacted by LTR5_Hs perturbation
top_VPR <- sort(target_genes_VPR[, "LTR5_Hs"],  decreasing = TRUE)
top_KRAB <- sort(target_genes_KRAB[, "LTR5_Hs"], decreasing = TRUE)

cat("\nTop 10 genes most dependent on LTR5_Hs (VPR):\n")
print(head(top_VPR, 10))

cat("\nTop 10 genes most dependent on LTR5_Hs (KRAB):\n")
print(head(top_KRAB, 10))

# Genes appearing in BOTH conditions = strongest candidates
shared <- intersect(names(top_VPR)[1:50], names(top_KRAB)[1:50])
cat("Genes regulated bidirectionally by LTR5_Hs:\n")
print(shared)


# ── 11. Target gene matrix ────────────────────────────────
# getTargetMatrix() estimates which genes depend most strongly on each TE subfamily 
# by computationally silencing TE subfamilies in the N matrix and measuring the loss in prediction accuracy.

# Computing target gene matrix (VPR)
target_genes_VPR  <- getTargetMatrix(signif_VPR,  preprocessed)

# Computing target gene matrix (KRAB)
target_genes_KRAB <- getTargetMatrix(signif_KRAB, preprocessed)

# Top genes most impacted by LTR5_Hs perturbation
top_VPR  <- sort(target_genes_VPR[, "LTR5_Hs"],  decreasing = TRUE)
top_KRAB <- sort(target_genes_KRAB[, "LTR5_Hs"], decreasing = TRUE)

top_VPR_df  <- annotate_targets(head(top_VPR, 10), ensembl_to_symbol)
top_KRAB_df <- annotate_targets(head(top_KRAB, 10), ensembl_to_symbol)

cat("\nTop 10 genes most dependent on LTR5_Hs (VPR):\n")
print(top_VPR_df)

cat("\nTop 10 genes most dependent on LTR5_Hs (KRAB):\n")
print(top_KRAB_df)

# Shared genes
shared <- intersect(
  names(top_VPR)[1:50],
  names(top_KRAB)[1:50]
)

shared_df <- data.frame(
  ensembl_gene_id = shared,
  hgnc_symbol = ensembl_to_symbol[shared],
  stringsAsFactors = FALSE
)
