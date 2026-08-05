# Local backward-compatible override for craTEs::getSignifActivities().
# Source after library(craTEs), then call getSignifActivities() without namespace.

.get_signif_delta_design <- function(E_centered,
                                     N,
                                     treatment_group,
                                     control_group,
                                     contrast_mode = c("auto", "paired", "group_mean", "all_pairs")) {
  contrast_mode <- match.arg(contrast_mode)

  if (identical(contrast_mode, "auto")) {
    contrast_mode <- if (length(treatment_group) == length(control_group)) {
      "paired"
    } else {
      "group_mean"
    }
  }

  if (identical(contrast_mode, "paired")) {
    if (length(treatment_group) != length(control_group)) {
      stop(
        "contrast_mode = 'paired' requires treatment_group and control_group to have the same length. ",
        "Use contrast_mode = 'group_mean' for unpaired or unbalanced groups.",
        call. = FALSE
      )
    }

    control <- c(sapply(control_group, function(x) E_centered[, x]))
    treatment <- c(sapply(treatment_group, function(x) E_centered[, x]))
    delta <- treatment - control
    N_augmented <- do.call(rbind, replicate(length(control_group), N, simplify = FALSE))
  } else if (identical(contrast_mode, "group_mean")) {
    treatment <- rowMeans(E_centered[, treatment_group, drop = FALSE])
    control <- rowMeans(E_centered[, control_group, drop = FALSE])
    delta <- treatment - control
    N_augmented <- N
  } else {
    delta_parts <- list()
    idx <- 1L
    for (treatment_sample in treatment_group) {
      for (control_sample in control_group) {
        delta_parts[[idx]] <- E_centered[, treatment_sample] - E_centered[, control_sample]
        idx <- idx + 1L
      }
    }
    delta <- unlist(delta_parts, use.names = FALSE)
    N_augmented <- do.call(rbind, replicate(length(delta_parts), N, simplify = FALSE))
  }

  if (length(delta) != nrow(N_augmented)) {
    stop(
      sprintf(
        "Internal contrast construction failed: delta has length %d but N_augmented has %d rows.",
        length(delta),
        nrow(N_augmented)
      ),
      call. = FALSE
    )
  }

  list(
    delta = delta,
    N_augmented = N_augmented,
    contrast_mode = contrast_mode
  )
}

getSignifActivities <- function(E_centered,
                                N,
                                treatment_group,
                                control_group,
                                exclude_genes = NULL,
                                contrast_mode = c("auto", "paired", "group_mean", "all_pairs")) {
  print("Starting new sample")
  contrast_mode <- match.arg(contrast_mode)

  stopifnot(dim(E_centered)[1] == dim(N)[1])

  stopifnot(all(sapply(treatment_group, function(x) x %in% colnames(E_centered))))
  stopifnot(all(sapply(control_group, function(x) x %in% colnames(E_centered))))
  stopifnot(all(sapply(control_group, function(x) !x %in% treatment_group)))

  genes_excluded <- c()
  if (!is.null(exclude_genes)) {
    genes_excluded <- exclude_genes[exclude_genes %in% rownames(E_centered)]
    if (length(genes_excluded) > 0) {
      to_delete <- which(rownames(E_centered) %in% genes_excluded)
      E_centered <- E_centered[-to_delete, , drop = FALSE]
      N <- N[-to_delete, , drop = FALSE]
    }
  }
  stopifnot(all(rownames(E_centered) == rownames(N)))

  contrast <- .get_signif_delta_design(
    E_centered = E_centered,
    N = N,
    treatment_group = treatment_group,
    control_group = control_group,
    contrast_mode = contrast_mode
  )

  design_matrix <- as.data.frame(cbind(delta = contrast$delta, contrast$N_augmented))
  fit <- lm(delta ~ ., data = design_matrix)
  coefs <- as.data.frame(summary(fit)$coefficients)

  p_adj <- p.adjust(coefs[-1, 4], method = "BH")
  p_adj <- append(p_adj, NA, 0)
  coefs$p_adj <- p_adj

  fit_summary <- summary(fit)

  res <- list()
  res$coefs <- coefs
  res$control_group <- control_group
  res$treatment_group <- treatment_group
  res$r.squared <- fit_summary$r.squared
  res$ajd.r.squared <- fit_summary$adj.r.squared
  res$adj.r.squared <- fit_summary$adj.r.squared
  res$aliased <- fit_summary$aliased
  res$fstatistic <- fit_summary$fstatistic
  res$rss <- deviance(fit)
  res$p <- ncol(N)
  res$N <- nrow(E_centered)
  res$genes_excluded <- genes_excluded
  res$contrast_mode <- contrast$contrast_mode
  res$n_treatment <- length(treatment_group)
  res$n_control <- length(control_group)
  res$n_delta_observations <- length(contrast$delta)
  res
}
