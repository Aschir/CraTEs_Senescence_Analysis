testthat::local_edition(3)

project_dir <- Sys.getenv("CRATES_TEST_PROJECT_DIR")
if (!nzchar(project_dir)) {
  project_dir <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

old_wd <- setwd(project_dir)
on.exit(setwd(old_wd), add = TRUE)

source(file.path(project_dir, "scripts", "craTEs_own_function.R"))
source(file.path(project_dir, "scripts", "most_active_te_subfamilies_f.R"))

testthat::test_that("senescence condition groups are inferred from sample names", {
  groups <- infer_senescence_condition_groups(c("PRO1", "PRO2", "ESEN1", "LSEN1"))

  testthat::expect_identical(names(groups), c("PRO", "ESEN", "LSEN"))
  testthat::expect_identical(groups$PRO, c("PRO1", "PRO2"))
  testthat::expect_identical(groups$ESEN, "ESEN1")
  testthat::expect_identical(groups$LSEN, "LSEN1")
})

testthat::test_that("condition means are computed from named sample groups", {
  counts <- matrix(
    c(
      2, 4, 10, 14, 100,
      6, 8, 20, 24, 200
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("TE_A", "TE_B"), c("PRO1", "PRO2", "ESEN1", "ESEN2", "LSEN1"))
  )
  groups <- list(PRO = c("PRO1", "PRO2"), ESEN = c("ESEN1", "ESEN2"), LSEN = "LSEN1")

  means <- condition_mean_matrix(counts, groups)

  testthat::expect_equal(unname(means["TE_A", ]), c(3, 12, 100))
  testthat::expect_equal(unname(means["TE_B", ]), c(7, 22, 200))
})

testthat::test_that("repeat locus IDs are parsed to TE subfamilies", {
  parsed <- extract_te_subfamily_from_feature_id(c(
    "AluJr_7q11.23_852",
    "LTR5_Hs_19p13.3_4",
    "(A)n_10p12.1_51",
    "L1PA2"
  ))

  testthat::expect_identical(parsed, c("AluJr", "LTR5_Hs", "(A)n", "L1PA2"))
})

testthat::test_that("activity ranking defaults to sum abs mean activity", {
  activity_scores <- matrix(
    c(
      5, 5, 5, 5, 5, 5,
      1, 1, 4, 4, 1, 1,
      0, 0, 0, 0, 9, 9,
      2, 2, 2, 2, 2, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("TE_A", "TE_B", "TE_C", "TE_D"),
      c("PRO1", "PRO2", "ESEN1", "ESEN2", "LSEN1", "LSEN2")
    )
  )
  groups <- infer_senescence_condition_groups(colnames(activity_scores))

  ranked <- rank_activity_te_subfamilies(
    activity_scores = activity_scores,
    condition_groups = groups,
    k_values = c(3, 5)
  )

  testthat::expect_identical(ranked$te_subfamily[1:3], c("TE_A", "TE_C", "TE_B"))
  testthat::expect_equal(ranked$ranking_score[ranked$te_subfamily == "TE_A"], 15)
  testthat::expect_equal(ranked$ranking_score[ranked$te_subfamily == "TE_C"], 9)
})

testthat::test_that("activity ranking supports condition-specific scopes", {
  activity_scores <- matrix(
    c(
      1, 1, -8, -8, 0, 0,
      4, 4, 3, 3, 1, 1,
      2, 2, 1, 1, 9, 9
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("TE_A", "TE_B", "TE_C"),
      c("PRO1", "PRO2", "ESEN1", "ESEN2", "LSEN1", "LSEN2")
    )
  )
  groups <- infer_senescence_condition_groups(colnames(activity_scores))

  ranked <- rank_activity_te_subfamilies(
    activity_scores = activity_scores,
    condition_groups = groups,
    k_values = 3,
    scope = "ESEN"
  )

  testthat::expect_identical(ranked$te_subfamily[1:2], c("TE_A", "TE_B"))
  testthat::expect_equal(ranked$ranking_score[ranked$te_subfamily == "TE_A"], 8)
  testthat::expect_equal(ranked$signed_activity_score[ranked$te_subfamily == "TE_A"], -8)
})

testthat::test_that("activity ranking can still call the legacy max-sum helper explicitly", {
  activity_scores <- matrix(
    c(
      5, 5, 5, 5, 5, 5,
      1, 1, 4, 4, 1, 1,
      0, 0, 0, 0, 9, 9,
      2, 2, 2, 2, 2, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("TE_A", "TE_B", "TE_C", "TE_D"),
      c("PRO1", "PRO2", "ESEN1", "ESEN2", "LSEN1", "LSEN2")
    )
  )
  groups <- infer_senescence_condition_groups(colnames(activity_scores))

  ranked <- rank_activity_te_subfamilies(
    activity_scores = activity_scores,
    condition_groups = groups,
    k_values = c(3, 5),
    ranking_metric = "max_sum_abs_activity_score"
  )

  testthat::expect_identical(ranked$te_subfamily[1:3], c("TE_C", "TE_A", "TE_B"))
  testthat::expect_equal(ranked$ranking_score[ranked$te_subfamily == "TE_C"], 18)
  testthat::expect_equal(ranked$ranking_score[ranked$te_subfamily == "TE_A"], 10)
})

testthat::test_that("raw TE extraction keeps TE rows and resolves duplicates by max sum", {
  fixture <- data.frame(
    feature_id = c("GENE1", "TE_B_1p36.33_1", "TE_B_2q21.1_2", "TE_A_3p14.2_3", "TE_Z"),
    gene_or_te = c("gene", "repeat", "repeat", "repeat", "repeat"),
    GSM1 = c(100, 1, 7, 5, 4),
    GSM2 = c(100, 1, 3, 5, 5),
    check.names = FALSE
  )
  fixture_path <- tempfile(fileext = ".tsv")
  utils::write.table(
    fixture,
    file = fixture_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  raw_te <- load_raw_te_subfamily_counts(
    E_path = fixture_path,
    reference_te_subfamilies = c("TE_A", "TE_B"),
    duplicate_method = "Max Sum",
    verbose = FALSE
  )

  testthat::expect_identical(sort(rownames(raw_te$counts)), c("TE_A", "TE_B", "TE_Z"))
  testthat::expect_equal(unname(raw_te$counts["TE_B", "GSM1"]), 7)
  testthat::expect_equal(unname(raw_te$counts["TE_B", "GSM2"]), 3)
  testthat::expect_equal(raw_te$selection_summary$duplicate_groups_before, 1)
  testthat::expect_equal(raw_te$selection_summary$duplicate_groups_after, 0)
})

testthat::test_that("TE-Seq expression table averages conditions before duplicate handling", {
  fixture <- data.frame(
    gene_id = c("GENE1", "TE_A_1p36.33_1", "TE_A_2q21.1_2", "TE_B_3p14.2_3"),
    gene_or_te = c("gene", "repeat", "repeat", "repeat"),
    PRO1 = c(999, 30, 0, 1),
    ESEN1 = c(999, 0, 12, 3),
    ESEN2 = c(999, 0, 12, 5),
    ESEN3 = c(999, 0, 12, 7),
    LSEN1 = c(999, 0, 0, 9),
    check.names = FALSE
  )
  fixture_path <- tempfile(fileext = ".tsv")
  utils::write.table(
    fixture,
    file = fixture_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  te_expression <- load_te_expression_condition_means(
    E_path = fixture_path,
    duplicate_method = "Max Sum",
    condition_groups = list(PRO = "PRO1", ESEN = c("ESEN1", "ESEN2", "ESEN3"), LSEN = "LSEN1"),
    verbose = FALSE
  )

  testthat::expect_identical(sort(rownames(te_expression$condition_means)), c("TE_A", "TE_B"))
  testthat::expect_equal(unname(te_expression$condition_means["TE_A", ]), c(30, 0, 0))
  testthat::expect_equal(unname(te_expression$condition_means["TE_B", ]), c(1, 5, 9))
  testthat::expect_equal(te_expression$selection_summary$expression_rows_kept_as_te, 3)
  testthat::expect_equal(te_expression$selection_summary$duplicate_groups_before, 1)
  testthat::expect_equal(te_expression$selection_summary$duplicate_groups_after, 0)
})

testthat::test_that("TE expression condition means are ranked globally and by condition", {
  condition_means <- matrix(
    c(
      100, 0, 0,
      0, 50, 0,
      10, 10, 10
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("TE_A", "TE_B", "TE_C"), c("PRO", "ESEN", "LSEN"))
  )

  global_ranked <- rank_te_expression_condition_means(condition_means, k_values = 3, scope = "global")
  pro_ranked <- rank_te_expression_condition_means(condition_means, k_values = 3, scope = "PRO")

  testthat::expect_identical(global_ranked$te_subfamily, c("TE_A", "TE_B", "TE_C"))
  testthat::expect_equal(global_ranked$ranking_score, c(100, 50, 30))
  testthat::expect_identical(pro_ranked$te_subfamily, c("TE_A", "TE_C", "TE_B"))
  testthat::expect_equal(pro_ranked$ranking_score, c(100, 10, 0))
})

testthat::test_that("top-k intersections compare ranked activity and raw lists", {
  activity_ranked <- data.frame(
    rank = 1:4,
    te_subfamily = c("TE_A", "TE_B", "TE_C", "TE_D"),
    ranking_score = c(10, 8, 6, 1),
    stringsAsFactors = FALSE
  )
  raw_ranked <- data.frame(
    rank = 1:4,
    te_subfamily = c("TE_B", "TE_Z", "TE_A", "TE_D"),
    raw_expression_sum = c(50, 40, 30, 5),
    stringsAsFactors = FALSE
  )

  intersections <- build_top_k_intersections(
    activity_ranked = activity_ranked,
    raw_ranked = raw_ranked,
    k_values = c(2, 3)
  )

  testthat::expect_equal(intersections$summary$n_intersection, c(1, 2))
  testthat::expect_identical(intersections$summary$scope, c("global", "global"))
  testthat::expect_identical(
    intersections$members$te_subfamily[intersections$members$k == 3],
    c("TE_A", "TE_B")
  )
})

testthat::test_that("global intersection top n uses combined ranks", {
  activity_ranked <- data.frame(
    rank = 1:5,
    te_subfamily = c("TE_A", "TE_B", "TE_C", "TE_D", "TE_E"),
    ranking_score = c(10, 8, 6, 4, 2),
    stringsAsFactors = FALSE
  )
  raw_ranked <- data.frame(
    rank = 1:5,
    te_subfamily = c("TE_D", "TE_A", "TE_E", "TE_B", "TE_Z"),
    raw_expression_sum = c(50, 40, 30, 20, 10),
    stringsAsFactors = FALSE
  )

  global_top <- build_global_intersection_top_n(
    activity_ranked = activity_ranked,
    raw_ranked = raw_ranked,
    top_n = 3
  )

  testthat::expect_identical(global_top$te_subfamily, c("TE_A", "TE_D", "TE_B"))
  testthat::expect_identical(global_top$combined_rank_sum, c(3L, 5L, 6L))
})

testthat::test_that("top-k intersections can compare condition-specific expression ranks", {
  activity_ranked <- data.frame(
    scope = "PRO",
    rank = 1:3,
    te_subfamily = c("TE_A", "TE_B", "TE_C"),
    ranking_score = c(10, 8, 2),
    stringsAsFactors = FALSE
  )
  expression_ranked <- data.frame(
    scope = "PRO",
    rank = 1:3,
    te_subfamily = c("TE_B", "TE_Z", "TE_A"),
    ranking_score = c(50, 40, 30),
    stringsAsFactors = FALSE
  )

  intersections <- build_top_k_intersections(
    activity_ranked = activity_ranked,
    raw_ranked = expression_ranked,
    k_values = c(2, 3),
    raw_score_col = "ranking_score",
    scope = "PRO"
  )

  testthat::expect_equal(intersections$summary$n_intersection, c(1, 2))
  testthat::expect_true(all(intersections$summary$scope == "PRO"))
  testthat::expect_identical(
    intersections$members$te_expression_rank[intersections$members$k == 3],
    c(3L, 1L)
  )
})
