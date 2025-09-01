context("filter function tests with scope and decomposition")

library(tidySummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(airway)

# Helper: airway-based test SE (subset for speed)
create_airway_test_se <- function() {
    data(airway)
    se <- airway
    se <- se[1:200, ]
    se
}

test_that("filter handles colData-only predicates (samples)", {
    se <- create_airway_test_se()

    # Reference using colData
    keep_samples <- which(colData(se)$dex == "trt")

    result <- se %>% filter(dex == "trt")

    expect_equal(nrow(result), nrow(se))
    expect_equal(ncol(result), length(keep_samples))
    expect_equal(colnames(result), colnames(se)[keep_samples])
    expect_equal(S4Vectors::metadata(result)$latest_filter_scope_report$scope, "coldata_only")
})

test_that("filter handles rowData-only predicates (features)", {
    se <- create_airway_test_se()

    # Reference using rowData
    keep_rows <- which(rowData(se)$gene_biotype == "protein_coding")

    result <- se %>% filter(gene_biotype == "protein_coding")

    expect_equal(ncol(result), ncol(se))
    expect_equal(nrow(result), length(keep_rows))
    expect_equal(rownames(result), rownames(se)[keep_rows])
    expect_equal(S4Vectors::metadata(result)$latest_filter_scope_report$scope, "rowdata_only")
})

## Mixed predicates across slots are handled by plyxp engine, which may enforce
## explicit rows()/cols() wrappers; covered elsewhere. We skip here to avoid flakiness.

test_that("is_composed detects multiple filter predicates", {
    se <- create_airway_test_se()
    expect_true(tidySummarizedExperiment:::is_composed("filter", dex == "trt", cell == cell[1]))
    expect_false(tidySummarizedExperiment:::is_composed("filter", dex == "trt"))
})

test_that("filter metadata stores latest scope report", {
    se <- create_airway_test_se()

    res1 <- se %>% filter(dex == "trt")
    expect_true(!is.null(S4Vectors::metadata(res1)$latest_filter_scope_report))
    expect_equal(S4Vectors::metadata(res1)$latest_filter_scope_report$scope, "coldata_only")

    res2 <- se %>% filter(gene_biotype == "protein_coding")
    expect_true(!is.null(S4Vectors::metadata(res2)$latest_filter_scope_report))
    expect_equal(S4Vectors::metadata(res2)$latest_filter_scope_report$scope, "rowdata_only")
})

test_that("internal analyze_query_scope_filter detects scopes", {
    se <- create_airway_test_se()

    scope <- tidySummarizedExperiment:::analyze_query_scope_filter(se, dex == "trt")
    expect_equal(scope$scope, "coldata_only")

    scope <- tidySummarizedExperiment:::analyze_query_scope_filter(se, gene_biotype == "protein_coding")
    expect_equal(scope$scope, "rowdata_only")

    scope <- tidySummarizedExperiment:::analyze_query_scope_filter(se, gene_biotype == "protein_coding", dex == "trt")
    expect_true(scope$scope %in% c("mixed", "unknown"))
})


test_that("assay-only filter yields non-rectangular tibble fallback", {
    se <- create_airway_test_se()

    # Reference tibble with non-rectangular selection using assay values
    res_tbl <- se %>%
        as_tibble(skip_GRanges = TRUE) %>%
        dplyr::filter(counts > 0)

    # Apply through SE interface; dependency analysis marks assay_only and we fallback
    res <- se %>% filter(counts > 0)

    expect_true(tibble::is_tibble(res))
    expect_identical(colnames(res_tbl), colnames(res))
    expect_identical(nrow(res_tbl), nrow(res))
})


