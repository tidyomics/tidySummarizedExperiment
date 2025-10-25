context("select function tests with scope analysis and fallback paths")

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

test_that("select coldata-only returns tibble expanded by rows and stores metadata", {
    se <- create_airway_test_se()

    res <- se %>% select(.sample, dex)

    expect_true(tibble::is_tibble(res))
    # Expect nrow expanded by number of features
    expect_equal(nrow(res), nrow(se) * ncol(se))
    expect_true(all(c(".sample", "dex") %in% colnames(res)))
})

test_that("select rowdata-only returns tibble expanded by columns and stores metadata", {
    se <- create_airway_test_se()

    res <- se %>% select(.feature, gene_biotype)

    expect_true(tibble::is_tibble(res))
    # Expect nrow expanded by number of samples
    expect_equal(nrow(res), ncol(se) * nrow(se))
    expect_true(all(c(".feature", "gene_biotype") %in% colnames(res)))
})

test_that("select assays-only without key columns falls back to tibble", {
    se <- create_airway_test_se()

    # Select only the assay 'counts' leaving key columns implicit
    res <- se %>% select(counts)

    expect_true(tibble::is_tibble(res))
    expect_true("counts" %in% colnames(res))
})

test_that("select mixed scope returns SE with selected row/col/assay pieces and metadata", {
    se <- create_airway_test_se()

    res <- se %>% select(.sample, .feature, counts)

    expect_s4_class(res, "SummarizedExperiment")
    expect_true("counts" %in% names(assays(res)))
    meta <- S4Vectors::metadata(res)
    expect_true(!is.null(meta$latest_select_scope_report))
    expect_true(meta$latest_select_scope_report$scope %in% c("mixed", "unknown"))
})

test_that("fallback when key columns missing converts to tibble and selects via dplyr", {
    se <- create_airway_test_se()

    # Force missing of key columns by selecting a non-key from colData only, but not .sample/.feature
    res <- se %>% select(dex)

    expect_true(tibble::is_tibble(res))
    expect_true("dex" %in% colnames(res))
})


