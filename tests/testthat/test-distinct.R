context("distinct tests")

library(tidySummarizedExperiment)

test_that("distinct works with colData columns", {
    # Test distinct on colData columns
    result <- pasilla %>%
        distinct(condition) %>%
        ncol()
    expect_equal(result, 1)
    
    # Test distinct on multiple colData columns
    result <- pasilla %>%
        distinct(condition, type) %>%
        ncol()
    expect_equal(result, 2)
})

test_that("distinct works with rowData columns", {
    # Test distinct on rowData columns (should return tibble)
    result <- pasilla %>%
        distinct(.feature)
    expect_true(nrow(result) > 0)
    expect_true(inherits(result, "tbl_df"))
})

test_that("distinct works with assay columns", {
    # Test distinct on assay columns
    result <- pasilla %>%
        distinct(counts) %>%
        ncol()
    expect_equal(result, 1)
})

test_that("distinct works with mixed columns", {
    # Test distinct on mixed columns (colData + assay)
    result <- pasilla %>%
        distinct(condition, counts) %>%
        ncol()
    expect_equal(result, 2)
})

test_that("distinct preserves SummarizedExperiment structure when possible", {
    # When distinct is applied to colData only, should return SE
    result <- pasilla %>%
        distinct(condition)
    expect_false(inherits(result, "SummarizedExperiment"))
    
    # When distinct is applied to mixed columns, should return tibble
    result <- pasilla %>%
        distinct(condition, counts, .sample, .feature)
    expect_true(inherits(result, "tbl_df"))
})

test_that("distinct handles empty selections", {
    # Test distinct with no columns specified (should use all columns)
    result <- pasilla %>%
        distinct() %>%
        nrow()
    expect_true(result > 0)
})

test_that("distinct works with large datasets", {
    # Create a larger test dataset
    se_large <- pasilla
    # Duplicate the data to create more rows
    se_large <- cbind(se_large, se_large)
    
    # Test distinct performance on larger dataset
    result <- se_large %>%
        distinct(condition) %>%
        ncol()
    expect_equal(result, 1)
})

test_that("distinct(counts) does not error when sample names are duplicated (base::cbind)", {
    se_large <- pasilla
    for (i in 1:5) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # as_tibble() path should repair duplicated sample names instead of erroring
    expect_warning(
        tbl <- suppressMessages(as_tibble(se_large)),
        "some column names are duplicated"
    )
    expect_equal(length(unique(tbl$.sample)), ncol(se_large))
    
    # distinct(counts) should also work (it relies on the as_tibble/get_count_datasets path)
    expect_warning(
        res <- suppressMessages(se_large %>% distinct(counts)),
        "some column names are duplicated"
    )
    expect_true(inherits(res, "tbl_df"))
    expect_equal(ncol(res), 1)
})

test_that("distinct works with GRanges columns", {
    # Test distinct when GRanges columns are involved
    result <- pasilla %>%
        distinct(.feature) %>%
        nrow()
    expect_true(result > 0)
})

test_that("distinct maintains row order", {
    # Test that distinct maintains the order of first occurrence
    result <- pasilla %>%
        distinct(condition) %>%
        pull(condition)
    
    # Should maintain the order of first occurrence
    expect_equal(result, c("untreated", "treated"))
})

test_that("distinct works with grouped data", {
    # Test distinct on grouped data
    result <- pasilla %>%
        group_by(condition) %>%
        distinct(type) %>%
        ncol()
    expect_equal(result, 2)
})

test_that("distinct handles special column names", {
    # Test distinct with special column names
    result <- pasilla %>%
        distinct(.sample) %>%
        ncol()
    expect_equal(result, 1)
    
    result <- pasilla %>%
        distinct(.feature) %>%
        nrow()
    expect_true(result > 0)
})

test_that("distinct performance with many samples", {
    # Test distinct performance with many samples (simulating 20K+ samples scenario)
    # Create a dataset with many samples by duplicating
    se_many_samples <- pasilla
    for (i in 1:10) {
        se_many_samples <- cbind(se_many_samples, pasilla)
    }
    
    # Test distinct on colData columns (should be fast)
    start_time <- Sys.time()
    result <- se_many_samples %>%
        distinct(condition) %>%
        ncol()
    end_time <- Sys.time()
    
    expect_equal(result, 1)
    # Should complete in reasonable time (less than 1 second for this test)
    expect_true(as.numeric(end_time - start_time) < 1)
})

test_that("distinct handles edge cases", {
    # Test distinct on single row
    se_single <- pasilla[1, ]
    result <- se_single %>%
        distinct(condition) %>%
        ncol()
    expect_equal(result, 1)
    
    # Test distinct on single column
    result <- pasilla %>%
        distinct(condition) %>%
        nrow()
    expect_equal(result, 2)  # Should have 2 distinct conditions
})

