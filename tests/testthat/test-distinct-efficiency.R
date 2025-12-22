context("distinct efficiency tests")

library(tidySummarizedExperiment)

test_that("distinct efficiency with large datasets", {
    # Create a large dataset by duplicating the pasilla dataset
    # This simulates the 20K+ samples scenario mentioned in the GitHub issue
    se_large <- pasilla
    
    # Duplicate the dataset multiple times to create many samples
    for (i in 1:20) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Test distinct on colData columns (should be efficient)
    start_time <- Sys.time()
    result_coldata <- se_large %>%
        distinct(condition) %>%
        ncol()
    end_time <- Sys.time()
    
    expect_equal(result_coldata, 1)
    expect_true(as.numeric(end_time - start_time) < 2)  # Should complete in reasonable time
    
    # Test distinct on multiple colData columns
    start_time <- Sys.time()
    result_multiple <- se_large %>%
        distinct(condition, type) %>%
        ncol()
    end_time <- Sys.time()
    
    expect_equal(result_multiple, 2)
    expect_true(as.numeric(end_time - start_time) < 2)
    
})

test_that("distinct performance comparison with select", {
    # Create a large dataset
    se_large <- pasilla
    for (i in 1:15) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Test distinct performance
    start_time_distinct <- Sys.time()
    result_distinct <- se_large %>%
        distinct(condition) %>%
        ncol()
    end_time_distinct <- Sys.time()
    time_distinct <- as.numeric(end_time_distinct - start_time_distinct)
    
    # Test select performance (should be faster due to optimization)
    start_time_select <- Sys.time()
    result_select <- se_large %>%
        select(condition) %>%
        ncol()
    end_time_select <- Sys.time()
    time_select <- as.numeric(end_time_select - start_time_select)
    
    expect_equal(result_distinct, 1)
    expect_equal(result_select, 1)
    
    # Both should complete in reasonable time
    expect_true(time_distinct < 2)
    expect_true(time_select < 1)
    
    # Note: distinct might be slower due to the current implementation
    # This test documents the current behavior and can be used to verify
    # improvements when distinct is optimized like select
})

test_that("distinct with mixed columns on large dataset", {
    # Create a large dataset
    se_large <- pasilla
    for (i in 1:10) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Test distinct with mixed columns (colData + assay)
    start_time <- Sys.time()
    result <- se_large %>%
        distinct(condition, counts) %>%
        ncol()
    end_time <- Sys.time()
    
    expect_equal(result, 2)
    expect_true(as.numeric(end_time - start_time) < 3)
})

test_that("distinct handles very large datasets gracefully", {
    # Create an even larger dataset
    se_very_large <- pasilla
    for (i in 1:50) {
        se_very_large <- cbind(se_very_large, pasilla)
    }
    
    # Test distinct on colData only (should be most efficient)
    start_time <- Sys.time()
    result <- se_very_large %>%
        distinct(condition) %>%
        ncol()
    end_time <- Sys.time()
    
    expect_equal(result, 1)
    # Should still complete in reasonable time even with very large dataset
    expect_true(as.numeric(end_time - start_time) < 5)
})

test_that("distinct memory usage with large datasets", {
    # Test that distinct doesn't cause memory issues
    se_large <- pasilla
    for (i in 1:25) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Test multiple distinct operations
    results <- list()
    for (i in 1:5) {
        result <- se_large %>%
            distinct(condition) %>%
            ncol()
        results[[i]] <- result
    }
    
    # All results should be the same
    expect_true(all(sapply(results, function(x) x == 1)))
})

test_that("distinct performance with grouped data on large dataset", {
    # Create a large dataset
    se_large <- pasilla
    for (i in 1:20) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Test distinct on grouped data
    start_time <- Sys.time()
    result <- se_large %>%
        group_by(condition) %>%
        distinct(type) %>%
        ncol()
    end_time <- Sys.time()
    
    expect_equal(result, 2)
    expect_true(as.numeric(end_time - start_time) < 3)
})

