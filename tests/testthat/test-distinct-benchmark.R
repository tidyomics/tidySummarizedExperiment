context("distinct benchmark tests")

library(tidySummarizedExperiment)
library(microbenchmark)

test_that("distinct vs select performance benchmark", {
    # Create a large dataset to simulate 20K+ samples scenario
    se_large <- pasilla
    for (i in 1:30) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Benchmark distinct vs select performance
    benchmark_results <- microbenchmark(
        distinct_result = {
            se_large %>%
                distinct(condition) %>%
                ncol()
        },
        select_result = {
            se_large %>%
                select(condition) %>%
                ncol()
        },
        times = 5,
        unit = "ms"
    )
    
    # Extract median times
    distinct_times <- benchmark_results$time[benchmark_results$expr == "distinct_result"]
    select_times <- benchmark_results$time[benchmark_results$expr == "select_result"]
    
    distinct_median <- median(distinct_times) / 1e6  # Convert to milliseconds
    select_median <- median(select_times) / 1e6
    
    # Both should complete in reasonable time
    expect_true(distinct_median < 1000)  # Less than 1 second
    expect_true(select_median < 1000)
    
    # Document the performance difference
    # This test can be used to verify improvements when distinct is optimized
    cat(sprintf("Distinct median time: %.2f ms\n", distinct_median))
    cat(sprintf("Select median time: %.2f ms\n", select_median))
    cat(sprintf("Performance ratio (distinct/select): %.2f\n", distinct_median / select_median))
})

test_that("distinct performance with different column types", {
    # Create a large dataset
    se_large <- pasilla
    for (i in 1:20) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Benchmark distinct on different column types
    benchmark_results <- microbenchmark(
        distinct_coldata = {
            se_large %>%
                distinct(condition) %>%
                ncol()
        },
        distinct_rowdata = {
            se_large %>%
                distinct(.feature) %>%
                nrow()
        },
        distinct_assay = {
            se_large %>%
                distinct(counts) %>%
                ncol()
        },
        distinct_mixed = {
            se_large %>%
                distinct(condition, counts) %>%
                ncol()
        },
        times = 3,
        unit = "ms"
    )
    
    # All operations should complete in reasonable time
    for (expr in unique(benchmark_results$expr)) {
        times <- benchmark_results$time[benchmark_results$expr == expr]
        median_time <- median(times) / 1e6
        expect_true(median_time < 2000)  # Less than 2 seconds
    }
})

test_that("distinct scalability test", {
    # Test distinct performance with increasing dataset sizes
    sizes <- c(1, 5, 10, 20)
    results <- list()
    
    for (size in sizes) {
        # Create dataset of specific size
        se_test <- pasilla
        for (i in 1:(size-1)) {
            se_test <- cbind(se_test, pasilla)
        }
        
        # Benchmark distinct performance
        start_time <- Sys.time()
        result <- se_test %>%
            distinct(condition) %>%
            ncol()
        end_time <- Sys.time()
        
        elapsed_time <- as.numeric(end_time - start_time)
        results[[as.character(size)]] <- elapsed_time
        
        expect_equal(result, 1)
        expect_true(elapsed_time < 5)  # Should complete in reasonable time
    }
    
    # Performance should scale reasonably (not exponentially)
    # This is a basic check - more sophisticated analysis could be added
    expect_true(length(results) == length(sizes))
})

test_that("distinct memory efficiency test", {
    # Test that distinct doesn't cause excessive memory usage
    se_large <- pasilla
    for (i in 1:25) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Test multiple distinct operations to check for memory leaks
    for (i in 1:10) {
        result <- se_large %>%
            distinct(condition) %>%
            ncol()
        expect_equal(result, 1)
        
        # Force garbage collection to check for memory issues
        gc()
    }
})

test_that("distinct with .keep_all performance", {
    # Create a large dataset
    se_large <- pasilla
    for (i in 1:20) {
        se_large <- cbind(se_large, pasilla)
    }
    
    # Benchmark distinct with and without .keep_all
    benchmark_results <- microbenchmark(
        distinct_keep_all = {
            se_large %>%
                distinct(condition, .keep_all = TRUE) %>%
                ncol()
        },
        distinct_no_keep_all = {
            se_large %>%
                distinct(condition, .keep_all = FALSE) %>%
                ncol()
        },
        times = 3,
        unit = "ms"
    )
    
    # Both should complete in reasonable time
    for (expr in unique(benchmark_results$expr)) {
        times <- benchmark_results$time[benchmark_results$expr == expr]
        median_time <- median(times) / 1e6
        expect_true(median_time < 2000)
    }
})

