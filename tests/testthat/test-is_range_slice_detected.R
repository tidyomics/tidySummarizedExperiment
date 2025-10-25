test_that("is_range_slice_ungrouped_detected works correctly", {
  local_edition(3)
  
  # Load required libraries
  library(rlang)
  library(tidySummarizedExperiment)
  
  # Access the internal function
  is_range_slice_ungrouped_detected <- tidySummarizedExperiment:::is_range_slice_ungrouped_detected
  
  # Load test data
  data(pasilla)
  
  # Test with range slice (1:3) on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1:3), quo(NULL)))
  
  # Test with single number on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(1), quo(NULL)))
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(5), quo(NULL)))
  
  # Test with multiple single numbers on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(1, 2, 3), quo(NULL)))
  
  # Test with mixed arguments (range and single) on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1:3, 5), quo(NULL)))
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1, 2:4), quo(NULL)))
  
  # Test with character arguments on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list("a", "b"), quo(NULL)))
  
  # Test with mixed numeric and character on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(1, "a"), quo(NULL)))
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1:3, "a"), quo(NULL)))
  
  # Test with empty list on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(), quo(NULL)))
  
  # Test with NULL on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(quo(NULL)), quo(NULL)))
  
  # Test with logical on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(TRUE, FALSE), quo(NULL)))
  
  # Test with complex range on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(seq(1, 10, 2)), quo(NULL)))
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(c(1, 3, 5, 7)), quo(NULL)))
  
  # Test with negative ranges on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(-1:-3), quo(NULL)))
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1:-1), quo(NULL)))
  
  # Test with decimal ranges on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(seq(1, 3, 0.5)), quo(NULL)))
})

test_that("is_range_slice_ungrouped_detected handles edge cases", {
  local_edition(3)
  
  # Load required libraries
  library(rlang)
  library(tidySummarizedExperiment)
  
  # Access the internal function
  is_range_slice_ungrouped_detected <- tidySummarizedExperiment:::is_range_slice_ungrouped_detected
  
  # Load test data
  data(pasilla)
  
  # Test with very large ranges on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1:1000), quo(NULL)))
  
  # Test with single element vector on ungrouped data (should be FALSE)
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(c(1)), quo(NULL)))
  
  # Test with named arguments on ungrouped data
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(x = 1:3), quo(NULL)))
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(x = 1), quo(NULL)))
  
  # Test with nested lists on ungrouped data
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(list(1:3)), quo(NULL)))
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(list(1)), quo(NULL)))
})

test_that("is_range_slice_ungrouped_detected handles grouping correctly", {
  local_edition(3)
  
  # Load required libraries
  library(rlang)
  library(tidySummarizedExperiment)
  
  # Access the internal function
  is_range_slice_ungrouped_detected <- tidySummarizedExperiment:::is_range_slice_ungrouped_detected
  
  # Load test data
  data(pasilla)
  
  # Test ungrouped data with range slice
  expect_true(is_range_slice_ungrouped_detected(pasilla, list(1:3), quo(NULL)))
  
  # Test grouped data with range slice
  grouped_data <- pasilla |> group_by(.sample)
  expect_false(is_range_slice_ungrouped_detected(grouped_data, list(1:3), quo(NULL)))
  
  # Test with .by parameter
  expect_false(is_range_slice_ungrouped_detected(pasilla, list(1:3), quo(.sample)))
})

test_that("slice_optimised handles .by parameter correctly", {
  local_edition(3)
  
  # Access the internal function
  slice_optimised <- tidySummarizedExperiment:::slice_optimised
  
  # Load test data
  library(tidySummarizedExperiment)
  data(pasilla)
  
  # Test that .by parameter is passed through correctly
  # This is more of an integration test to ensure the parameter is handled
  expect_error(
    slice_optimised(pasilla, 1:3, .by = NULL, .preserve = FALSE),
    "slice using a range doesn't work on ungrouped data"
  )
  
  # Test with single slice (should work)
  result <- slice_optimised(pasilla, 1, .by = NULL, .preserve = FALSE)
  expect_s4_class(result, "SummarizedExperiment")
})

test_that("slice method works correctly with range detection", {
  local_edition(3)
  
  # Load test data
  library(tidySummarizedExperiment)
  data(pasilla)
  
  # Test range slice on ungrouped data - should throw error
  expect_error(
    pasilla |> dplyr::slice(1:3),
    "slice using a range doesn't work on ungrouped data"
  )
  
  # Test single slice on ungrouped data - should work
  result <- pasilla |> dplyr::slice(1)
  expect_s4_class(result, "SummarizedExperiment")
  
  # Test range slice on grouped data - should work
  result_grouped <- pasilla |> 
    group_by(.sample) |> 
    dplyr::slice(1:3)
  expect_s3_class(result_grouped, "grouped_df")
  
  # Test range slice with .by parameter - should work
  result_by <- pasilla |> 
    dplyr::slice(1:3, .by = .sample)
  expect_s4_class(result_by, "SummarizedExperiment")
})
