context("dplyr test")

library(tidySummarizedExperiment)


test_that("append_samples", {
    pasilla_bind <- append_samples(pasilla, pasilla)

    # Check that the combined object has the expected number of samples
    expect_equal(ncol(pasilla_bind), ncol(pasilla) * 2)
    
    # Check that it's still a SummarizedExperiment
    expect_true(inherits(pasilla_bind, "SummarizedExperiment"))
})

test_that("distinct", {
    pasilla %>%
        distinct(condition) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("filter", {
    pasilla %>%
        filter(condition == "untreated") %>%
        nrow() %>%
        expect_equal(14599)
})

test_that("group_by", {
    pasilla %>%
        group_by(condition) %>%
        ncol() %>%
        expect_equal(5)
})

test_that("summarise", {
    pasilla %>%
        summarise(mean(counts)) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("mutate", {
    pasilla %>%
        mutate(condition = 1) %>%
        distinct(condition) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("rename", {
    pasilla %>%
        rename(groups = condition, type_2 = type) %>%
        select(groups) %>%
        ncol() %>%
        expect_equal(1)
})



test_that("slice", {
    pasilla %>%
        slice(1) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("select", {
    pasilla %>%
        select(-condition) %>%
        class() %>%
        as.character() %>%
        expect_equal("SummarizedExperiment")

    pasilla %>%
        select(condition) %>%
        class() %>%
        as.character() %>%
        .[1] %>%
        expect_equal("tbl_df")
})

test_that("sample_n", {
    pasilla %>%
        sample_n(50) %>%
        nrow() %>%
        expect_equal(50)
})

test_that("sample_frac", {
    pasilla %>%
        sample_frac(0.1) %>%
        nrow() %>%
        expect_equal(10219)
})

test_that("count", {
    pasilla %>%
        count(condition) %>%
        nrow() %>%
        expect_equal(2)
})


test_that("mutate counts", {
  
  se = tidySummarizedExperiment::pasilla |> mutate(counts_2 = counts) 

  se |> 
    pull(counts) |> 
    expect_equal(
      se |> pull(counts_2)
    )
  
  se = tidySummarizedExperiment::pasilla 
  assays(se, withDimnames = FALSE)$counts_2 = assays(se)$counts[,7:1]
  
  se |> 
    pull(counts) |> 
    expect_equal(
      se |> pull(counts_2)
    )
  
  se |> 
  tidySummarizedExperiment:::check_if_assays_are_NOT_overlapped(dim = "cols") |> 
    expect_equal(FALSE)
  
  se[,1] |> 
    tidySummarizedExperiment:::check_if_assays_are_NOT_overlapped(dim = "cols") |> 
    expect_equal(TRUE)
  
  })


test_that("group_split splits character columns", {
  data(pasilla)
  pasilla |> 
    group_split(condition) |> 
    length() |> 
    expect_equal(2)
})

test_that("group_split splits logical comparisons", {
  data(pasilla)
  pasilla |> 
    group_split(counts > 0) |> 
    length() |> 
    expect_equal(2)
})

test_that("group_split splits with mutliple arguments", {
  data(pasilla)
  pasilla |> 
    group_split(condition, counts > 0) |> 
    length() |> 
    expect_equal(4)
})

test_that("mutate features", {
  pasilla %>%
    mutate_features(new = 1:nrow(pasilla)) %>%
    rowData() %>%
    as_tibble() %>%
    pull(new) %>%
    expect_equal(1:nrow(pasilla))
})


test_that("mutate samples", {
  pasilla %>%
    mutate_samples(new = 1:ncol(pasilla)) %>%
    colData() %>%
    as_tibble() %>%
    pull(new) %>%
    expect_equal(1:ncol(pasilla))
})
