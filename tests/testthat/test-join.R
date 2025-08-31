context("join tests")

library(tidySummarizedExperiment)


test_that("left_join", {
    expect_equal(
        pasilla %>%
            left_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2)) %>%
            colData() %>%
            ncol(),
        pasilla %>%
            colData() %>%
            ncol() %>%
            sum(1)
    )
})

test_that("left_join 0 samples", {
 
    pasilla[0,] %>%
      left_join(pasilla %>%
                  distinct(condition) %>%
                  mutate(new_column = 1)) |> 
    as_tibble() |> 
      pull(new_column) %>%
      unique() |> 
      expect_equal(1)
  
})

test_that("inner_join", {
    pasilla %>% inner_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1)) %>%
        ncol() %>%
        expect_equal(4)
})

test_that("right_join", {
    pasilla %>% right_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1)) %>%
        ncol() %>%
        expect_equal(4)
})

test_that("full_join", {
    pasilla %>%
        full_join(tibble::tibble(condition = "A",     other = 1:4)) %>% nrow() %>%
        expect_equal(102197)
})
