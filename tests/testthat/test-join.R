context("join tests")

library(tidySummarizedExperiment)


test_that("left_join", {
    res_left <- pasilla %>%
        left_join(pasilla %>%
                      distinct(condition) %>%
                      mutate(new_column = 1:2), by = "condition")

    expect_equal(
        res_left %>%
            colData() %>%
            ncol(),
        pasilla %>%
            colData() %>%
            ncol() %>%
            sum(1)
    )

    # Metadata: scope should be coldata_only when joining by a colData key (condition)
    meta <- S4Vectors::metadata(res_left)
    expect_true(!is.null(meta$latest_join_scope_report))
    expect_equal(meta$latest_join_scope_report$scope, "coldata_only")
})

test_that("left_join 0 samples", {
 
    res_zero <- pasilla[0,] %>%
      left_join(pasilla %>%
                  distinct(condition) %>%
                  mutate(new_column = 1), by = "condition")

    # still SE; metadata should be set
    expect_s4_class(res_zero, "SummarizedExperiment")
    meta0 <- S4Vectors::metadata(res_zero)
    expect_true(!is.null(meta0$latest_join_scope_report))
    expect_equal(meta0$latest_join_scope_report$scope, "coldata_only")

    res_zero %>%
      as_tibble() |> 
      pull(new_column) %>%
      unique() |> 
      expect_equal(1)
  
})

test_that("inner_join", {
    res_inner <- pasilla %>% inner_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1), by = "condition")
    res_inner %>% ncol() %>% expect_equal(4)

    meta <- S4Vectors::metadata(res_inner)
    expect_true(!is.null(meta$latest_join_scope_report))
    expect_equal(meta$latest_join_scope_report$scope, "coldata_only")
})

test_that("right_join", {
    res_right <- pasilla %>% right_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1), by = "condition")
    res_right %>% ncol() %>% expect_equal(4)

    meta <- S4Vectors::metadata(res_right)
    expect_true(!is.null(meta$latest_join_scope_report))
    expect_equal(meta$latest_join_scope_report$scope, "coldata_only")
})

test_that("full_join", {
    res_full <- pasilla %>%
        full_join(tibble::tibble(condition = "A",     other = 1:4), by = "condition")
    res_full %>% nrow() %>% expect_equal(102197)

    if (methods::is(res_full, "SummarizedExperiment")) {
      meta <- S4Vectors::metadata(res_full)
      expect_true(!is.null(meta$latest_join_scope_report))
      expect_equal(meta$latest_join_scope_report$scope, "coldata_only")
    }
})
