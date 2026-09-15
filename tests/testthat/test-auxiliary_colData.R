context("auxiliary_colData")

library(tidySummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)

create_small_se <- function() {
    se <- SummarizedExperiment(
        assays = list(counts = matrix(seq_len(12), nrow = 3, ncol = 4)),
        colData = DataFrame(condition = c("a", "a", "b", "b"))
    )
    rownames(se) <- paste0("G", seq_len(3))
    colnames(se) <- paste0("S", seq_len(4))
    se
}

test_that("auxiliary_colData default is empty and column-aligned", {
    se <- create_small_se()
    aux <- auxiliary_colData(se)
    expect_s4_class(aux, "DataFrame")
    expect_equal(nrow(aux), ncol(se))
    expect_equal(ncol(aux), 0L)
    expect_equal(rownames(aux), colnames(se))
})

test_that("get_colnames_col is unchanged when auxiliary_colData is empty", {
    se <- create_small_se()
    expect_equal(
        tidySummarizedExperiment:::get_colnames_col(se),
        c(colnames(colData(se)), ".sample")
    )
})

test_that("join_auxiliary_colData is a no-op when auxiliary columns are empty", {
    se <- create_small_se()
    sample_tbl <- tibble::tibble(
        .sample = colnames(se),
        condition = colData(se)$condition
    )
    expect_identical(
        tidySummarizedExperiment:::join_auxiliary_colData(sample_tbl, se),
        sample_tbl
    )
})

test_that("join_auxiliary_colData joins by sample id and skips existing columns", {
    se <- create_small_se()
    sample_tbl <- tibble::tibble(
        .sample = rev(colnames(se)),
        condition = rev(as.character(colData(se)$condition))
    )
    aux <- DataFrame(
        UMAP1 = seq_len(ncol(se)),
        condition = "ignored",
        row.names = colnames(se)
    )

    out <- tidySummarizedExperiment:::join_auxiliary_colData(sample_tbl, se, aux = aux)

    expect_equal(out$.sample, rev(colnames(se)))
    expect_equal(out$UMAP1, rev(seq_len(ncol(se))))
    expect_equal(out$condition, rev(as.character(colData(se)$condition)))
    expect_false(any(grepl("condition\\.", names(out))))
})

test_that("tbl_with_joined_auxiliary_columns reports no extra columns by default", {
    se <- create_small_se()
    out <- tidySummarizedExperiment:::tbl_with_joined_auxiliary_columns(se)
    expect_equal(out$extra, character(0))
    expect_equal(nrow(out$tbl), ncol(se))
})
