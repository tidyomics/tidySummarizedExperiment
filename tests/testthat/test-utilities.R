test_that("tidySummarizedExperiment:::get_count_datasets works", {
    # Use testthat edition 3
    local_edition(3)
    
    # Generate test SE
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    rownames(se) <- paste0("G", seq_len(3))
    colnames(se) <- paste0("S", seq_len(3))
    SummarizedExperiment::assays(se) <- lapply(
        SummarizedExperiment::assays(se),
        function(x) {
            rownames(x) <- rownames(se)
            colnames(x) <- colnames(se)
            x
        }
    )
    
    # Check that dimnames are as expected
    expect_equal(rownames(se), rownames(assay(se, "mat1", withDimnames = FALSE)))
    expect_equal(rownames(se), rownames(assay(se, "mat2", withDimnames = FALSE)))
    expect_equal(rownames(se), rownames(assay(se, "mat3", withDimnames = FALSE)))
    expect_equal(colnames(se), colnames(assay(se, "mat1", withDimnames = FALSE)))
    expect_equal(colnames(se), colnames(assay(se, "mat2", withDimnames = FALSE)))
    expect_equal(colnames(se), colnames(assay(se, "mat3", withDimnames = FALSE)))
    # Check that matrix values are as expected
    expect_equal(assay(se, "mat1")[, 1], c(1, 2, 3), ignore_attr = TRUE)
    expect_equal(assay(se, "mat2")[, 2], c(13, 14, 15), ignore_attr = TRUE)
    expect_equal(assay(se, "mat3")[, 3], c(25, 26, 27), ignore_attr = TRUE)
    
    # SE and all assays have the same dimnames
    cds <- tidySummarizedExperiment:::get_count_datasets(se)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, all assays have the same
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    expect_equal(colnames(assay(se1, "mat1", withDimnames = FALSE)), paste0("S", seq_len(3)))
    expect_equal(rownames(assay(se1, "mat1", withDimnames = FALSE)), paste0("G", seq_len(3)))
    expect_null(colnames(se1))
    expect_null(rownames(se1))
    expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "have column names, but the SummarizedExperiment does not"), "have row names, but the SummarizedExperiment does not")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, assays 1 and 3 have the same, assay 2 does not have
    # SE dimnames will be set to those of assay 1, then assay 2 dimnames to those of the SE
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- 
        colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    expect_warning(expect_warning(expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names"))))
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, assays 2 and 3 have the same, assay 1 does not have
    # SE dimnames will be set to those of assay 2, then assay 1 dimnames to those of the SE
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- 
        colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    expect_warning(expect_warning(expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names"))))
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, assays have the same but in different order
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- c("S2", "S3", "S1")
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- c("G2", "G3", "G1")
    expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1)))
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, c(16, 17, 18, 10, 11, 12, 13, 14, 15))
    expect_equal(cds$mat3, c(21, 19, 20, 24, 22, 23, 27, 25, 26))
    
    # SE does not have dimnames, assays have nonoverlapping dimnames
    ## ...rownames
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- paste0("A", seq_len(3))
    expect_warning(expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have row names")))
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 18L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, c(rep(paste0("G", seq_len(3)), 3), 
                                 rep(paste0("A", seq_len(3)), 3)))
    expect_equal(cds$.sample, rep(rep(paste0("S", seq_len(3)), each = 3), 2))
    expect_equal(cds$mat1, c(seq(1, 9), rep(NA, 9)))
    expect_equal(cds$mat2, c(rep(NA, 9), seq(10, 18)))
    expect_equal(cds$mat3, c(seq(19, 27), rep(NA, 9)))
    
    ## ...colnames
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- paste0("A", seq_len(3))
    expect_warning(expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names")))
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 18L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(rep(paste0("G", seq_len(3)), 3), 2))
    expect_equal(cds$.sample, c(rep(paste0("S", seq_len(3)), each = 3),
                                rep(paste0("A", seq_len(3)), each = 3)))
    expect_equal(cds$mat1, c(seq(1, 9), rep(NA, 9)))
    expect_equal(cds$mat2, c(rep(NA, 9), seq(10, 18)))
    expect_equal(cds$mat3, c(seq(19, 27), rep(NA, 9)))
    
    # Neither SE nor assays have column names
    se1 <- se
    colnames(se1) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- tidySummarizedExperiment:::get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(as.character(seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # Neither SE nor assays have row names
    se1 <- se
    rownames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- tidySummarizedExperiment:::get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(as.character(seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # Neither SE nor assays have any dimnames
    se1 <- se
    rownames(se1) <- NULL
    colnames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- tidySummarizedExperiment:::get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(as.character(seq_len(3)), 3))
    expect_equal(cds$.sample, rep(as.character(seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE has dimnames, assays have the same dimnames but not overlapping with those of the SE
    se1 <- se
    rownames(se1) <- colnames(se1) <- paste0("A", seq_len(3))
    expect_warning(expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "don't agree with the column names of the SummarizedExperiment"))
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    se1 <- se
    rownames(se1) <- paste0("A", seq_len(3))
    expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "don't agree with the row names of the SummarizedExperiment")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE has dimnames, none of the assays have
    se1 <- se
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- tidySummarizedExperiment:::get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, one assay has duplicated colnames, one has no colnames
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE))[2] <- 
        colnames(assay(se1, "mat1", withDimnames = FALSE))[1]
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    expect_equal(colnames(assay(se1, "mat1", withDimnames = FALSE)), paste0("S", c(1, 1, 3)))
    expect_equal(rownames(assay(se1, "mat1", withDimnames = FALSE)), paste0("G", seq_len(3)))
    expect_null(colnames(assay(se1, "mat2", withDimnames = FALSE)))
    expect_equal(colnames(assay(se1, "mat3", withDimnames = FALSE)), paste0("S", seq_len(3)))
    expect_null(colnames(se1))
    expect_null(rownames(se1))
    expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "some column names are duplicated")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    
    # SE does not have dimnames, one assay has duplicated rownames, one has no rownames
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE))[2:3] <- 
        rownames(assay(se1, "mat1", withDimnames = FALSE))[1]
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    expect_equal(rownames(assay(se1, "mat1", withDimnames = FALSE)), paste0("G", c(1, 1, 1)))
    expect_equal(colnames(assay(se1, "mat1", withDimnames = FALSE)), paste0("S", seq_len(3)))
    expect_null(rownames(assay(se1, "mat2", withDimnames = FALSE)))
    expect_equal(rownames(assay(se1, "mat3", withDimnames = FALSE)), paste0("G", seq_len(3)))
    expect_null(colnames(se1))
    expect_null(rownames(se1))
    expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "some row names are duplicated")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    
    # SE has duplicated colnames
    se1 <- se
    colnames(se1) <- paste0("S", c(1, 1, 1))
    expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "some column names are duplicated")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_true(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    
    # SE has duplicated rownames
    se1 <- se
    rownames(se1) <- paste0("G", c(1, 2, 1))
    expect_warning(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "some row names are duplicated")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_true(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    
    # All assays + SE have duplicated colnames
    se1 <- se
    colnames(se1)[2] <- 
        colnames(assay(se1, "mat1", withDimnames = FALSE))[2] <- 
        colnames(assay(se1, "mat2", withDimnames = FALSE))[2] <- 
        colnames(assay(se1, "mat3", withDimnames = FALSE))[2] <- "S1"
    expect_true(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    expect_false(check_if_assays_are_NOT_overlapped(se1, dim = "cols"))
    expect_false(check_if_assays_are_NOT_overlapped(se1, dim = "rows"))
    
    # Two assays + SE have duplicated colnames
    se1 <- se
    colnames(se1)[2] <- 
        colnames(assay(se1, "mat1", withDimnames = FALSE))[2] <- 
        colnames(assay(se1, "mat3", withDimnames = FALSE))[2] <- "S1"
    expect_true(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    expect_true(check_if_assays_are_NOT_overlapped(se1, dim = "cols"))
    expect_false(check_if_assays_are_NOT_overlapped(se1, dim = "rows"))
    
    # Assays have duplicated colnames in different ways
    se1 <- se
    assay(se1, "mat2") <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- c("S1", "S1", "S2")
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- c("S1", "S2", "S2")
    expect_true(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    expect_true(check_if_assays_are_NOT_overlapped(se1, dim = "cols"))
    expect_false(check_if_assays_are_NOT_overlapped(se1, dim = "rows"))

    # All dimnames are NULL - not duplicated
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- 
        colnames(assay(se1, "mat1", withDimnames = FALSE)) <- 
        rownames(assay(se1, "mat2", withDimnames = FALSE)) <- 
        colnames(assay(se1, "mat2", withDimnames = FALSE)) <- 
        rownames(assay(se1, "mat3", withDimnames = FALSE)) <- 
        colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    
    # Unnamed assay(s)
    # se1 <- SummarizedExperiment::SummarizedExperiment(
    #     assays = list(
    #         matrix(seq_len(9), nrow = 3),
    #         mat2 = matrix(seq(10, 18), nrow = 3),
    #         matrix(seq(19, 27), nrow = 3)
    #     )
    # )
    # expect_error(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "assays must be named")
    # 
    # se1 <- SummarizedExperiment::SummarizedExperiment(
    #     assays = list(
    #         matrix(seq_len(9), nrow = 3),
    #         matrix(seq(10, 18), nrow = 3),
    #         matrix(seq(19, 27), nrow = 3)
    #     )
    # )
    # expect_error(cds <- tidySummarizedExperiment:::get_count_datasets(se1), "assays must be named")
})
