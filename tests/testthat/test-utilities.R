test_that("get_count_datasets works", {
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
    cds <- get_count_datasets(se)
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
    expect_warning(expect_warning(cds <- get_count_datasets(se1), "have column names, but the SummarizedExperiment does not"), "have row names, but the SummarizedExperiment does not")
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
    expect_warning(expect_warning(expect_warning(expect_warning(cds <- get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names"))))
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
    expect_warning(expect_warning(expect_warning(expect_warning(cds <- get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names"))))
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
    expect_warning(expect_warning(cds <- get_count_datasets(se1)))
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
    expect_warning(expect_warning(expect_warning(cds <- get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have row names")))
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
    expect_warning(expect_warning(expect_warning(cds <- get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names")))
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
    cds <- get_count_datasets(se1)
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
    cds <- get_count_datasets(se1)
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
    cds <- get_count_datasets(se1)
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
    expect_warning(expect_warning(cds <- get_count_datasets(se1), "don't agree with the column names of the SummarizedExperiment"))
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
    expect_warning(cds <- get_count_datasets(se1), "don't agree with the row names of the SummarizedExperiment")
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
    cds <- get_count_datasets(se1)
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
    expect_error(cds <- get_count_datasets(se1), "some column names are duplicated")
    
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
    expect_error(cds <- get_count_datasets(se1), "some row names are duplicated")
    
    # SE has duplicated colnames
    se1 <- se
    colnames(se1) <- paste0("S", c(1, 1, 1))
    expect_error(cds <- get_count_datasets(se1), "some column names are duplicated")
    expect_true(check_if_any_dimnames_duplicated(se1, dim = "cols"))
    expect_false(check_if_any_dimnames_duplicated(se1, dim = "rows"))
    
    # SE has duplicated rownames
    se1 <- se
    rownames(se1) <- paste0("G", c(1, 2, 1))
    expect_error(cds <- get_count_datasets(se1), "some row names are duplicated")
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
    # expect_error(cds <- get_count_datasets(se1), "assays must be named")
    # 
    # se1 <- SummarizedExperiment::SummarizedExperiment(
    #     assays = list(
    #         matrix(seq_len(9), nrow = 3),
    #         matrix(seq(10, 18), nrow = 3),
    #         matrix(seq(19, 27), nrow = 3)
    #     )
    # )
    # expect_error(cds <- get_count_datasets(se1), "assays must be named")
})

test_that("check_if_assays_are_NOT_consistently_ordered works correctly", {
    # Use testthat edition 3
    local_edition(3)
    
    # Test 1: SE with no assays - should return FALSE
    se_empty <- SummarizedExperiment::SummarizedExperiment(
        assays = list()
    )
    expect_false(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_empty))
    
    # Test 2: SE with assays but no colnames - should return FALSE
    se_no_colnames <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3)
        )
    )
    expect_false(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_no_colnames))
    
    # Test 3: SE with assays, some with colnames, some without - should return FALSE
    se_mixed <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3)
        )
    )
    colnames(assay(se_mixed, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_false(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_mixed))
    
    # Test 4: SE with assays, all with same colnames in same order - should return FALSE
    se_consistent <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3)
        )
    )
    colnames(assay(se_consistent, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_consistent, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_false(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_consistent))
    
    # Test 5: SE with assays, all with same colnames but different order - should return TRUE
    se_inconsistent <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3)
        )
    )
    colnames(assay(se_inconsistent, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_inconsistent, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_inconsistent))
    
    # Test 6: SE with assays, all with same colnames but different order (reverse) - should return TRUE
    se_reverse <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3)
        )
    )
    colnames(assay(se_reverse, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_reverse, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 2, 1))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_reverse))
    
    # Test 7: SE with assays, all with same colnames but different order (middle swap) - should return TRUE
    se_middle_swap <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3)
        )
    )
    colnames(assay(se_middle_swap, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_middle_swap, "mat2", withDimnames = FALSE)) <- paste0("S", c(1, 3, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_middle_swap))
    
    # Test 8: SE with 3 assays, first two consistent, third different - should return TRUE
    se_three_assays <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_three_assays, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_three_assays, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_three_assays, "mat3", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_three_assays))
    
    # Test 9: SE with 3 assays, all consistent - should return FALSE
    se_three_consistent <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_three_consistent, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_three_consistent, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_three_consistent, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_false(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_three_consistent))
    
    # Test 10: SE with assays, all with same colnames but different order (complex permutation) - should return TRUE
    se_complex <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(12), nrow = 3),
            mat2 = matrix(seq(13, 24), nrow = 3)
        )
    )
    colnames(assay(se_complex, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(4))
    colnames(assay(se_complex, "mat2", withDimnames = FALSE)) <- paste0("S", c(4, 2, 1, 3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_complex))
    
    # Test 11: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_all_different, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different, "mat3", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different))
    
    # Test 12: SE with assays, all with same colnames but different order (pairwise different) - should return TRUE
    se_pairwise <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_pairwise, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_pairwise, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_pairwise, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_pairwise))
    
    # Test 13: SE with assays, all with same colnames but different order (first and last same) - should return TRUE
    se_first_last <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_first_last, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_first_last, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_first_last, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_first_last))
    
    # Test 14: SE with assays, all with same colnames but different order (middle different) - should return TRUE
    se_middle_different <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_middle_different, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_middle_different, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_middle_different, "mat3", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_middle_different))
    
    # Test 15: SE with assays, all with same colnames but different order (first different) - should return TRUE
    se_first_different <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_first_different, "mat1", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_first_different, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_first_different, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_first_different))
    
    # Test 16: SE with assays, all with same colnames but different order (last different) - should return TRUE
    se_last_different <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_last_different, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_last_different, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_last_different, "mat3", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_last_different))
    
    # Test 17: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_2 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_2, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_all_different_2, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different_2, "mat3", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_2))
    
    # Test 18: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_3 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_3, "mat1", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different_3, "mat2", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_3, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_3))
    
    # Test 19: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_4 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_4, "mat1", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_4, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_all_different_4, "mat3", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_4))
    
    # Test 20: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_5 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_5, "mat1", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_5, "mat2", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different_5, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_5))
    
    # Test 21: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_6 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_6, "mat1", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different_6, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_all_different_6, "mat3", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_6))
    
    # Test 22: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_7 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_7, "mat1", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_all_different_7, "mat2", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_7, "mat3", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_7))
    
    # Test 23: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_8 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_8, "mat1", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different_8, "mat2", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_8, "mat3", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_8))
    
    # Test 24: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_9 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_9, "mat1", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_9, "mat2", withDimnames = FALSE)) <- paste0("S", seq_len(3))
    colnames(assay(se_all_different_9, "mat3", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_9))
    
    # Test 25: SE with assays, all with same colnames but different order (all different) - should return TRUE
    se_all_different_10 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    colnames(assay(se_all_different_10, "mat1", withDimnames = FALSE)) <- paste0("S", c(3, 1, 2))
    colnames(assay(se_all_different_10, "mat2", withDimnames = FALSE)) <- paste0("S", c(2, 3, 1))
    colnames(assay(se_all_different_10, "mat3", withDimnames = FALSE)) <- paste0("S", c(1, 2, 3))
    expect_true(tidySummarizedExperiment:::check_if_assays_are_NOT_consistently_ordered(se_all_different_10))
})
