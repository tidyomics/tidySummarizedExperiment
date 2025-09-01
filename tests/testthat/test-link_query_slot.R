context("link_query_slot simple test")

library(tidySummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)

test_that("extract_col_data works correctly", {
    # Use testthat edition 3
    local_edition(3)
    
    # Create SE with realistic colData
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:12, nrow = 3)),
        rowData = DataFrame(gene_id = paste0("GENE", 1:3)),
        colData = DataFrame(
            sample_id = paste0("SAMPLE", 1:4),
            condition = c("A", "B", "A", "B")
        )
    )
    rownames(se) <- paste0("gene", 1:3)
    colnames(se) <- paste0("sample", 1:4)
    
    # Create realistic tibble from SE manipulation
    test_tibble <- tibble::tibble(
        .feature = rep(paste0("gene", 1:3), 4),
        .sample = rep(paste0("sample", 1:4), each = 3),
        counts = as.vector(matrix(1:12, nrow = 3)),
        sample_id = rep(paste0("SAMPLE", 1:4), each = 3),
        condition = rep(c("A", "B", "A", "B"), each = 3),
        new_sample_col = rep(c("X", "Y", "Z", "W"), each = 3)  # New sample data
    )
    
    # Test extract_col_data
    col_data <- extract_col_data(test_tibble, se)
    
    expect_s4_class(col_data, "DataFrame")
    expect_equal(rownames(col_data), paste0("sample", 1:4))
    expect_true("sample_id" %in% colnames(col_data))
    expect_true("condition" %in% colnames(col_data))
    expect_true("new_sample_col" %in% colnames(col_data))
})

test_that("extract_row_data works correctly", {
    local_edition(3)
    
    # Create SE with realistic rowData
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:12, nrow = 3)),
        rowData = DataFrame(
            gene_id = paste0("GENE", 1:3),
            gene_type = c("A", "B", "C")
        ),
        colData = DataFrame(sample_id = paste0("SAMPLE", 1:4))
    )
    rownames(se) <- paste0("gene", 1:3)
    colnames(se) <- paste0("sample", 1:4)
    
    # Create realistic tibble
    test_tibble <- tibble::tibble(
        .feature = rep(paste0("gene", 1:3), 4),
        .sample = rep(paste0("sample", 1:4), each = 3),
        counts = as.vector(matrix(1:12, nrow = 3)),
        gene_id = rep(paste0("GENE", 1:3), 4),
        gene_type = rep(c("A", "B", "C"), 4),
        new_gene_col = rep(c("X", "Y", "Z"), 4)  # New feature data
    )
    
    # Test extract_row_data
    row_data <- extract_row_data(test_tibble, se)
    
    expect_s4_class(row_data, "DataFrame")
    expect_equal(rownames(row_data), paste0("gene", 1:3))
    expect_true("gene_id" %in% colnames(row_data))
    expect_true("gene_type" %in% colnames(row_data))
    expect_true("new_gene_col" %in% colnames(row_data))
})

test_that("extract_assay_colnames works correctly", {
    local_edition(3)
    
    # Create minimal SE
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:12, nrow = 3)),
        rowData = DataFrame(gene_id = paste0("GENE", 1:3)),
        colData = DataFrame(sample_id = paste0("SAMPLE", 1:4))
    )
    rownames(se) <- paste0("gene", 1:3)
    colnames(se) <- paste0("sample", 1:4)
    
    # Create tibble with new assay data
    test_tibble <- tibble::tibble(
        .feature = rep(paste0("gene", 1:3), 4),
        .sample = rep(paste0("sample", 1:4), each = 3),
        counts = as.vector(matrix(1:12, nrow = 3)),  # Existing assay
        logcounts = as.vector(matrix(13:24, nrow = 3)),  # New assay
        sample_id = rep(paste0("SAMPLE", 1:4), each = 3),
        gene_id = rep(paste0("GENE", 1:3), 4)
    )
    
    # Extract other data first
    col_data <- extract_col_data(test_tibble, se)
    row_data <- extract_row_data(test_tibble, se, col_data = col_data)
    
    # Test extract_assay_colnames
    assay_colnames <- extract_assay_colnames(test_tibble, se, col_data, row_data)
    
    expect_type(assay_colnames, "character")
    expect_true("logcounts" %in% assay_colnames)
    expect_false("counts" %in% assay_colnames)  # Existing assay should be excluded
    expect_false(".feature" %in% assay_colnames)
    expect_false(".sample" %in% assay_colnames)
})

test_that("functions work together in realistic scenario", {
    local_edition(3)
    
    # Create SE similar to what tidySummarizedExperiment would have
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:12, nrow = 3)),
        rowData = DataFrame(
            gene_id = paste0("GENE", 1:3),
            chromosome = c("chr1", "chr2", "chr1")
        ),
        colData = DataFrame(
            sample_id = paste0("SAMPLE", 1:4),
            condition = c("treated", "untreated", "treated", "untreated")
        )
    )
    rownames(se) <- paste0("gene", 1:3)
    colnames(se) <- paste0("sample", 1:4)
    
    # Simulate what would come from a mutate operation
    mutated_tibble <- tibble::tibble(
        .feature = rep(paste0("gene", 1:3), 4),
        .sample = rep(paste0("sample", 1:4), each = 3),
        counts = as.vector(matrix(1:12, nrow = 3)),
        # Existing metadata
        gene_id = rep(paste0("GENE", 1:3), 4),
        chromosome = rep(c("chr1", "chr2", "chr1"), 4),
        sample_id = rep(paste0("SAMPLE", 1:4), each = 3),
        condition = rep(c("treated", "untreated", "treated", "untreated"), each = 3),
        # New columns from mutation
        log_counts = log2(as.vector(matrix(1:12, nrow = 3)) + 1),  # New assay
        batch = rep(c("A", "B", "A", "B"), each = 3),  # New colData
        gene_length = rep(c(1000, 2000, 1500), 4)  # New rowData
    )
    
    # Extract all data types
    col_data <- extract_col_data(mutated_tibble, se)
    row_data <- extract_row_data(mutated_tibble, se, col_data = col_data)
    assay_colnames <- extract_assay_colnames(mutated_tibble, se, col_data, row_data)
    
    # Verify correct separation
    expect_true("batch" %in% colnames(col_data))
    expect_false("gene_length" %in% colnames(col_data))
    
    expect_true("gene_length" %in% colnames(row_data))
    expect_false("batch" %in% colnames(row_data))
    
    expect_true("log_counts" %in% assay_colnames)
    expect_false("counts" %in% assay_colnames)
    expect_false("batch" %in% assay_colnames)
    expect_false("gene_length" %in% assay_colnames)
})

test_that("functions handle minimal SE objects", {
    local_edition(3)
    
    # Create minimal SE
    minimal_se <- SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2))
    )
    rownames(minimal_se) <- c("gene1", "gene2")
    colnames(minimal_se) <- c("sample1", "sample2", "sample3")
    
    # Create minimal tibble
    test_tibble <- tibble::tibble(
        .feature = rep(c("gene1", "gene2"), 3),
        .sample = rep(c("sample1", "sample2", "sample3"), each = 2),
        counts = 1:6,
        new_sample_data = rep(c("A", "B", "C"), each = 2)
    )
    
    # Test functions work even with minimal SE
    col_data <- extract_col_data(test_tibble, minimal_se)
    row_data <- extract_row_data(test_tibble, minimal_se)
    assay_colnames <- extract_assay_colnames(test_tibble, minimal_se, col_data, row_data)
    
    expect_s4_class(col_data, "DataFrame")
    expect_s4_class(row_data, "DataFrame")
    expect_type(assay_colnames, "character")
    
         # Should properly classify the new data
     expect_true("new_sample_data" %in% colnames(col_data))
     expect_equal(length(assay_colnames), 0)  # No new assays in this case
 })

test_that("update_SE_from_tibble works correctly", {
    local_edition(3)
    
    # Create original SE
    original_se <- SummarizedExperiment(
        assays = list(counts = matrix(1:12, nrow = 3)),
        rowData = DataFrame(
            gene_id = paste0("GENE", 1:3),
            gene_type = c("A", "B", "C")
        ),
        colData = DataFrame(
            sample_id = paste0("SAMPLE", 1:4),
            condition = c("treated", "untreated", "treated", "untreated")
        )
    )
    rownames(original_se) <- paste0("gene", 1:3)
    colnames(original_se) <- paste0("sample", 1:4)
    
    # Create mutated tibble with new data in all categories
    mutated_tibble <- tibble::tibble(
        .feature = rep(paste0("gene", 1:3), 4),
        .sample = rep(paste0("sample", 1:4), each = 3),
        counts = as.vector(matrix(1:12, nrow = 3)),
        # Existing metadata
        gene_id = rep(paste0("GENE", 1:3), 4),
        gene_type = rep(c("A", "B", "C"), 4),
        sample_id = rep(paste0("SAMPLE", 1:4), each = 3),
        condition = rep(c("treated", "untreated", "treated", "untreated"), each = 3),
        # New data - simpler approach
        new_assay = as.vector(matrix(13:24, nrow = 3)),  # New assay data
        batch = rep(c("X", "Y", "X", "Y"), each = 3),  # New colData
        gene_length = rep(c(1000, 2000, 1500), 4)  # New rowData
    )
    
    # Test update_SE_from_tibble
    updated_se <- update_SE_from_tibble(mutated_tibble, original_se)
    
    # Verify it returns a SummarizedExperiment
    expect_s4_class(updated_se, "SummarizedExperiment")
    
    # Verify dimensions unchanged
    expect_equal(dim(updated_se), dim(original_se))
    expect_equal(rownames(updated_se), rownames(original_se))
    expect_equal(colnames(updated_se), colnames(original_se))
    
    # Verify original assays preserved
    expect_true("counts" %in% assayNames(updated_se))
    expect_equal(assay(updated_se, "counts"), assay(original_se, "counts"))
    
    # Verify new assay added
    expect_true("new_assay" %in% assayNames(updated_se))
    expected_new_assay <- matrix(13:24, nrow = 3, 
                                dimnames = list(paste0("gene", 1:3), paste0("sample", 1:4)))
    expect_equal(assay(updated_se, "new_assay"), expected_new_assay)
    
    # Verify original colData preserved and new data added
    expect_true("sample_id" %in% colnames(colData(updated_se)))
    expect_true("condition" %in% colnames(colData(updated_se)))
    expect_true("batch" %in% colnames(colData(updated_se)))
    expect_equal(colData(updated_se)$batch, c("X", "Y", "X", "Y"))
    
    # Verify original rowData preserved and new data added
    expect_true("gene_id" %in% colnames(rowData(updated_se)))
    expect_true("gene_type" %in% colnames(rowData(updated_se)))
    expect_true("gene_length" %in% colnames(rowData(updated_se)))
    expect_equal(rowData(updated_se)$gene_length, c(1000, 2000, 1500))
})

test_that("update_SE_from_tibble handles edge cases", {
    local_edition(3)
    
    # Create minimal SE
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2))
    )
    rownames(se) <- c("gene1", "gene2")
    colnames(se) <- c("sample1", "sample2", "sample3")
    
    # Test with tibble that only adds assay data
    assay_only_tibble <- tibble::tibble(
        .feature = rep(c("gene1", "gene2"), 3),
        .sample = rep(c("sample1", "sample2", "sample3"), each = 2),
        counts = as.vector(matrix(1:6, nrow = 2)),
        new_assay = as.vector(matrix(7:12, nrow = 2))
    )
    
    updated_se <- update_SE_from_tibble(assay_only_tibble, se)
    
    # Should add new assay
    expect_true("new_assay" %in% assayNames(updated_se))
    expected_matrix <- matrix(7:12, nrow = 2,
                             dimnames = list(c("gene1", "gene2"), c("sample1", "sample2", "sample3")))
    expect_equal(assay(updated_se, "new_assay"), expected_matrix)
    
    # Should not affect colData/rowData (they should remain minimal)
    expect_equal(ncol(colData(updated_se)), ncol(colData(se)))
    expect_equal(ncol(rowData(updated_se)), ncol(rowData(se)))
})

test_that("update_SE_from_tibble works with column_belonging", {
    local_edition(3)
    
    # Create SE
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2)),
        colData = DataFrame(sample_id = c("S1", "S2", "S3"))
    )
    rownames(se) <- c("gene1", "gene2")
    colnames(se) <- c("sample1", "sample2", "sample3")
    
    # Create tibble with ambiguous column that should be forced to sample data
    test_tibble <- tibble::tibble(
        .feature = rep(c("gene1", "gene2"), 3),
        .sample = rep(c("sample1", "sample2", "sample3"), each = 2),
        counts = as.vector(matrix(1:6, nrow = 2)),
        sample_id = rep(c("S1", "S2", "S3"), each = 2),  # Include existing colData
        ambiguous_col = rep(c("A", "B", "C"), each = 2)
    )
    
    # Force ambiguous_col to be treated as sample data
    column_belonging <- c("ambiguous_col" = ".sample")
    
    updated_se <- update_SE_from_tibble(test_tibble, se, column_belonging = column_belonging)
    
    # ambiguous_col should be in colData
    expect_true("ambiguous_col" %in% colnames(colData(updated_se)))
    expect_equal(colData(updated_se)$ambiguous_col, c("A", "B", "C"))
    
    # Should not be in rowData
    expect_false("ambiguous_col" %in% colnames(rowData(updated_se)))
})

test_that("update_SE_from_tibble preserves SE class and metadata", {
    local_edition(3)
    
    # Create SE with metadata
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2)),
        colData = DataFrame(sample_id = c("S1", "S2", "S3")),
        rowData = DataFrame(gene_id = c("G1", "G2")),
        metadata = list(experiment_info = "test_experiment")
    )
    rownames(se) <- c("gene1", "gene2")
    colnames(se) <- c("sample1", "sample2", "sample3")
    
    # Create simple mutated tibble
    mutated_tibble <- tibble::tibble(
        .feature = rep(c("gene1", "gene2"), 3),
        .sample = rep(c("sample1", "sample2", "sample3"), each = 2),
        counts = as.vector(matrix(1:6, nrow = 2)),
        sample_id = rep(c("S1", "S2", "S3"), each = 2),  # Include existing colData
        gene_id = rep(c("G1", "G2"), 3),  # Include existing rowData
        new_data = as.vector(matrix(7:12, nrow = 2))
    )
    
    updated_se <- update_SE_from_tibble(mutated_tibble, se)
    
    # Should preserve original class
    expect_s4_class(updated_se, class(se)[1])
    
    # Should preserve metadata
    expect_equal(metadata(updated_se), metadata(se))
    
    # Should preserve original data structure
    expect_equal(rownames(updated_se), rownames(se))
    expect_equal(colnames(updated_se), colnames(se))
})
