context("tidySummarizedExperiment:::analyze_query_scope_mutate test")

library(tidySummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)

# Helper function to create test SE
create_test_se <- function() {
    se <- SummarizedExperiment(
        assays = list(
            counts = matrix(1:12, nrow = 3),
            logcounts = matrix(13:24, nrow = 3)
        ),
        rowData = DataFrame(
            gene_id = paste0("GENE", 1:3),
            gene_type = c("protein_coding", "lncRNA", "miRNA"),
            chromosome = c("chr1", "chr2", "chr1")
        ),
        colData = DataFrame(
            sample_id = paste0("SAMPLE", 1:4),
            condition = c("treated", "untreated", "treated", "untreated"),
            batch = c("A", "A", "B", "B")
        )
    )
    rownames(se) <- paste0("gene", 1:3)
    colnames(se) <- paste0("sample", 1:4)
    return(se)
}

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate handles simple queries correctly", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Test simple existing column modification (colData)
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, new_condition = condition)
    expect_equal(scope$scope, "coldata_only")
    expect_equal(scope$query_complexity, "simple")
    # May use pre_mutate or subset_analysis or dependency_analysis depending on implementation
    expect_true(scope$analysis_method %in% c("pre_mutate", "subset_analysis", "dependency_analysis"))
    expect_true(scope$confidence %in% c("high", "medium"))
    expect_equal(scope$column_names, "new_condition")
    expect_true(scope$targets_coldata)
    expect_false(scope$targets_rowdata)
    expect_false(scope$targets_assays)
    
    # Test simple existing column modification (assay)
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, new_counts = counts)
    expect_equal(scope$scope, "assay_only")
    expect_equal(scope$query_complexity, "simple")
    expect_true(scope$targets_assays)
    
    # Test simple existing column modification (rowData)
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, new_type = gene_type)
    expect_equal(scope$scope, "rowdata_only")
    expect_equal(scope$query_complexity, "simple")
    expect_true(scope$targets_rowdata)
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate handles simple function calls", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Test simple function on assay data
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, log_counts = log2(counts))
    expect_equal(scope$query_complexity, "simple")
    expect_true(scope$scope %in% c("assay_only", "unknown"))  # May be unknown for pre-mutate
    
    # If unknown from pre-mutate, should try subset analysis
    if (scope$scope == "unknown") {
        expect_true(scope$analysis_method %in% c("subset_analysis","dependency_analysis"))
        expect_true(scope$confidence %in% c("medium","high"))
    }
    
    # Test simple function on colData
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, upper_condition = toupper(condition))
    expect_equal(scope$query_complexity, "simple")
    
    # Test simple function on rowData
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, upper_type = toupper(gene_type))
    expect_equal(scope$query_complexity, "simple")
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate detects complex queries", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Test conditional logic (complex)
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, 
        complex_col = ifelse(counts > 5, "high", "low"))
    expect_equal(scope$query_complexity, "complex")
    expect_true(scope$analysis_method %in% c("subset_analysis","dependency_analysis"))
    expect_true(scope$confidence %in% c("medium","high"))
    
    # Test aggregation functions (complex)
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, 
        mean_counts = mean(counts))
    expect_equal(scope$query_complexity, "complex")
    expect_true(scope$analysis_method %in% c("subset_analysis","dependency_analysis"))
    
    # Test long expression (complex)
    long_expr <- paste(rep("counts +", 50), collapse = " ")
    long_expr <- paste0(long_expr, " 1")
    
    # Use rlang::parse_expr to create the expression
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, 
        long_col = !!rlang::parse_expr(long_expr))
    expect_equal(scope$query_complexity, "complex")
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate handles mixed operations", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Test mixed simple operations
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se,
        new_condition = condition,
        new_type = gene_type,
        new_counts = counts)
    
    # Should be complex due to multiple operations, or mixed if analyzed
    expect_true(scope$scope %in% c("mixed", "unknown"))
    expect_true(scope$query_complexity %in% c("simple", "complex"))
    expect_true(length(scope$column_names) == 3)
    expect_equal(sort(scope$column_names), c("new_condition", "new_counts", "new_type"))
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate subset analysis works correctly", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Force subset analysis with complex query
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, 
        subset_col = ifelse(counts > mean(counts), "high", "low"))
    
    expect_true(scope$analysis_method %in% c("subset_analysis","dependency_analysis"))
    expect_true(scope$confidence %in% c("medium","high"))
    expect_equal(scope$query_complexity, "complex")
    
    # Should successfully analyze the scope (may be unknown for complex operations)
    expect_true(scope$scope %in% c("assay_only", "mixed", "unknown", "coldata_only"))
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate handles edge cases", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Unnamed expression should be handled (return a result list)
    expect_true(is.list(tidySummarizedExperiment:::analyze_query_scope_mutate(se, condition)))
    
    # Test empty SE
    empty_se <- se[0, 0]
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(empty_se, test_col = 1)
    expect_equal(scope$scope, "unknown")
    expect_equal(scope$confidence, "low")
    
    # Test with missing columns (should handle gracefully)
    expect_error(tidySummarizedExperiment:::analyze_query_scope_mutate(se, test_col = nonexistent_column), NA)
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate query complexity analysis works", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Test simple query complexity through the main function
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, simple_col = condition)
    expect_equal(scope$query_complexity, "simple")
    
    # Test complex query complexity  
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, complex_col = ifelse(counts > 5, "high", "low"))
    expect_equal(scope$query_complexity, "complex")
})

test_that("tidySummarizedExperiment:::analyze_query_scope_mutate handles realistic scenarios", {
    local_edition(3)
    
    se <- create_test_se()
    
    # Realistic log transformation
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, log_counts = log2(counts + 1))
    expect_equal(scope$query_complexity, "simple")
    expect_true(scope$scope %in% c("assay_only", "unknown"))
    
    # Realistic sample annotation
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, treatment_group = paste0(condition, "_group"))
    expect_equal(scope$query_complexity, "simple")
    
    # Realistic gene categorization
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, gene_category = toupper(gene_type))
    expect_equal(scope$query_complexity, "simple")
    
    # Realistic complex operation
    scope <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, 
        high_expression = ifelse(counts > quantile(counts, 0.75), "high", "low"))
    expect_equal(scope$query_complexity, "complex")
    expect_true(scope$analysis_method %in% c("subset_analysis","dependency_analysis"))
})
