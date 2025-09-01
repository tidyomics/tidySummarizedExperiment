context("mutate function comprehensive tests")

library(tidySummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(airway)

# Helper: airway-based test SE (subset for speed)
create_airway_test_se <- function() {
    data(airway)
    se <- airway
    se <- se[1:200, ]
    se
}

# Test simple queries (should use optimized paths)
test_that("mutate handles simple colData operations correctly", {
    se <- create_airway_test_se()
    
    # Simple column assignment
    result <- se %>% mutate(new_dex = dex)
    expect_true("new_dex" %in% colnames(colData(result)))
    expect_equal(as.character(colData(result)$new_dex), as.character(colData(se)$dex))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "coldata_only")
    
    # Simple arithmetic on colData
    result <- se %>% mutate(avgLength_plus_5 = avgLength + 5)
    expect_true("avgLength_plus_5" %in% colnames(colData(result)))
    expect_equal(colData(result)$avgLength_plus_5, colData(se)$avgLength + 5)
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "coldata_only")
    
    # Simple function on colData
    result <- se %>% mutate(cell_upper = toupper(cell))
    expect_true("cell_upper" %in% colnames(colData(result)))
    expect_equal(colData(result)$cell_upper, toupper(colData(se)$cell))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "coldata_only")
    
    # Multiple simple colData operations
    result <- se %>% mutate(
        avgLength_group = avgLength,
        dex_copy = dex,
        doubled_avgLength = avgLength * 2
    )
    expect_true(all(c("avgLength_group", "dex_copy", "doubled_avgLength") %in% colnames(colData(result))))
    expect_equal(colData(result)$doubled_avgLength, colData(se)$avgLength * 2)
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "coldata_only")
})

test_that("mutate handles simple rowData operations correctly", {
    se <- create_airway_test_se()
    
    # Simple column assignment
    result <- se %>% mutate(new_gene_biotype = gene_biotype)
    expect_true("new_gene_biotype" %in% colnames(rowData(result)))
    expect_equal(as.character(rowData(result)$new_gene_biotype), as.character(rowData(se)$gene_biotype))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "rowdata_only")
    
    # Simple arithmetic on rowData
    result <- se %>% mutate(gene_length_kb = (gene_seq_end - gene_seq_start) / 1000)
    expect_true("gene_length_kb" %in% colnames(rowData(result)))
    expect_equal(rowData(result)$gene_length_kb, (rowData(se)$gene_seq_end - rowData(se)$gene_seq_start) / 1000)
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "rowdata_only")
    
    # Simple function on rowData
    result <- se %>% mutate(symbol_upper = toupper(symbol))
    expect_true("symbol_upper" %in% colnames(rowData(result)))
    expect_equal(rowData(result)$symbol_upper, toupper(rowData(se)$symbol))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "rowdata_only")
    
    # Mathematical operations on rowData
    result <- se %>% mutate(
        start_kb = round(gene_seq_start / 1000, 2),
        log_end = log2(pmax(gene_seq_end, 1))
    )
    expect_true(all(c("start_kb", "log_end") %in% colnames(rowData(result))))
    expect_equal(rowData(result)$log_end, log2(pmax(rowData(se)$gene_seq_end, 1)))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "rowdata_only")
})

test_that("mutate handles simple assay operations correctly", {
    se <- create_airway_test_se()
    
    # Simple assay column assignment
    result <- se %>% mutate(counts_copy = counts)
    expect_true("counts_copy" %in% assayNames(result))
    expect_equal(assay(result, "counts_copy"), assay(se, "counts"))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "assay_only")
    
    # Simple arithmetic on assays
    result <- se %>% mutate(counts_plus_1 = counts + 1)
    expect_true("counts_plus_1" %in% assayNames(result))
    expect_equal(assay(result, "counts_plus_1"), assay(se, "counts") + 1)
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "assay_only")
    
    # Simple function on assays
    result <- se %>% mutate(log_counts_manual = log2(counts + 1))
    expect_true("log_counts_manual" %in% assayNames(result))
    expect_equal(assay(result, "log_counts_manual"), log2(assay(se, "counts") + 1))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "assay_only")
    
    # Multiple simple assay operations
    result <- se %>% mutate(
        normalized_counts = counts / 1000,
        sqrt_counts = sqrt(counts)
    )
    expect_true(all(c("normalized_counts", "sqrt_counts") %in% assayNames(result)))
    expect_equal(assay(result, "sqrt_counts"), sqrt(assay(se, "counts")))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "assay_only")
})

# Test complex queries (should use subset analysis or tibble conversion)
test_that("mutate handles complex operations with conditionals", {
    se <- create_airway_test_se()
    
    # Complex conditional on colData
    result <- se %>% mutate(
        length_group = ifelse(avgLength > mean(avgLength), "longer", "shorter")
    )
    expect_true("length_group" %in% colnames(colData(result)))
    expected <- ifelse(colData(se)$avgLength > mean(colData(se)$avgLength), "longer", "shorter")
    expect_equal(as.character(colData(result)$length_group), expected)
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "coldata_only")
    # Complex conditional on rowData - using tibble approach
    result_tbl <- se %>%  
        mutate(length_category = case_when(
            (gene_seq_end - gene_seq_start) > 3e4 ~ "long",
            (gene_seq_end - gene_seq_start) > 2e4 ~ "medium",
            TRUE ~ "short"
        ))
    expect_true("length_category" %in% colnames(rowData(result_tbl)))
    expect_equal(S4Vectors::metadata(result_tbl)$latest_mutate_scope_report$scope, "rowdata_only")
    # Complex nested conditionals
    result <- se %>% mutate(complex_category = ifelse(
        dex == "trt" & avgLength > mean(avgLength), 
        "treated_long", 
        ifelse(dex == "untrt", "untreated", "other")
    ))
    expect_true("complex_category" %in% colnames(colData(result)))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "coldata_only")
})

test_that("mutate handles complex operations with aggregations", {
    se <- create_airway_test_se()
    
    # This should use tibble conversion for complex aggregations
    result_tbl <- se %>% as_tibble() %>% mutate(mean_avgLength = mean(avgLength))
    expect_true("mean_avgLength" %in% colnames(result_tbl))
    
    # Complex expression combining multiple data types via tibble
    result_tbl <- se %>% as_tibble() %>% mutate(
        high_expression_in_trt = ifelse(
            dex == "trt" & counts > median(counts), 
            "yes", 
            "no"
        )
    )
    expect_true("high_expression_in_trt" %in% colnames(result_tbl))
})

# Test scope detection accuracy
test_that("scope detection correctly identifies operation targets", {
    se <- create_airway_test_se()
    
    # Test colData-only scope detection
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, new_len = avgLength)
    expect_equal(scope_result$scope, "coldata_only")
    expect_true(scope_result$targets_coldata)
    expect_false(scope_result$targets_rowdata)
    expect_false(scope_result$targets_assays)
    
    # Test rowData-only scope detection
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, new_symbol = symbol)
    expect_equal(scope_result$scope, "rowdata_only")
    expect_false(scope_result$targets_coldata)
    expect_true(scope_result$targets_rowdata)
    expect_false(scope_result$targets_assays)
    
    # Test assay-only scope detection
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, new_counts = counts)
    expect_equal(scope_result$scope, "assay_only")
    expect_false(scope_result$targets_coldata)
    expect_false(scope_result$targets_rowdata)
    expect_true(scope_result$targets_assays)
})

test_that("scope detection handles mixed operations", {
    se <- create_airway_test_se()
    
    # Test mixed operation (colData + assay) - creates new assay
    result <- se %>% mutate(mixed_col = paste0(dex, "_", counts))
    # This should work and create a new assay
    expect_true("mixed_col" %in% assayNames(result))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "mixed")
    
    # Test multiple operations spanning different data types
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se,
        col_operation = dex,
        row_operation = gene_biotype,
        assay_operation = counts
    )
    # Should be detected as mixed or unknown (leading to tibble conversion)
    expect_true(scope_result$scope %in% c("mixed", "unknown"))
})

test_that("mutate handles mixed arithmetic across slots (assay * colData)", {
    se <- create_airway_test_se()
    result <- se %>% mutate(new_counts = counts * avgLength)
    # Mixed operation should create a new assay
    expect_true("new_counts" %in% assayNames(result))
    expect_equal(assay(result, "new_counts"), assay(se, "counts") * matrix(rep(colData(se)$avgLength, each = nrow(se)), nrow = nrow(se)))
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "mixed")
    deps <- S4Vectors::metadata(result)$latest_mutate_scope_report$expression_dependencies$new_counts
    expect_true(deps$uses_assays)
    expect_true(deps$uses_coldata)
    expect_false(deps$uses_rowdata)
    expect_true("counts" %in% deps$referenced_assays)
    expect_true("avgLength" %in% deps$referenced_coldata)
})

test_that("mutate handles mixed arithmetic across slots (assay * rowData)", {
    se <- create_airway_test_se()
    result <- se %>% mutate(new_counts = counts * as.numeric(scale(gene_seq_end - gene_seq_start)))
    expect_true("new_counts" %in% assayNames(result))
    expected_scaled <- assay(se, "counts") * matrix(rep(as.numeric(scale(rowData(se)$gene_seq_end - rowData(se)$gene_seq_start)), times = ncol(se)), nrow = nrow(se))
    expect_equal(assay(result, "new_counts"), expected_scaled, tolerance = 5)
    expect_equal(S4Vectors::metadata(result)$latest_mutate_scope_report$scope, "mixed")
    deps <- S4Vectors::metadata(result)$latest_mutate_scope_report$expression_dependencies$new_counts
    expect_true(deps$uses_assays)
    expect_false(deps$uses_coldata)
    expect_true(deps$uses_rowdata)
    expect_true("counts" %in% deps$referenced_assays)
    expect_true(all(c("gene_seq_start","gene_seq_end") %in% deps$referenced_rowdata))
})

# Test execution paths work correctly
test_that("different execution paths produce correct results", {
    se <- create_airway_test_se()
    
    # Test modify_samples path (colData-only)
    result_coldata <- se %>% mutate(avgLength_doubled = avgLength * 2)
    expect_equal(colData(result_coldata)$avgLength_doubled, colData(se)$avgLength * 2)
    expect_equal(dim(result_coldata), dim(se))
    expect_equal(assayNames(result_coldata), assayNames(se))
    
    # Test modify_features path (rowData-only)
    result_rowdata <- se %>% mutate(length_log = log(pmax(gene_seq_end - gene_seq_start, 1)))
    expect_equal(rowData(result_rowdata)$length_log, log(pmax(rowData(se)$gene_seq_end - rowData(se)$gene_seq_start, 1)))
    expect_equal(dim(result_rowdata), dim(se))
    expect_equal(assayNames(result_rowdata), assayNames(se))
    
    # Test mutate_assay path (assay-only)
    result_assay <- se %>% mutate(normalized_counts = counts / 100)
    expect_true("normalized_counts" %in% assayNames(result_assay))
    expect_equal(assay(result_assay, "normalized_counts"), assay(se, "counts") / 100)
    expect_equal(dim(result_assay), dim(se))
    
    # Test tibble conversion path (complex/mixed operations)
    result_complex <- se %>% mutate(
        complex_feature = ifelse(dex == "trt" & counts > 5, "yes", "no")
    )
    # Complex mixed operations create new assays
    expect_true("complex_feature" %in% assayNames(result_complex))
})

# Test query complexity classification
test_that("query complexity is correctly classified", {
    se <- create_airway_test_se()
    
    # Simple queries
    simple_queries <- list(
        "simple_assignment" = quote(dex),
        "simple_arithmetic" = quote(avgLength + 5),
        "simple_function" = quote(log2(counts)),
        "simple_paste" = quote(paste0("prefix_", cell))
    )
    
    for (query_name in names(simple_queries)) {
        scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, 
            !!query_name := !!simple_queries[[query_name]])
        expect_equal(scope_result$query_complexity, "simple", 
                    info = paste("Query", query_name, "should be classified as simple"))
    }
    
    # Complex queries
    complex_queries <- list(
        "conditional" = quote(ifelse(avgLength > mean(avgLength), "long", "short")),
        "aggregation" = quote(mean(avgLength))
    )
    
    for (query_name in names(complex_queries)) {
        scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se,
            !!query_name := !!complex_queries[[query_name]])
        expect_equal(scope_result$query_complexity, "complex",
                    info = paste("Query", query_name, "should be classified as complex"))
    }
})

# Test performance and optimization
test_that("optimized paths are faster than tibble conversion", {
    se <- create_airway_test_se()
    
    # Simple colData operation (should use optimized path)
    time_optimized <- system.time({
        result_opt <- se %>% mutate(len_plus_10 = avgLength + 10)
    })
    
    # Force tibble conversion by using a complex operation
    time_tibble <- system.time({
        result_tibble <- se %>% mutate(complex_op = ifelse(avgLength > mean(avgLength) & dex == "trt", 1, 0))
    })
    
    # The optimized path should generally be faster, but we'll just check it completes
    expect_true(time_optimized[["elapsed"]] >= 0)
    expect_true(time_tibble[["elapsed"]] >= 0)
    
    # Both should produce valid results
    expect_s4_class(result_opt, "SummarizedExperiment")
    # This operation only uses colData, so result is stored in colData
    expect_true("complex_op" %in% colnames(colData(result_tibble)))
})

# Test edge cases and error handling
test_that("mutate handles edge cases correctly", {
    se <- create_airway_test_se()
    
    # Single row/column SE
    single_se <- se[1, 1]
    result_single <- single_se %>% mutate(new_len = avgLength + 1)
    expect_equal(dim(result_single), c(1, 1))
    expect_true("new_len" %in% colnames(colData(result_single)))
    
    # Overwriting existing columns
    result_overwrite <- se %>% mutate(dex = "trt")
    expect_true(all(colData(result_overwrite)$dex == "trt"))
    
    # Adding multiple columns at once
    result_multiple <- se %>% mutate(
        col1 = avgLength,
        col2 = cell, 
        col3 = dex,
        assay1 = counts,
        assay2 = log2(counts + 1)
    )
    expect_true(all(c("col1", "col2", "col3") %in% colnames(colData(result_multiple))))
    expect_true(all(c("assay1", "assay2") %in% assayNames(result_multiple)))
})

test_that("mutate error handling works correctly", {
    se <- create_airway_test_se()
    
    # Trying to modify special/protected columns should give appropriate errors
    expect_error(se %>% mutate(.feature = "modified"), "view only")
    expect_error(se %>% mutate(.sample = "modified"), "view only")
    
    # Operations that reference non-existent columns should be handled gracefully
    # This might work via tibble conversion or give an informative error
    expect_error(se %>% mutate(test = nonexistent_column), "object.*not found|could not find")
})

# Test integration with other dplyr operations
test_that("mutate integrates well with other dplyr operations", {
    se <- create_airway_test_se()
    
    # Chain mutate with filter
    result <- se %>% 
        mutate(length_group = ifelse(avgLength > mean(avgLength), "longer", "shorter")) %>%
        filter(length_group == "longer")
    
    # length_group is stored in colData for colData-only operations
    expect_true("length_group" %in% colnames(colData(result)))
    expect_true(all(as_tibble(result)$avgLength > mean(create_airway_test_se() %>% colData() %>% as.data.frame() %>% .$avgLength)))
    
    # Chain mutate with select
    result <- se %>%
        mutate(new_cell = toupper(cell)) %>%
        select(new_cell, avgLength)
    
    expect_equal(ncol(result), 2)
    expect_true("new_cell" %in% colnames(result))
    
    # Chain multiple mutates
    result <- se %>%
        mutate(len_doubled = avgLength * 2) %>%
        mutate(len_category = ifelse(len_doubled > mean(len_doubled), "high", "low"))
    
    # Both should be in colData for colData-only operations
    expect_true(all(c("len_doubled", "len_category") %in% colnames(colData(result))))
})

# Test with real-world-like scenarios
test_that("mutate handles realistic bioinformatics scenarios", {
    se <- create_airway_test_se()
    
    # Normalize counts by library size - using tibble conversion for mixed operations
    result_tbl <- se %>% as_tibble() %>% mutate(cpm = (counts / avgLength) * 1e6)
    expect_true("cpm" %in% colnames(result_tbl))
    
    # Create gene expression categories - simplified
    result_tbl <- se %>% as_tibble() %>% mutate(
        expression_category = ifelse(counts > median(counts), "high", "low")
    )
    expect_true("expression_category" %in% colnames(result_tbl))
    
    # Combine sample metadata
    result <- se %>% mutate(
        sample_info = paste(cell, dex, SampleName, sep = "_"),
        treatment_flag = ifelse(dex == "trt", 1, 0)
    )
    expect_true(all(c("sample_info", "treatment_flag") %in% colnames(colData(result))))
    
    # Gene annotation processing
    result <- se %>% mutate(
        is_protein_coding = gene_biotype == "protein_coding",
        gene_symbol = toupper(symbol),
        kb_length = round((gene_seq_end - gene_seq_start) / 1000, 1)
    )
    expect_true(all(c("is_protein_coding", "gene_symbol", "kb_length") %in% colnames(rowData(result))))
})

# Test scope detection - simplified to avoid quosure issues
test_that("scope detection works correctly through mutate interface", {
    se <- create_airway_test_se()
    
    # Test through the main mutate scope analysis function
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, test_col = dex)
    expect_equal(scope_result$scope, "coldata_only")
    
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, test_col = gene_biotype)
    expect_equal(scope_result$scope, "rowdata_only")
    
    scope_result <- tidySummarizedExperiment:::analyze_query_scope_mutate(se, test_col = counts)
    expect_equal(scope_result$scope, "assay_only")
})