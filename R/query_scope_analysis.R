

#' Analyze query complexity to determine analysis strategy
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param dots List of quosures from mutate expressions
#' @return List with is_simple (logical) and complexity_reasons (character vector)
#' @noRd
analyze_query_complexity <- function(se, dots) {
  
  complexity_reasons <- character(0)
  
  # Check for complex expressions
  for (i in seq_along(dots)) {
    expr_text <- rlang::quo_text(dots[[i]])
    
    # Simple assignment from existing column (variable name only)
    if (grepl("^[a-zA-Z][a-zA-Z0-9_\\.]*$", expr_text)) {
      next  # Simple column reference
    }
    
    # Simple arithmetic with constants
    if (grepl("^[a-zA-Z][a-zA-Z0-9_\\.]*\\s*[+\\-\\*/]\\s*[0-9.]+$", expr_text)) {
      next  # Simple arithmetic
    }
    
    # Simple function calls that are likely fast
    simple_functions <- c("log", "log2", "log10", "sqrt", "abs", "round", "ceiling", "floor",
                         "toupper", "tolower", "paste0", "paste", "as.character", "as.numeric")
    simple_pattern <- paste0("^(", paste(simple_functions, collapse = "|"), ")\\(")
    
    if (grepl(simple_pattern, expr_text)) {
      next  # Simple function call
    }
    
    # Check for complexity indicators
    if (grepl("ifelse|case_when|switch", expr_text)) {
      complexity_reasons <- c(complexity_reasons, "conditional_logic")
    }
    
    if (grepl("sum\\(|mean\\(|median\\(|max\\(|min\\(|sd\\(|var\\(", expr_text)) {
      complexity_reasons <- c(complexity_reasons, "aggregation_functions")
    }
    
    if (grepl("group_by|summarise|summarize", expr_text)) {
      complexity_reasons <- c(complexity_reasons, "grouping_operations")
    }
    
    if (grepl("\\[|\\$", expr_text)) {
      complexity_reasons <- c(complexity_reasons, "subsetting_operations")
    }
    
    # Long expressions (likely complex)
    if (nchar(expr_text) > 100) {
      complexity_reasons <- c(complexity_reasons, "long_expression")
    }
  }
  
  is_simple <- length(complexity_reasons) == 0
  
  return(list(
    is_simple = is_simple,
    complexity_reasons = complexity_reasons
  ))
}



#' Analyze expression dependencies to detect mixed scope operations
#'
#' This function examines the variables referenced in mutate expressions
#' to determine if they span multiple data types (assays, colData, rowData).
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param dots List of quosures from mutate expressions
#' @return List with dependency analysis results
#' @noRd
analyze_expression_dependencies <- function(se, dots) {
  
  # Get available column names from each data type
  coldata_cols <- colnames(colData(se))
  rowdata_cols <- if (.hasSlot(se, "rowData") || .hasSlot(se, "elementMetadata")) {
    colnames(rowData(se))
  } else {
    character(0)
  }
  assay_cols <- assayNames(se)
  
  # Analyze each expression
  expression_deps <- list()
  
  for (i in seq_along(dots)) {
    expr_name <- names(dots)[i]
    if (is.null(expr_name) || expr_name == "") {
      expr_name <- paste0("expr_", i)
    }
    
    # Extract variable names from the expression
    expr_text <- rlang::quo_text(dots[[i]])
    
    # Find which variables are referenced
    # Use a simple approach: check which column names appear as whole words
    uses_coldata <- if (length(coldata_cols) > 0) {
      any(sapply(coldata_cols, function(col) {
        grepl(paste0("\\b", col, "\\b"), expr_text)
      }, USE.NAMES = FALSE))
    } else {
      FALSE
    }
    
    uses_rowdata <- if (length(rowdata_cols) > 0) {
      any(sapply(rowdata_cols, function(col) {
        grepl(paste0("\\b", col, "\\b"), expr_text)
      }, USE.NAMES = FALSE))
    } else {
      FALSE
    }
    
    uses_assays <- if (length(assay_cols) > 0) {
      any(sapply(assay_cols, function(col) {
        grepl(paste0("\\b", col, "\\b"), expr_text)
      }, USE.NAMES = FALSE))
    } else {
      FALSE
    }
    
    # Count dependency types
    n_dep_types <- sum(uses_coldata, uses_rowdata, uses_assays)
    
    # Determine scope for this expression
    expr_scope <- if (n_dep_types == 0) {
      "unknown"
    } else if (n_dep_types == 1) {
      if (uses_coldata) "coldata_only"
      else if (uses_rowdata) "rowdata_only"
      else "assay_only"
    } else {
      "mixed"
    }
    
    expression_deps[[expr_name]] <- list(
      scope = expr_scope,
      uses_coldata = uses_coldata,
      uses_rowdata = uses_rowdata,
      uses_assays = uses_assays,
      expression_text = expr_text,
      referenced_coldata = if (length(coldata_cols) > 0) coldata_cols[sapply(coldata_cols, function(col) grepl(paste0("\\b", col, "\\b"), expr_text), USE.NAMES = FALSE)] else character(0),
      referenced_rowdata = if (length(rowdata_cols) > 0) rowdata_cols[sapply(rowdata_cols, function(col) grepl(paste0("\\b", col, "\\b"), expr_text), USE.NAMES = FALSE)] else character(0),
      referenced_assays = if (length(assay_cols) > 0) assay_cols[sapply(assay_cols, function(col) grepl(paste0("\\b", col, "\\b"), expr_text), USE.NAMES = FALSE)] else character(0)
    )
  }
  
  # Determine overall scope based on all expressions
  all_scopes <- sapply(expression_deps, function(x) x$scope)
  has_mixed <- any(all_scopes == "mixed")
  unique_scopes <- unique(all_scopes[all_scopes != "unknown"])
  
  overall_scope <- if (has_mixed || length(unique_scopes) > 1) {
    "mixed"
  } else if (length(unique_scopes) == 1) {
    unique_scopes[1]
  } else {
    "unknown"
  }
  
  return(list(
    overall_scope = overall_scope,
    expression_dependencies = expression_deps,
    analysis_method = "dependency_analysis",
    confidence = "high"
  ))
}

#' Analyze query scope by running on a subset of the data
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param dots List of quosures from mutate expressions
#' @return Analysis result similar to other analyze functions
#' @noRd
analyze_query_scope_subset <- function(se, dots) {
  
  # Create a small subset for fast analysis
  # Take first few rows and columns to keep it fast
  n_rows <- min(10, nrow(se))
  n_cols <- min(5, ncol(se))
  
  if (n_rows == 0 || n_cols == 0) {
    # Empty SE - cannot analyze
    return(list(
      scope = "unknown",
      targets_coldata = FALSE,
      targets_rowdata = FALSE,
      targets_assays = FALSE,
      new_coldata_cols = character(0),
      new_rowdata_cols = character(0),
      new_assay_cols = character(0),
      analysis_method = "subset_analysis",
      confidence = "low"
    ))
  }
  
  se_subset <- se[1:n_rows, 1:n_cols]
  
  # Convert to tibble and try the mutate operation
  tryCatch({
    # Convert subset to tibble
    tibble_subset <- se_subset %>% as_tibble()
    
    # Apply mutate expressions
    mutated_subset <- tibble_subset %>% dplyr::mutate(!!!dots)
    
    # Extract column names that were created/modified
    .cols <- names(dots)
    
    # Use post-operation analysis on the result
    result <- analyze_query_scope_post_operation(se_subset, .cols, mutated_subset, "mutate")
    result$analysis_method <- "subset_analysis"
    result$confidence <- "medium"  # Medium confidence due to subset
    
    return(result)
    
  }, error = function(e) {
    # If mutate fails on subset, return unknown
    return(list(
      scope = "unknown",
      targets_coldata = FALSE,
      targets_rowdata = FALSE,
      targets_assays = FALSE,
      new_coldata_cols = character(0),
      new_rowdata_cols = character(0),
      new_assay_cols = character(0),
      analysis_method = "subset_analysis",
      confidence = "low",
      error_message = as.character(e$message)
    ))
  })
}


#' Pre-mutate analysis: Analyze query scope based on column names only
#'
#' @keywords internal
#' @noRd
analyze_query_scope_pre_mutate <- function(se, .cols) {
  
  # Get existing column categories
  existing_coldata_cols <- colnames(colData(se))
  existing_rowdata_cols <- if (.hasSlot(se, "rowData") || .hasSlot(se, "elementMetadata")) {
    colnames(rowData(se))
  } else {
    c()
  }
  existing_assay_cols <- names(assays(se))
  
  # Analyze which existing categories are being modified
  modifies_coldata <- any(.cols %in% existing_coldata_cols)
  modifies_rowdata <- any(.cols %in% existing_rowdata_cols)
  modifies_assays <- any(.cols %in% existing_assay_cols)
  
  # Count modifications
  n_modifications <- sum(modifies_coldata, modifies_rowdata, modifies_assays)
  
  # For new columns, we cannot reliably determine their target without actual data
  new_cols <- setdiff(.cols, c(existing_coldata_cols, existing_rowdata_cols, existing_assay_cols))
  
  # We can only be certain about existing columns being modified
  # New columns are ambiguous without seeing the actual data
  ambiguous_cols <- new_cols
  
  # Determine scope - only based on existing column modifications
  targets_coldata <- modifies_coldata
  targets_rowdata <- modifies_rowdata  
  targets_assays <- modifies_assays
  
  # Count target types
  n_targets <- sum(targets_coldata, targets_rowdata, targets_assays)
  
  scope <- if (n_targets == 0 && length(ambiguous_cols) > 0) {
    "unknown"
  } else if (n_targets == 1) {
    if (targets_coldata) "coldata_only"
    else if (targets_rowdata) "rowdata_only"
    else "assay_only"
  } else if (n_targets > 1) {
    "mixed"
  } else {
    "unknown"
  }
  
  return(list(
    scope = scope,
    targets_coldata = targets_coldata,
    targets_rowdata = targets_rowdata,
    targets_assays = targets_assays,
    new_coldata_cols = character(0),  # Cannot determine without actual data
    new_rowdata_cols = character(0),  # Cannot determine without actual data
    new_assay_cols = character(0),    # Cannot determine without actual data
    ambiguous_cols = ambiguous_cols,
    analysis_method = "pre_mutate",
    confidence = if (length(ambiguous_cols) == 0) "high" else "low"
  ))
}

#' Analyze query scope based on modified tibble (works for any dplyr operation)
#'
#' This function can analyze the scope of any dplyr operation that modifies
#' columns by comparing the original SE with the resulting tibble.
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param .cols Character vector of column names that were modified/created
#' @param .data_modified The resulting tibble after dplyr operation
#' @param operation Character indicating the operation type (for documentation)
#' @noRd
analyze_query_scope_post_operation <- function(se, .cols, .data_modified, operation = "unknown") {
  
  # Use the existing functions to classify columns
  col_data <- extract_col_data(.data_modified, se)
  row_data <- extract_row_data(.data_modified, se, col_data = col_data)
  colnames_assay <- extract_assay_colnames(.data_modified, se, col_data, row_data)
  
  # Identify which of the specified columns fall into each category
  new_coldata_cols <- intersect(.cols, colnames(col_data))
  new_rowdata_cols <- intersect(.cols, colnames(row_data))
  new_assay_cols <- intersect(.cols, colnames_assay)
  
  # Determine targets
  targets_coldata <- length(new_coldata_cols) > 0
  targets_rowdata <- length(new_rowdata_cols) > 0
  targets_assays <- length(new_assay_cols) > 0
  
  # Count target types
  n_targets <- sum(targets_coldata, targets_rowdata, targets_assays)
  
  # Determine scope
  scope <- if (n_targets == 0) {
    "unknown"
  } else if (n_targets == 1) {
    if (targets_coldata) "coldata_only"
    else if (targets_rowdata) "rowdata_only"
    else "assay_only"
  } else {
    "mixed"
  }
  
  return(list(
    scope = scope,
    targets_coldata = targets_coldata,
    targets_rowdata = targets_rowdata,
    targets_assays = targets_assays,
    new_coldata_cols = new_coldata_cols,
    new_rowdata_cols = new_rowdata_cols,
    new_assay_cols = new_assay_cols,
    analysis_method = paste0("post_", operation),
    operation = operation,
    confidence = "high"
  ))
}




#' Apply any dplyr operation to feature metadata (rowData)
#'
#' This is a general function that allows applying any dplyr operation
#' directly to the rowData of a SummarizedExperiment object, providing
#' better performance when operations only target feature metadata.
#'
#' @param .data A SummarizedExperiment object
#' @param operation Character string specifying the dplyr operation (e.g., "mutate", "filter", "select")
#' @param ... Arguments passed to the specified dplyr operation
#'
#' @return A SummarizedExperiment with modified rowData
#'
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' 
#' # Mutate operation on features
#' modify_features(airway, "mutate", gene_length_kb = gene_length / 1000)
#' 
#' # Filter operation on features  
#' modify_features(airway, "filter", gene_biotype == "protein_coding")
#' 
#' # Select operation on features
#' modify_features(airway, "select", gene_name, gene_biotype)
#' 
#' # Arrange operation on features
#' modify_features(airway, "arrange", gene_name)
#' }
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @noRd
modify_features <- function(.data, operation, ...) {
  
  # Get the dplyr function
  dplyr_fn <- get(operation, envir = asNamespace("dplyr"))
  
  # Apply the operation to rowData
  modified_rowdata <- rowData(.data) |>
    tibble::as_tibble(rownames = "rowname__") |>
    dplyr_fn(...) |>
    DataFrame_rownames()
  
  # Handle operations that might change the number of features
  if (operation %in% c("filter", "slice", "sample_n", "sample_frac", "distinct")) {
    # These operations can change the number of rows
    # We need to subset the entire SE object accordingly
    
    # Get the row indices that remain after the operation
    original_rowdata <- rowData(.data) |> tibble::as_tibble(rownames = "rowname__")
    original_rowdata$.original_index <- seq_len(nrow(original_rowdata))
    
    filtered_with_index <- original_rowdata |>
      dplyr_fn(...) 
    
    # Subset the entire SE object (works for zero or more features)
    remaining_indices <- filtered_with_index$.original_index
    result_se <- .data[remaining_indices, ]
    
    # Update rowData with the modified version (without the index column)
    rowData(result_se) <- filtered_with_index |>
      dplyr::select(-.original_index) |>
      DataFrame_rownames()
    
    return(result_se)
    
  } else {
    # Operations that don't change the number of features
    rowData(.data) <- modified_rowdata
    return(.data)
  }
}

#' Apply any dplyr operation to sample metadata (colData)
#'
#' This is a general function that allows applying any dplyr operation
#' directly to the colData of a SummarizedExperiment object, providing
#' better performance when operations only target sample metadata.
#'
#' @param .data A SummarizedExperiment object
#' @param operation Character string specifying the dplyr operation (e.g., "mutate", "filter", "select")
#' @param ... Arguments passed to the specified dplyr operation
#'
#' @return A SummarizedExperiment with modified colData
#'
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' 
#' # Mutate operation on samples
#' modify_samples(airway, "mutate", new_dex = paste0("treatment_", dex))
#' 
#' # Filter operation on samples  
#' modify_samples(airway, "filter", dex == "trt")
#' 
#' # Select operation on samples
#' modify_samples(airway, "select", dex, cell)
#' 
#' # Arrange operation on samples
#' modify_samples(airway, "arrange", dex, cell)
#' }
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @noRd
modify_samples <- function(.data, operation, ...) {
  
  # Get the dplyr function
  dplyr_fn <- get(operation, envir = asNamespace("dplyr"))
  
  # Apply the operation to colData
  modified_coldata <- colData(.data) |>
    tibble::as_tibble(rownames = "rowname__") |>
    dplyr_fn(...) |>
    DataFrame_rownames()
  
  # Handle operations that might change the number of samples
  if (operation %in% c("filter", "slice", "sample_n", "sample_frac", "distinct")) {
    # These operations can change the number of rows
    # We need to subset the entire SE object accordingly
    
    # Get the row indices that remain after the operation
    original_coldata <- colData(.data) |> tibble::as_tibble(rownames = "rowname__")
    original_coldata$.original_index <- seq_len(nrow(original_coldata))
    
    filtered_with_index <- original_coldata |>
      dplyr_fn(...) 
    
    # Subset the entire SE object (works for zero or more samples)
    remaining_indices <- filtered_with_index$.original_index
    result_se <- .data[, remaining_indices]
    
    # Update colData with the modified version (without the index column)
    colData(result_se) <- filtered_with_index |>
      dplyr::select(-.original_index) |>
      DataFrame_rownames() 

    
    return(result_se)
    
  } else {
    # Operations that don't change the number of samples
    colData(.data) <- modified_coldata
    return(.data)
  }
}
