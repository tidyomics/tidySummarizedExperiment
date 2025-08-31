#' Analyze the scope of a mutate query to determine target data types
#'
#' This function analyzes a mutate query to determine whether it targets:
#' - Only colData (sample metadata)
#' - Only rowData (feature metadata) 
#' - Only assays (expression/count data)
#' - Mixed operations (combination of the above)
#'
#' @keywords internal
#'
#' @param se A SummarizedExperiment object
#' @param ... Mutate expressions (same as dplyr::mutate)
#'
#' @return A list containing:
#'   \item{scope}{Character: "coldata_only", "rowdata_only", "assay_only", "mixed", or "unknown"}
#'   \item{targets_coldata}{Logical: TRUE if query affects colData}
#'   \item{targets_rowdata}{Logical: TRUE if query affects rowData}
#'   \item{targets_assays}{Logical: TRUE if query affects assays}
#'   \item{new_coldata_cols}{Character vector of new colData columns}
#'   \item{new_rowdata_cols}{Character vector of new rowData columns}
#'   \item{new_assay_cols}{Character vector of new assay columns}
#'   \item{analysis_method}{Character: "pre_mutate", "post_mutate", or "subset_analysis"}
#'   \item{confidence}{Character: "high", "medium", or "low"}
#'   \item{query_complexity}{Character: "simple" or "complex"}
#'   \item{column_names}{Character vector of column names being created/modified}
#'
#' @examples
#' \dontrun{
#' # Simple query analysis
#' scope <- analyze_query_scope_mutate(se, new_col = condition)
#' 
#' # Complex query analysis  
#' scope <- analyze_query_scope_mutate(se, 
#'   log_counts = log2(counts + 1),
#'   batch_info = paste0("batch_", sample_id),
#'   gene_category = ifelse(counts > 100, "high", "low"))
#' }
#'
#' @noRd
analyze_query_scope_mutate <- function(se, ...) {
  
  # Capture the mutate expressions
  dots <- rlang::enquos(...)
  
  # Extract column names being created/modified
  .cols <- names(dots)
  
  # Handle unnamed expressions (e.g., from group_split operations)
  if (is.null(.cols) || any(.cols == "")) {
    # For unnamed expressions, we can't do scope analysis
    # Return unknown scope to fall back to general tibble conversion
    return(list(
      scope = "unknown",
      query_complexity = "unknown",
      column_names = character(0),
      analysis_method = "fallback_unnamed"
    ))
  }
  
  # Analyze query complexity
  complexity_analysis <- analyze_query_complexity(se, dots)
  
  # Try dependency analysis first - this can detect mixed operations
  dependency_result <- analyze_expression_dependencies(se, dots)
  
  # If dependency analysis detects mixed scope, use it directly
  if (dependency_result$overall_scope == "mixed") {
    return(list(
      scope = "mixed",
      targets_coldata = any(sapply(dependency_result$expression_dependencies, function(x) x$uses_coldata)),
      targets_rowdata = any(sapply(dependency_result$expression_dependencies, function(x) x$uses_rowdata)),
      targets_assays = any(sapply(dependency_result$expression_dependencies, function(x) x$uses_assays)),
      new_coldata_cols = character(0),  # Cannot determine without execution
      new_rowdata_cols = character(0),  # Cannot determine without execution
      new_assay_cols = character(0),    # Cannot determine without execution
      analysis_method = "dependency_analysis",
      confidence = "high",
      query_complexity = if (complexity_analysis$is_simple) "simple" else "complex",
      column_names = .cols,
      expression_dependencies = dependency_result$expression_dependencies
    ))
  }
  
  # If dependency analysis gives a clear single scope, use it
  if (dependency_result$overall_scope %in% c("coldata_only", "rowdata_only", "assay_only")) {
    return(list(
      scope = dependency_result$overall_scope,
      targets_coldata = dependency_result$overall_scope == "coldata_only",
      targets_rowdata = dependency_result$overall_scope == "rowdata_only",
      targets_assays = dependency_result$overall_scope == "assay_only",
      new_coldata_cols = character(0),  # Cannot determine without execution
      new_rowdata_cols = character(0),  # Cannot determine without execution
      new_assay_cols = character(0),    # Cannot determine without execution
      analysis_method = "dependency_analysis",
      confidence = "high",
      query_complexity = if (complexity_analysis$is_simple) "simple" else "complex",
      column_names = .cols,
      expression_dependencies = dependency_result$expression_dependencies
    ))
  }
  
  # Fall back to existing analysis methods if dependency analysis is inconclusive
  # Choose analysis strategy based on complexity
  if (complexity_analysis$is_simple) {
    # Simple query - try pre-mutate analysis first
    result <- analyze_query_scope_pre_mutate(se, .cols)
    result$query_complexity <- "simple"
    result$column_names <- .cols
    
    # If pre-mutate gives unknown result for simple query, try subset analysis
    if (result$scope == "unknown" && result$confidence == "low") {
      result <- analyze_query_scope_subset(se, dots)
      result$query_complexity <- "simple"
      result$column_names <- .cols
    }
  } else {
    # Complex query - run on subset for speed
    result <- analyze_query_scope_subset(se, dots)
    result$query_complexity <- "complex"
    result$column_names <- .cols
  }
  
  return(result)
}

#' @name mutate
#' @rdname mutate
#' @inherit dplyr::mutate
#' @family single table verbs
#'
#' @examples
#' data(pasilla)
#' pasilla |> mutate(logcounts=log2(counts))
#'
#' @importFrom rlang enquos
#' @importFrom dplyr mutate
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @importFrom purrr map
#' @export
mutate.SummarizedExperiment <- function(.data, ...) {

       # Check if query is composed (multiple expressions)
    if (is_composed("mutate", ...)) return(decompose_tidy_operation("mutate", ...)(.data))

        # Check for scope and dispatch elegantly
        scope_report <- analyze_query_scope_mutate(.data, ...)
        scope <- scope_report$scope

        if(scope == "coldata_only") modify_samples(.data, "mutate", ...)
        else if(scope == "rowdata_only") modify_features(.data, "mutate", ...)
        else if(scope == "assay_only") mutate_assay(.data, ...)
        else if(scope == "mixed") modify_se_plyxp(.data, "mutate", scope_report, ...)
        else mutate_via_tibble(.data, ...)


}

# Simple helper function for assay-only mutations
mutate_assay <- function(.data, ...) {
    # Get the mutation expressions
    dots <- enquos(...)
    
    # For each new assay column, add it directly to assays
    for (i in seq_along(dots)) {
        new_assay_name <- names(dots)[i]
        expr <- dots[[i]]
        
        # Evaluate the expression to create the new assay
        # This assumes the expression references existing assay names
        new_assay_data <- rlang::eval_tidy(expr, data = as.list(assays(.data)))
        
        # Add the new assay
        assays(.data)[[new_assay_name]] <- new_assay_data
    }
    
    return(.data)
}

# Internal helper: perform mutate via tibble conversion and map back to SE
mutate_via_tibble <- function(.data, ...) {
    # Check that we are not modifying a key column
    cols <- enquos(...) |> names()
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    special_columns <- get_special_columns(
        # Decrease the size of the dataset
        .data[1:min(100, nrow(.data)), 1:min(20, ncol(.data))]
    ) |> c(get_needed_columns(.data))
    
    tst <-
        intersect(
            cols,
            special_columns
        ) |> 
        length() |>
        gt(0)
    
    if (tst) {
        columns <-
            special_columns |>
                paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }
    
    # If Ranges column not in query perform fast as_tibble
    skip_GRanges <-
        get_GRanges_colnames() %in% 
        cols |>
        not()
    
    return(.data |>
        as_tibble(skip_GRanges=skip_GRanges) |>
        dplyr::mutate(...) |>
        update_SE_from_tibble(.data))
}