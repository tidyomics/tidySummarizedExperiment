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



#' Check if a dplyr query is composed (has multiple expressions)
#'
#' This function determines whether a dplyr operation contains multiple expressions
#' that could benefit from decomposition and individual analysis.
#'
#' @param operation Character string specifying the dplyr operation (e.g., "mutate", "filter")
#' @param ... The expressions passed to the dplyr function
#' @return Logical indicating whether the query is composed (TRUE) or simple (FALSE)
#'
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' 
#' # Simple query - not composed
#' is_composed("mutate", new_col = dex)  # FALSE
#' 
#' # Composed query - multiple expressions
#' is_composed("mutate", col1 = dex, col2 = cell)  # TRUE
#' 
#' # Filter with multiple conditions
#' is_composed("filter", dex == "trt", cell == "N61311")  # TRUE
#' }
#'
#' @keywords internal
#' @noRd
is_composed <- function(operation, ...) {
  
  # Capture the expressions
  dots <- rlang::enquos(...)
  
  # If no expressions, not composed
  if (length(dots) == 0) {
    return(FALSE)
  }
  
  # If only one expression, not composed
  if (length(dots) == 1) {
    return(FALSE)
  }
  
  # Multiple expressions - composed
  return(TRUE)
}

#' Function factory for substitute-based query decomposition
#'
#' Creates a decomposition function that uses substitute() to transform:
#' fx(a = x, b = y, c = z, .preserve = FALSE) into: 
#' fx(a = x, .preserve = FALSE) |> fx(b = y, .preserve = FALSE) |> fx(c = z, .preserve = FALSE)
#' 
#' Handles additional arguments (like .preserve, .by, etc.) by passing them to each decomposed step.
#' 
#' Can be used in two ways:
#' 1. Get decomposition info: substitute_decompose("mutate", x = a, y = b)
#' 2. Create executable function: substitute_decompose("mutate", x = a, y = b)(se)
#'
#' @param fx_name Character name of the function (e.g., "mutate", "filter", "select")
#' @param ... Main expressions to decompose AND additional named arguments
#' @return Function that can be applied to a SummarizedExperiment, or if no SE provided, decomposition info
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(airway)
#' data(airway)
#' 
#' # Basic mutate decomposition
#' mutate_fn <- substitute_decompose("mutate", new_dex = dex, new_cell = cell)
#' attr(mutate_fn, "pipeline_text")
#' # "mutate(new_dex = dex) |> mutate(new_cell = cell)"
#' 
#' # Execute the decomposed pipeline
#' result_se <- mutate_fn(airway)
#' 
#' # One-liner execution (your brilliant syntax!)
#' se_result <- substitute_decompose("mutate", dex_1 = dex, dex_2 = dex)(airway)
#' 
#' # With additional arguments
#' filter_preserve <- substitute_decompose("filter", 
#'                                        dex == "trt", 
#'                                        cell == "N61311", 
#'                                        .preserve = FALSE)
#' attr(filter_preserve, "pipeline_text")
#' # "filter(dex == \"trt\", .preserve = FALSE) |> filter(cell == \"N61311\", .preserve = FALSE)"
#' 
#' # Multiple additional arguments
#' mutate_with_args <- substitute_decompose("mutate", 
#'                                         log_counts = log2(counts + 1),
#'                                         is_treated = dex == "trt",
#'                                         .keep = "used")
#' 
#' # Universal pattern works with any dplyr operation
#' select_fn <- substitute_decompose("select", dex, cell, counts)
#' group_by_fn <- substitute_decompose("group_by", dex, cell)
#' summarise_fn <- substitute_decompose("summarise", 
#'                                     mean_counts = mean(counts),
#'                                     n_samples = n(),
#'                                     .groups = "drop")
#' 
#' # Access metadata
#' attr(mutate_fn, "decomposed")           # TRUE
#' attr(mutate_fn, "function_name")        # "mutate"
#' attr(mutate_fn, "individual_calls")     # Character vector of each step
#' attr(mutate_fn, "additional_args")      # List of .* arguments
#' }
#'
#' @keywords internal
#' @noRd
substitute_decompose <- function(fx_name, ...) {
  
  # Capture the expressions
  dots <- rlang::enquos(...)
  
  # Separate main expressions from additional arguments (like .preserve, .by, etc.)
  dot_names <- names(dots)
  additional_args <- list()
  main_expressions <- list()
  
  for (i in seq_along(dots)) {
    name <- dot_names[i]
    if (!is.null(name) && startsWith(name, ".")) {
      # Additional argument (starts with .)
      additional_args[[name]] <- dots[[i]]
    } else {
      # Main expression to decompose
      main_expressions <- append(main_expressions, dots[i])
    }
  }
  
  # If no names, all are main expressions
  if (is.null(dot_names)) {
    main_expressions <- dots
  }
  
  # Convert back to quosures for consistency
  main_expressions <- rlang::as_quosures(main_expressions, env = parent.frame())
  
  if (length(main_expressions) <= 1) {
    # Single expression - no decomposition needed, but still return executable function
    single_function <- function(se) {
      if (fx_name == "mutate") {
        return(se %>% mutate(!!!main_expressions, !!!additional_args))
      } else if (fx_name == "filter") {
        return(se %>% filter(!!!main_expressions, !!!additional_args))
      } else if (fx_name == "select") {
        return(se %>% select(!!!main_expressions, !!!additional_args))
      } else {
        # Generic execution
        op_call <- rlang::call2(fx_name, se, !!!main_expressions, !!!additional_args)
        return(rlang::eval_tidy(op_call, env = parent.frame()))
      }
    }
    
    attr(single_function, "decomposed") <- FALSE
    attr(single_function, "function_name") <- fx_name
    attr(single_function, "expressions") <- main_expressions
    attr(single_function, "additional_args") <- additional_args
    
    return(single_function)
  }
  
  # Create individual function calls as text
  individual_calls <- vector("character", length(main_expressions))
  individual_expressions <- vector("list", length(main_expressions))
  
  # Create additional args text for inclusion in each call
  additional_args_text <- ""
  if (length(additional_args) > 0) {
    additional_parts <- character(length(additional_args))
    for (j in seq_along(additional_args)) {
      arg_name <- names(additional_args)[j]
      arg_value <- rlang::quo_text(additional_args[[j]])
      additional_parts[j] <- paste0(arg_name, " = ", arg_value)
    }
    additional_args_text <- paste0(", ", paste(additional_parts, collapse = ", "))
  }
  
  for (i in seq_along(main_expressions)) {
    # Check if expression is named
    expr_name <- names(main_expressions)[i]
    is_named <- !is.null(expr_name) && expr_name != ""
    
    # For unnamed expressions, use a placeholder name for indexing but keep original for execution
    display_name <- if (!is_named) paste0("expr_", i) else expr_name
    
    single_expr <- main_expressions[i]
    # Only set names for named expressions
    if (is_named) {
      names(single_expr) <- expr_name
    }
    
    # Create the function call as text
    expr_text <- rlang::quo_text(single_expr[[1]])
    if (is_named) {
      # Named expression
      individual_calls[i] <- paste0(fx_name, "(", expr_name, " = ", expr_text, additional_args_text, ")")
    } else {
      # Unnamed expression
      individual_calls[i] <- paste0(fx_name, "(", expr_text, additional_args_text, ")")
    }
    
    individual_expressions[[i]] <- single_expr
    names(individual_calls)[i] <- display_name
    names(individual_expressions)[i] <- display_name
  }
  
  # Create pipeline representation as text
  pipeline_text <- paste(individual_calls, collapse = " |> ")
  
  # Create the executable function
  pipeline_function <- function(se) {
    # Execute each step sequentially
    result_se <- se
    
    for (i in seq_along(individual_expressions)) {
      step_expr <- individual_expressions[[i]]
      
      # Execute this step using the appropriate dplyr function, including additional args
      if (fx_name == "mutate") {
        result_se <- result_se %>% mutate(!!!step_expr, !!!additional_args)
      } else if (fx_name == "filter") {
        result_se <- result_se %>% filter(!!!step_expr, !!!additional_args)
      } else if (fx_name == "select") {
        result_se <- result_se %>% select(!!!step_expr, !!!additional_args)
      } else {
        # Generic execution for other operations
        op_call <- rlang::call2(fx_name, result_se, !!!step_expr, !!!additional_args)
        result_se <- rlang::eval_tidy(op_call, env = parent.frame())
      }
    }
    
    return(result_se)
  }
  
  # Add metadata to the function
  attr(pipeline_function, "decomposed") <- TRUE
  attr(pipeline_function, "function_name") <- fx_name
  attr(pipeline_function, "pipeline_text") <- pipeline_text
  attr(pipeline_function, "individual_calls") <- individual_calls
  attr(pipeline_function, "individual_expressions") <- individual_expressions
  attr(pipeline_function, "additional_args") <- additional_args
  
  return(pipeline_function)
}



#' Decompose composite dplyr query into pipeline using substitute()
#'
#' This elegant approach uses substitute() to transform:
#' mutate(a = x, b = y, c = z) 
#' into: 
#' mutate(a = x) |> mutate(b = y) |> mutate(c = z)
#' 
#' Then we can analyze each step individually with existing functions!
#'
#' @param se A SummarizedExperiment object
#' @param operation Character: "mutate", "filter", "select", etc.
#' @param ... Arguments for the dplyr operation
#' @return List with decomposed pipeline and individual analyses
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(airway)
#' data(airway)
#' 
#' # Decompose and analyze a complex mutate operation
#' result <- decompose_dplyr_query(airway, "mutate", 
#'                                new_condition = dex, 
#'                                new_type = cell, 
#'                                new_counts = counts)
#' 
#' # Access decomposition info
#' result$decomposition$pipeline_text
#' result$summary$optimization_potential
#' 
#' # See which steps target which data slots
#' result$summary$coldata_steps
#' result$summary$assay_steps
#' }
#'
#' @keywords internal
#' @noRd
decompose_dplyr_query <- function(se, operation, ...) {
  
  # Capture the expressions for later use
  dots <- rlang::enquos(...)
  
  # Use the function factory for clean decomposition (pass ... directly)
  decomposition_fn <- substitute_decompose(operation, ...)
  
  if (!attr(decomposition_fn, "decomposed")) {
    # Single expression - no decomposition needed
    return(list(
      decomposed = FALSE,
      single_analysis = if (operation == "mutate") {
        analyze_query_scope_mutate(se, ...)
      } else {
        analyze_dplyr_operation_scope(se, operation, ...)
      },
      decomposition = decomposition_fn
    ))
  }
  
  # For multiple expressions, analyze each individual operation
  individual_analyses <- list()
  
  for (i in seq_along(attr(decomposition_fn, "individual_expressions"))) {
    expr_name <- names(attr(decomposition_fn, "individual_expressions"))[i]
    single_expr <- attr(decomposition_fn, "individual_expressions")[[i]]
    
    # Analyze this individual expression
    tryCatch({
      if (operation == "mutate") {
        individual_analyses[[expr_name]] <- analyze_query_scope_mutate(se, !!!single_expr)
      } else {
        # For other operations, use the general analyzer
        individual_analyses[[expr_name]] <- analyze_dplyr_operation_scope(se, operation, !!!single_expr)
      }
    }, error = function(e) {
      individual_analyses[[expr_name]] <- list(
        scope = "error",
        error_message = as.character(e$message),
        expression_name = expr_name
      )
    })
  }
  
  # Summarize results
  scopes <- sapply(individual_analyses, function(x) x$scope)
  coldata_steps <- names(scopes)[scopes == "coldata_only"]
  rowdata_steps <- names(scopes)[scopes == "rowdata_only"]
  assay_steps <- names(scopes)[scopes == "assay_only"]
  unknown_steps <- names(scopes)[scopes %in% c("unknown", "error")]
  
  return(list(
    decomposed = TRUE,
    original_query = list(operation = operation, expressions = dots),
    decomposition = decomposition_fn,
    individual_analyses = individual_analyses,
    summary = list(
      total_steps = length(dots),
      coldata_steps = coldata_steps,
      rowdata_steps = rowdata_steps,
      assay_steps = assay_steps,
      unknown_steps = unknown_steps,
      can_optimize = length(unknown_steps) == 0,
      optimization_potential = paste0(
        "Can separate into ", length(coldata_steps), " colData, ",
        length(rowdata_steps), " rowData, and ",
        length(assay_steps), " assay operations"
      )
    )
  ))
}



#' General dplyr operation decomposer - separates composite queries into components
#'
#' This is the master function that can decompose ANY dplyr operation
#' (mutate, filter, select, etc.) into separate components based on scope.
#' This makes optimization available for ALL dplyr operations automatically.
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param operation Character: "mutate", "filter", "select", etc.
#' @param ... Arguments for the dplyr operation
#' @return List with decomposed components and execution strategy
#' @noRd
analyze_dplyr_operation_components <- function(se, operation, ...) {
  
  # This is the MASTER function that makes everything work for free!
  
  if (operation == "mutate") {
    # Use the universal decomposer which already handles mutate properly
    return(decompose_dplyr_query(se, operation, ...))
    
  } else if (operation == "filter") {
    # Implement filter decomposition
    return(analyze_filter_components(se, ...))
    
  } else if (operation == "select") {
    # Implement select decomposition  
    return(analyze_select_components(se, ...))
    
  } else if (operation %in% c("left_join", "right_join", "inner_join", "full_join")) {
    # Join operations not yet implemented - use generic approach
    return(analyze_generic_operation_components(se, operation, ...))
    
  } else {
    # Generic decomposition for other operations
    return(analyze_generic_operation_components(se, operation, ...))
  }
}

#' Generic operation component analyzer
#'
#' @keywords internal
#' @noRd
analyze_generic_operation_components <- function(se, operation, ...) {
  
  # For operations we don't have specialized decomposers yet,
  # use the general approach of running on subset and analyzing result
  
  result <- analyze_dplyr_operation_scope(se, operation, ...)
  
  # Wrap in component-style output for consistency
  return(list(
    overall_analysis = result,
    component_analyses = list(), # No individual components for generic operations
    summary = list(
      total_components = 1,
      coldata_components = if (result$scope == "coldata_only") operation else character(0),
      rowdata_components = if (result$scope == "rowdata_only") operation else character(0),
      assay_components = if (result$scope == "assay_only") operation else character(0),
      unknown_components = if (result$scope == "unknown") operation else character(0),
      can_optimize = result$scope != "unknown",
      optimization_strategy = paste("Single", operation, "operation -", result$scope)
    )
  ))
}

#' Filter operation component analyzer
#'
#' @keywords internal  
#' @noRd
analyze_filter_components <- function(se, ...) {
  
  # Capture filter expressions
  dots <- rlang::enquos(...)
  
  # Analyze what columns each filter condition uses
  # This is complex because filter expressions can reference multiple columns
  
  # For now, use the existing logic pattern from tidySE filter method
  # Check if filter uses only colData, only rowData, or mixed
  
  coldata_cols <- colnames(colData(se))
  rowdata_cols <- colnames(rowData(se))
  assay_cols <- assayNames(se)
  
  # Simple heuristic: check which column names appear in the filter expressions
  filter_text <- sapply(dots, rlang::quo_text)
  combined_filter_text <- paste(filter_text, collapse = " ")
  
  uses_coldata <- any(sapply(coldata_cols, function(col) grepl(paste0("\\b", col, "\\b"), combined_filter_text)))
  uses_rowdata <- any(sapply(rowdata_cols, function(col) grepl(paste0("\\b", col, "\\b"), combined_filter_text)))
  uses_assays <- any(sapply(assay_cols, function(col) grepl(paste0("\\b", col, "\\b"), combined_filter_text)))
  
  # Determine scope
  n_targets <- sum(uses_coldata, uses_rowdata, uses_assays)
  
  scope <- if (n_targets == 0) {
    "unknown"
  } else if (n_targets == 1) {
    if (uses_coldata) "coldata_only"
    else if (uses_rowdata) "rowdata_only" 
    else "assay_only"
  } else {
    "mixed"
  }
  
  # Create component analysis similar to mutate
  return(list(
    overall_analysis = list(
      scope = scope,
      operation = "filter",
      analysis_method = "component_filter",
      uses_coldata = uses_coldata,
      uses_rowdata = uses_rowdata,
      uses_assays = uses_assays
    ),
    component_analyses = list(), # Filter is typically one logical expression
    summary = list(
      total_components = 1,
      coldata_components = if (uses_coldata) "filter" else character(0),
      rowdata_components = if (uses_rowdata) "filter" else character(0),
      assay_components = if (uses_assays) "filter" else character(0),
      unknown_components = if (scope == "unknown") "filter" else character(0),
      can_optimize = scope %in% c("coldata_only", "rowdata_only"),
      optimization_strategy = switch(scope,
        "coldata_only" = "Use se[, filtered_samples] - sample subset only",
        "rowdata_only" = "Use se[filtered_features, ] - feature subset only", 
        "mixed" = "Requires full tibble conversion",
        "unknown" = "Cannot determine filter scope"
      )
    )
  ))
}

#' Select operation component analyzer
#'
#' @keywords internal
#' @noRd  
analyze_select_components <- function(se, ...) {
  
  # For select, we need to determine which columns are being selected
  # and classify them by scope
  
  # This is simpler than filter because select just chooses columns
  dots <- rlang::enquos(...)
  
  # Use dplyr to determine which columns would be selected
  # by running on a subset
  tryCatch({
    test_tibble <- se %>% as_tibble() %>% head(1)
    selected_cols <- test_tibble %>% dplyr::select(...) %>% colnames()
    
    # Classify selected columns
    coldata_cols <- intersect(selected_cols, colnames(colData(se)))
    rowdata_cols <- intersect(selected_cols, colnames(rowData(se)))
    assay_cols <- intersect(selected_cols, assayNames(se))
    special_cols <- intersect(selected_cols, c(".sample", ".feature"))
    
    # Determine scope
    n_types <- sum(length(coldata_cols) > 0, length(rowdata_cols) > 0, length(assay_cols) > 0)
    
    scope <- if (n_types == 0 && length(special_cols) > 0) {
      "metadata_only"  # Just .sample/.feature columns
    } else if (n_types == 1) {
      if (length(coldata_cols) > 0) "coldata_only"
      else if (length(rowdata_cols) > 0) "rowdata_only"
      else "assay_only"
    } else if (n_types > 1) {
      "mixed"
    } else {
      "unknown"
    }
    
    return(list(
      overall_analysis = list(
        scope = scope,
        operation = "select",
        analysis_method = "component_select",
        selected_columns = selected_cols
      ),
      component_analyses = list(
        coldata = list(columns = coldata_cols, scope = if(length(coldata_cols) > 0) "coldata_only" else "none"),
        rowdata = list(columns = rowdata_cols, scope = if(length(rowdata_cols) > 0) "rowdata_only" else "none"),
        assays = list(columns = assay_cols, scope = if(length(assay_cols) > 0) "assay_only" else "none")
      ),
      summary = list(
        total_components = length(selected_cols),
        coldata_components = coldata_cols,
        rowdata_components = rowdata_cols,
        assay_components = assay_cols,
        unknown_components = character(0),
        can_optimize = scope != "mixed",
        optimization_strategy = switch(scope,
          "coldata_only" = "Extract colData only - very fast",
          "rowdata_only" = "Extract rowData only - very fast",
          "assay_only" = "Extract assay data only - fast",
          "metadata_only" = "Extract identifiers only - instant",
          "mixed" = "Requires full tibble conversion",
          "unknown" = "Cannot optimize"
        )
      )
    ))
    
  }, error = function(e) {
    return(list(
      overall_analysis = list(scope = "error", error_message = e$message),
      component_analyses = list(),
      summary = list(can_optimize = FALSE, optimization_strategy = "Error in select analysis")
    ))
  })
}

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

#' Backward compatibility function - analyze query scope with old interface
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param .cols Character vector of column names being created/modified in mutate
#' @param .data_mutated Optional. The resulting tibble after mutate (for post-analysis)
#' @noRd
analyze_query_scope <- function(se, .cols = NULL, .data_mutated = NULL) {
  
  # Determine analysis method
  analysis_method <- if (is.null(.data_mutated)) "pre_mutate" else "post_mutate"
  
  if (analysis_method == "pre_mutate") {
    return(analyze_query_scope_pre_mutate(se, .cols))
  } else {
    return(analyze_query_scope_post_operation(se, .cols, .data_mutated, "mutate"))
  }
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



#' Analyze the scope of any dplyr operation on SummarizedExperiment
#'
#' This function can analyze the scope of various dplyr operations by running them
#' on a subset of the data and analyzing the result.
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param operation Character: "mutate", "select", "rename", "filter", etc.
#' @param ... Arguments passed to the dplyr operation
#' @noRd
analyze_dplyr_operation_scope <- function(se, operation, ...) {
  
  # Create a small subset for fast analysis
  n_rows <- min(10, nrow(se))
  n_cols <- min(5, ncol(se))
  
  if (n_rows == 0 || n_cols == 0) {
    return(list(
      scope = "unknown",
      targets_coldata = FALSE,
      targets_rowdata = FALSE,
      targets_assays = FALSE,
      new_coldata_cols = character(0),
      new_rowdata_cols = character(0),
      new_assay_cols = character(0),
      analysis_method = paste0("subset_", operation),
      operation = operation,
      confidence = "low"
    ))
  }
  
  se_subset <- se[1:n_rows, 1:n_cols]
  
  tryCatch({
    # Convert subset to tibble
    tibble_subset <- se_subset %>% as_tibble()
    
    # Apply the dplyr operation
    result_tibble <- switch(operation,
      "mutate" = tibble_subset %>% dplyr::mutate(...),
      "select" = tibble_subset %>% dplyr::select(...),
      "rename" = tibble_subset %>% dplyr::rename(...),
      "filter" = tibble_subset %>% dplyr::filter(...),
      "summarise" = tibble_subset %>% dplyr::summarise(...),
      "summarize" = tibble_subset %>% dplyr::summarize(...),
      "transmute" = tibble_subset %>% dplyr::transmute(...),
      # Add more operations as needed
      stop(paste("Unsupported operation:", operation))
    )
    
    # Determine which columns were modified/created
    original_cols <- colnames(tibble_subset)
    result_cols <- colnames(result_tibble)
    
    # For operations that modify columns
    if (operation %in% c("mutate", "transmute", "summarise", "summarize")) {
      # Columns that are new or different
      .cols <- setdiff(result_cols, original_cols)
    } else if (operation == "rename") {
      # For rename, we need to capture the mapping - this is more complex
      # For now, assume all columns could be affected
      .cols <- result_cols
    } else if (operation == "select") {
      # For select, the selected columns are the result
      .cols <- result_cols
    } else if (operation == "filter") {
      # Filter doesn't create new columns, so no column changes
      .cols <- character(0)
    } else {
      .cols <- character(0)
    }
    
    # Use post-operation analysis
    result <- analyze_query_scope_post_operation(se_subset, .cols, result_tibble, operation)
    result$analysis_method <- paste0("subset_", operation)
    result$confidence <- "medium"  # Medium confidence due to subset
    
    return(result)
    
  }, error = function(e) {
    return(list(
      scope = "unknown",
      targets_coldata = FALSE,
      targets_rowdata = FALSE,
      targets_assays = FALSE,
      new_coldata_cols = character(0),
      new_rowdata_cols = character(0),
      new_assay_cols = character(0),
      analysis_method = paste0("subset_", operation),
      operation = operation,
      confidence = "low",
      error_message = as.character(e$message)
    ))
  })
}




# Helper function to convert tibble back to DataFrame with rownames
DataFrame_rownames <- function(data, rownames_col = "rowname__") {
  if (rownames_col %in% colnames(data)) {
    # Extract rownames and remove the rownames column
    row_names <- data[[rownames_col]]
    data_without_rownames <- data[, colnames(data) != rownames_col, drop = FALSE]
    
    # Convert to DataFrame with rownames
    result <- DataFrame(data_without_rownames, row.names = row_names)
    return(result)
  } else {
    # No rownames column, just convert to DataFrame
    return(as(data, "DataFrame"))
  }
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
