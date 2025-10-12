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
#' 1. Get decomposition info: decompose_tidy_operation("mutate", x = a, y = b)
#' 2. Create executable function: decompose_tidy_operation("mutate", x = a, y = b)(se)
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
#' mutate_fn <- decompose_tidy_operation("mutate", new_dex = dex, new_cell = cell)
#' attr(mutate_fn, "pipeline_text")
#' # "mutate(new_dex = dex) |> mutate(new_cell = cell)"
#' 
#' # Execute the decomposed pipeline
#' result_se <- mutate_fn(airway)
#' 
#' # One-liner execution (your brilliant syntax!)
#' se_result <- decompose_tidy_operation("mutate", dex_1 = dex, dex_2 = dex)(airway)
#' 
#' # With additional arguments
#' filter_preserve <- decompose_tidy_operation("filter", 
#'                                        dex == "trt", 
#'                                        cell == "N61311", 
#'                                        .preserve = FALSE)
#' attr(filter_preserve, "pipeline_text")
#' # "filter(dex == \"trt\", .preserve = FALSE) |> filter(cell == \"N61311\", .preserve = FALSE)"
#' 
#' # Multiple additional arguments
#' mutate_with_args <- decompose_tidy_operation("mutate", 
#'                                         log_counts = log2(counts + 1),
#'                                         is_treated = dex == "trt",
#'                                         .keep = "used")
#' 
#' # Universal pattern works with any dplyr operation
#' select_fn <- decompose_tidy_operation("select", dex, cell, counts)
#' group_by_fn <- decompose_tidy_operation("group_by", dex, cell)
#' summarise_fn <- decompose_tidy_operation("summarise", 
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
decompose_tidy_operation <- function(fx_name, ...) {
  
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