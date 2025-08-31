#' Convert mixed scope expressions to plyxp form
#'
#' This function converts expressions that reference multiple data types
#' (assays, colData, rowData) into plyxp form where colData and rowData
#' references are prefixed with .cols$ and .rows$ respectively.
#'
#' @keywords internal
#' @param expression_text Character string of the expression
#' @param referenced_coldata Character vector of colData columns referenced
#' @param referenced_rowdata Character vector of rowData columns referenced
#' @param referenced_assays Character vector of assay columns referenced
#' @return Character string of the converted expression
#' @noRd
convert_to_plyxp_form <- function(expression_text, referenced_coldata, referenced_rowdata, referenced_assays) {
  
  converted_expr <- expression_text
  
  # Convert colData references to .cols$column_name
  if (length(referenced_coldata) > 0) {
    for (col in referenced_coldata) {
      # Use word boundaries to ensure we only replace whole words
      pattern <- paste0("\\b", col, "\\b")
      replacement <- paste0(".cols$", col)
      converted_expr <- gsub(pattern, replacement, converted_expr)
    }
  }
  
  # Convert rowData references to .rows$column_name  
  if (length(referenced_rowdata) > 0) {
    for (col in referenced_rowdata) {
      # Use word boundaries to ensure we only replace whole words
      pattern <- paste0("\\b", col, "\\b")
      replacement <- paste0(".rows$", col)
      converted_expr <- gsub(pattern, replacement, converted_expr)
    }
  }
  
  # Assay references remain unchanged as they're the base context
  
  return(converted_expr)
}

#' Create a plyxp-style operation function for mixed scope operations
#'
#' This function creates a function that can efficiently execute mixed scope
#' operations using plyxp-style syntax with .cols and .rows references.
#'
#' @keywords internal
#' @param operation Character string specifying the operation (e.g., "mutate")
#' @param scope_report List from analyze_query_scope_mutate
#' @param ... Original expressions from the operation
#' @return Function that takes a SummarizedExperiment and executes the operation
#' @noRd
modify_se_plyxp <- function(.data, operation, scope_report, ...) {
  
  # Convert expressions to plyxp form (e.g., data * x -> data * .cols$x)
  dots <- rlang::enquos(...)
  converted_dots <- list()
  
  for (i in seq_along(dots)) {
    expr_name <- names(dots)[i]
    if (is.null(expr_name) || expr_name == "") {
      expr_name <- paste0("expr_", i)
    }
    
    # Get dependency info for this expression
    if (expr_name %in% names(scope_report$expression_dependencies)) {
      dep_info <- scope_report$expression_dependencies[[expr_name]]
      
      
      # Convert to plyxp form
      converted_text <- convert_to_plyxp_form(
        dep_info$expression_text,
        dep_info$referenced_coldata,
        dep_info$referenced_rowdata,
        dep_info$referenced_assays
      )
      
      # Create new quosure with converted expression
      converted_dots[[i]] <- rlang::quo(!!rlang::parse_expr(converted_text))
      names(converted_dots)[i] <- expr_name
      
    } else {
      # No dependency info, keep original
      converted_dots[[i]] <- dots[[i]]
      names(converted_dots)[i] <- expr_name
    }
  }
  
  # Convert back to quosures list
  converted_dots = rlang::as_quosures(converted_dots)

  # 1) Wrap within plyxp context
  .data <- plyxp::new_plyxp(.data)

  # 2) Execute the operation call constructed via call2
  #    e.g., mutate(.data, newdata = data * .cols$x)
  op_call <- rlang::call2(rlang::sym(operation), rlang::sym(".data"), !!!converted_dots)
  .data <- rlang::eval_tidy(op_call, env = rlang::env(.data = .data))

  # 3) Convert back to SummarizedExperiment
  .data <- plyxp::se(.data)

  return(.data)
}

