
#' Analyze the scope of a filter query to determine target data types
#'
#' This function analyzes a filter query to determine whether it targets:
#' - Only colData (sample metadata)
#' - Only rowData (feature metadata)
#' - Only assays (expression/count data)
#' - Mixed operations (combination of the above)
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param ... Filter expressions (same as dplyr::filter)
#' @return A list containing scope information and expression dependencies
#' @noRd
analyze_query_scope_filter <- function(se, ...) {
  
  # Capture the filter expressions
  dots <- rlang::enquos(...)
  
  # If no predicates provided, scope is unknown (no-op)
  if (length(dots) == 0) {
    return(list(
      scope = "unknown",
      targets_coldata = FALSE,
      targets_rowdata = FALSE,
      targets_assays = FALSE,
      analysis_method = "no_predicates",
      confidence = "low",
      expression_dependencies = list()
    ))
  }
  
  # Use dependency analysis to infer scope
  dependency_result <- analyze_expression_dependencies(se, dots)
  
  return(list(
    scope = dependency_result$overall_scope,
    targets_coldata = any(sapply(dependency_result$expression_dependencies, function(x) x$uses_coldata)),
    targets_rowdata = any(sapply(dependency_result$expression_dependencies, function(x) x$uses_rowdata)),
    targets_assays = any(sapply(dependency_result$expression_dependencies, function(x) x$uses_assays)),
    analysis_method = dependency_result$analysis_method,
    confidence = dependency_result$confidence,
    expression_dependencies = dependency_result$expression_dependencies
  ))
}


#' @name filter
#' @rdname filter
#' @inherit dplyr::filter
#' 
#' @examples
#' data(pasilla)
#' pasilla |>  filter(.sample == "untrt1")
#'
#' # Learn more in ?dplyr_tidy_eval
#' 
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
filter.SummarizedExperiment <- function(.data, ..., .preserve=FALSE) {
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }

    # Decompose composed queries (multiple predicates)
    if (is_composed("filter", ...)) {
        return(decompose_tidy_operation("filter", ..., .preserve = .preserve)(.data))
    }

    # Analyze scope and dispatch
    scope_report <- analyze_query_scope_filter(.data, ...)
    scope <- scope_report$scope

    result <-
    if (scope == "coldata_only") modify_samples(.data, "filter", ..., .preserve = .preserve)
    else if (scope == "rowdata_only") modify_features(.data, "filter", ..., .preserve = .preserve)
    else if (scope == "mixed") modify_se_plyxp(.data, "filter", scope_report, ..., .preserve = .preserve)
    else {
        # Fallback: perform via tibble and map back if rectangular
        new_meta <- .data |>
            as_tibble(skip_GRanges=TRUE) |>
            dplyr::filter(..., .preserve=.preserve)

        if (!is_rectangular(new_meta, .data)) {
            message("tidySummarizedExperiment says:",
                " The resulting data frame is not rectangular",
                " (all genes for all samples), a tibble is returned",
                " for independent data analysis.")
            return(new_meta)
        } else {
            .data[  
                unique(pull(new_meta,!!f_(.data)$symbol)), 
                unique(pull(new_meta,!!s_(.data)$symbol)) 
            ]
        }
    }

    # Record latest filter scope into metadata for testing/introspection
    meta <- S4Vectors::metadata(result)
    if (is.null(meta)) meta <- list()
    meta$latest_filter_scope_report <- scope_report
    S4Vectors::metadata(result) <- meta

    return(result)
}