#' @importFrom purrr when
#' Analyze scope of join keys (sample vs feature) for joins
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param join_keys Character vector of join key column names (already resolved)
#' @return List with scope ("coldata_only", "rowdata_only", "mixed", or "unknown") and details
#' @noRd
analyze_query_scope_join <- function(se, join_keys) {
  col_keys <- intersect(join_keys, colnames(colData(se)))
  row_keys <- if (.hasSlot(se, "rowData") || .hasSlot(se, "elementMetadata")) {
    intersect(join_keys, colnames(rowData(se)))
  } else character(0)

  uses_col <- length(col_keys) > 0
  uses_row <- length(row_keys) > 0

  scope <- if (uses_col && uses_row) {
    "mixed"
  } else if (uses_col) {
    "coldata_only"
  } else if (uses_row) {
    "rowdata_only"
  } else {
    "unknown"
  }

  list(
    scope = scope,
    targets_coldata = uses_col,
    targets_rowdata = uses_row,
    join_keys_coldata = col_keys,
    join_keys_rowdata = row_keys
  )
}

#' @importFrom purrr when
join_efficient_for_SE <- function(x, y, by = NULL, copy = FALSE, 
                                  suffix = c(".x", ".y"), join_function, 
                                  force_tibble_route = FALSE,
                                  ...) {
  
  # Comply to CRAN notes 
  . <- NULL 
  
  # Deprecation of special column names
  if (is_sample_feature_deprecated_used(x, when(by, !is.null(.) ~ by, ~ colnames(y)))) {
    x <- ping_old_special_column_into_metadata(x)
  }
  
  # Get the colnames of samples and feature datasets
  colnames_col <- get_colnames_col(x)
  colnames_row <- get_rownames_col(x)
  
  # See if join done by sample, feature or both
  columns_query <- by %>% when(
    !is.null(.) ~ choose_name_if_present(.), 
    ~ colnames(y) %>% intersect(c(colnames_col, colnames_row))
  )
  
  # Analyze join scope
  scope_report <- analyze_query_scope_join(x, columns_query)
  
  if (
    # Complex/mixed scope joins go via tibble route
    scope_report$scope %in% c("mixed", "unknown") |
    
    # If join is with something else, the inefficient generic solution might work, 
    # or delegate the proper error downstream
    (!any(columns_query %in% colnames_row) & !any(columns_query %in% colnames_col)) |
    
    # Needed for internal recurrence if outcome is not valid
    force_tibble_route) {
    
    # If I have a big dataset
    if (ncol(x) > 100) message("tidySummarizedExperiment says: if you are joining a dataframe both sample-wise and feature-wise, for efficiency (until further development), it is better to separate your joins and join datasets sample-wise OR feature-wise.")
    
    out <- x %>%
      as_tibble(skip_GRanges = TRUE) %>%
      join_function(y, by = by, copy = copy, suffix = suffix, ...) %>%
      when(
        
        # If duplicated sample-feature pair returns tibble
        !is_not_duplicated(., x) | !is_rectangular(., x) ~ {
          message(duplicated_cell_names)
          message(data_frame_returned_message)
          (.)
        },
        
        # Otherwise return updated tidySummarizedExperiment
        ~ update_SE_from_tibble(., x)
      )
    
    # Attach metadata if we returned an SE
    if (methods::is(out, "SummarizedExperiment")) {
      meta <- S4Vectors::metadata(out)
      if (is.null(meta)) meta <- list()
      meta$latest_join_scope_report <- scope_report
      S4Vectors::metadata(out) <- meta
    }
    
    return(out)
  }
  
  # Join only feature-wise
  else if (scope_report$scope == "rowdata_only") {
    
    row_data_tibble <-  
      rowData(x) %>% 
      as_tibble(rownames = f_(x)$name) %>%  
      join_function(y, by = by, copy = copy, suffix = suffix, ...) 
    
    # Check if the result is not SE then take the tibble route
    if (
      is.na(pull(row_data_tibble, !!f_(x)$symbol)) %>% any | 
      duplicated(pull(row_data_tibble, !!f_(x)$symbol)) %>% any |
      pull(row_data_tibble, !!f_(x)$symbol) %>% setdiff(rownames(colData(x))) %>% length() %>% gt(0)
    ) return(join_efficient_for_SE(x, y, by = by, copy = copy, suffix = suffix, 
                                   join_function, force_tibble_route = TRUE, ...))
    
    row_data <- 
      row_data_tibble %>% 
      data.frame(row.names = pull(., !!f_(x)$symbol)) %>%
      select(-!!f_(x)$symbol) %>%
      DataFrame()
    
    # Subset in case of an inner join, or a right join
    x <- x[rownames(row_data),]  
    
    # Tranfer annotation
    rowData(x) <- row_data
    
    # Attach latest join scope metadata
    meta <- S4Vectors::metadata(x)
    if (is.null(meta)) meta <- list()
    meta$latest_join_scope_report <- scope_report
    S4Vectors::metadata(x) <- meta
    
    # Return
    x
  }
  
  # Join only sample-wise
  else if (scope_report$scope == "coldata_only") {
    
    col_data_tibble <- 
      colData(x) %>% 
      as_tibble(rownames = s_(x)$name) %>%  
      join_function(y, by = by, copy = copy, suffix = suffix, ...)
    
    # Check if the result is not SE then take the tibble route
    if (
      is.na(pull(col_data_tibble, !!s_(x)$symbol)) %>% any | 
      duplicated(pull(col_data_tibble, !!s_(x)$symbol)) %>% any |
      pull(col_data_tibble, !!s_(x)$symbol) %>% setdiff(rownames(colData(x))) %>% length() %>% gt(0)
    ) return(join_efficient_for_SE(x, y, by = by, copy = copy, suffix = suffix, 
                                   join_function, force_tibble_route = TRUE, ...))
    
    col_data <- 
      col_data_tibble %>% 
      data.frame(row.names = pull(., !!s_(x)$symbol)) %>%
      select(-!!s_(x)$symbol) %>%
      DataFrame()
    
    # Subset in case of an inner join, or a right join
    x <- x[,rownames(col_data)]  
    
    # Transfer annotation
    colData(x) <- col_data
    
    # Attach latest join scope metadata
    meta <- S4Vectors::metadata(x)
    if (is.null(meta)) meta <- list()
    meta$latest_join_scope_report <- scope_report
    S4Vectors::metadata(x) <- meta
    
    # Return
    x
  }
  
  else stop("tidySummarizedExperiment says: ERROR FOR DEVELOPERS: this option should not exist. In join utility.")
}

#' @name left_join
#' @rdname left_join
#' @inherit dplyr::left_join
#'
#' @examples
#' data(pasilla)
#' 
#' tt <- pasilla 
#' tt |> left_join(tt |>
#'     distinct(condition) |>
#'     mutate(new_column=1:2) |>
#'     slice(1))
#'
#' @importFrom dplyr left_join
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
left_join.SummarizedExperiment <- function(x, y, by = NULL,
    copy = FALSE, suffix = c(".x", ".y"), ...) {
  
    join_efficient_for_SE(x, y, by = by, copy = copy,
        suffix = suffix, dplyr::left_join, ...)

}

#' @name inner_join
#' @rdname inner_join
#' @inherit dplyr::inner_join
#'
#' @examples
#' data(pasilla)
#' 
#' tt <- pasilla 
#' tt |> inner_join(tt |>
#'     distinct(condition) |>
#'     mutate(new_column=1:2) |>
#'     slice(1))
#'
#' @importFrom dplyr inner_join
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
inner_join.SummarizedExperiment <- function(x, y, by=NULL,
    copy=FALSE, suffix=c(".x", ".y"), ...) {
  
    join_efficient_for_SE(x, y, by=by, copy=copy,
        suffix=suffix, dplyr::inner_join, ...)

}

#' @name right_join
#' @rdname right_join
#' @inherit dplyr::right_join
#'
#' @examples
#' data(pasilla)
#' 
#' tt <- pasilla
#' tt |> right_join(tt |>
#'     distinct(condition) |>
#'     mutate(new_column=1:2) |>
#'     slice(1))
#'
#' @importFrom dplyr right_join
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
right_join.SummarizedExperiment <- function(x, y, by=NULL,
    copy=FALSE, suffix=c(".x", ".y"), ...) {
  
    join_efficient_for_SE(x, y, by=by, copy=copy,
        suffix=suffix, dplyr::right_join, ...)
}

#' @name full_join
#' @rdname full_join
#' @inherit dplyr::full_join
#'
#' @examples
#' data(pasilla)
#' 
#' tt <- pasilla
#' tt |> full_join(tibble::tibble(condition="treated", dose=10))
#'
#' @importFrom dplyr full_join
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
full_join.SummarizedExperiment <- function(x, y, by=NULL,
    copy=FALSE, suffix=c(".x", ".y"), ...) {

    join_efficient_for_SE(x, y, by=by, copy=copy,
        suffix=suffix, dplyr::full_join, ...)
}