
#' Analyze scope of a select operation based on selected columns
#'
#' @keywords internal
#' @param se A SummarizedExperiment object
#' @param columns_query Character vector of selected column names
#' @return List with scope ("coldata_only", "rowdata_only", "assay_only", "mixed", or "unknown") and details
#' @noRd
analyze_query_scope_select <- function(se, ...) {
  
  # Derive selected columns from a tiny preview to avoid heavy work
  preview <- {
    if (ncol(se) > 0) se[1, 1, drop = FALSE] else se
  }
  tbl <- tibble::as_tibble(preview)
  loc <- tidyselect::eval_select(rlang::expr(c(...)), tbl)
  columns_query <- colnames(dplyr::select(tbl, loc))
  
  # Identify presence per domain (including special id columns)
  col_keys <- intersect(columns_query, get_colnames_col(se))
  row_keys <- intersect(columns_query, get_rownames_col(se))
  assay_keys <- intersect(columns_query, assayNames(se))
  
  targets_coldata <- length(col_keys) > 0
  targets_rowdata <- length(row_keys) > 0
  targets_assays <- length(assay_keys) > 0
  
  n_targets <- sum(targets_coldata, targets_rowdata, targets_assays)
  scope <- if (n_targets == 0) {
    "unknown"
  } else if (n_targets == 1) {
    if (targets_coldata) "coldata_only"
    else if (targets_rowdata) "rowdata_only"
    else "assay_only"
  } else {
    "mixed"
  }
  
  list(
    scope = scope,
    targets_coldata = targets_coldata,
    targets_rowdata = targets_rowdata,
    targets_assays = targets_assays,
    selected_coldata_cols = col_keys,
    selected_rowdata_cols = row_keys,
    selected_assay_cols = assay_keys,
    selected_columns = columns_query,
    analysis_method = "select_columns",
    confidence = "high"
  )
}



#' @name select
#' @rdname select
#' @inherit dplyr::select
#'
#' @examples
#' data(pasilla)
#' pasilla |> select(.sample, .feature, counts)
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr select
#' @export
select.SummarizedExperiment <- function(.data, ...) {
   
    # colnames_col <- get_colnames_col(.data)
    # colnames_row <- get_rownames_col(.data)
    # colnames_assay = .data@assays@data |> names()
    
    . <- NULL
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
  
    # Warning if column names of assays do not overlap
    if (check_if_assays_are_NOT_consistently_ordered(.data)) {
    
        warning(
            "tidySummarizedExperiment says:",
            " the assays in your SummarizedExperiment have column names,",
            " but their order is not the same. Assays were internally",
            " reordered to be consistent with each other.",
            " To avoid unwanted behaviour it is highly reccomended",
            " to have assays with the same order of colnames and rownames."
        )
        
        # reorder assay colnames before printing
        # Rearrange if assays has colnames and rownames
        .data <- order_assays_internally_to_be_consistent(.data)
    
    }

    # Determine what domains (rowData/colData/assays) the select touches
    scope_report <- analyze_query_scope_select(.data, ...)
  
    # Early fast paths: if selection is limited to a single metadata slot, return a tibble
    # expanded to the full rectangular shape (repeat rows or columns accordingly)
    if (scope_report$scope == "coldata_only") {
        col_data_tibble <- 
            colData(.data) |> 
            as_tibble(rownames = s_(.data)$name)
        message(
            "tidySummarizedExperiment says:",
            " Key columns are missing.",
            " A data frame is returned for independent data analysis."
        )
        return(
            col_data_tibble |> 
                dplyr::select(tidyselect::eval_select(rlang::expr(c(...)), col_data_tibble)) |> 
                slice(rep(1:n(), each=nrow(!!.data)))
        )
    }
    if (scope_report$scope == "rowdata_only") {
        row_data_tibble <-
            rowData(.data) |> 
            as_tibble(rownames=f_(.data)$name)
        message("tidySummarizedExperiment says:",
            " Key columns are missing.",
            " A data frame is returned for independent data analysis.")
        return(
            row_data_tibble |> 
                dplyr::select(tidyselect::eval_select(rlang::expr(c(...)), row_data_tibble)) |> 
                slice(rep(1:n(), ncol(!!.data)  ))
        )
    }

    # Mixed or assay selections: build only the pieces required by scope
    build_row <- scope_report$targets_rowdata || scope_report$scope == "mixed"
    build_col <- scope_report$targets_coldata || scope_report$scope == "mixed"
    selected_row_cols <- scope_report$selected_rowdata_cols
    selected_col_cols <- scope_report$selected_coldata_cols
    selected_assay_cols <- scope_report$selected_assay_cols

    # Build rowData selection only when rowData is targeted
    if (build_row) {
        row_data_tibble <-
            rowData(.data) |> 
            as_tibble(rownames=f_(.data)$name)
        row_data_DF <-
            row_data_tibble |> 
            select(any_of(c(selected_row_cols)), !!f_(.data)$symbol) |>
            suppressWarnings() %>% 
            data.frame(row.names=pull(., !!f_(.data)$symbol)) |>
            select(-!!f_(.data)$symbol) |>
            DataFrame()
        # If SE does not have rownames, align
        if(rownames(.data) |> is.null()) rownames(row_data_DF)  = NULL
    } else {
        row_data_tibble <- tibble::tibble(!!f_(.data)$name := character())
        row_data_DF <- S4Vectors::DataFrame()
    }
    
    # Build colData selection only when colData is targeted
    if (build_col) {
        col_data_tibble <- 
            colData(.data) |> 
            as_tibble(rownames = s_(.data)$name)
        col_data_DF <-
            col_data_tibble |>  
            select(any_of(c(selected_col_cols)), !!s_(.data)$symbol) |>
            data.frame(row.names=pull(col_data_tibble, !!s_(.data)$symbol)) |>
            select(-!!s_(.data)$symbol) |>
            DataFrame()
        # If SE does not have colnames, align
        if(colnames(.data) |> is.null()) rownames(col_data_DF)  = NULL
    } else {
        col_data_tibble <- tibble::tibble(!!s_(.data)$name := character())
        col_data_DF <- S4Vectors::DataFrame()
    }
    
    # Subset assays list to only the selected assay columns (if any)
    count_data <-
        assays(.data)@listData %>%
            .[names(assays(.data)@listData) %in% selected_assay_cols]
  
    # Fallback: if required key columns are not all present in the result, 
    # convert to tibble and delegate selection to dplyr
    if (!all(c(get_needed_columns(.data)) %in% scope_report$selected_columns)) {
        if (ncol(.data)>100) {
            message("tidySummarizedExperiment says:",
                " You are doing a complex selection both sample-wise",
                " and feature-wise. In the latter case, for efficiency",
                " (until further development), it is better to separate",
                " your selects sample-wise OR feature-wise.")
        }
        message("tidySummarizedExperiment says:",
            " Key columns are missing.",
            " A data frame is returned for independent data analysis.")
    
        .data |>
            as_tibble(skip_GRanges=TRUE) |>
            select_helper(...) 
    } else {
        # Single or multi-slot SE return: set rowData/colData/assays that were built above
        rowData(.data) <- row_data_DF
        colData(.data) <- col_data_DF
        assays(.data) <- count_data
        
        # Attach latest select scope metadata for testing/introspection
        meta <- S4Vectors::metadata(.data)
        if (is.null(meta)) meta <- list()
        meta$latest_select_scope_report <- scope_report
        S4Vectors::metadata(.data) <- meta
        .data
    }
}