#' @name rename
#' @rdname rename
#' @inherit dplyr::rename
#' @family single table verbs
#'
#' @examples
#' data(pasilla)
#' pasilla |> rename(cond=condition)
#'
#' @importFrom tidyselect eval_select
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom dplyr rename
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
rename.SummarizedExperiment <- function(.data, ...) {

    # Check if query is composed (multiple expressions)
    if (is_composed("rename", ...)) {
        return(decompose_tidy_operation("rename", ...)(.data))
    }

    # Analyze the scope of the rename operation
    dots <- enquos(...)
    # In rename, the old names are the values (what we're renaming FROM)
    old_names <- sapply(dots, quo_name)
    
    # Determine which domain the columns belong to
    coldata_cols <- old_names[old_names %in% colnames(colData(.data))]
    rowdata_cols <- old_names[old_names %in% colnames(rowData(.data))]
    

    
    # Check for special/view-only columns
    special_columns <- get_special_columns(
        # Decrease the size of the dataset
        .data[1:min(100, nrow(.data)), 1:min(20, ncol(.data))]
    ) |> c(get_needed_columns(.data))
    
    if (any(old_names %in% special_columns)) {
        columns <- intersect(old_names, special_columns) |> paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only: ",
            columns,
            " (it is not present in the colData or rowData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }
    
    # Dispatch to appropriate domain-specific function
    if (length(coldata_cols) > 0) {
        # ColData-only rename
        return(modify_samples(.data, "rename", ...))
    } else if (length(rowdata_cols) > 0) {
        # RowData-only rename  
        return(modify_features(.data, "rename", ...))
    } else {
        # No matching columns found
        stop("tidySummarizedExperiment says:",
            " the columns you are trying to rename (",
            paste(old_names, collapse=", "),
            ") are not found in colData or rowData.")
    }
}