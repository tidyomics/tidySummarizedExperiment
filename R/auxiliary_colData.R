#' Column-aligned auxiliary metadata
#'
#' Extension point for extra sample-aligned (column-aligned) variables that
#' should be visible to tidy verbs but are not stored in `colData()`.
#' Downstream packages can implement a method for their class; for example
#' tidySingleCellExperiment can expose reduced dimensions as view-only columns.
#'
#' The default method returns an empty [S4Vectors::DataFrame] with one row per
#' column of `x`. Auxiliary columns are treated as colData-aligned for query
#' scope (e.g. `filter()`) and as view-only when writing back to the object.
#'
#' @param x A `SummarizedExperiment` (or subclass).
#' @param ... Additional arguments passed to methods.
#'
#' @return A [S4Vectors::DataFrame] with `nrow` equal to `ncol(x)` and row
#'   names matching `colnames(x)` when those exist.
#'
#' @examples
#' data(pasilla)
#' auxiliary_colData(pasilla)
#'
#' @export
auxiliary_colData <- function(x, ...) {
    UseMethod("auxiliary_colData")
}

#' @rdname auxiliary_colData
#' @export
auxiliary_colData.default <- function(x, ...) {
    .auxiliary_colData_empty(x)
}

#' @rdname auxiliary_colData
#' @export
auxiliary_colData.SummarizedExperiment <- function(x, ...) {
    .auxiliary_colData_empty(x)
}

#' @importFrom S4Vectors DataFrame
.auxiliary_colData_empty <- function(x) {
    nms <- colnames(x)
    if (is.null(nms)) nms <- as.character(seq_len(ncol(x)))
    DataFrame(row.names = nms)
}

colnames_auxiliary_colData <- function(x) {
    colnames(auxiliary_colData(x))
}

#' Join `auxiliary_colData()` onto a sample-aligned tibble by sample id.
#' Columns already present on `sample_tbl` are skipped.
#' @noRd
join_auxiliary_colData <- function(sample_tbl, x, aux = auxiliary_colData(x)) {
    extra <- setdiff(colnames(aux), colnames(sample_tbl))
    if (!length(extra)) {
        return(sample_tbl)
    }

    key <- s_(x)$name
    aux_tbl <- tibble::as_tibble(
        as.data.frame(aux[, extra, drop = FALSE]),
        rownames = key
    )
    dplyr::left_join(sample_tbl, aux_tbl, by = key)
}

#' Bind auxiliary columns onto a sample tibble for dplyr evaluation.
#' Returns the tibble and the names that must not be written back to colData.
#' @noRd
tbl_with_joined_auxiliary_columns <- function(.data) {
    sample_tbl <- subset_sample_data(.data)
    aux <- auxiliary_colData(.data)
    extra <- setdiff(colnames(aux), colnames(sample_tbl))
    list(
        tbl = join_auxiliary_colData(sample_tbl, .data, aux = aux),
        extra = extra
    )
}
