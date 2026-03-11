#' Extract feature-wise (transcript-wise) information
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description subset_feature_data() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with only feature-related columns
#'
#'
#'
#' @name subset_feature_data
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#'
#' @param join_GRanges If `TRUE` (default), join range info from `rowRanges` for `RangedSummarizedExperiment`. If `FALSE`, return only rowData with `.feature` as the key column.
#'
#' @details This function extracts only feature-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` with feature-related information
#'
#'
#'
#'
#' @examples
#' ## Load airway dataset for examples
#'
#'   data('airway', package = 'airway')
#'   # Ensure a 'condition' column exists for examples expecting it
#'
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#' library(airway)
#' data(airway)
#' airway <- airway[1:100, 1:5]
#' 
#' 	subset_feature_data(airway)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @docType methods
#' @rdname subset_feature_data-methods
#' @export
#'
#'
setGeneric("subset_feature_data", function(.data, join_GRanges = TRUE)
  standardGeneric("subset_feature_data"))

#' Extract sample-wise information
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description subset_sample_data() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with only sample-related columns
#'
#'
#'
#' @name subset_sample_data
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#'
#'
#' @details This function extracts only sample-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` with sample-related information
#'
#'
#'
#'
#' @examples
#' ## Load airway dataset for examples
#'
#'   data('airway', package = 'airway')
#'   # Ensure a 'condition' column exists for examples expecting it
#'
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#' library(airway)
#' data(airway)
#' airway <- airway[1:100, 1:5]
#' 
#' 	subset_sample_data(airway)
#'
#'
#' @docType methods
#' @rdname subset_sample_data-methods
#' @export
#'
#'
setGeneric("subset_sample_data", function(.data)
  standardGeneric("subset_sample_data"))



# Set internal
.subset_sample_data = function(.data) {

  colData(.data) |>
    # If reserved column names are present add .x
    setNames(
      colnames(colData(.data)) |>
        str_replace("^sample$", "sample.x")
    ) |>
    # Convert to tibble
    tibble::as_tibble(rownames = sample__$name)
}

#' subset_sample_data
#'
#' @docType methods
#' @rdname subset_sample_data-methods
#' @export
#'
#' @return A consistent object (to the input)
setMethod("subset_sample_data",
          "SummarizedExperiment",
          .subset_sample_data)

#' subset_sample_data
#'
#' @docType methods
#' @rdname subset_sample_data-methods
#' @export
#'
#' @importFrom stringr str_replace
#' @importFrom rlang enquo quo_name
#' @importFrom dplyr select left_join
#' @importFrom SummarizedExperiment colData rowData
#'
#' @return A consistent object (to the input)
setMethod("subset_sample_data",
          "RangedSummarizedExperiment",
          .subset_sample_data)



# Set internal
.subset_feature_data = function(.data, join_GRanges = TRUE) {

  # Fix NOTEs
  . = NULL

  # Get rowData and fix column names
  row_data <- rowData(.data)
  col_names <- colnames(row_data)
  fixed_col_names <- str_replace(col_names, "^feature$", "feature.x")
  row_data_fixed <- setNames(row_data, fixed_col_names)

  gene_info <-
    row_data_fixed |>
    tibble::as_tibble(rownames = feature__$name)

  if (!join_GRanges) {
    return(gene_info)
  }

  range_info <-
    get_special_datasets(.data) |>
    reduce(left_join, by = feature__$name)

  if (nrow(range_info) > 0) {
    gene_info |> left_join(range_info, by = feature__$name)
  } else {
    gene_info
  }
}

#' subset_feature_data
#'
#' @docType methods
#' @rdname subset_feature_data-methods
#' @export
#'
#' @return A consistent object (to the input)
setMethod("subset_feature_data",
          "SummarizedExperiment",
          .subset_feature_data)

#' subset_feature_data
#'
#' @docType methods
#' @rdname subset_feature_data-methods
#' @export
#'
#' @return A consistent object (to the input)
setMethod("subset_feature_data",
          "RangedSummarizedExperiment",
          .subset_feature_data)

