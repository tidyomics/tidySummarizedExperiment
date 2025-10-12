#' Extract column data from mutated tibble
#'
#' @keywords internal
#'
#' @param .data_mutated A tibble
#' @param se A tidySummarizedExperiment
#' @param column_belonging Named vector indicating column source
#' @param colnames_col Character vector of column data column names
#' @param special_columns Character vector of special columns
#' @param colnames_row Character vector of row data column names
#'
#' @return DataFrame for colData
#'
#' @noRd
extract_col_data <- function(.data_mutated, se, column_belonging = NULL, 
                             colnames_col = NULL, special_columns = NULL, 
                             colnames_row = NULL) {
  
  # Comply to CRAN notes 
  . <- NULL 
  
  # Get colnames_col if not provided
  if (is.null(colnames_col)) {
    colnames_col <- 
      colnames(colData(se)) %>% 
      c(s_(se)$name) %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging == s_(se)$name]))
  }
  
  # Get special_columns if not provided
  if (is.null(special_columns)) {
    special_columns <- get_special_columns(
      # Decrease the size of the dataset
      se[1:min(100, nrow(se)), min(1, ncol(se)):min(20, ncol(se))]
    ) 
  }
  
  # Get colnames_row if not provided
  if (is.null(colnames_row)) {
    colnames_row <- se %>%
      when(
        .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
        TRUE ~ c()
      ) %>% 
      c(f_(se)$name) %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging == f_(se)$name]))
  }
  
  # This is used if I have one column with one value that can be mapped to rows and columns
  new_colnames_col = c()
  
  # This condiction is because if I don't have any samples, the new column 
  # could be mapped to samples and return NA in the final SE
  if(ncol(se) > 0)
    new_colnames_col <- 
    .data_mutated %>%
    select_if(!colnames(.) %in% setdiff(colnames_col, s_(se)$name)) %>% 
    
    # Eliminate special columns that are read only. Assays
    select_if(!colnames(.) %in% special_columns) %>%
    select_if(!colnames(.) %in% colnames_row) %>%
    # Replace for subset
    select(!!s_(se)$symbol, get_subset_columns(., !!s_(se)$symbol)) %>% 
    colnames()
  
  col_data <-
    .data_mutated %>%
    
    select(c(colnames_col, new_colnames_col)) %>%
    
    
    # Filter Sample is NA from SE that have 0 samples
    filter(!is.na(!!s_(se)$symbol)) 
  
  # This works even if I have 0 samples
  duplicated_samples = col_data |> pull(!!s_(se)$symbol) %>% duplicated() 
  if(duplicated_samples |> which() |>  length() > 0)
    # Make fast distinct()
    col_data = col_data |> filter(!duplicated_samples)
  
  col_data = 
    col_data |> 
    
    # In case unitary SE subset does not work
    data.frame(row.names = pull(col_data, !!s_(se)$symbol), check.names = FALSE) %>%
    select(-!!s_(se)$symbol) %>%
    DataFrame(check.names = FALSE)
  
  return(col_data)
}

#' Extract row data from mutated tibble
#'
#' @keywords internal
#'
#' @param .data_mutated A tibble
#' @param se A tidySummarizedExperiment
#' @param special_columns Character vector of special columns
#' @param col_data DataFrame of column data (for exclusion)
#'
#' @return DataFrame for rowData
#'
#' @noRd
extract_row_data <- function(.data_mutated, se, special_columns = NULL, 
                             col_data = NULL) {
  
  # Comply to CRAN notes 
  . <- NULL 
  
  # Get special_columns if not provided
  if (is.null(special_columns)) {
    special_columns <- get_special_columns(
      # Decrease the size of the dataset
      se[1:min(100, nrow(se)), min(1, ncol(se)):min(20, ncol(se))]
    ) 
  }
  
  # Get col_data column names for exclusion
  col_data_colnames <- if (!is.null(col_data)) colnames(col_data) else c()
  
  row_data <-
    .data_mutated %>%
    
    # Eliminate special columns that are read only 
    select_if(!colnames(.) %in% special_columns) %>%
    
    #eliminate sample columns directly
    select_if(!colnames(.) %in% c(s_(se)$name, col_data_colnames)) %>%
    
    # select(one_of(colnames(rowData(se))))
    # Replace for subset
    select(!!f_(se)$symbol,  get_subset_columns(., !!f_(se)$symbol) ) %>%
    
    # Make fast distinct()
    filter(pull(., !!f_(se)$symbol) %>% duplicated() %>% not()) %>% 
    
    # In case unitary SE subset does not work because all same
    data.frame(row.names = pull(.,f_(se)$symbol), check.names = FALSE) %>%
    select(-!!f_(se)$symbol) %>%
    DataFrame(check.names = FALSE)
  
  return(row_data)
}

#' Extract assay column names from mutated tibble
#'
#' @keywords internal
#'
#' @param .data_mutated A tibble
#' @param se A tidySummarizedExperiment
#' @param col_data DataFrame of column data
#' @param row_data DataFrame of row data
#'
#' @return Character vector of assay column names
#'
#' @noRd
extract_assay_colnames <- function(.data_mutated, se, col_data = NULL, 
                                  row_data = NULL) {
  
  # Count-like data that is NOT in the assay slot already 
  colnames_assay <-
    colnames(.data_mutated) %>% 
    setdiff(c(f_(se)$name, s_(se)$name, colnames(as.data.frame(head(rowRanges(se), 1))) )) %>%
    setdiff(if (!is.null(col_data)) colnames(col_data) else c()) %>% 
    setdiff(if (!is.null(row_data)) colnames(row_data) else c()) %>%
    setdiff(assays(se) %>% names)
  
  return(colnames_assay)
}
