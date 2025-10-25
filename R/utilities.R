#' Drop attribute to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the attribute
#' 
#' @noRd
#'
#' @return A tibble with an additional attribute
drop_attr <- function(var, name) {
  attr(var, name) <- NULL
  var
}

#' Drop all attribute to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @noRd
#' 
#' @return A tibble with an additional attribute
drop_all_attr <- function(var) {
  
  N <- var %>% attributes() %>% names()
  
  for (n in N) {
    var <- var %>% drop_attr(n)
  }
  
  var
}

#' Get matrix from tibble
#'
#' @keywords internal
#' 
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' as_matrix(head(dplyr::select(tidybulk::counts_mini, feature, count)), rownames=feature)
#'
#' @noRd
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames <- enquo(rownames)
  
  tbl %>%
    # Throw warning if data frame is not numerical beside the rownames column (if present)
    when(
      do_check &&
        (.) %>%
        # If rownames defined eliminate it from the data frame
        when(!quo_is_null(rownames) ~ (.)[,-1], ~ (.)) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any() ~ {
          warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
          (.)
        },
      ~ (.)
    ) %>%
    as.data.frame() %>%
    
    # Deal with rownames column if present
    when(
      !quo_is_null(rownames) ~ (.) %>%
        magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
        select(-1),
      ~ (.)
    ) %>%
    
    # Convert to matrix
    as.matrix()
}

#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal
#'
#' @param .data A tidySummarizedExperiment
#'
#' @noRd
to_tib <- function(.data) {
  colData(.data) %>%
    as.data.frame() %>%
    as_tibble(rownames = "cell")
}

# Greater than
gt <- function(a, b) {
  a > b
}

# Smaller than
st <- function(a, b) {
  a < b
}

# Negation
not <- function(is) {
  !is
}

# Raise to the power
pow <- function(a, b) {
  a^b
}

# Equals
eq <- function(a, b) {
  a == b
}

prepend <- function(x, values, before = 1) {
  n <- length(x)
  stopifnot(before > 0 && before <= n)
  if (before == 1) {
    c(values, x)
  }
  else {
    c(x[seq_len(before - 1)], values, x[before:n])
  }
}

#' Add class to abject
#'
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
#' 
#' @noRd
add_class <- function(var, name) {
  if (!name %in% class(var)) {
    class(var) <- prepend(class(var), name)
  }
  var
}

#' Remove class to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' 
#' @noRd
drop_class <- function(var, name) {
  class(var) <- class(var)[!class(var) %in% name]
  var
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySummarizedExperiment
#' @param transcripts A character
#' @param all A boolean
#'
#'
#' @return A tidySummarizedExperiment object
#'
#'
#' @noRd
get_abundance_sc_wide <- function(.data, transcripts = NULL, all = FALSE) {
  
  # Solve CRAN warnings
  . <- NULL
  
  # For SCE there is not filed for variable features
  variable_feature <- c()
  
  # Check if output would be too big without forcing
  if (
    length(variable_feature) == 0 &&
    is.null(transcripts) &&
    !all
  ) {
    stop("
                Your object does not contain variable feature labels,
                feature argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of feature names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
  }
  
  # Get variable features if existing
  if (
    length(variable_feature) > 0 &&
    is.null(transcripts) &&
    !all
  ) {
    variable_genes <- variable_feature
  } else {
    variable_genes <- NULL
  }
  
  # Just grub last assay
  assays(.data) %>%
    as.list() %>%
    tail(1) %>%
    .[[1]] %>%
    when(
      variable_genes %>% is.null() %>% `!`() ~ (.)[variable_genes, , drop = FALSE],
      transcripts %>% is.null() %>% `!`() ~ (.)[transcripts, , drop = FALSE],
      ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
    ) %>%
    as.matrix() %>%
    t() %>%
    as_tibble(rownames = "cell")
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom purrr map2
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySummarizedExperiment
#' @param transcripts A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#'
#' @return A tidySummarizedExperiment object
#'
#'
#' @noRd
get_abundance_sc_long <- function(.data, transcripts = NULL, all = FALSE,
                                  exclude_zeros = FALSE) {
  
  # Solve CRAN warnings
  . <- NULL
  
  # For SCE there is not filed for variable features
  variable_feature <- c()
  
  # Check if output would be too big without forcing
  if (
    length(variable_feature) == 0 &&
    is.null(transcripts) &&
    !all
  ) {
    stop("
                Your object does not contain variable feature labels,
                feature argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of feature names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
  }
  
  
  # Get variable features if existing
  if (
    length(variable_feature) > 0 &&
    is.null(transcripts) &&
    !all
  ) {
    variable_genes <- variable_feature
  } else {
    variable_genes <- NULL
  }
  
  assay_names <- assays(.data) %>% names()
  
  
  assays(.data) %>%
    as.list() %>%
    
    # Take active assay
    map2(
      assay_names,
      
      ~ .x %>%
        when(
          variable_genes %>% is.null() %>% `!`() ~ .x[variable_genes, , drop = FALSE],
          transcripts %>% is.null() %>% `!`() ~ .x[toupper(rownames(.x)) %in% toupper(transcripts), , drop=FALSE],
          all ~ .x,
          ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
        ) %>%
        
        # Replace 0 with NA
        when(exclude_zeros ~ (.) %>% {
          x <- (.)
          x[x == 0] <- NA
          x
        }, ~ (.)) %>%
        as.matrix() %>%
        data.frame() %>%
        as_tibble(rownames=f_(.data)$name) %>%
        tidyr::pivot_longer(
          cols = -!!f_(.data)$symbol,
          names_to = "cell",
          values_to = "abundance" %>% paste(.y, sep = "_"),
          values_drop_na = TRUE
        )
      # %>%
      # mutate_if(is.character, as.factor) %>%
    ) %>%
    Reduce(function(...) left_join(..., by = c(f_(.data)$name, "cell")), .)
}


#' @importFrom methods .hasSlot
#' @importFrom S4Vectors DataFrame
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom S4Vectors head 
#' @importFrom rlang .data
#'
#' @keywords internal
#'
#' @param .data_mutated A tibble
#' @param se A tidySummarizedExperiment
#'
#' @noRd
update_SE_from_tibble <- function(.data_mutated, se, column_belonging = NULL) {
  
  # Comply to CRAN notes 
  . <- NULL 
  
  # Extract column data
  col_data <- extract_col_data(.data_mutated, se, column_belonging)
  
  # This to avoid the mismatch between missing column names for counts 
  # and numerical row names for colData
  row_names_col = 
    col_data %>%
    rownames() %>%
    when(
      colnames(se) %>% is.null ~ as.integer(.),
      ~ (.)
    )
  
  # Extract row data
  row_data <- extract_row_data(.data_mutated, se, col_data = col_data)
  
  # This to avoid the mismatch between missing row names for counts 
  # and numerical row names for rowData
  row_names_row <- 
    row_data %>%
    rownames() %>%
    when(
      rownames(se) %>% is.null ~ as.integer(.),
      ~ (.)
    )
  
  # Subset if needed. This function is used by many dplyr utilities
  se <- se[row_names_row, row_names_col]
  
  # Update
  colData(se) <- col_data
  rowData(se) <- row_data
  
  # Extract assay column names
  colnames_assay <- extract_assay_colnames(.data_mutated, se, col_data, row_data)
  
  if (length(colnames_assay) > 0)
    assays(se, withDimnames = FALSE) <-
    assays(se, withDimnames = FALSE) %>% c(
      .data_mutated %>% 
        
        # Select assays column
        select(!!f_(se)$symbol, !!s_(se)$symbol, colnames_assay) %>% 
        
        # Pivot for generalising to many assays
        pivot_longer(cols = -c(!!s_(se)$symbol, !!f_(se)$symbol)) %>%
        nest(data___ = c(!!s_(se)$symbol, !!f_(se)$symbol, value)) %>%
        
        # Convert to matrix and to named list
        mutate(data___ = map2(
          data___, name,
          ~ {
            .x = 
              .x %>%
              spread(!!s_(se)$symbol, value) %>% 
              as_matrix(rownames = !!f_(se)$symbol)  %>% 
              suppressWarnings()
            
            # Rearrange if assays has colnames and rownames
            if(!is.null(rownames(se)) & !is.null(rownames(.x))) .x = .x[rownames(se),,drop=FALSE]
            if(!is.null(colnames(se)) & !is.null(colnames(.x))) .x = .x[,colnames(se),drop=FALSE]
            
            
            .x %>%
              
              list() %>%
              setNames(.y)
          }
        )) %>%
        
        # Create correct list
        pull(data___) %>%
        reduce(c)
    )
  
  # return
  se
}
#' @importFrom methods is
#' @keywords internal
#' @param slice_args List of arguments passed to slice
#' @param .data The data object to check for grouping
#' @param .by The .by argument from slice
#' @return Logical indicating if this is an ungrouped range slice
#' @noRd
is_range_slice_ungrouped_detected <- function(slice_args, .data, .by= NULL) {
  # First check if this is a range slice
  is_range_slice <- any(sapply(slice_args, function(x) {
    if (is.numeric(x)) {
      length(x) > 1
    } else {
      FALSE
    }
  }))
  
  # If not a range slice, return FALSE (no error needed)
  if (!is_range_slice) {
    return(FALSE)
  }
  
  # If it is a range slice, check if the data is grouped
  # For SummarizedExperiment objects, only check .by parameter
  is_grouped <- !quo_is_null(.by)
  
  # Return TRUE if ungrouped (should throw error)
  !is_grouped
}

slice_optimised <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  
  . <- NULL
  
  .by <- enquo(.by)

  # For slice operations that might select multiple rows, we need to handle the case where
  # the result cannot maintain the SummarizedExperiment structure
  # In this case, we should return a tibble with the equivalent data
  # Check if this is an ungrouped range slice
  slice_args <- list(...)
  if (is_range_slice_ungrouped_detected(slice_args, .data, .by)) {
    # For range slices on ungrouped data, throw an error with a helpful message
    stop("tidySummarizedExperiment says: slice using a range doesn't work on ungrouped data. Please use group_by() first or convert to tibble with as_tibble() before using slice with ranges.")
  }
  
  
  # For single row slices, use the original optimized approach
  # This simulated tibble only gets samples and features so we know those that have been completely omitted already
  # In order to save time for the as_tibble conversion
  simulated_slice = 
    simulate_feature_sample_from_tibble(.data, !!.by) %>% 
    dplyr::slice(..., .by = !!.by, .preserve = .preserve)
  
  # Remove .by column only if it's not a special column
  if (!quo_name(.by) %in% c(f_(.data)$name, s_(.data)$name)) {
    simulated_slice <- simulated_slice %>% select(-!!.by)
  }
 
  .data = 
    .data %>%
    
    # Subset the object for samples and features present in the simulated data
    .[rownames(.data) %in% simulated_slice[,f_(.data)$name][[1]], 
      colnames(.data) %in% simulated_slice[,s_(.data)$name][[1]]
      ] 
  .data = .data %>% 
    inner_join(simulated_slice, by = c(f_(.data)$name, s_(.data)$name)) 
  
  # If the result is already a tibble (due to mixed scope join), return it directly
  if (.data %>% is("tbl")) {
    return(.data)
  }
  
  # If order do not match with the one proposed by slice convert to tibble
  if (.data %>% is("tbl") %>% not()) {
    .data_for_match <- .data %>% select(!!f_(.data)$symbol, !!s_(.data)$symbol) %>% as_tibble() 
    
    x <- c(pull(.data_for_match, !!f_(.data)$symbol), pull(.data_for_match, !!s_(.data)$symbol))
    y <- c(pull(simulated_slice, !!f_(.data)$symbol), pull(simulated_slice, !!s_(.data)$symbol))
    
    if (identical(x, y) %>% not()) {
      left_join(simulated_slice, as_tibble(.data), by = c(f_(.data)$name, s_(.data)$name))
    } else {
      .data
    }
  }
}


#' @importFrom purrr map_chr
#' @importFrom dplyr select_if
#'
#' @keywords internal
#'
#' @param se A tidySummarizedExperiment
#'
#' @noRd
#'
get_special_columns <- function(.data) {
  colnames_special <-
    get_special_datasets(.data) %>%
    
    # In case any of those have feature of sample in column names
    map(
      ~ .x %>%
        select_if(!colnames(.) %in% get_needed_columns(.data)) %>%
        colnames()
    ) %>%
    unlist() %>%
    as.character()
  
  colnames_counts <-
    get_count_datasets(.data) %>%
    select(-!!f_(.data)$symbol, -!!s_(.data)$symbol) %>%
    colnames()
  
  colnames_special %>% c(colnames_counts)
}

#' @importFrom dplyr select
#' @importFrom tidyselect one_of
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom tibble rowid_to_column
#' 
#' @noRd
get_special_datasets <- function(se) {
  
  rr =  se %>%
    rowRanges() 
  
  rr %>%
    when( 
      
      # If no ranges
      as.data.frame(.) %>%
        nrow() %>%
        equals(0) ~ tibble(),
      
      # If it is a range list (multiple rows per feature)
      is(., "CompressedGRangesList") ~ {
        
        # If GRanges does not have row names
        if (is.null(rr@partitioning@NAMES)) {
          rr@partitioning@NAMES <- as.character(1:nrow(se))
        }
        
        tibble::as_tibble(rr) %>%
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(se) %>%
          nest(GRangesList = -group_name) %>%
          rename(!!f_(se)$symbol := group_name)
        
      },
      
      # If standard GRanges (one feature per line)
      ~ {
        
        # If GRanges does not have row names
        if (is.null(rr@ranges@NAMES)) {
          rr@ranges@NAMES <- as.character(1:nrow(se))
        }
        
        tibble::as_tibble(rr) %>% 
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(se) %>% 
          mutate(!!f_(se)$symbol := rr@ranges@NAMES) 
      }
      
    ) %>%
    list()
}

check_se_dimnames <- function(se) {

    # Stop if any column or row names are duplicated
    if (check_if_any_dimnames_duplicated(se, dim = "cols")) {
        stop("tidySummarizedExperiment says: some column names are duplicated")
    }
    if (check_if_any_dimnames_duplicated(se, dim = "rows")) {
        stop("tidySummarizedExperiment says: some row names are duplicated")
    }

    # Stop if column names of assays do not overlap, or if some assays have 
    # column names and others don't
    if (check_if_assays_are_NOT_overlapped(se, dim = "cols")) { 
        warning( 
            "tidySummarizedExperiment says: at least one of the assays in your SummarizedExperiment have column names, but they don't completely overlap between assays. It is strongly recommended to make the assays consistent, to avoid erroneous matching of samples." 
        )
    }
    # Same for row names
    if (check_if_assays_are_NOT_overlapped(se, dim = "rows")) { 
        warning( 
            "tidySummarizedExperiment says: at least one of the assays in your SummarizedExperiment have row names, but they don't completely overlap between assays. It is strongly recommended to make the assays consistent, to avoid erroneous matching of features." 
        )
    }
    
    # If the assays have dimnames but the SE does not, throw a warning and set 
    # the dimnames of the SE to those of the first assay with dimnames.
    # (At this point we know that all assays have the same dimnames (could be 
    # NULL), but they could be in different order)
    if (is.null(colnames(se)) && 
        length(assays(se)) > 0) {
        cn <- vapply(assays(se, withDimnames = FALSE), function(x) !is.null(colnames(x)), FALSE)
        if (any(cn)) {
            idx <- which(cn)[1]
            warning(
                "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names, but the SummarizedExperiment does not. Setting colnames(se) to column names of first assay with column names (assay ", idx, ")."
            )
            colnames(se) <- colnames(assays(se, withDimnames = FALSE)[[idx]])
        }
    }
    if (is.null(rownames(se)) && 
        length(assays(se)) > 0) {
        rn <- vapply(assays(se, withDimnames = FALSE), function(x) !is.null(rownames(x)), FALSE)
        if (any(rn)) {
            idx <- which(rn)[1]
            warning(
                "tidySummarizedExperiment says: the assays in your SummarizedExperiment have row names, but the SummarizedExperiment does not. Setting rownames(se) to row names of first assay with row names (assay ", idx, ")."
            )
            rownames(se) <- rownames(assays(se, withDimnames = FALSE)[[idx]])
        }
    }
    
    # If the assays as well as the SE have dimnames, but they don't overlap 
    # (they may be in different order), throw an error.
    if (!is.null(colnames(se)) &&
        length(assays(se)) > 0 && 
        !is.null(colnames(assays(se, withDimnames = FALSE)[[1]])) && 
        !all(colnames(assays(se, withDimnames = FALSE)[[1]]) %in% colnames(se))) {
        warning(
            "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names, but they don't agree with the column names of the SummarizedExperiment object itself. It is strongly recommended to make the assays consistent, to avoid erroneous matching of samples."
        )

    }
  
  if (is.null(rownames(se)) && 
      length(assays(se)) > 0) {
    rn <- vapply(assays(se, withDimnames = FALSE), function(x) !is.null(rownames(x)), FALSE)
    if (any(rn)) {
      idx <- which(rn)[1]
      warning(
        "tidySummarizedExperiment says: the assays in your SummarizedExperiment have row names, but the SummarizedExperiment does not. Setting rownames(se) to row names of first assay with row names (assay ", idx, ")."
      )
      rownames(se) <- rownames(assays(se, withDimnames = FALSE)[[idx]])
    }
  }
  
  # If the assays as well as the SE have dimnames, but they don't overlap 
  # (they may be in different order), throw an error.
  if (!is.null(colnames(se)) &&
      length(assays(se)) > 0 && 
      !is.null(colnames(assays(se, withDimnames = FALSE)[[1]])) && 
      !all(colnames(assays(se, withDimnames = FALSE)[[1]]) %in% colnames(se))) {
    warning(
      "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names, but they don't agree with the column names of the SummarizedExperiment object itself. It is strongly recommended to make the assays consistent, to avoid erroneous matching of samples."
    )
  }
  if (!is.null(rownames(se)) &&
      length(assays(se)) > 0 && 
      !is.null(rownames(assays(se, withDimnames = FALSE)[[1]])) && 
      !all(rownames(assays(se, withDimnames = FALSE)[[1]]) %in% rownames(se))) {
    warning(
      "tidySummarizedExperiment says: the assays in your SummarizedExperiment have row names, but they don't agree with the row names of the SummarizedExperiment object itself. It is strongly recommended to make the assays consistent, to avoid erroneous matching of features."
    )
  }
  
  se
}

#' @importFrom tidyr gather
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @importFrom SummarizedExperiment assays
#' @importFrom magrittr equals
#' 
#' @noRd
get_count_datasets <- function(se) { 
  # Check that dimnames are consistent
  se <- check_se_dimnames(se)
  
  # Join assays
  list_assays = 
    map2( 
    assays(se, withDimnames = FALSE) %>% as.list(),
    names(assays(se)),
    ~ {
      
      # If the counts are in a sparse matrix convert to a matrix
      # This might happen because the user loaded tidySummarizedExperiment and is 
      # print a SingleCellExperiment
      if (is(.x, "dgCMatrix") | is(.x, "DelayedArray")) {
        .x <- as.matrix(.x) 
      }
      
      # Rearrange if assays has colnames and rownames
      if (!is.null(rownames(se)) && !is.null(rownames(.x)) && 
          all(rownames(se) %in% rownames(.x))) {
        .x = .x[rownames(se), , drop = FALSE]
      }
      if (!is.null(colnames(se)) && !is.null(colnames(.x)) && 
          all(colnames(se) %in% colnames(.x))) {
        .x = .x[, colnames(se), drop = FALSE]
      }
      
      # If I don't have assay colnames and rownames add them
      if (!is.null(rownames(se)) && is.null(rownames(.x))) rownames(.x) = rownames(se) 
      if (!is.null(colnames(se)) && is.null(colnames(.x))) colnames(.x) = colnames(se) 
      
      .x = 
        .x %>%
        # matrix() %>%
        # as.data.frame() %>% 
        
        # In case I have a sparse matrix
        as.matrix() |> 
        tibble::as_tibble(rownames = f_(se)$name, .name_repair = "minimal") %>%
        
        # If the matrix does not have sample names, fix column names
        when(colnames(.x) %>% is.null() & is.null(colnames(se)) ~ setNames(., c(
          f_(se)$name,  seq_len(ncol(.x)) 
        )),
        ~ (.)
        ) 
      
      # Avoid dug if SE if completely empty, no rows, no columns
      if(.x |> select(-!!f_(se)$symbol) |> ncol() == 0) return(.x)
      
      .x |> 
        pivot_longer(names_to = s_(se)$name, values_to = .y, cols = -!!f_(se)$symbol, cols_vary = "slowest")
      
      #%>%
      #  rename(!!.y := count)
    }) %>%
    
    # Add dummy sample or feature if we have empty assay. 
    # This is needed for a correct visualisation of the tibble form
    map(~when(
      f_(se)$name %in% colnames(.x) %>% not ~ mutate(.x, !!f_(se)$symbol := as.character(NA)),
      s_(se)$name %in% colnames(.x) %>% not ~ mutate(.x, !!s_(se)$symbol := as.character(NA)),
      ~ .x
    )) 
  
  # If assays is non empty 
  if(list_assays |> length() > 0)
    list_assays |> 
    reduce(full_join, by = c(f_(se)$name, s_(se)$name))
  
  # If assays is empty 
  else {
    
    # If I don't have row column names
    if(se |> rownames() |> is.null()) rn = nrow(se) |> seq_len() |> as.character()
    else rn = rownames(se)
    if(se |> colnames() |> is.null()) cn = ncol(se) |> seq_len() |> as.character()
    else cn = colnames(se)
    
    expand.grid(  rn, cn  ) |> 
             setNames(c(f_(se)$name, s_(se)$name)) |> 
             tibble::as_tibble()
  }
   
  
}

get_needed_columns <- function(.data) {
  c(f_(.data)$name, s_(.data)$name)
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
#' 
#' @noRd
quo_names <- function(v) {
  v <- quo_name(quo_squash(v))
  gsub("^c\\(|`|\\)$", "", v) %>%
    strsplit(", ") %>%
    unlist()
}

#' @importFrom dplyr select
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
#'
#' @noRd
select_helper <- function(.data, ...) {
  loc <- tidyselect::eval_select(expr(c(...)), .data)
  
  dplyr::select(.data, loc)
}

outersect <- function(x, y) {
  sort(c(
    setdiff(x, y),
    setdiff(y, x)
  ))
}

#' @importFrom dplyr distinct_at
#' @importFrom dplyr vars
#' @importFrom purrr map
#' @importFrom magrittr equals
#' 
#' @noRd
get_subset_columns <- function(.data, .col) {
  
  # Comply with CRAN NOTES
  . <- NULL
  
  # Make col names
  .col <- enquo(.col)
  
  # x-annotation df
  n_x <- .data %>%
    distinct_at(vars(!!.col)) %>%
    nrow()
  
  # element wise columns
  .data %>%
    select(-!!.col) %>%
    colnames() %>%
    map(
      ~
        .x %>%
        when(
          .data %>%
            distinct_at(vars(!!.col, .x)) %>%
            nrow() %>%
            equals(n_x) ~ (.),
          ~NULL
        )
    ) %>%
    
    # Drop NULL
    {
      (.)[lengths((.)) != 0]
    } %>%
    unlist()
}

#' @importFrom purrr map_int
#' @importFrom purrr map
is_split_by_transcript <- function(.my_data) {
  
  se <- .my_data[[1]]
  tot_length <- .my_data %>% map(~ pull(.x, !!f_(se)$symbol) ) %>% unlist %>% unique() %>% length
  all_lengths <- .my_data %>% map_int(~ pull(.x, !!f_(se)$symbol) %>% unique() %>% length) 
  
  all_lengths %>% unique %>% length() %>% gt(1) |
    (all_lengths != tot_length) %>% any()
}

is_split_by_sample <- function(.my_data) {
  
  se <- .my_data[[1]]
  tot_length <- .my_data %>% map(~ pull(.x, !!s_(se)$symbol) ) %>% unlist %>% unique() %>% length
  all_lengths <- .my_data %>% map_int(~ pull(.x, !!s_(se)$symbol) %>% unique() %>% length) 
  
  all_lengths %>% unique %>% length() %>% gt(1) |
    (all_lengths != tot_length) %>% any()
}


get_GRanges_colnames <- function() {
  "GenomicRanges"
}


eliminate_GRanges_metadata_columns_also_present_in_Rowdata <- function(.my_data, se) {
  .my_data %>%
    select(-one_of(colnames(rowData(se)))) %>%
    
    # In case there is not metadata column
    suppressWarnings() 
}

subset_tibble_output <- function(.data, count_info, sample_info, gene_info, range_info, .subset) {
  # This function outputs a tibble after subsetting the columns
  .subset <- enquo(.subset)
  
  # Build template of the output
  output_colnames <- 
    slice(count_info, 0) %>%
    left_join(slice(sample_info, 0), by = s_(.data)$name) %>%
    left_join(slice(gene_info, 0), by = f_(.data)$name) %>%
    when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by = f_(.data)$name), ~ (.)) %>%
    select(!!.subset) %>%
    colnames()
  
  
  # Sample table
  sample_info <- 
    sample_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(s_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  range_info <- 
    range_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(f_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  gene_info <- 
    gene_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(f_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  count_info <- 
    count_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(f_(.data)$name, s_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  if (
    !is.null(count_info) & 
    (
      !is.null(sample_info) & !is.null(gene_info) | 
      
      # Make exception for weirs cases (e.g. c(sample, counts))
      (colnames(count_info) %>% outersect(c(f_(.data)$name, s_(.data)$name)) %>% 
       length() %>% gt(0))
    )
  ) {
    output_df <- 
      count_info %>%
      when(!is.null(sample_info) ~ (.) %>% left_join(sample_info, by=s_(.data)$name), ~ (.)) %>%
      when(!is.null(gene_info) ~ (.) %>% left_join(gene_info, by=f_(.data)$name), ~ (.)) %>%
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.))
  } else if (!is.null(sample_info) ) {
    output_df <- sample_info
  } else if (!is.null(gene_info)) {
    output_df <- gene_info %>%
      
      # If present join GRanges
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.))
  }
  
  output_df %>%
    
    # Cleanup
    select(one_of(output_colnames)) %>%
    suppressWarnings()
}

#' @importFrom stringr str_replace
change_reserved_column_names <- function(col_data, .data) {
  
  # Fix  NOTEs
  . = NULL
  
  col_data %>%
    
    setNames(
      colnames(.) %>% 
        sapply(function(x) if (x == f_(.data)$name) sprintf("%s.x", f_(.data)$name) else x) %>% 
        sapply(function(x) if (x == s_(.data)$name) sprintf("%s.x", s_(.data)$name) else x) %>% 
        str_replace("^coordinate$", "coordinate.x")
    ) 
  
}

choose_name_if_present <- function(x) {
  columns_query <- c()
  for (i in 1:length(x)) {
    if (is.null(names(x[i]))) columns_query[i] <- x[i]
    else columns_query[i] <- names(x[i])
  }
  
  columns_query
}


get_ellipse_colnames <- function(...) {
  (enquos(..., .ignore_empty = "all") %>% map(~ quo_name(.x)) %>% unlist)
}

get_colnames_col <- function(x) {
  colnames(colData(x)) %>% 
    c(s_(x)$name) 
}

get_rownames_col <- function(x) {
  x %>%
    when(
      .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
      TRUE ~ c()
    ) %>% 
    c(f_(x)$name) 
}

# This function is used for the change of special sample column to .sample
# Check if "sample" is included in the query and is not part of any other existing annotation
#' @importFrom stringr str_detect
#' @importFrom stringr regex
is_sample_feature_deprecated_used <- function(.data, user_columns, use_old_special_names = FALSE) {
  
  old_standard_is_used_for_sample <- 
    (
      ( any(str_detect(user_columns, regex("\\bsample\\b"))) & 
          !any(str_detect(user_columns, regex("\\W*(\\.sample)\\W*")))  ) |
        "sample" %in% user_columns 
    ) & 
    !"sample" %in% c(colnames(rowData(.data)), colnames(colData(.data)))
  
  old_standard_is_used_for_feature <- 
    (
      ( any(str_detect(user_columns, regex("\\bfeature\\b"))) & 
          !any(str_detect(user_columns, regex("\\W*(\\.feature)\\W*")))  ) |
        "feature" %in% user_columns 
    ) & 
    !"feature" %in% c(colnames(rowData(.data)), colnames(colData(.data)))
  
  old_standard_is_used <- old_standard_is_used_for_sample | old_standard_is_used_for_feature
  
  if (old_standard_is_used) {
    warning("tidySummarizedExperiment says: from version 1.3.1, the special columns including sample/feature id (colnames(se), rownames(se)) has changed to \".sample\" and \".feature\". This dataset is returned with the old-style vocabulary (feature and sample), however we suggest to update your workflow to reflect the new vocabulary (.feature, .sample)")
    
    use_old_special_names <- TRUE
  }
  
  use_old_special_names
}

data_frame_returned_message <- "tidySummarizedExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names <- "tidySummarizedExperiment says: This operation lead to duplicated feature names. A data frame is returned for independent data analysis."

# Key column names
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
ping_old_special_column_into_metadata <- function(.data) {
  
  metadata(.data)$feature__ <- get_special_column_name_symbol("feature")
  metadata(.data)$sample__ <- get_special_column_name_symbol("sample")
  
  .data
}

get_special_column_name_symbol <- function(name) {
  list(name = name, symbol = as.symbol(name))
}

# This function produce artificially the feature ans sample column, 
# to make optimisation before as_tibble is called
# for big datasets
simulate_feature_sample_from_tibble <- function(.data, ...) {

    . <- NULL
    r <- rownames(.data) %>% .[rep(1:length(.), ncol(.data) )]
    c <- colnames(.data) %>% .[rep(1:length(.), each = nrow(.data) )]
    
    # Create base tibble with feature and sample columns
    result_tibble <- tibble(!!f_(.data)$symbol := r,  !!s_(.data)$symbol := c)
    
    # BELOW IS THE MORE COMPLEX GROUP CASE
    # THAT SLOWS DOWS 10x
    # IN THE FUTURE THIS COULD BE OPTIMISED
    
    # Store original dimensions for use in rowData/colData processing
    n_features_orig <- nrow(.data)
    n_samples_orig <- ncol(.data)
    feature_names_orig <- rownames(.data)
    sample_names_orig <- colnames(.data)
    
    # Store column names to avoid issues with .data references
    feature_col_name <- f_(.data)$name
    sample_col_name <- s_(.data)$name
    feature_col_symbol <- f_(.data)$symbol
    sample_col_symbol <- s_(.data)$symbol
    
    # Check if ... parameters reference columns in rowData or colData
    # and attach only those columns to the output
    
    # Extract column names from ... parameters
    slice_colnames <- get_ellipse_colnames(...)
    
    # Filter out any non-character elements (like quosures)
    slice_colnames <- slice_colnames[is.character(slice_colnames)]
    
    # Only proceed if there are slice parameters
    if (length(slice_colnames) > 0) {
        # Get available rowData and colData column names
        rowdata_colnames <- if (.hasSlot(.data, "rowData") | .hasSlot(.data, "elementMetadata")) {
            colnames(rowData(.data))
        } else {
            character(0)
        }
        
        coldata_colnames <- colnames(colData(.data))
        
        # Check which slice columns are in rowData
        rowdata_matches <- intersect(slice_colnames, rowdata_colnames)
        
        # Check which slice columns are in colData  
        coldata_matches <- intersect(slice_colnames, coldata_colnames)
        
        # Add rowData columns if they are referenced in ...
        if (length(rowdata_matches) > 0) {
            # Create rowData tibble with repeated rows for each sample
            rowdata_df <- as.data.frame(rowData(.data)[, rowdata_matches, drop = FALSE])
            rowdata_df[[feature_col_name]] <- rownames(rowData(.data))
            
            # Create expanded tibble more efficiently
            rowdata_tibble <- tibble(
                !!feature_col_symbol := rep(rowdata_df[[feature_col_name]], each = n_samples_orig),
                !!sample_col_symbol := rep(sample_names_orig, n_features_orig)
            )
            
            # Add the rowData columns
            for (col in rowdata_matches) {
                rowdata_tibble[[col]] <- rep(rowdata_df[[col]], each = n_samples_orig)
            }
            
            # Join with result tibble
            result_tibble <- result_tibble %>%
                left_join(rowdata_tibble, by = c(feature_col_name, sample_col_name))
        }
        
        # Add colData columns if they are referenced in ...
        if (length(coldata_matches) > 0) {
            # Create colData tibble with repeated columns for each feature
            coldata_df <- as.data.frame(colData(.data)[, coldata_matches, drop = FALSE])
            coldata_df[[sample_col_name]] <- rownames(colData(.data))
            
            # Create expanded tibble more efficiently
            coldata_tibble <- tibble(
                !!feature_col_symbol := rep(feature_names_orig, each = n_samples_orig),
                !!sample_col_symbol := rep(coldata_df[[sample_col_name]], n_features_orig)
            )
            
            # Add the colData columns
            for (col in coldata_matches) {
                coldata_tibble[[col]] <- rep(coldata_df[[col]], n_features_orig)
            }
            
            # Join with result tibble
            result_tibble <- result_tibble %>%
                left_join(coldata_tibble, by = c(feature_col_name, sample_col_name))
        }
    }
    
    return(result_tibble)
}

feature__ <- get_special_column_name_symbol(".feature")
sample__ <- get_special_column_name_symbol(".sample")

#' @importFrom S4Vectors metadata
f_ <- function(x) {
  # Check if old deprecated columns are used
  if ("feature__" %in% names(metadata(x))) feature__ = metadata(x)$feature__
  return(feature__)
}

#' @importFrom S4Vectors metadata
s_ <- function(x) {
  if ("sample__" %in% names(metadata(x))) sample__ = metadata(x)$sample__
  return(sample__)
}

split_SummarizedExperiment_by_feature_to_list <- function(.data) {
  if (nrow(.data) > 1000)
    message("tidySummarizedExperiment says: grouping a SummarizedExperiment by feature takes 1 minute for ~ 10,000 features.")
  map(1:nrow(.data), ~ .data[.x,])
}

is_filer_columns_in_column_selection <- function(.data, ...) {
  # columns = enquos(columns)
  tryCatch({
    .data |>
      slice(0) |>
      dplyr::filter(...)
    TRUE
  },
  error = function(e) FALSE)
}

check_if_assays_are_NOT_consistently_ordered <- function(se) {
  
  # If I have any assay at all
  assays(se) |> length() |> gt(0) &&
    
    # If I have more than one assay with colnames
    Filter(
      Negate(is.null),
      assays(se, withDimnames = FALSE) |>  
        as.list() |> 
        map(colnames)
    ) |> 
    length() |>
    gt(0) &&
    
    # If I have lack of consistency
    se |> 
    assays(withDimnames = FALSE) |>
    as.list() |>
    purrr::map_dfr(colnames) |>
    apply(1, function(x) x |> unique() |> length()) |>
    equals(1) |>
    all() |>  
    not()
}

check_if_any_dimnames_duplicated <- function(se, dim = "cols") {
    stopifnot(dim %in% c("rows", "cols"))
    if (dim == "rows") {
        dimnames_function <- rownames
        nbr_unique_dimnames_function <- function(x) length(unique(rownames(x)))
        length_function <- nrow
    } else {
        dimnames_function <- colnames
        nbr_unique_dimnames_function <- function(x) length(unique(colnames(x)))
        length_function <- ncol
    }
    
    # Check assays
    # If I have any assay at all
    assays_check <- assays(se) |> length() |> gt(0) &&
        
        # If I have at least one assay with dimnames
        Filter(
            Negate(is.null),
            assays(se, withDimnames = FALSE) |>  
                as.list() |> 
                map(dimnames_function)
        ) |> 
        length() |>
        gt(0) &&
        
        # If any named assay have fewer unique names than expected
        assays(se, withDimnames = FALSE) |>  
        as.list() |> 
        map(dimnames_function) |>
        Filter(Negate(is.null), x = _) |>
        map(unique) |> 
        map(length) |>
        reduce(min) |> 
        equals(length_function(se)) |> 
        not()
    
    # Check SE object
    se_check <- !is.null(dimnames_function(se)) &&
        nbr_unique_dimnames_function(se) != length_function(se)
    
    # Return TRUE if either of the two checks return TRUE
    assays_check || se_check
}

check_if_assays_are_NOT_overlapped <- function(se, dim = "cols") {

    stopifnot(dim %in% c("rows", "cols"))
    if (dim == "rows") {
        dimnames_function <- rownames
        length_function <- nrow
    } else {
        dimnames_function <- colnames
        length_function <- ncol
    }
    is_identical_for_reduce <- function(x,y) if (identical(x,y)) x else FALSE
    
    # If I have any assay at all
    assays(se) |> length() |> gt(0) &&
        
        # If I have at least one assay with dimnames
        Filter(
            Negate(is.null),
            assays(se, withDimnames = FALSE) |>  
                as.list() |> 
                map(dimnames_function)
        ) |> 
        length() |>
        gt(0) &&
        
        # If I have lack of consistency
        # This will be TRUE also if some assays have dimnames and other don't
        # For each assay, sort the dimnames, then check that they are all the 
        # same. Can't check for unique length, since some names may be repeated
        # If they're not all the same, the reduce() step will return FALSE; 
        # otherwise, returns the (shared) dimnames
        assays(se, withDimnames = FALSE) |>  
        as.list() |> 
        map(dimnames_function) |> 
        map(sort) |>
        reduce(is_identical_for_reduce) |> 
        is.logical()

}

order_assays_internally_to_be_consistent <- function(se) {
  se <- check_se_dimnames(se)
  se |> 
    assays(withDimnames = FALSE) =
    map2(
      assays(se, withDimnames = FALSE) %>% as.list(),
      names(assays(se)),
      ~ {
        if (!is.null(rownames(se)) && !is.null(rownames(.x)) && 
            all(rownames(se) %in% rownames(.x))) {
          .x = .x[rownames(se), , drop = FALSE]
        }
        if (!is.null(colnames(se)) && !is.null(colnames(.x)) && 
            all(colnames(se) %in% colnames(.x))) {
          .x = .x[, colnames(se), drop = FALSE]
        }
        .x
        
      })
  
  se
}

#' @importFrom SummarizedExperiment cbind
reduce_cbind_se <- function(se_list){
  
  do.call(cbind, se_list)
}

#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom SummarizedExperiment rbind
#' @importFrom SummarizedExperiment elementMetadata
#' @importFrom SummarizedExperiment elementMetadata<-
reduce_rbind_se <- function(se_list){
  
  # rbind does not accept elementMetadata so I merge and take it off
  element_metadata = se_list %>% map(elementMetadata) |> reduce(rbind)
  
  # Drop elementMetadata
  se_list = se_list |> map(~{
    elementMetadata(.x) = NULL
    .x
  })
  
  # Bind
  se = do.call(rbind, se_list)
  rm(se_list)
  
  # Put elementMetadata back in - THE (safe) ASSUMPTION IS THAT THE ORDER DOES NOT CHANGE
  elementMetadata(se)  = element_metadata
  
  # Return
  se
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