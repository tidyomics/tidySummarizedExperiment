#' @name bind_rows
#' @rdname bind_rows
#' @inherit ttservice::bind_rows
#' @param add.cell.ids Appends the corresponding values to
#' @noRd
#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom SummarizedExperiment cbind
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' @importFrom S4Vectors SimpleList
#' @importFrom ttservice bind_rows
#' @importFrom lifecycle deprecate_warn
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
bind_rows.SummarizedExperiment <- function(..., .id=NULL, add.cell.ids=NULL) {
    lifecycle::deprecate_warn(
        when = "1.19.5",
        what = "bind_rows()",
        with = "append_samples()",
        details = "bind_rows is not a generic method in dplyr and may cause conflicts. Use append_samples from ttservice instead."
    )
    tts <- flatten_if(dots_values(...), is_spliced)

    new_obj <- 
        tts %>%
            when(
                is_split_by_sample(.) & is_split_by_transcript(.) ~ 
                    stop("tidySummarizedExperiment says:",
                        " bind_rows cannot be applied to splits both by",
                        " sample- and feature-wise information"),
                is_split_by_sample(.) ~ cbind(.[[1]], .[[2]]) ,
                is_split_by_transcript(.) ~ rbind(.[[1]], .[[2]]),
            
                # If there is not split, then bind the samples
                ~ cbind(.[[1]], .[[2]])
            )
    

    # If duplicated sample names
    if (new_obj |> colnames() |> duplicated() |> which() |> length() |> gt(0)) {
        warning("tidySummarizedExperiment says:",
            " you have duplicated sample names, they will be made unique.")
    }
    unique_colnames <- make.unique(colnames(new_obj), sep="_")

    colnames(new_obj) <- unique_colnames

    # Change also all assays colnames
    assays(new_obj) <- assays(new_obj)@listData |> map(~ {
        colnames(.x) <- unique_colnames
        .x
    }) |> 
        SimpleList()

    new_obj
}

#' @name append_samples
#' @rdname append_samples
#' @title Append samples from multiple SummarizedExperiment objects
#' 
#' @description
#' Append samples from multiple SummarizedExperiment objects by column-binding them.
#' This function is equivalent to `cbind` but provides a tidyverse-like interface.
#' 
#' @param x First SummarizedExperiment object to combine
#' @param ... Additional SummarizedExperiment objects to combine by samples
#' @param .id Object identifier (currently not used)
#' 
#' @return A combined SummarizedExperiment object
#' 
#' @examples
#' data(se)
#' append_samples(se, se)
#' 
#' @importFrom ttservice append_samples
#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom SummarizedExperiment cbind
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' @importFrom S4Vectors SimpleList
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
append_samples.SummarizedExperiment <- function(x, ..., .id = NULL) {
    # Combine all arguments into a list
    tts <- flatten_if(list(x, ...), is_spliced)
    new_obj <- do.call(cbind, tts)

    # If duplicated sample names
    if (any(duplicated(colnames(new_obj)))) {
        warning("tidySummarizedExperiment says:",
                " you have duplicated sample names, they will be made unique.")
        unique_colnames <- make.unique(colnames(new_obj), sep = "_")
        colnames(new_obj) <- unique_colnames

        # Change also all assays colnames
        assays(new_obj) <- assays(new_obj)@listData |> lapply(function(.x) {
            colnames(.x) <- unique_colnames
            .x
        }) |> SimpleList()
    }

    new_obj
}

#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom rlang dots_values
#' @importFrom ttservice bind_cols
bind_cols_internal <- function(..., .id=NULL, column_belonging=NULL) {
    tts <- tts <- flatten_if(dots_values(...), is_spliced)


    # If I have column corresponding bind directly
    # Without tranformation to tibble
    if(!is.null(column_belonging)){
      
      # For colData
      colData_additions = column_belonging[column_belonging==s_(tts[[1]] )$name] |> names()
      
      data_frame_to_attach = 
        tts[[1]] |> 
        select(!!s_(tts[[1]] )$symbol) |> 
        suppressMessages() |> 
        bind_cols(tts[[2]] |> select(all_of(colData_additions))) |> 
        distinct() 
      
      # Set row names
      data_frame_to_attach = 
        data_frame_to_attach |> 
        select(-1) |> 
        DataFrame(row.names = data_frame_to_attach |> pull(1)) 
      
      # Reorder
      data_frame_to_attach = data_frame_to_attach[match(rownames(data_frame_to_attach), colnames(tts[[1]])), , drop=FALSE]
      
      # Attach
      colData(tts[[1]]) = cbind(colData(tts[[1]]), data_frame_to_attach)
      
      # For rowData
      rowData_additions = column_belonging[column_belonging==f_(tts[[1]] )$name] |> names()
      
      data_frame_to_attach = 
        tts[[1]] |> 
        select(!!f_(tts[[1]] )$symbol) |> 
        suppressMessages() |> 
        bind_cols(tts[[2]] |> select(all_of(rowData_additions))) |> 
        distinct() 
      
      # Set row names
      data_frame_to_attach = 
        data_frame_to_attach |> 
        select(-1) |> 
        DataFrame(row.names = data_frame_to_attach |> pull(1)) 
      
      # Reorder
      data_frame_to_attach = data_frame_to_attach[match(rownames(data_frame_to_attach), rownames(tts[[1]])), , drop=FALSE]
      
      # Attach
      rowData(tts[[1]]) = cbind(rowData(tts[[1]]), data_frame_to_attach)
      
      tts[[1]]
    }
    
    # If I DON'T have column corresponding go through tibble
    else
      tts[[1]] |> 
          as_tibble(skip_GRanges = TRUE) |>
          dplyr::bind_cols(tts[[2]], .id=.id) %>%
          when(
  
              # If the column added are not sample-wise or feature-wise return tibble
              (colnames(tts[[2]]) %in% c(
                  get_subset_columns(., !!s_(tts[[1]])$symbol),
                  get_subset_columns(., !!f_(tts[[1]])$symbol)
              )
              ) |> all() ~ update_SE_from_tibble(., tts[[1]], column_belonging = column_belonging),
  
              # Return tiblle
              ~ {
                warning("tidySummarizedExperiment says:",
                        " The new columns do not include pure sample-wise",
                        " or feature-wise. A data frame is returned for",
                        " independent data analysis.")
                (.)
              }
          )

}

bind_cols_ <- function(..., .id=NULL) { bind_cols_internal(..., .id=NULL) }

#' @rdname bind_rows
#' @aliases bind_cols
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
bind_cols.SummarizedExperiment <- bind_cols_

#' @rdname bind_rows
#' @aliases bind_cols
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
bind_cols.RangedSummarizedExperiment <- bind_cols_

#' @name distinct
#' @rdname distinct
#' @inherit dplyr::distinct
#'
#' @examples
#' data(pasilla)
#' pasilla |> distinct(.sample)
#'
#' @importFrom dplyr distinct
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
distinct.SummarizedExperiment <- function(.data, ..., .keep_all=FALSE) {
  
    # message(data_frame_returned_message)

    distinct_columns <- (enquos(..., .ignore_empty="all") |>
        map(~ quo_name(.x)) |>
        unlist())
  
    # If Ranges column not in query perform fast as_tibble
    skip_GRanges <- get_GRanges_colnames() %in% 
        distinct_columns |>
        not()
  
    # Deprecation of special column names
    if (is_sample_feature_deprecated_used(.data, distinct_columns)) {
        .data= ping_old_special_column_into_metadata(.data)
    }
  
    .data |>
        select(...) |> 
        as_tibble(skip_GRanges=skip_GRanges) |>
        dplyr::distinct(..., .keep_all=.keep_all)
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
  
    # Understand what the filtering is about
    is_filter_from_samples <- .data |>
        colData() |>
        as_tibble() |> 
        is_filer_columns_in_column_selection(...)

    is_filter_from_features <- .data |>
        rowData() |>
        as_tibble() |> 
        is_filer_columns_in_column_selection(...)
  
    # Get the simpler route if filter is only on samples
    if (is_filter_from_samples & !is_filter_from_features){
        filtered_samples <- 
            colData(.data) |> 
            as_tibble(rownames=s_(.data)$name) |> 
            dplyr::filter(..., .preserve=.preserve) |> 
            pull(!!s_(.data)$symbol)
    
        return(.data[,filtered_samples])
    } else if (!is_filter_from_samples & is_filter_from_features) {
    # Get the simpler route if filter is only on features
        filtered_features <- rowData(.data) |> 
            as_tibble(rownames=f_(.data)$name) |> 
            dplyr::filter(..., .preserve=.preserve) |> 
            pull(!!f_(.data)$symbol)
    
        return(.data[filtered_features,])
    } else {
    # If filtering is based on both features and samples
    # Do filtering
    new_meta <- .data |>
        as_tibble(skip_GRanges=TRUE) |>
        dplyr::filter(..., .preserve=.preserve)
    
        # If data cannot be a SummarizedExperiment
        if (!is_rectangular(new_meta, .data)) {
            message("tidySummarizedExperiment says:",
                " The resulting data frame is not rectangular",
                " (all genes for all samples), a tibble is returned",
                " for independent data analysis.")
            return(new_meta)
        } else {
            return(.data[  
                unique(pull(new_meta,!!f_(.data)$symbol)), 
                unique(pull(new_meta,!!s_(.data)$symbol)) 
            ])
        }
    }  
}

#' @name group_by
#' @rdname group_by
#' @inherit dplyr::group_by
#'
#' @examples
#' data(pasilla)
#' pasilla  |> group_by(.sample)
#'     
#' @importFrom dplyr group_by_drop_default
#' @importFrom dplyr group_by
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
group_by.SummarizedExperiment <- function(.data, ...,
    .add=FALSE, .drop=group_by_drop_default(.data)) {
    message(data_frame_returned_message)

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
  
    .data |>
        as_tibble() |>
        dplyr::group_by(..., .add=.add, .drop=.drop)
}

#' @name summarise
#' @aliases summarize
#' @inherit dplyr::summarise
#' @family single table verbs
#' 
#' @examples
#' data(pasilla)
#' pasilla |> summarise(mean(counts))
#'
#' @importFrom dplyr summarise
#' @importFrom purrr map
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
summarise.SummarizedExperiment <- function(.data, ...) {
    message(data_frame_returned_message)

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
  
    # If Ranges column not in query perform fast as_tibble
    skip_GRanges <-
        get_GRanges_colnames() %in% 
            (enquos(..., .ignore_empty="all") |>
            map(~ quo_name(.x)) |>
            unlist()) |>
            not()

    .data |>
        as_tibble(skip_GRanges=skip_GRanges) |>
        dplyr::summarise(...)
}

#' @name summarise
#' @rdname summarise
#' @importFrom dplyr summarize
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
summarize.SummarizedExperiment <- summarise.SummarizedExperiment

#' @name mutate
#' @rdname mutate
#' @inherit dplyr::mutate
#' @family single table verbs
#'
#' @examples
#' data(pasilla)
#' pasilla |> mutate(logcounts=log2(counts))
#'
#' @importFrom rlang enquos
#' @importFrom dplyr mutate
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @importFrom purrr map
#' @export
mutate.SummarizedExperiment <- function(.data, ...) {
    # Check that we are not modifying a key column
    cols <- enquos(...) |> names()
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    special_columns <- get_special_columns(
        # Decrease the size of the dataset
        .data[1:min(100, nrow(.data)), 1:min(20, ncol(.data))]
    ) |> c(get_needed_columns(.data))
    
    tst <-
        intersect(
            cols,
            special_columns
        ) |> 
        length() |>
        gt(0)


    if (tst) {
        columns <-
            special_columns |>
                paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }

    # If Ranges column not in query perform fast as_tibble
    skip_GRanges <-
        get_GRanges_colnames() %in% 
        cols |>
        not()
    
    .data |>
        as_tibble(skip_GRanges=skip_GRanges) |>
        dplyr::mutate(...) |>
        update_SE_from_tibble(.data)
}

#' Mutate features
#'
#' Allows mutate call on features (rowData)
#' of a SummarizedExperiment
#'
#' @param .data a SummarizedExperiment
#' @param ... extra arguments passed to dplyr::mutate
#'
#' @return a SummarizedExperiment with modified rowData
#' 
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' @export
mutate_features <- function(.data, ...) {
  feature_info <- rowData(.data) |>
    tibble::as_tibble() |>
    dplyr::mutate(...) |>
    as("DataFrame")
  rowData(.data) <- feature_info
  return(.data)
}

#' Mutate samples
#'
#' Allows mutate call on samples (colData)
#' of a SummarizedExperiment
#'
#' @param .data a SummarizedExperiment
#' @param ... extra arguments passed to dplyr::mutate
#'
#' @return a SummarizedExperiment with modified colData
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @export
mutate_samples <- function(.data, ...) {
  sample_info <- colData(.data) |>
    tibble::as_tibble() |>
    dplyr::mutate(...) |>
    as("DataFrame")
  colData(.data) <- sample_info
  return(.data)
}

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

    # Cols data frame
    cols_data_frame <- 
      bind_cols(
      colData(.data) |> as_tibble() |> slice(0),
      rowData(.data) |> as_tibble() |> slice(0)
    )
  
    # Check that we are not modifying a key column
    cols <- tidyselect::eval_select(expr(c(...)), cols_data_frame)

    # Check if column is row-wise of column-wise
    old_names <- cols_data_frame[,cols] |> colnames()
    new_names <- cols |> names()
    
    # If renaming col and row data at the same time,
    # it is too complicate, so error
    if(
        old_names %in% colnames(colData(.data)) |> any() &
        old_names %in% colnames(rowData(.data)) |> any()
    )
      stop("tidySummarizedExperiment says:",
        " renaming columns from both colData and rowData at the same time",
        " is an unfeasible abstraction using dplyr.",
        " Please run two `rename` commands for",
        " sample-wise and feature-wise columns.")
    
    special_columns <- get_special_columns(
        # Decrease the size of the dataset
        .data[1:min(100, nrow(.data)), 1:min(20, ncol(.data))]
    ) |> c(get_needed_columns(.data))
    
    tst <-
        intersect(
            cols |> names(),
            special_columns
        ) |>
        length() |>
        gt(0)

    # If column in view-only columns stop
    if (tst) {
        columns <-
            special_columns |>
            paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }

    # Rename sample annotation
    if(old_names %in% colnames(colData(.data)) |> any())
        colData(.data) <- dplyr::rename(colData(.data) |>
            as.data.frame(), ...) |>
            DataFrame()

    # Rename gene annotation
    if(old_names %in% colnames(rowData(.data)) |> any())
        rowData(.data) <- dplyr::rename(rowData(.data) |>
            as.data.frame(), ...) |>
            DataFrame()

    .data
}


#' @name rowwise
#' @rdname rowwise
#' @inherit dplyr::rowwise
#'
#' @examples
#' # TODO
#'
#' @importFrom dplyr rowwise
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
rowwise.SummarizedExperiment <- function(data, ...) {
    message(data_frame_returned_message)

    data |>
        as_tibble() |>
        dplyr::rowwise()
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
#'     mutate(new_column=1:2))
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr left_join
#' @importFrom dplyr count
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
left_join.SummarizedExperiment <- function(x, y, by=NULL,
    copy=FALSE, suffix=c(".x", ".y"), ...) {
  
    join_efficient_for_SE(x, y, by=by, copy=copy,
        suffix=suffix, dplyr::left_join, ...)

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

#' @name slice
#' @rdname slice
#' @aliases slice_head slice_tail 
#'   slice_sample slice_min slice_max
#' @inherit dplyr::slice
#' @family single table verbs
#' 
#' @examples
#' data(pasilla)
#' pasilla |> slice(1)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr slice
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
slice.SummarizedExperiment <- function(.data, ..., .preserve=FALSE) {
  
    slice_optimised(.data, ..., .preserve=.preserve)
    
    # .data |>
    #     as_tibble(skip_GRanges = T) |>
    #     dplyr::slice(..., .preserve=.preserve) %>%
    #     when(
    # 
    #         # If duplicated sample-feature pair returns tibble
    #         !is_not_duplicated(., .data) | 
    #         !is_rectangular(., .data) ~ {
    #             message(duplicated_cell_names)
    #             (.)
    #         },
    # 
    #         # Otherwise return updated tidySummarizedExperiment
    #         ~ update_SE_from_tibble(., .data)
    #     )
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
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
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

    # See if join done by sample, feature or both
    columns_query <-
        .data %>%
        {
            if(ncol(.) > 0) .[1,1, drop=FALSE]
            else (.)
        } |>
        as_tibble() |> 
        select_helper(...) |> 
        colnames()
  
    row_data_tibble <-
        rowData(.data) |> 
        as_tibble(rownames=f_(.data)$name) 
  
    row_data_DF <-
        row_data_tibble |> 
        select(one_of(columns_query), !!f_(.data)$symbol) |>
        suppressWarnings() %>% 
        data.frame(row.names=pull(., !!f_(.data)$symbol)) |>
        select(-!!f_(.data)$symbol) |>
        DataFrame()
    
    # If SE does not have rownames, 
    # I have to take them our of here, otherwise count integration, 
    # which is a matrix and behaved differently from DataFrame fails
    if(rownames(.data) |> is.null()) rownames(row_data_DF)  = NULL
    
    col_data_tibble <- 
        colData(.data) |> 
        as_tibble(rownames = s_(.data)$name) 
  
    col_data_DF <-
        col_data_tibble |>  
        select(any_of(columns_query), !!s_(.data)$symbol) |>
        data.frame(row.names=pull(col_data_tibble, !!s_(.data)$symbol)) |>
        select(-!!s_(.data)$symbol) |>
        DataFrame()
  
    # If SE does not have rownames, 
    # I have to take them our of here, otherwise count integration, 
    # which is a matrix and behaved differently from DataFrame fails
    if(colnames(.data) |> is.null()) rownames(col_data_DF)  = NULL
    
    count_data <-
        assays(.data)@listData %>%
            .[names(assays(.data)@listData) %in% columns_query]
  
    # If it's just col data
    if (ncol(row_data_DF) == 0 &
        !f_(.data)$name %in% columns_query &
        length(count_data) == 0 &
        (ncol(col_data_DF) > 0 | s_(.data)$name %in% columns_query)) {
        message(
            "tidySummarizedExperiment says:",
            " Key columns are missing.",
            " A data frame is returned for independent data analysis."
        )
    
        col_data_tibble |> 
            select_helper(...) |> 
            slice(rep(1:n(), each=nrow(!!.data)))
    }

    # If it's just row data
    else if (ncol(col_data_DF) == 0 &
        !s_(.data)$name %in% columns_query &
        length(count_data) == 0 &
        (ncol(row_data_DF) > 0 | f_(.data)$name %in% columns_query)){
        message("tidySummarizedExperiment says:",
            " Key columns are missing.",
            " A data frame is returned for independent data analysis.")
    
        row_data_tibble |> 
            select_helper(...) |> 
            slice(rep(1:n(), ncol(!!.data)  ))
    }
  
    else if (!all(c(get_needed_columns(.data)) %in% columns_query)) {
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
        rowData(.data) <- row_data_DF
        colData(.data) <- col_data_DF
        assays(.data) <- count_data
        .data
    }
}

#' @name sample_n
#' @rdname sample_n
#' @aliases sample_frac
#' @inherit dplyr::sample_n
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' data(pasilla)
#' pasilla |> sample_n(50)
#' pasilla |> sample_frac(0.1)
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr sample_n
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
sample_n.SummarizedExperiment <- function(tbl, size, replace=FALSE,
    weight=NULL, .env=NULL, ...) {
    lifecycle::signal_superseded("1.0.0", "sample_n()", "slice_sample()")

    message(data_frame_returned_message)

    tbl |>
        as_tibble() |>
        dplyr::sample_n(size, replace=replace, weight=weight, .env=.env, ...)
}

#' @rdname sample_n
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr sample_frac
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
sample_frac.SummarizedExperiment <- function(tbl, size=1, replace=FALSE,
    weight=NULL, .env=NULL, ...) {
    lifecycle::signal_superseded("1.0.0", "sample_frac()", "slice_sample()")

    message(data_frame_returned_message)

    tbl |>
        as_tibble() |>
        dplyr::sample_frac(size, replace=replace, weight=weight, .env=.env, ...)
}

#' @name count
#' @rdname count
#' @inherit dplyr::count
#' 
#' @examples
#' data(se)
#' se |> count(dex)
#'     
#' @importFrom dplyr count
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
count.SummarizedExperiment <- function(x, ..., wt=NULL,
    sort=FALSE, name=NULL, .drop=group_by_drop_default(x)) {
    message(data_frame_returned_message)

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
  
  
    # If Ranges column not in query perform fast as_tibble
    skip_GRanges <-
        get_GRanges_colnames() %in% 
        (enquos(..., .ignore_empty = "all") |>
            map(~ quo_name(.x)) |> unlist()) |>
        not()
  
    x |>
        as_tibble(skip_GRanges=skip_GRanges) |>
        dplyr::count(..., wt=!!enquo(wt), sort=sort, name=name, .drop=.drop)
}

#' @name pull
#' @rdname pull
#' @inherit dplyr::pull
#' 
#' @examples
#' data(pasilla)
#' pasilla |> pull(feature)
#'     
#' @importFrom ellipsis check_dots_used
#' @importFrom dplyr pull
#' @importFrom SummarizedExperiment assay
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
pull.SummarizedExperiment <- function(.data, var=-1, name=NULL, ...) {
  
    # Fix NOTEs
    . <- NULL
  
    var <- enquo(var)
    name <- enquo(name)

    quo_name_name <- name %>% when(quo_is_null(.) ~ NULL, quo_name(name))
    
    # Deprecation of special column names
    if (is_sample_feature_deprecated_used(
        .data, 
        quo_name(var)
    )) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    # If Ranges column not in query perform fast as_tibble
    skip_GRanges <- 
        get_GRanges_colnames() %in% 
        quo_name(var) |>
        not()
    
    # Subset column annotation
    if (all(c(quo_names(var), quo_name_name) %in% colnames(colData(.data)))) {
        return(colData(.data)[, quo_names(var)] %>%
            .[rep(1:length(.), each=nrow(.data) )])
    }
    
    # Subset row annotation
    if(all(c(quo_names(var), quo_name_name) %in% colnames(rowData(.data)))) {
        return(colData(.data)[,quo_names(var)] %>%
            .[rep(1:length(.), ncol(.data) )])
    }
    
    # This returns a vector column wise. With the first sample and all features, 
    # second sample and all features, etc..
    if(all(c(quo_names(var), quo_name_name) %in% names(.data@assays@data))){

        # Warning if column names of assays do not overlap
        if (check_if_assays_are_NOT_consistently_ordered(.data)) {
            warning(
                "tidySummarizedExperiment says:",
                " the assays in your SummarizedExperiment have column names, ",
                "but their order is not the same. Pulling assays can return ",
                "data in a order you don't expect. To avoid unwanted behaviour",
                " it is highly recommended to have assays with the same order",
                " of colnames and rownames" 
            )
            
            # reorder assay colnames before printing
            # Rearrange if assays has colnames and rownames
            .data <- order_assays_internally_to_be_consistent(.data)
        
        }
        return(assay(.data, quo_names(var)) |> as.matrix() |> as.vector())       
    }
    
    # Subset rowranges
    if (all(c(quo_names(var), quo_name_name) %in%
            colnames(as.data.frame(rowRanges(.data))))) {
        return( as.data.frame(rowRanges(.data))[,quo_names(var)] %>%
            .[rep(1:length(.), ncol(.data) )])        
    }
    
    # Otherwise (SHOULD NOT HAPPEN) use the long general procedure
    colData(.data) <- colData(.data)[,colnames(colData(.data)) %in%
        c(quo_names(var), quo_name_name), drop=FALSE ]
    rowData(.data) <- rowData(.data)[,colnames(rowData(.data)) %in%
        c(quo_names(var), quo_name_name), drop=FALSE ]
  
    .data |>
        as_tibble(skip_GRanges=skip_GRanges) |>
        dplyr::pull(var=!!var, name=!!name, ...)
}

#' @name group_split
#' @rdname group_split
#' @inherit dplyr::group_split
#' 
#' @examples
#' data(pasilla, package = "tidySummarizedExperiment")
#' pasilla |> group_split(condition)
#' pasilla |> group_split(counts > 0)
#' pasilla |> group_split(condition, counts > 0)
#' 
#' @importFrom ellipsis check_dots_used
#' @importFrom dplyr group_split
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 2.1.4, https://CRAN.R-project.org/package=dplyr
#' @export
group_split.SummarizedExperiment <- function(.tbl, ..., .keep = TRUE) {
  
  var_list <- enquos(...)
  data_nested <- NULL
  
  nested <- .tbl |> 
    mutate(!!!var_list) |> 
    nest(data_nested = -(substring(as.character(var_list), 2)))
  
  if(.keep) {
    grouped_data <- nested |> 
      pull(data_nested)
    
    grouping_cols <- nested |> 
      select(substring(as.character(var_list), 2))
    
    for(i in 1:length(grouped_data)) {
      grouped_data[[i]] <- grouped_data[[i]] |> 
        mutate(grouping_cols[i,])
    }
    
    grouped_data
    
  } else {
    nested |> 
      pull(data_nested)
  }
    
}


