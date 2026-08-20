#' @name bind_rows
#' @rdname bind_rows
#' @inherit ttservice::bind_rows
#' @param add.cell.ids Appends the corresponding values to
#' @importFrom rlang dots_values
#' @importFrom tidyselect all_of
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom SummarizedExperiment cbind
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' @importFrom S4Vectors SimpleList
#' @importFrom ttservice bind_rows
#' @importFrom tidyprint tidy_message
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

    if (is_split_by_sample(tts) && is_split_by_transcript(tts)) {
        tidy_stop("bind_rows cannot be applied to splits both by sample- and feature-wise information")
    }
    new_obj <- if (is_split_by_sample(tts)) {
        cbind(tts[[1]], tts[[2]])
    } else if (is_split_by_transcript(tts)) {
        rbind(tts[[1]], tts[[2]])
    } else {
        cbind(tts[[1]], tts[[2]])
    }
    

    # If duplicated sample names
    if (new_obj |> colnames() |> duplicated() |> which() |> length() |> gt(0)) {
        tidy_warning("you have duplicated sample names, they will be made unique.")
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
        tidy_warning("you have duplicated sample names, they will be made unique.")
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
    else {
      bound <- tts[[1]] |>
          as_tibble(skip_GRanges = TRUE) |>
          dplyr::bind_cols(tts[[2]], .id = .id)
      if (all(colnames(tts[[2]]) %in% c(
          get_subset_columns(bound, !!s_(tts[[1]])$symbol),
          get_subset_columns(bound, !!f_(tts[[1]])$symbol)
      ))) {
        update_SE_from_tibble(bound, tts[[1]], column_belonging = column_belonging)
      } else {
        tidy_warning("The new columns do not include pure sample-wise or feature-wise. A data frame is returned for independent data analysis.")
        bound
      }
    }

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
    tidy_message(data_frame_returned_message)

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
    tidy_message(data_frame_returned_message)

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
    tidy_message(data_frame_returned_message)

    data |>
        as_tibble() |>
        dplyr::rowwise()
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
slice.SummarizedExperiment <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  
  .by <- enquo(.by)
  
    slice_optimised(.data, ..., .by = !!.by, .preserve = .preserve)
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

    tidy_message(data_frame_returned_message)

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

    tidy_message(data_frame_returned_message)

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
    tidy_message(data_frame_returned_message)

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

    quo_name_name <- if (quo_is_null(name)) NULL else quo_name(name)
    
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
            tidy_warning(paste0(
                "the assays in your SummarizedExperiment have column names, ",
                "but their order is not the same. Pulling assays can return ",
                "data in a order you don't expect. To avoid unwanted behaviour ",
                "it is highly recommended to have assays with the same order ",
                "of colnames and rownames."
            ))
            
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


