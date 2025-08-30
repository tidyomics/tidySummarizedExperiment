#' @name unnest
#' @rdname unnest
#' @inherit tidyr::unnest
#' @aliases unnest_summarized_experiment
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' tidySummarizedExperiment::pasilla |>
#'     nest(data=-condition) |>
#'     unnest(data)
#' 
#' @importFrom tidyr unnest
#' @importFrom rlang quo_name
#' @importFrom purrr imap
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
unnest.tidySummarizedExperiment_nested <- function(data, cols, ...,
    keep_empty=FALSE, ptype=NULL, names_sep=NULL,
    names_repair="check_unique", .drop, .id, .sep, .preserve) {
    
    cols <- enquo(cols)
    
    unnest_summarized_experiment(data, !!cols, ...,
        keep_empty=keep_empty, ptype=ptype,
        names_sep=names_sep, names_repair=names_repair)
   
}


#' @rdname unnest
#' @examples
#' tidySummarizedExperiment::pasilla |>
#'     nest(data=-condition) |>
#'     unnest_summarized_experiment(data)
#' 
#' @importFrom tidyr unnest
#' @importFrom purrr when
#' @importFrom rlang quo_name
#' @importFrom purrr imap
#' @importFrom purrr map_lgl
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
unnest_summarized_experiment <- function(data, cols, ...,
    keep_empty=FALSE, ptype=NULL, names_sep=NULL,
    names_repair="check_unique", .drop, .id, .sep, .preserve) {
    . <- NULL

    # Need this otherwise crashes map
    .data_ <- data
  
    cols <- enquo(cols)
  
    # If the column is not SE do normal stuff
    if (
        data %>%
        pull(!!cols) %>%
        .[[1]] %>%
        class() %>%
        as.character() %in%
        c("SummarizedExperiment",
            "RangedSummarizedExperiment") %>%
        all() %>%
        not()
    ) {
        return(
            data %>%
            drop_class("tidySummarizedExperiment_nested") %>%
            tidyr::unnest(!!cols, ..., keep_empty=keep_empty,
                ptype=ptype, names_sep=names_sep,
                names_repair=names_repair) %>%
            add_class("tidySummarizedExperiment_nested")
        )
    }


    # If both nested by transcript and sample
    if (s_(se)$name %in% colnames(data) &
        f_(se)$name %in% colnames(data) ) {
        stop("tidySummarizedExperiment says:",
            " for the moment nesting both by sample- and feature-wise",
            " information is not possible. Please ask this feature",
            " to github/stemangiola/tidySummarizedExperiment")
    }
  
    # If both nested not by transcript nor sample
    if(!s_(se)$name %in% colnames(data) &
        !f_(se)$name %in% colnames(data)) {
    
        my_se <- pull(.data_, !!cols) %>% .[[1]] 

    
        # Mark if columns belong to feature or sample
        my_unnested_tibble =
          data |> 
          mutate(!!cols := map(!!cols, ~ as_tibble(.x))) |>
          select(-any_of(c(s_(my_se)$name, f_(my_se)$name))) |> 
          unnest(!!cols)
    
        # Get which column is relative to feature or sample
        sample_columns <- my_unnested_tibble %>%
            get_subset_columns(!!s_(my_se)$symbol)
        transcript_columns <- my_unnested_tibble %>%
            get_subset_columns(!!f_(my_se)$symbol)
    
        source_column <-
            c(
                rep(s_(my_se)$name,
                    length(sample_columns)) %>%
                    setNames(sample_columns),
                rep(f_(my_se)$name,
                    length(transcript_columns)) %>%
                    setNames(transcript_columns)
            )
    
        # Drop if SE is null
        if (data |> filter(map_lgl(!!cols, is.null)) |> nrow() > 0) {
            warning("tidySummarizedcExperiment says:",
                " some SummarizedExperiment objects to",
                " unnest were <NULL>, and were elminated")
            data <- data |> filter(!map_lgl(!!cols, is.null))
        }
    
        # Do my trick to unnest
    data = 
      data |>
      
      # I have to use this because imap behave strangely
      rowid_to_column(var = "i___") |> 
      mutate(!!cols := map2(
        !!cols, i___, ~ .x %>% 
          bind_cols_internal(
            
            # Attach back the columns used for nesting
            .data_ %>%
              select(-!!cols, -any_of(c(s_(my_se)$name, f_(my_se)$name))) %>%
              slice(rep(as.integer(.y), ncol(.x) * nrow(.x))),
            
            # Column sample-wise or feature-wise
            column_belonging =
              source_column[
                .data_ %>%
                  select(-!!cols, -any_of(c(s_(my_se)$name, f_(my_se)$name))) %>%
                  colnames()
              ]
          )
      )) |> 
      
      # I have to use this because imap behave strangely
      select(-i___)
    
        # Understand if split was done feature 
        if(identical(
            data |> pull(!!cols) |> magrittr::extract2(1) |>
                colnames() |> sort(),
            data |> pull(!!cols) |> map(colnames) |> 
                reduce(intersect) |> sort()
        )) {
            return(data |> pull(!!cols) |> reduce_rbind_se())
        } 
        # Understand if split was done sample 
        else if (identical(
            data |> pull(!!cols) |> magrittr::extract2(1) |>
                rownames() |> sort(),
            data |> pull(!!cols) |> map(rownames) |>
                reduce(intersect) |> sort()
        )) {
            return(data |> pull(!!cols) |> reduce_cbind_se())
        }
        # If neither there is something wrong
        else {
            stop("tidySummarizedcExperiment says: not the sample names",
                " nor the feature names overlap through your nesting.",
                " The nesting (due to the underlying",
                " SummarizedExperiment::cbind and",
                " SummarizedExperiment::rbind requirements)",
                " needs to be rectangular.)")
        }

    }
  
    # If column is SE nd only feature
    if (f_(se)$name %in% colnames(data)) {
        se <- do.call(SummarizedExperiment::rbind, pull(data, !!cols))
        rowData(se) <- cbind(rowData(se),
            data %>% select(-!!cols, -!!f_(se)$symbol))
        return(se)
    }
  
    # If column is SE nd only sample
    if (s_(se)$name %in% colnames(data)) {
        se <- data |> pull(!!cols) |> reduce_cbind_se()
        colData(se) <- cbind(colData(se),
            data %>% select(-!!cols, -!!s_(se)$symbol))
        return(se)
    }
}

#' @name nest
#' @rdname nest
#' @inherit tidyr::nest
#' @return `tidySummarizedExperiment_nested`
#'
#' @examples
#' tidySummarizedExperiment::pasilla |>
#'     nest(data=-condition)
#'     
#' @importFrom rlang enquos
#' @importFrom rlang :=
#' @importFrom purrr when
#' @importFrom purrr pmap
#' @importFrom tidyr nest
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
nest.SummarizedExperiment <- function(.data, ..., .names_sep=NULL) {
    cols <- enquos(...)
    col_name_data <- names(cols)

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    # if the data does not have raw, names or clumn names
    # we have to add them, otherwise the nesting and a nesting will get confused with 
    # the link between Sample wise, columns, and Sample IDs
    if(rownames(.data) |> is.null() | colnames(.data) |> is.null() )
      warning("tidySummarizedExperiment says: the nesting and unnesting operations require row names and column names to avoid side-effects. Therefore, doors will be added as \"1\", \"2\", \"3\".")
    if(rownames(.data) |> is.null() ) rownames(.data) = .data |> nrow() |> seq_len() |> as.character()
    if(colnames(.data) |> is.null() ) colnames(.data) = .data |> ncol() |> seq_len() |> as.character()
    
    my_data__ <- .data 
    
    # Names
    sample_name <- s_(my_data__)$name
    feature_name <- f_(my_data__)$name
    sample_symbol <- s_(my_data__)$symbol
    feature_symbol <- f_(my_data__)$symbol
    
    # Check if the nesting is too complicated
    # for the moment without optimisation
    my_test_nest <- 
        my_data__[min(1, nrow(my_data__)),min(1, ncol(my_data__))] %>%
        as_tibble() %>%
        tidyr::nest(...) 

    # Understand what the nesting is about
    my_nesting_column <- my_test_nest |>
        select(-!!as.symbol(col_name_data)) |>
        colnames()

    # Check that sample or feature are in the nesting
    if(
        # Check column intersection
        c(f_(.data)$name) %>%
            intersect(colnames(my_test_nest)) %>%
            length() %>% `>` (0) &
      
        # Check that other column are there
        length(colnames(my_test_nest)) > 2
    ) {
        stop("tidySummarizedExperiment says:",
            " You cannot have the columns feature among the nesting",
            " mixed with other nesting for efficiency reasons.",
            " Please consider to convert to_tibble() first.",
            " We are working for optimising a generalised solution of nest().")
    }

    # If I nest only for .feature -> THIS WORKS ONLY WITH THE CHECK ABOVE
    if (feature_name %in% colnames(my_test_nest)) {
        return(
            my_data__ %>%
            # This is needed otherwise nest goes into loop and fails
            as_tibble() %>%
            tidyr::nest(...) %>% 
            mutate(
                !!as.symbol(col_name_data) := 
                 split_SummarizedExperiment_by_feature_to_list(!!.data)
            ) %>%
            # Coerce to tidySummarizedExperiment_nested for unnesting
            add_class("tidySummarizedExperiment_nested")
        )
    } 
    

    my_data__ %>%
        select(!!sample_symbol, !!feature_symbol, all_of(my_nesting_column)) |> 
        as_tibble() %>%
        tidyr::nest(...) |> 

          mutate(
            !!as.symbol(col_name_data) := pmap(

              # Add sample feature to map if nesting by those
                list(!!as.symbol(col_name_data)) %>%

                  # Check if nested by sample
                  when(sample_name %in% colnames(my_test_nest) ~ c(., list(!!sample_symbol)), ~ (.)) %>%

                  # Check if nested by feature
                  when(feature_name %in% colnames(my_test_nest) ~ c(., list(!!feature_symbol)), ~ (.)) , ~ { 
                    
                    # VERY COMPLICATE WAY TO DO THIS. SIMPLIFY IN THE FUTURE
                    
                    # Check if nested by sample
                    if(sample_name %in% colnames(my_test_nest)) { my_samples=..2 }
                    
                    # Here I am filtering because if I have 0 samples this leads to failure
                    else my_samples= ..1 |> filter(!is.na(!!sample_symbol)) |> pull(!!sample_symbol)

                    # Check if nested by feature and sample
                    if(sample_name %in% colnames(my_test_nest) & feature_name %in% colnames(my_test_nest)) {my_transcripts=..3}
                    else if(feature_name %in% colnames(my_test_nest)) my_transcripts=..2
                    
                    # Here I am filtering because if I have 0 features this leads to failure
                    else my_transcripts= ..1 |> filter(!is.na(!!feature_symbol)) |>  pull(!!feature_symbol)
                    
                    # if the summarised experiment does not have feature ID or Sample ID, 
                    # convert back those to integers
                    # Note to self: this is convoluted, because if feature ID and sample ID 
                    # were to stay integers all along, these will not be needed
                    if(rownames(my_data__) |> is.null()) my_transcripts = as.integer(my_transcripts)
                    if(colnames(my_data__) |> is.null()) my_samples = as.integer(my_samples)
                    
                    my_data__[unique(my_transcripts),unique(my_samples)] |>
                      select(-any_of(
                        my_nesting_column |> 
                          setdiff(c(sample_name, feature_name))
                      )) 
                 
                }
            )
        ) %>%
        # Coerce to tidySummarizedExperiment_nested for unnesting
        add_class("tidySummarizedExperiment_nested")
}

#' @name extract
#' @rdname extract
#' @inherit tidyr::extract
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' tidySummarizedExperiment::pasilla |>
#'     extract(type, into="sequencing", regex="([a-z]*)_end", convert=TRUE)
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom tidyr extract
#' @importFrom rlang enquo
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
extract.SummarizedExperiment <- function(data, col,
    into, regex="([[:alnum:]]+)", remove=TRUE,
    convert=FALSE, ...) {

    . <- NULL
    se <- tidySummarizedExperiment::se
    col <- enquo(col)

    # Deprecation of special column names
    if (is_sample_feature_deprecated_used(
        data, c(quo_name(col), into)
    )) {
        data <- ping_old_special_column_into_metadata(data)
    }

    
    special_columns <- get_special_columns(  
        # Decrease the size of the dataset
        data[1:min(100, nrow(data)), 1:min(20, ncol(data))]
    ) |> 
        c(get_needed_columns(data))
    
    tst <- intersect(quo_names(into),  special_columns) %>%
        length() %>%
        gt(0) & remove


    if (tst) {
        columns <- special_columns |>  paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }

    # Subset column annotation
    if (
        (
            all(quo_names(col) %in% colnames(colData(data))) |
            (quo_name(col) == s_(se)$name & !remove)
        ) & 
        !s_(se)$name %in% into
    ) {
        colData(data) <- colData(data) %>% 
            as.data.frame() %>% 
            as_tibble(rownames=s_(data)$name) %>% 
            tidyr::extract(col=!!col, into=into, regex=regex,
                remove=remove, convert=convert, ...) %>% 
            data.frame(row.names=pull(., !!s_(se)$symbol),
                check.names=FALSE) %>%
            select(-!!s_(se)$symbol) %>%
            DataFrame(check.names=FALSE)
      
        return(data)
    }
     
    # Subset row annotation
    if (
        (
            all( quo_names(col) %in% colnames(rowData(data))) | 
            (quo_name(col) == f_(se)$name & !remove)
        ) & 
        !f_(se)$name %in% into
    ) {
        rowData(data) <- rowData(data) %>% 
            as.data.frame() %>% 
            as_tibble(rownames=f_(data)$name) %>% 
            tidyr::extract(col=!!col, into=into, regex=regex,
                remove=remove, convert=convert, ...) %>% 
            data.frame(row.names=pull(., !!f_(se)$symbol),
                check.names=FALSE) %>%
            select(-!!f_(se)$symbol) %>%
            DataFrame(check.names=FALSE)
      
        return(data)
    }
    
    data %>%
        as_tibble(skip_GRanges=TRUE) %>%
        tidyr::extract(col=!!col, into=into, regex=regex,
            remove=remove, convert=convert, ...) %>%
        update_SE_from_tibble(data)
}

#' @name pivot_longer
#' @rdname pivot_longer
#' @inherit tidyr::pivot_longer
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' # See vignette("pivot") for examples and explanation
#' library(dplyr)
#' tidySummarizedExperiment::pasilla %>%
#'     pivot_longer(c(condition, type),
#'         names_to="name", values_to="value")
#' 
#' @importFrom tidyr pivot_longer
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
pivot_longer.SummarizedExperiment <- function(data,
    cols, ..., cols_vary = "fastest", names_to = "name", 
    names_prefix = NULL, names_sep = NULL, names_pattern = NULL, 
    names_ptypes = NULL, names_transform = NULL, names_repair = "check_unique", 
    values_to = "value", values_drop_na = FALSE, values_ptypes = NULL, 
    values_transform = NULL) {

    cols <- enquo(cols)

    message(data_frame_returned_message)

    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
        data, 
        c(quo_names(cols))
    )) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    data %>%
        as_tibble(skip_GRanges=TRUE) %>%
        tidyr::pivot_longer(!!cols, ..., cols_vary=cols_vary, 
            names_to=names_to, names_prefix=names_prefix, 
            names_sep=names_sep, names_pattern=names_pattern, 
            names_ptypes=names_ptypes, names_transform=names_transform,
            names_repair=names_repair, values_to=values_to,
            values_drop_na=values_drop_na, values_ptypes=values_ptypes, 
            values_transform=values_transform
        )
}

#' @name pivot_wider
#' @rdname pivot_wider
#' @inherit tidyr::pivot_wider
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' # See vignette("pivot") for examples and explanation
#' library(dplyr)
#' tidySummarizedExperiment::pasilla %>%
#'     pivot_wider(names_from=feature, values_from=counts)
#' 
#' @importFrom tidyr pivot_wider
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
pivot_wider.SummarizedExperiment <- function(data,
    ..., id_cols = NULL, id_expand = FALSE, names_from = name,
    names_prefix = "", names_sep = "_", names_glue = NULL,
    names_sort = FALSE, names_vary = "fastest", names_expand = FALSE,
    names_repair = "check_unique", values_from = value,
    values_fill = NULL, values_fn = NULL, unused_fn = NULL) {
    id_cols <- enquo(id_cols)
    name <- enquo(names_from)
    value <- enquo(values_from)

    message(data_frame_returned_message)

    # Deprecation of special column names
    if (is_sample_feature_deprecated_used(
        data, 
        c(quo_names(id_cols), quo_names(name), quo_names(value))
    )) {
        data <- ping_old_special_column_into_metadata(data)
    }
  
    data %>%
        as_tibble(skip_GRanges=TRUE) %>%
        tidyr::pivot_wider(..., id_cols=!!id_cols, id_expand=id_expand,
            names_from=!!name, names_prefix=names_prefix,
            names_sep=names_sep, names_glue=names_glue,
            names_sort=names_sort, names_vary=names_vary,
            names_expand=names_expand, names_repair=names_repair,
            values_from=!!value, values_fill=values_fill,
            values_fn=values_fn, unused_fn=unused_fn
        )
}

#' @name unite
#' @rdname unite
#' @inherit tidyr::unite
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' tidySummarizedExperiment::pasilla |>
#'     unite("group", c(condition, type))
#'     
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom rlang enquo enquos quo_name
#' @importFrom tidyr unite
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
unite.SummarizedExperiment <- function(data, col, ...,
    sep="_", remove=TRUE, na.rm=FALSE) {
    
    se <- tidySummarizedExperiment::se
    . <- NULL
    
    # Check that we are not modifying a key column
    cols <- enquo(col)

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
  
    special_columns <- get_special_columns(
        # Decrease the size of the dataset
        data[1:min(100, nrow(data)), 1:min(20, ncol(data))]
    ) |> 
        c(get_needed_columns(data))
    
    tst <-
        intersect(
            cols %>% quo_names(),
            special_columns
        ) %>%
        length() %>%
        gt(0) & remove

    if (tst) {
        columns <- special_columns |>  paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }

    columns_to_unite <- data[1,1] %>% select(...) %>% colnames()
    
    # Subset column annotation
    if (all(columns_to_unite %in% colnames(colData(data))) &
        !s_(se)$name %in% col) {
        colData(data) <-
            colData(data) %>% 
            as.data.frame() %>% 
            as_tibble(rownames=s_(data)$name) %>% 
            tidyr::unite(!!cols, ..., sep=sep,
                remove=remove, na.rm=na.rm) %>%
            data.frame(row.names=pull(., !!s_(se)$symbol),
                check.names=FALSE) %>%
            select(-!!s_(se)$symbol) %>%
            DataFrame(check.names=FALSE)

        return(data)
    }
    
    # Subset row annotation
    if (all(columns_to_unite %in% colnames(rowData(data))) &
        !f_(se)$name %in% col) {
        rowData(data) <-
            rowData(data) %>% 
            as.data.frame() %>% 
            as_tibble(rownames=f_(data)$name) %>% 
            tidyr::unite(!!cols, ..., sep=sep,
                remove=remove, na.rm=na.rm) %>%
            data.frame(row.names=pull(., !!f_(se)$symbol),
                check.names=FALSE) %>%
            select(-!!f_(se)$symbol) %>%
            DataFrame(check.names = FALSE)

        return(data)
    }

    # Otherwise go simple and slow
    data %>%
        as_tibble(skip_GRanges=TRUE) %>%
        tidyr::unite(!!cols, ..., sep=sep, remove=remove, na.rm=na.rm) %>%
        update_SE_from_tibble(data)
}

#' @name separate
#' @rdname separate
#' @inherit tidyr::separate
#' @return `tidySummarizedExperiment`
#' 
#' @examples
#' un <- tidySummarizedExperiment::pasilla |>
#'     unite("group", c(condition, type))
#' un |> separate(col=group, into=c("condition", "type"))
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom tidyr separate
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H., Vaughan, D. (2023). tidyr: Tidy Messy Data. R package version 2.0.0, https://CRAN.R-project.org/package=tidyr
#' @export
separate.SummarizedExperiment <- function(data, col,
    into, sep="[^[:alnum:]]+", remove=TRUE,
    convert=FALSE, extra="warn", fill="warn", ...) {

    # Fix NOTEs
    . <- NULL
    se <- tidySummarizedExperiment::se
  
    # Check that we are not modifying a key column
    cols <- enquo(col)

    # Deprecation of special column names
    if (is_sample_feature_deprecated_used(
        data, 
        c(quo_names(cols))
    )) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    special_columns <- get_special_columns(
        # Decrease the size of the dataset
        data[1:min(100, nrow(data)), 1:min(20, ncol(data))]
    ) |> 
        c(get_needed_columns(data))
    
    tst <-
        intersect(
            cols %>% quo_names(),
            special_columns
        ) %>%
        length() %>%
        gt(0) & remove
 
    if (tst) {
        columns <- special_columns |>  paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says:",
            " you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one."
        )
    }


    columns_to_unite <- data[1,1] %>% select(!!cols) %>%
        suppressMessages() %>% colnames()
    
    # Subset column annotation
    if (all(columns_to_unite %in% colnames(colData(data))) &
        (!s_(data)$name %in% into)) {
        colData(data) <-
            colData(data) %>% 
            as.data.frame() %>% 
            as_tibble(rownames=s_(data)$name) %>% 
            tidyr::separate(!!cols, into=into, sep=sep,
                remove=remove, convert=convert,
                extra=extra, fill=fill, ...) %>%
            data.frame(row.names=pull(., !!s_(se)$symbol),
                check.names=FALSE) %>%
            select(-!!s_(se)$symbol) %>%
            DataFrame(check.names=FALSE)
      
        return(data)
    }
    
    # Subset row annotation
    if (all(columns_to_unite %in% colnames(rowData(data))) &
        (!f_(data)$name %in% into)) {
        rowData(data) <-
            rowData(data) %>% 
            as.data.frame() %>% 
            as_tibble(rownames = f_(data)$name) %>% 
            tidyr::separate(!!cols, into=into, sep=sep,
                remove=remove, convert=convert,
                extra=extra, fill=fill, ...) %>%
            data.frame(row.names=pull(., !!f_(se)$symbol),
                check.names=FALSE) %>%
            select(-!!f_(data)$symbol) %>%
            DataFrame(check.names=FALSE)

        return(data)
    }
    
    # Otherwise go simple and slow
    data %>%
        as_tibble(skip_GRanges=TRUE) %>%
        tidyr::separate(!!cols, into=into, sep=sep,
            remove=remove, convert=convert,
            extra=extra, fill=fill, ...) %>%
        update_SE_from_tibble(data)
}
